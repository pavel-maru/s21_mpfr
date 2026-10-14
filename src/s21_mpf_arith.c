#include <stdlib.h>
#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

/* ============================================================
   Сложение / вычитание: raw
   ============================================================ */

/* Вспомогательная: буфер для сдвинутого операнда.
   Если count лимбов влезает в S21_MPF_STACK_LIMBS — используем
   стек, иначе calloc. Возвращает указатель и признак владения. */
static uint64_t *shifted_alloc(size_t count, uint64_t *stack_buf,
                               int *owned) {
  if (count <= S21_MPF_STACK_LIMBS) {
    memset(stack_buf, 0, count * sizeof(uint64_t));
    *owned = 0;
    return stack_buf;
  }
  *owned = 1;
  return calloc(count, sizeof(uint64_t));
}

void s21_mpf_add_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->sign == y->sign) {
    size_t count = s21_mpf_limbs_for_prec(res->prec);
    const s21_mpf_t *a = x, *b = y;
    if (a->exp < b->exp) {
      const s21_mpf_t *t = a;
      a = b;
      b = t;
    }

    int64_t shift = a->exp - b->exp;
    uint64_t stack_buf[S21_MPF_STACK_LIMBS];
    int owned = 0;
    uint64_t *b_shifted = shifted_alloc(count, stack_buf, &owned);
    if (shift < (int64_t)count * 64)
      s21_mpf_shift_right_into(b_shifted, count, b->limbs, count, (int)shift);

    uint64_t carry = 0;
    for (size_t i = 0; i < count; i++) {
      uint64_t sum;
      uint64_t c1 = __builtin_add_overflow(a->limbs[i], b_shifted[i], &sum);
      uint64_t c2 = __builtin_add_overflow(sum, carry, &sum);
      res->limbs[i] = sum;
      carry = c1 | c2;
    }
    if (owned) free(b_shifted);

    res->exp = a->exp;
    res->sign = x->sign;
    res->kind = S21_MPF_NORMAL;

    if (carry) {
      uint64_t bit_in = carry;
      for (int i = (int)count - 1; i >= 0; i--) {
        uint64_t nb = res->limbs[i] & 1ULL;
        res->limbs[i] = (res->limbs[i] >> 1) | (bit_in << 63);
        bit_in = nb;
      }
      res->exp += 1;
    }
    s21_mpf_normalize(res);
    return;
  }

  int cmp = s21_mpf_cmp_abs(x, y);
  if (cmp == 0) {
    s21_mpf_set_zero(res, 0);
    return;
  }

  const s21_mpf_t *a = (cmp > 0) ? x : y;
  const s21_mpf_t *b = (cmp > 0) ? y : x;
  size_t count = s21_mpf_limbs_for_prec(res->prec);

  int64_t shift = a->exp - b->exp;
  uint64_t stack_buf[S21_MPF_STACK_LIMBS];
  int owned = 0;
  uint64_t *b_shifted = shifted_alloc(count, stack_buf, &owned);
  if (shift < (int64_t)count * 64)
    s21_mpf_shift_right_into(b_shifted, count, b->limbs, count, (int)shift);

  uint64_t borrow = 0;
  for (size_t i = 0; i < count; i++) {
    uint64_t diff;
    uint64_t b1 = __builtin_sub_overflow(a->limbs[i], b_shifted[i], &diff);
    uint64_t b2 = __builtin_sub_overflow(diff, borrow, &diff);
    res->limbs[i] = diff;
    borrow = b1 | b2;
  }
  if (owned) free(b_shifted);

  res->exp = a->exp;
  res->sign = a->sign;
  res->kind = S21_MPF_NORMAL;
  s21_mpf_normalize(res);
}

void s21_mpf_sub_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  /* x - y = x + (-y). */
  s21_mpf_stack_t neg_y_s;
  s21_mpf_stack_init(&neg_y_s, y->prec);
  s21_mpf_neg_raw(&neg_y_s.mpf, y);
  s21_mpf_add_raw(res, x, &neg_y_s.mpf);
  s21_mpf_stack_clear(&neg_y_s);
}

void s21_mpf_neg_raw(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_set(res, x);
  if (res->kind != S21_MPF_NAN) res->sign = !res->sign;
}

/* ============================================================
   Публичные add / sub
   ============================================================ */

static int add_special(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 1;
  }
  if (x->kind == S21_MPF_INF && y->kind == S21_MPF_INF) {
    if (x->sign == y->sign)
      s21_mpf_set_inf(res, x->sign);
    else
      s21_mpf_set_nan(res);
    return 1;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, x->sign);
    return 1;
  }
  if (y->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, y->sign);
    return 1;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set(res, y);
    return 1;
  }
  if (y->kind == S21_MPF_ZERO) {
    s21_mpf_set(res, x);
    return 1;
  }
  return 0;
}

int s21_mpf_add(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (s21_mpf_check_binary(res, x, y)) return -1;

  if (!add_special(res, x, y)) {
    uint32_t wp = res->prec + 64;
    s21_mpf_stack_t xw_s, yw_s, rw_s;
    s21_mpf_stack_init(&xw_s, wp);
    s21_mpf_stack_init(&yw_s, wp);
    s21_mpf_stack_init(&rw_s, wp);
    s21_mpf_set(&xw_s.mpf, x);
    s21_mpf_set(&yw_s.mpf, y);

    s21_mpf_add_raw(&rw_s.mpf, &xw_s.mpf, &yw_s.mpf);
    s21_mpf_set(res, &rw_s.mpf);

    s21_mpf_stack_clear(&xw_s);
    s21_mpf_stack_clear(&yw_s);
    s21_mpf_stack_clear(&rw_s);
  }
  return 0;
}

int s21_mpf_sub(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (s21_mpf_check_binary(res, x, y)) return -1;

  s21_mpf_stack_t neg_y_s;
  s21_mpf_stack_init(&neg_y_s, y->prec);
  s21_mpf_neg_raw(&neg_y_s.mpf, y);
  int rc = s21_mpf_add(res, x, &neg_y_s.mpf);
  s21_mpf_stack_clear(&neg_y_s);
  return rc;
}

/* ============================================================
   Умножение с округлением RNDN
   ============================================================ */

/* ============================================================
   Хелперы для Toom-Cook 3-way.

   Все операции беззнаковые; промежуточные значения гарантированно
   неотрицательны (проверено алгебраически на формулах интерполяции).
   ============================================================ */

/* s[0..len) *= m, m ∈ {1, 2, 3}.  Возвращает carry — верхний
   лимб, который caller должен записать в s[len] (гарантируя
   что там есть место).  Для m == 1 carry всегда 0. */
static uint64_t mulsmall_carry(uint64_t *s, size_t len, uint64_t m) {
  if (m == 1) return 0;
  uint64_t carry = 0;
  for (size_t i = 0; i < len; i++) {
    unsigned __int128 p = (unsigned __int128)s[i] * m + carry;
    s[i] = (uint64_t)p;
    carry = (uint64_t)(p >> 64);
  }
  return carry;
}

/* s[0..s_len) += a[0..a_len).  Carry распространяется вверх,
   обрезается на s_len. */
static void add_into(uint64_t *s, size_t s_len, const uint64_t *a,
                     size_t a_len) {
  uint64_t carry = 0;
  size_t i = 0;
  for (; i < a_len && i < s_len; i++) {
    uint64_t t1, t2, c1, c2;
    c1 = __builtin_add_overflow(s[i], a[i], &t1);
    c2 = __builtin_add_overflow(t1, carry, &t2);
    s[i] = t2;
    carry = c1 | c2;
  }
  while (carry && i < s_len) {
    s[i]++;
    if (s[i] != 0) carry = 0;
    i++;
  }
}

/* s[0..len) <<= shift, shift ∈ [1, 63].  Верхние shift бит
   теряются; caller гарантирует что в них нули. */
static void lshift_inplace(uint64_t *s, size_t len, int shift) {
  if (shift <= 0) return;
  int bshift = shift;
  for (int i = (int)len - 1; i > 0; i--) {
    s[i] = (s[i] << bshift) | (s[i - 1] >> (64 - bshift));
  }
  s[0] <<= bshift;
}

/* s[0..len) /= m, m ∈ {2, 3, 6}.  Деление точное. */
static void div_small_raw(uint64_t *s, size_t len, uint64_t m) {
  uint64_t rem = 0;
  for (int i = (int)len - 1; i >= 0; i--) {
    unsigned __int128 v = ((unsigned __int128)rem << 64) | s[i];
    s[i] = (uint64_t)(v / m);
    rem = (uint64_t)(v % m);
  }
}

/* s[0..k+2) = a0 + m·a1 + m²·a2, m ∈ {1, 2, 3}.
   a0, a1, a2 — по k лимбов.  Результат занимает k+1 лимбов
   (старший может быть ненулевым), s[k+1] = 0.

   Horner: s = a2; s = s·m + a1; s = s·m + a0. */
static void eval_poly(uint64_t *s, const uint64_t *a0, const uint64_t *a1,
                      const uint64_t *a2, size_t k, uint64_t m) {
  memcpy(s, a2, k * sizeof(uint64_t));
  s[k]     = mulsmall_carry(s, k, m);
  s[k + 1] = 0;
  add_into(s, k + 2, a1, k);
  s[k + 1] += mulsmall_carry(s, k + 1, m);
  add_into(s, k + 2, a0, k);
}

/* ------------------------------------------------------------
   Умножение мантисс: школьный алгоритм как base case для
   Карацубы.

   prod = a * b, a и b по n лимбов, результат — 2n лимбов.
   prod должен быть уже выделен под 2n лимбов.
   ------------------------------------------------------------ */

void s21_mpf_mul_mant_schoolbook(uint64_t *prod, const uint64_t *a,
                                 const uint64_t *b, size_t n) {
  memset(prod, 0, 2 * n * sizeof(uint64_t));

  for (size_t i = 0; i < n; i++) {
    uint64_t carry = 0;
    for (size_t j = 0; j < n; j++) {
      unsigned __int128 p =
          (unsigned __int128)a[i] * b[j] + prod[i + j] + carry;
      prod[i + j] = (uint64_t)p;
      carry = (uint64_t)(p >> 64);
    }
    size_t k = i + n;
    while (carry != 0 && k < 2 * n) {
      unsigned __int128 s = (unsigned __int128)prod[k] + carry;
      prod[k] = (uint64_t)s;
      carry = (uint64_t)(s >> 64);
      k++;
    }
  }
}

/* ------------------------------------------------------------
   Карацуба.

   Для a = a0 + a1·B^h, b = b0 + b1·B^h, где h = n/2 (floor),
   B = 2^64:
     z0 = a0·b0
     z2 = a1·b1
     z1 = (a0 + a1)·(b0 + b1) - z0 - z2
     a·b = z0 + z1·B^h + z2·B^(2h)

   Три подпроизведения размера ~n/2 вместо четырёх в schoolbook.
   Рекурсия идёт до S21_MPF_KARATSUBA_THRESHOLD, ниже — schoolbook.

   scratch — единый пул для текущего уровня и всей рекурсии.
   На каждом уровне фиксированный layout:

     [0 .. 2h)              z0
     [2h .. 2h + 2m)        z2
     [2h+2m .. 2h+2m+2m+2)  z1 = s·t
     [.. + m+1)             s = a0 + a1
     [.. + m+1)             t = b0 + b1
     [..  ]                 sub — вложенная рекурсия

   где h = n/2 (floor), m = n - h.  Верхняя оценка полного
   размера scratch: ~8n лимбов, с запасом аллоцируем 10n + 256.
   ------------------------------------------------------------ */

/* dst (len limbs) += src (src_len limbs) на позиции offset.
   Carry распространяется вверх по dst. */
static void add_shifted(uint64_t *dst, const uint64_t *src, size_t src_len,
                        size_t offset, size_t dst_len) {
  uint64_t carry = 0;
  size_t i = 0;
  for (; i < src_len && offset + i < dst_len; i++) {
    uint64_t t1, t2, c1, c2;
    c1 = __builtin_add_overflow(dst[offset + i], src[i], &t1);
    c2 = __builtin_add_overflow(t1, carry, &t2);
    dst[offset + i] = t2;
    carry = c1 | c2;
  }
  size_t j = offset + i;
  while (carry && j < dst_len) {
    dst[j]++;
    if (dst[j] != 0) carry = 0;
    j++;
  }
}

/* a (a_len limbs) -= b (b_len limbs, дополняется нулями до a_len).
   Borrow распространяется вверх. */
static void sub_in_place_padded(uint64_t *a, const uint64_t *b, size_t b_len,
                                size_t a_len) {
  uint64_t borrow = 0;
  for (size_t i = 0; i < a_len; i++) {
    uint64_t bb = (i < b_len) ? b[i] : 0;
    uint64_t t1, t2, b1, b2;
    b1 = __builtin_sub_overflow(a[i], bb, &t1);
    b2 = __builtin_sub_overflow(t1, borrow, &t2);
    a[i] = t2;
    borrow = b1 | b2;
  }
}

/* s (m+1 limbs) = a0 (h limbs) + a1 (m limbs).  h <= m. */
static void add_split(uint64_t *s, const uint64_t *a0, size_t h,
                      const uint64_t *a1, size_t m) {
  uint64_t carry = 0;
  for (size_t i = 0; i < h; i++) {
    uint64_t t1, t2, c1, c2;
    c1 = __builtin_add_overflow(a0[i], a1[i], &t1);
    c2 = __builtin_add_overflow(t1, carry, &t2);
    s[i] = t2;
    carry = c1 | c2;
  }
  for (size_t i = h; i < m; i++) {
    uint64_t t1, c1;
    c1 = __builtin_add_overflow(a1[i], carry, &t1);
    s[i] = t1;
    carry = c1;
  }
  s[m] = carry;
}

static void mul_mant_rec(uint64_t *prod, const uint64_t *a, const uint64_t *b,
                         size_t n, uint64_t *scratch) {
  if (n < S21_MPF_KARATSUBA_THRESHOLD) {
    s21_mpf_mul_mant_schoolbook(prod, a, b, n);
    return;
  }

  size_t h = n / 2;
  size_t m = n - h; /* m == h или m == h+1 */

  uint64_t *z0  = scratch;
  uint64_t *z2  = z0 + 2 * h;
  uint64_t *z1  = z2 + 2 * m;
  uint64_t *s   = z1 + 2 * m + 2;
  uint64_t *t   = s + (m + 1);
  uint64_t *sub = t + (m + 1);

  add_split(s, a, h, a + h, m);
  add_split(t, b, h, b + h, m);

  mul_mant_rec(z0, a, b, h, sub);
  mul_mant_rec(z2, a + h, b + h, m, sub);
  mul_mant_rec(z1, s, t, m + 1, sub);

  /* z1 -= z0, z1 -= z2 (z0, z2 дополняются нулями до 2m+2). */
  sub_in_place_padded(z1, z0, 2 * h, 2 * m + 2);
  sub_in_place_padded(z1, z2, 2 * m, 2 * m + 2);

  /* Сборка: prod = z0 + z1·B^h + z2·B^(2h). */
  memset(prod, 0, 2 * n * sizeof(uint64_t));
  add_shifted(prod, z0, 2 * h, 0, 2 * n);
  add_shifted(prod, z1, 2 * m + 2, h, 2 * n);
  add_shifted(prod, z2, 2 * m, 2 * h, 2 * n);
}

/* Точка входа.  Экспонирована в заголовке как _karatsuba — имя
   подчёркивает, что ниже порога функция автоматически переходит
   на schoolbook, а вызов из production-кода (mul_raw) идёт
   именно сюда. */
void s21_mpf_mul_mant_karatsuba(uint64_t *prod, const uint64_t *a,
                                const uint64_t *b, size_t n) {
  if (n < S21_MPF_KARATSUBA_THRESHOLD) {
    s21_mpf_mul_mant_schoolbook(prod, a, b, n);
    return;
  }
  if (n >= S21_MPF_TOOMCOOK_THRESHOLD && n % 3 == 0) {
    s21_mpf_mul_mant_toomcook(prod, a, b, n);
    return;
  }
  /* Карацуба. */
  size_t scratch_len = 10 * n + 256;
  uint64_t *scratch = malloc(scratch_len * sizeof(uint64_t));
  if (scratch == NULL) {
    s21_mpf_mul_mant_schoolbook(prod, a, b, n);
    return;
  }
  mul_mant_rec(prod, a, b, n, scratch);
  free(scratch);
}

/* ============================================================
   Toom-Cook 3-way (Stage 2: n кратно 3).

   Для a = a0 + a1·B^k + a2·B^(2k), b = b0 + b1·B^k + b2·B^(2k),
   B = 2^64, k = n/3:

   C(x) = a(x)·b(x), где a(x) = a0 + a1·x + a2·x².
   Вычисляем C в точках x ∈ {0, 1, 2, 3, ∞}:
     w0 = C(0)      = a0·b0
     w1 = C(1)      = (a0+a1+a2)·(b0+b1+b2)
     w2 = C(2)      = (a0+2a1+4a2)·(b0+2b1+4b2)
     w3 = C(3)      = (a0+3a1+9a2)·(b0+3b1+9b2)
     w4 = C(∞)      = a2·b2

   Коэффициенты c0..c4 восстанавливаются формулами, использующими
   только сложение, вычитание и деление на 2/3/6.  Промежуточные
   значения неотрицательны.

   Scratch layout на уровень:

     wa1, wa2, wa3:  3·(k+2)
     wb1, wb2, wb3:  3·(k+2)
     w0, w1, w2, w3, w4:  5·(2k+4)
     tmp:            2k+4

   Итого 18k + 36 ≈ 6n + 36 лимбов.  Округляем до 6n + 64. */
static void mul_mant_toom_rec(uint64_t *prod, const uint64_t *a,
                              const uint64_t *b, size_t n, uint64_t *scratch) {
  size_t k = n / 3;
  size_t L = 2 * k + 2;

  uint64_t *wa1 = scratch;
  uint64_t *wa2 = wa1 + (k + 2);
  uint64_t *wa3 = wa2 + (k + 2);
  uint64_t *wb1 = wa3 + (k + 2);
  uint64_t *wb2 = wb1 + (k + 2);
  uint64_t *wb3 = wb2 + (k + 2);
  uint64_t *w0  = wb3 + (k + 2);
  uint64_t *w1  = w0 + (2 * k + 4);
  uint64_t *w2  = w1 + (2 * k + 4);
  uint64_t *w3  = w2 + (2 * k + 4);
  uint64_t *w4  = w3 + (2 * k + 4);
  uint64_t *tmp = w4 + (2 * k + 4);

  /* Оценки в точках 1, 2, 3.  Результат k+1 лимбов. */
  eval_poly(wa1, a, a + k, a + 2 * k, k, 1);
  eval_poly(wa2, a, a + k, a + 2 * k, k, 2);
  eval_poly(wa3, a, a + k, a + 2 * k, k, 3);
  eval_poly(wb1, b, b + k, b + 2 * k, k, 1);
  eval_poly(wb2, b, b + k, b + 2 * k, k, 2);
  eval_poly(wb3, b, b + k, b + 2 * k, k, 3);

  /* Рекурсивные умножения.  w0, w4 — 2k лимбов; w1, w2, w3 — 2k+2. */
  s21_mpf_mul_mant_karatsuba(w0, a, b, k);
  memset(w0 + 2 * k, 0, 4 * sizeof(uint64_t));
  s21_mpf_mul_mant_karatsuba(w4, a + 2 * k, b + 2 * k, k);
  memset(w4 + 2 * k, 0, 4 * sizeof(uint64_t));
  s21_mpf_mul_mant_karatsuba(w1, wa1, wb1, k + 1);
  s21_mpf_mul_mant_karatsuba(w2, wa2, wb2, k + 1);
  s21_mpf_mul_mant_karatsuba(w3, wa3, wb3, k + 1);

  /* u = w1 − w0 − w4 в w1. */
  sub_in_place_padded(w1, w0, 2 * k, L);
  sub_in_place_padded(w1, w4, 2 * k, L);

  /* v = w2 − w0 − 16·w4 в w2. */
  sub_in_place_padded(w2, w0, 2 * k, L);
  memcpy(tmp, w4, 2 * k * sizeof(uint64_t));
  tmp[2 * k]     = 0;
  tmp[2 * k + 1] = 0;
  lshift_inplace(tmp, L, 4);         /* tmp = 16·w4 */
  sub_in_place_padded(w2, tmp, L, L);

  /* w = w3 − w0 − 81·w4 в w3. */
  sub_in_place_padded(w3, w0, 2 * k, L);
  memcpy(tmp, w4, 2 * k * sizeof(uint64_t));
  tmp[2 * k]     = mulsmall_carry(tmp, 2 * k, 81);
  tmp[2 * k + 1] = 0;
  sub_in_place_padded(w3, tmp, L, L);

  /* c3 = (w + 3u − 3v) / 6 в w3.  w3 (бывшее w) + 3·w1 (u) − 3·w2 (v). */
  memcpy(tmp, w1, L * sizeof(uint64_t));
  mulsmall_carry(tmp, L - 1, 3);      /* tmp = 3u (L-1 лимбов) */
  /* tmp занимает L лимбов максимум, L-й уже 0.  Дополним до L. */
  /* Здесь tmp[L-1] уже 0, что корректно. */
  add_into(w3, L, tmp, L);
  memcpy(tmp, w2, L * sizeof(uint64_t));
  mulsmall_carry(tmp, L - 1, 3);
  sub_in_place_padded(w3, tmp, L, L);
  div_small_raw(w3, L, 6);            /* w3 = c3 */

  /* c2 = (v − 2u − 6·c3) / 2 в w2. */
  memcpy(tmp, w1, L * sizeof(uint64_t));
  mulsmall_carry(tmp, L - 1, 2);
  sub_in_place_padded(w2, tmp, L, L);
  memcpy(tmp, w3, L * sizeof(uint64_t));
  mulsmall_carry(tmp, L - 1, 6);
  sub_in_place_padded(w2, tmp, L, L);
  div_small_raw(w2, L, 2);            /* w2 = c2 */

  /* c1 = u − c2 − c3 в w1. */
  sub_in_place_padded(w1, w2, L, L);
  sub_in_place_padded(w1, w3, L, L);

  /* Реконструкция: prod = c0 + c1·B^k + c2·B^(2k) + c3·B^(3k) + c4·B^(4k). */
  memset(prod, 0, 6 * k * sizeof(uint64_t));
  add_shifted(prod, w0, L, 0,       6 * k);
  add_shifted(prod, w1, L, k,       6 * k);
  add_shifted(prod, w2, L, 2 * k,   6 * k);
  add_shifted(prod, w3, L, 3 * k,   6 * k);
  add_shifted(prod, w4, L, 4 * k,   6 * k);
}

/* Публичная точка входа Toom-Cook.  Если n % 3 != 0 или
   n < S21_MPF_TOOMCOOK_THRESHOLD — откатываемся на Карацубу. */
void s21_mpf_mul_mant_toomcook(uint64_t *prod, const uint64_t *a,
                               const uint64_t *b, size_t n) {
  if (n < S21_MPF_TOOMCOOK_THRESHOLD || n % 3 != 0) {
    s21_mpf_mul_mant_karatsuba(prod, a, b, n);
    return;
  }
  size_t scratch_len = 6 * n + 64;
  uint64_t *scratch = malloc(scratch_len * sizeof(uint64_t));
  if (scratch == NULL) {
    s21_mpf_mul_mant_schoolbook(prod, a, b, n);
    return;
  }
  mul_mant_toom_rec(prod, a, b, n, scratch);
  free(scratch);
}

/* Операнды NORMAL, все одинаковой точности. Внутри — рабочие
   guard bits wp = res->prec + 64; произведение умещается в 2*wp
   бит, финальное округление до res->prec — RNDN через set_round.

   xw, yw, temp размещаются в едином скрэтч-буфере:
     - при малых prec (total <= S21_MPF_MUL_SCRATCH_LIMBS) — стек;
     - иначе — один calloc вместо трёх init2 + отдельного prod.
   temp.limbs служит и приёмником mul_mant, отдельный prod
   не нужен. */
int s21_mpf_mul_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  uint32_t wp = res->prec + 64;
  uint32_t prod_prec = 2 * wp;
  size_t l_wp = s21_mpf_limbs_for_prec(wp);
  size_t l_prod = s21_mpf_limbs_for_prec(prod_prec);
  size_t total = 2 * l_wp + l_prod;

  uint64_t stack_buf[S21_MPF_MUL_SCRATCH_LIMBS];
  uint64_t *base;
  int on_heap = 0;
  if (total <= S21_MPF_MUL_SCRATCH_LIMBS) {
    memset(stack_buf, 0, total * sizeof(uint64_t));
    base = stack_buf;
  } else {
    base = calloc(total, sizeof(uint64_t));
    if (base == NULL) return -1;
    on_heap = 1;
  }

  /* Собираем три mpf вручную: limbs указывают внутрь base,
     s21_mpf_init2/clear не вызываются — освобождать нечего. */
  s21_mpf_t xw = { .limbs = base,
                   .exp   = 0,
                   .prec  = wp,
                   .sign  = 0,
                   .kind  = S21_MPF_ZERO };
  s21_mpf_t yw = { .limbs = base + l_wp,
                   .exp   = 0,
                   .prec  = wp,
                   .sign  = 0,
                   .kind  = S21_MPF_ZERO };
  s21_mpf_t temp = { .limbs = base + 2 * l_wp,
                     .exp   = 0,
                     .prec  = prod_prec,
                     .sign  = 0,
                     .kind  = S21_MPF_ZERO };

  s21_mpf_set(&xw, x);
  s21_mpf_set(&yw, y);

  s21_mpf_mul_mant_karatsuba(temp.limbs, xw.limbs, yw.limbs, l_wp);

  temp.exp  = xw.exp + yw.exp;
  temp.sign = xw.sign ^ yw.sign;
  temp.kind = S21_MPF_NORMAL;
  s21_mpf_normalize(&temp);

  s21_mpf_set(res, &temp);

  if (on_heap) free(base);
  return 0;
}

static int mul_special(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 1;
  }
  if (x->kind == S21_MPF_INF || y->kind == S21_MPF_INF) {
    int other_zero = (x->kind == S21_MPF_INF && y->kind == S21_MPF_ZERO) ||
                     (y->kind == S21_MPF_INF && x->kind == S21_MPF_ZERO);
    if (other_zero)
      s21_mpf_set_nan(res);
    else
      s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 1;
  }
  if (x->kind == S21_MPF_ZERO || y->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign ^ y->sign);
    return 1;
  }
  return 0;
}

int s21_mpf_mul(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (s21_mpf_check_binary(res, x, y)) return -1;

  int rc = 0;
  if (!mul_special(res, x, y)) {
    rc = s21_mpf_mul_raw(res, x, y);
  }
  return rc;
}

/* ============================================================
   Деление: алгоритм D Кнута (TAOCP 4.3.1)
   ============================================================ */

static uint32_t div_norm_shift(uint32_t prec, size_t n) {
  uint32_t msb_in_top = prec - 1 - (uint32_t)(64 * (n - 1));
  return 63 - msb_in_top;
}

static void div_shift_left_into(uint64_t *V, const uint64_t *y_limbs, size_t n,
                                uint32_t s) {
  if (s == 0) {
    memcpy(V, y_limbs, n * sizeof(uint64_t));
    return;
  }
  V[0] = y_limbs[0] << s;
  for (size_t i = 1; i < n; i++) {
    V[i] = (y_limbs[i] << s) | (y_limbs[i - 1] >> (64 - s));
  }
}

static void div_prepare_dividend(uint64_t *U, const uint64_t *x_limbs,
                                 size_t n) {
  memcpy(U + n, x_limbs, n * sizeof(uint64_t));
}

static void div_knuth(uint64_t *U, const uint64_t *V, size_t n, uint64_t *Q) {
  for (int j = (int)n; j >= 0; j--) {
    unsigned __int128 num = ((unsigned __int128)U[j + n] << 64) | U[j + n - 1];
    unsigned __int128 qhat_wide = num / V[n - 1];
    uint64_t qhat = (qhat_wide > UINT64_MAX) ? UINT64_MAX : (uint64_t)qhat_wide;
    unsigned __int128 rhat = num - (unsigned __int128)qhat * V[n - 1];

    if (n >= 2) {
      while ((unsigned __int128)qhat * V[n - 2] > (rhat << 64) + U[j + n - 2]) {
        qhat--;
        rhat += V[n - 1];
        if (rhat >> 64) break;
      }
    }

    uint64_t borrow = 0;
    uint64_t carry = 0;
    for (size_t i = 0; i < n; i++) {
      unsigned __int128 p = (unsigned __int128)qhat * V[i] + carry;
      carry = (uint64_t)(p >> 64);
      uint64_t p_lo = (uint64_t)p;

      uint64_t t1, t2, c1, c2;
      c1 = __builtin_sub_overflow(U[j + i], p_lo, &t1);
      c2 = __builtin_sub_overflow(t1, borrow, &t2);
      U[j + i] = t2;
      borrow = c1 | c2;
    }

    uint64_t t1, t2, c1, c2;
    c1 = __builtin_sub_overflow(U[j + n], carry, &t1);
    c2 = __builtin_sub_overflow(t1, borrow, &t2);
    U[j + n] = t2;
    uint64_t final_borrow = c1 | c2;

    if (final_borrow) {
      qhat--;
      uint64_t carry2 = 0;
      for (size_t i = 0; i < n; i++) {
        uint64_t s1, s2, cb1, cb2;
        cb1 = __builtin_add_overflow(U[j + i], V[i], &s1);
        cb2 = __builtin_add_overflow(s1, carry2, &s2);
        U[j + i] = s2;
        carry2 = cb1 | cb2;
      }
      U[j + n] += carry2;
    }

    Q[j] = qhat;
  }
}

static int quotient_has_overflow(const uint64_t *Q, size_t q_len,
                                 uint32_t prec) {
  uint32_t pw = prec / 64;
  uint32_t pb = prec % 64;
  if (pw >= q_len) return 0;
  return (int)((Q[pw] >> pb) & 1ULL);
}

static void quotient_shift_right_one(uint64_t *Q, size_t len) {
  for (size_t w = 0; w < len; w++) {
    uint64_t high = (w + 1 < len) ? (Q[w + 1] << 63) : 0;
    Q[w] = (Q[w] >> 1) | high;
  }
}

/* Операнды NORMAL — спецзначения обрабатываются в div_special.

   Рабочие буферы U (2n+1), V (n), Q (n+1) — в одном стековом
   массиве, если 4n+2 <= 4 * S21_MPF_STACK_LIMBS. Иначе calloc.
   Зануляется только используемая часть (total лимбов), не весь
   стековый буфер. */
int s21_mpf_div_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  uint32_t prec = res->prec;
  size_t n = s21_mpf_limbs_for_prec(prec);

  size_t u_len = 2 * n + 1;
  size_t v_len = n;
  size_t q_len = n + 1;
  size_t total = u_len + v_len + q_len;

  uint64_t stack_buf[4 * S21_MPF_STACK_LIMBS];
  int owned = 0;
  uint64_t *base;
  if (total <= 4 * S21_MPF_STACK_LIMBS) {
    memset(stack_buf, 0, total * sizeof(uint64_t));
    base = stack_buf;
  } else {
    base = calloc(total, sizeof(uint64_t));
    owned = 1;
  }
  uint64_t *U = base;
  uint64_t *V = base + u_len;
  uint64_t *Q = base + u_len + v_len;

  uint32_t s = div_norm_shift(prec, n);
  div_shift_left_into(V, y->limbs, n, s);
  div_prepare_dividend(U, x->limbs, n);
  div_knuth(U, V, n, Q);

  int overflow = quotient_has_overflow(Q, q_len, prec);
  if (overflow) quotient_shift_right_one(Q, q_len);

  memcpy(res->limbs, Q, n * sizeof(uint64_t));
  if (owned) free(base);

  uint32_t mask_bits = prec % 64;
  if (mask_bits != 0) res->limbs[n - 1] &= (1ULL << mask_bits) - 1;

  res->exp = x->exp - y->exp + (overflow ? 1 : 0);
  res->sign = x->sign ^ y->sign;
  res->kind = S21_MPF_NORMAL;
  s21_mpf_normalize(res);
  return 0;
}

static int div_special(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 1;
  }
  if (x->kind == S21_MPF_INF && y->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 1;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 1;
  }
  if (y->kind == S21_MPF_INF) {
    s21_mpf_set_zero(res, x->sign ^ y->sign);
    return 1;
  }
  if (y->kind == S21_MPF_ZERO) {
    if (x->kind == S21_MPF_ZERO)
      s21_mpf_set_nan(res);
    else
      s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 1;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign ^ y->sign);
    return 1;
  }
  return 0;
}

int s21_mpf_div(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (s21_mpf_check_binary(res, x, y)) return -1;

  int rc = 0;
  if (!div_special(res, x, y)) {
    uint32_t wp = res->prec + 64;
    s21_mpf_stack_t xw_s, yw_s, rw_s;
    s21_mpf_stack_init(&xw_s, wp);
    s21_mpf_stack_init(&yw_s, wp);
    s21_mpf_stack_init(&rw_s, wp);
    s21_mpf_set(&xw_s.mpf, x);
    s21_mpf_set(&yw_s.mpf, y);

    rc = s21_mpf_div_raw(&rw_s.mpf, &xw_s.mpf, &yw_s.mpf);
    if (rc == 0) s21_mpf_set(res, &rw_s.mpf);

    s21_mpf_stack_clear(&xw_s);
    s21_mpf_stack_clear(&yw_s);
    s21_mpf_stack_clear(&rw_s);
  }
  return rc;
}

void s21_mpf_div_small(s21_mpf_t *x, uint32_t n) {
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  uint64_t rem = 0;
  for (int i = (int)count - 1; i >= 0; i--) {
    unsigned __int128 v = ((unsigned __int128)rem << 64) | x->limbs[i];
    x->limbs[i] = (uint64_t)(v / (uint32_t)n);
    rem = (uint64_t)(v % (uint32_t)n);
  }
  s21_mpf_normalize(x);
}

/* ============================================================
   Публичные neg / abs
   ============================================================ */

void s21_mpf_neg(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return;
  s21_mpf_neg_raw(res, x);
}

void s21_mpf_abs(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return;
  s21_mpf_set(res, x);
  res->sign = 0;
}
