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

static void s21_mpf_mul_mant(uint64_t *prod, const uint64_t *a,
                             const uint64_t *b, size_t count) {
  memset(prod, 0, 2 * count * sizeof(uint64_t));

  for (size_t i = 0; i < count; i++) {
    uint64_t carry = 0;
    for (size_t j = 0; j < count; j++) {
      unsigned __int128 p =
          (unsigned __int128)a[i] * b[j] + prod[i + j] + carry;
      prod[i + j] = (uint64_t)p;
      carry = (uint64_t)(p >> 64);
    }
    size_t k = i + count;
    while (carry != 0 && k < 2 * count) {
      unsigned __int128 s = (unsigned __int128)prod[k] + carry;
      prod[k] = (uint64_t)s;
      carry = (uint64_t)(s >> 64);
      k++;
    }
  }
}

/* Операнды NORMAL, все одинаковой точности. Внутри — рабочие
   guard bits wp = res->prec + 64; произведение умещается в 2*wp
   бит, финальное округление до res->prec — RNDN через set_round.

   NOTE: здесь пока не применена схема со стековыми mpf: temp
   требует 2*wp лимбов, что не влезает в S21_MPF_STACK_LIMBS при
   больших prec. Отдельная задача. */
int s21_mpf_mul_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  uint32_t wp = res->prec + 64;

  s21_mpf_t xw, yw, temp;
  s21_mpf_init2(&xw, wp);
  s21_mpf_init2(&yw, wp);
  s21_mpf_init2(&temp, 2 * wp);

  s21_mpf_set(&xw, x);
  s21_mpf_set(&yw, y);

  size_t wcount = s21_mpf_limbs_for_prec(wp);
  uint64_t *prod = malloc(2 * wcount * sizeof(uint64_t));
  s21_mpf_mul_mant(prod, xw.limbs, yw.limbs, wcount);
  memcpy(temp.limbs, prod, s21_mpf_limbs_for_prec(2 * wp) * sizeof(uint64_t));
  free(prod);

  temp.exp = xw.exp + yw.exp;
  temp.sign = xw.sign ^ yw.sign;
  temp.kind = S21_MPF_NORMAL;
  s21_mpf_normalize(&temp);

  s21_mpf_set(res, &temp);

  s21_mpf_clear(&temp);
  s21_mpf_clear(&xw);
  s21_mpf_clear(&yw);
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

  if (!mul_special(res, x, y)) {
    s21_mpf_mul_raw(res, x, y);
  }
  return 0;
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

  if (!div_special(res, x, y)) {
    uint32_t wp = res->prec + 64;
    s21_mpf_stack_t xw_s, yw_s, rw_s;
    s21_mpf_stack_init(&xw_s, wp);
    s21_mpf_stack_init(&yw_s, wp);
    s21_mpf_stack_init(&rw_s, wp);
    s21_mpf_set(&xw_s.mpf, x);
    s21_mpf_set(&yw_s.mpf, y);

    s21_mpf_div_raw(&rw_s.mpf, &xw_s.mpf, &yw_s.mpf);
    s21_mpf_set(res, &rw_s.mpf);

    s21_mpf_stack_clear(&xw_s);
    s21_mpf_stack_clear(&yw_s);
    s21_mpf_stack_clear(&rw_s);
  }
  return 0;
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
