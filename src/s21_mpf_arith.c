#include <stdlib.h>
#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

/* ============================================================
   Сложение / вычитание: raw
   ============================================================ */

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
    uint64_t *b_shifted = calloc(count, sizeof(uint64_t));
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
    free(b_shifted);

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
  uint64_t *b_shifted = calloc(count, sizeof(uint64_t));
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
  free(b_shifted);

  res->exp = a->exp;
  res->sign = a->sign;
  res->kind = S21_MPF_NORMAL;
  s21_mpf_normalize(res);
}

void s21_mpf_sub_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  /* x - y = x + (-y). Устраняет дублирование ветки разных знаков. */
  s21_mpf_t neg_y;
  s21_mpf_init2(&neg_y, y->prec);
  s21_mpf_neg_raw(&neg_y, y);
  s21_mpf_add_raw(res, x, &neg_y);
  s21_mpf_clear(&neg_y);
}

void s21_mpf_neg_raw(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_set(res, x);
  if (res->kind != S21_MPF_NAN) res->sign = !res->sign;
}

/* ============================================================
   Публичные add / sub
   ============================================================ */

/* Возвращает 1, если результат уже определён (спецслучай), 0 — надо считать. */
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
    s21_mpf_t xw, yw, rw;
    s21_mpf_init2(&xw, wp);
    s21_mpf_init2(&yw, wp);
    s21_mpf_init2(&rw, wp);
    s21_mpf_set(&xw, x);
    s21_mpf_set(&yw, y);

    s21_mpf_add_raw(&rw, &xw, &yw);
    s21_mpf_set(res, &rw);

    s21_mpf_clear(&xw);
    s21_mpf_clear(&yw);
    s21_mpf_clear(&rw);
  }
  return 0;
}

int s21_mpf_sub(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (s21_mpf_check_binary(res, x, y)) return -1;

  s21_mpf_t neg_y;
  s21_mpf_init2(&neg_y, y->prec);
  s21_mpf_neg_raw(&neg_y, y);
  int rc = s21_mpf_add(res, x, &neg_y);
  s21_mpf_clear(&neg_y);
  return rc;
}

/* ============================================================
   Умножение: raw + публичная
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

int s21_mpf_mul_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_INF || y->kind == S21_MPF_INF) {
    int other_zero = (x->kind == S21_MPF_INF && y->kind == S21_MPF_ZERO) ||
                     (y->kind == S21_MPF_INF && x->kind == S21_MPF_ZERO);
    if (other_zero)
      s21_mpf_set_nan(res);
    else
      s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO || y->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign ^ y->sign);
    return 0;
  }

  size_t count = s21_mpf_limbs_for_prec(res->prec);
  uint64_t *prod = malloc(2 * count * sizeof(uint64_t));
  s21_mpf_mul_mant(prod, x->limbs, y->limbs, count);

  uint32_t word_shift = res->prec / 64;
  uint32_t bit_shift = res->prec % 64;

  for (size_t i = 0; i < count; i++) {
    int src = (int)i + (int)word_shift;
    if (src >= (int)(2 * count)) {
      res->limbs[i] = 0;
    } else {
      uint64_t v = prod[src];
      if (bit_shift && src + 1 < (int)(2 * count))
        v = (v >> bit_shift) | (prod[src + 1] << (64 - bit_shift));
      else if (bit_shift)
        v >>= bit_shift;
      res->limbs[i] = v;
    }
  }
  free(prod);

  uint32_t mask_bits = res->prec % 64;
  if (mask_bits != 0) {
    uint64_t mask = (1ULL << mask_bits) - 1;
    res->limbs[count - 1] &= mask;
  }

  res->exp = x->exp + y->exp;
  res->sign = x->sign ^ y->sign;
  res->kind = S21_MPF_NORMAL;
  s21_mpf_normalize(res);
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
    /* Работаем в точности wp = 2 * prec результата. */
    uint32_t wp = 2 * res->prec;
    s21_mpf_t xw, yw;
    s21_mpf_init2(&xw, wp);
    s21_mpf_init2(&yw, wp);
    s21_mpf_set(&xw, x);
    s21_mpf_set(&yw, y);

    size_t wcount = s21_mpf_limbs_for_prec(wp);
    uint64_t *prod = malloc(2 * wcount * sizeof(uint64_t));
    s21_mpf_mul_mant(prod, xw.limbs, yw.limbs, wcount);

    /* temp.prec равен реальной точности произведения 2 * wp. */
    uint32_t prod_prec = 2 * wp;
    s21_mpf_t temp;
    s21_mpf_init2(&temp, prod_prec);
    memcpy(temp.limbs, prod,
           s21_mpf_limbs_for_prec(prod_prec) * sizeof(uint64_t));
    free(prod);

    temp.exp = xw.exp + yw.exp;
    temp.sign = xw.sign ^ yw.sign;
    temp.kind = S21_MPF_NORMAL;
    s21_mpf_normalize(&temp);

    s21_mpf_set(res, &temp);

    s21_mpf_clear(&temp);
    s21_mpf_clear(&xw);
    s21_mpf_clear(&yw);
  }
  return 0;
}

/* ============================================================
   Деление: длинное деление побитово
   ============================================================ */

/* rem <<= 1 (len лимбов, перенос между лимбами внутри len). */
static void rem_shift_left_one(uint64_t *rem, size_t len) {
  uint64_t carry = 0;
  for (size_t w = 0; w < len; w++) {
    uint64_t nc = rem[w] >> 63;
    rem[w] = (rem[w] << 1) | carry;
    carry = nc;
  }
}

/* rem >= y? У rem на один лимб больше, чем у y; старший лимб rem —
   перенос из предыдущего сдвига. */
static int rem_ge_y(const uint64_t *rem, size_t y_len, const uint64_t *y) {
  if (rem[y_len] != 0) return 1;
  for (int w = (int)y_len - 1; w >= 0; w--) {
    if (rem[w] < y[w]) return 0;
    if (rem[w] > y[w]) return 1;
  }
  return 1;
}

/* rem -= y (y_len лимбов; borrow выходит в rem[y_len]). */
static void rem_sub_y(uint64_t *rem, size_t y_len, const uint64_t *y) {
  uint64_t borrow = 0;
  for (size_t w = 0; w < y_len; w++) {
    uint64_t diff;
    uint64_t b1 = __builtin_sub_overflow(rem[w], y[w], &diff);
    uint64_t b2 = __builtin_sub_overflow(diff, borrow, &diff);
    rem[w] = diff;
    borrow = b1 | b2;
  }
  rem[y_len] -= borrow;
}

/* Установлен ли в q бит prec? Если да, частное >= 2^prec и нужна
   нормализация сдвигом на 1 бит вправо с инкрементом экспоненты. */
static int quotient_has_overflow(const uint64_t *q, size_t q_len,
                                 uint32_t prec) {
  uint32_t pw = prec / 64;
  uint32_t pb = prec % 64;
  if (pw >= q_len) return 0;
  return (int)((q[pw] >> pb) & 1ULL);
}

/* q >>= 1 (перенос между лимбами внутри len). */
static void quotient_shift_right_one(uint64_t *q, size_t len) {
  for (size_t w = 0; w < len; w++) {
    uint64_t high_bit = (w + 1 < len) ? (q[w + 1] << 63) : 0;
    q[w] = (q[w] >> 1) | high_bit;
  }
}

int s21_mpf_div_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  /* Операнды NORMAL — спецзначения обработаны в div_special. */
  uint32_t prec = res->prec;
  size_t count = s21_mpf_limbs_for_prec(prec);
  size_t rem_len = count + 1;

  uint64_t *rem = calloc(rem_len, sizeof(uint64_t));
  uint64_t *q = calloc(rem_len, sizeof(uint64_t));

  /* Длинное деление: dividend = x.mant · 2^prec (2·prec бит, MSB-первым),
     divisor = y.mant (prec бит). Частное — не более prec+1 бит. */
  for (int i = 2 * (int)prec - 1; i >= 0; i--) {
    rem_shift_left_one(rem, rem_len);
    if (i >= (int)prec) {
      int src = i - (int)prec;
      rem[0] |= (x->limbs[src / 64] >> (src % 64)) & 1ULL;
    }
    if (rem_ge_y(rem, count, y->limbs)) {
      rem_sub_y(rem, count, y->limbs);
      q[i / 64] |= 1ULL << (i % 64);
    }
  }
  free(rem);

  int overflow = quotient_has_overflow(q, rem_len, prec);
  if (overflow) quotient_shift_right_one(q, rem_len);

  memcpy(res->limbs, q, count * sizeof(uint64_t));
  free(q);

  uint32_t mask_bits = prec % 64;
  if (mask_bits != 0) res->limbs[count - 1] &= (1ULL << mask_bits) - 1;

  res->exp = x->exp - y->exp + (overflow ? 1 : 0);
  res->sign = x->sign ^ y->sign;
  res->kind = S21_MPF_NORMAL;
  s21_mpf_normalize(res);
  return 0;
}

/* Возвращает 1, если результат уже определён (спецслучай), 0 — надо считать.
   Покрывает все пары x, y, где хотя бы один не NORMAL. */
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
    s21_mpf_t xw, yw, rw;
    s21_mpf_init2(&xw, wp);
    s21_mpf_init2(&yw, wp);
    s21_mpf_init2(&rw, wp);
    s21_mpf_set(&xw, x);
    s21_mpf_set(&yw, y);

    s21_mpf_div_raw(&rw, &xw, &yw);
    s21_mpf_set(res, &rw);

    s21_mpf_clear(&xw);
    s21_mpf_clear(&yw);
    s21_mpf_clear(&rw);
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
  /* abs(-0) = +0, abs(-inf) = +inf, abs(-x) = +x. NaN остаётся NaN. */
  res->sign = 0;
}
