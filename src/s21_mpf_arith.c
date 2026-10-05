#include "s21_mpf.h"
#include "s21_mpf_internal.h"

#include <stdlib.h>
#include <string.h>

/* ============================================================
   Сдвиг вправо (внутренний, только для add_raw/sub_raw)
   ============================================================ */

static void s21_mpf_shift_right_into(uint64_t *dst, const uint64_t *src,
                                     size_t count, int shift) {
  if (shift == 0) {
    memcpy(dst, src, count * sizeof(uint64_t));
    return;
  }
  int word_shift = shift / 64;
  int bit_shift = shift % 64;

  for (size_t i = 0; i < count; i++) {
    int src_idx = (int)i + word_shift;
    if (src_idx >= (int)count) {
      dst[i] = 0;
    } else {
      uint64_t v = src[src_idx];
      if (bit_shift && src_idx + 1 < (int)count) {
        v = (v >> bit_shift) | (src[src_idx + 1] << (64 - bit_shift));
      } else if (bit_shift) {
        v >>= bit_shift;
      }
      dst[i] = v;
    }
  }
}

/* ============================================================
   Сложение / вычитание: raw
   ============================================================ */

void s21_mpf_add_raw(s21_mpf_t *res, const s21_mpf_t *x,
                     const s21_mpf_t *y) {
  if (x->sign == y->sign) {
    size_t count = s21_mpf_limbs_for_prec(res->prec);
    const s21_mpf_t *a = x, *b = y;
    if (a->exp < b->exp) { const s21_mpf_t *t = a; a = b; b = t; }

    int64_t shift = a->exp - b->exp;
    uint64_t *b_shifted = calloc(count, sizeof(uint64_t));
    if (shift < (int64_t)count * 64) {
      s21_mpf_shift_right_into(b_shifted, b->limbs, count, (int)shift);
    }

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
  } else {
    int cmp = s21_mpf_cmp_abs(x, y);
    if (cmp == 0) {
      s21_mpf_set_zero(res, 0);
    } else {
      const s21_mpf_t *a = (cmp > 0) ? x : y;
      const s21_mpf_t *b = (cmp > 0) ? y : x;

      size_t count = s21_mpf_limbs_for_prec(res->prec);
      int64_t shift = a->exp - b->exp;
      uint64_t *b_shifted = calloc(count, sizeof(uint64_t));
      if (shift < (int64_t)count * 64) {
        s21_mpf_shift_right_into(b_shifted, b->limbs, count, (int)shift);
      }

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
  }
}

void s21_mpf_sub_raw(s21_mpf_t *res, const s21_mpf_t *x,
                     const s21_mpf_t *y) {
  const s21_mpf_t *a = x, *b = y;
  if (a->exp < b->exp) { const s21_mpf_t *t = a; a = b; b = t; }
  if (a->exp == b->exp && s21_mpf_cmp_abs(a, b) < 0) {
    const s21_mpf_t *t = a; a = b; b = t;
  }

  int64_t shift = a->exp - b->exp;
  size_t count = s21_mpf_limbs_for_prec(a->prec);

  uint64_t *b_shifted = calloc(count, sizeof(uint64_t));
  if (shift < (int64_t)count * 64) {
    s21_mpf_shift_right_into(b_shifted, b->limbs, count, (int)shift);
  }

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
  res->sign = 0;
  res->kind = S21_MPF_NORMAL;
  s21_mpf_normalize(res);
}

void s21_mpf_neg_raw(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_set(res, x);
  if (res->kind != S21_MPF_NAN) res->sign = !res->sign;
}

/* ============================================================
   Публичные add / sub
   ============================================================ */

int s21_mpf_add(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;
  if (res->prec != x->prec || x->prec != y->prec) return -1;

  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF && y->kind == S21_MPF_INF) {
    if (x->sign == y->sign) s21_mpf_set_inf(res, x->sign);
    else s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_INF) { s21_mpf_set_inf(res, x->sign); return 0; }
  if (y->kind == S21_MPF_INF) { s21_mpf_set_inf(res, y->sign); return 0; }
  if (x->kind == S21_MPF_ZERO) { s21_mpf_set(res, y); return 0; }
  if (y->kind == S21_MPF_ZERO) { s21_mpf_set(res, x); return 0; }

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
  return 0;
}

int s21_mpf_sub(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;
  if (res->prec != x->prec || x->prec != y->prec) return -1;

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

int s21_mpf_mul_raw(s21_mpf_t *res, const s21_mpf_t *x,
                    const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF || y->kind == S21_MPF_INF) {
    int other_zero = (x->kind == S21_MPF_INF && y->kind == S21_MPF_ZERO) ||
                     (y->kind == S21_MPF_INF && x->kind == S21_MPF_ZERO);
    if (other_zero) s21_mpf_set_nan(res);
    else s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO || y->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign ^ y->sign); return 0;
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
      if (bit_shift && src + 1 < (int)(2 * count)) {
        v = (v >> bit_shift) | (prod[src + 1] << (64 - bit_shift));
      } else if (bit_shift) {
        v >>= bit_shift;
      }
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

int s21_mpf_mul(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;
  if (res->prec != x->prec || x->prec != y->prec) return -1;

  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF || y->kind == S21_MPF_INF) {
    int other_zero = (x->kind == S21_MPF_INF && y->kind == S21_MPF_ZERO) ||
                     (y->kind == S21_MPF_INF && x->kind == S21_MPF_ZERO);
    if (other_zero) s21_mpf_set_nan(res);
    else s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO || y->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign ^ y->sign); return 0;
  }

  /* Работаем в точности wp = 2 * prec результата: мантиссы
     операндов поднимаются в wp бит, произведение умещается
     в 2 * wp бит, старшие биты после приведения дают точность
     prec с запасом на округление. */
  uint32_t wp = 2 * res->prec;
  s21_mpf_t xw, yw;
  s21_mpf_init2(&xw, wp);
  s21_mpf_init2(&yw, wp);
  s21_mpf_set(&xw, x);
  s21_mpf_set(&yw, y);

  size_t wcount = s21_mpf_limbs_for_prec(wp);
  uint64_t *prod = malloc(2 * wcount * sizeof(uint64_t));
  s21_mpf_mul_mant(prod, xw.limbs, yw.limbs, wcount);

  /* Раньше temp.prec выставлялся как 2 * wcount * 64 — округление
     до границы лимба. Это совпадает с 2 * wp только когда prec
     кратно 32; иначе normalize() уводил MSB выше реального
     значения, а exp не получал поправки, и результат уезжал
     на (2 * wcount * 64 - 2 * wp) бит. Теперь temp.prec равен
     реальной точности произведения 2 * wp, копируются только
     значащие лимбы. */
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
  return 0;
}

/* ============================================================
   Деление: raw + публичная + div_small
   ============================================================ */

int s21_mpf_div_raw(s21_mpf_t *res, const s21_mpf_t *x,
                    const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF && y->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF) { s21_mpf_set_inf(res, x->sign ^ y->sign); return 0; }
  if (y->kind == S21_MPF_INF) { s21_mpf_set_zero(res, x->sign ^ y->sign); return 0; }
  if (y->kind == S21_MPF_ZERO) {
    if (x->kind == S21_MPF_ZERO) s21_mpf_set_nan(res);
    else s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) { s21_mpf_set_zero(res, x->sign ^ y->sign); return 0; }

  uint32_t prec = res->prec;
  size_t count = s21_mpf_limbs_for_prec(prec);

  uint64_t *rem = calloc(count + 1, sizeof(uint64_t));
  uint64_t *q = calloc(count + 1, sizeof(uint64_t));

  for (int i = 2 * (int)prec - 1; i >= 0; i--) {
    uint64_t carry = 0;
    for (size_t w = 0; w < count + 1; w++) {
      uint64_t nc = rem[w] >> 63;
      rem[w] = (rem[w] << 1) | carry;
      carry = nc;
    }
    if (i >= (int)prec) {
      int src = i - (int)prec;
      uint64_t b = (x->limbs[src / 64] >> (src % 64)) & 1ULL;
      rem[0] |= b;
    }

    int geq = (rem[count] != 0);
    if (!geq) {
      geq = 1;
      for (int w = (int)count - 1; w >= 0; w--) {
        if (rem[w] < y->limbs[w]) { geq = 0; break; }
        if (rem[w] > y->limbs[w]) { geq = 1; break; }
      }
    }

    if (geq) {
      uint64_t borrow = 0;
      for (size_t w = 0; w < count; w++) {
        uint64_t diff;
        uint64_t b1 = __builtin_sub_overflow(rem[w], y->limbs[w], &diff);
        uint64_t b2 = __builtin_sub_overflow(diff, borrow, &diff);
        rem[w] = diff;
        borrow = b1 | b2;
      }
      rem[count] -= borrow;
      q[i / 64] |= (1ULL << (i % 64));
    }
  }
  free(rem);

  uint32_t pw = prec / 64;
  uint32_t pb = prec % 64;
  int overflow = 0;
  if (pw < count + 1) {
    if ((q[pw] >> pb) & 1ULL) overflow = 1;
  }
  if (overflow) {
    for (size_t w = 0; w < count + 1; w++) {
      uint64_t high_bit = (w + 1 < count + 1) ? (q[w + 1] << 63) : 0;
      q[w] = (q[w] >> 1) | high_bit;
    }
  }

  memcpy(res->limbs, q, count * sizeof(uint64_t));
  free(q);

  uint32_t mask_bits = prec % 64;
  if (mask_bits != 0) {
    res->limbs[count - 1] &= (1ULL << mask_bits) - 1;
  }

  res->exp = x->exp - y->exp + (overflow ? 1 : 0);
  res->sign = x->sign ^ y->sign;
  res->kind = S21_MPF_NORMAL;
  s21_mpf_normalize(res);
  return 0;
}

int s21_mpf_div(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;
  if (res->prec != x->prec || x->prec != y->prec) return -1;

  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF || y->kind == S21_MPF_INF || y->kind == S21_MPF_ZERO ||
      x->kind == S21_MPF_ZERO) {
    return s21_mpf_div_raw(res, x, y);
  }

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
  if (res->kind == S21_MPF_NORMAL || res->kind == S21_MPF_INF) {
    res->sign = 0;
  }
}
