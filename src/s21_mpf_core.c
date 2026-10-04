#include "s21_mpf.h"

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================
   Инициализация и освобождение
   ============================================================ */

void s21_mpf_init2(s21_mpf_t *x, uint32_t prec) {
  if (prec < 2) prec = 2;
  x->prec = prec;
  x->limbs = calloc(s21_mpf_limbs_for_prec(prec), sizeof(uint64_t));
  x->exp = 0;
  x->sign = 0;
  x->kind = S21_MPF_ZERO;
}

void s21_mpf_init(s21_mpf_t *x) { s21_mpf_init2(x, 256); }

void s21_mpf_clear(s21_mpf_t *x) {
  if (x->limbs) {
    free(x->limbs);
    x->limbs = NULL;
  }
  x->prec = 0;
  x->exp = 0;
  x->sign = 0;
  x->kind = S21_MPF_ZERO;
}

void s21_mpf_set_prec(s21_mpf_t *x, uint32_t prec) {
  if (prec < 2) prec = 2;
  if (prec == x->prec) return;

  uint64_t *new_limbs = calloc(s21_mpf_limbs_for_prec(prec), sizeof(uint64_t));

  size_t old_count = s21_mpf_limbs_for_prec(x->prec);
  size_t new_count = s21_mpf_limbs_for_prec(prec);
  size_t n = old_count < new_count ? old_count : new_count;
  memcpy(new_limbs, x->limbs, n * sizeof(uint64_t));

  free(x->limbs);
  x->limbs = new_limbs;
  x->prec = prec;
}

/* ============================================================
   Низкоуровневые операции над битами
   ============================================================ */

int s21_mpf_get_bit(const s21_mpf_t *x, uint32_t pos) {
  if (pos >= x->prec) return 0;
  uint32_t word = pos / 64;
  uint32_t bit = pos % 64;
  return (int)((x->limbs[word] >> bit) & 1ULL);
}

void s21_mpf_set_bit(s21_mpf_t *x, uint32_t pos, int value) {
  if (x == NULL || pos >= x->prec) return;
  uint32_t word = pos / 64;
  uint32_t bit = pos % 64;
  if (value) {
    x->limbs[word] |= (1ULL << bit);
  } else {
    x->limbs[word] &= ~(1ULL << bit);
  }
}

int s21_mpf_msb(const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) return -1;
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    if (x->limbs[i] != 0) {
      return i * 64 + (63 - __builtin_clzll(x->limbs[i]));
    }
  }
  return -1;
}

int s21_mpf_lsb(const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) return -1;
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (size_t i = 0; i < count; i++) {
    if (x->limbs[i] != 0) {
      return (int)(i * 64 + __builtin_ctzll(x->limbs[i]));
    }
  }
  return -1;
}

void s21_mpf_shift_left_into(uint64_t *dst, const uint64_t *src,
                             size_t count, int shift) {
  if (shift == 0) {
    memcpy(dst, src, count * sizeof(uint64_t));
    return;
  }
  int word_shift = shift / 64;
  int bit_shift = shift % 64;

  for (int i = (int)count - 1; i >= 0; i--) {
    int src_idx = i - word_shift;
    if (src_idx < 0) {
      dst[i] = 0;
    } else {
      uint64_t v = src[src_idx];
      if (bit_shift && src_idx > 0) {
        v = (v << bit_shift) | (src[src_idx - 1] >> (64 - bit_shift));
      } else if (bit_shift) {
        v <<= bit_shift;
      }
      dst[i] = v;
    }
  }
}

/* ============================================================
   Нормализация мантиссы
   ============================================================ */

void s21_mpf_normalize(s21_mpf_t *x) {
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  int last = (int)count - 1;
  while (last >= 0 && x->limbs[last] == 0) last--;

  if (last < 0) {
    x->kind = S21_MPF_ZERO;
    x->exp = 0;
    return;
  }

  int top_bit_global = last * 64 + (63 - __builtin_clzll(x->limbs[last]));
  int shift = (int)(x->prec - 1) - top_bit_global;

  if (shift > 0) {
    int word_shift = shift / 64;
    int bit_shift = shift % 64;

    for (int i = (int)count - 1; i >= 0; i--) {
      int src = i - word_shift;
      if (src < 0) {
        x->limbs[i] = 0;
      } else {
        uint64_t v = x->limbs[src];
        if (bit_shift && src > 0) {
          v = (v << bit_shift) | (x->limbs[src - 1] >> (64 - bit_shift));
        } else if (bit_shift) {
          v <<= bit_shift;
        }
        x->limbs[i] = v;
      }
    }
    x->exp -= shift;
  } else if (shift < 0) {
    int rshift = -shift;
    int word_shift = rshift / 64;
    int bit_shift = rshift % 64;

    for (size_t i = 0; i < count; i++) {
      int src = (int)i + word_shift;
      if (src >= (int)count) {
        x->limbs[i] = 0;
      } else {
        uint64_t v = x->limbs[src];
        if (bit_shift && src + 1 < (int)count) {
          v = (v >> bit_shift) | (x->limbs[src + 1] << (64 - bit_shift));
        } else if (bit_shift) {
          v >>= bit_shift;
        }
        x->limbs[i] = v;
      }
    }
    x->exp += rshift;
  }

  uint32_t mask_bits = x->prec % 64;
  if (mask_bits != 0) {
    uint64_t mask = (1ULL << mask_bits) - 1;
    x->limbs[count - 1] &= mask;
  }

  x->kind = S21_MPF_NORMAL;
}

/* ============================================================
   Присваивание
   ============================================================ */

void s21_mpf_set_zero(s21_mpf_t *x, int sign) {
  memset(x->limbs, 0, s21_mpf_limbs_for_prec(x->prec) * sizeof(uint64_t));
  x->exp = 0;
  x->sign = sign ? 1 : 0;
  x->kind = S21_MPF_ZERO;
}

void s21_mpf_set_nan(s21_mpf_t *x) {
  memset(x->limbs, 0, s21_mpf_limbs_for_prec(x->prec) * sizeof(uint64_t));
  x->exp = 0;
  x->sign = 0;
  x->kind = S21_MPF_NAN;
}

void s21_mpf_set_inf(s21_mpf_t *x, int sign) {
  memset(x->limbs, 0, s21_mpf_limbs_for_prec(x->prec) * sizeof(uint64_t));
  x->exp = 0;
  x->sign = sign ? 1 : 0;
  x->kind = S21_MPF_INF;
}

void s21_mpf_set_ui(s21_mpf_t *x, unsigned long v) {
  if (v == 0) {
    s21_mpf_set_zero(x, 0);
    return;
  }
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  memset(x->limbs, 0, count * sizeof(uint64_t));
  x->limbs[0] = v;
  x->sign = 0;
  x->exp = (int64_t)x->prec;
  s21_mpf_normalize(x);
}

void s21_mpf_set_si(s21_mpf_t *x, long v) {
  if (v < 0) {
    s21_mpf_set_ui(x, (unsigned long)(-(v + 1)) + 1UL);
    x->sign = 1;
  } else {
    s21_mpf_set_ui(x, (unsigned long)v);
  }
}

static void s21_decompose_double(double v, uint64_t *mant, int64_t *exp,
                                  int *sign) {
  union {
    double d;
    uint64_t u;
  } u;
  u.d = v;

  *sign = (int)(u.u >> 63);
  int64_t e_raw = (int64_t)((u.u >> 52) & 0x7FF);
  uint64_t frac = u.u & 0xFFFFFFFFFFFFFULL;

  if (e_raw == 0) {
    *mant = frac;
    *exp = 1 - 1023;
  } else {
    *mant = (1ULL << 52) | frac;
    *exp = e_raw - 1022;
  }
}

void s21_mpf_set_d(s21_mpf_t *x, double v) {
  if (v != v) { s21_mpf_set_nan(x); return; }
  if (v == 1.0 / 0.0) { s21_mpf_set_inf(x, 0); return; }
  if (v == -1.0 / 0.0) { s21_mpf_set_inf(x, 1); return; }
  if (v == 0.0) { s21_mpf_set_zero(x, (1.0 / v < 0) ? 1 : 0); return; }

  uint64_t mant_d;
  int64_t exp_d;
  int sign_d;
  s21_decompose_double(v, &mant_d, &exp_d, &sign_d);

  size_t count = s21_mpf_limbs_for_prec(x->prec);
  memset(x->limbs, 0, count * sizeof(uint64_t));

  x->limbs[0] = mant_d;
  x->sign = sign_d;
  x->exp = exp_d - 53 + (int64_t)x->prec;

  s21_mpf_normalize(x);
}

/* Присваивание по умолчанию: RNDN. */
void s21_mpf_set(s21_mpf_t *dst, const s21_mpf_t *src) {
  s21_mpf_set_round(dst, src, S21_MPF_RNDN);
}

/* Присваивание с явным режимом округления.

   При dst->prec >= src->prec операция точная (расширение).
   При dst->prec < src->prec отбрасываемые биты анализируются:
     round_bit = старший отбрасываемый бит
     sticky    = любые остальные отбрасываемые биты
     lsb       = младший бит сохранённой части

   RNDN: если round_bit=0 — не округляем; если round_bit=1 и sticky=1 —
         округляем; если round_bit=1 и sticky=0 (точная середина) —
         округляем к чётному (по lsb).
   RNDZ: всегда вниз (к нулю).
   RNDU: вверх, если результат положительный.
   RNDD: вниз, если результат отрицательный. */
int s21_mpf_set_round(s21_mpf_t *dst, const s21_mpf_t *src, s21_mpf_rnd_t rnd) {
  if (dst == NULL || src == NULL) return -1;
  if (dst == src) return 0;

  size_t dst_count = s21_mpf_limbs_for_prec(dst->prec);
  memset(dst->limbs, 0, dst_count * sizeof(uint64_t));

  /* Особые значения или совпадающая точность: копирование как есть. */
  if (src->kind != S21_MPF_NORMAL || dst->prec == src->prec) {
    size_t src_count = s21_mpf_limbs_for_prec(src->prec);
    size_t n = dst_count < src_count ? dst_count : src_count;
    memcpy(dst->limbs, src->limbs, n * sizeof(uint64_t));
    dst->exp = src->exp;
    dst->sign = src->sign;
    dst->kind = src->kind;
    return 0;
  }

  size_t src_count = s21_mpf_limbs_for_prec(src->prec);

  /* Расширение точности: сдвиг влево, значение не меняется. */
  if (dst->prec > src->prec) {
    int shift = (int)(dst->prec - src->prec);
    s21_mpf_shift_left_into(dst->limbs, src->limbs, dst_count, shift);
    dst->exp = src->exp;
    dst->sign = src->sign;
    dst->kind = S21_MPF_NORMAL;
    return 0;
  }

  /* Сужение точности: анализируем отбрасываемые биты. */
  int shift = (int)(src->prec - dst->prec);
  int word_shift = shift / 64;
  int bit_shift = shift % 64;

  /* round_bit = старший отбрасываемый бит (позиция shift - 1). */
  int round_bit = 0;
  {
    int pos = shift - 1;
    round_bit = (int)((src->limbs[pos / 64] >> (pos % 64)) & 1ULL);
  }

  /* sticky = любой ненулевой бит среди младших (shift - 1). */
  int sticky = 0;
  {
    int top = shift - 1;
    int full_words = top / 64;
    int rem_bits = top % 64;
    for (int i = 0; i < full_words; i++) {
      if (src->limbs[i] != 0) { sticky = 1; break; }
    }
    if (!sticky && rem_bits > 0) {
      uint64_t mask = (1ULL << rem_bits) - 1;
      if ((src->limbs[full_words] & mask) != 0) sticky = 1;
    }
  }

  /* Сдвиг src вправо на shift бит → dst. */
  for (size_t i = 0; i < dst_count; i++) {
    int src_idx = (int)i + word_shift;
    if (src_idx >= (int)src_count) {
      dst->limbs[i] = 0;
    } else {
      uint64_t v = src->limbs[src_idx];
      if (bit_shift && src_idx + 1 < (int)src_count) {
        v = (v >> bit_shift) | (src->limbs[src_idx + 1] << (64 - bit_shift));
      } else if (bit_shift) {
        v >>= bit_shift;
      }
      dst->limbs[i] = v;
    }
  }

  /* Отрезать биты выше dst->prec. */
  uint32_t mask_bits = dst->prec % 64;
  if (mask_bits != 0 && dst_count > 0) {
    dst->limbs[dst_count - 1] &= (1ULL << mask_bits) - 1;
  }

  dst->exp = src->exp;
  dst->sign = src->sign;
  dst->kind = S21_MPF_NORMAL;

  /* Решение об округлении вверх. */
  int round_up = 0;
  if (round_bit || sticky) {
    switch (rnd) {
      case S21_MPF_RNDZ:
        round_up = 0;
        break;
      case S21_MPF_RNDU:
        round_up = (src->sign == 0);
        break;
      case S21_MPF_RNDD:
        round_up = (src->sign == 1);
        break;
      case S21_MPF_RNDN:
      default:
        if (round_bit == 0) {
          round_up = 0;
        } else if (sticky) {
          round_up = 1;
        } else {
          round_up = (int)(dst->limbs[0] & 1ULL);
        }
        break;
    }
  }

  if (round_up) {
    uint64_t carry = 1;
    for (size_t i = 0; i < dst_count && carry; i++) {
      dst->limbs[i]++;
      if (dst->limbs[i] != 0) carry = 0;
    }
    if (carry) {
      /* Переполнение за пределы prec бит: мантисса стала 2^prec,
         что представляется как 2^(prec-1) с exp+1. */
      memset(dst->limbs, 0, dst_count * sizeof(uint64_t));
      uint32_t top = dst->prec - 1;
      dst->limbs[top / 64] = 1ULL << (top % 64);
      dst->exp += 1;
    }
  }

  return 0;
}

/* ============================================================
   Утилиты
   ============================================================ */

int s21_mpf_is_nan(const s21_mpf_t *x) { return x->kind == S21_MPF_NAN; }
int s21_mpf_is_inf(const s21_mpf_t *x) { return x->kind == S21_MPF_INF; }
int s21_mpf_is_zero(const s21_mpf_t *x) { return x->kind == S21_MPF_ZERO; }
int s21_mpf_is_normal(const s21_mpf_t *x) { return x->kind == S21_MPF_NORMAL; }

int s21_mpf_sign(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_ZERO || x->kind == S21_MPF_NAN) return 0;
  return x->sign ? -1 : 1;
}

uint32_t s21_mpf_get_prec(const s21_mpf_t *x) { return x->prec; }

/* ============================================================
   Сравнение
   ============================================================ */

static int s21_mpf_cmp_mant(const s21_mpf_t *x, const s21_mpf_t *y) {
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    if (x->limbs[i] < y->limbs[i]) return -1;
    if (x->limbs[i] > y->limbs[i]) return 1;
  }
  return 0;
}

int s21_mpf_cmp_abs(const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) return 0;

  int xinf = (x->kind == S21_MPF_INF);
  int yinf = (y->kind == S21_MPF_INF);
  if (xinf && yinf) return 0;
  if (xinf) return 1;
  if (yinf) return -1;

  int xzero = (x->kind == S21_MPF_ZERO);
  int yzero = (y->kind == S21_MPF_ZERO);
  if (xzero && yzero) return 0;
  if (xzero) return -1;
  if (yzero) return 1;

  if (x->exp < y->exp) return -1;
  if (x->exp > y->exp) return 1;

  if (x->prec == y->prec) return s21_mpf_cmp_mant(x, y);

  double xd = 0, yd = 0;
  {
    size_t xc = s21_mpf_limbs_for_prec(x->prec);
    for (int i = (int)xc - 1; i >= 0; i--)
      xd = xd * 18446744073709551616.0 + (double)x->limbs[i];
    xd = ldexp(xd, (int)(x->exp - (int64_t)x->prec));

    size_t yc = s21_mpf_limbs_for_prec(y->prec);
    for (int i = (int)yc - 1; i >= 0; i--)
      yd = yd * 18446744073709551616.0 + (double)y->limbs[i];
    yd = ldexp(yd, (int)(y->exp - (int64_t)y->prec));
  }
  if (xd < yd) return -1;
  if (xd > yd) return 1;
  return 0;
}

int s21_mpf_cmp(const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) return 0;

  int xs = x->kind == S21_MPF_ZERO ? 0 : (x->sign ? -1 : 1);
  int ys = y->kind == S21_MPF_ZERO ? 0 : (y->sign ? -1 : 1);

  if (x->kind != S21_MPF_ZERO && y->kind != S21_MPF_ZERO && xs != ys) {
    return xs < ys ? -1 : 1;
  }

  int cmp_abs = s21_mpf_cmp_abs(x, y);
  if (xs < 0) return -cmp_abs;
  return cmp_abs;
}

int s21_mpf_equal(const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) return 0;
  if (x->kind != y->kind) return 0;
  if (x->kind == S21_MPF_ZERO) return 1;
  if (x->kind == S21_MPF_INF) return x->sign == y->sign;
  if (x->sign != y->sign) return 0;
  if (x->prec != y->prec) return 0;
  if (x->exp != y->exp) return 0;
  return s21_mpf_cmp_mant(x, y) == 0;
}

int s21_mpf_zero_p(const s21_mpf_t *x) { return x->kind == S21_MPF_ZERO; }

int s21_mpf_integer_p(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_ZERO) return 1;
  if (x->kind != S21_MPF_NORMAL) return 0;
  if (x->exp >= (int64_t)x->prec) return 1;

  int64_t fractional_bits = (int64_t)x->prec - x->exp;
  if (fractional_bits <= 0) return 1;
  if (fractional_bits >= (int64_t)x->prec) return 0;

  uint32_t fb = (uint32_t)fractional_bits;
  uint32_t full_words = fb / 64;
  uint32_t rem_bits = fb % 64;

  for (uint32_t i = 0; i < full_words; i++) {
    if (x->limbs[i] != 0) return 0;
  }
  if (rem_bits) {
    uint64_t mask = (1ULL << rem_bits) - 1;
    if ((x->limbs[full_words] & mask) != 0) return 0;
  }
  return 1;
}

/* ============================================================
   Вспомогательные функции для арифметики
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

static void s21_mpf_add_magnitude(s21_mpf_t *res, const s21_mpf_t *x,
                                  const s21_mpf_t *y) {
  const s21_mpf_t *a = x, *b = y;
  if (a->exp < b->exp) {
    const s21_mpf_t *t = a;
    a = b;
    b = t;
  }

  int64_t shift = a->exp - b->exp;
  size_t count = s21_mpf_limbs_for_prec(a->prec);

  uint64_t *b_shifted = calloc(count, sizeof(uint64_t));
  if (shift < (int64_t)a->prec) {
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
  res->sign = 0;
  res->kind = S21_MPF_NORMAL;

  if (carry) {
    uint64_t bit_in = carry;
    for (int i = (int)count - 1; i >= 0; i--) {
      uint64_t new_bit_in = res->limbs[i] & 1ULL;
      res->limbs[i] = (res->limbs[i] >> 1) | (bit_in << 63);
      bit_in = new_bit_in;
    }
    res->exp += 1;
  }

  s21_mpf_normalize(res);
}

static void s21_mpf_sub_magnitude(s21_mpf_t *res, const s21_mpf_t *x,
                                  const s21_mpf_t *y) {
  const s21_mpf_t *a = x, *b = y;
  if (a->exp < b->exp) {
    const s21_mpf_t *t = a;
    a = b;
    b = t;
  }
  if (a->exp == b->exp && s21_mpf_cmp_mant(a, b) < 0) {
    const s21_mpf_t *t = a;
    a = b;
    b = t;
  }

  int64_t shift = a->exp - b->exp;
  size_t count = s21_mpf_limbs_for_prec(a->prec);

  uint64_t *b_shifted = calloc(count, sizeof(uint64_t));
  if (shift < (int64_t)a->prec) {
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

/* ============================================================
   Арифметика
   ============================================================ */

void s21_mpf_neg(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return;
  s21_mpf_set(res, x);
  if (res->kind != S21_MPF_NAN) {
    res->sign = !res->sign;
  }
}

void s21_mpf_abs(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return;
  s21_mpf_set(res, x);
  if (res->kind == S21_MPF_NORMAL || res->kind == S21_MPF_INF) {
    res->sign = 0;
  }
}

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

  s21_mpf_t tmp;
  s21_mpf_init2(&tmp, x->prec);

  if (x->sign == y->sign) {
    s21_mpf_add_magnitude(&tmp, x, y);
    tmp.sign = x->sign;
  } else {
    int cmp = s21_mpf_cmp_abs(x, y);
    if (cmp == 0) {
      s21_mpf_set_zero(&tmp, 0);
    } else if (cmp > 0) {
      s21_mpf_sub_magnitude(&tmp, x, y);
      tmp.sign = x->sign;
    } else {
      s21_mpf_sub_magnitude(&tmp, y, x);
      tmp.sign = y->sign;
    }
  }

  s21_mpf_set(res, &tmp);
  s21_mpf_clear(&tmp);
  return 0;
}

int s21_mpf_sub(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;
  if (res->prec != x->prec || x->prec != y->prec) return -1;

  s21_mpf_t neg_y;
  s21_mpf_init2(&neg_y, y->prec);
  s21_mpf_neg(&neg_y, y);
  int rc = s21_mpf_add(res, x, &neg_y);
  s21_mpf_clear(&neg_y);
  return rc;
}

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

int s21_mpf_div(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;
  if (res->prec != x->prec || x->prec != y->prec) return -1;

  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF && y->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, x->sign ^ y->sign); return 0;
  }
  if (y->kind == S21_MPF_INF) {
    s21_mpf_set_zero(res, x->sign ^ y->sign); return 0;
  }
  if (y->kind == S21_MPF_ZERO) {
    if (x->kind == S21_MPF_ZERO) s21_mpf_set_nan(res);
    else s21_mpf_set_inf(res, x->sign ^ y->sign);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign ^ y->sign); return 0;
  }

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

/* ============================================================
   Деление на маленькое целое
   ============================================================ */

static void s21_mpf_div_small(s21_mpf_t *x, uint32_t n) {
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
   Квадратный корень
   ============================================================ */

static int s21_mpf_round_to_prec(s21_mpf_t *res, const s21_mpf_t *src,
                                 uint32_t target_prec) {
  if (res->prec != target_prec) return -1;
  if (src->prec == target_prec) {
    s21_mpf_set(res, src);
    return 0;
  }
  if (src->prec < target_prec) return -1;
  s21_mpf_set(res, src);
  return 0;
}

int s21_mpf_sqrt(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;

  if (x->kind == S21_MPF_NAN) { s21_mpf_set_nan(res); return 0; }
  if (x->kind == S21_MPF_ZERO) { s21_mpf_set_zero(res, 0); return 0; }
  if (x->sign == 1) { s21_mpf_set_nan(res); return 0; }
  if (x->kind == S21_MPF_INF) { s21_mpf_set_inf(res, 0); return 0; }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, x_cur, x_next, tmp;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&x_cur, work_prec);
  s21_mpf_init2(&x_next, work_prec);
  s21_mpf_init2(&tmp, work_prec);

  s21_mpf_set(&a, x);

  int64_t half_exp = a.exp / 2;
  s21_mpf_set_ui(&x_cur, 1);
  x_cur.exp = half_exp + 1;
  s21_mpf_normalize(&x_cur);

  for (int i = 0; i < 200; i++) {
    if (s21_mpf_div(&tmp, &a, &x_cur) != 0) {
      s21_mpf_set_nan(res);
      goto cleanup;
    }
    if (s21_mpf_add(&x_next, &x_cur, &tmp) != 0) {
      s21_mpf_set_nan(res);
      goto cleanup;
    }
    x_next.exp -= 1;
    s21_mpf_normalize(&x_next);

    if (s21_mpf_cmp(&x_next, &x_cur) == 0) break;
    s21_mpf_set(&x_cur, &x_next);
  }

  if (res->prec == work_prec) s21_mpf_set(res, &x_cur);
  else s21_mpf_round_to_prec(res, &x_cur, res->prec);

cleanup:
  s21_mpf_clear(&a);
  s21_mpf_clear(&x_cur);
  s21_mpf_clear(&x_next);
  s21_mpf_clear(&tmp);
  return 0;
}

/* ============================================================
   Экспонента
   ============================================================ */

int s21_mpf_exp(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;

  if (x->kind == S21_MPF_NAN) { s21_mpf_set_nan(res); return 0; }
  if (x->kind == S21_MPF_ZERO) { s21_mpf_set_ui(res, 1); return 0; }
  if (x->kind == S21_MPF_INF) {
    if (x->sign == 0) s21_mpf_set_inf(res, 0);
    else s21_mpf_set_zero(res, 0);
    return 0;
  }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, y, sum, term;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&y, work_prec);
  s21_mpf_init2(&sum, work_prec);
  s21_mpf_init2(&term, work_prec);

  s21_mpf_set(&a, x);
  s21_mpf_set(&y, &a);

  int s = 0;
  while (y.exp > -1 && s < 128) {
    y.exp -= 1;
    s++;
  }

  s21_mpf_set_ui(&sum, 1);
  s21_mpf_set_ui(&term, 1);

  for (int n = 1; n < 300; n++) {
    s21_mpf_mul(&term, &term, &y);
    s21_mpf_div_small(&term, (uint32_t)n);
    s21_mpf_add(&sum, &sum, &term);

    if (term.kind == S21_MPF_ZERO) break;
    if (term.exp < sum.exp - (int64_t)work_prec - 4) break;
  }

  for (int i = 0; i < s; i++) {
    s21_mpf_mul(&sum, &sum, &sum);
  }

  s21_mpf_set(res, &sum);

  s21_mpf_clear(&a);
  s21_mpf_clear(&y);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&term);
  return 0;
}

/* ============================================================
   Ряд atanh(z)
   ============================================================ */

static void s21_mpf_atanh_series(s21_mpf_t *out, const s21_mpf_t *z,
                                  uint32_t work_prec) {
  s21_mpf_t z2, term, sum, tmp;
  s21_mpf_init2(&z2, work_prec);
  s21_mpf_init2(&term, work_prec);
  s21_mpf_init2(&sum, work_prec);
  s21_mpf_init2(&tmp, work_prec);

  s21_mpf_mul(&z2, z, z);
  s21_mpf_set(&term, z);
  s21_mpf_set(&sum, z);

  for (int n = 3; n < 500; n += 2) {
    s21_mpf_mul(&term, &term, &z2);
    s21_mpf_set(&tmp, &term);
    s21_mpf_div_small(&tmp, (uint32_t)n);
    s21_mpf_add(&sum, &sum, &tmp);

    if (term.kind == S21_MPF_ZERO) break;
    if (term.exp < sum.exp - (int64_t)work_prec - 4) break;
  }

  s21_mpf_set(out, &sum);
  s21_mpf_clear(&z2);
  s21_mpf_clear(&term);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&tmp);
}

/* ============================================================
   Логарифм
   ============================================================ */

int s21_mpf_log(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;

  if (x->kind == S21_MPF_NAN) { s21_mpf_set_nan(res); return 0; }
  if (x->kind == S21_MPF_ZERO) { s21_mpf_set_inf(res, 1); return 0; }
  if (x->sign == 1) { s21_mpf_set_nan(res); return 0; }
  if (x->kind == S21_MPF_INF) { s21_mpf_set_inf(res, 0); return 0; }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, y, z, num, den, log_y, ln2, k_ln2, result, one, three, kk;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&y, work_prec);
  s21_mpf_init2(&z, work_prec);
  s21_mpf_init2(&num, work_prec);
  s21_mpf_init2(&den, work_prec);
  s21_mpf_init2(&log_y, work_prec);
  s21_mpf_init2(&ln2, work_prec);
  s21_mpf_init2(&k_ln2, work_prec);
  s21_mpf_init2(&result, work_prec);
  s21_mpf_init2(&one, work_prec);
  s21_mpf_init2(&three, work_prec);
  s21_mpf_init2(&kk, work_prec);

  s21_mpf_set(&a, x);
  s21_mpf_set_ui(&one, 1);

  if (s21_mpf_cmp(&a, &one) == 0) {
    s21_mpf_set_zero(res, 0);
    goto cleanup;
  }

  int64_t k = a.exp - 1;
  s21_mpf_set(&y, &a);
  y.exp = 1;
  s21_mpf_normalize(&y);

  s21_mpf_sub(&num, &y, &one);
  s21_mpf_add(&den, &y, &one);
  s21_mpf_div(&z, &num, &den);

  s21_mpf_atanh_series(&log_y, &z, work_prec);
  log_y.exp += 1;
  s21_mpf_normalize(&log_y);

  if (k == 0) {
    s21_mpf_set(res, &log_y);
    goto cleanup;
  }

  s21_mpf_set_ui(&three, 3);
  s21_mpf_div(&z, &one, &three);
  s21_mpf_atanh_series(&ln2, &z, work_prec);
  ln2.exp += 1;
  s21_mpf_normalize(&ln2);

  s21_mpf_set_si(&kk, (long)k);
  s21_mpf_mul(&k_ln2, &kk, &ln2);

  s21_mpf_add(&result, &log_y, &k_ln2);
  s21_mpf_set(res, &result);

cleanup:
  s21_mpf_clear(&a);
  s21_mpf_clear(&y);
  s21_mpf_clear(&z);
  s21_mpf_clear(&num);
  s21_mpf_clear(&den);
  s21_mpf_clear(&log_y);
  s21_mpf_clear(&ln2);
  s21_mpf_clear(&k_ln2);
  s21_mpf_clear(&result);
  s21_mpf_clear(&one);
  s21_mpf_clear(&three);
  s21_mpf_clear(&kk);
  return 0;
}

/* ============================================================
   Утилиты тригонометрии
   ============================================================ */

static void s21_mpf_trunc(s21_mpf_t *res, const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) { s21_mpf_set(res, x); return; }
  if (x->exp <= 0) { s21_mpf_set_zero(res, x->sign); return; }

  s21_mpf_set(res, x);
  if (res->exp >= (int64_t)res->prec) return;

  uint32_t keep = (uint32_t)res->exp;
  for (uint32_t i = 0; i < res->prec - keep; i++) {
    s21_mpf_set_bit(res, i, 0);
  }
  s21_mpf_normalize(res);
}

static void s21_mpf_round_nearest(s21_mpf_t *res, const s21_mpf_t *x) {
  if (x->kind == S21_MPF_ZERO) { s21_mpf_set_zero(res, x->sign); return; }

  s21_mpf_t half, shifted;
  s21_mpf_init2(&half, x->prec);
  s21_mpf_init2(&shifted, x->prec);

  s21_mpf_set_d(&half, 0.5);
  if (x->sign) s21_mpf_neg(&half, &half);

  s21_mpf_add(&shifted, x, &half);
  s21_mpf_trunc(res, &shifted);

  s21_mpf_clear(&half);
  s21_mpf_clear(&shifted);
}

static void s21_mpf_atan_small_series(s21_mpf_t *out, const s21_mpf_t *z,
                                       uint32_t work_prec) {
  s21_mpf_t z2, term, sum, tmp;
  s21_mpf_init2(&z2, work_prec);
  s21_mpf_init2(&term, work_prec);
  s21_mpf_init2(&sum, work_prec);
  s21_mpf_init2(&tmp, work_prec);

  s21_mpf_mul(&z2, z, z);
  s21_mpf_set(&term, z);
  s21_mpf_set(&sum, z);

  for (int n = 3; n < 4000; n += 2) {
    s21_mpf_mul(&term, &term, &z2);
    s21_mpf_neg(&term, &term);

    s21_mpf_set(&tmp, &term);
    s21_mpf_div_small(&tmp, (uint32_t)n);
    s21_mpf_add(&sum, &sum, &tmp);

    if (term.kind == S21_MPF_ZERO) break;
    if (term.exp < sum.exp - (int64_t)work_prec - 4) break;
  }

  s21_mpf_set(out, &sum);
  s21_mpf_clear(&z2);
  s21_mpf_clear(&term);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&tmp);
}

static void s21_mpf_compute_pi(s21_mpf_t *pi, uint32_t work_prec) {
  s21_mpf_t one, five, k239, z, a, b, tmp;
  s21_mpf_init2(&one, work_prec);
  s21_mpf_init2(&five, work_prec);
  s21_mpf_init2(&k239, work_prec);
  s21_mpf_init2(&z, work_prec);
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&b, work_prec);
  s21_mpf_init2(&tmp, work_prec);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_set_ui(&five, 5);
  s21_mpf_set_ui(&k239, 239);

  s21_mpf_div(&z, &one, &five);
  s21_mpf_atan_small_series(&a, &z, work_prec);
  s21_mpf_set_ui(&tmp, 16);
  s21_mpf_mul(&a, &a, &tmp);

  s21_mpf_div(&z, &one, &k239);
  s21_mpf_atan_small_series(&b, &z, work_prec);
  s21_mpf_set_ui(&tmp, 4);
  s21_mpf_mul(&b, &b, &tmp);

  s21_mpf_sub(pi, &a, &b);

  s21_mpf_clear(&one);
  s21_mpf_clear(&five);
  s21_mpf_clear(&k239);
  s21_mpf_clear(&z);
  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&tmp);
}

/* ============================================================
   Константа π
   ============================================================ */

int s21_mpf_pi(s21_mpf_t *res) {
  if (res == NULL) return -1;
  uint32_t work_prec = res->prec + 64;
  s21_mpf_t pi_w;
  s21_mpf_init2(&pi_w, work_prec);
  s21_mpf_compute_pi(&pi_w, work_prec);
  s21_mpf_set(res, &pi_w);
  s21_mpf_clear(&pi_w);
  return 0;
}

/* ============================================================
   Тригонометрия: прямые
   ============================================================ */

int s21_mpf_sin(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign); return 0;
  }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, pi, two_pi, q, q_r, r, term, sum, x2, tmp;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&pi, work_prec);
  s21_mpf_init2(&two_pi, work_prec);
  s21_mpf_init2(&q, work_prec);
  s21_mpf_init2(&q_r, work_prec);
  s21_mpf_init2(&r, work_prec);
  s21_mpf_init2(&term, work_prec);
  s21_mpf_init2(&sum, work_prec);
  s21_mpf_init2(&x2, work_prec);
  s21_mpf_init2(&tmp, work_prec);

  s21_mpf_set(&a, x);

  s21_mpf_compute_pi(&pi, work_prec);
  s21_mpf_set(&two_pi, &pi);
  two_pi.exp += 1;
  s21_mpf_normalize(&two_pi);

  s21_mpf_div(&q, &a, &two_pi);
  s21_mpf_round_nearest(&q_r, &q);
  s21_mpf_mul(&tmp, &q_r, &two_pi);
  s21_mpf_sub(&r, &a, &tmp);

  s21_mpf_set(&term, &r);
  s21_mpf_set(&sum, &r);
  s21_mpf_mul(&x2, &r, &r);

  for (int n = 1; n < 400; n++) {
    s21_mpf_mul(&term, &term, &x2);
    s21_mpf_neg(&term, &term);
    s21_mpf_div_small(&term, (uint32_t)(2 * n));
    s21_mpf_div_small(&term, (uint32_t)(2 * n + 1));
    s21_mpf_add(&sum, &sum, &term);

    if (term.kind == S21_MPF_ZERO) break;
    if (term.exp < sum.exp - (int64_t)work_prec - 4) break;
  }

  s21_mpf_set(res, &sum);

  s21_mpf_clear(&a);
  s21_mpf_clear(&pi);
  s21_mpf_clear(&two_pi);
  s21_mpf_clear(&q);
  s21_mpf_clear(&q_r);
  s21_mpf_clear(&r);
  s21_mpf_clear(&term);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&x2);
  s21_mpf_clear(&tmp);
  return 0;
}

int s21_mpf_cos(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_ui(res, 1); return 0;
  }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, pi, two_pi, q, q_r, r, term, sum, x2, tmp, one;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&pi, work_prec);
  s21_mpf_init2(&two_pi, work_prec);
  s21_mpf_init2(&q, work_prec);
  s21_mpf_init2(&q_r, work_prec);
  s21_mpf_init2(&r, work_prec);
  s21_mpf_init2(&term, work_prec);
  s21_mpf_init2(&sum, work_prec);
  s21_mpf_init2(&x2, work_prec);
  s21_mpf_init2(&tmp, work_prec);
  s21_mpf_init2(&one, work_prec);

  s21_mpf_set(&a, x);

  s21_mpf_compute_pi(&pi, work_prec);
  s21_mpf_set(&two_pi, &pi);
  two_pi.exp += 1;
  s21_mpf_normalize(&two_pi);

  s21_mpf_div(&q, &a, &two_pi);
  s21_mpf_round_nearest(&q_r, &q);
  s21_mpf_mul(&tmp, &q_r, &two_pi);
  s21_mpf_sub(&r, &a, &tmp);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_set(&term, &one);
  s21_mpf_set(&sum, &one);
  s21_mpf_mul(&x2, &r, &r);

  for (int n = 1; n < 400; n++) {
    s21_mpf_mul(&term, &term, &x2);
    s21_mpf_neg(&term, &term);
    s21_mpf_div_small(&term, (uint32_t)(2 * n - 1));
    s21_mpf_div_small(&term, (uint32_t)(2 * n));
    s21_mpf_add(&sum, &sum, &term);

    if (term.kind == S21_MPF_ZERO) break;
    if (term.exp < sum.exp - (int64_t)work_prec - 4) break;
  }

  s21_mpf_set(res, &sum);

  s21_mpf_clear(&a);
  s21_mpf_clear(&pi);
  s21_mpf_clear(&two_pi);
  s21_mpf_clear(&q);
  s21_mpf_clear(&q_r);
  s21_mpf_clear(&r);
  s21_mpf_clear(&term);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&x2);
  s21_mpf_clear(&tmp);
  s21_mpf_clear(&one);
  return 0;
}

int s21_mpf_tan(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res); return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign); return 0;
  }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, s, c, result;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&s, work_prec);
  s21_mpf_init2(&c, work_prec);
  s21_mpf_init2(&result, work_prec);

  s21_mpf_set(&a, x);
  s21_mpf_sin(&s, &a);
  s21_mpf_cos(&c, &a);

  if (s21_mpf_is_zero(&c)) {
    s21_mpf_set_inf(res, 0);
    goto cleanup;
  }

  s21_mpf_div(&result, &s, &c);
  s21_mpf_set(res, &result);

cleanup:
  s21_mpf_clear(&a);
  s21_mpf_clear(&s);
  s21_mpf_clear(&c);
  s21_mpf_clear(&result);
  return 0;
}

/* ============================================================
   Тригонометрия: обратные
   ============================================================ */

int s21_mpf_atan(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) { s21_mpf_set_nan(res); return 0; }
  if (x->kind == S21_MPF_ZERO) { s21_mpf_set_zero(res, x->sign); return 0; }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t one, pi, pi_2, a, tmp, sq;
  s21_mpf_init2(&one, work_prec);
  s21_mpf_init2(&pi, work_prec);
  s21_mpf_init2(&pi_2, work_prec);
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&tmp, work_prec);
  s21_mpf_init2(&sq, work_prec);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_compute_pi(&pi, work_prec);

  if (x->kind == S21_MPF_INF) {
    s21_mpf_set(&pi_2, &pi);
    pi_2.exp -= 1;
    s21_mpf_normalize(&pi_2);
    pi_2.sign = x->sign;
    s21_mpf_set(res, &pi_2);
    goto cleanup;
  }

  int sign = x->sign;
  s21_mpf_abs(&a, x);

  int invert = 0;
  if (s21_mpf_cmp(&a, &one) > 0) {
    s21_mpf_div(&a, &one, &a);
    invert = 1;
  }

  int reductions = 0;
  while (a.exp > -2 && reductions < 16) {
    s21_mpf_mul(&sq, &a, &a);
    s21_mpf_add(&sq, &sq, &one);
    s21_mpf_sqrt(&sq, &sq);
    s21_mpf_add(&sq, &sq, &one);
    s21_mpf_div(&a, &a, &sq);
    reductions++;
  }

  s21_mpf_atan_small_series(&a, &a, work_prec);

  a.exp += reductions;
  s21_mpf_normalize(&a);

  if (invert) {
    s21_mpf_set(&pi_2, &pi);
    pi_2.exp -= 1;
    s21_mpf_normalize(&pi_2);
    s21_mpf_sub(&a, &pi_2, &a);
  }

  if (sign) a.sign = 1;
  s21_mpf_set(res, &a);

cleanup:
  s21_mpf_clear(&one);
  s21_mpf_clear(&pi);
  s21_mpf_clear(&pi_2);
  s21_mpf_clear(&a);
  s21_mpf_clear(&tmp);
  s21_mpf_clear(&sq);
  return 0;
}

int s21_mpf_asin(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) { s21_mpf_set_nan(res); return 0; }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, one, sq, num, den, at, result;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&one, work_prec);
  s21_mpf_init2(&sq, work_prec);
  s21_mpf_init2(&num, work_prec);
  s21_mpf_init2(&den, work_prec);
  s21_mpf_init2(&at, work_prec);
  s21_mpf_init2(&result, work_prec);

  s21_mpf_set(&a, x);

  if (a.kind == S21_MPF_INF) {
    s21_mpf_set_nan(res); goto cleanup;
  }
  s21_mpf_abs(&a, &a);
  s21_mpf_set_ui(&one, 1);
  if (s21_mpf_cmp(&a, &one) > 0) {
    s21_mpf_set_nan(res); goto cleanup;
  }

  s21_mpf_set(&a, x);
  s21_mpf_mul(&sq, &a, &a);
  s21_mpf_sub(&den, &one, &sq);
  s21_mpf_sqrt(&den, &den);
  s21_mpf_div(&num, &a, &den);

  s21_mpf_atan(&at, &num);
  s21_mpf_set(res, &at);

cleanup:
  s21_mpf_clear(&a);
  s21_mpf_clear(&one);
  s21_mpf_clear(&sq);
  s21_mpf_clear(&num);
  s21_mpf_clear(&den);
  s21_mpf_clear(&at);
  s21_mpf_clear(&result);
  return 0;
}

int s21_mpf_acos(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) { s21_mpf_set_nan(res); return 0; }

  uint32_t work_prec = res->prec + 64;

  s21_mpf_t a, one, as, pi_2, result;
  s21_mpf_init2(&a, work_prec);
  s21_mpf_init2(&one, work_prec);
  s21_mpf_init2(&as, work_prec);
  s21_mpf_init2(&pi_2, work_prec);
  s21_mpf_init2(&result, work_prec);

  s21_mpf_set(&a, x);
  s21_mpf_set_ui(&one, 1);
  s21_mpf_abs(&a, &a);

  if (a.kind == S21_MPF_INF || s21_mpf_cmp(&a, &one) > 0) {
    s21_mpf_set_nan(res); goto cleanup;
  }

  if (s21_mpf_cmp(&a, &one) == 0) {
    if (x->sign == 0) {
      s21_mpf_set_zero(res, 0);
    } else {
      s21_mpf_compute_pi(&result, work_prec);
      s21_mpf_set(res, &result);
    }
    goto cleanup;
  }

  s21_mpf_asin(&as, x);

  s21_mpf_compute_pi(&pi_2, work_prec);
  pi_2.exp -= 1;
  s21_mpf_normalize(&pi_2);

  s21_mpf_sub(&result, &pi_2, &as);
  s21_mpf_set(res, &result);

cleanup:
  s21_mpf_clear(&a);
  s21_mpf_clear(&one);
  s21_mpf_clear(&as);
  s21_mpf_clear(&pi_2);
  s21_mpf_clear(&result);
  return 0;
}

/* ============================================================
   Печать
   ============================================================ */

void s21_mpf_print(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_NAN) { printf("NaN\n"); return; }
  if (x->kind == S21_MPF_INF) {
    printf("%sinf\n", x->sign ? "-" : "+"); return;
  }
  if (x->kind == S21_MPF_ZERO) {
    printf("%s0\n", x->sign ? "-" : "+"); return;
  }
  printf("%s0x", x->sign ? "-" : "+");
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    printf("%016" PRIx64, x->limbs[i]);
  }
  printf("p%" PRId64 " [prec=%u]\n", x->exp, x->prec);
}

void s21_mpf_print_d(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_NAN) { printf("nan\n"); return; }
  if (x->kind == S21_MPF_INF) {
    printf("%sinf\n", x->sign ? "-" : "+"); return;
  }
  if (x->kind == S21_MPF_ZERO) {
    printf("%s0\n", x->sign ? "-" : "+"); return;
  }

  size_t count = s21_mpf_limbs_for_prec(x->prec);
  double val = 0.0;
  for (int i = (int)count - 1; i >= 0; i--) {
    val = val * 18446744073709551616.0 + (double)x->limbs[i];
  }
  int64_t e = x->exp - (int64_t)x->prec;
  double result = ldexp(val, (int)e);
  printf("%s%.17g\n", x->sign ? "-" : "+", result);
}
