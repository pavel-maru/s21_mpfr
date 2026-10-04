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

/* Смена точности. TODO: правильное округление при уменьшении prec. */
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

/* ============================================================
   Нормализация мантиссы
   ============================================================

   Приводит число к каноническому виду:
     - старший бит мантиссы стоит на позиции (prec - 1)
     - все биты выше prec нулевые
     - exp скорректирован так, что значение не изменилось
     - kind установлен в NORMAL (или ZERO, если всё нулевое) */
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

/* Разбор double на компоненты:
   v = (-1)^sign * mant * 2^(exp - 53), где mant — 53-битное
   целое со старшим установленным битом. */
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
  if (v != v) {
    s21_mpf_set_nan(x);
    return;
  }
  if (v == 1.0 / 0.0) {
    s21_mpf_set_inf(x, 0);
    return;
  }
  if (v == -1.0 / 0.0) {
    s21_mpf_set_inf(x, 1);
    return;
  }
  if (v == 0.0) {
    s21_mpf_set_zero(x, (1.0 / v < 0) ? 1 : 0);
    return;
  }

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

void s21_mpf_set(s21_mpf_t *dst, const s21_mpf_t *src) {
  if (dst == src) return;
  size_t dst_count = s21_mpf_limbs_for_prec(dst->prec);
  size_t src_count = s21_mpf_limbs_for_prec(src->prec);
  memset(dst->limbs, 0, dst_count * sizeof(uint64_t));
  size_t n = dst_count < src_count ? dst_count : src_count;
  memcpy(dst->limbs, src->limbs, n * sizeof(uint64_t));
  dst->exp = src->exp;
  dst->sign = src->sign;
  dst->kind = src->kind;
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

/* Сравнение мантисс побитово, начиная со старшего бита.
   Предполагается, что x->prec == y->prec и обе мантиссы
   нормализованы (старший значащий бит на позиции prec-1). */
static int s21_mpf_cmp_mant(const s21_mpf_t *x, const s21_mpf_t *y) {
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    if (x->limbs[i] < y->limbs[i]) return -1;
    if (x->limbs[i] > y->limbs[i]) return 1;
  }
  return 0;
}

/* Сравнение по абсолютной величине. Требует, чтобы x и y имели
   одинаковую точность. Для разных точностей — временно через double. */
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

  if (x->prec == y->prec) {
    return s21_mpf_cmp_mant(x, y);
  }

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

/* Сдвиг мантиссы вправо на shift бит: dst = src >> shift.
   Старшие биты заполняются нулями. */
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

/* |x| + |y|. Оба NORMAL. Результат в res (sign=0). */
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

/* |x| - |y|. Требует |x| >= |y|, оба NORMAL. Результат в res (sign=0). */
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
  s21_mpf_set(res, x);
  if (res->kind != S21_MPF_NAN) {
    res->sign = !res->sign;
  }
}

void s21_mpf_abs(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_set(res, x);
  if (res->kind == S21_MPF_NORMAL || res->kind == S21_MPF_INF) {
    res->sign = 0;
  }
}

int s21_mpf_add(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res->prec != x->prec || x->prec != y->prec) return -1;

  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_INF && y->kind == S21_MPF_INF) {
    if (x->sign == y->sign) {
      s21_mpf_set_inf(res, x->sign);
    } else {
      s21_mpf_set_nan(res);
    }
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, x->sign);
    return 0;
  }
  if (y->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, y->sign);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set(res, y);
    return 0;
  }
  if (y->kind == S21_MPF_ZERO) {
    s21_mpf_set(res, x);
    return 0;
  }

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
  if (res->prec != x->prec || x->prec != y->prec) return -1;

  s21_mpf_t neg_y;
  s21_mpf_init2(&neg_y, y->prec);
  s21_mpf_neg(&neg_y, y);
  int rc = s21_mpf_add(res, x, &neg_y);
  s21_mpf_clear(&neg_y);
  return rc;
}

/* ============================================================
   Печать
   ============================================================ */

void s21_mpf_print(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_NAN) {
    printf("NaN\n");
    return;
  }
  if (x->kind == S21_MPF_INF) {
    printf("%sinf\n", x->sign ? "-" : "+");
    return;
  }
  if (x->kind == S21_MPF_ZERO) {
    printf("%s0\n", x->sign ? "-" : "+");
    return;
  }
  printf("%s0x", x->sign ? "-" : "+");
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    printf("%016" PRIx64, x->limbs[i]);
  }
  printf("p%" PRId64 " [prec=%u]\n", x->exp, x->prec);
}

void s21_mpf_print_d(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_NAN) {
    printf("nan\n");
    return;
  }
  if (x->kind == S21_MPF_INF) {
    printf("%sinf\n", x->sign ? "-" : "+");
    return;
  }
  if (x->kind == S21_MPF_ZERO) {
    printf("%s0\n", x->sign ? "-" : "+");
    return;
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
