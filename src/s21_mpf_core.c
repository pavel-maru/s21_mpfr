#include "s21_mpf.h"

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
  /* NB: семантика set_prec будет отдельно ревизоваться (не сохраняет
     инвариант MSB@prec-1). Сейчас поведение оставлено как было. */
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

void s21_mpf_set(s21_mpf_t *dst, const s21_mpf_t *src) {
  s21_mpf_set_round(dst, src, S21_MPF_RNDN);
}

int s21_mpf_set_round(s21_mpf_t *dst, const s21_mpf_t *src, s21_mpf_rnd_t rnd) {
  if (dst == NULL || src == NULL) return -1;
  if (dst == src) return 0;

  size_t dst_count = s21_mpf_limbs_for_prec(dst->prec);
  memset(dst->limbs, 0, dst_count * sizeof(uint64_t));

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

  if (dst->prec > src->prec) {
    int shift = (int)(dst->prec - src->prec);
    s21_mpf_shift_left_into(dst->limbs, src->limbs, dst_count, shift);
    dst->exp = src->exp;
    dst->sign = src->sign;
    dst->kind = S21_MPF_NORMAL;
    return 0;
  }

  int shift = (int)(src->prec - dst->prec);
  int word_shift = shift / 64;
  int bit_shift = shift % 64;

  int round_bit = 0;
  {
    int pos = shift - 1;
    round_bit = (int)((src->limbs[pos / 64] >> (pos % 64)) & 1ULL);
  }

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

  uint32_t mask_bits = dst->prec % 64;
  if (mask_bits != 0 && dst_count > 0) {
    dst->limbs[dst_count - 1] &= (1ULL << mask_bits) - 1;
  }

  dst->exp = src->exp;
  dst->sign = src->sign;
  dst->kind = S21_MPF_NORMAL;

  int round_up = 0;
  if (round_bit || sticky) {
    switch (rnd) {
      case S21_MPF_RNDZ: round_up = 0; break;
      case S21_MPF_RNDU: round_up = (src->sign == 0); break;
      case S21_MPF_RNDD: round_up = (src->sign == 1); break;
      case S21_MPF_RNDN:
      default:
        if (round_bit == 0) round_up = 0;
        else if (sticky) round_up = 1;
        else round_up = (int)(dst->limbs[0] & 1ULL);
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
