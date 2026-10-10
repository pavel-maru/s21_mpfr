#include <stdlib.h>
#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

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
    s21_mpf_shift_left_into(x->limbs, count, x->limbs, count, shift);
    x->exp -= shift;
  } else if (shift < 0) {
    s21_mpf_shift_right_into(x->limbs, count, x->limbs, count, -shift);
    x->exp += -shift;
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
    /* Субнормальные: value = frac * 2^(-1074).
       Вес LSB мантиссы — 2^(-1074), значит exp_d - 53 = -1074,
       exp_d = -1021. */
    *mant = frac;
    *exp = 1 - 1022;
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
  s21_mpf_set_round(dst, src, S21_MPF_RNDN);
}

/* ============================================================
   Подпрограммы для управления точностью
   ============================================================ */

/* Есть ли единичные биты ниже позиции (shift - 1)? */
static int mpf_sticky_below(const s21_mpf_t *src, int shift) {
  int top = shift - 1;
  int full_word = top / 64;
  int rem_bits = top % 64;
  for (int i = 0; i < full_word; i++) {
    if (src->limbs[i] != 0) return 1;
  }
  if (rem_bits > 0) {
    uint64_t mask = (1ULL << rem_bits) - 1;
    if ((src->limbs[full_word] & mask) != 0) return 1;
  }
  return 0;
}

/* Решение об округлении вверх. */
static int mpf_round_up(s21_mpf_rnd_t rnd, int round_bit, int sticky, int sign,
                        uint64_t lsb) {
  if (!round_bit && !sticky) return 0;
  switch (rnd) {
    case S21_MPF_RNDZ:
      return 0;
    case S21_MPF_RNDU:
      return sign == 0;
    case S21_MPF_RNDD:
      return sign == 1;
    case S21_MPF_RNDN:
    default:
      if (!round_bit) return 0;
      if (sticky) return 1;
      return (int)(lsb & 1ULL);
  }
}

/* Прибавить 1 к мантиссе. При переполнении prec-битной сетки —
   нормализовать к 2^(prec-1) и увеличить exp. Два случая:
     (а) +1 вынеслось за все лимбы — обрабатывается carry-веткой;
     (б) +1 перешло в бит prec внутри последнего лимба (prec % 64 != 0) —
         обрабатывается вторым блоком. */
static void mpf_add_one(uint64_t *limbs, size_t count, uint32_t prec,
                        int64_t *exp) {
  uint64_t carry = 1;
  for (size_t i = 0; i < count && carry; i++) {
    limbs[i]++;
    if (limbs[i] != 0) carry = 0;
  }

  if (carry) {
    /* Случай (а). */
    memset(limbs, 0, count * sizeof(uint64_t));
    uint32_t top = prec - 1;
    limbs[top / 64] = 1ULL << (top % 64);
    *exp += 1;
    return;
  }

  /* Случай (б). Проверяем бит на позиции prec. */
  uint32_t pw = prec / 64;
  uint32_t pb = prec % 64;
  if (pw < count && ((limbs[pw] >> pb) & 1ULL)) {
    for (size_t i = 0; i < count; i++) {
      uint64_t high = (i + 1 < count) ? limbs[i + 1] : 0;
      limbs[i] = (limbs[i] >> 1) | (high << 63);
    }
    *exp += 1;
  }
}

int s21_mpf_set_round(s21_mpf_t *dst, const s21_mpf_t *src, s21_mpf_rnd_t rnd) {
  if (dst == NULL || src == NULL) return -1;
  if (dst == src) return 0;

  size_t dst_count = s21_mpf_limbs_for_prec(dst->prec);
  memset(dst->limbs, 0, dst_count * sizeof(uint64_t));

  /* Спецзначения и равные точности — копирование как есть. */
  if (src->kind != S21_MPF_NORMAL || dst->prec == src->prec) {
    size_t src_count = s21_mpf_limbs_for_prec(src->prec);
    size_t n = dst_count < src_count ? dst_count : src_count;
    memcpy(dst->limbs, src->limbs, n * sizeof(uint64_t));
    dst->exp = src->exp;
    dst->sign = src->sign;
    dst->kind = src->kind;
    return 0;
  }

  dst->exp = src->exp;
  dst->sign = src->sign;
  dst->kind = S21_MPF_NORMAL;

  size_t src_count = s21_mpf_limbs_for_prec(src->prec);

  /* Расширение — всегда точно. */
  if (dst->prec > src->prec) {
    s21_mpf_shift_left_into(dst->limbs, dst_count, src->limbs, src_count,
                            (int)(dst->prec - src->prec));
    return 0;
  }

  /* Сужение — округляем. */
  int shift = (int)(src->prec - dst->prec);
  int round_bit =
      (int)((src->limbs[(shift - 1) / 64] >> ((shift - 1) % 64)) & 1ULL);
  int sticky = mpf_sticky_below(src, shift);

  s21_mpf_shift_right_into(dst->limbs, dst_count, src->limbs, src_count, shift);

  uint32_t mask_bits = dst->prec % 64;
  if (mask_bits != 0 && dst_count > 0)
    dst->limbs[dst_count - 1] &= (1ULL << mask_bits) - 1;

  if (mpf_round_up(rnd, round_bit, sticky, src->sign, dst->limbs[0]))
    mpf_add_one(dst->limbs, dst_count, dst->prec, &dst->exp);

  return 0;
}

/* Смена точности на месте, сохраняя значение.
   Расширение — точно.  Сужение — округление по RNDN.
   Спецзначения (ZERO / INF / NAN) переносятся как есть.

   Опорный инвариант: exp не зависит от prec и равен «показателю
   MSB + 1», то есть вес старшего бита мантиссы — это 2^(exp - 1).
   При смене prec мантисса сдвигается так, чтобы этот вес не
   изменился, а exp остаётся на месте (кроме случая округления
   вверх через границу, где exp инкрементируется). */
int s21_mpf_set_prec(s21_mpf_t *x, uint32_t new_prec) {
  if (x == NULL) return -1;
  if (new_prec < 2) new_prec = 2;
  if (new_prec == x->prec) return 0;

  size_t new_count = s21_mpf_limbs_for_prec(new_prec);
  uint64_t *new_limbs = calloc(new_count, sizeof(uint64_t));
  if (new_limbs == NULL) return -1;

  /* Спецзначения: буфер перераспределяем, поля не трогаем. */
  if (x->kind != S21_MPF_NORMAL) {
    free(x->limbs);
    x->limbs = new_limbs;
    x->prec = new_prec;
    return 0;
  }

  size_t old_count = s21_mpf_limbs_for_prec(x->prec);
  int64_t new_exp = x->exp;

  if (new_prec > x->prec) {
    /* Расширение — точно. */
    s21_mpf_shift_left_into(new_limbs, new_count, x->limbs, old_count,
                            (int)(new_prec - x->prec));
  } else {
    /* Сужение — округление RNDN. */
    int delta = (int)(x->prec - new_prec);
    int round_bit =
        (int)((x->limbs[(delta - 1) / 64] >> ((delta - 1) % 64)) & 1ULL);
    int sticky = mpf_sticky_below(x, delta);

    s21_mpf_shift_right_into(new_limbs, new_count, x->limbs, old_count, delta);

    uint32_t mask_bits = new_prec % 64;
    if (mask_bits != 0) new_limbs[new_count - 1] &= (1ULL << mask_bits) - 1;

    if (mpf_round_up(S21_MPF_RNDN, round_bit, sticky, x->sign, new_limbs[0]))
      mpf_add_one(new_limbs, new_count, new_prec, &new_exp);
  }

  /* Коммит: только здесь объект меняется. */
  free(x->limbs);
  x->limbs = new_limbs;
  x->prec = new_prec;
  x->exp = new_exp;
  return 0;
}

/* ============================================================
   Временные mpf на стеке
   ============================================================ */

void s21_mpf_stack_init(s21_mpf_stack_t *st, uint32_t prec) {
  if (prec < 2) prec = 2;
  size_t count = s21_mpf_limbs_for_prec(prec);

  st->mpf.prec = prec;
  st->mpf.exp = 0;
  st->mpf.sign = 0;
  st->mpf.kind = S21_MPF_ZERO;

  if (count <= S21_MPF_STACK_LIMBS) {
    /* Зануляем только используемую часть буфера (count младших
       лимбов). Остальные стековые лимбы не читаются — это
       эквивалент calloc'а, но без лишней работы. */
    memset(st->stack_buf, 0, count * sizeof(uint64_t));
    st->mpf.limbs = st->stack_buf;
    st->heap_limbs = NULL;
  } else {
    st->mpf.limbs = calloc(count, sizeof(uint64_t));
    st->heap_limbs = st->mpf.limbs;
  }
}

void s21_mpf_stack_clear(s21_mpf_stack_t *st) {
  if (st->heap_limbs != NULL) free(st->heap_limbs);
  st->heap_limbs = NULL;
  st->mpf.limbs = NULL;
  st->mpf.prec = 0;
  st->mpf.exp = 0;
  st->mpf.sign = 0;
  st->mpf.kind = S21_MPF_ZERO;
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
