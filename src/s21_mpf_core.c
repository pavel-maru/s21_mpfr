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
     - kind установлен в NORMAL (или ZERO, если всё нулевое)

   Формула: значение x = mant * 2^(exp - prec). После сдвига
   мантиссы влево на k бит значение mant умножается на 2^k, значит
   exp надо уменьшить на k, чтобы x не изменилось. */
static void s21_mpf_normalize(s21_mpf_t *x) {
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  int last = (int)count - 1;
  while (last >= 0 && x->limbs[last] == 0) last--;

  if (last < 0) {
    x->kind = S21_MPF_ZERO;
    x->exp = 0;
    return;
  }

  /* Позиция старшего значащего бита во всём массиве */
  int top_bit_global = last * 64 + (63 - __builtin_clzll(x->limbs[last]));

  /* Сколько бит нужно сдвинуть влево, чтобы старший бит оказался
     на позиции (prec - 1) */
  int shift = (int)(x->prec - 1) - top_bit_global;

  if (shift > 0) {
    /* Сдвиг влево на shift бит */
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
    /* Сдвиг вправо на -shift бит. Для set_* это не происходит
       (там все биты влезают), но функция универсальна. */
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

  /* Замаскировать биты выше prec, если prec не кратен 64 */
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
  /* До нормализации считаем, что мантисса занимает prec бит
     с ведущими нулями. Значение x = v = mant * 2^(exp - prec),
     где mant = v, значит exp = prec. */
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
    /* Денормализованное. Для MVP грубо обрабатываем. */
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

  /* Помещаем mant_d (53 бита) в мантиссу mpf как есть.
     Значение x = mant_d * 2^(exp_d - 53). В формате mpf:
     x = mant_mpf * 2^(exp_mpf - prec), где mant_mpf = mant_d.
     Значит exp_mpf = exp_d - 53 + prec. */
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

  /* Собираем значение через длинное деление:
     mant = сумма limbs[i] * 2^(64*i). Мы знаем exp. Хотим
     получить значение в double. mant может быть очень большим,
     поэтому собираем пошагово с масштабированием. */
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  double val = 0.0;
  for (int i = (int)count - 1; i >= 0; i--) {
    val = val * 18446744073709551616.0 + (double)x->limbs[i];
  }
  /* Теперь val = mant (в double). Нужно умножить на 2^(exp - prec). */
  int64_t e = x->exp - (int64_t)x->prec;
  /* Масштабируем через ldexp, который работает до 2^1023 */
  double result = ldexp(val, (int)e);
  printf("%s%.17g\n", x->sign ? "-" : "+", result);
}
