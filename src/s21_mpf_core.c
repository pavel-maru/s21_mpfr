#include "s21_mpf.h"

#include <inttypes.h>
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

/* Установить точность. TODO: правильно округлять при уменьшении. */
void s21_mpf_set_prec(s21_mpf_t *x, uint32_t prec) {
  if (prec < 2) prec = 2;
  if (prec == x->prec) return;

  uint64_t *new_limbs = calloc(s21_mpf_limbs_for_prec(prec), sizeof(uint64_t));

  /* Пока без округления: копируем что влезает */
  size_t old_count = s21_mpf_limbs_for_prec(x->prec);
  size_t new_count = s21_mpf_limbs_for_prec(prec);
  size_t n = old_count < new_count ? old_count : new_count;
  memcpy(new_limbs, x->limbs, n * sizeof(uint64_t));

  free(x->limbs);
  x->limbs = new_limbs;
  x->prec = prec;
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

/* Нормализовать: mant должен занимать ровно prec бит, старший бит установлен.
   mant передаётся как массив лимбов, содержащий значащие биты в младших
   позициях. Функция сдвигает влево так, чтобы старший бит оказался на
   позиции prec-1, и соответствующим образом корректирует exp.
   Возвращает: 0 = OK, иначе — число потерянных бит. */
static int s21_mpf_normalize(s21_mpf_t *x) {
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  size_t last = count - 1;

  /* Считаем ведущие нули в старшем лимбе */
  if (x->limbs[last] == 0) {
    /* Всё число — ноль */
    x->kind = S21_MPF_ZERO;
    return 0;
  }

  int clz = __builtin_clzll(x->limbs[last]);
  /* Хотим, чтобы старший бит был на позиции prec-1.
     prec может быть не кратен 64. */
  int target_shift = clz - (int)(count * 64 - x->prec);
  /* target_shift > 0 — сдвигаем влево, < 0 — вправо */

  if (target_shift > 0) {
    /* Сдвиг влево на target_shift бит */
    int word_shift = target_shift / 64;
    int bit_shift = target_shift % 64;
    for (size_t i = last; i < count; i--) {  /* безопасно, size_t обёртка */
      /* обходим — используем индексы int64 */
    }
    /* Реализуем сдвиг через простой цикл (для MVP) */
    for (int i = (int)count - 1; i >= word_shift; i--) {
      int src = i - word_shift;
      uint64_t v = x->limbs[src];
      if (bit_shift && src > 0) {
        v = (v << bit_shift) | (x->limbs[src - 1] >> (64 - bit_shift));
      } else if (bit_shift) {
        v <<= bit_shift;
      }
      x->limbs[i] = v;
    }
    for (int i = 0; i < word_shift; i++) x->limbs[i] = 0;
    x->exp -= target_shift;
  } else if (target_shift < 0) {
    /* Сдвиг вправо — потеря битов, для нормализации set_* это не нужно,
       но оставим на будущее */
    /* (пока пропускаем) */
  }

  x->kind = S21_MPF_NORMAL;
  return 0;
}

void s21_mpf_set_ui(s21_mpf_t *x, unsigned long v) {
  if (v == 0) {
    s21_mpf_set_zero(x, 0);
    return;
  }
  memset(x->limbs, 0, s21_mpf_limbs_for_prec(x->prec) * sizeof(uint64_t));
  x->limbs[0] = v;
  x->sign = 0;
  x->exp = 64;  /* x = v * 2^(64 - prec)... скорректируется в normalize */
  s21_mpf_normalize(x);
}

void s21_mpf_set_si(s21_mpf_t *x, long v) {
  if (v < 0) {
    s21_mpf_set_ui(x, (unsigned long)(-v));
    x->sign = 1;
  } else {
    s21_mpf_set_ui(x, (unsigned long)v);
  }
}

/* Разбор double на компоненты.
   Возвращает мантиссу как uint64_t (со скрытой единицей), 
   экспоненту в формате "x = mant * 2^(exp - 53)", знак. */
static void s21_decompose_double(double v, uint64_t *mant, int64_t *exp,
                                  int *sign) {
  union { double d; uint64_t u; } u;
  u.d = v;

  *sign = (int)(u.u >> 63);
  int64_t e_raw = (int64_t)((u.u >> 52) & 0x7FF);
  uint64_t frac = u.u & 0xFFFFFFFFFFFFFULL;  /* младшие 52 бита */

  if (e_raw == 0) {
    /* Денормализованное или ноль (TODO: обработать денормализованные) */
    *mant = frac;
    *exp = 1 - 1023 - 52;  /* упрощённо */
  } else {
    *mant = (1ULL << 52) | frac;
    *exp = e_raw - 1023 - 52 + 53;  /* = e_raw - 1022 */
  }
}

void s21_mpf_set_d(s21_mpf_t *x, double v) {
  /* Обработка особых случаев */
  if (v != v) { s21_mpf_set_nan(x); return; }
  if (v == 1.0 / 0.0) { s21_mpf_set_inf(x, 0); return; }
  if (v == -1.0 / 0.0) { s21_mpf_set_inf(x, 1); return; }
  if (v == 0.0) {
    s21_mpf_set_zero(x, (1.0 / v < 0) ? 1 : 0);
    return;
  }

  uint64_t mant_d;
  int64_t exp_d;
  int sign_d;
  s21_decompose_double(v, &mant_d, &exp_d, &sign_d);

  /* mant_d занимает 53 бита (старший установлен).
     Переводим в mpf с prec битами: сдвигаем влево на (prec - 53).
     exp корректируем так: x = mant_d * 2^(exp_d - 53) = mant_mpf * 2^(exp_mpf - prec)
     mant_mpf = mant_d << (prec - 53)
     => exp_mpf = exp_d - 53 + (prec - 53) ... не совсем.
     Проще: mant_mpf * 2^(exp_mpf - prec) = mant_d * 2^(exp_d - 53)
     mant_d << (prec - 53) * 2^(exp_mpf - prec) = mant_d * 2^(exp_d - 53)
     2^(prec - 53) * 2^(exp_mpf - prec) = 2^(exp_d - 53)
     exp_mpf = exp_d
  */
  memset(x->limbs, 0, s21_mpf_limbs_for_prec(x->prec) * sizeof(uint64_t));

  if (x->prec >= 53) {
    x->limbs[0] = mant_d << (x->prec >= 64 ? 0 : (53 - x->prec));
    /* Более корректно: помещаем mant_d в старшие 53 бита общего prec-битного поля */
    /* Пока упрощённо, сдвиг влево на (prec - 53) при prec >= 53 */
    uint32_t shift = x->prec - 53;
    int word = shift / 64;
    int bit = shift % 64;
    memset(x->limbs, 0, s21_mpf_limbs_for_prec(x->prec) * sizeof(uint64_t));
    x->limbs[word] = mant_d << bit;
    if (bit && word + 1 < (int)s21_mpf_limbs_for_prec(x->prec)) {
      x->limbs[word + 1] = mant_d >> (64 - bit);
    }
    x->exp = exp_d;
  } else {
    /* Уменьшаем точность — теряем биты (TODO: округление) */
    uint32_t drop = 53 - x->prec;
    x->limbs[0] = mant_d >> drop;
    x->exp = exp_d;
  }

  x->sign = sign_d;
  x->kind = S21_MPF_NORMAL;
}

void s21_mpf_set(s21_mpf_t *dst, const s21_mpf_t *src) {
  if (dst == src) return;
  if (dst->prec != src->prec) {
    /* Разные точности: пока просто пересоздаём с той же точностью, что у dst,
       копируем что влезает (TODO: правильное округление) */
    size_t dst_count = s21_mpf_limbs_for_prec(dst->prec);
    size_t src_count = s21_mpf_limbs_for_prec(src->prec);
    memset(dst->limbs, 0, dst_count * sizeof(uint64_t));
    size_t n = dst_count < src_count ? dst_count : src_count;
    memcpy(dst->limbs, src->limbs, n * sizeof(uint64_t));
  } else {
    memcpy(dst->limbs, src->limbs,
           s21_mpf_limbs_for_prec(dst->prec) * sizeof(uint64_t));
  }
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

int s21_mpf_sign(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_ZERO || x->kind == S21_MPF_NAN) return 0;
  return x->sign ? -1 : 1;
}

uint32_t s21_mpf_get_prec(const s21_mpf_t *x) { return x->prec; }

void s21_mpf_print(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_NAN) { printf("NaN\n"); return; }
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
    printf("%016lx", x->limbs[i]);
  }
  printf("p%" PRId64 " [prec=%u]\n", x->exp, x->prec);
}

void s21_mpf_print_d(const s21_mpf_t *x) {
  /* Временная реализация: собираем double из старших лимбов */
  if (x->kind == S21_MPF_NAN) { printf("nan\n"); return; }
  if (x->kind == S21_MPF_INF) {
    printf("%sinf\n", x->sign ? "-" : "+");
    return;
  }
  if (x->kind == S21_MPF_ZERO) {
    printf("%s0\n", x->sign ? "-" : "+");
    return;
  }
  /* Упрощённо: берём старший лимб и собираем приблизительное double.
     Для отладки достаточно. */
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  uint64_t hi = x->limbs[count - 1];
  int clz = __builtin_clzll(hi);
  int64_t real_exp = x->exp - (int64_t)x->prec + (int64_t)(count * 64) - clz;
  double val = (double)hi;
  for (int i = (int)count - 2; i >= 0; i--) {
    val = val * 18446744073709551616.0 + (double)x->limbs[i];
    if (val > 1e300) break;
  }
  double scale = 1.0;
  int64_t e = real_exp;
  while (e > 0) { scale *= 2.0; e--; }
  while (e < 0) { scale /= 2.0; e++; }
  printf("%s%g\n", x->sign ? "-" : "+", val * scale);
}
