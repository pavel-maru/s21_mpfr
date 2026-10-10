#ifndef S21_MPF_H
#define S21_MPF_H

#include <stddef.h>
#include <stdint.h>

/* Режимы округления (совместимы с MPFR) */
typedef enum {
  S21_MPF_RNDN = 0, /* к ближайшему, ties to even */
  S21_MPF_RNDZ,     /* к нулю */
  S21_MPF_RNDU,     /* вверх (+inf) */
  S21_MPF_RNDD      /* вниз (-inf) */
} s21_mpf_rnd_t;

/* Тип значения */
typedef enum {
  S21_MPF_ZERO = 0,
  S21_MPF_NORMAL,
  S21_MPF_INF,
  S21_MPF_NAN
} s21_mpf_kind_t;

/* Число произвольной точности.
   Инварианты для NORMAL:
     - limbs — little-endian массив, limbs[0] — младшие 64 бита
     - старший значащий бит находится на позиции (prec - 1)
     - x = (-1)^sign * mant * 2^(exp - prec), где mant = sum(limbs[i] *
   2^(64*i))
     - биты выше prec нулевые */
typedef struct {
  uint64_t *limbs; /* мантисса, little-endian */
  int64_t exp;     /* x = mant * 2^(exp - prec) */
  uint32_t prec;   /* точность в битах */
  int sign;        /* 0 = +, 1 = - */
  int kind;        /* ZERO / NORMAL / INF / NAN */
} s21_mpf_t;

/* Количество 64-битных лимбов для точности prec */
static inline size_t s21_mpf_limbs_for_prec(uint32_t prec) {
  return (prec + 63) / 64;
}

/* ============== Инициализация ============== */

void s21_mpf_init(s21_mpf_t *x); /* prec = 256 по умолчанию */
void s21_mpf_init2(s21_mpf_t *x, uint32_t prec);
void s21_mpf_clear(s21_mpf_t *x);

/* Смена точности на месте, сохраняя значение.
   Расширение — точно.  Сужение — округление по RNDN.
   0 — успех, -1 — ошибка (x == NULL или ошибка аллокации). */
int s21_mpf_set_prec(s21_mpf_t *x, uint32_t prec);

/* ============== Присваивание ============== */

void s21_mpf_set_zero(s21_mpf_t *x, int sign);
void s21_mpf_set_nan(s21_mpf_t *x);
void s21_mpf_set_inf(s21_mpf_t *x, int sign);
void s21_mpf_set_ui(s21_mpf_t *x, unsigned long v);
void s21_mpf_set_si(s21_mpf_t *x, long v);
void s21_mpf_set_d(s21_mpf_t *x, double v);

/* Присваивание по умолчанию (RNDN при смене точности). */
void s21_mpf_set(s21_mpf_t *dst, const s21_mpf_t *src);

/* Присваивание с указанным режимом округления при смене точности.
   Всегда возвращает 0 (или -1 при NULL). */
int s21_mpf_set_round(s21_mpf_t *dst, const s21_mpf_t *src, s21_mpf_rnd_t rnd);

/* ============== Утилиты ============== */

int s21_mpf_is_nan(const s21_mpf_t *x);
int s21_mpf_is_inf(const s21_mpf_t *x);
int s21_mpf_is_zero(const s21_mpf_t *x);
int s21_mpf_is_normal(const s21_mpf_t *x);
int s21_mpf_sign(const s21_mpf_t *x);
uint32_t s21_mpf_get_prec(const s21_mpf_t *x);

/* Низкоуровневые операции над битами мантиссы.
   Сдвиги объявлены во внутреннем заголовке s21_mpf_internal.h. */
int s21_mpf_get_bit(const s21_mpf_t *x, uint32_t pos);
int s21_mpf_msb(const s21_mpf_t *x);
int s21_mpf_lsb(const s21_mpf_t *x);
void s21_mpf_set_bit(s21_mpf_t *x, uint32_t pos, int value);

/* ============== Сравнение ============== */

/* Возвращает -1 / 0 / +1. Для NaN — 0 (как в MPFR). */
int s21_mpf_cmp(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_cmp_abs(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_equal(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_zero_p(const s21_mpf_t *x);
int s21_mpf_integer_p(const s21_mpf_t *x);

/* ============== Арифметика ============== */

int s21_mpf_add(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_sub(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_mul(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_div(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

void s21_mpf_neg(s21_mpf_t *res, const s21_mpf_t *x);
void s21_mpf_abs(s21_mpf_t *res, const s21_mpf_t *x);
void s21_mpf_normalize(s21_mpf_t *x);

/* ============== Элементарные функции ============== */

int s21_mpf_sqrt(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_exp(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_log(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== Константы ============== */

int s21_mpf_pi(s21_mpf_t *res);

/* ============== Тригонометрия: прямые ============== */

int s21_mpf_sin(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_cos(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_tan(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== Тригонометрия: обратные ============== */

int s21_mpf_atan(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_asin(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_acos(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== I/O ============== */

/* Разбор строки в mpf.  Формат: [+-]?(digits[.digits]?|.digits)
   ([eE][+-]?digits)?  Спецзначения: inf, infinity, nan
   (регистронезависимо).  base — только 10.  Возвращает 0 при
   успехе, -1 при ошибке. */
int s21_mpf_set_str(s21_mpf_t *x, const char *str, int base);

/* Конвертирует x в десятичную строку из n_digits значащих цифр
   (с округлением по rnd) и записывает в buf.

   Формат: знак (если отрицательное), затем n_digits цифр без
   десятичной точки.  Позиция точки относительно начала строки
   цифр возвращается в *exp_out: value ≈ 0.DDDD... × 10^exp_out
   (конвенция MPFR).

   Спецзначения: "nan", "inf", "-inf", "0", "-0"; для них
   *exp_out = 0.

   n_digits == 0 → значение по умолчанию (20).

   0 при успехе, -1 при ошибке (NULL, base != 10, буфер мал). */
int s21_mpf_get_str(char *buf, size_t bufsize, int *exp_out,
                    int base, uint32_t n_digits,
                    const s21_mpf_t *x, s21_mpf_rnd_t rnd);

/* ============== Печать ============== */

void s21_mpf_print(const s21_mpf_t *x);
void s21_mpf_print_d(const s21_mpf_t *x);

#endif
