#ifndef S21_MPF_H
#define S21_MPF_H

#include <stddef.h>
#include <stdint.h>

/* Режимы округления (совместимы с MPFR) */
typedef enum {
  S21_MPF_RNDN = 0,  /* к ближайшему, ties to even */
  S21_MPF_RNDZ,      /* к нулю */
  S21_MPF_RNDU,      /* вверх (+inf) */
  S21_MPF_RNDD       /* вниз (-inf) */
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
     - x = (-1)^sign * mant * 2^(exp - prec), где mant = sum(limbs[i] * 2^(64*i))
     - биты выше prec нулевые */
typedef struct {
  uint64_t *limbs;   /* мантисса, little-endian */
  int64_t   exp;     /* x = mant * 2^(exp - prec) */
  uint32_t  prec;    /* точность в битах */
  int       sign;    /* 0 = +, 1 = - */
  int       kind;    /* ZERO / NORMAL / INF / NAN */
} s21_mpf_t;

/* Количество 64-битных лимбов для точности prec */
static inline size_t s21_mpf_limbs_for_prec(uint32_t prec) {
  return (prec + 63) / 64;
}

/* ============== Инициализация ============== */

void s21_mpf_init(s21_mpf_t *x);                /* prec = 256 по умолчанию */
void s21_mpf_init2(s21_mpf_t *x, uint32_t prec);
void s21_mpf_clear(s21_mpf_t *x);
void s21_mpf_set_prec(s21_mpf_t *x, uint32_t prec);

/* ============== Присваивание ============== */

void s21_mpf_set_zero(s21_mpf_t *x, int sign);
void s21_mpf_set_nan(s21_mpf_t *x);
void s21_mpf_set_inf(s21_mpf_t *x, int sign);
void s21_mpf_set_ui(s21_mpf_t *x, unsigned long v);
void s21_mpf_set_si(s21_mpf_t *x, long v);
void s21_mpf_set_d(s21_mpf_t *x, double v);
void s21_mpf_set(s21_mpf_t *dst, const s21_mpf_t *src);

/* ============== Утилиты ============== */

int      s21_mpf_is_nan(const s21_mpf_t *x);
int      s21_mpf_is_inf(const s21_mpf_t *x);
int      s21_mpf_is_zero(const s21_mpf_t *x);
int      s21_mpf_is_normal(const s21_mpf_t *x);
int      s21_mpf_sign(const s21_mpf_t *x);
uint32_t s21_mpf_get_prec(const s21_mpf_t *x);

/* Низкоуровневые операции над битами мантиссы */
int  s21_mpf_get_bit(const s21_mpf_t *x, uint32_t pos);
int  s21_mpf_msb(const s21_mpf_t *x);
int  s21_mpf_lsb(const s21_mpf_t *x);
void s21_mpf_set_bit(s21_mpf_t *x, uint32_t pos, int value);
void s21_mpf_shift_left_into(uint64_t *dst, const uint64_t *src,
                             size_t count, int shift);

/* ============== Сравнение ============== */

/* Возвращает -1 / 0 / +1. Для NaN — 0 (как в MPFR). */
int s21_mpf_cmp(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_cmp_abs(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_equal(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_zero_p(const s21_mpf_t *x);
int s21_mpf_integer_p(const s21_mpf_t *x);

/* ============== Арифметика ============== */

/* res = x + y. Все три должны иметь одинаковую точность. */
int s21_mpf_add(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_sub(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_mul(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_div(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

/* res = -x, res = |x| */
void s21_mpf_neg(s21_mpf_t *res, const s21_mpf_t *x);
void s21_mpf_abs(s21_mpf_t *res, const s21_mpf_t *x);

/* Нормализация (обычно вызывается автоматически) */
void s21_mpf_normalize(s21_mpf_t *x);

/* ============== Элементарные функции ============== */

int s21_mpf_sqrt(s21_mpf_t *res, const s21_mpf_t *x);   /* x >= 0 */
int s21_mpf_exp(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_log(s21_mpf_t *res, const s21_mpf_t *x);    /* x > 0 */

/* ============== Константы ============== */

/* res = π с точностью res->prec */
int s21_mpf_pi(s21_mpf_t *res);

/* ============== Тригонометрия: прямые ============== */

int s21_mpf_sin(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_cos(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_tan(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== Тригонометрия: обратные ============== */

int s21_mpf_atan(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_asin(s21_mpf_t *res, const s21_mpf_t *x);   /* |x| <= 1 */
int s21_mpf_acos(s21_mpf_t *res, const s21_mpf_t *x);   /* |x| <= 1 */

/* ============== Печать ============== */

void s21_mpf_print(const s21_mpf_t *x);
void s21_mpf_print_d(const s21_mpf_t *x);

#endif
