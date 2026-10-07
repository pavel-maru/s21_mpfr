#ifndef S21_MPF_INTERNAL_H
#define S21_MPF_INTERNAL_H

#include "s21_mpf.h"

/* ------------------------------------------------------------
   Единые сдвиги мантиссы.

   dst_count и src_count — размеры массивов в лимбах, могут
   различаться. dst и src могут указывать на один буфер при
   dst_count == src_count. Любой индекс за пределами src
   трактуется как ноль.
   ------------------------------------------------------------ */

void s21_mpf_shift_left_into(uint64_t *dst, size_t dst_count,
                             const uint64_t *src, size_t src_count,
                             int shift);
void s21_mpf_shift_right_into(uint64_t *dst, size_t dst_count,
                              const uint64_t *src, size_t src_count,
                              int shift);

/* ------------------------------------------------------------
   Сырая арифметика.

   Все эти функции работают в предположении, что операнды уже
   приведены к одной точности (= prec результата) и что специальные
   значения (NaN / Inf / Zero) обработаны на уровне публичных
   обёрток. Никаких проверок на NULL и prec.
   ------------------------------------------------------------ */

void s21_mpf_add_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
void s21_mpf_sub_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
void s21_mpf_neg_raw(s21_mpf_t *res, const s21_mpf_t *x);
int  s21_mpf_mul_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int  s21_mpf_div_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

void s21_mpf_div_small(s21_mpf_t *x, uint32_t n);

/* ------------------------------------------------------------
   Проверки аргументов. 0 — успех, -1 — ошибка.
   ------------------------------------------------------------ */

static inline int s21_mpf_check_binary(const s21_mpf_t *res,
                                       const s21_mpf_t *x,
                                       const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;
  if (res->prec != x->prec || x->prec != y->prec) return -1;
  return 0;
}

static inline int s21_mpf_check_unary(const s21_mpf_t *res,
                                      const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  return 0;
}

/* ------------------------------------------------------------
   Внутренние _impl-варианты тригонометрии.

   pi_opt — уже вычисленная π (может быть NULL). Если NULL и π
   действительно нужна операции — она вычисляется лениво внутри.
   Спецзначения (NaN / Inf / Zero) обрабатываются в публичных
   обёртках, до вызова _impl.
   ------------------------------------------------------------ */

void s21_mpf_sin_impl(s21_mpf_t *res, const s21_mpf_t *x,
                      const s21_mpf_t *pi_opt);
void s21_mpf_cos_impl(s21_mpf_t *res, const s21_mpf_t *x,
                      const s21_mpf_t *pi_opt);
void s21_mpf_atan_impl(s21_mpf_t *res, const s21_mpf_t *x,
                       const s21_mpf_t *pi_opt);
void s21_mpf_asin_impl(s21_mpf_t *res, const s21_mpf_t *x,
                       const s21_mpf_t *pi_opt);

#endif
