#ifndef S21_MPF_INTERNAL_H
#define S21_MPF_INTERNAL_H

#include "s21_mpf.h"

/* ------------------------------------------------------------
   Сырая арифметика.

   Все эти функции работают в предположении, что операнды уже
   приведены к одной точности (= prec результата) и что специальные
   значения (NaN / Inf / Zero) обработаны на уровне публичных
   обёрток в s21_mpf_arith.c. Никаких проверок на NULL и prec.
   ------------------------------------------------------------ */

void s21_mpf_add_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
void s21_mpf_sub_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
void s21_mpf_neg_raw(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_mul_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_div_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

void s21_mpf_div_small(s21_mpf_t *x, uint32_t n);

#endif
