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
int s21_mpf_mul_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_div_raw(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

void s21_mpf_div_small(s21_mpf_t *x, uint32_t n);

/* ------------------------------------------------------------
   Временные mpf-объекты на стеке.

   Идея: большинство вызовов арифметики работает на малых prec,
   где 3-5 malloc/calloc на операцию перевешивают саму арифметику.
   Обёртка s21_mpf_stack_t содержит сам mpf и стековый буфер; при
   prec, влезающих в S21_MPF_STACK_LIMBS лимбов, limbs указывает
   на stack_buf, иначе — на calloc'нутый heap_limbs.

   Использование:
     s21_mpf_stack_t tmp;
     s21_mpf_stack_init(&tmp, prec);
     ...работать с tmp.mpf...
     s21_mpf_stack_clear(&tmp);

   ВАЖНО: такие объекты нельзя освобождать через s21_mpf_clear —
   free() попытается освободить указатель на стек. Только через
   s21_mpf_stack_clear.
   ------------------------------------------------------------ */

#define S21_MPF_STACK_LIMBS 32

typedef struct {
  s21_mpf_t mpf;
  uint64_t *heap_limbs; /* NULL, если limbs указывает на stack_buf */
  uint64_t stack_buf[S21_MPF_STACK_LIMBS];
} s21_mpf_stack_t;

void s21_mpf_stack_init(s21_mpf_stack_t *st, uint32_t prec);
void s21_mpf_stack_clear(s21_mpf_stack_t *st);

/* ------------------------------------------------------------
   Проверки аргументов. 0 — успех, -1 — ошибка.
   ------------------------------------------------------------ */

static inline int s21_mpf_check_binary(const s21_mpf_t *res, const s21_mpf_t *x,
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
