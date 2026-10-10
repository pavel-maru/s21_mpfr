#include <stdlib.h>
#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

/* ============================================================
   Квадратный корень
   ============================================================ */

int s21_mpf_sqrt(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, 0);
    return 0;
  }
  if (x->sign == 1) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, 0);
    return 0;
  }

  uint32_t wp = res->prec + 64;
  s21_mpf_t a, x_cur, x_next, tmp;
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&x_cur, wp);
  s21_mpf_init2(&x_next, wp);
  s21_mpf_init2(&tmp, wp);

  s21_mpf_set(&a, x);

  int64_t half_exp = a.exp / 2;
  s21_mpf_set_ui(&x_cur, 1);
  x_cur.exp = half_exp + 1;
  s21_mpf_normalize(&x_cur);

  for (int i = 0; i < 200; i++) {
    s21_mpf_div_raw(&tmp, &a, &x_cur);
    s21_mpf_add_raw(&x_next, &x_cur, &tmp);
    x_next.exp -= 1;
    s21_mpf_normalize(&x_next);
    if (s21_mpf_cmp(&x_next, &x_cur) == 0) break;
    s21_mpf_set(&x_cur, &x_next);
  }

  s21_mpf_set(res, &x_cur);

  s21_mpf_clear(&a);
  s21_mpf_clear(&x_cur);
  s21_mpf_clear(&x_next);
  s21_mpf_clear(&tmp);
  return 0;
}

/* ============================================================
   Экспонента
   ============================================================ */

int s21_mpf_exp(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_ui(res, 1);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    if (x->sign == 0)
      s21_mpf_set_inf(res, 0);
    else
      s21_mpf_set_zero(res, 0);
    return 0;
  }

  uint32_t wp = res->prec + 64;
  s21_mpf_t a, y, sum, term;
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&y, wp);
  s21_mpf_init2(&sum, wp);
  s21_mpf_init2(&term, wp);

  s21_mpf_set(&a, x);
  s21_mpf_set(&y, &a);

  int s = 0;
  while (y.exp > -1 && s < 128) {
    y.exp -= 1;
    s++;
  }

  s21_mpf_set_ui(&sum, 1);
  s21_mpf_set_ui(&term, 1);

  for (int n = 1; n < 300; n++) {
    s21_mpf_mul_raw(&term, &term, &y);
    s21_mpf_div_small(&term, (uint32_t)n);
    s21_mpf_add_raw(&sum, &sum, &term);

    if (term.kind == S21_MPF_ZERO || term.exp < sum.exp - (int64_t)wp - 4)
      break;
  }

  for (int i = 0; i < s; i++) s21_mpf_mul_raw(&sum, &sum, &sum);

  s21_mpf_set(res, &sum);

  s21_mpf_clear(&a);
  s21_mpf_clear(&y);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&term);
  return 0;
}

/* ============================================================
   Серия для atanh / atan — одна реализация на оба случая.

   alternating == 0 → atanh(z) = z + z³/3 + z⁵/5 + ...
   alternating == 1 → atan(z)  = z - z³/3 + z⁵/5 - ...
   ============================================================ */

static void s21_mpf_atan_series(s21_mpf_t *out, const s21_mpf_t *z, uint32_t wp,
                                int alternating) {
  s21_mpf_t z2, term, sum, tmp;
  s21_mpf_init2(&z2, wp);
  s21_mpf_init2(&term, wp);
  s21_mpf_init2(&sum, wp);
  s21_mpf_init2(&tmp, wp);

  s21_mpf_mul_raw(&z2, z, z);
  s21_mpf_set(&term, z);
  s21_mpf_set(&sum, z);

  for (int n = 3; n < 4000; n += 2) {
    s21_mpf_mul_raw(&term, &term, &z2);
    if (alternating) s21_mpf_neg_raw(&term, &term);

    s21_mpf_set(&tmp, &term);
    s21_mpf_div_small(&tmp, (uint32_t)n);
    s21_mpf_add_raw(&sum, &sum, &tmp);

    if (term.kind == S21_MPF_ZERO || term.exp < sum.exp - (int64_t)wp - 4)
      break;
  }

  s21_mpf_set(out, &sum);
  s21_mpf_clear(&z2);
  s21_mpf_clear(&term);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&tmp);
}

/* ============================================================
   Натуральный логарифм

   log(a) = log(y) + k·ln2, где a = y · 2^k, y ∈ [1, 2).
   log(y) = 2 · atanh((y-1)/(y+1)),
   ln2    = 2 · atanh(1/3).
   ============================================================ */

/* Разложение a = y · 2^k с приведением y к [1, 2). */
static void log_reduce_to_unit(s21_mpf_t *y, int64_t *k, const s21_mpf_t *a) {
  *k = a->exp - 1;
  s21_mpf_set(y, a);
  y->exp = 1;
  s21_mpf_normalize(y);
}

/* log(y) = 2 · atanh((y-1)/(y+1)) для y ∈ [1, 2). */
static void log_unit_interval(s21_mpf_t *log_y, const s21_mpf_t *y,
                              uint32_t wp) {
  s21_mpf_t one, num, den, z;
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&num, wp);
  s21_mpf_init2(&den, wp);
  s21_mpf_init2(&z, wp);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_sub_raw(&num, y, &one);
  s21_mpf_add_raw(&den, y, &one);
  s21_mpf_div_raw(&z, &num, &den);

  s21_mpf_atan_series(log_y, &z, wp, 0);
  log_y->exp += 1;
  s21_mpf_normalize(log_y);

  s21_mpf_clear(&one);
  s21_mpf_clear(&num);
  s21_mpf_clear(&den);
  s21_mpf_clear(&z);
}

/* ln2 = 2 · atanh(1/3). */
static void log_compute_ln2(s21_mpf_t *ln2, uint32_t wp) {
  s21_mpf_t one, three, z;
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&three, wp);
  s21_mpf_init2(&z, wp);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_set_ui(&three, 3);
  s21_mpf_div_raw(&z, &one, &three);
  s21_mpf_atan_series(ln2, &z, wp, 0);
  ln2->exp += 1;
  s21_mpf_normalize(ln2);

  s21_mpf_clear(&one);
  s21_mpf_clear(&three);
  s21_mpf_clear(&z);
}

int s21_mpf_log(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_inf(res, 1);
    return 0;
  }
  if (x->sign == 1) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, 0);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  int64_t k;
  s21_mpf_t y, log_y;
  s21_mpf_init2(&y, wp);
  s21_mpf_init2(&log_y, wp);

  log_reduce_to_unit(&y, &k, x);
  log_unit_interval(&log_y, &y, wp);

  if (k == 0) {
    s21_mpf_set(res, &log_y);
  } else {
    s21_mpf_t ln2, kk, k_ln2, result;
    s21_mpf_init2(&ln2, wp);
    s21_mpf_init2(&kk, wp);
    s21_mpf_init2(&k_ln2, wp);
    s21_mpf_init2(&result, wp);

    log_compute_ln2(&ln2, wp);
    s21_mpf_set_si(&kk, (long)k);
    s21_mpf_mul_raw(&k_ln2, &kk, &ln2);
    s21_mpf_add_raw(&result, &log_y, &k_ln2);
    s21_mpf_set(res, &result);

    s21_mpf_clear(&ln2);
    s21_mpf_clear(&kk);
    s21_mpf_clear(&k_ln2);
    s21_mpf_clear(&result);
  }

  s21_mpf_clear(&y);
  s21_mpf_clear(&log_y);
  return 0;
}

/* ============================================================
   Утилиты тригонометрии
   ============================================================ */

static void s21_mpf_trunc(s21_mpf_t *res, const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) {
    s21_mpf_set(res, x);
    return;
  }
  if (x->exp <= 0) {
    s21_mpf_set_zero(res, x->sign);
    return;
  }

  s21_mpf_set(res, x);
  if (res->exp >= (int64_t)res->prec) return;

  uint32_t keep = (uint32_t)res->exp;
  for (uint32_t i = 0; i < res->prec - keep; i++) s21_mpf_set_bit(res, i, 0);
  s21_mpf_normalize(res);
}

static void s21_mpf_round_nearest(s21_mpf_t *res, const s21_mpf_t *x) {
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return;
  }

  s21_mpf_t half, shifted;
  s21_mpf_init2(&half, x->prec);
  s21_mpf_init2(&shifted, x->prec);

  s21_mpf_set_d(&half, 0.5);
  if (x->sign) s21_mpf_neg_raw(&half, &half);

  s21_mpf_add_raw(&shifted, x, &half);
  s21_mpf_trunc(res, &shifted);

  s21_mpf_clear(&half);
  s21_mpf_clear(&shifted);
}

/* ============================================================
   Константа π (формула Мэчина)
   ============================================================ */

static void s21_mpf_compute_pi(s21_mpf_t *pi, uint32_t wp) {
  s21_mpf_t one, five, k239, z, a, b, tmp;
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&five, wp);
  s21_mpf_init2(&k239, wp);
  s21_mpf_init2(&z, wp);
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&b, wp);
  s21_mpf_init2(&tmp, wp);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_set_ui(&five, 5);
  s21_mpf_set_ui(&k239, 239);

  s21_mpf_div_raw(&z, &one, &five);
  s21_mpf_atan_series(&a, &z, wp, 1);
  s21_mpf_set_ui(&tmp, 16);
  s21_mpf_mul_raw(&a, &a, &tmp);

  s21_mpf_div_raw(&z, &one, &k239);
  s21_mpf_atan_series(&b, &z, wp, 1);
  s21_mpf_set_ui(&tmp, 4);
  s21_mpf_mul_raw(&b, &b, &tmp);

  s21_mpf_sub_raw(pi, &a, &b);

  s21_mpf_clear(&one);
  s21_mpf_clear(&five);
  s21_mpf_clear(&k239);
  s21_mpf_clear(&z);
  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&tmp);
}

int s21_mpf_pi(s21_mpf_t *res) {
  if (res == NULL) return -1;
  uint32_t wp = res->prec + 64;
  s21_mpf_t pi_w;
  s21_mpf_init2(&pi_w, wp);
  s21_mpf_compute_pi(&pi_w, wp);
  s21_mpf_set(res, &pi_w);
  s21_mpf_clear(&pi_w);
  return 0;
}

/* ============================================================
   Ленивое получение π

   Если pi_in != NULL — возвращаем его, ничего не выделяем.
   Иначе — при первом запросе выделяем pi_local, считаем π и
   помечаем *owned = 1. Повторные запросы с уже установленным
   *owned возвращают ранее вычисленную π без пересчёта.
   ============================================================ */

static const s21_mpf_t *mpf_get_pi(const s21_mpf_t *pi_in, s21_mpf_t *pi_local,
                                   uint32_t wp, int *owned) {
  if (pi_in != NULL) return pi_in;
  if (*owned) return pi_local;
  s21_mpf_init2(pi_local, wp);
  s21_mpf_compute_pi(pi_local, wp);
  *owned = 1;
  return pi_local;
}

/* ============================================================
   Reduction по модулю 2π: r = a - round(a / 2π) · 2π.
   Сводит a к (-π, π].  pi должно быть на точности не меньше wp.
   ============================================================ */

static void reduce_2pi(s21_mpf_t *r, const s21_mpf_t *a, const s21_mpf_t *pi,
                       uint32_t wp) {
  s21_mpf_t two_pi, q, q_r, tmp;
  s21_mpf_init2(&two_pi, wp);
  s21_mpf_init2(&q, wp);
  s21_mpf_init2(&q_r, wp);
  s21_mpf_init2(&tmp, wp);

  s21_mpf_set(&two_pi, pi);
  two_pi.exp += 1;
  s21_mpf_normalize(&two_pi);

  s21_mpf_div_raw(&q, a, &two_pi);
  s21_mpf_round_nearest(&q_r, &q);
  s21_mpf_mul_raw(&tmp, &q_r, &two_pi);
  s21_mpf_neg_raw(&tmp, &tmp);
  s21_mpf_add_raw(r, a, &tmp);

  s21_mpf_clear(&two_pi);
  s21_mpf_clear(&q);
  s21_mpf_clear(&q_r);
  s21_mpf_clear(&tmp);
}

/* ============================================================
   sin — публичный + _impl
   ============================================================ */

void s21_mpf_sin_impl(s21_mpf_t *res, const s21_mpf_t *x,
                      const s21_mpf_t *pi_in) {
  uint32_t wp = res->prec + 64;
  s21_mpf_t pi_local;
  int pi_owned = 0;
  const s21_mpf_t *pi = mpf_get_pi(pi_in, &pi_local, wp, &pi_owned);

  s21_mpf_t a, r, term, sum, x2;
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&r, wp);
  s21_mpf_init2(&term, wp);
  s21_mpf_init2(&sum, wp);
  s21_mpf_init2(&x2, wp);

  s21_mpf_set(&a, x);
  reduce_2pi(&r, &a, pi, wp);

  s21_mpf_set(&term, &r);
  s21_mpf_set(&sum, &r);
  s21_mpf_mul_raw(&x2, &r, &r);

  for (int n = 1; n < 400; n++) {
    s21_mpf_mul_raw(&term, &term, &x2);
    s21_mpf_neg_raw(&term, &term);
    s21_mpf_div_small(&term, (uint32_t)(2 * n));
    s21_mpf_div_small(&term, (uint32_t)(2 * n + 1));
    s21_mpf_add_raw(&sum, &sum, &term);

    if (term.kind == S21_MPF_ZERO || term.exp < sum.exp - (int64_t)wp - 4)
      break;
  }

  s21_mpf_set(res, &sum);

  s21_mpf_clear(&a);
  s21_mpf_clear(&r);
  s21_mpf_clear(&term);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&x2);

  if (pi_owned) s21_mpf_clear(&pi_local);
}

int s21_mpf_sin(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }
  s21_mpf_sin_impl(res, x, NULL);
  return 0;
}

/* ============================================================
   cos — публичный + _impl
   ============================================================ */

void s21_mpf_cos_impl(s21_mpf_t *res, const s21_mpf_t *x,
                      const s21_mpf_t *pi_in) {
  uint32_t wp = res->prec + 64;
  s21_mpf_t pi_local;
  int pi_owned = 0;
  const s21_mpf_t *pi = mpf_get_pi(pi_in, &pi_local, wp, &pi_owned);

  s21_mpf_t a, r, term, sum, x2, one;
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&r, wp);
  s21_mpf_init2(&term, wp);
  s21_mpf_init2(&sum, wp);
  s21_mpf_init2(&x2, wp);
  s21_mpf_init2(&one, wp);

  s21_mpf_set(&a, x);
  reduce_2pi(&r, &a, pi, wp);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_set(&term, &one);
  s21_mpf_set(&sum, &one);
  s21_mpf_mul_raw(&x2, &r, &r);

  for (int n = 1; n < 400; n++) {
    s21_mpf_mul_raw(&term, &term, &x2);
    s21_mpf_neg_raw(&term, &term);
    s21_mpf_div_small(&term, (uint32_t)(2 * n - 1));
    s21_mpf_div_small(&term, (uint32_t)(2 * n));
    s21_mpf_add_raw(&sum, &sum, &term);

    if (term.kind == S21_MPF_ZERO || term.exp < sum.exp - (int64_t)wp - 4)
      break;
  }

  s21_mpf_set(res, &sum);

  s21_mpf_clear(&a);
  s21_mpf_clear(&r);
  s21_mpf_clear(&term);
  s21_mpf_clear(&sum);
  s21_mpf_clear(&x2);
  s21_mpf_clear(&one);

  if (pi_owned) s21_mpf_clear(&pi_local);
}

int s21_mpf_cos(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_ui(res, 1);
    return 0;
  }
  s21_mpf_cos_impl(res, x, NULL);
  return 0;
}

/* ============================================================
   tan — одна π на sin и cos
   ============================================================ */

int s21_mpf_tan(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }

  uint32_t wp = res->prec + 64;
  /* sin_impl и cos_impl работают на wp' = s.prec + 64. Чтобы π
     хватило на весь расчёт без потери guard-битов, берём её на
     prec = res->prec + 128. */
  uint32_t pi_prec = res->prec + 128;

  s21_mpf_t pi, s, c, result;
  s21_mpf_init2(&pi, pi_prec);
  s21_mpf_init2(&s, wp);
  s21_mpf_init2(&c, wp);
  s21_mpf_init2(&result, wp);

  s21_mpf_compute_pi(&pi, pi_prec);
  s21_mpf_sin_impl(&s, x, &pi);
  s21_mpf_cos_impl(&c, x, &pi);

  if (s21_mpf_is_zero(&c)) {
    /* tan(π/2 + kπ) = ±inf, знак берём из sin. */
    s21_mpf_set_inf(res, s.sign);
  } else {
    s21_mpf_div_raw(&result, &s, &c);
    s21_mpf_set(res, &result);
  }

  s21_mpf_clear(&pi);
  s21_mpf_clear(&s);
  s21_mpf_clear(&c);
  s21_mpf_clear(&result);
  return 0;
}

/* ============================================================
   atan — публичный + _impl

   _impl принимает только NORMAL. π нужна только для |x| > 1;
   для |x| ≤ 1 не вычисляется вообще.
   ============================================================ */

void s21_mpf_atan_impl(s21_mpf_t *res, const s21_mpf_t *x,
                       const s21_mpf_t *pi_in) {
  uint32_t wp = res->prec + 64;
  s21_mpf_t pi_local;
  int pi_owned = 0;

  s21_mpf_t one, a, sq;
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&sq, wp);
  s21_mpf_set_ui(&one, 1);

  int sign = x->sign;
  s21_mpf_abs(&a, x);

  int invert = 0;
  const s21_mpf_t *pi = NULL;
  if (s21_mpf_cmp(&a, &one) > 0) {
    /* π нужна только здесь. */
    pi = mpf_get_pi(pi_in, &pi_local, wp, &pi_owned);
    s21_mpf_div_raw(&a, &one, &a);
    invert = 1;
  }

  int reductions = 0;
  while (a.exp > -2 && reductions < 16) {
    s21_mpf_mul_raw(&sq, &a, &a);
    s21_mpf_add_raw(&sq, &sq, &one);
    s21_mpf_sqrt(&sq, &sq);
    s21_mpf_add_raw(&sq, &sq, &one);
    s21_mpf_div_raw(&a, &a, &sq);
    reductions++;
  }

  s21_mpf_atan_series(&a, &a, wp, 1);
  a.exp += reductions;
  s21_mpf_normalize(&a);

  if (invert) {
    s21_mpf_t pi_2;
    s21_mpf_init2(&pi_2, wp);
    s21_mpf_set(&pi_2, pi);
    pi_2.exp -= 1;
    s21_mpf_normalize(&pi_2);
    s21_mpf_sub_raw(&a, &pi_2, &a);
    s21_mpf_clear(&pi_2);
  }

  if (sign) a.sign = 1;
  s21_mpf_set(res, &a);

  s21_mpf_clear(&one);
  s21_mpf_clear(&a);
  s21_mpf_clear(&sq);
  if (pi_owned) s21_mpf_clear(&pi_local);
}

int s21_mpf_atan(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    /* atan(±inf) = ±π/2 */
    uint32_t wp = res->prec + 64;
    s21_mpf_t pi;
    s21_mpf_init2(&pi, wp);
    s21_mpf_compute_pi(&pi, wp);
    pi.exp -= 1;
    s21_mpf_normalize(&pi);
    pi.sign = x->sign;
    s21_mpf_set(res, &pi);
    s21_mpf_clear(&pi);
    return 0;
  }
  s21_mpf_atan_impl(res, x, NULL);
  return 0;
}

/* ============================================================
   asin — публичный + _impl

   _impl принимает только NORMAL, |x| < 1.  Для |x| = 1 обёртка
   отдельно возвращает ±π/2.
   ============================================================ */

void s21_mpf_asin_impl(s21_mpf_t *res, const s21_mpf_t *x,
                       const s21_mpf_t *pi_in) {
  uint32_t wp = res->prec + 64;
  s21_mpf_t a, one, sq, num, den, at;
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&sq, wp);
  s21_mpf_init2(&num, wp);
  s21_mpf_init2(&den, wp);
  s21_mpf_init2(&at, wp);

  s21_mpf_set(&a, x);
  s21_mpf_set_ui(&one, 1);

  s21_mpf_mul_raw(&sq, &a, &a);    /* x² */
  s21_mpf_sub_raw(&sq, &one, &sq); /* 1 - x² */
  s21_mpf_sqrt(&den, &sq);         /* √(1 - x²) */
  s21_mpf_div_raw(&num, &a, &den); /* x / √(1 - x²) */
  s21_mpf_atan_impl(&at, &num, pi_in);
  s21_mpf_set(res, &at);

  s21_mpf_clear(&a);
  s21_mpf_clear(&one);
  s21_mpf_clear(&sq);
  s21_mpf_clear(&num);
  s21_mpf_clear(&den);
  s21_mpf_clear(&at);
}

int s21_mpf_asin(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }

  /* NORMAL: проверяем |x| ≤ 1. */
  s21_mpf_t one;
  s21_mpf_init2(&one, x->prec);
  s21_mpf_set_ui(&one, 1);
  int cmp = s21_mpf_cmp_abs(x, &one);
  s21_mpf_clear(&one);

  if (cmp > 0) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (cmp == 0) {
    /* x = ±1 → asin = ±π/2 */
    uint32_t wp = res->prec + 64;
    s21_mpf_t pi;
    s21_mpf_init2(&pi, wp);
    s21_mpf_compute_pi(&pi, wp);
    pi.exp -= 1;
    s21_mpf_normalize(&pi);
    pi.sign = x->sign;
    s21_mpf_set(res, &pi);
    s21_mpf_clear(&pi);
    return 0;
  }
  s21_mpf_asin_impl(res, x, NULL);
  return 0;
}

/* ============================================================
   acos — одна π на весь вызов
   ============================================================ */

int s21_mpf_acos(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    /* acos(0) = π/2 */
    uint32_t wp = res->prec + 64;
    s21_mpf_t pi;
    s21_mpf_init2(&pi, wp);
    s21_mpf_compute_pi(&pi, wp);
    pi.exp -= 1;
    s21_mpf_normalize(&pi);
    s21_mpf_set(res, &pi);
    s21_mpf_clear(&pi);
    return 0;
  }

  /* NORMAL */
  uint32_t wp = res->prec + 64;
  s21_mpf_t pi, a, one;
  s21_mpf_init2(&pi, wp);
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&one, wp);
  s21_mpf_set_ui(&one, 1);
  s21_mpf_abs(&a, x);

  if (s21_mpf_cmp(&a, &one) > 0) {
    s21_mpf_set_nan(res);
  } else {
    /* Одна π на всё: и для ветки |x| = 1, и для передачи в asin_impl. */
    s21_mpf_compute_pi(&pi, wp);

    if (s21_mpf_cmp(&a, &one) == 0) {
      if (x->sign)
        s21_mpf_set(res, &pi);
      else
        s21_mpf_set_zero(res, 0);
    } else {
      s21_mpf_t as, pi_2, result;
      s21_mpf_init2(&as, wp);
      s21_mpf_init2(&pi_2, wp);
      s21_mpf_init2(&result, wp);

      s21_mpf_asin_impl(&as, x, &pi); /* π идёт вниз по цепочке */

      s21_mpf_set(&pi_2, &pi);
      pi_2.exp -= 1;
      s21_mpf_normalize(&pi_2);
      s21_mpf_sub_raw(&result, &pi_2, &as);
      s21_mpf_set(res, &result);

      s21_mpf_clear(&as);
      s21_mpf_clear(&pi_2);
      s21_mpf_clear(&result);
    }
  }

  s21_mpf_clear(&pi);
  s21_mpf_clear(&a);
  s21_mpf_clear(&one);
  return 0;
}
