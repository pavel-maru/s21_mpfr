#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

/* ============================================================
   Квадратный корень
   ============================================================ */

/* Начальное приближение для sqrt из double.

   a — NORMAL и положительное.  seed — mpf (уже инициализирован
   на нужную точность).

   Возвращает 0 при успехе, -1 если double-представление
   неприменимо (a.exp слишком экстремален).  При -1 caller
   использует fallback 2^(a.exp/2 + 1). */
static int sqrt_seed_from_double(s21_mpf_t *seed, const s21_mpf_t *a) {
  size_t count = s21_mpf_limbs_for_prec(a->prec);
  int last = (int)count - 1;
  while (last > 0 && a->limbs[last] == 0) last--;
  if (a->limbs[last] == 0) return -1;

  /* top64 — 64 старших бита мантиссы, MSB на позиции 63. */
  int msb_local = 63 - __builtin_clzll(a->limbs[last]);
  int shift = 63 - msb_local;
  uint64_t top64 = a->limbs[last] << shift;

  /* m_norm = top64 / 2^63 ∈ [1, 2). */
  double m_norm = (double)top64 / 9223372036854775808.0;

  /* a = m_norm · 2^(exp-1).  sqrt(a) = sqrt(m_norm·2^rem) · 2^half,
     где E = exp-1, half = floor(E/2), rem = E - 2·half ∈ {0,1}. */
  int64_t E = a->exp - 1;
  int64_t half = E / 2;
  if (E < 0 && (E & 1)) half--;
  int64_t rem = E - 2 * half;

  double scaled = m_norm * (rem ? 2.0 : 1.0);
  double root = sqrt(scaled);

  if (half > INT32_MAX || half < INT32_MIN) return -1;
  double result = ldexp(root, (int)half);

  if (!isfinite(result) || result <= 0.0) return -1;

  s21_mpf_set_d(seed, result);
  return 0;
}

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

  if (sqrt_seed_from_double(&x_cur, &a) != 0) {
    /* Fallback: 2^(a.exp/2 + 1). */
    int64_t half_exp = a.exp / 2;
    s21_mpf_set_ui(&x_cur, 1);
    x_cur.exp = half_exp + 1;
    s21_mpf_normalize(&x_cur);
  }

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

  /* Halving до |y| ≤ 2^-5.  При |y| = 2^-5 ряд exp(y) сходится
     за ~wp/5 членов вместо ~wp при |y| < 1.  Дополнительные 4
     halving (по одной mul каждый) компенсируются многократным
     сокращением ряда. */
  int s = 0;
  while (y.exp > -5 && s < 128) {
    y.exp -= 1;
    s++;
  }

  s21_mpf_set_ui(&sum, 1);
  s21_mpf_set_ui(&term, 1);

  /* |y| ≤ 2^-5 после halving, поэтому term_n ≤ 2^-5n и ряд
     гарантированно сходится за wp/5 + небольшой запас членов.
     Значение 4000 соответствует wp до 20000, то есть prec до
     ~20000 бит.  Больше не нужно: при wp → 4000 нам хватает
     ~800 членов с огромным запасом. */
  int max_terms = (int)wp + 10;
  for (int n = 1; n < max_terms; n++) {
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

static void log_reduce_to_unit(s21_mpf_t *y, int64_t *k, const s21_mpf_t *a) {
  *k = a->exp - 1;
  s21_mpf_set(y, a);
  y->exp = 1;
  s21_mpf_normalize(y);
}

/* log(y) = 2^m · log(y^(1/2^m)) для y ∈ [1, 2).

   После m итераций sqrt аргумент atanh-ряда уменьшается в 2^(m+1)
   раз, а z² — в 2^(2m+2) раз.  Ряд сходится за ~wp / (2m+2) членов
   вместо ~wp / log2(9) при m = 0.

   При m = 5 для prec = 4096 ряд сокращается с ~1300 до ~350 членов.
   Стоимость 5 sqrt (~500 мкс суммарно с v1.9.0 double seed) окупается
   многократно: ~1.5 мс чистой экономии.

   Ошибка округления от m sqrt накапливается на уровне 2^-wp (каждый
   sqrt даёт результат с точностью wp), затем масштабируется сдвигом
   экспоненты на m бит — точность сохраняется. */
static void log_unit_interval(s21_mpf_t *log_y, const s21_mpf_t *y,
                              uint32_t wp) {
  const int m = 5;

  s21_mpf_t ym, one, num, den, z;
  s21_mpf_init2(&ym, wp);
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&num, wp);
  s21_mpf_init2(&den, wp);
  s21_mpf_init2(&z, wp);

  s21_mpf_set(&ym, y);
  for (int i = 0; i < m; i++) {
    s21_mpf_sqrt(&ym, &ym);
  }

  s21_mpf_set_ui(&one, 1);
  s21_mpf_sub_raw(&num, &ym, &one);
  s21_mpf_add_raw(&den, &ym, &one);
  s21_mpf_div_raw(&z, &num, &den);

  s21_mpf_atan_series(log_y, &z, wp, 0);
  log_y->exp += 1; /* ×2 из atanh */
  log_y->exp += m; /* ×2^m из sqrt reduction */
  s21_mpf_normalize(log_y);

  s21_mpf_clear(&ym);
  s21_mpf_clear(&one);
  s21_mpf_clear(&num);
  s21_mpf_clear(&den);
  s21_mpf_clear(&z);
}

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
   tan — общая часть tan_impl + публичные обёртки
   ============================================================ */

/* Общая часть: pi уже дано, спецзначения уже обработаны. */
static int tan_impl(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *pi) {
  uint32_t wp = res->prec + 64;
  s21_mpf_t s, c, result;
  s21_mpf_init2(&s, wp);
  s21_mpf_init2(&c, wp);
  s21_mpf_init2(&result, wp);

  s21_mpf_sin_impl(&s, x, pi);
  s21_mpf_cos_impl(&c, x, pi);

  if (s21_mpf_is_zero(&c)) {
    /* tan(π/2 + kπ) = ±inf, знак берём из sin. */
    s21_mpf_set_inf(res, s.sign);
  } else {
    s21_mpf_div_raw(&result, &s, &c);
    s21_mpf_set(res, &result);
  }

  s21_mpf_clear(&s);
  s21_mpf_clear(&c);
  s21_mpf_clear(&result);
  return 0;
}

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

  /* sin_impl и cos_impl работают на wp' = res->prec + 128. */
  uint32_t pi_prec = res->prec + 128;
  s21_mpf_t pi;
  s21_mpf_init2(&pi, pi_prec);
  s21_mpf_compute_pi(&pi, pi_prec);
  int rc = tan_impl(res, x, &pi);
  s21_mpf_clear(&pi);
  return rc;
}

/* ============================================================
   atan — публичный + _impl
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
   acos — публичный, использует asin_impl внутри
   ============================================================ */

int s21_mpf_acos(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
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

      s21_mpf_asin_impl(&as, x, &pi);

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

/* ============================================================
   Публичный API с явной π.

   Контракт: pi->prec >= res->prec + 128.  Меньшая точность π
   даёт -1.  Единое требование для всех шести функций, потому
   что tan_with_pi и asin_with_pi внутри проходят две ступени
   guard-битов (res->prec + 64, затем ещё +64).
   ============================================================ */

/* Общая проверка аргументов для _with_pi. */
static int with_pi_check(const s21_mpf_t *res, const s21_mpf_t *x,
                         const s21_mpf_t *pi) {
  if (res == NULL || x == NULL || pi == NULL) return -1;
  if (pi->prec < res->prec + 128) return -1;
  return 0;
}

int s21_mpf_sin_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                        const s21_mpf_t *pi) {
  if (with_pi_check(res, x, pi)) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }
  s21_mpf_sin_impl(res, x, pi);
  return 0;
}

int s21_mpf_cos_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                        const s21_mpf_t *pi) {
  if (with_pi_check(res, x, pi)) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_ui(res, 1);
    return 0;
  }
  s21_mpf_cos_impl(res, x, pi);
  return 0;
}

int s21_mpf_tan_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                        const s21_mpf_t *pi) {
  if (with_pi_check(res, x, pi)) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }
  return tan_impl(res, x, pi);
}

int s21_mpf_atan_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                         const s21_mpf_t *pi) {
  if (with_pi_check(res, x, pi)) return -1;
  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    /* atan(±inf) = ±π/2, используем переданную π. */
    uint32_t wp = res->prec + 64;
    s21_mpf_t pi_2;
    s21_mpf_init2(&pi_2, wp);
    s21_mpf_set(&pi_2, pi);
    pi_2.exp -= 1;
    s21_mpf_normalize(&pi_2);
    pi_2.sign = x->sign;
    s21_mpf_set(res, &pi_2);
    s21_mpf_clear(&pi_2);
    return 0;
  }
  s21_mpf_atan_impl(res, x, pi);
  return 0;
}

int s21_mpf_asin_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                         const s21_mpf_t *pi) {
  if (with_pi_check(res, x, pi)) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }

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
    s21_mpf_t pi_2;
    s21_mpf_init2(&pi_2, wp);
    s21_mpf_set(&pi_2, pi);
    pi_2.exp -= 1;
    s21_mpf_normalize(&pi_2);
    pi_2.sign = x->sign;
    s21_mpf_set(res, &pi_2);
    s21_mpf_clear(&pi_2);
    return 0;
  }
  s21_mpf_asin_impl(res, x, pi);
  return 0;
}

int s21_mpf_acos_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                         const s21_mpf_t *pi) {
  if (with_pi_check(res, x, pi)) return -1;
  if (x->kind == S21_MPF_NAN || x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }

  uint32_t wp = res->prec + 64;
  s21_mpf_t a, one, pi_2;
  s21_mpf_init2(&a, wp);
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&pi_2, wp);
  s21_mpf_set_ui(&one, 1);

  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set(&pi_2, pi);
    pi_2.exp -= 1;
    s21_mpf_normalize(&pi_2);
    s21_mpf_set(res, &pi_2);
    s21_mpf_clear(&a);
    s21_mpf_clear(&one);
    s21_mpf_clear(&pi_2);
    return 0;
  }

  s21_mpf_abs(&a, x);
  int cmp = s21_mpf_cmp(&a, &one);

  if (cmp > 0) {
    s21_mpf_set_nan(res);
  } else if (cmp == 0) {
    if (x->sign)
      s21_mpf_set(res, pi);
    else
      s21_mpf_set_zero(res, 0);
  } else {
    s21_mpf_t as, result;
    s21_mpf_init2(&as, wp);
    s21_mpf_init2(&result, wp);

    s21_mpf_asin_impl(&as, x, pi);

    s21_mpf_set(&pi_2, pi);
    pi_2.exp -= 1;
    s21_mpf_normalize(&pi_2);
    s21_mpf_sub_raw(&result, &pi_2, &as);
    s21_mpf_set(res, &result);

    s21_mpf_clear(&as);
    s21_mpf_clear(&result);
  }

  s21_mpf_clear(&a);
  s21_mpf_clear(&one);
  s21_mpf_clear(&pi_2);
  return 0;
}
