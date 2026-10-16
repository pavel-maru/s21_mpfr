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
   Евклидова норма hypot(x, y) = √(x² + y²).

   Прямое вычисление x² + y² переполняется при |x| ~ 2^(prec/2)
   и теряет точность при разнице порядков.  Масштабируем на
   m = max(|x|, |y|):

     hypot(x, y) = m · √((x/m)² + (y/m)²)

   При m = 0 возвращаем +0.  Спецзначения — по C99:
   hypot(±inf, anything) = +inf (даже если второй аргумент NaN);
   hypot(NaN, y ≠ ±inf) = NaN. */
int s21_mpf_hypot(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;

  /* ±inf доминирует над NaN — это стандартное поведение hypot. */
  if (x->kind == S21_MPF_INF || y->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, 0);
    return 0;
  }
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO && y->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, 0);
    return 0;
  }
  /* Один из аргументов — ноль: hypot(x, 0) = |x|, точно. */
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_abs(res, y);
    return 0;
  }
  if (y->kind == S21_MPF_ZERO) {
    s21_mpf_abs(res, x);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t ax_s, ay_s, m_s, xm_s, ym_s, sum_s, sq_s, root_s;
  s21_mpf_stack_init(&ax_s,   wp);
  s21_mpf_stack_init(&ay_s,   wp);
  s21_mpf_stack_init(&m_s,    wp);
  s21_mpf_stack_init(&xm_s,   wp);
  s21_mpf_stack_init(&ym_s,   wp);
  s21_mpf_stack_init(&sum_s,  wp);
  s21_mpf_stack_init(&sq_s,   wp);
  s21_mpf_stack_init(&root_s, wp);

  s21_mpf_abs(&ax_s.mpf, x);
  s21_mpf_abs(&ay_s.mpf, y);

  if (s21_mpf_cmp(&ax_s.mpf, &ay_s.mpf) >= 0)
    s21_mpf_set(&m_s.mpf, &ax_s.mpf);
  else
    s21_mpf_set(&m_s.mpf, &ay_s.mpf);

  s21_mpf_div_raw(&xm_s.mpf, &ax_s.mpf, &m_s.mpf);
  s21_mpf_div_raw(&ym_s.mpf, &ay_s.mpf, &m_s.mpf);
  s21_mpf_mul_raw(&sq_s.mpf,  &xm_s.mpf, &xm_s.mpf);
  s21_mpf_mul_raw(&sum_s.mpf, &ym_s.mpf, &ym_s.mpf);
  s21_mpf_add_raw(&sum_s.mpf, &sum_s.mpf, &sq_s.mpf);
  s21_mpf_sqrt(&root_s.mpf, &sum_s.mpf);
  s21_mpf_mul_raw(&sum_s.mpf, &root_s.mpf, &m_s.mpf);
  s21_mpf_set(res, &sum_s.mpf);

  s21_mpf_stack_clear(&ax_s);
  s21_mpf_stack_clear(&ay_s);
  s21_mpf_stack_clear(&m_s);
  s21_mpf_stack_clear(&xm_s);
  s21_mpf_stack_clear(&ym_s);
  s21_mpf_stack_clear(&sum_s);
  s21_mpf_stack_clear(&sq_s);
  s21_mpf_stack_clear(&root_s);
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
   log2(x) через reduction по степеням двойки.

   x = 2^k · y,  y ∈ [1, 2),  log2(x) = k + log2(y).

   Разложение по 2^k — точное, потому что 2^k представимо
   в бинарном mpf.  Для x = 2^m: k = m, y = 1, log2(1) = 0,
   результат = m точно (без округления).

   log2(y) = ln(y) / ln(2), где ln(y) — из log_unit_interval
   (тот же atanh-ряд с sqrt-reduction, что и в log).

   Стоимость: один log_unit_interval + ln2 + одно деление.
   Это дешевле, чем полный log(x) / ln(2), потому что
   log_reduce_to_unit уже сделан. */
int s21_mpf_log2(s21_mpf_t *res, const s21_mpf_t *x) {
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
  s21_mpf_t y, log_y, ln2, quot, kk;
  s21_mpf_init2(&y, wp);
  s21_mpf_init2(&log_y, wp);
  s21_mpf_init2(&ln2, wp);
  s21_mpf_init2(&quot, wp);
  s21_mpf_init2(&kk, wp);

  log_reduce_to_unit(&y, &k, x);
  log_unit_interval(&log_y, &y, wp);
  log_compute_ln2(&ln2, wp);
  s21_mpf_div_raw(&quot, &log_y, &ln2);
  s21_mpf_set_si(&kk, (long)k);
  s21_mpf_add_raw(&quot, &quot, &kk);
  s21_mpf_set(res, &quot);

  s21_mpf_clear(&y);
  s21_mpf_clear(&log_y);
  s21_mpf_clear(&ln2);
  s21_mpf_clear(&quot);
  s21_mpf_clear(&kk);
  return 0;
}

/* ============================================================
   log10(x) = log(x) / ln(10).

   Reduction по 10^k невозможен точно: 10^k содержит множитель
   5^k и не представимо в бинарном mpf.  Поэтому — прямое
   деление; ln(10) вычисляется через s21_mpf_log от x = 10. */
int s21_mpf_log10(s21_mpf_t *res, const s21_mpf_t *x) {
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

  s21_mpf_t ten, log_x, log_ten, quot;
  s21_mpf_init2(&ten, wp);
  s21_mpf_init2(&log_x, wp);
  s21_mpf_init2(&log_ten, wp);
  s21_mpf_init2(&quot, wp);

  s21_mpf_set_ui(&ten, 10);
  s21_mpf_log(&log_x, x);
  s21_mpf_log(&log_ten, &ten);
  s21_mpf_div_raw(&quot, &log_x, &log_ten);
  s21_mpf_set(res, &quot);

  s21_mpf_clear(&ten);
  s21_mpf_clear(&log_x);
  s21_mpf_clear(&log_ten);
  s21_mpf_clear(&quot);
  return 0;
}

/* ============================================================
   log(1 + x) и exp(x) - 1 с повышенной точностью у нуля.

   Обе функции страдают от cancelation: при малых |x|
   1 + x или exp(x) - 1 теряют ~log2(1/|x|) бит.  Добавляем
   guard-биты, пропорциональные -x->exp.

   Для |x| < 2^(-prec/2) ряд Тейлора отклоняется от x меньше
   чем на prec бит — возвращаем x напрямую.
   ============================================================ */

int s21_mpf_log1p(s21_mpf_t *res, const s21_mpf_t *x) {
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
    if (x->sign) s21_mpf_set_nan(res);
    else         s21_mpf_set_inf(res, 0);
    return 0;
  }

  /* Домен: x > -1.  Проверка через cmp_abs против 1. */
  if (x->sign) {
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
      s21_mpf_set_inf(res, 1);
      return 0;
    }
  }

  /* |x| < 2^(-prec/2): log1p(x) = x с точностью выше prec. */
  if (x->exp < -(int64_t)res->prec / 2) {
    s21_mpf_set(res, x);
    return 0;
  }

  int64_t guard = (x->exp < 0) ? -x->exp : 0;
  if (guard > S21_MPF_MAX_ARG_BITS) {
    s21_mpf_set(res, x);
    return 0;
  }
  uint32_t wp = res->prec + 64 + (uint32_t)guard;

  s21_mpf_stack_t t_s, one_s;
  s21_mpf_stack_init(&t_s, wp);
  s21_mpf_stack_init(&one_s, wp);

  s21_mpf_set(&t_s.mpf, x);
  s21_mpf_set_ui(&one_s.mpf, 1);
  s21_mpf_add_raw(&t_s.mpf, &t_s.mpf, &one_s.mpf); /* 1 + x */
  s21_mpf_log(res, &t_s.mpf);

  s21_mpf_stack_clear(&t_s);
  s21_mpf_stack_clear(&one_s);
  return 0;
}

int s21_mpf_expm1(s21_mpf_t *res, const s21_mpf_t *x) {
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
    if (x->sign) s21_mpf_set_si(res, -1);
    else         s21_mpf_set_inf(res, 0);
    return 0;
  }

  /* |x| < 2^(-prec/2): expm1(x) = x. */
  if (x->exp < -(int64_t)res->prec / 2) {
    s21_mpf_set(res, x);
    return 0;
  }

  int64_t guard = (x->exp < 0) ? -x->exp : 0;
  if (guard > S21_MPF_MAX_ARG_BITS) {
    s21_mpf_set(res, x);
    return 0;
  }
  uint32_t wp = res->prec + 64 + (uint32_t)guard;

  s21_mpf_stack_t t_s, one_s;
  s21_mpf_stack_init(&t_s, wp);
  s21_mpf_stack_init(&one_s, wp);

  s21_mpf_set(&t_s.mpf, x);
  s21_mpf_exp(&t_s.mpf, &t_s.mpf);
  s21_mpf_set_ui(&one_s.mpf, 1);
  s21_mpf_sub_raw(&t_s.mpf, &t_s.mpf, &one_s.mpf);
  s21_mpf_set(res, &t_s.mpf);

  s21_mpf_stack_clear(&t_s);
  s21_mpf_stack_clear(&one_s);
  return 0;
}

/* ============================================================
   Кубический корень через метод Ньютона.

   Итерация:  y_{n+1} = (2·y_n + x / y_n²) / 3.

   Сходимость квадратичная при y_0, близком к ∛x.  Seed из
   double cbrt(x), если применим; иначе fallback 2^(a.exp/3).

   Работает для обоих знаков: ∛x для x < 0 — отрицательное
   число, Newton сохраняет знак в силу 2·y + x/y² (при y < 0
   и x < 0 оба слагаемых отрицательны).
   ============================================================ */

/* Начальное приближение ∛a из double.
   a — NORMAL, a != 0.  Возвращает 0 при успехе, -1 если
   double-приближение неприменимо (a.exp слишком экстремален). */
static int cbrt_seed_from_double(s21_mpf_t *seed, const s21_mpf_t *a) {
  size_t count = s21_mpf_limbs_for_prec(a->prec);
  int last = (int)count - 1;
  while (last > 0 && a->limbs[last] == 0) last--;
  if (a->limbs[last] == 0) return -1;

  int msb_local = 63 - __builtin_clzll(a->limbs[last]);
  int shift = 63 - msb_local;
  uint64_t top64 = a->limbs[last] << shift;

  /* m_norm = top64 / 2^63 ∈ [1, 2), знак отдельно. */
  double m_norm = (double)top64 / 9223372036854775808.0;

  /* a = ± m_norm · 2^(exp-1).  ∛a = ∛m_norm · 2^((exp-1)/3).
     E = exp - 1,  E = 3·k + r,  r ∈ {0, 1, 2}. */
  int64_t E = a->exp - 1;
  int64_t k = E / 3;
  int64_t r = E - 3 * k;
  /* При E < 0 C-деление округляет к нулю; нормализуем остаток. */
  if (r < 0) {
    r += 3;
    k -= 1;
  }

  double scaled = m_norm;
  if (r == 1) scaled *= 2.0;
  else if (r == 2) scaled *= 4.0;

  double root = cbrt(scaled);
  if (root <= 0.0) return -1;
  if (a->sign) root = -root;

  if (k > INT32_MAX || k < INT32_MIN) return -1;
  double result = ldexp(root, (int)k);

  if (!isfinite(result) || result == 0.0) return -1;

  s21_mpf_set_d(seed, result);
  return 0;
}

int s21_mpf_cbrt(s21_mpf_t *res, const s21_mpf_t *x) {
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
    s21_mpf_set_inf(res, x->sign);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t a_s, y_s, y2_s, tmp_s, three_s, two_s, next_s;
  s21_mpf_stack_init(&a_s,     wp);
  s21_mpf_stack_init(&y_s,     wp);
  s21_mpf_stack_init(&y2_s,    wp);
  s21_mpf_stack_init(&tmp_s,   wp);
  s21_mpf_stack_init(&three_s, wp);
  s21_mpf_stack_init(&two_s,   wp);
  s21_mpf_stack_init(&next_s,  wp);

  s21_mpf_set(&a_s.mpf, x);

  if (cbrt_seed_from_double(&y_s.mpf, &a_s.mpf) != 0) {
    /* Fallback: 2^(exp/3 + 1), знак из a. */
    int64_t k = a_s.mpf.exp / 3;
    s21_mpf_set_ui(&y_s.mpf, 1);
    y_s.mpf.exp = k + 1;
    y_s.mpf.sign = a_s.mpf.sign;
    s21_mpf_normalize(&y_s.mpf);
  }

  s21_mpf_set_ui(&two_s.mpf, 2);
  s21_mpf_set_ui(&three_s.mpf, 3);

  for (int i = 0; i < 200; i++) {
    /* y2 = y² */
    s21_mpf_mul_raw(&y2_s.mpf, &y_s.mpf, &y_s.mpf);
    /* tmp = a / y² */
    s21_mpf_div_raw(&tmp_s.mpf, &a_s.mpf, &y2_s.mpf);
    /* next = 2y + tmp */
    s21_mpf_mul_raw(&next_s.mpf, &y_s.mpf, &two_s.mpf);
    s21_mpf_add_raw(&next_s.mpf, &next_s.mpf, &tmp_s.mpf);
    /* next = next / 3 */
    s21_mpf_div_raw(&next_s.mpf, &next_s.mpf, &three_s.mpf);

    if (s21_mpf_cmp(&next_s.mpf, &y_s.mpf) == 0) break;
    s21_mpf_set(&y_s.mpf, &next_s.mpf);
  }

  s21_mpf_set(res, &y_s.mpf);

  s21_mpf_stack_clear(&a_s);
  s21_mpf_stack_clear(&y_s);
  s21_mpf_stack_clear(&y2_s);
  s21_mpf_stack_clear(&tmp_s);
  s21_mpf_stack_clear(&three_s);
  s21_mpf_stack_clear(&two_s);
  s21_mpf_stack_clear(&next_s);
  return 0;
}

/* ============================================================
   Округление и остаток.

   s21_mpf_trunc — публичная, экспонирована в s21_mpf.h.
   s21_mpf_round_nearest — внутренняя, используется в reduce_2pi
   (тригонометрия).
   ============================================================ */

int s21_mpf_trunc(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind != S21_MPF_NORMAL) {
    s21_mpf_set(res, x);
    return 0;
  }
  if (x->exp <= 0) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }

  s21_mpf_set(res, x);
  if (res->exp >= (int64_t)res->prec) return 0;

  uint32_t keep = (uint32_t)res->exp;
  for (uint32_t i = 0; i < res->prec - keep; i++) s21_mpf_set_bit(res, i, 0);
  s21_mpf_normalize(res);
  return 0;
}

/* Есть ли у x ненулевая дробная часть?
   Предполагается x->kind == S21_MPF_NORMAL. */
static int mpf_has_fraction(const s21_mpf_t *x) {
  if (x->exp <= 0) return 1;
  if (x->exp >= (int64_t)x->prec) return 0;
  uint32_t keep = (uint32_t)x->exp;
  for (uint32_t i = 0; i < x->prec - keep; i++) {
    if (s21_mpf_get_bit(x, i)) return 1;
  }
  return 0;
}

int s21_mpf_floor(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind != S21_MPF_NORMAL) {
    s21_mpf_set(res, x);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t t_s, one_s;
  s21_mpf_stack_init(&t_s, wp);
  s21_mpf_stack_init(&one_s, wp);

  s21_mpf_set(&t_s.mpf, x);
  int has_frac = mpf_has_fraction(&t_s.mpf);

  s21_mpf_trunc(&t_s.mpf, &t_s.mpf);

  if (has_frac && x->sign) {
    s21_mpf_set_ui(&one_s.mpf, 1);
    s21_mpf_sub_raw(&t_s.mpf, &t_s.mpf, &one_s.mpf);
  }
  s21_mpf_set(res, &t_s.mpf);

  s21_mpf_stack_clear(&t_s);
  s21_mpf_stack_clear(&one_s);
  return 0;
}

int s21_mpf_ceil(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind != S21_MPF_NORMAL) {
    s21_mpf_set(res, x);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t t_s, one_s;
  s21_mpf_stack_init(&t_s, wp);
  s21_mpf_stack_init(&one_s, wp);

  s21_mpf_set(&t_s.mpf, x);
  int has_frac = mpf_has_fraction(&t_s.mpf);

  s21_mpf_trunc(&t_s.mpf, &t_s.mpf);

  if (has_frac && !x->sign) {
    s21_mpf_set_ui(&one_s.mpf, 1);
    s21_mpf_add_raw(&t_s.mpf, &t_s.mpf, &one_s.mpf);
  }
  s21_mpf_set(res, &t_s.mpf);

  s21_mpf_stack_clear(&t_s);
  s21_mpf_stack_clear(&one_s);
  return 0;
}

int s21_mpf_round(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind != S21_MPF_NORMAL) {
    s21_mpf_set(res, x);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t t_s, half_s;
  s21_mpf_stack_init(&t_s, wp);
  s21_mpf_stack_init(&half_s, wp);

  s21_mpf_set(&t_s.mpf, x);
  s21_mpf_set_d(&half_s.mpf, 0.5);
  if (x->sign) half_s.mpf.sign = 1;

  s21_mpf_add_raw(&t_s.mpf, &t_s.mpf, &half_s.mpf);
  s21_mpf_trunc(&t_s.mpf, &t_s.mpf);
  s21_mpf_set(res, &t_s.mpf);

  s21_mpf_stack_clear(&t_s);
  s21_mpf_stack_clear(&half_s);
  return 0;
}

int s21_mpf_fmod(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;

  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (y->kind == S21_MPF_INF) {
    s21_mpf_set(res, x);
    return 0;
  }
  if (y->kind == S21_MPF_ZERO) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    s21_mpf_set_zero(res, x->sign);
    return 0;
  }

  /* Оба NORMAL.  Если |x| < |y|, fmod = x. */
  if (s21_mpf_cmp_abs(x, y) < 0) {
    s21_mpf_set(res, x);
    return 0;
  }

  /* Guard-точность на cancelation в x - trunc(x/y)*y.
     При x/y ≈ 2^E теряем E бит; добавляем E к рабочей точности. */
  int64_t E = (int64_t)x->exp - (int64_t)y->exp;
  if (E < 0) E = 0;
  if (E > S21_MPF_MAX_ARG_BITS) {
    s21_mpf_set_nan(res);
    return 0;
  }
  uint32_t wp = res->prec + 64 + (uint32_t)E;

  s21_mpf_stack_t xw_s, yw_s, q_s, n_s, prod_s, r_s;
  s21_mpf_stack_init(&xw_s,   wp);
  s21_mpf_stack_init(&yw_s,   wp);
  s21_mpf_stack_init(&q_s,    wp);
  s21_mpf_stack_init(&n_s,    wp);
  s21_mpf_stack_init(&prod_s, wp);
  s21_mpf_stack_init(&r_s,    wp);

  s21_mpf_set(&xw_s.mpf, x);
  s21_mpf_set(&yw_s.mpf, y);

  s21_mpf_div_raw(&q_s.mpf, &xw_s.mpf, &yw_s.mpf);
  s21_mpf_trunc(&n_s.mpf, &q_s.mpf);
  s21_mpf_mul_raw(&prod_s.mpf, &n_s.mpf, &yw_s.mpf);
  s21_mpf_sub_raw(&r_s.mpf, &xw_s.mpf, &prod_s.mpf);

  /* C99: знак fmod = знак x, включая случай r = 0. */
  if (s21_mpf_is_zero(&r_s.mpf)) r_s.mpf.sign = x->sign;

  s21_mpf_set(res, &r_s.mpf);

  s21_mpf_stack_clear(&xw_s);
  s21_mpf_stack_clear(&yw_s);
  s21_mpf_stack_clear(&q_s);
  s21_mpf_stack_clear(&n_s);
  s21_mpf_stack_clear(&prod_s);
  s21_mpf_stack_clear(&r_s);
  return 0;
}

/* ============================================================
   Fused multiply-add: x·y + z с одним округлением.

   Точность.  Стандартный подход — считать на wp = 2·prec + 64:
   произведение двух prec-битных чисел точно укладывается
   в 2·prec бит, сложение с z — до 2·prec + 1 бит, и вся
   накопленная точность затем округляется один раз до prec.

   Спецзначения — по C99 F.10.10.1.  Порядок проверок важен:
   NaN-протокол, затем inf в x/y, затем inf в z, затем
   x·y = ±0, и только в конце общий случай.

   Знак нуля.  Когда x·y = ±0 и z = ±0, результат — по
   правилу IEEE 754 для суммы двух нулей: sign = prod_sign
   AND z_sign.  То есть fma(+0, y, -0) = +0, fma(-0, y, +0)
   = +0, fma(-0, y, -0) = -0.
   ============================================================ */

int s21_mpf_fma(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y,
                const s21_mpf_t *z) {
  if (res == NULL || x == NULL || y == NULL || z == NULL) return -1;

  /* NaN-протокол: любой NaN → NaN. */
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN ||
      z->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }

  /* Тип произведения x·y. */
  int prod_inf = 0;
  int prod_zero = 0;

  if (x->kind == S21_MPF_INF || y->kind == S21_MPF_INF) {
    if (x->kind == S21_MPF_ZERO || y->kind == S21_MPF_ZERO) {
      /* inf · 0 = NaN. */
      s21_mpf_set_nan(res);
      return 0;
    }
    prod_inf = 1;
  } else if (x->kind == S21_MPF_ZERO || y->kind == S21_MPF_ZERO) {
    prod_zero = 1;
  }

  int prod_sign = x->sign ^ y->sign;

  if (prod_inf) {
    /* inf + z: конфликт знаков → NaN. */
    if (z->kind == S21_MPF_INF && z->sign != prod_sign) {
      s21_mpf_set_nan(res);
      return 0;
    }
    s21_mpf_set_inf(res, prod_sign);
    return 0;
  }

  if (z->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, z->sign);
    return 0;
  }

  if (prod_zero) {
    /* x·y = ±0.  Результат — z + (±0). */
    if (z->kind == S21_MPF_ZERO) {
      /* IEEE 754: сумма двух нулей даёт −0 только когда оба −0. */
      int zero_sign = (prod_sign && z->sign) ? 1 : 0;
      s21_mpf_set_zero(res, zero_sign);
      return 0;
    }
    s21_mpf_set(res, z);
    return 0;
  }

  /* Все три NORMAL, x·y ≠ 0.  Точное умножение на 2·prec,
     сложение с z на 2·prec + 64, одно округление до prec. */
  uint32_t wp = 2 * res->prec + 64;

  s21_mpf_stack_t xs, ys, zs, prod_s, sum_s;
  s21_mpf_stack_init(&xs, wp);
  s21_mpf_stack_init(&ys, wp);
  s21_mpf_stack_init(&zs, wp);
  s21_mpf_stack_init(&prod_s, wp);
  s21_mpf_stack_init(&sum_s, wp);

  s21_mpf_set(&xs.mpf, x);
  s21_mpf_set(&ys.mpf, y);
  s21_mpf_set(&zs.mpf, z);

  int rc = s21_mpf_mul_raw(&prod_s.mpf, &xs.mpf, &ys.mpf);
  if (rc == 0) {
    s21_mpf_add_raw(&sum_s.mpf, &prod_s.mpf, &zs.mpf);
    s21_mpf_set(res, &sum_s.mpf); /* сужение wp → prec, RNDN */
  }

  s21_mpf_stack_clear(&xs);
  s21_mpf_stack_clear(&ys);
  s21_mpf_stack_clear(&zs);
  s21_mpf_stack_clear(&prod_s);
  s21_mpf_stack_clear(&sum_s);
  return rc;
}

/* Внутренний аналог trunc для reduce_2pi: возвращает void,
   чтобы не переписывать вызов в round_nearest. */
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

void s21_mpf_compute_pi_machin(s21_mpf_t *pi, uint32_t wp) {
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

/* ============================================================
   π через Borwein / Gauss-Legendre AGM.

   Итерация:

     a_{n+1} = (aₙ + bₙ) / 2
     b_{n+1} = √(aₙ · bₙ)
     t_{n+1} = tₙ − pₙ · (aₙ − a_{n+1})²
     p_{n+1} = 2 · pₙ

   с начальными a₀ = 1, b₀ = 1/√2, t₀ = 1/4, p₀ = 1.

   Сходимость квадратичная: число правильных бит удваивается
   за итерацию.  Для wp = 4096 нужно ~12 итераций; берём
   log₂(wp) + 3 с запасом.

   Внутренняя точность wp + 64 — чтобы накопленная за 15
   итераций ошибка округления не проникла в значащую часть
   после финального приведения к wp. */

/* Инициализация начальных значений AGM:
   a₀ = 1, b₀ = 1/√2, t₀ = 1/4, p₀ = 1.
   Вспомогательные one / two / four создаются и освобождаются
   внутри. */
static void pi_agm_init(s21_mpf_t *a, s21_mpf_t *b, s21_mpf_t *t,
                        s21_mpf_t *p, uint32_t wp_int) {
  s21_mpf_t one, two, four;
  s21_mpf_init2(&one, wp_int);
  s21_mpf_init2(&two, wp_int);
  s21_mpf_init2(&four, wp_int);

  s21_mpf_set_ui(&one, 1);
  s21_mpf_set_ui(&two, 2);
  s21_mpf_set_ui(&four, 4);

  s21_mpf_set(a, &one);
  s21_mpf_set(p, &one);
  s21_mpf_div_raw(t, &one, &four);

  s21_mpf_div_raw(b, &one, &two);
  s21_mpf_sqrt(b, b);

  s21_mpf_clear(&one);
  s21_mpf_clear(&two);
  s21_mpf_clear(&four);
}

/* Итерация AGM:

     a_{n+1} = (aₙ + bₙ) / 2
     b_{n+1} = √(aₙ · bₙ)
     t_{n+1} = tₙ − pₙ · (aₙ − a_{n+1})²
     p_{n+1} = 2 · pₙ

   Сходимость квадратичная: при wp_int = 4096 достаточно ~12
   итераций.  Цикл останавливается либо по max_iter, либо
   раньше, если a и b совпали побитово. */
static void pi_agm_iterate(s21_mpf_t *a, s21_mpf_t *b, s21_mpf_t *t,
                           s21_mpf_t *p, int max_iter, uint32_t wp_int) {
  s21_mpf_t tmp1, tmp2, tmp3;
  s21_mpf_init2(&tmp1, wp_int);
  s21_mpf_init2(&tmp2, wp_int);
  s21_mpf_init2(&tmp3, wp_int);

  for (int i = 0; i < max_iter; i++) {
    /* tmp1 = (a + b) / 2 */
    s21_mpf_add_raw(&tmp1, a, b);
    s21_mpf_div_small(&tmp1, 2);

    /* tmp2 = √(a · b) */
    s21_mpf_mul_raw(&tmp2, a, b);
    s21_mpf_sqrt(&tmp2, &tmp2);

    /* tmp3 = p · (a − tmp1)², вклад в t */
    s21_mpf_sub_raw(&tmp3, a, &tmp1);
    s21_mpf_mul_raw(&tmp3, &tmp3, &tmp3);
    s21_mpf_mul_raw(&tmp3, &tmp3, p);

    s21_mpf_sub_raw(t, t, &tmp3);

    s21_mpf_set(a, &tmp1);
    s21_mpf_set(b, &tmp2);

    /* p *= 2 */
    s21_mpf_add_raw(p, p, p);

    if (s21_mpf_cmp(a, b) == 0) break;
  }

  s21_mpf_clear(&tmp1);
  s21_mpf_clear(&tmp2);
  s21_mpf_clear(&tmp3);
}

/* π = (a + b)² / (4 · t).

   Считается на wp_int, потом сужается до pi->prec через set_round.
   Прямое div_raw в pi сломало бы алгоритм: pi имеет prec = wp,
   а (a+b)² и 4t — wp_int = wp + 64; div_raw читает только
   limbs_for_prec(res->prec) младших лимбов операндов и потерял
   бы верхние 64 бита мантисс. */
static void pi_agm_finalize(s21_mpf_t *pi, const s21_mpf_t *a,
                            const s21_mpf_t *b, const s21_mpf_t *t,
                            uint32_t wp_int) {
  s21_mpf_t tmp1, tmp2, pi_wide;
  s21_mpf_init2(&tmp1, wp_int);
  s21_mpf_init2(&tmp2, wp_int);
  s21_mpf_init2(&pi_wide, wp_int);

  s21_mpf_add_raw(&tmp1, a, b);
  s21_mpf_mul_raw(&tmp1, &tmp1, &tmp1); /* (a+b)² */
  s21_mpf_add_raw(&tmp2, t, t);         /* 2t */
  s21_mpf_add_raw(&tmp2, &tmp2, &tmp2); /* 4t */
  s21_mpf_div_raw(&pi_wide, &tmp1, &tmp2);

  s21_mpf_set(pi, &pi_wide); /* сужение wp_int -> pi->prec, RNDN */

  s21_mpf_clear(&tmp1);
  s21_mpf_clear(&tmp2);
  s21_mpf_clear(&pi_wide);
}

void s21_mpf_compute_pi_agm(s21_mpf_t *pi, uint32_t wp) {
  uint32_t wp_int = wp + 64;

  /* Число итераций: log₂(wp) + 3.  Считаем через сдвиги,
     без libm-зависимости от формата float. */
  int max_iter = 3;
  for (uint32_t q = wp; q > 1; q >>= 1) max_iter++;

  s21_mpf_t a, b, t, p;
  s21_mpf_init2(&a, wp_int);
  s21_mpf_init2(&b, wp_int);
  s21_mpf_init2(&t, wp_int);
  s21_mpf_init2(&p, wp_int);

  pi_agm_init(&a, &b, &t, &p, wp_int);
  pi_agm_iterate(&a, &b, &t, &p, max_iter, wp_int);
  pi_agm_finalize(pi, &a, &b, &t, wp_int);

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&t);
  s21_mpf_clear(&p);
}

int s21_mpf_pi(s21_mpf_t *res) {
  if (res == NULL) return -1;
  uint32_t wp = res->prec + 64;
  s21_mpf_t pi_w;
  s21_mpf_init2(&pi_w, wp);
  s21_mpf_compute_pi_agm(&pi_w, wp);
  s21_mpf_set(res, &pi_w);
  s21_mpf_clear(&pi_w);
  return 0;
}

int s21_mpf_ln2(s21_mpf_t *res) {
  if (res == NULL) return -1;
  uint32_t wp = res->prec + 64;
  s21_mpf_t ln2_w;
  s21_mpf_init2(&ln2_w, wp);
  log_compute_ln2(&ln2_w, wp);
  s21_mpf_set(res, &ln2_w);
  s21_mpf_clear(&ln2_w);
  return 0;
}

int s21_mpf_e(s21_mpf_t *res) {
  if (res == NULL) return -1;
  s21_mpf_t one;
  s21_mpf_init2(&one, res->prec);
  s21_mpf_set_ui(&one, 1);
  int rc = s21_mpf_exp(res, &one);
  s21_mpf_clear(&one);
  return rc;
}

/* ============================================================
   Возведение в степень: pow(x, y) = x^y.

   Спецзначения — по C99 (Annex F) / MPFR:

     pow(any, ±0)  = 1            включая pow(NaN, 0)
     pow(+1, any)  = 1            включая pow(1, NaN), pow(1, ±inf)
     pow(-1, ±inf) = 1
     pow(NaN, y≠0) = NaN
     pow(x≠+1, NaN) = NaN
     pow(x<0, y нецелое) = NaN
     pow(±0, y>0)  = ±0           знак = sign_x & (y mod 2)
     pow(±0, y<0)  = ±inf         то же
     pow(±inf, y>0) = ±inf        знак = sign_x & (y mod 2)
     pow(±inf, y<0) = ±0          то же
     pow(x, ±inf)  = +inf / +0    если |x| > 1 / |x| < 1
                    1             если |x| = 1

   Общий случай:
     x > 0           r = exp(y · log x), знак +
     x < 0, y целое  r = exp(y · log |x|), знак = sign_x & (y mod 2)
   ============================================================ */

/* x — ровно +1 или −1?  Возвращает 1/0, is_neg — знак. */
static int mpf_is_unit(const s21_mpf_t *x, int *is_neg) {
  if (x->kind != S21_MPF_NORMAL) return 0;
  if (x->exp != 1) return 0;

  uint32_t top_word = (x->prec - 1) / 64;
  uint32_t top_bit = (x->prec - 1) % 64;
  size_t count = s21_mpf_limbs_for_prec(x->prec);

  if (x->limbs[top_word] != (1ULL << top_bit)) return 0;
  for (uint32_t i = 0; i < top_word; i++)
    if (x->limbs[i] != 0) return 0;
  for (size_t i = top_word + 1; i < count; i++)
    if (x->limbs[i] != 0) return 0;

  *is_neg = x->sign;
  return 1;
}

/* y — нечётное целое?  Предполагается y->kind == NORMAL, y != 0.

   Чётность определяется битом мантиссы на позиции (prec − exp):
   именно он имеет вес 2^0 в значении y.  При exp > prec младший
   значащий вес ≥ 2^1, значит y гарантированно чётно. */
static int mpf_is_odd_integer(const s21_mpf_t *y) {
  if (y->exp <= 0) return 0;
  if (y->exp > (int64_t)y->prec) return 0;
  uint32_t pos = (uint32_t)((int64_t)y->prec - y->exp);
  return s21_mpf_get_bit(y, pos);
}

int s21_mpf_pow(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y) {
  if (res == NULL || x == NULL || y == NULL) return -1;

  /* pow(anything, ±0) = 1 — раньше всего, включая pow(NaN, 0). */
  if (y->kind == S21_MPF_ZERO) {
    s21_mpf_set_ui(res, 1);
    return 0;
  }

  /* pow(+1, anything) = 1 — раньше NaN-протокола, включая pow(1, NaN). */
  int x_neg_unit;
  int x_unit = mpf_is_unit(x, &x_neg_unit);
  if (x_unit && !x_neg_unit) {
    s21_mpf_set_ui(res, 1);
    return 0;
  }

  /* NaN-протокол. */
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }

  int y_odd = (y->kind == S21_MPF_NORMAL) ? mpf_is_odd_integer(y) : 0;
  int y_int = s21_mpf_integer_p(y);

  /* y = ±inf. */
  if (y->kind == S21_MPF_INF) {
    /* |x| = 1 → 1 (включая x = -1). */
    if (x_unit) {
      s21_mpf_set_ui(res, 1);
      return 0;
    }
    s21_mpf_t one;
    s21_mpf_init2(&one, x->prec);
    s21_mpf_set_ui(&one, 1);
    int cmp = s21_mpf_cmp_abs(x, &one);
    s21_mpf_clear(&one);

    if (cmp > 0) {
      if (y->sign) s21_mpf_set_zero(res, 0);
      else         s21_mpf_set_inf(res, 0);
    } else {
      if (y->sign) s21_mpf_set_inf(res, 0);
      else         s21_mpf_set_zero(res, 0);
    }
    return 0;
  }

  /* x = ±inf. */
  if (x->kind == S21_MPF_INF) {
    int sign = (x->sign && y_odd) ? 1 : 0;
    if (y->sign) s21_mpf_set_zero(res, sign);
    else         s21_mpf_set_inf(res, sign);
    return 0;
  }

  /* x = ±0. */
  if (x->kind == S21_MPF_ZERO) {
    int sign = (x->sign && y_odd) ? 1 : 0;
    if (y->sign) s21_mpf_set_inf(res, sign);
    else         s21_mpf_set_zero(res, sign);
    return 0;
  }

  /* Оба NORMAL.  x < 0 и y нецелое → NaN. */
  if (x->sign && !y_int) {
    s21_mpf_set_nan(res);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t lg_s, yw_s, ylg_s, r_s;
  s21_mpf_stack_init(&lg_s,  wp);
  s21_mpf_stack_init(&yw_s,  wp);
  s21_mpf_stack_init(&ylg_s, wp);
  s21_mpf_stack_init(&r_s,   wp);

  if (x->sign) {
    s21_mpf_t xa;
    s21_mpf_init2(&xa, wp);
    s21_mpf_abs(&xa, x);
    s21_mpf_log(&lg_s.mpf, &xa);
    s21_mpf_clear(&xa);
  } else {
    s21_mpf_log(&lg_s.mpf, x);
  }

  s21_mpf_set(&yw_s.mpf, y);
  s21_mpf_mul_raw(&ylg_s.mpf, &yw_s.mpf, &lg_s.mpf);
  s21_mpf_exp(&r_s.mpf, &ylg_s.mpf);

  if (x->sign) r_s.mpf.sign = y_odd ? 1 : 0;

  s21_mpf_set(res, &r_s.mpf);

  s21_mpf_stack_clear(&lg_s);
  s21_mpf_stack_clear(&yw_s);
  s21_mpf_stack_clear(&ylg_s);
  s21_mpf_stack_clear(&r_s);
  return 0;
}

/* ============================================================
   Целочисленный показатель: pow_ui и pow_si.

   pow_ui использует бинарное возведение в степень:
     result = 1; base = x
     while n > 0:
       if n & 1: result *= base
       n >>= 1
       if n > 0: base *= base

   Точнее общего pow: нет exp(log) round-trip, нет guard-потерь.
   Для x = 2 результат — точно 2^n.

   pow_si для n < 0: 1 / pow_ui(x, |n|).  Для n = LONG_MIN
   используем безопасное взятие модуля через -(n + 1) + 1.
   ============================================================ */

int s21_mpf_pow_ui(s21_mpf_t *res, const s21_mpf_t *x, unsigned long n) {
  if (res == NULL || x == NULL) return -1;

  if (n == 0) {
    s21_mpf_set_ui(res, 1);
    return 0;
  }
  if (n == 1) {
    s21_mpf_set(res, x);
    return 0;
  }

  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }

  int odd_n = (int)(n & 1);

  if (x->kind == S21_MPF_ZERO) {
    int sign = (x->sign && odd_n) ? 1 : 0;
    s21_mpf_set_zero(res, sign);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    int sign = (x->sign && odd_n) ? 1 : 0;
    s21_mpf_set_inf(res, sign);
    return 0;
  }

  /* x — NORMAL.  Бинарное возведение. */
  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t base_s, result_s;
  s21_mpf_stack_init(&base_s,   wp);
  s21_mpf_stack_init(&result_s, wp);

  s21_mpf_set(&base_s.mpf, x);
  s21_mpf_set_ui(&result_s.mpf, 1);

  while (n > 0) {
    if (n & 1) s21_mpf_mul_raw(&result_s.mpf, &result_s.mpf, &base_s.mpf);
    n >>= 1;
    if (n > 0) s21_mpf_mul_raw(&base_s.mpf, &base_s.mpf, &base_s.mpf);
  }

  s21_mpf_set(res, &result_s.mpf);

  s21_mpf_stack_clear(&base_s);
  s21_mpf_stack_clear(&result_s);
  return 0;
}

int s21_mpf_pow_si(s21_mpf_t *res, const s21_mpf_t *x, long n) {
  if (res == NULL || x == NULL) return -1;

  if (n >= 0) return s21_mpf_pow_ui(res, x, (unsigned long)n);

  /* n < 0.  Безопасное |n| даже для LONG_MIN. */
  unsigned long abs_n = (unsigned long)(-(n + 1)) + 1UL;

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t p_s, one_s, r_s;
  s21_mpf_stack_init(&p_s,   wp);
  s21_mpf_stack_init(&one_s, wp);
  s21_mpf_stack_init(&r_s,   wp);

  int rc = s21_mpf_pow_ui(&p_s.mpf, x, abs_n);
  if (rc != 0) {
    s21_mpf_stack_clear(&p_s);
    s21_mpf_stack_clear(&one_s);
    s21_mpf_stack_clear(&r_s);
    return rc;
  }

  /* 1 / p_s: обработка спецзначений до div_raw.
     div_raw ожидает оба операнда NORMAL и совпадающие с res->prec. */
  if (p_s.mpf.kind != S21_MPF_NORMAL) {
    if (p_s.mpf.kind == S21_MPF_ZERO)
      s21_mpf_set_inf(res, p_s.mpf.sign);
    else if (p_s.mpf.kind == S21_MPF_INF)
      s21_mpf_set_zero(res, p_s.mpf.sign);
    else
      s21_mpf_set_nan(res);
    s21_mpf_stack_clear(&p_s);
    s21_mpf_stack_clear(&one_s);
    s21_mpf_stack_clear(&r_s);
    return 0;
  }

  s21_mpf_set_ui(&one_s.mpf, 1);
  s21_mpf_div_raw(&r_s.mpf, &one_s.mpf, &p_s.mpf); /* r_s.prec = wp */
  s21_mpf_set(res, &r_s.mpf);

  s21_mpf_stack_clear(&p_s);
  s21_mpf_stack_clear(&one_s);
  s21_mpf_stack_clear(&r_s);
  return 0;
}

/* ============================================================
   Корень n-й степени: root(x, n) = x^(1/n).

   Спецзначения:
     root(x, 1)          = x
     root(x, 2)          = sqrt(x)
     root(x, 0)          = NaN
     root(+0, n)         = +0
     root(-0, нечётное)  = -0
     root(-0, чётное)    = +0
     root(+inf, n)       = +inf
     root(-inf, нечёт.)  = -inf
     root(-inf, чётное)  = NaN
     root(x<0, чётное)   = NaN

   Общий случай: exp(log|x| / n), знак отдельно. */
int s21_mpf_root(s21_mpf_t *res, const s21_mpf_t *x, unsigned long n) {
  if (res == NULL || x == NULL) return -1;

  if (n == 0) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (n == 1) {
    s21_mpf_set(res, x);
    return 0;
  }
  if (n == 2) return s21_mpf_sqrt(res, x);

  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }

  int odd_n = (int)(n & 1);

  if (x->kind == S21_MPF_ZERO) {
    int sign = (x->sign && odd_n) ? 1 : 0;
    s21_mpf_set_zero(res, sign);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    if (x->sign && !odd_n) {
      s21_mpf_set_nan(res);
    } else {
      s21_mpf_set_inf(res, x->sign);
    }
    return 0;
  }

  /* x — NORMAL.  Отрицательный x требует нечётного n. */
  if (x->sign && !odd_n) {
    s21_mpf_set_nan(res);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t xa_s, log_a_s, n_s, quot_s, result_s;
  s21_mpf_stack_init(&xa_s,     wp);
  s21_mpf_stack_init(&log_a_s,  wp);
  s21_mpf_stack_init(&n_s,      wp);
  s21_mpf_stack_init(&quot_s,   wp);
  s21_mpf_stack_init(&result_s, wp);

  s21_mpf_abs(&xa_s.mpf, x);
  s21_mpf_log(&log_a_s.mpf, &xa_s.mpf);
  s21_mpf_set_ui(&n_s.mpf, n);
  s21_mpf_div_raw(&quot_s.mpf, &log_a_s.mpf, &n_s.mpf);
  s21_mpf_exp(&result_s.mpf, &quot_s.mpf);

  if (x->sign) result_s.mpf.sign = 1;
  s21_mpf_set(res, &result_s.mpf);

  s21_mpf_stack_clear(&xa_s);
  s21_mpf_stack_clear(&log_a_s);
  s21_mpf_stack_clear(&n_s);
  s21_mpf_stack_clear(&quot_s);
  s21_mpf_stack_clear(&result_s);
  return 0;
}

/* ============================================================
   Гиперболические функции.

     sinh(x) = (e^x − e^−x) / 2
     cosh(x) = (e^x + e^−x) / 2
     tanh(x) = sinh(x) / cosh(x)

   Cancelation есть только в sinh и tanh: при малых |x| вычитание
   двух почти равных чисел теряет ~log2(1/|x|) бит.  Компенсируем
   guard-битами, пропорциональными −exp(x).  cosh(x) — оба
   слагаемых одного знака, cancelation нет.

   Короткое замыкание при |x| < 2^(−prec/2): sinh(x) ≈ x,
   tanh(x) ≈ x с точностью выше prec — возвращаем x напрямую.
   ============================================================ */

/* Общий каркас sinh/cosh: r = (e^x ± e^−x) / 2.
   subtract == 1 → sinh, subtract == 0 → cosh. */
static int hyp_exp_pair(s21_mpf_t *res, const s21_mpf_t *x, int subtract) {
  uint32_t base_prec = res->prec;

  /* Короткое замыкание: |x| < 2^(−prec/2) → sinh(x) = x. */
  if (subtract && x->exp < -(int64_t)base_prec / 2) {
    s21_mpf_set(res, x);
    return 0;
  }

  int64_t guard = 0;
  if (subtract && x->exp < 0) {
    guard = -x->exp;
    if (guard > (int64_t)base_prec) guard = (int64_t)base_prec;
  }
  uint32_t wp = base_prec + 64 + (uint32_t)guard;

  s21_mpf_stack_t ep_s, em_s, sum_s;
  s21_mpf_stack_init(&ep_s,  wp);
  s21_mpf_stack_init(&em_s,  wp);
  s21_mpf_stack_init(&sum_s, wp);

  s21_mpf_set(&ep_s.mpf, x);
  s21_mpf_neg_raw(&em_s.mpf, x);

  s21_mpf_exp(&ep_s.mpf, &ep_s.mpf);
  s21_mpf_exp(&em_s.mpf, &em_s.mpf);

  if (subtract)
    s21_mpf_sub_raw(&sum_s.mpf, &ep_s.mpf, &em_s.mpf);
  else
    s21_mpf_add_raw(&sum_s.mpf, &ep_s.mpf, &em_s.mpf);

  sum_s.mpf.exp -= 1; /* × ½ */
  s21_mpf_normalize(&sum_s.mpf);
  s21_mpf_set(res, &sum_s.mpf);

  s21_mpf_stack_clear(&ep_s);
  s21_mpf_stack_clear(&em_s);
  s21_mpf_stack_clear(&sum_s);
  return 0;
}

int s21_mpf_sinh(s21_mpf_t *res, const s21_mpf_t *x) {
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
    s21_mpf_set_inf(res, x->sign);
    return 0;
  }
  return hyp_exp_pair(res, x, 1);
}

int s21_mpf_cosh(s21_mpf_t *res, const s21_mpf_t *x) {
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
    s21_mpf_set_inf(res, 0);
    return 0;
  }
  return hyp_exp_pair(res, x, 0);
}

int s21_mpf_tanh(s21_mpf_t *res, const s21_mpf_t *x) {
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
    s21_mpf_set_ui(res, 1);
    res->sign = x->sign;
    return 0;
  }

  /* Малые |x|: tanh(x) ≈ x с точностью выше prec. */
  if (x->exp < -(int64_t)res->prec / 2) {
    s21_mpf_set(res, x);
    return 0;
  }

  /* Большие |x|: tanh насыщается к ±1; при |x|·log₂(e) > prec
     младший бит единицы уже не различим. */
  if (x->exp > (int64_t)res->prec) {
    s21_mpf_set_ui(res, 1);
    res->sign = x->sign;
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t s_s, c_s, r_s;
  s21_mpf_stack_init(&s_s, wp);
  s21_mpf_stack_init(&c_s, wp);
  s21_mpf_stack_init(&r_s, wp);

  hyp_exp_pair(&s_s.mpf, x, 1);
  hyp_exp_pair(&c_s.mpf, x, 0);
  s21_mpf_div_raw(&r_s.mpf, &s_s.mpf, &c_s.mpf);
  s21_mpf_set(res, &r_s.mpf);

  s21_mpf_stack_clear(&s_s);
  s21_mpf_stack_clear(&c_s);
  s21_mpf_stack_clear(&r_s);
  return 0;
}

int s21_mpf_asinh(s21_mpf_t *res, const s21_mpf_t *x) {
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
    s21_mpf_set_inf(res, x->sign);
    return 0;
  }

  /* Малые |x|: asinh(x) = x − x³/6 + …  → x с точностью выше prec. */
  if (x->exp < -(int64_t)res->prec / 2) {
    s21_mpf_set(res, x);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t xa_s, sq_s, one_s, sum_s, root_s;
  s21_mpf_stack_init(&xa_s,   wp);
  s21_mpf_stack_init(&sq_s,   wp);
  s21_mpf_stack_init(&one_s,  wp);
  s21_mpf_stack_init(&sum_s,  wp);
  s21_mpf_stack_init(&root_s, wp);

  s21_mpf_abs(&xa_s.mpf, x);
  s21_mpf_set_ui(&one_s.mpf, 1);
  s21_mpf_mul_raw(&sq_s.mpf, &xa_s.mpf, &xa_s.mpf);
  s21_mpf_add_raw(&sum_s.mpf, &sq_s.mpf, &one_s.mpf);
  s21_mpf_sqrt(&root_s.mpf, &sum_s.mpf);
  s21_mpf_add_raw(&root_s.mpf, &root_s.mpf, &xa_s.mpf);
  s21_mpf_log(&sum_s.mpf, &root_s.mpf);

  if (x->sign) sum_s.mpf.sign = 1;
  s21_mpf_set(res, &sum_s.mpf);

  s21_mpf_stack_clear(&xa_s);
  s21_mpf_stack_clear(&sq_s);
  s21_mpf_stack_clear(&one_s);
  s21_mpf_stack_clear(&sum_s);
  s21_mpf_stack_clear(&root_s);
  return 0;
}

int s21_mpf_acosh(s21_mpf_t *res, const s21_mpf_t *x) {
  if (res == NULL || x == NULL) return -1;
  if (x->kind == S21_MPF_NAN) {
    s21_mpf_set_nan(res);
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    s21_mpf_set_inf(res, 0);
    return 0;
  }
  if (x->kind == S21_MPF_ZERO || x->sign) {
    s21_mpf_set_nan(res);
    return 0;
  }

  /* Проверка x vs 1. */
  s21_mpf_t one;
  s21_mpf_init2(&one, x->prec);
  s21_mpf_set_ui(&one, 1);
  int cmp = s21_mpf_cmp(x, &one);
  s21_mpf_clear(&one);
  if (cmp == 0) {
    s21_mpf_set_zero(res, 0);
    return 0;
  }
  if (cmp < 0) {
    s21_mpf_set_nan(res);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t xw_s, one_s, xm_s, xp_s, prod_s, root_s, sum_s;
  s21_mpf_stack_init(&xw_s,   wp);
  s21_mpf_stack_init(&one_s,  wp);
  s21_mpf_stack_init(&xm_s,   wp);
  s21_mpf_stack_init(&xp_s,   wp);
  s21_mpf_stack_init(&prod_s, wp);
  s21_mpf_stack_init(&root_s, wp);
  s21_mpf_stack_init(&sum_s,  wp);

  s21_mpf_set(&xw_s.mpf, x); /* приведение x к рабочей точности wp */
  s21_mpf_set_ui(&one_s.mpf, 1);
  s21_mpf_sub_raw(&xm_s.mpf, &xw_s.mpf, &one_s.mpf); /* x − 1 */
  s21_mpf_add_raw(&xp_s.mpf, &xw_s.mpf, &one_s.mpf); /* x + 1 */
  s21_mpf_mul_raw(&prod_s.mpf, &xm_s.mpf, &xp_s.mpf);
  s21_mpf_sqrt(&root_s.mpf, &prod_s.mpf);
  s21_mpf_add_raw(&sum_s.mpf, &xw_s.mpf, &root_s.mpf);
  s21_mpf_log(res, &sum_s.mpf);

  s21_mpf_stack_clear(&xw_s);
  s21_mpf_stack_clear(&one_s);
  s21_mpf_stack_clear(&xm_s);
  s21_mpf_stack_clear(&xp_s);
  s21_mpf_stack_clear(&prod_s);
  s21_mpf_stack_clear(&root_s);
  s21_mpf_stack_clear(&sum_s);
  return 0;
}

int s21_mpf_atanh(s21_mpf_t *res, const s21_mpf_t *x) {
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
    s21_mpf_set_nan(res);
    return 0;
  }

  /* |x| = 1 → ±inf, |x| > 1 → NaN. */
  s21_mpf_t one;
  s21_mpf_init2(&one, x->prec);
  s21_mpf_set_ui(&one, 1);
  int cmp = s21_mpf_cmp_abs(x, &one);
  s21_mpf_clear(&one);
  if (cmp == 0) {
    s21_mpf_set_inf(res, x->sign);
    return 0;
  }
  if (cmp > 0) {
    s21_mpf_set_nan(res);
    return 0;
  }

  /* Малые |x|: atanh(x) = x + x³/3 + … → x с точностью выше prec. */
  if (x->exp < -(int64_t)res->prec / 2) {
    s21_mpf_set(res, x);
    return 0;
  }

  uint32_t wp = res->prec + 64;

  s21_mpf_stack_t xw_s, one_s, num_s, den_s, ratio_s;
  s21_mpf_stack_init(&xw_s,    wp);
  s21_mpf_stack_init(&one_s,   wp);
  s21_mpf_stack_init(&num_s,   wp);
  s21_mpf_stack_init(&den_s,   wp);
  s21_mpf_stack_init(&ratio_s, wp);

  s21_mpf_set(&xw_s.mpf, x); /* приведение x к рабочей точности wp */
  s21_mpf_set_ui(&one_s.mpf, 1);
  s21_mpf_add_raw(&num_s.mpf, &one_s.mpf, &xw_s.mpf); /* 1 + x */
  s21_mpf_sub_raw(&den_s.mpf, &one_s.mpf, &xw_s.mpf); /* 1 − x */
  s21_mpf_div_raw(&ratio_s.mpf, &num_s.mpf, &den_s.mpf);
  s21_mpf_log(&num_s.mpf, &ratio_s.mpf);
  num_s.mpf.exp -= 1; /* × ½ */
  s21_mpf_normalize(&num_s.mpf);
  s21_mpf_set(res, &num_s.mpf);

  s21_mpf_stack_clear(&xw_s);
  s21_mpf_stack_clear(&one_s);
  s21_mpf_stack_clear(&num_s);
  s21_mpf_stack_clear(&den_s);
  s21_mpf_stack_clear(&ratio_s);
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
  s21_mpf_compute_pi_agm(pi_local, wp);
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

/* Сколько «лишних» бит нужно добавить к рабочей точности из-за
   cancelation при reduce_2pi.  Для типичных |x| < 2 это 0;
   растёт как log2(|x|). */
static int64_t trig_extra_bits(const s21_mpf_t *x) {
  int64_t E = x->exp - 1;
  return E > 0 ? E : 0;
}

/* ============================================================
   sin — публичный + _impl
   ============================================================ */

void s21_mpf_sin_impl(s21_mpf_t *res, const s21_mpf_t *x,
                      const s21_mpf_t *pi_in) {
  int64_t E = trig_extra_bits(x);
  if (E > S21_MPF_MAX_ARG_BITS) {
    s21_mpf_set_nan(res);
    return;
  }
  uint32_t wp = res->prec + 64 + (uint32_t)E;

  if (pi_in != NULL && pi_in->prec < wp) {
    s21_mpf_set_nan(res);
    return;
  }

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
  int64_t E = trig_extra_bits(x);
  if (E > S21_MPF_MAX_ARG_BITS) {
    s21_mpf_set_nan(res);
    return;
  }
  uint32_t wp = res->prec + 64 + (uint32_t)E;

  if (pi_in != NULL && pi_in->prec < wp) {
    s21_mpf_set_nan(res);
    return;
  }

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

  int64_t E = trig_extra_bits(x);
  if (E > S21_MPF_MAX_ARG_BITS) {
    s21_mpf_set_nan(res);
    return 0;
  }
  uint32_t pi_prec = res->prec + 128 + (uint32_t)E;

  s21_mpf_t pi;
  s21_mpf_init2(&pi, pi_prec);
  s21_mpf_compute_pi_agm(&pi, pi_prec);
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
    s21_mpf_compute_pi_agm(&pi, wp);
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
    s21_mpf_compute_pi_agm(&pi, wp);
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
    s21_mpf_compute_pi_agm(&pi, wp);
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
    s21_mpf_compute_pi_agm(&pi, wp);

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
  if (x->kind != S21_MPF_NORMAL) return 0; /* спецзначения — раньше */
  int64_t E = trig_extra_bits(x);
  if (E > S21_MPF_MAX_ARG_BITS) return -1;
  uint64_t need = (uint64_t)res->prec + 128 + (uint64_t)E;
  if (pi->prec < need) return -1;
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
