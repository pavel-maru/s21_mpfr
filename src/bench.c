/* ============================================================
   s21_mpfr — микробенчмарк публичных операций.

   Для каждой операции и каждой точности из фиксированного
   набора запускаем серию вызовов, замеряем суммарное время
   через CLOCK_MONOTONIC, берём лучший из TRIALS прогонов
   (защита от шума планировщика), печатаем ns/op.

   Отдельно показаны sin/cos/tan с явно переданной π — они
   не пересчитывают π на каждый вызов (формула Мэчина,
   ~1200 членов при prec = 4096), поэтому разница заметна
   на порядок.
   ============================================================ */

#define _POSIX_C_SOURCE 200809L

#include <locale.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "s21_mpf.h"

typedef int (*binary_fn_t)(s21_mpf_t *, const s21_mpf_t *, const s21_mpf_t *);
typedef int (*unary_fn_t)(s21_mpf_t *, const s21_mpf_t *);
typedef int (*unary_pi_fn_t)(s21_mpf_t *, const s21_mpf_t *,
                             const s21_mpf_t *);

#define TRIALS 3

static volatile uint64_t g_sink = 0;

static void consume(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_NORMAL && x->limbs != NULL) {
    g_sink ^= x->limbs[0];
  }
}

static double now_ns(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}

/* ------------------------------------------------------------
   Бинарные.
   ------------------------------------------------------------ */

static void bench_binary(const char *name, binary_fn_t fn, uint32_t prec,
                         int iters) {
  s21_mpf_t a, b, r;
  s21_mpf_init2(&a, prec);
  s21_mpf_init2(&b, prec);
  s21_mpf_init2(&r, prec);
  s21_mpf_set_d(&a, 1.5);
  s21_mpf_set_d(&b, 2.5);

  for (int i = 0; i < 5; i++) fn(&r, &a, &b);

  double best = 1e30;
  for (int t = 0; t < TRIALS; t++) {
    double t0 = now_ns();
    for (int i = 0; i < iters; i++) fn(&r, &a, &b);
    double dt = (now_ns() - t0) / (double)iters;
    if (dt < best) best = dt;
  }
  consume(&r);
  printf("  %-5s %6u  %12.3f   %8d iter\n", name, prec, best, iters);

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&r);
}

/* ------------------------------------------------------------
   Унарные (без π).
   ------------------------------------------------------------ */

static void bench_unary(const char *name, unary_fn_t fn, uint32_t prec,
                        int iters) {
  s21_mpf_t a, r;
  s21_mpf_init2(&a, prec);
  s21_mpf_init2(&r, prec);
  s21_mpf_set_d(&a, 1.5);

  for (int i = 0; i < 5; i++) fn(&r, &a);

  double best = 1e30;
  for (int t = 0; t < TRIALS; t++) {
    double t0 = now_ns();
    for (int i = 0; i < iters; i++) fn(&r, &a);
    double dt = (now_ns() - t0) / (double)iters;
    if (dt < best) best = dt;
  }
  consume(&r);
  printf("  %-5s %6u  %12.3f   %8d iter\n", name, prec, best, iters);

  s21_mpf_clear(&a);
  s21_mpf_clear(&r);
}

/* ------------------------------------------------------------
   Унарные с явной π.
   ------------------------------------------------------------ */

static void bench_unary_with_pi(const char *name, unary_pi_fn_t fn,
                                uint32_t prec, int iters) {
  s21_mpf_t a, r, pi;
  s21_mpf_init2(&a, prec);
  s21_mpf_init2(&r, prec);
  s21_mpf_init2(&pi, prec + 128);
  s21_mpf_set_d(&a, 1.5);
  s21_mpf_pi(&pi);

  for (int i = 0; i < 5; i++) fn(&r, &a, &pi);

  double best = 1e30;
  for (int t = 0; t < TRIALS; t++) {
    double t0 = now_ns();
    for (int i = 0; i < iters; i++) fn(&r, &a, &pi);
    double dt = (now_ns() - t0) / (double)iters;
    if (dt < best) best = dt;
  }
  consume(&r);
  printf("  %-5s %6u  %12.3f   %8d iter  (with pi)\n", name, prec, best, iters);

  s21_mpf_clear(&a);
  s21_mpf_clear(&r);
  s21_mpf_clear(&pi);
}

/* ------------------------------------------------------------
   Число итераций по (prec, operation).
   ------------------------------------------------------------ */

typedef struct {
  uint32_t prec;
  int addsub;
  int mul;
  int div;
  int sqrt;
  int transc;
} bench_config_t;

static const bench_config_t configs[] = {
    {  64, 200000, 100000, 50000,  5000, 1000},
    { 256,  50000,  20000,  5000,  1000,  300},
    {1024,   5000,   2000,   500,   200,   30},
    {4096,   2000,    500,    50,    30,    5},
    /* Два дополнительных масштаба — посмотреть поведение mul
       выше порога Карацубы и оценить crossover до Toom-Cook. */
    {8192,   1000,    200,    20,    10,    2},
    {16384,   500,     60,     6,     3,    1},
};

/* Обёртка для s21_mpf_pi — публичная функция принимает только
   res, а bench_unary работает с unary_fn_t (res, x). */
static int bench_pi_wrapper(s21_mpf_t *res, const s21_mpf_t *x) {
  (void)x;
  return s21_mpf_pi(res);
}

int main(void) {
  setlocale(LC_NUMERIC, "C");

  printf("s21_mpfr benchmark — best of %d trials, CLOCK_MONOTONIC\n", TRIALS);
  printf("  %-5s %6s  %12s   %8s\n", "op", "prec", "ns/op", "iters");
  printf("  ------------------------------------------------\n");

  const size_t nc = sizeof(configs) / sizeof(configs[0]);
  for (size_t i = 0; i < nc; i++) {
    const bench_config_t *c = &configs[i];
    printf("  --- prec = %u ---\n", c->prec);
    bench_binary("add",  s21_mpf_add,  c->prec, c->addsub);
    bench_binary("sub",  s21_mpf_sub,  c->prec, c->addsub);
    bench_binary("mul",  s21_mpf_mul,  c->prec, c->mul);
    bench_binary("div",  s21_mpf_div,  c->prec, c->div);
    bench_unary ("sqrt", s21_mpf_sqrt, c->prec, c->sqrt);
    bench_unary ("exp",  s21_mpf_exp,  c->prec, c->transc);
    bench_unary ("log",  s21_mpf_log,  c->prec, c->transc);
    printf("\n");
    bench_unary ("pi",   bench_pi_wrapper, c->prec, c->transc);
    bench_unary_with_pi("sin", s21_mpf_sin_with_pi, c->prec, c->transc);
    bench_unary_with_pi("cos", s21_mpf_cos_with_pi, c->prec, c->transc);
    bench_unary_with_pi("tan", s21_mpf_tan_with_pi, c->prec, c->transc);
    printf("\n");
    bench_unary ("sin", s21_mpf_sin, c->prec, c->transc);
    bench_unary ("cos", s21_mpf_cos, c->prec, c->transc);
    bench_unary ("tan", s21_mpf_tan, c->prec, c->transc);
    printf("\n");
  }

  printf("(checksum %llu)\n", (unsigned long long)g_sink);
  return 0;
}
