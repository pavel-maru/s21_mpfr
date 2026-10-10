/* ============================================================
   s21_mpfr — микробенчмарк публичных операций.

   Для каждой операции и каждой точности из фиксированного
   набора запускаем серию вызовов, замеряем суммарное время
   через CLOCK_MONOTONIC, берём лучший из TRIALS прогонов
   (защита от шума планировщика), печатаем ns/op.

   Число итераций подобрано так, чтобы каждая серия длилась
   ~10-200 мс: мелкие операции на малых точностях гоняем
   сотнями тысяч раз, тяжёлые трансцендентные на 4096 битах —
   единицами.

   Вызовы идут через публичный API — то, что видит пользователь.
   Время включает аллокации временных объектов внутри операций.
   ============================================================ */

#define _POSIX_C_SOURCE 200809L

#include <locale.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include "s21_mpf.h"

typedef int (*binary_fn_t)(s21_mpf_t *, const s21_mpf_t *, const s21_mpf_t *);
typedef int (*unary_fn_t)(s21_mpf_t *, const s21_mpf_t *);

#define TRIALS 3

/* Аккумулятор: компилятор не может выбросить вычисления, потому
   что их результат кто-то читает. Само значение неважно — важно,
   что чтение происходит. */
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
   Бинарные: add / sub / mul / div.
   ------------------------------------------------------------ */

static void bench_binary(const char *name, binary_fn_t fn, uint32_t prec,
                         int iters) {
  s21_mpf_t a, b, r;
  s21_mpf_init2(&a, prec);
  s21_mpf_init2(&b, prec);
  s21_mpf_init2(&r, prec);
  s21_mpf_set_d(&a, 1.5);
  s21_mpf_set_d(&b, 2.5);

  for (int i = 0; i < 5; i++) fn(&r, &a, &b); /* warmup */

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
   Унарные: sqrt / exp / log / sin.
   ------------------------------------------------------------ */

static void bench_unary(const char *name, unary_fn_t fn, uint32_t prec,
                        int iters) {
  s21_mpf_t a, r;
  s21_mpf_init2(&a, prec);
  s21_mpf_init2(&r, prec);
  s21_mpf_set_d(&a, 1.5);

  for (int i = 0; i < 5; i++) fn(&r, &a); /* warmup */

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
   Число итераций по (prec, operation). Подобрано эмпирически:
   add занимает ~10x меньше, чем mul, поэтому iters разные.
   Ориентир — каждая серия от 10 мс до 200 мс.
   ------------------------------------------------------------ */

typedef struct {
  uint32_t prec;
  int addsub;   /* add, sub */
  int mul;
  int div;
  int sqrt;
  int transc;   /* exp, log, sin */
} bench_config_t;

static const bench_config_t configs[] = {
    {  64, 200000, 100000, 50000,  5000, 1000},
    { 256,  50000,  20000,  5000,  1000,  300},
    {1024,   5000,   2000,   500,   200,   30},
    {4096,   2000,    500,    50,    30,    5},
};

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
    bench_unary ("sin",  s21_mpf_sin,  c->prec, c->transc);
    printf("\n");
  }

  printf("(checksum %llu)\n", (unsigned long long)g_sink);
  return 0;
}
