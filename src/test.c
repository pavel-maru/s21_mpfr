#include "s21_mpf.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ASSERT(cond, msg) do { \
  if (!(cond)) { \
    fprintf(stderr, "FAIL: %s at %s:%d\n", msg, __FILE__, __LINE__); \
    exit(1); \
  } \
} while (0)

static void test_init_clear(void) {
  s21_mpf_t x;
  s21_mpf_init2(&x, 128);
  ASSERT(x.prec == 128, "prec не установлен");
  ASSERT(x.kind == S21_MPF_ZERO, "не ZERO по умолчанию");
  s21_mpf_clear(&x);
  printf("[ok] init/clear\n");
}

static void test_set_ui(void) {
  s21_mpf_t x;
  s21_mpf_init2(&x, 256);

  s21_mpf_set_ui(&x, 42);
  ASSERT(x.kind == S21_MPF_NORMAL, "42 должен быть NORMAL");
  ASSERT(x.sign == 0, "42 положительный");
  /* После нормализации старший бит на позиции prec-1 = 255 */
  int msb = s21_mpf_msb(&x);
  ASSERT(msb == 255, "42: старший бит должен быть на позиции 255");
  /* Проверка значения через print_d */
  s21_mpf_print_d(&x);  /* ожидаем 42 */

  s21_mpf_set_ui(&x, 0);
  ASSERT(x.kind == S21_MPF_ZERO, "0 должен быть ZERO");

  s21_mpf_set_ui(&x, 1);
  ASSERT(x.kind == S21_MPF_NORMAL, "1 должен быть NORMAL");

  s21_mpf_clear(&x);
  printf("[ok] set_ui\n");
}

static void test_set_si(void) {
  s21_mpf_t x;
  s21_mpf_init2(&x, 256);

  s21_mpf_set_si(&x, -7);
  ASSERT(x.sign == 1, "-7 отрицательный");
  ASSERT(x.kind == S21_MPF_NORMAL, "-7 NORMAL");
  s21_mpf_print_d(&x);  /* ожидаем -7 */

  s21_mpf_set_si(&x, 0);
  ASSERT(x.kind == S21_MPF_ZERO, "0 (si) ZERO");

  s21_mpf_clear(&x);
  printf("[ok] set_si\n");
}

static void test_set_d(void) {
  s21_mpf_t x;
  s21_mpf_init2(&x, 256);

  s21_mpf_set_d(&x, 3.14);
  ASSERT(x.kind == S21_MPF_NORMAL, "3.14 NORMAL");
  s21_mpf_print_d(&x);  /* ожидаем 3.14... */

  s21_mpf_set_d(&x, -0.5);
  s21_mpf_print_d(&x);  /* ожидаем -0.5 */

  s21_mpf_set_d(&x, 1e10);
  s21_mpf_print_d(&x);  /* ожидаем 1e+10 */

  s21_mpf_clear(&x);
  printf("[ok] set_d\n");
}

static void test_special(void) {
  s21_mpf_t x;
  s21_mpf_init(&x);

  s21_mpf_set_nan(&x);
  ASSERT(s21_mpf_is_nan(&x), "NaN");

  s21_mpf_set_inf(&x, 0);
  ASSERT(s21_mpf_is_inf(&x) && s21_mpf_sign(&x) == 1, "+inf");

  s21_mpf_set_inf(&x, 1);
  ASSERT(s21_mpf_is_inf(&x) && s21_mpf_sign(&x) == -1, "-inf");

  s21_mpf_set_zero(&x, 1);
  ASSERT(s21_mpf_is_zero(&x) && s21_mpf_sign(&x) == 0, "-0");

  s21_mpf_clear(&x);
  printf("[ok] special values\n");
}

static void test_precision(void) {
  s21_mpf_t x;
  for (uint32_t prec = 2; prec <= 512; prec += 50) {
    s21_mpf_init2(&x, prec);
    s21_mpf_set_ui(&x, 12345);
    int msb = s21_mpf_msb(&x);
    ASSERT(msb == (int)prec - 1, "старший бит не на prec-1");
    s21_mpf_clear(&x);
  }
  printf("[ok] precisions 2..512\n");
}

int main(void) {
  printf("=== s21_mpf: базовые тесты ===\n");
  test_init_clear();
  test_set_ui();
  test_set_si();
  test_set_d();
  test_special();
  test_precision();
  printf("=== Все тесты прошли ===\n");
  return 0;
}
