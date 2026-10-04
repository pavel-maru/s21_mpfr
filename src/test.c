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
  int msb = s21_mpf_msb(&x);
  ASSERT(msb == 255, "42: старший бит должен быть на позиции 255");
  s21_mpf_print_d(&x);

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
  s21_mpf_print_d(&x);

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
  s21_mpf_print_d(&x);

  s21_mpf_set_d(&x, -0.5);
  s21_mpf_print_d(&x);

  s21_mpf_set_d(&x, 1e10);
  s21_mpf_print_d(&x);

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

static void test_cmp(void) {
  s21_mpf_t a, b;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&b, 256);

  s21_mpf_set_ui(&a, 42);
  s21_mpf_set_ui(&b, 42);
  ASSERT(s21_mpf_cmp(&a, &b) == 0, "42 == 42");
  ASSERT(s21_mpf_equal(&a, &b), "42 == 42 (equal)");

  s21_mpf_set_ui(&a, 42);
  s21_mpf_set_ui(&b, 43);
  ASSERT(s21_mpf_cmp(&a, &b) == -1, "42 < 43");
  ASSERT(s21_mpf_cmp(&b, &a) == 1, "43 > 42");

  s21_mpf_set_si(&a, -5);
  s21_mpf_set_ui(&b, 3);
  ASSERT(s21_mpf_cmp(&a, &b) == -1, "-5 < 3");

  s21_mpf_set_si(&a, -5);
  s21_mpf_set_si(&b, -3);
  ASSERT(s21_mpf_cmp(&a, &b) == -1, "-5 < -3");

  s21_mpf_set_zero(&a, 0);
  s21_mpf_set_ui(&b, 1);
  ASSERT(s21_mpf_cmp(&a, &b) == -1, "0 < 1");

  s21_mpf_set_ui(&a, 100);
  s21_mpf_set_ui(&b, 99);
  ASSERT(s21_mpf_cmp_abs(&a, &b) == 1, "|100| > |99|");

  s21_mpf_set_si(&a, -100);
  s21_mpf_set_ui(&b, 99);
  ASSERT(s21_mpf_cmp_abs(&a, &b) == 1, "|-100| > |99|");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_ui(&b, 999999);
  ASSERT(s21_mpf_cmp(&a, &b) == 1, "+inf > 999999");

  s21_mpf_set_inf(&a, 1);
  ASSERT(s21_mpf_cmp(&a, &b) == -1, "-inf < 999999");

  s21_mpf_set_ui(&a, 100);
  ASSERT(s21_mpf_integer_p(&a), "100 — целое");

  s21_mpf_set_ui(&a, 2);
  ASSERT(s21_mpf_integer_p(&a), "2 — целое");

  s21_mpf_set_d(&a, 0.5);
  ASSERT(!s21_mpf_integer_p(&a), "0.5 — не целое");

  s21_mpf_set_d(&a, 1.0);
  ASSERT(s21_mpf_integer_p(&a), "1.0 — целое");

  s21_mpf_set_d(&a, -3.0);
  ASSERT(s21_mpf_integer_p(&a), "-3.0 — целое");

  s21_mpf_set_d(&a, 1e10);
  ASSERT(s21_mpf_integer_p(&a), "1e10 — целое");

  s21_mpf_set_d(&a, 3.14);
  ASSERT(!s21_mpf_integer_p(&a), "3.14 — не целое");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  printf("[ok] cmp / cmp_abs / equal / integer_p\n");
}

static void test_add_sub(void) {
  s21_mpf_t a, b, c, expected;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&b, 256);
  s21_mpf_init2(&c, 256);
  s21_mpf_init2(&expected, 256);

  s21_mpf_set_ui(&a, 2);
  s21_mpf_set_ui(&b, 3);
  s21_mpf_add(&c, &a, &b);
  s21_mpf_set_ui(&expected, 5);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "2 + 3 = 5");

  s21_mpf_set_ui(&a, 10);
  s21_mpf_set_ui(&b, 4);
  s21_mpf_sub(&c, &a, &b);
  s21_mpf_set_ui(&expected, 6);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "10 - 4 = 6");

  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_si(&b, -3);
  s21_mpf_add(&c, &a, &b);
  s21_mpf_set_ui(&expected, 2);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "5 + (-3) = 2");

  s21_mpf_set_ui(&a, 3);
  s21_mpf_set_ui(&b, 5);
  s21_mpf_sub(&c, &a, &b);
  s21_mpf_set_si(&expected, -2);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "3 - 5 = -2");

  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_ui(&b, 5);
  s21_mpf_sub(&c, &a, &b);
  ASSERT(s21_mpf_is_zero(&c), "5 - 5 = 0");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_set_ui(&b, 42);
  s21_mpf_add(&c, &a, &b);
  s21_mpf_set_ui(&expected, 42);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "0 + 42 = 42");

  s21_mpf_set_ui(&a, 1000000);
  s21_mpf_set_ui(&b, 999999);
  s21_mpf_add(&c, &a, &b);
  s21_mpf_set_ui(&expected, 1999999);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "1000000 + 999999 = 1999999");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_ui(&b, 1);
  s21_mpf_add(&c, &a, &b);
  ASSERT(s21_mpf_is_inf(&c) && c.sign == 0, "inf + 1 = inf");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_inf(&b, 1);
  s21_mpf_add(&c, &a, &b);
  ASSERT(s21_mpf_is_nan(&c), "inf + (-inf) = NaN");

  s21_mpf_set_si(&a, -7);
  s21_mpf_neg(&b, &a);
  s21_mpf_set_ui(&expected, 7);
  ASSERT(s21_mpf_cmp(&b, &expected) == 0, "-(-7) = 7");

  s21_mpf_set_si(&a, -7);
  s21_mpf_abs(&b, &a);
  ASSERT(s21_mpf_cmp(&b, &expected) == 0, "|-7| = 7");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&expected);
  printf("[ok] add / sub / neg / abs\n");
}

static void test_mul(void) {
  s21_mpf_t a, b, c, expected;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&b, 256);
  s21_mpf_init2(&c, 256);
  s21_mpf_init2(&expected, 256);

  s21_mpf_set_ui(&a, 2);
  s21_mpf_set_ui(&b, 3);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 6);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "2 * 3 = 6");

  s21_mpf_set_ui(&a, 10);
  s21_mpf_set_ui(&b, 10);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 100);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "10 * 10 = 100");

  s21_mpf_set_si(&a, -3);
  s21_mpf_set_ui(&b, 4);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_si(&expected, -12);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "-3 * 4 = -12");

  s21_mpf_set_si(&a, -3);
  s21_mpf_set_si(&b, -4);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 12);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "-3 * (-4) = 12");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_set_ui(&b, 5);
  s21_mpf_mul(&c, &a, &b);
  ASSERT(s21_mpf_is_zero(&c), "0 * 5 = 0");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_set_ui(&b, 1);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 1);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "1 * 1 = 1");

  s21_mpf_set_ui(&a, 12345);
  s21_mpf_set_ui(&b, 6789);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 83810205);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "12345 * 6789 = 83810205");

  s21_mpf_set_ui(&a, 123456789);
  s21_mpf_set_ui(&b, 987654321);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 121932631112635269UL);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0,
         "123456789 * 987654321 = 121932631112635269");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_ui(&b, 5);
  s21_mpf_mul(&c, &a, &b);
  ASSERT(s21_mpf_is_inf(&c) && c.sign == 0, "inf * 5 = inf");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_ui(&b, 0);
  s21_mpf_mul(&c, &a, &b);
  ASSERT(s21_mpf_is_nan(&c), "inf * 0 = NaN");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_si(&b, -1);
  s21_mpf_mul(&c, &a, &b);
  ASSERT(s21_mpf_is_inf(&c) && c.sign == 1, "inf * (-1) = -inf");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&expected);
  printf("[ok] mul\n");
}

int main(void) {
  printf("=== s21_mpf: базовые тесты ===\n");
  test_init_clear();
  test_set_ui();
  test_set_si();
  test_set_d();
  test_special();
  test_precision();
  test_cmp();
  test_add_sub();
  test_mul();
  printf("=== Все тесты прошли ===\n");
  return 0;
}
