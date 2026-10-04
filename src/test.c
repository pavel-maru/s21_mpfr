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

static void test_div(void) {
  s21_mpf_t a, b, c, expected;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&b, 256);
  s21_mpf_init2(&c, 256);
  s21_mpf_init2(&expected, 256);

  s21_mpf_set_ui(&a, 6);
  s21_mpf_set_ui(&b, 3);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_ui(&expected, 2);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "6 / 3 = 2");

  s21_mpf_set_ui(&a, 100);
  s21_mpf_set_ui(&b, 25);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_ui(&expected, 4);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "100 / 25 = 4");

  s21_mpf_set_ui(&a, 1000000000);
  s21_mpf_set_ui(&b, 1000);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_ui(&expected, 1000000);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "1e9 / 1e3 = 1e6");

  s21_mpf_set_si(&a, -12);
  s21_mpf_set_ui(&b, 4);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_si(&expected, -3);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "-12 / 4 = -3");

  s21_mpf_set_si(&a, -12);
  s21_mpf_set_si(&b, -4);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_ui(&expected, 3);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "-12 / -4 = 3");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_set_ui(&b, 1);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_ui(&expected, 1);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "1 / 1 = 1");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_set_ui(&b, 5);
  s21_mpf_div(&c, &a, &b);
  ASSERT(s21_mpf_is_zero(&c), "0 / 5 = 0");

  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_inf(&b, 0);
  s21_mpf_div(&c, &a, &b);
  ASSERT(s21_mpf_is_zero(&c), "5 / inf = 0");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_ui(&b, 5);
  s21_mpf_div(&c, &a, &b);
  ASSERT(s21_mpf_is_inf(&c) && c.sign == 0, "inf / 5 = inf");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_set_inf(&b, 0);
  s21_mpf_div(&c, &a, &b);
  ASSERT(s21_mpf_is_nan(&c), "inf / inf = NaN");

  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_ui(&b, 0);
  s21_mpf_div(&c, &a, &b);
  ASSERT(s21_mpf_is_inf(&c) && c.sign == 0, "5 / 0 = inf");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_set_ui(&b, 0);
  s21_mpf_div(&c, &a, &b);
  ASSERT(s21_mpf_is_nan(&c), "0 / 0 = NaN");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_set_ui(&b, 3);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_t back;
  s21_mpf_init2(&back, 256);
  s21_mpf_mul(&back, &c, &b);
  s21_mpf_set_ui(&expected, 1);
  s21_mpf_t diff;
  s21_mpf_init2(&diff, 256);
  s21_mpf_sub(&diff, &back, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_t threshold;
  s21_mpf_init2(&threshold, 256);
  s21_mpf_set_d(&threshold, 1e-60);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "(1/3)*3 ≈ 1");
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  s21_mpf_clear(&back);

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&expected);
  printf("[ok] div\n");
}

static void test_null_checks(void) {
  s21_mpf_t a, b, c;
  s21_mpf_init2(&a, 64);
  s21_mpf_init2(&b, 64);
  s21_mpf_init2(&c, 64);

  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_ui(&b, 3);

  /* add */
  ASSERT(s21_mpf_add(NULL, &a, &b) == -1, "add: NULL res");
  ASSERT(s21_mpf_add(&c, NULL, &b) == -1, "add: NULL x");
  ASSERT(s21_mpf_add(&c, &a, NULL) == -1, "add: NULL y");

  /* sub */
  ASSERT(s21_mpf_sub(NULL, &a, &b) == -1, "sub: NULL res");
  ASSERT(s21_mpf_sub(&c, NULL, &b) == -1, "sub: NULL x");
  ASSERT(s21_mpf_sub(&c, &a, NULL) == -1, "sub: NULL y");

  /* mul */
  ASSERT(s21_mpf_mul(NULL, &a, &b) == -1, "mul: NULL res");
  ASSERT(s21_mpf_mul(&c, NULL, &b) == -1, "mul: NULL x");
  ASSERT(s21_mpf_mul(&c, &a, NULL) == -1, "mul: NULL y");

  /* div */
  ASSERT(s21_mpf_div(NULL, &a, &b) == -1, "div: NULL res");
  ASSERT(s21_mpf_div(&c, NULL, &b) == -1, "div: NULL x");
  ASSERT(s21_mpf_div(&c, &a, NULL) == -1, "div: NULL y");

  /* neg / abs — void, но не должны падать */
  s21_mpf_neg(NULL, &a);
  s21_mpf_neg(&c, NULL);
  s21_mpf_abs(NULL, &a);
  s21_mpf_abs(&c, NULL);

  /* Разные точности — тоже -1 */
  s21_mpf_t d;
  s21_mpf_init2(&d, 128);
  s21_mpf_set_ui(&d, 7);
  ASSERT(s21_mpf_add(&c, &a, &d) == -1, "add: prec mismatch");
  ASSERT(s21_mpf_mul(&c, &a, &d) == -1, "mul: prec mismatch");
  s21_mpf_clear(&d);

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  printf("[ok] null checks\n");
}

static void test_bit_utils(void) {
  s21_mpf_t x;
  s21_mpf_init2(&x, 128);
  s21_mpf_set_ui(&x, 1);

  /* set_bit / get_bit round-trip */
  for (uint32_t pos = 0; pos < 128; pos += 17) {
    s21_mpf_set_bit(&x, pos, 1);
    ASSERT(s21_mpf_get_bit(&x, pos) == 1, "set_bit 1 -> get_bit 1");
    s21_mpf_set_bit(&x, pos, 0);
    ASSERT(s21_mpf_get_bit(&x, pos) == 0, "set_bit 0 -> get_bit 0");
  }

  /* Бит за пределами prec игнорируется */
  s21_mpf_set_bit(&x, 200, 1);
  ASSERT(s21_mpf_get_bit(&x, 200) == 0, "bit above prec ignored");

  s21_mpf_clear(&x);

  /* shift_left_into */
  {
    uint64_t src[2] = {0x1ULL, 0};
    uint64_t dst[2] = {0};
    s21_mpf_shift_left_into(dst, src, 2, 1);
    ASSERT(dst[0] == 0x2ULL, "shift left by 1: low");
    ASSERT(dst[1] == 0, "shift left by 1: high");

    uint64_t src2[2] = {0x8000000000000000ULL, 0};
    uint64_t dst2[2] = {0};
    s21_mpf_shift_left_into(dst2, src2, 2, 1);
    ASSERT(dst2[0] == 0, "shift through boundary: low zero");
    ASSERT(dst2[1] == 0x1ULL, "shift through boundary: high 1");

    uint64_t src3[2] = {0xDEADBEEFULL, 0};
    uint64_t dst3[2] = {0};
    s21_mpf_shift_left_into(dst3, src3, 2, 64);
    ASSERT(dst3[0] == 0, "shift by 64: low zero");
    ASSERT(dst3[1] == 0xDEADBEEFULL, "shift by 64: high moved");

    uint64_t src4[2] = {0xCAFEULL, 0xF00DULL};
    uint64_t dst4[2] = {0};
    s21_mpf_shift_left_into(dst4, src4, 2, 0);
    ASSERT(dst4[0] == 0xCAFEULL && dst4[1] == 0xF00DULL, "shift by 0 copy");
  }

  printf("[ok] bit utilities\n");
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
  test_div();
  test_null_checks();
  test_bit_utils();
  printf("=== Все тесты прошли ===\n");
  return 0;
}
