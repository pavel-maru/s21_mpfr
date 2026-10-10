#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

#define ASSERT(cond, msg)                                              \
  do {                                                                 \
    if (!(cond)) {                                                     \
      fprintf(stderr, "FAIL: %s at %s:%d\n", msg, __FILE__, __LINE__); \
      exit(1);                                                         \
    }                                                                  \
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

static void test_small_prec(void) {
  s21_mpf_t x, y, z, expected;

  /* prec = 2: представимы 0, 1, 2, 3 */
  s21_mpf_init2(&x, 2);
  s21_mpf_init2(&y, 2);
  s21_mpf_init2(&z, 2);
  s21_mpf_init2(&expected, 2);

  s21_mpf_set_ui(&x, 1);
  ASSERT(x.kind == S21_MPF_NORMAL, "1 при prec=2 — NORMAL");
  ASSERT(s21_mpf_msb(&x) == 1, "1 при prec=2 — MSB на 1");

  s21_mpf_set_ui(&x, 3);
  ASSERT(s21_mpf_msb(&x) == 1, "3 при prec=2 — MSB на 1");

  s21_mpf_set_ui(&x, 1);
  s21_mpf_set_ui(&y, 2);
  s21_mpf_add(&z, &x, &y);
  s21_mpf_set_ui(&expected, 3);
  ASSERT(s21_mpf_cmp(&z, &expected) == 0, "1 + 2 = 3 при prec=2");

  s21_mpf_set_ui(&x, 3);
  s21_mpf_set_ui(&y, 1);
  s21_mpf_sub(&z, &x, &y);
  s21_mpf_set_ui(&expected, 2);
  ASSERT(s21_mpf_cmp(&z, &expected) == 0, "3 - 1 = 2 при prec=2");

  s21_mpf_clear(&x);
  s21_mpf_clear(&y);
  s21_mpf_clear(&z);
  s21_mpf_clear(&expected);

  /* prec = 3 */
  s21_mpf_init2(&x, 3);
  s21_mpf_init2(&y, 3);
  s21_mpf_init2(&z, 3);
  s21_mpf_init2(&expected, 3);

  s21_mpf_set_ui(&x, 5);
  s21_mpf_set_ui(&y, 2);
  s21_mpf_add(&z, &x, &y);
  s21_mpf_set_ui(&expected, 7);
  ASSERT(s21_mpf_cmp(&z, &expected) == 0, "5 + 2 = 7 при prec=3");

  s21_mpf_clear(&x);
  s21_mpf_clear(&y);
  s21_mpf_clear(&z);
  s21_mpf_clear(&expected);

  /* prec = 5 */
  s21_mpf_init2(&x, 5);
  s21_mpf_init2(&y, 5);
  s21_mpf_init2(&z, 5);
  s21_mpf_init2(&expected, 5);

  s21_mpf_set_ui(&x, 15);
  s21_mpf_set_ui(&y, 16);
  s21_mpf_add(&z, &x, &y);
  s21_mpf_set_ui(&expected, 31);
  ASSERT(s21_mpf_cmp(&z, &expected) == 0, "15 + 16 = 31 при prec=5");

  s21_mpf_clear(&x);
  s21_mpf_clear(&y);
  s21_mpf_clear(&z);
  s21_mpf_clear(&expected);

  /* prec = 8 */
  s21_mpf_init2(&x, 8);
  s21_mpf_init2(&y, 8);
  s21_mpf_init2(&z, 8);
  s21_mpf_init2(&expected, 8);

  s21_mpf_set_ui(&x, 100);
  s21_mpf_set_ui(&y, 27);
  s21_mpf_sub(&z, &x, &y);
  s21_mpf_set_ui(&expected, 73);
  ASSERT(s21_mpf_cmp(&z, &expected) == 0, "100 - 27 = 73 при prec=8");

  s21_mpf_clear(&x);
  s21_mpf_clear(&y);
  s21_mpf_clear(&z);
  s21_mpf_clear(&expected);

  printf("[ok] малые precisions (2, 3, 5, 8)\n");
}

/* Регрессия: при сужении 66 → 2 бита (word_shift = 1, dst_count = 1,
   src_count = 2) shift_right_into раньше сравнивал src_idx с dst_count
   вместо src_count и терял старший лимб. Именно поэтому падал
   "1 + 2 = 3 при prec=2". */
static void test_round_narrow_across_limb_boundary(void) {
  s21_mpf_t a, b, c, e;
  s21_mpf_init2(&a, 2);
  s21_mpf_init2(&b, 2);
  s21_mpf_init2(&c, 2);
  s21_mpf_init2(&e, 2);

  s21_mpf_set_ui(&a, 1);
  s21_mpf_set_ui(&b, 2);
  s21_mpf_add(&c, &a, &b);
  s21_mpf_set_ui(&e, 3);
  ASSERT(s21_mpf_cmp(&c, &e) == 0,
         "1 + 2 = 3 при prec=2 (regression, limb boundary)");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&e);
  printf("[ok] round narrow across limb boundary\n");
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

static void test_cmp_mixed_prec(void) {
  s21_mpf_t a256, b64, eps;
  s21_mpf_init2(&a256, 256);
  s21_mpf_init2(&b64, 64);
  s21_mpf_init2(&eps, 256);

  /* Равные значения в разных точностях */
  s21_mpf_set_ui(&a256, 1);
  s21_mpf_set_ui(&b64, 1);
  ASSERT(s21_mpf_cmp_abs(&a256, &b64) == 0, "1 (256) == 1 (64)");
  ASSERT(s21_mpf_cmp(&a256, &b64) == 0, "1 (256) == 1 (64)");

  s21_mpf_set_d(&a256, 0.5);
  s21_mpf_set_d(&b64, 0.5);
  ASSERT(s21_mpf_cmp_abs(&a256, &b64) == 0, "0.5 (256) == 0.5 (64)");

  /* a = 1 + 2^-100 в prec=256, b = 1 в prec=64.
     Разница 2^-100 << 2^-52, double её не видит. */
  s21_mpf_set_ui(&a256, 1);
  s21_mpf_set_ui(&eps, 1);
  eps.exp -= 100;
  s21_mpf_normalize(&eps);
  s21_mpf_add(&a256, &a256, &eps);

  s21_mpf_set_ui(&b64, 1);

  ASSERT(s21_mpf_cmp_abs(&a256, &b64) > 0,
         "1 + 2^-100 (256) > 1 (64) — double-приближение не видит разницы");
  ASSERT(s21_mpf_cmp_abs(&b64, &a256) < 0, "1 (64) < 1 + 2^-100 (256)");

  s21_mpf_clear(&a256);
  s21_mpf_clear(&b64);
  s21_mpf_clear(&eps);
  printf("[ok] cmp / cmp_abs (mixed prec)\n");
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

static void test_nonmult32_prec(void) {
  s21_mpf_t a, b, c, expected;

  /* prec = 33 — не кратно 32 */
  s21_mpf_init2(&a, 33);
  s21_mpf_init2(&b, 33);
  s21_mpf_init2(&c, 33);
  s21_mpf_init2(&expected, 33);

  s21_mpf_set_ui(&a, 3);
  s21_mpf_set_ui(&b, 5);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 15);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "3 * 5 = 15 при prec=33");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&expected);

  /* prec = 100 */
  s21_mpf_init2(&a, 100);
  s21_mpf_init2(&b, 100);
  s21_mpf_init2(&c, 100);
  s21_mpf_init2(&expected, 100);

  s21_mpf_set_ui(&a, 1234);
  s21_mpf_set_ui(&b, 5678);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 7006652);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "1234 * 5678 = 7006652 при prec=100");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&expected);

  /* prec = 130 */
  s21_mpf_init2(&a, 130);
  s21_mpf_init2(&b, 130);
  s21_mpf_init2(&c, 130);
  s21_mpf_init2(&expected, 130);

  s21_mpf_set_ui(&a, 1000000);
  s21_mpf_set_ui(&b, 1000000);
  s21_mpf_mul(&c, &a, &b);
  s21_mpf_set_ui(&expected, 1000000000000UL);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "1e6 * 1e6 = 1e12 при prec=130");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&expected);

  printf("[ok] prec, не кратные 32 (33, 100, 130)\n");
}

static void test_mul_many_prec(void) {
  const uint32_t precs[] = {4,   5,   7,   32,  33,  63,  64, 65,
                            100, 127, 128, 129, 200, 256, 500};
  for (size_t i = 0; i < sizeof(precs) / sizeof(precs[0]); i++) {
    uint32_t p = precs[i];
    s21_mpf_t a, b, c, e;
    s21_mpf_init2(&a, p);
    s21_mpf_init2(&b, p);
    s21_mpf_init2(&c, p);
    s21_mpf_init2(&e, p);

    s21_mpf_set_ui(&a, 3);
    s21_mpf_set_ui(&b, 5);
    s21_mpf_mul(&c, &a, &b);
    s21_mpf_set_ui(&e, 15);
    ASSERT(s21_mpf_cmp(&c, &e) == 0, "3 * 5 = 15 при любом p >= 4");

    s21_mpf_clear(&a);
    s21_mpf_clear(&b);
    s21_mpf_clear(&c);
    s21_mpf_clear(&e);
  }
  printf("[ok] mul at many precisions\n");
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

static void test_div_many_prec(void) {
  const uint32_t precs[] = {4,   5,   7,   32,  33,  63,  64, 65,
                            100, 127, 128, 129, 200, 256, 500};
  for (size_t i = 0; i < sizeof(precs) / sizeof(precs[0]); i++) {
    uint32_t p = precs[i];
    s21_mpf_t a, b, c, e;
    s21_mpf_init2(&a, p);
    s21_mpf_init2(&b, p);
    s21_mpf_init2(&c, p);
    s21_mpf_init2(&e, p);

    s21_mpf_set_ui(&a, 100);
    s21_mpf_set_ui(&b, 4);
    s21_mpf_div(&c, &a, &b);
    s21_mpf_set_ui(&e, 25);
    ASSERT(s21_mpf_cmp(&c, &e) == 0, "100 / 4 = 25 при любом p");

    s21_mpf_set_ui(&a, 1000);
    s21_mpf_set_ui(&b, 8);
    s21_mpf_div(&c, &a, &b);
    s21_mpf_set_ui(&e, 125);
    ASSERT(s21_mpf_cmp(&c, &e) == 0, "1000 / 8 = 125 при любом p");

    s21_mpf_clear(&a);
    s21_mpf_clear(&b);
    s21_mpf_clear(&c);
    s21_mpf_clear(&e);
  }
  printf("[ok] div at many precisions\n");
}

/* Точные степени двойки и round-trip div/mul. */
static void test_div_powers_of_two(void) {
  s21_mpf_t a, b, c, e;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&b, 256);
  s21_mpf_init2(&c, 256);
  s21_mpf_init2(&e, 256);

  /* 2^100 / 2^40 = 2^60 */
  s21_mpf_set_ui(&a, 1);
  a.exp = 101;
  s21_mpf_normalize(&a);
  s21_mpf_set_ui(&b, 1);
  b.exp = 41;
  s21_mpf_normalize(&b);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_ui(&e, 1);
  e.exp = 61;
  s21_mpf_normalize(&e);
  ASSERT(s21_mpf_cmp(&c, &e) == 0, "2^100 / 2^40 = 2^60");

  /* 2^200 / 2^199 = 2 (минимальный шаг экспоненты) */
  s21_mpf_set_ui(&a, 1);
  a.exp = 201;
  s21_mpf_normalize(&a);
  s21_mpf_set_ui(&b, 1);
  b.exp = 200;
  s21_mpf_normalize(&b);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_set_ui(&e, 2);
  ASSERT(s21_mpf_cmp(&c, &e) == 0, "2^200 / 2^199 = 2");

  /* Round-trip: (7/3)*3 ≈ 7.  Проверяет согласованность div/mul. */
  s21_mpf_set_ui(&a, 7);
  s21_mpf_set_ui(&b, 3);
  s21_mpf_div(&c, &a, &b);
  s21_mpf_mul(&c, &c, &b);
  s21_mpf_sub(&c, &c, &a);
  s21_mpf_abs(&c, &c);
  s21_mpf_t threshold;
  s21_mpf_init2(&threshold, 256);
  s21_mpf_set_d(&threshold, 1e-70);
  ASSERT(s21_mpf_cmp(&c, &threshold) < 0, "(7/3)*3 ≈ 7");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&c);
  s21_mpf_clear(&e);
  s21_mpf_clear(&threshold);
  printf("[ok] div: powers of two and round-trip\n");
}

static void test_null_checks(void) {
  s21_mpf_t a, b, c;
  s21_mpf_init2(&a, 64);
  s21_mpf_init2(&b, 64);
  s21_mpf_init2(&c, 64);

  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_ui(&b, 3);

  ASSERT(s21_mpf_add(NULL, &a, &b) == -1, "add: NULL res");
  ASSERT(s21_mpf_add(&c, NULL, &b) == -1, "add: NULL x");
  ASSERT(s21_mpf_add(&c, &a, NULL) == -1, "add: NULL y");

  ASSERT(s21_mpf_sub(NULL, &a, &b) == -1, "sub: NULL res");
  ASSERT(s21_mpf_sub(&c, NULL, &b) == -1, "sub: NULL x");
  ASSERT(s21_mpf_sub(&c, &a, NULL) == -1, "sub: NULL y");

  ASSERT(s21_mpf_mul(NULL, &a, &b) == -1, "mul: NULL res");
  ASSERT(s21_mpf_mul(&c, NULL, &b) == -1, "mul: NULL x");
  ASSERT(s21_mpf_mul(&c, &a, NULL) == -1, "mul: NULL y");

  ASSERT(s21_mpf_div(NULL, &a, &b) == -1, "div: NULL res");
  ASSERT(s21_mpf_div(&c, NULL, &b) == -1, "div: NULL x");
  ASSERT(s21_mpf_div(&c, &a, NULL) == -1, "div: NULL y");

  s21_mpf_neg(NULL, &a);
  s21_mpf_neg(&c, NULL);
  s21_mpf_abs(NULL, &a);
  s21_mpf_abs(&c, NULL);

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

  for (uint32_t pos = 0; pos < 128; pos += 17) {
    s21_mpf_set_bit(&x, pos, 1);
    ASSERT(s21_mpf_get_bit(&x, pos) == 1, "set_bit 1 -> get_bit 1");
    s21_mpf_set_bit(&x, pos, 0);
    ASSERT(s21_mpf_get_bit(&x, pos) == 0, "set_bit 0 -> get_bit 0");
  }

  s21_mpf_set_bit(&x, 200, 1);
  ASSERT(s21_mpf_get_bit(&x, 200) == 0, "bit above prec ignored");

  s21_mpf_clear(&x);

  {
    uint64_t src[2] = {0x1ULL, 0};
    uint64_t dst[2] = {0};
    s21_mpf_shift_left_into(dst, 2, src, 2, 1);
    ASSERT(dst[0] == 0x2ULL, "shift left by 1: low");
    ASSERT(dst[1] == 0, "shift left by 1: high");

    uint64_t src2[2] = {0x8000000000000000ULL, 0};
    uint64_t dst2[2] = {0};
    s21_mpf_shift_left_into(dst2, 2, src2, 2, 1);
    ASSERT(dst2[0] == 0, "shift through boundary: low zero");
    ASSERT(dst2[1] == 0x1ULL, "shift through boundary: high 1");

    uint64_t src3[2] = {0xDEADBEEFULL, 0};
    uint64_t dst3[2] = {0};
    s21_mpf_shift_left_into(dst3, 2, src3, 2, 64);
    ASSERT(dst3[0] == 0, "shift by 64: low zero");
    ASSERT(dst3[1] == 0xDEADBEEFULL, "shift by 64: high moved");

    uint64_t src4[2] = {0xCAFEULL, 0xF00DULL};
    uint64_t dst4[2] = {0};
    s21_mpf_shift_left_into(dst4, 2, src4, 2, 0);
    ASSERT(dst4[0] == 0xCAFEULL && dst4[1] == 0xF00DULL, "shift by 0 copy");

    /* Регрессия на разные размеры dst/src:
       src_count = 2, dst_count = 1, сдвиг = 64 — старший лимб должен
       переехать в младший, не потеряться. */
    uint64_t src5[2] = {0, 0x5ULL};
    uint64_t dst5[1] = {0};
    s21_mpf_shift_right_into(dst5, 1, src5, 2, 64);
    ASSERT(dst5[0] == 0x5ULL, "right shift 64 across limb boundary");
  }

  printf("[ok] bit utilities\n");
}

static void test_sqrt(void) {
  s21_mpf_t a, s, expected;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&s, 256);
  s21_mpf_init2(&expected, 256);

  s21_mpf_set_ui(&a, 4);
  s21_mpf_sqrt(&s, &a);
  s21_mpf_set_ui(&expected, 2);
  ASSERT(s21_mpf_cmp(&s, &expected) == 0, "sqrt(4) = 2");

  s21_mpf_set_ui(&a, 9);
  s21_mpf_sqrt(&s, &a);
  s21_mpf_set_ui(&expected, 3);
  ASSERT(s21_mpf_cmp(&s, &expected) == 0, "sqrt(9) = 3");

  s21_mpf_set_ui(&a, 16);
  s21_mpf_sqrt(&s, &a);
  s21_mpf_set_ui(&expected, 4);
  ASSERT(s21_mpf_cmp(&s, &expected) == 0, "sqrt(16) = 4");

  s21_mpf_set_ui(&a, 1000000);
  s21_mpf_sqrt(&s, &a);
  s21_mpf_set_ui(&expected, 1000);
  ASSERT(s21_mpf_cmp(&s, &expected) == 0, "sqrt(1e6) = 1000");

  s21_mpf_set_ui(&a, 144);
  s21_mpf_sqrt(&s, &a);
  s21_mpf_set_ui(&expected, 12);
  ASSERT(s21_mpf_cmp(&s, &expected) == 0, "sqrt(144) = 12");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_sqrt(&s, &a);
  s21_mpf_set_ui(&expected, 1);
  ASSERT(s21_mpf_cmp(&s, &expected) == 0, "sqrt(1) = 1");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_sqrt(&s, &a);
  ASSERT(s21_mpf_is_zero(&s), "sqrt(0) = 0");

  s21_mpf_set_si(&a, -1);
  s21_mpf_sqrt(&s, &a);
  ASSERT(s21_mpf_is_nan(&s), "sqrt(-1) = NaN");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_sqrt(&s, &a);
  ASSERT(s21_mpf_is_inf(&s) && s.sign == 0, "sqrt(+inf) = +inf");

  s21_mpf_set_nan(&a);
  s21_mpf_sqrt(&s, &a);
  ASSERT(s21_mpf_is_nan(&s), "sqrt(NaN) = NaN");

  s21_mpf_set_ui(&a, 2);
  s21_mpf_sqrt(&s, &a);
  s21_mpf_t sq;
  s21_mpf_init2(&sq, 256);
  s21_mpf_mul(&sq, &s, &s);
  s21_mpf_set_ui(&expected, 2);
  s21_mpf_t diff;
  s21_mpf_init2(&diff, 256);
  s21_mpf_sub(&diff, &sq, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_t threshold;
  s21_mpf_init2(&threshold, 256);
  s21_mpf_set_d(&threshold, 1e-60);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "sqrt(2)^2 ≈ 2");
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  s21_mpf_clear(&sq);

  ASSERT(s21_mpf_sqrt(NULL, &a) == -1, "sqrt: NULL res");
  ASSERT(s21_mpf_sqrt(&s, NULL) == -1, "sqrt: NULL x");

  s21_mpf_clear(&a);
  s21_mpf_clear(&s);
  s21_mpf_clear(&expected);
  printf("[ok] sqrt\n");
}

static void test_exp(void) {
  s21_mpf_t a, e, expected, diff, threshold;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&e, 256);
  s21_mpf_init2(&expected, 256);
  s21_mpf_init2(&diff, 256);
  s21_mpf_init2(&threshold, 256);

  s21_mpf_set_ui(&a, 0);
  s21_mpf_exp(&e, &a);
  s21_mpf_set_ui(&expected, 1);
  ASSERT(s21_mpf_cmp(&e, &expected) == 0, "exp(0) = 1");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_exp(&e, &a);
  s21_mpf_set_d(&expected, 2.718281828459045);
  s21_mpf_sub(&diff, &e, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-14);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "exp(1) ≈ e");

  s21_mpf_set_si(&a, -1);
  s21_mpf_exp(&e, &a);
  s21_mpf_set_d(&expected, 0.367879441171442);
  s21_mpf_sub(&diff, &e, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "exp(-1) ≈ 1/e");

  s21_mpf_set_ui(&a, 2);
  s21_mpf_exp(&e, &a);
  s21_mpf_set_d(&expected, 7.389056098930650);
  s21_mpf_sub(&diff, &e, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-13);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "exp(2) ≈ 7.389");

  s21_mpf_set_ui(&a, 10);
  s21_mpf_exp(&e, &a);
  s21_mpf_set_d(&expected, 22026.465794806718);
  s21_mpf_sub(&diff, &e, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-10);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "exp(10) ≈ 22026.5");

  s21_mpf_set_nan(&a);
  s21_mpf_exp(&e, &a);
  ASSERT(s21_mpf_is_nan(&e), "exp(NaN) = NaN");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_exp(&e, &a);
  ASSERT(s21_mpf_is_inf(&e) && e.sign == 0, "exp(+inf) = +inf");

  s21_mpf_set_inf(&a, 1);
  s21_mpf_exp(&e, &a);
  ASSERT(s21_mpf_is_zero(&e), "exp(-inf) = 0");

  ASSERT(s21_mpf_exp(NULL, &a) == -1, "exp: NULL res");
  ASSERT(s21_mpf_exp(&e, NULL) == -1, "exp: NULL x");

  s21_mpf_clear(&a);
  s21_mpf_clear(&e);
  s21_mpf_clear(&expected);
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  printf("[ok] exp\n");
}

static void test_log(void) {
  s21_mpf_t a, l, expected, diff, threshold;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&l, 256);
  s21_mpf_init2(&expected, 256);
  s21_mpf_init2(&diff, 256);
  s21_mpf_init2(&threshold, 256);

  s21_mpf_set_ui(&a, 1);
  s21_mpf_log(&l, &a);
  ASSERT(s21_mpf_is_zero(&l), "log(1) = 0");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_exp(&a, &a);
  s21_mpf_log(&l, &a);
  s21_mpf_set_ui(&expected, 1);
  s21_mpf_sub(&diff, &l, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-14);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "log(e) ≈ 1");

  s21_mpf_set_ui(&a, 2);
  s21_mpf_log(&l, &a);
  s21_mpf_set_d(&expected, 0.693147180559945);
  s21_mpf_sub(&diff, &l, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "log(2) ≈ ln2");

  s21_mpf_set_ui(&a, 10);
  s21_mpf_log(&l, &a);
  s21_mpf_set_d(&expected, 2.302585092994046);
  s21_mpf_sub(&diff, &l, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-13);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "log(10) ≈ 2.3026");

  s21_mpf_set_d(&a, 0.5);
  s21_mpf_log(&l, &a);
  s21_mpf_set_d(&expected, -0.693147180559945);
  s21_mpf_sub(&diff, &l, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-14);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "log(0.5) ≈ -ln2");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_log(&l, &a);
  ASSERT(s21_mpf_is_inf(&l) && l.sign == 1, "log(0) = -inf");

  s21_mpf_set_si(&a, -1);
  s21_mpf_log(&l, &a);
  ASSERT(s21_mpf_is_nan(&l), "log(-1) = NaN");

  s21_mpf_set_nan(&a);
  s21_mpf_log(&l, &a);
  ASSERT(s21_mpf_is_nan(&l), "log(NaN) = NaN");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_log(&l, &a);
  ASSERT(s21_mpf_is_inf(&l) && l.sign == 0, "log(+inf) = +inf");

  ASSERT(s21_mpf_log(NULL, &a) == -1, "log: NULL res");
  ASSERT(s21_mpf_log(&l, NULL) == -1, "log: NULL x");

  s21_mpf_clear(&a);
  s21_mpf_clear(&l);
  s21_mpf_clear(&expected);
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  printf("[ok] log\n");
}

static void test_pi(void) {
  s21_mpf_t pi, expected, diff, threshold;
  s21_mpf_init2(&pi, 256);
  s21_mpf_init2(&expected, 256);
  s21_mpf_init2(&diff, 256);
  s21_mpf_init2(&threshold, 256);

  s21_mpf_pi(&pi);
  s21_mpf_set_d(&expected, 3.141592653589793);
  s21_mpf_sub(&diff, &pi, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-15);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "π ≈ 3.141592653589793");

  ASSERT(s21_mpf_pi(NULL) == -1, "pi: NULL");

  s21_mpf_clear(&pi);
  s21_mpf_clear(&expected);
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  printf("[ok] pi\n");
}

static void test_sin_cos(void) {
  s21_mpf_t a, s, c, expected, diff, threshold, pi;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&s, 256);
  s21_mpf_init2(&c, 256);
  s21_mpf_init2(&expected, 256);
  s21_mpf_init2(&diff, 256);
  s21_mpf_init2(&threshold, 256);
  s21_mpf_init2(&pi, 256);

  s21_mpf_set_d(&threshold, 1e-14);

  s21_mpf_set_ui(&a, 0);
  s21_mpf_sin(&s, &a);
  ASSERT(s21_mpf_is_zero(&s), "sin(0) = 0");

  s21_mpf_cos(&c, &a);
  s21_mpf_set_ui(&expected, 1);
  ASSERT(s21_mpf_cmp(&c, &expected) == 0, "cos(0) = 1");

  s21_mpf_set_d(&a, 1.0);
  s21_mpf_sin(&s, &a);
  s21_mpf_set_d(&expected, 0.8414709848078965);
  s21_mpf_sub(&diff, &s, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "sin(1) ≈ 0.84147");

  s21_mpf_cos(&c, &a);
  s21_mpf_set_d(&expected, 0.5403023058681398);
  s21_mpf_sub(&diff, &c, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "cos(1) ≈ 0.54030");

  s21_mpf_pi(&pi);
  s21_mpf_set(&a, &pi);
  a.exp -= 1;
  s21_mpf_normalize(&a);
  s21_mpf_sin(&s, &a);
  s21_mpf_set_ui(&expected, 1);
  s21_mpf_sub(&diff, &s, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "sin(π/2) ≈ 1");

  s21_mpf_cos(&c, &pi);
  s21_mpf_set_si(&expected, -1);
  s21_mpf_sub(&diff, &c, &expected);
  s21_mpf_abs(&diff, &diff);
  s21_mpf_set_d(&threshold, 1e-13);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "cos(π) ≈ -1");

  s21_mpf_set_d(&a, 2.5);
  s21_mpf_sin(&s, &a);
  s21_mpf_neg(&a, &a);
  s21_mpf_sin(&expected, &a);
  s21_mpf_neg(&expected, &expected);
  ASSERT(s21_mpf_cmp(&s, &expected) == 0, "sin(-x) = -sin(x)");

  s21_mpf_set_nan(&a);
  s21_mpf_sin(&s, &a);
  ASSERT(s21_mpf_is_nan(&s), "sin(NaN) = NaN");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_sin(&s, &a);
  ASSERT(s21_mpf_is_nan(&s), "sin(+inf) = NaN");

  ASSERT(s21_mpf_sin(NULL, &a) == -1, "sin: NULL res");
  ASSERT(s21_mpf_sin(&s, NULL) == -1, "sin: NULL x");

  s21_mpf_clear(&a);
  s21_mpf_clear(&s);
  s21_mpf_clear(&c);
  s21_mpf_clear(&expected);
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  s21_mpf_clear(&pi);
  printf("[ok] sin / cos\n");
}

static void test_atan(void) {
  s21_mpf_t a, at, expected, diff, threshold;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&at, 256);
  s21_mpf_init2(&expected, 256);
  s21_mpf_init2(&diff, 256);
  s21_mpf_init2(&threshold, 256);

  s21_mpf_set_d(&threshold, 1e-14);

  s21_mpf_set_ui(&a, 0);
  s21_mpf_atan(&at, &a);
  ASSERT(s21_mpf_is_zero(&at), "atan(0) = 0");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_atan(&at, &a);
  s21_mpf_set_d(&expected, 0.7853981633974483);
  s21_mpf_sub(&diff, &at, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "atan(1) ≈ π/4");

  s21_mpf_set_d(&a, 0.5);
  s21_mpf_atan(&at, &a);
  s21_mpf_set_d(&expected, 0.4636476090008061);
  s21_mpf_sub(&diff, &at, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "atan(0.5) ≈ 0.46365");

  s21_mpf_set_ui(&a, 10);
  s21_mpf_atan(&at, &a);
  s21_mpf_set_d(&expected, 1.4711276743037347);
  s21_mpf_sub(&diff, &at, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "atan(10) ≈ 1.4711");

  s21_mpf_set_si(&a, -1);
  s21_mpf_atan(&at, &a);
  s21_mpf_set_d(&expected, -0.7853981633974483);
  s21_mpf_sub(&diff, &at, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "atan(-1) ≈ -π/4");

  s21_mpf_set_inf(&a, 0);
  s21_mpf_atan(&at, &a);
  s21_mpf_set_d(&expected, 1.5707963267948966);
  s21_mpf_sub(&diff, &at, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "atan(+inf) ≈ π/2");

  s21_mpf_set_inf(&a, 1);
  s21_mpf_atan(&at, &a);
  s21_mpf_set_d(&expected, -1.5707963267948966);
  s21_mpf_sub(&diff, &at, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "atan(-inf) ≈ -π/2");

  ASSERT(s21_mpf_atan(NULL, &a) == -1, "atan: NULL res");
  ASSERT(s21_mpf_atan(&at, NULL) == -1, "atan: NULL x");

  s21_mpf_clear(&a);
  s21_mpf_clear(&at);
  s21_mpf_clear(&expected);
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  printf("[ok] atan\n");
}

static void test_tan_asin_acos(void) {
  s21_mpf_t a, r, expected, diff, threshold;
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&r, 256);
  s21_mpf_init2(&expected, 256);
  s21_mpf_init2(&diff, 256);
  s21_mpf_init2(&threshold, 256);

  s21_mpf_set_d(&threshold, 1e-14);

  s21_mpf_set_ui(&a, 0);
  s21_mpf_tan(&r, &a);
  ASSERT(s21_mpf_is_zero(&r), "tan(0) = 0");

  s21_mpf_set_d(&a, 1.0);
  s21_mpf_tan(&r, &a);
  s21_mpf_set_d(&expected, 1.5574077246549023);
  s21_mpf_sub(&diff, &r, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "tan(1) ≈ 1.5574");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_asin(&r, &a);
  ASSERT(s21_mpf_is_zero(&r), "asin(0) = 0");

  s21_mpf_set_d(&a, 0.5);
  s21_mpf_asin(&r, &a);
  s21_mpf_set_d(&expected, 0.5235987755982989);
  s21_mpf_sub(&diff, &r, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "asin(0.5) ≈ π/6");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_asin(&r, &a);
  s21_mpf_set_d(&expected, 1.5707963267948966);
  s21_mpf_sub(&diff, &r, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "asin(1) ≈ π/2");

  s21_mpf_set_ui(&a, 2);
  s21_mpf_asin(&r, &a);
  ASSERT(s21_mpf_is_nan(&r), "asin(2) = NaN");

  s21_mpf_set_ui(&a, 0);
  s21_mpf_acos(&r, &a);
  s21_mpf_set_d(&expected, 1.5707963267948966);
  s21_mpf_sub(&diff, &r, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "acos(0) ≈ π/2");

  s21_mpf_set_ui(&a, 1);
  s21_mpf_acos(&r, &a);
  ASSERT(s21_mpf_is_zero(&r), "acos(1) = 0");

  s21_mpf_set_si(&a, -1);
  s21_mpf_acos(&r, &a);
  s21_mpf_set_d(&expected, 3.141592653589793);
  s21_mpf_sub(&diff, &r, &expected);
  s21_mpf_abs(&diff, &diff);
  ASSERT(s21_mpf_cmp(&diff, &threshold) < 0, "acos(-1) ≈ π");

  s21_mpf_set_ui(&a, 2);
  s21_mpf_acos(&r, &a);
  ASSERT(s21_mpf_is_nan(&r), "acos(2) = NaN");

  ASSERT(s21_mpf_tan(NULL, &a) == -1, "tan: NULL res");
  ASSERT(s21_mpf_tan(&r, NULL) == -1, "tan: NULL x");
  ASSERT(s21_mpf_asin(NULL, &a) == -1, "asin: NULL res");
  ASSERT(s21_mpf_asin(&r, NULL) == -1, "asin: NULL x");
  ASSERT(s21_mpf_acos(NULL, &a) == -1, "acos: NULL res");
  ASSERT(s21_mpf_acos(&r, NULL) == -1, "acos: NULL x");

  s21_mpf_clear(&a);
  s21_mpf_clear(&r);
  s21_mpf_clear(&expected);
  s21_mpf_clear(&diff);
  s21_mpf_clear(&threshold);
  printf("[ok] tan / asin / acos\n");
}

static void test_set_round(void) {
  /* x = 1 + 3·2⁻⁶⁴ при prec=256, приведённое к prec=64.
     Это ровно середина между 1+2⁻⁶³ и 1+2⁻⁶²:
       RNDN → 1+2⁻⁶² (ties-to-even)
       RNDZ → 1+2⁻⁶³
       RNDU → 1+2⁻⁶²
       RNDD → 1+2⁻⁶³ */
  s21_mpf_t hi, rndn, rndz, rndu, rndd;
  s21_mpf_init2(&hi, 256);
  s21_mpf_init2(&rndn, 64);
  s21_mpf_init2(&rndz, 64);
  s21_mpf_init2(&rndu, 64);
  s21_mpf_init2(&rndd, 64);

  s21_mpf_set_ui(&hi, 1);
  s21_mpf_t tiny;
  s21_mpf_init2(&tiny, 256);
  s21_mpf_set_ui(&tiny, 3);
  tiny.exp -= 64;
  s21_mpf_normalize(&tiny);
  s21_mpf_add(&hi, &hi, &tiny);
  s21_mpf_clear(&tiny);

  s21_mpf_set_round(&rndn, &hi, S21_MPF_RNDN);
  s21_mpf_set_round(&rndz, &hi, S21_MPF_RNDZ);
  s21_mpf_set_round(&rndu, &hi, S21_MPF_RNDU);
  s21_mpf_set_round(&rndd, &hi, S21_MPF_RNDD);

  s21_mpf_t lo62;
  s21_mpf_init2(&lo62, 64);
  s21_mpf_set_ui(&lo62, 1);
  s21_mpf_t t62;
  s21_mpf_init2(&t62, 64);
  s21_mpf_set_ui(&t62, 1);
  t62.exp -= 62;
  s21_mpf_normalize(&t62);
  s21_mpf_add(&lo62, &lo62, &t62);
  s21_mpf_clear(&t62);

  s21_mpf_t lo63;
  s21_mpf_init2(&lo63, 64);
  s21_mpf_set_ui(&lo63, 1);
  s21_mpf_t t63;
  s21_mpf_init2(&t63, 64);
  s21_mpf_set_ui(&t63, 1);
  t63.exp -= 63;
  s21_mpf_normalize(&t63);
  s21_mpf_add(&lo63, &lo63, &t63);
  s21_mpf_clear(&t63);

  ASSERT(s21_mpf_cmp(&rndn, &lo62) == 0,
         "RNDN(1 + 3·2⁻⁶⁴) = 1 + 2⁻⁶² (ties-to-even)");
  ASSERT(s21_mpf_cmp(&rndz, &lo63) == 0,
         "RNDZ(1 + 3·2⁻⁶⁴) = 1 + 2⁻⁶³ (truncate)");
  ASSERT(s21_mpf_cmp(&rndu, &lo62) == 0, "RNDU(1 + 3·2⁻⁶⁴) = 1 + 2⁻⁶² (up)");
  ASSERT(s21_mpf_cmp(&rndd, &lo63) == 0, "RNDD(1 + 3·2⁻⁶⁴) = 1 + 2⁻⁶³ (down)");

  /* Расширение точности — всегда точное. */
  s21_mpf_t widened;
  s21_mpf_init2(&widened, 128);
  s21_mpf_set_round(&widened, &lo62, S21_MPF_RNDZ);
  ASSERT(s21_mpf_cmp(&widened, &lo62) == 0, "widening is exact");

  ASSERT(s21_mpf_set_round(NULL, &hi, S21_MPF_RNDN) == -1, "NULL dst");
  ASSERT(s21_mpf_set_round(&rndn, NULL, S21_MPF_RNDN) == -1, "NULL src");

  s21_mpf_clear(&hi);
  s21_mpf_clear(&rndn);
  s21_mpf_clear(&rndz);
  s21_mpf_clear(&rndu);
  s21_mpf_clear(&rndd);
  s21_mpf_clear(&lo62);
  s21_mpf_clear(&lo63);
  s21_mpf_clear(&widened);
  printf("[ok] set_round (RNDN/RNDZ/RNDU/RNDD)\n");
}

/* Добавить 2^-e к x (положительное e). */
static void add_pow2(s21_mpf_t *x, int e) {
  s21_mpf_t t;
  s21_mpf_init2(&t, x->prec);
  s21_mpf_set_ui(&t, 1);
  t.exp -= e;
  s21_mpf_normalize(&t);
  s21_mpf_add(x, x, &t);
  s21_mpf_clear(&t);
}

static void test_set_round_nonmult64(void) {
  s21_mpf_t a, b, e;

  /* Сценарий 1: widening 33 -> 130, значение сохраняется. */
  s21_mpf_init2(&a, 33);
  s21_mpf_init2(&b, 130);
  s21_mpf_init2(&e, 130);
  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_ui(&e, 5);
  ASSERT(s21_mpf_set_round(&b, &a, S21_MPF_RNDN) == 0, "wide: rc");
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "wide 33->130: 5 preserved");
  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&e);

  /* Сценарий 2: narrow exact 130 -> 33. */
  s21_mpf_init2(&a, 130);
  s21_mpf_init2(&b, 33);
  s21_mpf_init2(&e, 33);
  s21_mpf_set_ui(&a, 5);
  s21_mpf_set_ui(&e, 5);
  ASSERT(s21_mpf_set_round(&b, &a, S21_MPF_RNDN) == 0, "narrow exact: rc");
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "narrow 130->33: 5 preserved");
  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&e);

  /* Сценарий 3: tie-to-even при некратном целевом prec.
     x = 1 + 2^-63 при prec=256, сужение до 63.
     round_bit = 1, sticky = 0, lsb(mant_after) = 0 → округление к 1. */
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&b, 63);
  s21_mpf_init2(&e, 63);
  s21_mpf_set_ui(&a, 1);
  add_pow2(&a, 63); /* a = 1 + 2^-63 */
  s21_mpf_set_ui(&e, 1);
  ASSERT(s21_mpf_set_round(&b, &a, S21_MPF_RNDN) == 0, "tie: rc");
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "tie 256->63: RNDN -> 1");
  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&e);

  /* Сценарий 4: round + sticky, все четыре режима.
     x = 1 + 3·2^-64 при prec=256, сужение до 63.
     Мантисса источника: биты 192 (2^-63) и 191 (2^-64).
     round_bit = бит 192 = 1, sticky = бит 191 = 1.
     Ожидания:
       RNDN → 1 + 2^-62
       RNDZ → 1
       RNDU → 1 + 2^-62
       RNDD → 1 */
  s21_mpf_init2(&a, 256);
  s21_mpf_init2(&b, 63);
  s21_mpf_init2(&e, 63);
  s21_mpf_set_ui(&a, 1);
  add_pow2(&a, 63);
  add_pow2(&a, 64); /* a = 1 + 3·2^-64 */

  s21_mpf_set_round(&b, &a, S21_MPF_RNDN);
  s21_mpf_set_ui(&e, 1);
  add_pow2(&e, 62); /* e = 1 + 2^-62 */
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "RNDN: 1 + 3·2^-64 -> 1 + 2^-62");

  s21_mpf_set_round(&b, &a, S21_MPF_RNDZ);
  s21_mpf_set_ui(&e, 1);
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "RNDZ: 1 + 3·2^-64 -> 1");

  s21_mpf_set_round(&b, &a, S21_MPF_RNDU);
  s21_mpf_set_ui(&e, 1);
  add_pow2(&e, 62);
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "RNDU: 1 + 3·2^-64 -> 1 + 2^-62");

  s21_mpf_set_round(&b, &a, S21_MPF_RNDD);
  s21_mpf_set_ui(&e, 1);
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "RNDD: 1 + 3·2^-64 -> 1");

  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&e);

  /* Сценарий 5: round-up через границу мантиссы (случай (б) в mpf_add_one).
     x = 2^64 - 1 при prec=64, сужение до 63.
     mant_after = 2^63 - 1, round_bit = 1, lsb = 1 → +1 → 2^63.
     Без нормализации получилось бы MSB на позиции 63 вместо 62.
     Ожидание: 2^64 (mant = 2^62, exp = 65). */
  s21_mpf_init2(&a, 64);
  s21_mpf_init2(&b, 63);
  s21_mpf_init2(&e, 63);
  s21_mpf_set_ui(&a, 0xFFFFFFFFFFFFFFFFUL); /* 2^64 - 1 */
  ASSERT(s21_mpf_set_round(&b, &a, S21_MPF_RNDN) == 0, "overflow: rc");
  s21_mpf_set_ui(&e, 1);
  e.exp += 64; /* e = 2^64 при prec=63: mant=2^62, exp=65 */
  ASSERT(s21_mpf_cmp(&b, &e) == 0, "overflow 64->63: 2^64-1 -> 2^64");
  s21_mpf_clear(&a);
  s21_mpf_clear(&b);
  s21_mpf_clear(&e);

  printf("[ok] set_round: prec not multiple of 64\n");
}

/* ============================================================
   ULP-метрика для трансцендентных функций.

   Считаем f(x) при prec = N и prec = 2N. Округляем 2N-результат
   к N битам, вычитаем, получаем |diff|. Нормируем на ulp(N) —
   получаем gap в ULP младшей точности. Здоровый алгоритм даёт
   gap 0 или 1; gap > 2 — сигнал о систематической потере
   точности (недостаточная рабочая точность, обрезание ряда,
   накопление ошибок).

   Идея самосогласованная, без внешнего эталона: если f(x)
   при N и 2N согласованы, значит вычисления сходятся.
   ============================================================ */

/* Приближение mpf через double. Точность ~53 бита — для
   диагностики достаточно, поскольку diff обычно много меньше
   единицы ulp и представим в double без потерь. */
static double mpf_to_double(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_ZERO) return 0.0;
  if (x->kind != S21_MPF_NORMAL) return 0.0; /* не ожидается */

  size_t count = s21_mpf_limbs_for_prec(x->prec);
  double val = 0.0;
  int taken = 0;
  int i = (int)count - 1;
  for (; i >= 0 && taken < 2; i--, taken++) {
    val = val * 18446744073709551616.0 + (double)x->limbs[i];
  }
  int64_t missed = (int64_t)(i + 1);
  int64_t e = x->exp - (int64_t)x->prec + 64 * missed;
  double result = ldexp(val, (int)e);
  return x->sign ? -result : result;
}

typedef int (*unary_fn_t)(s21_mpf_t *, const s21_mpf_t *);

/* Обёртка для s21_mpf_pi — она не принимает аргумента. */
static int wrap_pi(s21_mpf_t *res, const s21_mpf_t *x) {
  (void)x;
  return s21_mpf_pi(res);
}

/* Возвращает gap = |r_lo − round_to_N(r_hi)| / ulp(r_lo),
   где r_lo = f(arg) при prec N, r_hi = f(arg) при prec 2N. */
static double measure_ulp_gap(unary_fn_t fn, double arg, uint32_t N) {
  s21_mpf_t x_lo, x_hi, r_lo, r_hi, r_hi_N, diff;
  s21_mpf_init2(&x_lo, N);
  s21_mpf_init2(&x_hi, 2 * N);
  s21_mpf_init2(&r_lo, N);
  s21_mpf_init2(&r_hi, 2 * N);
  s21_mpf_init2(&r_hi_N, N);
  s21_mpf_init2(&diff, N);

  s21_mpf_set_d(&x_lo, arg);
  s21_mpf_set_d(&x_hi, arg);
  fn(&r_lo, &x_lo);
  fn(&r_hi, &x_hi);
  s21_mpf_set_round(&r_hi_N, &r_hi, S21_MPF_RNDN);

  s21_mpf_sub(&diff, &r_lo, &r_hi_N);
  s21_mpf_abs(&diff, &diff);

  double gap = 0.0;
  if (!s21_mpf_is_zero(&diff)) {
    double d = mpf_to_double(&diff);
    gap = ldexp(d, (int)N - (int)r_lo.exp);
  }

  s21_mpf_clear(&x_lo);
  s21_mpf_clear(&x_hi);
  s21_mpf_clear(&r_lo);
  s21_mpf_clear(&r_hi);
  s21_mpf_clear(&r_hi_N);
  s21_mpf_clear(&diff);
  return gap;
}

static void test_ulp_diagnostics(void) {
  struct {
    const char *name;
    unary_fn_t fn;
    double arg;
  } cases[] = {
      {"sqrt(2)", s21_mpf_sqrt, 2.0},
      {"sqrt(10)", s21_mpf_sqrt, 10.0},
      {"exp(1)", s21_mpf_exp, 1.0},
      {"exp(-2)", s21_mpf_exp, -2.0},
      {"log(2)", s21_mpf_log, 2.0},
      {"log(10)", s21_mpf_log, 10.0},
      {"sin(1)", s21_mpf_sin, 1.0},
      {"cos(1)", s21_mpf_cos, 1.0},
      {"atan(1)", s21_mpf_atan, 1.0},
      {"atan(10)", s21_mpf_atan, 10.0},
      {"asin(0.5)", s21_mpf_asin, 0.5},
      {"acos(0.5)", s21_mpf_acos, 0.5},
      {"pi", wrap_pi, 0.0},
  };
  const uint32_t Ns[] = {32, 64, 128, 256};
  const size_t nc = sizeof(cases) / sizeof(cases[0]);
  const size_t nn = sizeof(Ns) / sizeof(Ns[0]);

  printf("  %-11s", "function");
  for (size_t j = 0; j < nn; j++) printf("  N=%-3u", Ns[j]);
  printf("\n");

  for (size_t i = 0; i < nc; i++) {
    printf("  %-11s", cases[i].name);
    for (size_t j = 0; j < nn; j++) {
      double gap = measure_ulp_gap(cases[i].fn, cases[i].arg, Ns[j]);
      printf("  %5.1f", gap);
      ASSERT(gap <= 2.0, "ulp gap > 2");
    }
    printf("\n");
  }
  printf("[ok] ulp diagnostics\n");
}

static void test_set_prec(void) {
  s21_mpf_t x, expected;

  /* 1. Расширение сохраняет значение. */
  s21_mpf_init2(&x, 32);
  s21_mpf_init2(&expected, 128);
  s21_mpf_set_ui(&x, 42);
  s21_mpf_set_ui(&expected, 42);
  ASSERT(s21_mpf_set_prec(&x, 128) == 0, "expand: rc");
  ASSERT(x.prec == 128, "expand: prec");
  ASSERT(s21_mpf_cmp(&x, &expected) == 0, "expand: value preserved");
  s21_mpf_clear(&x);
  s21_mpf_clear(&expected);

  /* 2. Сужение точно (значение представимо). */
  s21_mpf_init2(&x, 128);
  s21_mpf_init2(&expected, 32);
  s21_mpf_set_ui(&x, 7);
  s21_mpf_set_ui(&expected, 7);
  ASSERT(s21_mpf_set_prec(&x, 32) == 0, "narrow exact: rc");
  ASSERT(x.prec == 32, "narrow exact: prec");
  ASSERT(s21_mpf_cmp(&x, &expected) == 0, "narrow exact: value");
  s21_mpf_clear(&x);
  s21_mpf_clear(&expected);

  /* 3. Сужение с округлением: 1 + 2^-5 при prec=128 → prec=4.
     Ближайшее 4-битное значение — 1.0 (расстояние 1/32),
     против 1.125 (расстояние 3/32). */
  s21_mpf_init2(&x, 128);
  s21_mpf_set_ui(&x, 1);
  {
    s21_mpf_t tiny;
    s21_mpf_init2(&tiny, 128);
    s21_mpf_set_ui(&tiny, 1);
    tiny.exp -= 5;
    s21_mpf_normalize(&tiny);
    s21_mpf_add(&x, &x, &tiny);
    s21_mpf_clear(&tiny);
  }
  ASSERT(s21_mpf_set_prec(&x, 4) == 0, "narrow round: rc");
  ASSERT(x.prec == 4, "narrow round: prec");
  s21_mpf_init2(&expected, 4);
  s21_mpf_set_ui(&expected, 1);
  ASSERT(s21_mpf_cmp(&x, &expected) == 0, "narrow round: → 1.0");
  s21_mpf_clear(&x);
  s21_mpf_clear(&expected);

  /* 4. Сужение с округлением вверх и переполнением мантиссы:
     1.111…1₂ → 2.0.  Проверяет, что mpf_add_one корректно
     обрабатывает случай выхода за границу prec. */
  s21_mpf_init2(&x, 128);
  s21_mpf_set_ui(&x, 1);
  {
    s21_mpf_t frac, pow;
    s21_mpf_init2(&frac, 128);
    s21_mpf_init2(&pow, 128);
    s21_mpf_set_ui(&pow, 1);
    /* Складываем 2^-1 + 2^-2 + ... + 2^-2, пока чуть-чуть не дойдём
       до 2. Проще: построить сумму 2 - 2^-100. */
    s21_mpf_t two;
    s21_mpf_init2(&two, 128);
    s21_mpf_set_ui(&two, 2);
    s21_mpf_set_ui(&pow, 1);
    pow.exp -= 100;
    s21_mpf_normalize(&pow);
    s21_mpf_sub(&x, &two, &pow); /* x = 2 - 2^-100, чуть меньше 2 */
    s21_mpf_clear(&frac);
    s21_mpf_clear(&pow);
    s21_mpf_clear(&two);
  }
  ASSERT(s21_mpf_set_prec(&x, 4) == 0, "narrow to 2.0: rc");
  ASSERT(x.prec == 4, "narrow to 2.0: prec");
  s21_mpf_init2(&expected, 4);
  s21_mpf_set_ui(&expected, 2);
  ASSERT(s21_mpf_cmp(&x, &expected) == 0, "narrow to 2.0: value");
  s21_mpf_clear(&x);
  s21_mpf_clear(&expected);

  /* 5. Same prec — no-op. */
  s21_mpf_init2(&x, 64);
  s21_mpf_set_ui(&x, 5);
  uint64_t saved_exp = (uint64_t)x.exp;
  ASSERT(s21_mpf_set_prec(&x, 64) == 0, "same: rc");
  ASSERT(x.prec == 64, "same: prec");
  ASSERT((uint64_t)x.exp == saved_exp, "same: exp");
  s21_mpf_clear(&x);

  /* 6. Clamp < 2. */
  s21_mpf_init2(&x, 32);
  s21_mpf_set_ui(&x, 1);
  ASSERT(s21_mpf_set_prec(&x, 0) == 0, "clamp 0: rc");
  ASSERT(x.prec == 2, "clamp 0: prec == 2");
  ASSERT(s21_mpf_set_prec(&x, 1) == 0, "clamp 1: rc");
  ASSERT(x.prec == 2, "clamp 1: prec == 2");
  s21_mpf_clear(&x);

  /* 7. Спецзначения. */
  s21_mpf_init2(&x, 32);
  s21_mpf_set_inf(&x, 1);
  ASSERT(s21_mpf_set_prec(&x, 128) == 0, "inf: rc");
  ASSERT(x.prec == 128, "inf: prec");
  ASSERT(s21_mpf_is_inf(&x) && x.sign == 1, "inf: value");
  s21_mpf_set_nan(&x);
  ASSERT(s21_mpf_set_prec(&x, 64) == 0, "nan: rc");
  ASSERT(s21_mpf_is_nan(&x), "nan: value");
  s21_mpf_set_zero(&x, 1);
  ASSERT(s21_mpf_set_prec(&x, 96) == 0, "zero: rc");
  ASSERT(s21_mpf_is_zero(&x) && x.sign == 1, "zero: value");
  s21_mpf_clear(&x);

  /* 8. NULL. */
  ASSERT(s21_mpf_set_prec(NULL, 128) == -1, "NULL");

  printf("[ok] set_prec (in-place, expand/narrow/RNDN, special)\n");
}

int main(void) {
  printf("=== s21_mpf: базовые тесты ===\n");
  test_init_clear();
  test_set_ui();
  test_set_si();
  test_set_d();
  test_special();
  test_precision();
  test_small_prec();
  test_round_narrow_across_limb_boundary();
  test_cmp();
  test_cmp_mixed_prec();
  test_add_sub();
  test_mul();
  test_nonmult32_prec();
  test_mul_many_prec();
  test_div();
  test_div_many_prec();
  test_div_powers_of_two();
  test_null_checks();
  test_bit_utils();
  test_sqrt();
  test_exp();
  test_log();
  test_pi();
  test_sin_cos();
  test_atan();
  test_tan_asin_acos();
  test_set_round();
  test_set_round_nonmult64();
  test_set_prec();
  test_ulp_diagnostics();
  printf("=== Все тесты прошли ===\n");
  return 0;
}
