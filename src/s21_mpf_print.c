#include <inttypes.h>
#include <math.h>
#include <stdio.h>

#include "s21_mpf.h"

void s21_mpf_fprint(FILE *stream, const s21_mpf_t *x) {
  if (stream == NULL || x == NULL) return;
  if (x->kind == S21_MPF_NAN) {
    fprintf(stream, "NaN\n");
    return;
  }
  if (x->kind == S21_MPF_INF) {
    fprintf(stream, "%sinf\n", x->sign ? "-" : "+");
    return;
  }
  if (x->kind == S21_MPF_ZERO) {
    fprintf(stream, "%s0\n", x->sign ? "-" : "+");
    return;
  }
  fprintf(stream, "%s0x", x->sign ? "-" : "+");
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--)
    fprintf(stream, "%016" PRIx64, x->limbs[i]);
  fprintf(stream, "p%" PRId64 " [prec=%u]\n", x->exp, x->prec);
}

void s21_mpf_fprint_d(FILE *stream, const s21_mpf_t *x) {
  if (stream == NULL || x == NULL) return;
  if (x->kind == S21_MPF_NAN) {
    fprintf(stream, "nan\n");
    return;
  }
  if (x->kind == S21_MPF_INF) {
    fprintf(stream, "%sinf\n", x->sign ? "-" : "+");
    return;
  }
  if (x->kind == S21_MPF_ZERO) {
    fprintf(stream, "%s0\n", x->sign ? "-" : "+");
    return;
  }

  size_t count = s21_mpf_limbs_for_prec(x->prec);
  double val = 0.0;
  int taken = 0;
  int i = (int)count - 1;
  for (; i >= 0 && taken < 2; i--, taken++) {
    val = val * 18446744073709551616.0 + (double)x->limbs[i];
  }
  int64_t missed_limbs = (int64_t)(i + 1);
  int64_t e = x->exp - (int64_t)x->prec + 64 * missed_limbs;
  double result = ldexp(val, (int)e);
  fprintf(stream, "%s%.17g\n", x->sign ? "-" : "+", result);
}

void s21_mpf_print(const s21_mpf_t *x) { s21_mpf_fprint(stdout, x); }
void s21_mpf_print_d(const s21_mpf_t *x) { s21_mpf_fprint_d(stdout, x); }