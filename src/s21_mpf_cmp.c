#include "s21_mpf.h"

/* Быстрое сравнение мантисс при равных prec. */
static int s21_mpf_cmp_mant(const s21_mpf_t *x, const s21_mpf_t *y) {
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    if (x->limbs[i] < y->limbs[i]) return -1;
    if (x->limbs[i] > y->limbs[i]) return 1;
  }
  return 0;
}

/* pos-й бит мантиссы v, расширенной до n бит (n >= v->prec).
   Биты младше (n - v->prec) равны нулю; остальные совпадают
   с мантиссой v. */
static int s21_mpf_mant_bit_padded(const s21_mpf_t *v, uint32_t pos,
                                   uint32_t n) {
  uint32_t shift = n - v->prec;
  if (pos < shift) return 0;
  uint32_t idx = pos - shift;
  return (int)((v->limbs[idx / 64] >> (idx % 64)) & 1ULL);
}

/* Сравнение мантисс при разных prec. Обе мантиссы нормализованы
   (MSB на позиции prec-1), при равных exp вес MSB одинаков, поэтому
   выравниваем их по общему размеру n = max(prec_x, prec_y) сдвигом
   влево и сравниваем побитово от старшего к младшему. */
static int s21_mpf_cmp_mant_padded(const s21_mpf_t *x, const s21_mpf_t *y) {
  uint32_t n = x->prec > y->prec ? x->prec : y->prec;
  for (int64_t i = (int64_t)n - 1; i >= 0; i--) {
    int bx = s21_mpf_mant_bit_padded(x, (uint32_t)i, n);
    int by = s21_mpf_mant_bit_padded(y, (uint32_t)i, n);
    if (bx != by) return bx < by ? -1 : 1;
  }
  return 0;
}

int s21_mpf_cmp_abs(const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) return 0;

  int xinf = (x->kind == S21_MPF_INF);
  int yinf = (y->kind == S21_MPF_INF);
  if (xinf && yinf) return 0;
  if (xinf) return 1;
  if (yinf) return -1;

  int xzero = (x->kind == S21_MPF_ZERO);
  int yzero = (y->kind == S21_MPF_ZERO);
  if (xzero && yzero) return 0;
  if (xzero) return -1;
  if (yzero) return 1;

  if (x->exp < y->exp) return -1;
  if (x->exp > y->exp) return 1;

  if (x->prec == y->prec) return s21_mpf_cmp_mant(x, y);
  return s21_mpf_cmp_mant_padded(x, y);
}

int s21_mpf_cmp(const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) return 0;

  int xs = x->kind == S21_MPF_ZERO ? 0 : (x->sign ? -1 : 1);
  int ys = y->kind == S21_MPF_ZERO ? 0 : (y->sign ? -1 : 1);

  if (x->kind != S21_MPF_ZERO && y->kind != S21_MPF_ZERO && xs != ys) {
    return xs < ys ? -1 : 1;
  }

  int cmp_abs = s21_mpf_cmp_abs(x, y);
  if (xs < 0) return -cmp_abs;
  return cmp_abs;
}

int s21_mpf_equal(const s21_mpf_t *x, const s21_mpf_t *y) {
  if (x->kind == S21_MPF_NAN || y->kind == S21_MPF_NAN) return 0;
  if (x->kind != y->kind) return 0;
  if (x->kind == S21_MPF_ZERO) return 1;
  if (x->kind == S21_MPF_INF) return x->sign == y->sign;
  if (x->sign != y->sign) return 0;
  if (x->prec != y->prec) return 0;
  if (x->exp != y->exp) return 0;
  return s21_mpf_cmp_mant(x, y) == 0;
}

int s21_mpf_zero_p(const s21_mpf_t *x) { return x->kind == S21_MPF_ZERO; }

int s21_mpf_integer_p(const s21_mpf_t *x) {
  if (x->kind == S21_MPF_ZERO) return 1;
  if (x->kind != S21_MPF_NORMAL) return 0;
  if (x->exp >= (int64_t)x->prec) return 1;

  int64_t fractional_bits = (int64_t)x->prec - x->exp;
  if (fractional_bits <= 0) return 1;
  if (fractional_bits >= (int64_t)x->prec) return 0;

  uint32_t fb = (uint32_t)fractional_bits;
  uint32_t full_words = fb / 64;
  uint32_t rem_bits = fb % 64;

  for (uint32_t i = 0; i < full_words; i++) {
    if (x->limbs[i] != 0) return 0;
  }
  if (rem_bits) {
    uint64_t mask = (1ULL << rem_bits) - 1;
    if ((x->limbs[full_words] & mask) != 0) return 0;
  }
  return 1;
}
