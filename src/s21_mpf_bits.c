#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

int s21_mpf_get_bit(const s21_mpf_t *x, uint32_t pos) {
  if (pos >= x->prec) return 0;
  return (int)((x->limbs[pos / 64] >> (pos % 64)) & 1ULL);
}

void s21_mpf_set_bit(s21_mpf_t *x, uint32_t pos, int value) {
  if (x == NULL || pos >= x->prec) return;
  uint32_t word = pos / 64;
  uint32_t bit = pos % 64;
  if (value)
    x->limbs[word] |= (1ULL << bit);
  else
    x->limbs[word] &= ~(1ULL << bit);
}

int s21_mpf_msb(const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) return -1;
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    if (x->limbs[i] != 0) return i * 64 + (63 - __builtin_clzll(x->limbs[i]));
  }
  return -1;
}

int s21_mpf_lsb(const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) return -1;
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (size_t i = 0; i < count; i++) {
    if (x->limbs[i] != 0) return (int)(i * 64 + __builtin_ctzll(x->limbs[i]));
  }
  return -1;
}

/* Чтение лимба с трактовкой выхода за границы как нуля. */
static inline uint64_t fetch_limb(const uint64_t *src, size_t src_count,
                                  int idx) {
  if (idx < 0 || idx >= (int)src_count) return 0;
  return src[idx];
}

/* Сдвиг влево на shift бит. dst и src могут совпадать при
   dst_count == src_count; обход от старших индексов к младшим
   делает операцию безопасной in-place. */
void s21_mpf_shift_left_into(uint64_t *dst, size_t dst_count,
                             const uint64_t *src, size_t src_count, int shift) {
  if (shift == 0) {
    size_t n = dst_count < src_count ? dst_count : src_count;
    memcpy(dst, src, n * sizeof(uint64_t));
    for (size_t i = n; i < dst_count; i++) dst[i] = 0;
    return;
  }
  int word_shift = shift / 64;
  int bit_shift = shift % 64;

  for (int i = (int)dst_count - 1; i >= 0; i--) {
    int idx_lo = i - word_shift;
    int idx_hi = idx_lo - 1;
    uint64_t lo = fetch_limb(src, src_count, idx_lo);
    uint64_t hi = fetch_limb(src, src_count, idx_hi);
    uint64_t v = lo;
    if (bit_shift) v = (v << bit_shift) | (hi >> (64 - bit_shift));
    dst[i] = v;
  }
}

/* Сдвиг вправо на shift бит. Обход от младших индексов к старшим
   делает операцию безопасной in-place при dst_count == src_count. */
void s21_mpf_shift_right_into(uint64_t *dst, size_t dst_count,
                              const uint64_t *src, size_t src_count,
                              int shift) {
  if (shift == 0) {
    size_t n = dst_count < src_count ? dst_count : src_count;
    memcpy(dst, src, n * sizeof(uint64_t));
    for (size_t i = n; i < dst_count; i++) dst[i] = 0;
    return;
  }
  int word_shift = shift / 64;
  int bit_shift = shift % 64;

  for (size_t i = 0; i < dst_count; i++) {
    int idx_lo = (int)i + word_shift;
    int idx_hi = idx_lo + 1;
    uint64_t lo = fetch_limb(src, src_count, idx_lo);
    uint64_t hi = fetch_limb(src, src_count, idx_hi);
    uint64_t v = lo;
    if (bit_shift) v = (v >> bit_shift) | (hi << (64 - bit_shift));
    dst[i] = v;
  }
}
