#include <string.h>

#include "s21_mpf.h"

int s21_mpf_get_bit(const s21_mpf_t *x, uint32_t pos) {
  if (pos >= x->prec) return 0;
  uint32_t word = pos / 64;
  uint32_t bit = pos % 64;
  return (int)((x->limbs[word] >> bit) & 1ULL);
}

void s21_mpf_set_bit(s21_mpf_t *x, uint32_t pos, int value) {
  if (x == NULL || pos >= x->prec) return;
  uint32_t word = pos / 64;
  uint32_t bit = pos % 64;
  if (value) {
    x->limbs[word] |= (1ULL << bit);
  } else {
    x->limbs[word] &= ~(1ULL << bit);
  }
}

int s21_mpf_msb(const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) return -1;
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (int i = (int)count - 1; i >= 0; i--) {
    if (x->limbs[i] != 0) {
      return i * 64 + (63 - __builtin_clzll(x->limbs[i]));
    }
  }
  return -1;
}

int s21_mpf_lsb(const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) return -1;
  size_t count = s21_mpf_limbs_for_prec(x->prec);
  for (size_t i = 0; i < count; i++) {
    if (x->limbs[i] != 0) {
      return (int)(i * 64 + __builtin_ctzll(x->limbs[i]));
    }
  }
  return -1;
}

void s21_mpf_shift_left_into(uint64_t *dst, const uint64_t *src, size_t count,
                             int shift) {
  if (shift == 0) {
    memcpy(dst, src, count * sizeof(uint64_t));
    return;
  }
  int word_shift = shift / 64;
  int bit_shift = shift % 64;

  for (int i = (int)count - 1; i >= 0; i--) {
    int src_idx = i - word_shift;
    if (src_idx < 0) {
      dst[i] = 0;
    } else {
      uint64_t v = src[src_idx];
      if (bit_shift && src_idx > 0) {
        v = (v << bit_shift) | (src[src_idx - 1] >> (64 - bit_shift));
      } else if (bit_shift) {
        v <<= bit_shift;
      }
      dst[i] = v;
    }
  }
}
