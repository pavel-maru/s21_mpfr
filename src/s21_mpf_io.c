/* ============================================================
   s21_mpfr — I/O: парсинг строк (set_str) и вывод (get_str).
   ============================================================ */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

/* ============================================================
   set_str
   ============================================================ */

static int parse_special(const char *s, int sign, s21_mpf_t *x) {
  if (tolower((unsigned char)s[0]) == 'i' &&
      tolower((unsigned char)s[1]) == 'n' &&
      tolower((unsigned char)s[2]) == 'f') {
    size_t len = 3;
    if (tolower((unsigned char)s[3]) == 'i' &&
        tolower((unsigned char)s[4]) == 'n' &&
        tolower((unsigned char)s[5]) == 'i' &&
        tolower((unsigned char)s[6]) == 't' &&
        tolower((unsigned char)s[7]) == 'y') {
      len = 8;
    }
    if (s[len] == '\0') {
      s21_mpf_set_inf(x, sign);
      return 1;
    }
  }
  if (tolower((unsigned char)s[0]) == 'n' &&
      tolower((unsigned char)s[1]) == 'a' &&
      tolower((unsigned char)s[2]) == 'n' && s[3] == '\0') {
    s21_mpf_set_nan(x);
    return 1;
  }
  return 0;
}

/* res = 10^k для k >= 0. res->prec должен быть >= wp. */
static int pow10_mpf(s21_mpf_t *res, int64_t k, uint32_t wp) {
  s21_mpf_set_ui(res, 1);
  if (k == 0) return 0;

  s21_mpf_t base;
  s21_mpf_init2(&base, wp);
  s21_mpf_set_ui(&base, 10);

  while (k > 0) {
    if (k & 1) s21_mpf_mul_raw(res, res, &base);
    k >>= 1;
    if (k > 0) s21_mpf_mul_raw(&base, &base, &base);
  }

  s21_mpf_clear(&base);
  return 0;
}

int s21_mpf_set_str(s21_mpf_t *x, const char *str, int base) {
  if (x == NULL || str == NULL) return -1;
  if (base != 10) return -1;

  const char *p = str;

  while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;

  int sign = 0;
  if (*p == '+') {
    p++;
  } else if (*p == '-') {
    sign = 1;
    p++;
  }

  if (parse_special(p, sign, x)) return 0;

  uint32_t wp = x->prec + 64;

  s21_mpf_t M;
  s21_mpf_init2(&M, wp);
  s21_mpf_set_ui(&M, 0);

  int has_digit = 0;
  int64_t frac_digits = 0;
  int seen_nonzero = 0;

  while (isdigit((unsigned char)*p)) {
    has_digit = 1;
    int d = *p - '0';
    if (d != 0) seen_nonzero = 1;
    if (seen_nonzero) {
      s21_mpf_t ten, dval;
      s21_mpf_init2(&ten, wp);
      s21_mpf_init2(&dval, wp);
      s21_mpf_set_ui(&ten, 10);
      s21_mpf_mul_raw(&M, &M, &ten);
      if (d != 0) {
        s21_mpf_set_ui(&dval, (unsigned long)d);
        s21_mpf_add_raw(&M, &M, &dval);
      }
      s21_mpf_clear(&ten);
      s21_mpf_clear(&dval);
    }
    p++;
  }

  if (*p == '.') {
    p++;
    while (isdigit((unsigned char)*p)) {
      has_digit = 1;
      int d = *p - '0';
      if (d != 0) seen_nonzero = 1;
      if (seen_nonzero) {
        s21_mpf_t ten, dval;
        s21_mpf_init2(&ten, wp);
        s21_mpf_init2(&dval, wp);
        s21_mpf_set_ui(&ten, 10);
        s21_mpf_mul_raw(&M, &M, &ten);
        if (d != 0) {
          s21_mpf_set_ui(&dval, (unsigned long)d);
          s21_mpf_add_raw(&M, &M, &dval);
        }
        s21_mpf_clear(&ten);
        s21_mpf_clear(&dval);
      }
      frac_digits++;
      p++;
    }
  }

  if (!has_digit) {
    s21_mpf_clear(&M);
    return -1;
  }

  if (!seen_nonzero) {
    s21_mpf_clear(&M);
    s21_mpf_set_zero(x, sign);
    return 0;
  }

  int64_t e_part = 0;
  if (*p == 'e' || *p == 'E') {
    p++;
    int e_sign = 0;
    if (*p == '+') {
      p++;
    } else if (*p == '-') {
      e_sign = 1;
      p++;
    }
    int has_exp_digits = 0;
    while (isdigit((unsigned char)*p)) {
      has_exp_digits = 1;
      e_part = e_part * 10 + (*p - '0');
      p++;
    }
    if (!has_exp_digits) {
      s21_mpf_clear(&M);
      return -1;
    }
    if (e_sign) e_part = -e_part;
  }

  if (*p != '\0') {
    s21_mpf_clear(&M);
    return -1;
  }

  int64_t E = e_part - frac_digits;

  if (E != 0) {
    s21_mpf_t p10;
    s21_mpf_init2(&p10, wp);
    int64_t absE = E < 0 ? -E : E;
    pow10_mpf(&p10, absE, wp);
    if (E > 0)
      s21_mpf_mul_raw(&M, &M, &p10);
    else
      s21_mpf_div_raw(&M, &M, &p10);
    s21_mpf_clear(&p10);
  }

  M.sign = sign;
  s21_mpf_set_round(x, &M, S21_MPF_RNDN);

  s21_mpf_clear(&M);
  return 0;
}

/* ============================================================
   get_str
   ============================================================ */

int s21_mpf_get_str(char *buf, size_t bufsize, int *exp_out,
                    int base, uint32_t n_digits,
                    const s21_mpf_t *x, s21_mpf_rnd_t rnd) {
  if (buf == NULL || x == NULL || exp_out == NULL) return -1;
  if (base != 10) return -1;
  if (n_digits == 0) n_digits = 20;

  /* Спецзначения */
  if (x->kind == S21_MPF_NAN) {
    if (bufsize < 4) return -1;
    strcpy(buf, "nan");
    *exp_out = 0;
    return 0;
  }
  if (x->kind == S21_MPF_INF) {
    if (bufsize < 5) return -1;
    strcpy(buf, x->sign ? "-inf" : "inf");
    *exp_out = 0;
    return 0;
  }
  if (x->kind == S21_MPF_ZERO) {
    if (bufsize < 3) return -1;
    strcpy(buf, x->sign ? "-0" : "0");
    *exp_out = 0;
    return 0;
  }

  size_t need = (size_t)n_digits + 2; /* sign + digits + '\0' */
  if (bufsize < need) return -1;

  uint32_t wp = n_digits * 4 + 192;

  size_t digits_start = 0;
  if (x->sign) buf[digits_start++] = '-';

  /* E_est ≈ floor(log10(x)); после нормализации y ∈ [1, 10), E — точное. */
  int64_t E = (int64_t)((double)x->exp * 0.30102999566398119521);

  s21_mpf_t xw, pow10, y;
  s21_mpf_init2(&xw, wp);
  s21_mpf_init2(&pow10, wp);
  s21_mpf_init2(&y, wp);
  s21_mpf_set(&xw, x);
  /* Знак уже учтён в buf (см. digits_start выше).  Для дальнейшей
     работы нужен модуль: цикл нормализации y в [1, 10) не
     завершается для отрицательного y (y < 1 всегда истинно). */
  xw.sign = 0;

  if (E >= 0) {
    pow10_mpf(&pow10, E, wp);
    s21_mpf_div_raw(&y, &xw, &pow10);
  } else {
    pow10_mpf(&pow10, -E, wp);
    s21_mpf_mul_raw(&y, &xw, &pow10);
  }

  s21_mpf_t one, ten;
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&ten, wp);
  s21_mpf_set_ui(&one, 1);
  s21_mpf_set_ui(&ten, 10);

  while (s21_mpf_cmp(&y, &ten) >= 0) {
    s21_mpf_div_raw(&y, &y, &ten);
    E++;
  }
  while (s21_mpf_cmp(&y, &one) < 0) {
    s21_mpf_mul_raw(&y, &y, &ten);
    E--;
  }

  /* Извлечение n_digits цифр: d = floor(y), y = (y - d) * 10. */
  s21_mpf_stack_t tmp_s;
  s21_mpf_stack_init(&tmp_s, wp);

  for (uint32_t i = 0; i < n_digits; i++) {
    int d = 0;
    for (int k = 9; k >= 0; k--) {
      s21_mpf_set_ui(&tmp_s.mpf, (unsigned long)k);
      if (s21_mpf_cmp(&y, &tmp_s.mpf) >= 0) {
        d = k;
        break;
      }
    }
    buf[digits_start + i] = (char)('0' + d);
    if (d > 0) {
      s21_mpf_set_ui(&tmp_s.mpf, (unsigned long)d);
      s21_mpf_sub_raw(&y, &y, &tmp_s.mpf);
    }
    s21_mpf_mul_raw(&y, &y, &ten);
  }
  s21_mpf_stack_clear(&tmp_s);

  /* Округление: y ∈ [0, 10) — остаток. */
  s21_mpf_t five;
  s21_mpf_init2(&five, wp);
  s21_mpf_set_ui(&five, 5);
  int cmp5 = s21_mpf_cmp(&y, &five);
  s21_mpf_clear(&five);

  int round_up = 0;
  int nonzero_rem = !s21_mpf_is_zero(&y);

  switch (rnd) {
    case S21_MPF_RNDN:
      if (cmp5 > 0)
        round_up = 1;
      else if (cmp5 == 0)
        round_up = ((buf[digits_start + n_digits - 1] - '0') & 1);
      break;
    case S21_MPF_RNDZ:
      round_up = 0;
      break;
    case S21_MPF_RNDU:
      round_up = (nonzero_rem && !x->sign);
      break;
    case S21_MPF_RNDD:
      round_up = (nonzero_rem && x->sign);
      break;
  }

  if (round_up) {
    int i = (int)n_digits - 1;
    while (i >= 0) {
      size_t idx = digits_start + (size_t)i;
      if (buf[idx] < '9') {
        buf[idx]++;
        break;
      }
      buf[idx] = '0';
      i--;
    }
    if (i < 0) {
      /* Перенос через все 9: 999… → 100…0, E++. */
      for (uint32_t j = 0; j < n_digits; j++) {
        buf[digits_start + j] = '0';
      }
      buf[digits_start] = '1';
      E++;
    }
  }

  buf[digits_start + n_digits] = '\0';
  *exp_out = (int)(E + 1);

  s21_mpf_clear(&xw);
  s21_mpf_clear(&pow10);
  s21_mpf_clear(&y);
  s21_mpf_clear(&one);
  s21_mpf_clear(&ten);
  return 0;
}
