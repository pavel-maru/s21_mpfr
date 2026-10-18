/* ============================================================
   s21_mpfr — I/O: парсинг строк (set_str) и вывод (get_str).
   ============================================================ */

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"

/* ============================================================
   Цифры в произвольном основании.

   base ∈ [2, 36], алфавит 0-9a-z (регистронезависимо на входе,
   строчные на выходе).  Для base ≤ 10 заглавные буквы не
   являются цифрами (digit_value('A') = 10 ≥ base).
   ============================================================ */

static int digit_value(char c) {
  unsigned char uc = (unsigned char)c;
  if (uc >= '0' && uc <= '9') return uc - '0';
  if (uc >= 'a' && uc <= 'z') return uc - 'a' + 10;
  if (uc >= 'A' && uc <= 'Z') return uc - 'A' + 10;
  return -1;
}

static char digit_char(int d) {
  return (char)(d < 10 ? '0' + d : 'a' + d - 10);
}

static int is_digit_in_base(char c, int base) {
  int d = digit_value(c);
  return d >= 0 && d < base;
}

/* ============================================================
   set_str
   ============================================================ */

static int parse_special(const char *s, int sign, s21_mpf_t *x) {
  size_t slen = strlen(s);

  if (slen >= 3 && tolower((unsigned char)s[0]) == 'i' &&
      tolower((unsigned char)s[1]) == 'n' &&
      tolower((unsigned char)s[2]) == 'f') {
    size_t len = 3;
    if (slen >= 8 && tolower((unsigned char)s[3]) == 'i' &&
        tolower((unsigned char)s[4]) == 'n' &&
        tolower((unsigned char)s[5]) == 'i' &&
        tolower((unsigned char)s[6]) == 't' &&
        tolower((unsigned char)s[7]) == 'y') {
      len = 8;
    }
    if (slen == len) {
      s21_mpf_set_inf(x, sign);
      return 1;
    }
  }
  if (slen == 3 && tolower((unsigned char)s[0]) == 'n' &&
      tolower((unsigned char)s[1]) == 'a' &&
      tolower((unsigned char)s[2]) == 'n') {
    s21_mpf_set_nan(x);
    return 1;
  }
  return 0;
}

/* res = base^k для k >= 0.  res->prec должен быть >= wp. */
static int pow_base_mpf(s21_mpf_t *res, unsigned base, int64_t k, uint32_t wp) {
  s21_mpf_set_ui(res, 1);
  if (k == 0) return 0;

  s21_mpf_t b;
  s21_mpf_init2(&b, wp);
  s21_mpf_set_ui(&b, (unsigned long)base);

  while (k > 0) {
    if (k & 1) s21_mpf_mul_raw(res, res, &b);
    k >>= 1;
    if (k > 0) s21_mpf_mul_raw(&b, &b, &b);
  }

  s21_mpf_clear(&b);
  return 0;
}

int s21_mpf_set_str(s21_mpf_t *x, const char *str, int base) {
  if (x == NULL || str == NULL) return -1;
  if (base < 2 || base > 36) return -1;

  const char *p = str;

  while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;

  int sign = 0;
  if (*p == '+') {
    p++;
  } else if (*p == '-') {
    sign = 1;
    p++;
  }

  /* C99 hexfloat prefix 0x / 0X -- only for base 16. */
  if (base == 16 && p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) p += 2;

  if (parse_special(p, sign, x)) return 0;

  uint32_t wp = x->prec + 64;

  s21_mpf_t M, bval;
  s21_mpf_init2(&M, wp);
  s21_mpf_init2(&bval, wp);
  s21_mpf_set_ui(&M, 0);
  s21_mpf_set_ui(&bval, (unsigned long)base);

  int has_digit = 0;
  int64_t frac_digits = 0;
  int seen_nonzero = 0;

  while (is_digit_in_base(*p, base)) {
    has_digit = 1;
    int d = digit_value(*p);
    if (d != 0) seen_nonzero = 1;
    if (seen_nonzero) {
      s21_mpf_t dval;
      s21_mpf_init2(&dval, wp);
      s21_mpf_mul_raw(&M, &M, &bval);
      if (d != 0) {
        s21_mpf_set_ui(&dval, (unsigned long)d);
        s21_mpf_add_raw(&M, &M, &dval);
      }
      s21_mpf_clear(&dval);
    }
    p++;
  }

  if (*p == '.') {
    p++;
    while (is_digit_in_base(*p, base)) {
      has_digit = 1;
      int d = digit_value(*p);
      if (d != 0) seen_nonzero = 1;
      if (seen_nonzero) {
        s21_mpf_t dval;
        s21_mpf_init2(&dval, wp);
        s21_mpf_mul_raw(&M, &M, &bval);
        if (d != 0) {
          s21_mpf_set_ui(&dval, (unsigned long)d);
          s21_mpf_add_raw(&M, &M, &dval);
        }
        s21_mpf_clear(&dval);
      }
      frac_digits++;
      p++;
    }
  }

  if (!has_digit) {
    s21_mpf_clear(&M);
    s21_mpf_clear(&bval);
    return -1;
  }

  /* C99 hexfloat: p / P-экспонента -- десятичное число, задающее
     степень двойки.  Только для base 16.  Значение итога:
     M · base^(-frac_digits) · 2^(p_exp).  Применяется коррекцией
     M.exp в самом конце -- умножение на степень двойки точно. */
  int64_t p_exp = 0;
  int has_p_exp = 0;
  if (base == 16 && (*p == 'p' || *p == 'P')) {
    p++;
    int pe_sign = 0;
    if (*p == '+') {
      p++;
    } else if (*p == '-') {
      pe_sign = 1;
      p++;
    }
    int has_pe_digits = 0;
    while (isdigit((unsigned char)*p)) {
      has_pe_digits = 1;
      if (p_exp > (INT64_MAX - (*p - '0')) / 10) {
        s21_mpf_clear(&M);
        s21_mpf_clear(&bval);
        return -1;
      }
      p_exp = p_exp * 10 + (*p - '0');
      p++;
    }
    if (!has_pe_digits) {
      s21_mpf_clear(&M);
      s21_mpf_clear(&bval);
      return -1;
    }
    if (pe_sign) p_exp = -p_exp;
    has_p_exp = 1;
  }

  /* Десятичная экспонента e / E -- только для base 10. */
  int64_t e_part = 0;
  if (base == 10 && (*p == 'e' || *p == 'E')) {
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
      s21_mpf_clear(&bval);
      return -1;
    }
    if (e_sign) e_part = -e_part;
  }

  if (*p != '\0') {
    s21_mpf_clear(&M);
    s21_mpf_clear(&bval);
    return -1;
  }

  if (!seen_nonzero) {
    s21_mpf_clear(&M);
    s21_mpf_clear(&bval);
    s21_mpf_set_zero(x, sign);
    return 0;
  }

  int64_t E = e_part - frac_digits;

  if (E != 0) {
    s21_mpf_t pb;
    s21_mpf_init2(&pb, wp);
    int64_t absE = E < 0 ? -E : E;
    pow_base_mpf(&pb, (unsigned)base, absE, wp);
    if (E > 0)
      s21_mpf_mul_raw(&M, &M, &pb);
    else
      s21_mpf_div_raw(&M, &M, &pb);
    s21_mpf_clear(&pb);
  }

  /* Применяем hexfloat-экспоненту: M *= 2^(p_exp).  Сумма M.exp
     и p_exp может переполнить int64_t -- обрабатываем явно. */
  if (has_p_exp) {
    int64_t new_exp;
    if (__builtin_add_overflow(M.exp, p_exp, &new_exp)) {
      s21_mpf_clear(&M);
      s21_mpf_clear(&bval);
      return -1;
    }
    M.exp = new_exp;
  }

  M.sign = sign;
  s21_mpf_set_round(x, &M, S21_MPF_RNDN);

  s21_mpf_clear(&M);
  s21_mpf_clear(&bval);
  return 0;
}

/* ============================================================
   get_str
   ============================================================ */

/* Обработка NAN / INF / ZERO.  Возвращает 1, если результат
   уже записан, 0 — если x нормальное. */
static int get_str_special(char *buf, size_t bufsize, int *exp_out,
                           const s21_mpf_t *x) {
  if (x->kind == S21_MPF_NAN) {
    if (bufsize < 4) return -1;
    strcpy(buf, "nan");
    *exp_out = 0;
    return 1;
  }
  if (x->kind == S21_MPF_INF) {
    size_t need = x->sign ? 5u : 4u;
    if (bufsize < need) return -1;
    strcpy(buf, x->sign ? "-inf" : "inf");
    *exp_out = 0;
    return 1;
  }
  if (x->kind == S21_MPF_ZERO) {
    if (bufsize < 3) return -1;
    strcpy(buf, x->sign ? "-0" : "0");
    *exp_out = 0;
    return 1;
  }
  return 0;
}

/* Извлечение n_digits цифр по основанию base из y ∈ [1, base)
   в buf[digits_start..).  y изменяется: на выходе это остаток
   y ∈ [0, base) после последнего умножения на base. */
static void get_str_extract(char *buf, size_t digits_start, uint32_t n_digits,
                            s21_mpf_t *y, unsigned base, uint32_t wp) {
  s21_mpf_t base_val;
  s21_mpf_stack_t tmp_s;
  s21_mpf_init2(&base_val, wp);
  s21_mpf_stack_init(&tmp_s, wp);
  s21_mpf_set_ui(&base_val, (unsigned long)base);

  for (uint32_t i = 0; i < n_digits; i++) {
    int d = 0;
    for (int k = (int)base - 1; k >= 0; k--) {
      s21_mpf_set_ui(&tmp_s.mpf, (unsigned long)k);
      if (s21_mpf_cmp(y, &tmp_s.mpf) >= 0) {
        d = k;
        break;
      }
    }
    buf[digits_start + i] = digit_char(d);
    if (d > 0) {
      s21_mpf_set_ui(&tmp_s.mpf, (unsigned long)d);
      s21_mpf_sub_raw(y, y, &tmp_s.mpf);
    }
    s21_mpf_mul_raw(y, y, &base_val);
  }

  s21_mpf_stack_clear(&tmp_s);
  s21_mpf_clear(&base_val);
}

/* Решение об округлении вверх.  rem = y после извлечения цифр,
   лежит в [0, base).  last_digit_val — значение последней
   записанной цифры (для ties-to-even). */
static int get_str_should_round_up(const s21_mpf_t *rem, int last_digit_val,
                                   int sign, s21_mpf_rnd_t rnd, unsigned base,
                                   uint32_t wp) {
  s21_mpf_t half;
  s21_mpf_init2(&half, wp);
  s21_mpf_set_ui(&half, (unsigned long)base);
  half.exp -= 1; /* base / 2 */
  s21_mpf_normalize(&half);
  int cmp_half = s21_mpf_cmp(rem, &half);
  s21_mpf_clear(&half);

  int nonzero_rem = !s21_mpf_is_zero(rem);

  switch (rnd) {
    case S21_MPF_RNDN:
      if (cmp_half > 0) return 1;
      if (cmp_half == 0) return (last_digit_val & 1);
      return 0;
    case S21_MPF_RNDZ:
      return 0;
    case S21_MPF_RNDU:
      return nonzero_rem && !sign;
    case S21_MPF_RNDD:
      return nonzero_rem && sign;
  }
  return 0;
}

/* Инкремент последней цифры с переносом через все (base−1).
   При переносе через всю строку ((base−1)(base−1)… → 100…0)
   увеличивает *E. */
static void get_str_round_up(char *buf, uint32_t n_digits, unsigned base,
                             int64_t *E) {
  int i = (int)n_digits - 1;
  while (i >= 0) {
    int d = digit_value(buf[i]);
    if (d < (int)base - 1) {
      buf[i] = digit_char(d + 1);
      return;
    }
    buf[i] = '0';
    i--;
  }
  for (uint32_t j = 0; j < n_digits; j++) buf[j] = '0';
  buf[0] = '1';
  (*E)++;
}

int s21_mpf_get_str(char *buf, size_t bufsize, int *exp_out, int base,
                    uint32_t n_digits, const s21_mpf_t *x, s21_mpf_rnd_t rnd) {
  if (buf == NULL || x == NULL || exp_out == NULL) return -1;
  if (base < 2 || base > 36) return -1;
  if (n_digits == 0) n_digits = 20;

  int rc_special = get_str_special(buf, bufsize, exp_out, x);
  if (rc_special > 0) return 0;
  if (rc_special < 0) return -1;

  size_t need = (size_t)n_digits + 2;
  if (bufsize < need) return -1;

  /* Рабочая точность: log2(base) бит на цифру плюс 192 запаса. */
  double log2_base = log((double)base) / log(2.0);
  uint32_t wp = (uint32_t)((double)n_digits * log2_base) + 192;

  size_t digits_start = 0;
  if (x->sign) buf[digits_start++] = '-';

  /* E_est ≈ floor(log_base(x)) = x.exp / log2(base). */
  int64_t E = (int64_t)((double)x->exp / log2_base);

  s21_mpf_t xw, powb, y, one, base_val;
  s21_mpf_init2(&xw, wp);
  s21_mpf_init2(&powb, wp);
  s21_mpf_init2(&y, wp);
  s21_mpf_init2(&one, wp);
  s21_mpf_init2(&base_val, wp);
  s21_mpf_set(&xw, x);
  /* Знак уже учтён в buf.  Для дальнейшей работы нужен модуль:
     цикл нормализации y ∈ [1, base) не завершается для y < 0. */
  xw.sign = 0;

  if (E >= 0) {
    pow_base_mpf(&powb, (unsigned)base, E, wp);
    s21_mpf_div_raw(&y, &xw, &powb);
  } else {
    pow_base_mpf(&powb, (unsigned)base, -E, wp);
    s21_mpf_mul_raw(&y, &xw, &powb);
  }

  s21_mpf_set_ui(&one, 1);
  s21_mpf_set_ui(&base_val, (unsigned long)base);
  while (s21_mpf_cmp(&y, &base_val) >= 0) {
    s21_mpf_div_raw(&y, &y, &base_val);
    E++;
  }
  while (s21_mpf_cmp(&y, &one) < 0) {
    s21_mpf_mul_raw(&y, &y, &base_val);
    E--;
  }

  get_str_extract(buf, digits_start, n_digits, &y, (unsigned)base, wp);

  int last_digit_val = digit_value(buf[digits_start + n_digits - 1]);
  if (get_str_should_round_up(&y, last_digit_val, x->sign, rnd, (unsigned)base,
                              wp)) {
    get_str_round_up(buf + digits_start, n_digits, (unsigned)base, &E);
  }

  buf[digits_start + n_digits] = '\0';
  *exp_out = (int)(E + 1);

  s21_mpf_clear(&xw);
  s21_mpf_clear(&powb);
  s21_mpf_clear(&y);
  s21_mpf_clear(&one);
  s21_mpf_clear(&base_val);
  return 0;
}
