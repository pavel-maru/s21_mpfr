/* ============================================================
   s21_mpfr — интерактивная утилита (REPL).

   Формат ввода:
     <func> [arg1] [arg2] [prec]
     <number>            — то же, что 'show <number>'

   Если arg1/arg2/prec пропущены — утилита спросит их отдельно.
   prec — точность в битах; по умолчанию 256.  Суффикс 'd' у
   числа означает десятичные цифры: 77d == 256 бит.

   Аргументы парсятся через s21_mpf_set_str, а не через strtod:
   никакого ограничения 53 битами, 3.14 — это ровно 3.14.

   Тригонометрия использует кэш π: при первом вызове sin/cos/tan
   считается π на текущую точность, дальше переиспользуется.
   Кэш сбрасывается при смене точности.
   ============================================================ */

#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <locale.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "s21_mpf.h"
#include "s21_mpf_internal.h"   /* S21_MPF_MAX_ARG_BITS */

#define S21_LOG2_10 3.3219280948873626
#define S21_LOG10_2 0.30102999566398120

/* ------------------------------------------------------------
   Кэш π для тригонометрических функций.

   Хранит π на точности (последний запрошенный prec + 128).
   При запросе другой точности — пересчитывает.
   ------------------------------------------------------------ */

static s21_mpf_t g_pi_cache;
static uint32_t g_pi_cache_prec = 0;
static bool g_pi_cache_valid = false;

static void pi_cache_invalidate(void) {
  if (g_pi_cache_valid) {
    s21_mpf_clear(&g_pi_cache);
    g_pi_cache_valid = false;
    g_pi_cache_prec = 0;
  }
}

/* Сколько «лишних» бит нужно добавить к точности π из-за cancelation
   в reduce_2pi.  Для |x| < 2 это 0; растёт как log2(|x|).  То же,
   что trig_extra_bits в s21_mpf_func.c. */
static int64_t arg_extra_bits(const s21_mpf_t *x) {
  if (x->kind != S21_MPF_NORMAL) return 0;
  int64_t E = x->exp - 1;
  return E > 0 ? E : 0;
}

/* Возвращает π на точности не меньше res_prec + 128 + E.
   При E > S21_MPF_MAX_ARG_BITS возвращает NULL — вызывающий
   должен откатиться на обычную функцию без явной π. */
static const s21_mpf_t *pi_cache_get(uint32_t res_prec, int64_t E) {
  if (E > S21_MPF_MAX_ARG_BITS) return NULL;
  uint32_t need = res_prec + 128 + (uint32_t)E;
  if (!g_pi_cache_valid || g_pi_cache_prec < need) {
    pi_cache_invalidate();
    s21_mpf_init2(&g_pi_cache, need);
    s21_mpf_pi(&g_pi_cache);
    g_pi_cache_prec = need;
    g_pi_cache_valid = true;
  }
  return &g_pi_cache;
}

/* ------------------------------------------------------------
   Обёртки тригонометрии через кэш π.
   ------------------------------------------------------------ */

static int wrap_sin(s21_mpf_t *r, const s21_mpf_t *x) {
  const s21_mpf_t *pi = pi_cache_get(r->prec, arg_extra_bits(x));
  if (pi == NULL) return s21_mpf_sin(r, x);
  return s21_mpf_sin_with_pi(r, x, pi);
}

static int wrap_cos(s21_mpf_t *r, const s21_mpf_t *x) {
  const s21_mpf_t *pi = pi_cache_get(r->prec, arg_extra_bits(x));
  if (pi == NULL) return s21_mpf_cos(r, x);
  return s21_mpf_cos_with_pi(r, x, pi);
}

static int wrap_tan(s21_mpf_t *r, const s21_mpf_t *x) {
  const s21_mpf_t *pi = pi_cache_get(r->prec, arg_extra_bits(x));
  if (pi == NULL) return s21_mpf_tan(r, x);
  return s21_mpf_tan_with_pi(r, x, pi);
}

static int wrap_atan(s21_mpf_t *r, const s21_mpf_t *x) {
  const s21_mpf_t *pi = pi_cache_get(r->prec, arg_extra_bits(x));
  if (pi == NULL) return s21_mpf_atan(r, x);
  return s21_mpf_atan_with_pi(r, x, pi);
}

static int wrap_asin(s21_mpf_t *r, const s21_mpf_t *x) {
  const s21_mpf_t *pi = pi_cache_get(r->prec, arg_extra_bits(x));
  if (pi == NULL) return s21_mpf_asin(r, x);
  return s21_mpf_asin_with_pi(r, x, pi);
}

static int wrap_acos(s21_mpf_t *r, const s21_mpf_t *x) {
  const s21_mpf_t *pi = pi_cache_get(r->prec, arg_extra_bits(x));
  if (pi == NULL) return s21_mpf_acos(r, x);
  return s21_mpf_acos_with_pi(r, x, pi);
}

static int wrap_abs(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_abs(res, x);
  return 0;
}

static int wrap_neg(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_neg(res, x);
  return 0;
}

/* ------------------------------------------------------------
   Таблицы функций.
   ------------------------------------------------------------ */

typedef int (*s21_unary_fn)(s21_mpf_t *, const s21_mpf_t *);
typedef int (*s21_binary_fn)(s21_mpf_t *, const s21_mpf_t *, const s21_mpf_t *);

typedef struct {
  const char *name;
  s21_unary_fn fn;
} unary_entry_t;

typedef struct {
  const char *name;
  s21_binary_fn fn;
} binary_entry_t;

static const unary_entry_t unary_table[] = {
    { "sqrt", s21_mpf_sqrt },
    { "exp",  s21_mpf_exp  },
    { "log",  s21_mpf_log  },
    { "sin",  wrap_sin     },
    { "cos",  wrap_cos     },
    { "tan",  wrap_tan     },
    { "atan", wrap_atan    },
    { "asin", wrap_asin    },
    { "acos", wrap_acos    },
    { "abs",  wrap_abs     },
    { "neg",  wrap_neg     },
};
static const size_t unary_count = sizeof(unary_table) / sizeof(unary_table[0]);

static const binary_entry_t binary_table[] = {
    { "add", s21_mpf_add },
    { "sub", s21_mpf_sub },
    { "mul", s21_mpf_mul },
    { "div", s21_mpf_div },
};
static const size_t binary_count =
    sizeof(binary_table) / sizeof(binary_table[0]);

/* ------------------------------------------------------------
   Глобальное состояние REPL.
   ------------------------------------------------------------ */

static uint32_t g_default_prec = 256;
static uint32_t g_dec_digits = 20;
static bool g_hex_output = false;
static bool g_interactive = false;

/* ------------------------------------------------------------
   Мелкие утилиты.
   ------------------------------------------------------------ */

static s21_unary_fn find_unary(const char *name) {
  for (size_t i = 0; i < unary_count; i++)
    if (strcmp(unary_table[i].name, name) == 0) return unary_table[i].fn;
  return NULL;
}

static s21_binary_fn find_binary(const char *name) {
  for (size_t i = 0; i < binary_count; i++)
    if (strcmp(binary_table[i].name, name) == 0) return binary_table[i].fn;
  return NULL;
}

static void print_dec(const s21_mpf_t *x, uint32_t n_digits) {
  if (n_digits == 0) n_digits = 20;
  if (n_digits > 1000) n_digits = 1000;

  char buf[1024];
  int e = 0;
  int rc = s21_mpf_get_str(buf, sizeof(buf), &e, 10, n_digits, x, S21_MPF_RNDN);
  if (rc != 0) {
    printf("(get_str error)\n");
    return;
  }

  if (x->kind != S21_MPF_NORMAL) {
    printf("%s\n", buf);
    return;
  }

  const char *p = buf;
  const char *sign = "";
  if (*p == '-') {
    sign = "-";
    p++;
  }
  const char *digits = p;

  if (e > 30 || e < -10) {
    int last = (int)n_digits - 1;
    while (last > 0 && digits[last] == '0') last--;
    printf("%s%c", sign, digits[0]);
    if (last > 0) {
      printf(".");
      fwrite(digits + 1, 1, (size_t)last, stdout);
    }
    printf("e%d\n", e - 1);
    return;
  }

  char out[1100];
  int pos = 0;
  if (*sign) out[pos++] = '-';

  if (e > 0) {
    uint32_t int_len = (uint32_t)e;
    if (int_len <= n_digits) {
      memcpy(out + pos, digits, int_len);
      pos += (int)int_len;
      out[pos++] = '.';
      memcpy(out + pos, digits + int_len, n_digits - int_len);
      pos += (int)(n_digits - int_len);
    } else {
      memcpy(out + pos, digits, n_digits);
      pos += (int)n_digits;
      for (uint32_t i = 0; i < int_len - n_digits; i++) out[pos++] = '0';
    }
  } else if (e == 0) {
    out[pos++] = '0';
    out[pos++] = '.';
    memcpy(out + pos, digits, n_digits);
    pos += (int)n_digits;
  } else {
    int zeros = -e;
    out[pos++] = '0';
    out[pos++] = '.';
    for (int i = 0; i < zeros; i++) out[pos++] = '0';
    memcpy(out + pos, digits, n_digits);
    pos += (int)n_digits;
  }
  out[pos] = '\0';

  char *dot = strchr(out, '.');
  if (dot) {
    int q = pos - 1;
    while (q > (int)(dot - out) && out[q] == '0') q--;
    if (q == (int)(dot - out)) q--;
    out[q + 1] = '\0';
  }

  printf("%s\n", out);
}

static void print_value(const s21_mpf_t *x) {
  if (g_hex_output)
    s21_mpf_print(x);
  else
    print_dec(x, g_dec_digits);
}

static void prompt(const char *text) {
  if (g_interactive) {
    fputs(text, stdout);
    fflush(stdout);
  }
}

static bool read_line(char *buf, size_t maxlen) {
  if (fgets(buf, (int)maxlen, stdin) == NULL) return false;
  buf[strcspn(buf, "\r\n")] = '\0';
  return true;
}

static bool is_number_token(const char *s) {
  if (s == NULL || *s == '\0') return false;
  const char *p = s;
  if (*p == '+' || *p == '-') p++;
  if (*p == '\0') return false;
  if (*p == '.') p++;
  return isdigit((unsigned char)*p) != 0;
}

static bool parse_uint(const char *s, uint32_t *out) {
  if (s == NULL || *s == '\0') return false;
  errno = 0;
  char *end = NULL;
  unsigned long v = strtoul(s, &end, 10);
  if (end == s || *end != '\0') return false;
  if (errno == ERANGE || v > UINT32_MAX) return false;
  *out = (uint32_t)v;
  return true;
}

static bool parse_prec(const char *s, uint32_t *out_bits) {
  if (s == NULL || *s == '\0') return false;

  size_t len = strlen(s);
  bool is_dec = false;
  if (s[len - 1] == 'd' || s[len - 1] == 'D') {
    is_dec = true;
    len--;
  }
  if (len == 0 || len >= 32) return false;

  char buf[32];
  memcpy(buf, s, len);
  buf[len] = '\0';

  errno = 0;
  char *end = NULL;
  unsigned long v = strtoul(buf, &end, 10);
  if (end == buf || *end != '\0') return false;
  if (errno == ERANGE || v > UINT32_MAX) return false;

  uint32_t bits;
  if (is_dec) {
    if (v == 0) return false;
    double b = (double)v * S21_LOG2_10;
    if (b > (double)UINT32_MAX - 2.0) return false;
    bits = (uint32_t)b;
    if ((double)bits < b) bits++;
    if (bits < 2) bits = 2;
  } else {
    if (v < 2) return false;
    bits = (uint32_t)v;
  }
  *out_bits = bits;
  return true;
}

static void print_prec_info(uint32_t bits) {
  uint32_t dec = (uint32_t)((double)bits * S21_LOG10_2);
  printf("prec = %u bits (~%u decimal digits)\n", bits, dec);
}

static bool ask_prec(const char *label, uint32_t *out) {
  char line[128];
  for (;;) {
    prompt(label);
    if (!read_line(line, sizeof(line))) return false;
    if (*line == '\0') return true;
    if (!parse_prec(line, out)) {
      fputs("(ожидается целое >= 2, или <n>d для десятичных цифр)\n", stderr);
      continue;
    }
    return true;
  }
}

static bool ask_arg(const char *label, char *buf, size_t maxlen) {
  for (;;) {
    prompt(label);
    if (!read_line(buf, maxlen)) return false;
    if (*buf == '\0') {
      fputs("(ожидается число)\n", stderr);
      continue;
    }
    return true;
  }
}

/* ------------------------------------------------------------
   Команды.
   ------------------------------------------------------------ */

static void cmd_help(void) {
  printf("Унарные:  ");
  for (size_t i = 0; i < unary_count; i++)
    printf("%s%s", unary_table[i].name, i + 1 < unary_count ? " " : "\n");
  printf("Бинарные: ");
  for (size_t i = 0; i < binary_count; i++)
    printf("%s%s", binary_table[i].name, i + 1 < binary_count ? " " : "\n");
  printf("Константы: pi [prec]\n");
  printf("Прочее:    setprec <n>[d], prec, show <num>, hex, dec [n],\n");
  printf("           quit/exit/q, help/?\n");
  printf("Формат:    <func> [arg1] [arg2] [prec]\n");
  printf("           <number>          == show <number>\n");
  printf("Точность:  число бит (например 256), либо <n>d — десятичных\n");
  printf("           цифр (например 77d == 256 бит).\n");
  printf("Аргументы: парсятся через s21_mpf_set_str (без double),\n");
  printf("           поддерживают '3.14', '1.5e-3', 'inf', 'nan'.\n");
  printf("Вывод:     dec      — десятичный через s21_mpf_get_str,\n");
  printf("                      по умолчанию 20 значащих цифр;\n");
  printf("           dec <n>  — установить число значащих цифр (1..1000);\n");
  printf("           hex      — полная мантисса в hex.\n");
  printf("π:         sin/cos/tan/atan/asin/acos кэшируют π на текущую\n");
  printf("           точность; при setprec кэш сбрасывается.\n");
  printf("Примеры:   sin 1.5 256\n");
  printf("           show 3.14159265358979323846\n");
  printf("           1e-30\n");
}

static void cmd_pi(const char *prec_tok) {
  uint32_t prec = g_default_prec;
  if (prec_tok != NULL) {
    if (!parse_prec(prec_tok, &prec)) {
      fprintf(stderr, "pi: некорректная точность '%s'\n", prec_tok);
      return;
    }
  }
  s21_mpf_t r;
  s21_mpf_init2(&r, prec);
  int rc = s21_mpf_pi(&r);
  if (rc == 0) print_value(&r);
  s21_mpf_clear(&r);
}

static void cmd_setprec(const char *tok) {
  uint32_t p;
  if (!parse_prec(tok, &p)) {
    fprintf(stderr, "setprec: ожидается целое >= 2, или <n>d\n");
    return;
  }
  g_default_prec = p;
  pi_cache_invalidate();
  print_prec_info(p);
}

static void cmd_show(const char *tok) {
  s21_mpf_t x;
  s21_mpf_init2(&x, g_default_prec);
  if (s21_mpf_set_str(&x, tok, 10) != 0) {
    fprintf(stderr, "show: не удалось разобрать '%s'\n", tok);
    s21_mpf_clear(&x);
    return;
  }
  print_value(&x);
  s21_mpf_clear(&x);
}

static void cmd_unary(const char *name, const char *arg_tok,
                      const char *prec_tok) {
  s21_unary_fn fn = find_unary(name);
  if (fn == NULL) {
    fprintf(stderr, "неизвестная унарная функция '%s'\n", name);
    return;
  }

  char arg_buf[256];
  const char *arg_str = arg_tok;
  if (arg_str == NULL) {
    char label[64];
    snprintf(label, sizeof(label), "arg (%s): ", name);
    if (!ask_arg(label, arg_buf, sizeof(arg_buf))) return;
    arg_str = arg_buf;
  }

  uint32_t prec = g_default_prec;
  if (prec_tok != NULL) {
    if (!parse_prec(prec_tok, &prec)) {
      fprintf(stderr, "%s: некорректная точность '%s'\n", name, prec_tok);
      return;
    }
  } else if (g_interactive && arg_tok == NULL) {
    char label[64];
    snprintf(label, sizeof(label), "prec [%u]: ", g_default_prec);
    if (!ask_prec(label, &prec)) return;
  }

  s21_mpf_t x, r;
  s21_mpf_init2(&x, prec);
  s21_mpf_init2(&r, prec);
  if (s21_mpf_set_str(&x, arg_str, 10) != 0) {
    fprintf(stderr, "%s: не удалось разобрать '%s'\n", name, arg_str);
    s21_mpf_clear(&x);
    s21_mpf_clear(&r);
    return;
  }
  int rc = fn(&r, &x);
  if (rc == 0) print_value(&r);
  else fprintf(stderr, "%s: ошибка\n", name);
  s21_mpf_clear(&x);
  s21_mpf_clear(&r);
}

static void cmd_binary(const char *name, const char *a_tok,
                       const char *b_tok, const char *prec_tok) {
  s21_binary_fn fn = find_binary(name);
  if (fn == NULL) {
    fprintf(stderr, "неизвестная бинарная функция '%s'\n", name);
    return;
  }

  char a_buf[256], b_buf[256];
  const char *a_str = a_tok;
  const char *b_str = b_tok;

  if (a_str == NULL) {
    char label[64];
    snprintf(label, sizeof(label), "arg1 (%s): ", name);
    if (!ask_arg(label, a_buf, sizeof(a_buf))) return;
    a_str = a_buf;
  }
  if (b_str == NULL) {
    char label[64];
    snprintf(label, sizeof(label), "arg2 (%s): ", name);
    if (!ask_arg(label, b_buf, sizeof(b_buf))) return;
    b_str = b_buf;
  }

  uint32_t prec = g_default_prec;
  if (prec_tok != NULL) {
    if (!parse_prec(prec_tok, &prec)) {
      fprintf(stderr, "%s: некорректная точность '%s'\n", name, prec_tok);
      return;
    }
  } else if (g_interactive && a_tok == NULL) {
    char label[64];
    snprintf(label, sizeof(label), "prec [%u]: ", g_default_prec);
    if (!ask_prec(label, &prec)) return;
  }

  s21_mpf_t x, y, r;
  s21_mpf_init2(&x, prec);
  s21_mpf_init2(&y, prec);
  s21_mpf_init2(&r, prec);

  if (s21_mpf_set_str(&x, a_str, 10) != 0) {
    fprintf(stderr, "%s: не удалось разобрать '%s'\n", name, a_str);
    goto cleanup;
  }
  if (s21_mpf_set_str(&y, b_str, 10) != 0) {
    fprintf(stderr, "%s: не удалось разобрать '%s'\n", name, b_str);
    goto cleanup;
  }

  {
    int rc = fn(&r, &x, &y);
    if (rc == 0) print_value(&r);
    else fprintf(stderr, "%s: ошибка\n", name);
  }

cleanup:
  s21_mpf_clear(&x);
  s21_mpf_clear(&y);
  s21_mpf_clear(&r);
}

static bool process_line(char *line) {
  const char *tok[4] = { NULL, NULL, NULL, NULL };
  size_t n = 0;
  for (char *p = strtok(line, " \t"); p != NULL && n < 4;
       p = strtok(NULL, " \t")) {
    tok[n++] = p;
  }
  if (n == 0) return true;

  const char *cmd = tok[0];

  if (strcmp(cmd, "quit") == 0 || strcmp(cmd, "exit") == 0 ||
      strcmp(cmd, "q") == 0) {
    return false;
  }

  if (strcmp(cmd, "help") == 0 || strcmp(cmd, "?") == 0) {
    cmd_help();
    return true;
  }

  if (strcmp(cmd, "prec") == 0) {
    print_prec_info(g_default_prec);
    return true;
  }

  if (strcmp(cmd, "setprec") == 0) {
    if (n < 2) {
      fprintf(stderr, "setprec: нужен аргумент\n");
      return true;
    }
    cmd_setprec(tok[1]);
    return true;
  }

  if (strcmp(cmd, "hex") == 0) {
    g_hex_output = true;
    printf("вывод: hex\n");
    return true;
  }

  if (strcmp(cmd, "dec") == 0) {
    if (n >= 2) {
      uint32_t d;
      if (!parse_uint(tok[1], &d) || d < 1 || d > 1000) {
        fprintf(stderr, "dec: ожидается число 1..1000\n");
        return true;
      }
      g_dec_digits = d;
    }
    g_hex_output = false;
    printf("вывод: dec (%u digits)\n", g_dec_digits);
    return true;
  }

  if (strcmp(cmd, "pi") == 0) {
    cmd_pi(n >= 2 ? tok[1] : NULL);
    return true;
  }

  if (strcmp(cmd, "show") == 0) {
    if (n < 2) {
      fprintf(stderr, "show: нужен аргумент\n");
      return true;
    }
    cmd_show(tok[1]);
    return true;
  }

  if (find_unary(cmd) != NULL) {
    cmd_unary(cmd, n >= 2 ? tok[1] : NULL, n >= 3 ? tok[2] : NULL);
    return true;
  }

  if (find_binary(cmd) != NULL) {
    cmd_binary(cmd, n >= 2 ? tok[1] : NULL, n >= 3 ? tok[2] : NULL,
               n >= 4 ? tok[3] : NULL);
    return true;
  }

  if (is_number_token(cmd)) {
    cmd_show(cmd);
    return true;
  }

  fprintf(stderr, "неизвестная команда '%s' (help — справка)\n", cmd);
  return true;
}

int main(void) {
  setlocale(LC_NUMERIC, "C");

  g_interactive = isatty(STDIN_FILENO) != 0;

  if (g_interactive) {
    printf("s21_mpfr REPL. Команды: help, quit. ");
    printf("Точность по умолчанию: %u бит. ", g_default_prec);
    printf("Вывод: dec (%u digits).\n", g_dec_digits);
  }

  char line[512];
  for (;;) {
    prompt("> ");
    if (!read_line(line, sizeof(line))) break;
    if (!process_line(line)) break;
  }

  pi_cache_invalidate();
  if (g_interactive) printf("\nдо встречи.\n");
  return 0;
}
