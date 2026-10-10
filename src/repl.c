/* ============================================================
   s21_mpfr — интерактивная утилита (REPL).

   Формат ввода:
     <func> [arg1] [arg2] [prec]

   Если arg1/arg2/prec пропущены — утилита спросит их отдельно.
   prec — точность в битах; по умолчанию 256.  Суффикс 'd' у
   числа означает десятичные цифры: 77d == 256 бит.

   Специальные команды:
     help, ?            — список функций
     quit, exit, q      — выход
     pi [prec]          — константа π
     prec               — показать текущую точность
     setprec <n>[d]     — изменить точность по умолчанию
     hex                — переключить вывод на полную мантиссу
     dec                — переключить вывод на double-приближение
   ============================================================ */

#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <locale.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "s21_mpf.h"

/* log2(10) и log10(2) с запасом точности, чтобы не тянуть
   их каждый раз через log(). */
#define S21_LOG2_10 3.3219280948873626
#define S21_LOG10_2 0.30102999566398120

/* ------------------------------------------------------------
   Таблицы функций.
   ------------------------------------------------------------ */

typedef int (*s21_unary_fn)(s21_mpf_t *, const s21_mpf_t *);
typedef int (*s21_binary_fn)(s21_mpf_t *, const s21_mpf_t *, const s21_mpf_t *);

/* Обёртки для функций, возвращающих void. */
static int wrap_abs(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_abs(res, x);
  return 0;
}

static int wrap_neg(s21_mpf_t *res, const s21_mpf_t *x) {
  s21_mpf_neg(res, x);
  return 0;
}

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
  { "sin",  s21_mpf_sin  },
  { "cos",  s21_mpf_cos  },
  { "tan",  s21_mpf_tan  },
  { "atan", s21_mpf_atan },
  { "asin", s21_mpf_asin },
  { "acos", s21_mpf_acos },
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

static void print_value(const s21_mpf_t *x) {
  if (g_hex_output)
    s21_mpf_print(x);      /* полная мантисса, p<exp> [prec=N] */
  else
    s21_mpf_print_d(x);    /* double-приближение, %.17g */
}

/* Печать приглашения — только в интерактивном режиме. */
static void prompt(const char *text) {
  if (g_interactive) {
    fputs(text, stdout);
    fflush(stdout);
  }
}

/* Прочитать строку из stdin. Возвращает false на EOF.
   Строка сохраняется без trailing newline. */
static bool read_line(char *buf, size_t maxlen) {
  if (fgets(buf, (int)maxlen, stdin) == NULL) return false;
  buf[strcspn(buf, "\r\n")] = '\0';
  return true;
}

/* Разбор double из строки. Возвращает true при успехе. */
static bool parse_double(const char *s, double *out) {
  if (s == NULL || *s == '\0') return false;
  errno = 0;
  char *end = NULL;
  double v = strtod(s, &end);
  if (end == s || *end != '\0') return false;
  if (errno == ERANGE && (v == 0.0 || v == 1.0 / 0.0 || v == -1.0 / 0.0))
    return false;
  *out = v;
  return true;
}

/* Разбор точности. Принимает:
     "256"   — 256 бит (минимум 2)
     "77d"   — 77 десятичных цифр, что даёт ceil(77·log2 10) = 256 бит
   Возвращает биты. */
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
    if ((double)bits < b) bits++; /* ceil */
    if (bits < 2) bits = 2;
  } else {
    if (v < 2) return false;
    bits = (uint32_t)v;
  }
  *out_bits = bits;
  return true;
}

/* Показать точность в битах и приблизительное число десятичных цифр. */
static void print_prec_info(uint32_t bits) {
  uint32_t dec = (uint32_t)((double)bits * S21_LOG10_2);
  printf("prec = %u bits (~%u decimal digits)\n", bits, dec);
}

/* ------------------------------------------------------------
   Ввод параметров с приглашением (используется, когда токенов
   в основной строке не хватает).
   ------------------------------------------------------------ */

static bool ask_double(const char *label, double *out) {
  char line[128];
  for (;;) {
    prompt(label);
    if (!read_line(line, sizeof(line))) return false;
    if (*line == '\0') {
      fputs("(ожидается число)\n", stderr);
      continue;
    }
    if (!parse_double(line, out)) {
      fputs("(не удалось разобрать число)\n", stderr);
      continue;
    }
    return true;
  }
}

static bool ask_prec(const char *label, uint32_t *out) {
  char line[128];
  for (;;) {
    prompt(label);
    if (!read_line(line, sizeof(line))) return false;
    if (*line == '\0') {
      /* пустой ввод — оставить значение по умолчанию */
      return true;
    }
    if (!parse_prec(line, out)) {
      fputs("(ожидается целое >= 2, или <n>d для десятичных цифр)\n", stderr);
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
  printf("Прочее:    setprec <n>[d], prec, hex, dec, quit/exit/q, help/?\n");
  printf("Формат:    <func> [arg1] [arg2] [prec]\n");
  printf("Точность:  число бит (например 256), либо <n>d — десятичных\n");
  printf("           цифр (например 77d == 256 бит).\n");
  printf("Пример:    sin 1.5 256\n");
  printf("           sin 1.5 100d\n");
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
  print_prec_info(p);
}

static void cmd_unary(const char *name, const char *arg_tok,
                      const char *prec_tok) {
  s21_unary_fn fn = find_unary(name);
  if (fn == NULL) {
    fprintf(stderr, "неизвестная унарная функция '%s'\n", name);
    return;
  }

  double d;
  if (arg_tok != NULL) {
    if (!parse_double(arg_tok, &d)) {
      fprintf(stderr, "%s: некорректный аргумент '%s'\n", name, arg_tok);
      return;
    }
  } else {
    char label[64];
    snprintf(label, sizeof(label), "arg (%s): ", name);
    if (!ask_double(label, &d)) return;
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
  s21_mpf_set_d(&x, d);
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

  double da, db;
  if (a_tok != NULL) {
    if (!parse_double(a_tok, &da)) {
      fprintf(stderr, "%s: некорректный первый аргумент '%s'\n", name, a_tok);
      return;
    }
  } else {
    char label[64];
    snprintf(label, sizeof(label), "arg1 (%s): ", name);
    if (!ask_double(label, &da)) return;
  }
  if (b_tok != NULL) {
    if (!parse_double(b_tok, &db)) {
      fprintf(stderr, "%s: некорректный второй аргумент '%s'\n", name, b_tok);
      return;
    }
  } else {
    char label[64];
    snprintf(label, sizeof(label), "arg2 (%s): ", name);
    if (!ask_double(label, &db)) return;
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
  s21_mpf_set_d(&x, da);
  s21_mpf_set_d(&y, db);
  int rc = fn(&r, &x, &y);
  if (rc == 0) print_value(&r);
  else fprintf(stderr, "%s: ошибка\n", name);
  s21_mpf_clear(&x);
  s21_mpf_clear(&y);
  s21_mpf_clear(&r);
}

/* ------------------------------------------------------------
   Разбор одной строки и диспетчеризация.
   Возвращает false, если пользователь попросил выйти.
   ------------------------------------------------------------ */

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
    g_hex_output = false;
    printf("вывод: dec\n");
    return true;
  }

  if (strcmp(cmd, "pi") == 0) {
    cmd_pi(n >= 2 ? tok[1] : NULL);
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

  fprintf(stderr, "неизвестная команда '%s' (help — справка)\n", cmd);
  return true;
}

/* ------------------------------------------------------------
   main.
   ------------------------------------------------------------ */

int main(void) {
  setlocale(LC_NUMERIC, "C");

  g_interactive = isatty(STDIN_FILENO) != 0;

  if (g_interactive) {
    printf("s21_mpfr REPL. Команды: help, quit. ");
    printf("Точность по умолчанию: %u бит.\n", g_default_prec);
  }

  char line[512];
  for (;;) {
    prompt("> ");
    if (!read_line(line, sizeof(line))) break;
    if (!process_line(line)) break;
  }

  if (g_interactive) printf("\nдо встречи.\n");
  return 0;
}
