#ifndef S21_MPF_H
#define S21_MPF_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/* Режимы округления (совместимы с MPFR) */
typedef enum {
  S21_MPF_RNDN = 0, /* к ближайшему, ties to even */
  S21_MPF_RNDZ,     /* к нулю */
  S21_MPF_RNDU,     /* вверх (+inf) */
  S21_MPF_RNDD      /* вниз (-inf) */
} s21_mpf_rnd_t;

/* Тип значения */
typedef enum {
  S21_MPF_ZERO = 0,
  S21_MPF_NORMAL,
  S21_MPF_INF,
  S21_MPF_NAN
} s21_mpf_kind_t;

/* Число произвольной точности.
   Инварианты для NORMAL:
     - limbs — little-endian массив, limbs[0] — младшие 64 бита
     - старший значащий бит находится на позиции (prec - 1)
     - x = (-1)^sign * mant * 2^(exp - prec), где mant = sum(limbs[i] *
   2^(64*i))
     - биты выше prec нулевые */
typedef struct {
  uint64_t *limbs; /* мантисса, little-endian */
  int64_t exp;     /* x = mant * 2^(exp - prec) */
  uint32_t prec;   /* точность в битах */
  int sign;        /* 0 = +, 1 = - */
  int kind;        /* ZERO / NORMAL / INF / NAN */
} s21_mpf_t;

/* Количество 64-битных лимбов для точности prec */
static inline size_t s21_mpf_limbs_for_prec(uint32_t prec) {
  return (prec + 63) / 64;
}

/* ============== Инициализация ============== */

void s21_mpf_init(s21_mpf_t *x); /* prec = 256 по умолчанию */
/* x не должен быть NULL.  prec зажимается до >= 2.  При неудаче
   calloc (практически только на prec в миллионах) x->limbs
   остаётся NULL, kind = ZERO — объект непригоден до clear и
   повторного init2.  Публичные функции эту ситуацию не
   проверяют: на реальных точностях calloc не отказывает. */
void s21_mpf_init2(s21_mpf_t *x, uint32_t prec);
void s21_mpf_clear(s21_mpf_t *x);

/* Смена точности на месте, сохраняя значение.
   Расширение — точно.  Сужение — округление по RNDN.
   0 — успех, -1 — ошибка (x == NULL или ошибка аллокации). */
int s21_mpf_set_prec(s21_mpf_t *x, uint32_t prec);

/* Смена точности на месте с указанным режимом округления.
   Расширение — точно при любом rnd.  Сужение — rnd.
   0 — успех, -1 — ошибка. */
int s21_mpf_set_prec_round(s21_mpf_t *x, uint32_t prec, s21_mpf_rnd_t rnd);

/* ============== Присваивание ============== */

void s21_mpf_set_zero(s21_mpf_t *x, int sign);
void s21_mpf_set_nan(s21_mpf_t *x);
void s21_mpf_set_inf(s21_mpf_t *x, int sign);
void s21_mpf_set_ui(s21_mpf_t *x, unsigned long v);
void s21_mpf_set_si(s21_mpf_t *x, long v);
void s21_mpf_set_d(s21_mpf_t *x, double v);

/* Присваивание по умолчанию (RNDN при смене точности). */
void s21_mpf_set(s21_mpf_t *dst, const s21_mpf_t *src);

/* Присваивание с указанным режимом округления при смене точности.
   Всегда возвращает 0 (или -1 при NULL). */
int s21_mpf_set_round(s21_mpf_t *dst, const s21_mpf_t *src, s21_mpf_rnd_t rnd);

/* ============== Утилиты ============== */

int s21_mpf_is_nan(const s21_mpf_t *x);
int s21_mpf_is_inf(const s21_mpf_t *x);
int s21_mpf_is_zero(const s21_mpf_t *x);
int s21_mpf_is_normal(const s21_mpf_t *x);
int s21_mpf_sign(const s21_mpf_t *x);
uint32_t s21_mpf_get_prec(const s21_mpf_t *x);

/* Низкоуровневые операции над битами мантиссы.
   Сдвиги объявлены во внутреннем заголовке s21_mpf_internal.h. */
int s21_mpf_get_bit(const s21_mpf_t *x, uint32_t pos);
int s21_mpf_msb(const s21_mpf_t *x);
int s21_mpf_lsb(const s21_mpf_t *x);
void s21_mpf_set_bit(s21_mpf_t *x, uint32_t pos, int value);

/* ============== Сравнение ============== */

/* Возвращает -1 / 0 / +1. Для NaN — 0 (как в MPFR). */
int s21_mpf_cmp(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_cmp_abs(const s21_mpf_t *x, const s21_mpf_t *y);
/* Побитовое сравнение представления: возвращает 1, только когда
   совпадают kind, sign, exp, prec и мантисса.  Для одного и
   того же значения при разных точностях — 0:
   equal(1@64, 1@128) == 0, хотя cmp(1@64, 1@128) == 0.
   Для числового равенства используйте s21_mpf_cmp(x, y) == 0. */
int s21_mpf_equal(const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_zero_p(const s21_mpf_t *x);
int s21_mpf_integer_p(const s21_mpf_t *x);

/* ============== Округление к целому ============== */

/* Усечение к нулю: отбрасывает дробную часть.
   trunc(±0) = ±0, trunc(±inf) = ±inf, trunc(NaN) = NaN. */
int s21_mpf_trunc(s21_mpf_t *res, const s21_mpf_t *x);

/* Округление к -inf.  floor(±0) = ±0, floor(±inf) = ±inf. */
int s21_mpf_floor(s21_mpf_t *res, const s21_mpf_t *x);

/* Округление к +inf.  ceil(±0) = ±0, ceil(±inf) = ±inf. */
int s21_mpf_ceil(s21_mpf_t *res, const s21_mpf_t *x);

/* Округление к ближайшему целому, ties away from zero (C99 round).
   round(±0) = ±0, round(±inf) = ±inf, round(NaN) = NaN. */
int s21_mpf_round(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== Остаток ============== */

/* fmod(x, y) = x - trunc(x/y) * y.  Знак результата = знак x.
   Спецзначения по C99:
     fmod(x, ±inf) = x       (для конечного x)
     fmod(±inf, y) = NaN
     fmod(x, 0)    = NaN
     fmod(±0, y)   = ±0      (знак x)
   При |x| >> |y| добавляется guard-точность ~ log2(|x/y|);
   выше S21_MPF_MAX_ARG_BITS — NaN. */
int s21_mpf_fmod(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

/* ============== Арифметика ============== */

int s21_mpf_add(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_sub(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_mul(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);
int s21_mpf_div(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

/* ============== Fused multiply-add ============== */

/* res = x·y + z с одним округлением в конце (RNDN).
   Промежуточное произведение x·y не округляется — это
   отличает fma от mul + add.  Классический пример: при
   x = 1 + 2^-50, y = 1 - 2^-50, z = -1 и prec = 100 fma
   даёт -2^-100, тогда как mul+add возвращает 0 (mul теряет
   младшие биты при округлении до prec).

   Спецзначения — по C99 / IEEE 754-2008 F.10.10.1:
     fma(NaN, y, z)          = NaN
     fma(±inf, 0, z)         = NaN
     fma(±inf, y≠0, z)       = ±inf (знак = sign(inf)·sign(y))
     fma(±inf, y≠0, z=∓inf)  = NaN (конфликт знаков)
     fma(x, y, ±inf)         = ±inf
     fma(±0, y, z)           = z  (±0 сохраняется по RNDN,
                                    если z = ∓0)

   0 при успехе, -1 при NULL. */
int s21_mpf_fma(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y,
                const s21_mpf_t *z);

void s21_mpf_neg(s21_mpf_t *res, const s21_mpf_t *x);
void s21_mpf_abs(s21_mpf_t *res, const s21_mpf_t *x);
void s21_mpf_normalize(s21_mpf_t *x);

/* ============== Возведение в степень ============== */

/* x^y.  Семантика — C99 pow / MPFR (Annex F):
     pow(any, ±0)         = 1     (включая pow(NaN, 0))
     pow(+1, any)         = 1     (включая pow(1, NaN), pow(1, ±inf))
     pow(-1, ±inf)        = 1
     pow(x<0, y нецелое)  = NaN
     pow(±0, y>0)         = ±0    знак = sign_x & (y mod 2)
     pow(±0, y<0)         = ±inf  то же
     pow(±inf, y>0)       = ±inf  то же
     pow(±inf, y<0)       = ±0    то же
     pow(|x|>1, ±inf)     = ±inf / ±0
     pow(|x|<1, ±inf)     = ±0 / ±inf

   Общий случай: x > 0 → exp(y·log x).
   x < 0, y целое → exp(y·log|x|), знак = sign_x & (y mod 2).

   0 при успехе, -1 при NULL. */
int s21_mpf_pow(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

/* Целочисленный показатель: x^n, n ≥ 0.
   Бинарное возведение в степень, без exp/log — точнее и быстрее
   общего pow.  Спецзначения — по общим правилам pow. */
int s21_mpf_pow_ui(s21_mpf_t *res, const s21_mpf_t *x, unsigned long n);

/* Целочисленный знаковый показатель.  Для n < 0 —
   1 / pow_ui(x, -n) с корректной обработкой x = ±0 / ±inf. */
int s21_mpf_pow_si(s21_mpf_t *res, const s21_mpf_t *x, long n);

/* Корень n-й степени: x^(1/n), n ≥ 1.
     root(x, 1) = x
     root(x, 2) = sqrt(x)
     root(x, 0) = NaN
     root(x<0, чётное n) = NaN
     root(±0, нечётное n) = ±0
   Общий случай: exp(log|x| / n), знак отдельно. */
int s21_mpf_root(s21_mpf_t *res, const s21_mpf_t *x, unsigned long n);

/* ============== Элементарные функции ============== */

int s21_mpf_sqrt(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_exp(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_log(s21_mpf_t *res, const s21_mpf_t *x);

/* log(1 + x) с повышенной точностью при малых |x|.
   Для малых |x| добавка guard-битов, пропорциональная -x->exp,
   компенсирует cancelation в 1 + x.  Для |x| < 2^(-prec/2)
   результат равен x с точностью выше prec.

   Домен x > -1:  log1p(-1) = -inf, log1p(x<-1) = NaN,
   log1p(+inf) = +inf, log1p(-inf) = NaN. */
int s21_mpf_log1p(s21_mpf_t *res, const s21_mpf_t *x);

/* exp(x) - 1 с повышенной точностью при малых |x|.
   Guard-биты компенсируют cancelation в exp(x) - 1.
   Для |x| < 2^(-prec/2) результат равен x.

   expm1(-inf) = -1, expm1(+inf) = +inf. */
int s21_mpf_expm1(s21_mpf_t *res, const s21_mpf_t *x);

/* Кубический корень через метод Ньютона с seed из double.
   Определён для любого знака: cbrt(-8) = -2.
   cbrt(±0) = ±0, cbrt(±inf) = ±inf. */
int s21_mpf_cbrt(s21_mpf_t *res, const s21_mpf_t *x);

/* log2(x) через reduction: x = 2^k · y, y ∈ [1, 2),
   log2(x) = k + log2(y).  Для x = 2^m результат точен
   (k = m, log2(1) = 0). */
int s21_mpf_log2(s21_mpf_t *res, const s21_mpf_t *x);

/* log10(x) = log(x) / ln(10).
   10^k не представимо точно в бинарном mpf, поэтому
   reduction 10^k не даёт выигрыша против прямого деления. */
int s21_mpf_log10(s21_mpf_t *res, const s21_mpf_t *x);

/* hypot(x, y) = √(x² + y²) со масштабированием по max(|x|, |y|).
   Без масштабирования x² + y² теряет точность при разнице
   порядков и переполняется при |x| ~ 2^(prec/2).  Спецзначения
   по C99: hypot(±inf, NaN) = +inf, hypot(NaN, non-inf) = NaN. */
int s21_mpf_hypot(s21_mpf_t *res, const s21_mpf_t *x, const s21_mpf_t *y);

/* ============== Константы ============== */

int s21_mpf_pi(s21_mpf_t *res);
int s21_mpf_ln2(s21_mpf_t *res);
int s21_mpf_e(s21_mpf_t *res);

/* ============== Тригонометрия: прямые ============== */

int s21_mpf_sin(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_cos(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_tan(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== Тригонометрия: обратные ============== */

int s21_mpf_atan(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_asin(s21_mpf_t *res, const s21_mpf_t *x);
int s21_mpf_acos(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== Тригонометрия с явной π ==============

   Варианты sin / cos / tan / atan / asin / acos, принимающие
   уже вычисленную π параметром.  Позволяет не пересчитывать π
   на каждый вызов — при prec = 4096 только на π уходит ~10 мс
   (формула Мэчина, ~1200 членов ряда).

   Типичное использование:

     s21_mpf_t pi;
     s21_mpf_init2(&pi, prec + 128);
     s21_mpf_pi(&pi);

     for (...) {
       s21_mpf_sin_with_pi(&r, &x, &pi);   // без пересчёта π
     }

   Контракт: pi->prec >= res->prec + 128 + E, где E = max(0, exp-1)
   при x->exp = exp.  Для |x| < 2 требование сводится к прежнему
   +128; с ростом |x| растёт, потому что аргумент-редукция
   reduce_2pi теряет ~log2(|x|) бит точности.  Меньшая точность π
   даёт -1.  При E > S21_MPF_MAX_ARG_BITS функция возвращает NaN
   (аргумент слишком велик для данного уровня точности).

   Спецзначения x (NaN / ±inf / ±0) обрабатываются как в
   одноимённых функциях без π.

   0 при успехе, -1 при ошибке (NULL, точность π мала). */
int s21_mpf_sin_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                        const s21_mpf_t *pi);
int s21_mpf_cos_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                        const s21_mpf_t *pi);
int s21_mpf_tan_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                        const s21_mpf_t *pi);
int s21_mpf_atan_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                         const s21_mpf_t *pi);
int s21_mpf_asin_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                         const s21_mpf_t *pi);
int s21_mpf_acos_with_pi(s21_mpf_t *res, const s21_mpf_t *x,
                         const s21_mpf_t *pi);

/* ============== Гиперболические функции ============== */

/* sinh(x) = (e^x − e^−x) / 2.  Нечётная.  sinh(±inf) = ±inf. */
int s21_mpf_sinh(s21_mpf_t *res, const s21_mpf_t *x);

/* cosh(x) = (e^x + e^−x) / 2.  Чётная, cosh(x) ≥ 1.  cosh(±inf) = +inf. */
int s21_mpf_cosh(s21_mpf_t *res, const s21_mpf_t *x);

/* tanh(x) = sinh(x) / cosh(x).  Нечётная, |tanh| < 1.  tanh(±inf) = ±1. */
int s21_mpf_tanh(s21_mpf_t *res, const s21_mpf_t *x);

/* asinh(x) = log(x + √(x²+1)).  Нечётная, определена везде. */
int s21_mpf_asinh(s21_mpf_t *res, const s21_mpf_t *x);

/* acosh(x) = log(x + √((x−1)(x+1))).  Домен x ≥ 1.
   acosh(1) = 0, acosh(x<1) = NaN, acosh(+inf) = +inf. */
int s21_mpf_acosh(s21_mpf_t *res, const s21_mpf_t *x);

/* atanh(x) = ½·log((1+x)/(1−x)).  Домен |x| < 1.
   atanh(±1) = ±inf, atanh(|x|>1) = NaN. */
int s21_mpf_atanh(s21_mpf_t *res, const s21_mpf_t *x);

/* ============== I/O ============== */

/* Разбор строки в mpf.  base ∈ [2, 36], алфавит цифр 0-9a-z
   (регистронезависимо на входе).  Спецзначения: inf, infinity,
   nan (регистронезависимо).

   Формат:
     base 10:  [+-]? digits [. digits]?  [eE] [+-]? dec_digits
     base 16:  [+-]? [0x] hex_digits [. hex_digits]?
                     [pP] [+-]? dec_digits
     прочие:   [+-]? digits [. digits]?

   Для base 16 поддерживается C99 hexfloat: префикс 0x / 0X
   (опциональный), p / P-экспонента — десятичное число, задающее
   степень двойки.  Значение = mant × 2^(p_exp).  Без p —
   обычный hex без сдвига экспоненты.

   Ведущие пробелы (space, tab, \n, \r) пропускаются.  Хвостовые
   не допускаются: set_str(&x, "42 ", 10) возвращает -1.  Это
   намеренно — строгий разбор ловит случайные лишние символы,
   но легко споткнуться при чтении из файлов.

   0 при успехе, -1 при ошибке (NULL, base вне [2, 36],
   неразобранный хвост, экспонента без цифр). */
int s21_mpf_set_str(s21_mpf_t *x, const char *str, int base);

/* Конвертирует x в строку из n_digits значащих цифр по основанию
   base ∈ [2, 36] (цифры 0-9a-z строчные) с округлением по rnd
   и записывает в buf.

   Формат: знак (если отрицательное), затем n_digits цифр без
   разделителя.  Позиция точки относительно начала строки цифр
   возвращается в *exp_out: value ≈ 0.DDDD... × base^exp_out
   (конвенция MPFR).

   Спецзначения: "nan", "inf", "-inf", "0", "-0"; для них
   *exp_out = 0 (base-агностично).

   n_digits == 0 → значение по умолчанию (20).

   0 при успехе, -1 при ошибке (NULL, base вне [2, 36], буфер мал). */
int s21_mpf_get_str(char *buf, size_t bufsize, int *exp_out, int base,
                    uint32_t n_digits, const s21_mpf_t *x, s21_mpf_rnd_t rnd);

/* ============== Печать ============== */

/* Вывод в произвольный поток.  Формат идентичен print / print_d;
   fprint(stdout, x) эквивалентен print(x).  NULL stream или x —
   no-op. */
void s21_mpf_fprint(FILE *stream, const s21_mpf_t *x);
void s21_mpf_fprint_d(FILE *stream, const s21_mpf_t *x);

/* Тонкие обёртки: print(x) == fprint(stdout, x);
   print_d(x) == fprint_d(stdout, x). */
void s21_mpf_print(const s21_mpf_t *x);
void s21_mpf_print_d(const s21_mpf_t *x);

#endif
