# ============================================================
# s21_mpfr — произвольная точность на C11
# ============================================================
# Цели:
#   make                     — статическая библиотека libs21_mpf.a
#   make test                — clean + сборка + запуск unit-тестов
#   make repl                — интерактивная утилита s21_mpf_repl
#   make bench               — микробенчмарк операций
#   make clang-format        — применить clang-format ко всем исходникам
#   make clang-format-check  — dry-run с --Werror (для CI)
#   make valgrind            — утечки памяти (memcheck)
#   make asan                — сборка с AddressSanitizer
#   make ubsan               — сборка с UndefinedBehaviorSanitizer
#   make gcov_report         — покрытие (lcov + genhtml)
#   make clean               — удалить артефакты
# ============================================================

CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11 -O2
LDFLAGS = -lm
CLANG_FORMAT ?= clang-format

SRC = s21_mpf_core.c \
      s21_mpf_bits.c \
      s21_mpf_cmp.c \
      s21_mpf_arith.c \
      s21_mpf_func.c \
      s21_mpf_print.c
OBJ = $(SRC:.c=.o)
LIB = libs21_mpf.a

TEST_SRC = test.c
TEST = test_mpf

REPL_SRC = repl.c
REPL = s21_mpf_repl

BENCH_SRC = bench.c
BENCH = s21_mpf_bench

HDRS = s21_mpf.h s21_mpf_internal.h
FORMAT_FILES = $(SRC) $(HDRS) $(TEST_SRC) $(REPL_SRC) $(BENCH_SRC)

# ============================================================
# Основное: библиотека
# ============================================================

all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^

%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) -I. -c $< -o $@

# ============================================================
# Тесты
# ============================================================

test:
	$(MAKE) clean
	$(MAKE) $(TEST)
	./$(TEST)

$(TEST): $(TEST_SRC) $(LIB)
	$(CC) $(CFLAGS) -I. $(TEST_SRC) $(LIB) $(LDFLAGS) -o $@

# ============================================================
# REPL
# ============================================================

repl: $(REPL)
	./$(REPL)

$(REPL): $(REPL_SRC) $(LIB)
	$(CC) $(CFLAGS) -I. $(REPL_SRC) $(LIB) $(LDFLAGS) -o $@

# ============================================================
# Бенчмарк
# ============================================================

bench: $(BENCH)
	./$(BENCH)

$(BENCH): $(BENCH_SRC) $(LIB)
	$(CC) $(CFLAGS) -I. $(BENCH_SRC) $(LIB) $(LDFLAGS) -o $@

# ============================================================
# Форматирование
# ============================================================

clang-format:
	$(CLANG_FORMAT) -i $(FORMAT_FILES)

clang-format-check:
	$(CLANG_FORMAT) --dry-run --Werror $(FORMAT_FILES)

# ============================================================
# Проверки памяти и неопределённого поведения
# ============================================================

valgrind: test
	valgrind --tool=memcheck --leak-check=yes --show-leak-kinds=definite ./$(TEST)

asan: clean
	$(CC) $(CFLAGS) -fsanitize=address -g -I. $(TEST_SRC) $(SRC) $(LDFLAGS) -o $(TEST)
	./$(TEST)

ubsan: clean
	$(CC) $(CFLAGS) -fsanitize=undefined -g -I. $(TEST_SRC) $(SRC) $(LDFLAGS) -o $(TEST)
	./$(TEST)

# ============================================================
# Покрытие
# ============================================================

gcov_report: clean
	$(CC) $(CFLAGS) --coverage -I. $(TEST_SRC) $(SRC) $(LDFLAGS) -o $(TEST)
	./$(TEST)
	lcov -t "s21_mpfr" -o s21_mpfr.info -c -d . --no-external
	genhtml -o report s21_mpfr.info

# ============================================================
# Очистка
# ============================================================

clean:
	rm -f $(OBJ) $(LIB) $(TEST) $(REPL) $(BENCH) *.gcda *.gcno *.info
	rm -rf report

.PHONY: all test repl bench clang-format clang-format-check valgrind asan ubsan \
        gcov_report clean
