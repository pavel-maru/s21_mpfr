CC = gcc
CFLAGS = -Wall -Wextra -Werror -std=c11 -O2
LDFLAGS = -lm

SRC = src/s21_mpf_core.c
OBJ = $(SRC:.c=.o)
LIB = libs21_mpf.a

TEST_SRC = src/test.c
TEST = test_mpf

all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^

%.o: %.c src/s21_mpf.h
	$(CC) $(CFLAGS) -Isrc -c $< -o $@

test: $(LIB)
	$(CC) $(CFLAGS) -Isrc $(TEST_SRC) $(LIB) $(LDFLAGS) -o $(TEST)
	./$(TEST)

clean:
	rm -f $(OBJ) $(LIB) $(TEST)

.PHONY: all clean test
