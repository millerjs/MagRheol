DEPS = $(OBJS:.o=.d)
CC = gcc
CFLAGS = -g3 -Wall -fpic -std=gnu99 -MMD -MP -DDEBUG -D_DEBUG -O3 -Werror
LDLIBS = -lm -lrt -pthread  -D_POSIX_C_SOURCE=200112L -D_GNU_SOURCE

##### Objects
BIN = libjosh.a
OBJS = libjosh.o threadpool.o cll.o

all: $(BIN)

$(BIN): $(OBJS)
	ar  rcs $(BIN) $(OBJS)

%.d: %.c

check-syntax:
	$(CC) -o nul -S ${CHK_SOURCES}

clean:
	rm -f $(OBJS) $(BIN) sample *.d
