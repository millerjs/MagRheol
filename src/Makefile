DEPS = $(OBJS:.o=.d)
BIN = ./magrheol

CC = gcc

CFLAGS = -ggdb -g -Wall -fpic  -MMD -MP -DDEBUG -O3 -lm -std=c99 -Werror
LDLIBS = -lm -lrt -pthread -Llibjosh/ -ljosh

OBJS = maghreol.o \
	particles.o \
	domain.o

all: $(BIN)

$(BIN): $(OBJS)
	$(MAKE) -C libjosh/	
	$(CC) $(DFLAGS) $(OBJS) $(LIB) $(LDLIBS) -o ../$(BIN)

%.d: %.c

check-syntax:
	$(CC) -o nul -S ${CHK_SOURCES}

clean:
	$(MAKE) clean -C libjosh/ 
	rm -f $(OBJS) ../$(BIN) *.d
