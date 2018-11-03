all: hello

IDIR_GSL_HEADERS = libs/gsl-1.16
LIB_GSL=libs/gsl-1.16/.libs/libgsl.a

INCLUDES=-I libs/seq_file -I libs/string_buffer -I libs/cJSON -I $(IDIR_GSL_HEADERS)
LIBDIRS=-Llibs/string_buffer
CFLAGS=-Wall -Wextra
LINK=-lstrbuf -lz -lm -lpthread 
CC=gcc

SRCS=$(wildcard *.c) libs/cJSON/cJSON.c
HDRS=$(wildcard *.h)

ifdef VERBOSE
	CFLAGS := $(CFLAGS) -DVERBOSE=1
endif

ifdef DEBUG
	OPT = -O0 -Wstack-protector -fstack-protector
	DEBUG_ARGS = -g -ggdb3 -DDEBUG=1
else
	OPT = -O3 # -DNDEBUG=1
	DEBUG_ARGS = 
endif

hello: src/tools/hello.c $(SRCS) $(HDRS)
	$(CC) $(OPT) $(DEBUG_ARGS) $(CFLAGS) $(INCLUDES) $(LIBDIRS) -I . -o hello src/tools/hello.c $(SRCS) $(LIB_GSL) $(LINK)

clean:
	rm -rf hello *.dSYM *.greg

.PHONY: all clean
