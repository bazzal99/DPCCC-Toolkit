CC      = gcc
CFLAGS  = -O2 -std=c99 -Wall -Wextra
LDFLAGS = -lm
TARGET  = dpccc_decoder_free
SRCS    = src/main.c src/designer.c src/rtz_search.c \
          ../../shared/encoder.c ../../shared/interleaver.c \
          ../../shared/utils.c

.PHONY: all clean run
all: $(TARGET)
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
run: all
	./$(TARGET)
clean:
	rm -f $(TARGET) result.txt temps.txt _dedup_tmp.txt
