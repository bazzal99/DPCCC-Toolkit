CC      = gcc
CFLAGS  = -O2 -std=c99 -Wall -Wextra
LDFLAGS = -lm -lpthread
TARGET  = dpccc_impulse
SRCS    = src/main.c src/designer.c \
          ../../shared/encoder.c ../../shared/interleaver.c \
          ../../shared/decoder.c ../../shared/utils.c

.PHONY: all clean run
all: $(TARGET)
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
run: all
	./$(TARGET)
clean:
	rm -f $(TARGET) result.txt temps.txt _dedup_tmp.txt repeates.txt
