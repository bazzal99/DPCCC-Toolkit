CC      = gcc
CFLAGS  = -O2 -std=c99 -Wall -Wextra
LDFLAGS = -lm
TARGET  = dpccc_sim
SRCS    = src/main.c \
          ../shared/encoder.c ../shared/interleaver.c \
          ../shared/decoder.c ../shared/utils.c

.PHONY: all clean run
all: $(TARGET)
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
run: all
	./$(TARGET)
clean:
	rm -f $(TARGET) result.txt
