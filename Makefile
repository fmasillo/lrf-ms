CC = g++
CFLAGS = -Wall -Wextra -O3 -std=c++20 -march=native -fopenmp

SRCS = matching_statistics.cpp libsais/src/libsais.c
OBJS = matching_statistics.o libsais/src/libsais.o

TARGET = lrf_ms

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f $(TARGET) $(OBJS)
