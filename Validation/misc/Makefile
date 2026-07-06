CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c11
LDFLAGS = -lz
TARGET = lca_cmp
SRC = lca_cmp.c
TARGET2 = getrank_count
SRC2 = getrank_count.c

.PHONY: all clean

all: $(TARGET) $(TARGET2)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

$(TARGET2): $(SRC2)
	$(CC) $(CFLAGS) -o $(TARGET2) $(SRC2) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(TARGET2)
