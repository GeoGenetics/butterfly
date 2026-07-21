CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c11
LDFLAGS = -lz
TARGET = lca_cmp
SRC = lca_cmp.c
TARGET2 = getrank_count
SRC2 = getrank_count.c
TARGET3 = bam_cmp
SRC3 = bam_cmp.c
.PHONY: all clean

all: $(TARGET) $(TARGET2) $(TARGET3)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

$(TARGET2): $(SRC2)
	$(CC) $(CFLAGS) -o $(TARGET2) $(SRC2) $(LDFLAGS)

$(TARGET3): $(SRC3)
	$(CC) $(CFLAGS) -o $(TARGET3) $(SRC3) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(TARGET2) $(TARGET3)
