CC ?= cc
CFLAGS ?= -O2 -Wall -Wextra -std=c11

pavages: pavages.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f pavages

.PHONY: clean
