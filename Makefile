.DEFAULT_GOAL := interactive

include config.mk

SRC = main.cpp

clean:
	rm -f main

interactive:
	${CPPC} ${WFLAGS} ${CPPFLAGSINTERACTIVE} ${SRC} -DINTERACTIVE -o main

non-interactive:
	${CPPC} ${WFLAGS} ${CPPFLAGS} ${SRC} -o main

.PHONY: clean interactive non-interactive
