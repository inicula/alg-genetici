#flags
CPPSTD = -std=c++20
WFLAGS = -Wall -Wextra -Wpedantic
CPPFLAGSINTERACTIVE = ${CPPSTD} -Os -march=native -flto -fno-exceptions -fno-rtti
CPPFLAGS = ${CPPSTD} -O3 -march=native -flto -fno-exceptions -fno-rtti

#compiler
CPPC = g++
