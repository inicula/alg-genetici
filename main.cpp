#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <random>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/color.h>
#include "util.h"

#define ECHO_INPUT
#define STDERR_COLOR fg(fmt::terminal_color::bright_red),

/* using declarations */
using i64 = std::int64_t;
using u64 = std::uint64_t;
using chromosome_t = std::vector<bool>;
/* using declarations */

/* global variables */
static double a;
static double b;
static i64 num_unique_chromosomes;
static u64 chromosome_length;
static i64 population_size;
static std::array<double, 3> coefficients;
static i64 precision;
static double crossover_prob;
static double mutation_prob;
static i64 num_generations;
/* global variables */

/* function declarations */
static void read_input();
static i64 chromosome_to_integer(const chromosome_t&);
static double fitness(const double);
static chromosome_t generate_random_chromosome();
static double integer_to_domain(double);
/* function declarations */

/* function definitions */
void read_input()
{
        read(a, b);
        read(population_size);
        read(coefficients[0], coefficients[1], coefficients[2]);
        read(precision, crossover_prob, mutation_prob, num_generations);

        num_unique_chromosomes = (b - a) * std::pow(10, precision);
        chromosome_length = std::ceil(std::log2(num_unique_chromosomes));
}

i64 chromosome_to_integer(const chromosome_t& c)
{
        i64 result = 0;
        for(std::size_t i = 0; i < chromosome_length; ++i)
        {
                i64 mask = c[i];
                mask = mask << i;

                result |= mask;
        }

        return result;
}

double fitness(const double x)
{
        return coefficients[2] * x * x + coefficients[1] * x + coefficients[0];
}

chromosome_t generate_random_chromosome()
{
        static std::mt19937 rng(std::random_device{}());
        static std::uniform_int_distribution<i64> dist(0, num_unique_chromosomes);

        chromosome_t res;

        i64 x = dist(rng);
        for(std::size_t i = 0; i < chromosome_length; ++i)
        {
                i64 mask = 1l << i;

                res.push_back(x & mask);
        }

        return res;
}

double integer_to_domain(const double x)
{
        double total = num_unique_chromosomes;

        return (x / total) * (b - a) + a;
}
/* function definitions */

int main()
{
        read_input();

#ifdef ECHO_INPUT
        fmt::print(stderr, STDERR_COLOR "Function domain interval: [{}, {})\n", a, b);
        fmt::print(stderr, STDERR_COLOR "Population size: {}\n", population_size);
        fmt::print(stderr, STDERR_COLOR "Coefficients: {}\n", fmt::join(coefficients, ", "));
        fmt::print(stderr, STDERR_COLOR "Precision: {}\n", precision);
        fmt::print(stderr, STDERR_COLOR "Crossover probability: {}\n", crossover_prob);
        fmt::print(stderr, STDERR_COLOR "Mutation probability: {}\n", mutation_prob);
        fmt::print(stderr, STDERR_COLOR "Number of generations: {}\n", num_generations);
#endif

        const auto t = generate_random_chromosome();
        const auto i = chromosome_to_integer(t);
        const auto d = integer_to_domain(i);
        const auto f = fitness(d);
        fmt::print("{}\n{}\n{}\n{}\n", t, i, d, f);
}
