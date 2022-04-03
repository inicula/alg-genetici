#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <random>
#include <algorithm>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/color.h>
#include "util.h"

/* macros */
#define ECHO_INPUT
#define STDERR_COLOR fg(fmt::terminal_color::bright_red),

/* using declarations */
using i64 = std::int64_t;
using u64 = std::uint64_t;
using chromosome_t = std::vector<bool>;

/* struct definitions */
struct generation_t
{
        std::vector<chromosome_t> chromosomes;
        std::vector<i64> integer_reps;
        std::vector<double> domain_values;
        std::vector<double> fitness_values;
        std::vector<double> selection_probs;
        std::vector<double> interval_points;
        double fitness_sum;
};

/* global variables */
static std::mt19937 rng(std::random_device{}());

/* from input */
static double a;
static double b;
static i64 num_unique_chromosomes;
static u64 chromosome_length;
static u64 population_size;
static std::array<double, 3> coefficients;
static i64 precision;
static double crossover_prob;
static double mutation_prob;
static u64 num_generations;

/* function declarations */
static void read_input();
static i64 chromosome_to_integer(const chromosome_t&);
static double fitness(const double);
static chromosome_t generate_random_chromosome();
static double integer_to_domain(double);
static generation_t generation_from_chromosomes(generation_t&&);
static generation_t op_select(generation_t);
static generation_t op_cross(generation_t);
static generation_t op_mutate(generation_t);

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

generation_t generation_from_chromosomes(generation_t&& g)
{
        /* clear memory from previous generation */
        g.fitness_sum = 0;
        [](auto&&... vecs)
        {
                (
                    [](auto& vec)
                    {
                            vec.clear();
                    }(vecs),
                    ...);
        }(g.integer_reps, g.domain_values, g.fitness_values, g.selection_probs,
          g.interval_points);

        for(std::size_t i = 0; i < population_size; ++i)
        {
                g.integer_reps.push_back(chromosome_to_integer(g.chromosomes[i]));
                g.domain_values.push_back(integer_to_domain(g.integer_reps[i]));
                g.fitness_values.push_back(fitness(g.domain_values[i]));
                g.fitness_sum += g.fitness_values[i];
        }

        g.interval_points.push_back(0);
        for(std::size_t i = 0; i < population_size; ++i)
        {
                g.selection_probs.push_back(g.fitness_values[i] / g.fitness_sum);
                g.interval_points.push_back(g.interval_points.back() + g.selection_probs[i]);
        }

        return g;
}

generation_t random_generation()
{
        generation_t g;
        g.chromosomes.resize(population_size);
        std::generate(g.chromosomes.begin(), g.chromosomes.end(), generate_random_chromosome);

        return generation_from_chromosomes(std::move(g));
}

generation_t op_select(generation_t g)
{
        generation_t next;
        return next;
}

generation_t op_cross(generation_t g)
{
        return g;
}

generation_t op_mutate(generation_t g)
{
        return g;
}

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

        generation_t g = random_generation();

        for(u64 i = 0; i < num_generations; ++i)
        {
                /* select */
                g = op_select(std::move(g));

                /* cross */
                g = op_cross(std::move(g));

                /* mutate */
                g = op_mutate(std::move(g));
        }
}
