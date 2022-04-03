#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <cassert>

#define FMT_HEADER_ONLY
#include "third-party/fmt-8.1.0/include/fmt/core.h"
#include "third-party/fmt-8.1.0/include/fmt/ranges.h"
#include "third-party/fmt-8.1.0/include/fmt/color.h"

/* macros */
#define ECHO_INPUT
#define STDERR_COLOR fg(fmt::terminal_color::bright_red),
#define GREEN fmt::fg(fmt::terminal_color::bright_green)
#define GRAY fmt::fg(fmt::terminal_color::bright_black)

/* using declarations */
using i64 = std::int64_t;
using u64 = std::uint64_t;
using chromosome_t = std::vector<bool>;

/* struct definitions */
struct generation_t
{
        double max_fitness() const;
        double avg_fitness() const;

        std::vector<chromosome_t> chromosomes;
        std::vector<i64> integer_reps;
        std::vector<double> domain_values;
        std::vector<double> fitness_values;
        std::vector<double> selection_probs;
        std::vector<double> interval_points;
        double fitness_sum;
};

/* global variables */
static FILE* fptr = std::fopen("report.txt", "w");
static int phase = 0;
static std::vector<std::pair<double, double>> max_avgs;
static std::ifstream finput;
static std::mt19937 rng(std::random_device{}());
static std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
#include "util.h"

/* from input */
static double a;
static double b;
static u64 chromosome_length;
static u64 population_size;
static std::array<double, 3> coefficients;
static i64 precision;
static double crossover_prob;
static double mutation_prob;
static u64 num_generations;

/* function declarations */
static void clear_screen();
static void read_input();
static i64 chromosome_to_integer(const chromosome_t&);
static double fitness(const double);
static chromosome_t generate_random_chromosome();
static double integer_to_domain(double);
static generation_t generation_from_chromosomes(generation_t&&);
static generation_t random_generation();
static generation_t op_selection(generation_t&&);
static generation_t op_crossover(generation_t&&, std::unordered_set<chromosome_t>&);
static generation_t op_mutations(generation_t&&);
static std::string chromosome_to_str(const chromosome_t&);
static std::string chromosome_to_str_diff(const chromosome_t&, const chromosome_t&);
static void print_normal(const generation_t&, const std::string_view, FILE*);
static void print_selected(const generation_t&, const generation_t&);
static void print_crossover(const generation_t&, const std::unordered_set<chromosome_t>&);
static void print_mutations(const generation_t&, const generation_t&);
static void menu_pause();
static void cprint(auto&&...);

/* function definitions */
double generation_t::max_fitness() const
{
        return *std::max_element(fitness_values.begin(), fitness_values.end());
}

double generation_t::avg_fitness() const
{
        return fitness_sum / (double)population_size;
}

void clear_screen()
{
        std::system("clear");
}

void read_input()
{
        read(a, b);
        read(population_size);
        read(coefficients[0], coefficients[1], coefficients[2]);
        read(precision, crossover_prob, mutation_prob, num_generations);

        const u64 n_unique_chromosomes = (b - a) * std::pow(10, precision);
        chromosome_length = std::ceil(std::log2(n_unique_chromosomes));
}

i64 chromosome_to_integer(const chromosome_t& c)
{
        i64 result = 0;
        for(std::size_t i = 0; i < chromosome_length; ++i)
        {
                assert(i < c.size());
                i64 mask = c[i];
                mask = mask << i;

                result |= mask;
        }

        return result;
}

double fitness(const double x)
{
        const double res = coefficients[2] * x * x + coefficients[1] * x + coefficients[0];
        return res;
}

chromosome_t generate_random_chromosome()
{
        static std::uniform_int_distribution<i64> bit_distrib(0, 1);

        chromosome_t res(chromosome_length);
        for(std::size_t i = 0; i < chromosome_length; ++i)
        {
                res[i] = bit_distrib(rng);
        }

        return res;
}

double integer_to_domain(const double x)
{
        const double dom = ((b - a) / (std::pow(2, chromosome_length) - 1)) * x + a;
        return dom;
}

generation_t generation_from_chromosomes(generation_t&& g)
{
        /* clear memory from previous generation */
        g.integer_reps.clear();
        g.domain_values.clear();
        g.fitness_values.clear();
        g.selection_probs.clear();
        g.interval_points.clear();
        g.fitness_sum = 0;

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

generation_t op_selection(generation_t&& g)
{
        std::vector<chromosome_t> chromosomes_after_selection;
        for(std::size_t i = 0; i < population_size; ++i)
        {
                const double pos = prob_distrib(rng);

                const auto it =
                    std::upper_bound(g.interval_points.begin(), g.interval_points.end(), pos);

                const std::size_t idx = it - g.interval_points.begin();
                chromosomes_after_selection.push_back(g.chromosomes[idx - 1]);

                cprint(fptr, fmt::runtime("{}\n"), std::make_pair(idx, pos));
        }

        g.chromosomes = std::move(chromosomes_after_selection);
        return generation_from_chromosomes(std::move(g));
}

generation_t op_crossover(generation_t&& g,
                          std::unordered_set<chromosome_t>& selected_for_cross)
{
        static std::uniform_int_distribution<u64> cut_distrib(0, (u64)chromosome_length);

        std::vector<chromosome_t> new_chromosomes;

        std::vector<std::size_t> idx_participants;
        for(std::size_t i = 0; i < population_size; ++i)
        {
                const double x = prob_distrib(rng);
                if(x < crossover_prob)
                {
                        selected_for_cross.insert(g.chromosomes[i]);
                        idx_participants.push_back(i);
                }
                else
                {
                        new_chromosomes.push_back(g.chromosomes[i]);
                }

                cprint(fptr, fmt::runtime("{:<3}: {}, u = {:>18.16f} {}\n"), i + 1,
                       chromosome_to_str(g.chromosomes[i]), x,
                       x < crossover_prob ? "part of crossover" : "not part of crossover");
        }

        cprint(fptr, fmt::runtime("\n"));

        if(idx_participants.size() % 2 == 1)
        {
                new_chromosomes.push_back(g.chromosomes[idx_participants.back()]);
                idx_participants.pop_back();
        }

        if(idx_participants.empty())
        {
                return g;
        }

        /* shuffle crossover pairs */
        std::shuffle(idx_participants.begin(), idx_participants.end(), rng);

        const i64 len = idx_participants.size();
        for(i64 k = 0; k < len - 1; k += 2)
        {
                const std::size_t i = idx_participants[k];
                const std::size_t j = idx_participants[k + 1];

                chromosome_t c1 = g.chromosomes[i];
                chromosome_t c2 = g.chromosomes[j];

                const u64 cut_point = cut_distrib(rng);

                cprint(fptr, fmt::runtime("Crossover between {} and {}:\n"), i + 1, j + 1);
                cprint(fptr, fmt::runtime("{} X {}, cut point: {}\n"), chromosome_to_str(c1),
                       chromosome_to_str(c2), cut_point);

                for(std::size_t pos = 0; pos < cut_point; ++pos)
                {
                        std::swap(c1[pos], c2[pos]);
                }

                cprint(fptr, fmt::runtime("Results: {}, {}\n"), chromosome_to_str(c1),
                       chromosome_to_str(c2));

                new_chromosomes.push_back(std::move(c1));
                new_chromosomes.push_back(std::move(c2));
        }

        g.chromosomes = std::move(new_chromosomes);

        return generation_from_chromosomes(std::move(g));
}

generation_t op_mutations(generation_t&& g)
{
        bool mutated = false;

        for(std::size_t j = 0; j < population_size; ++j)
        {
                auto& chromosome = g.chromosomes[j];

                for(std::size_t i = 0; i < chromosome_length; ++i)
                {
                        const double pos = prob_distrib(rng);
                        if(pos < mutation_prob)
                        {
                                cprint(fptr,
                                       fmt::runtime("Chromosome {}, gene {} was mutated\n"),
                                       j + 1, i + 1);

                                chromosome[i] = !chromosome[i];
                                mutated = true;
                        }
                }
        }

        if(!mutated)
        {
                return g;
        }

        return generation_from_chromosomes(std::move(g));
}

std::string chromosome_to_str(const chromosome_t& c)
{
        static constexpr std::array<char, 2> chars = {'0', '1'};

        std::string str;
        str.reserve(chromosome_length);

        for(auto it = c.rbegin(); it != c.rend(); ++it)
        {
                str.push_back(chars[*it]);
        }

        return str;
}

std::string chromosome_to_str_diff(const chromosome_t& c1, const chromosome_t& c2)
{
        const auto str1 = chromosome_to_str(c1);
        const auto str2 = chromosome_to_str(c2);

        std::string res;
        for(std::size_t i = 0; i < str1.size(); ++i)
        {
                if(str1[i] == str2[i])
                {
                        res += fmt::format(GRAY, "{}", str2[i]);
                }
                else
                {
                        res += fmt::format(GREEN, "{}", str2[i]);
                }
        }

        return res;
}

void print_normal(const generation_t& g, const std::string_view header = "",
                  FILE* fptr = stdout)
{
        fmt::print(fptr, "{}\n\n", header);

        for(std::size_t i = 0; i < population_size; ++i)
        {
                fmt::print(fptr, "{:<3}: {}, x = {:>10.6f}, f = {:>18.16f}\n", i + 1,
                           chromosome_to_str(g.chromosomes[i]), g.domain_values[i],
                           g.fitness_values[i]);
        }

        fmt::print(fptr, "\nMaximum fitness: {}, Average fitness: {}\n", g.max_fitness(),
                   g.fitness_sum / (double)population_size);
}

void print_selected(const generation_t& prev_g, const generation_t& g)
{
        fmt::print("SELECTED\n\n");

        std::unordered_set ch_after(g.chromosomes.begin(), g.chromosomes.end());
        for(std::size_t i = 0; i < population_size; ++i)
        {
                fmt::print(ch_after.contains(prev_g.chromosomes[i]) ? GREEN : GRAY,
                           "{:<3}: {}, x = {:>10.6f}, f = {:>18.16f}\n", i + 1,
                           chromosome_to_str(prev_g.chromosomes[i]), prev_g.domain_values[i],
                           prev_g.fitness_values[i]);
        }

        fmt::print("\nMaximum fitness: {}, Average fitness: {}\n", prev_g.max_fitness(),
                   prev_g.fitness_sum / (double)population_size);
}

void print_crossover(const generation_t& prev_g,
                     const std::unordered_set<chromosome_t>& picked)
{

        fmt::print("SELECTED FOR CROSSOVER\n\n");

        for(std::size_t i = 0; i < population_size; ++i)
        {
                fmt::print(picked.contains(prev_g.chromosomes[i]) ? GREEN : GRAY,
                           "{:<3}: {}, x = {:>10.6f}, f = {:>18.16f}\n", i + 1,
                           chromosome_to_str(prev_g.chromosomes[i]), prev_g.domain_values[i],
                           prev_g.fitness_values[i]);
        }

        fmt::print("\nMaximum fitness: {}, Average fitness: {}\n", prev_g.max_fitness(),
                   prev_g.fitness_sum / (double)population_size);
}

void print_mutations(const generation_t& prev_g, const generation_t& g)
{
        fmt::print("MUTATIONS (EACH GENE HAS p = {})\n\n", mutation_prob);

        std::unordered_set ch_after(g.chromosomes.begin(), g.chromosomes.end());
        for(std::size_t i = 0; i < population_size; ++i)
        {
                fmt::print("{:<3}: {}, x = {:>10.6f}, f = {:>18.16f}\n", i + 1,
                           chromosome_to_str_diff(prev_g.chromosomes[i], g.chromosomes[i]),
                           g.domain_values[i], g.fitness_values[i]);
        }

        fmt::print("\nMaximum fitness: {}, Average fitness: {}\n", g.max_fitness(),
                   g.fitness_sum / (double)population_size);
}

void menu_pause()
{
#ifdef INTERACTIVE
        fmt::print(stderr, "Press RETURN to go to next step!\n");
        std::cin.get();
        std::system("clear");
#endif
}

void cprint(auto&&... args)
{
        if(phase == 0)
        {
                fmt::print(std::forward<decltype(args)>(args)...);
        }
}

int main(const int argc, const char* argv[])
{
        if(argc != 2)
        {
                fmt::print(
                    stderr,
                    "usage:\ninteractive mode: make interactive && ./main "
                    "<input-file>\nnon-interactive "
                    "mode: make non-interactive && ./main <input-file> && cat report.txt\n");
                return EXIT_FAILURE;
        }

#ifndef INTERACTIVE
        std::freopen("/dev/null", "w", stdout);
#endif

        finput = std::ifstream(argv[1]);
        read_input();

#ifndef ECHO_INPUT
        fmt::print(stderr, STDERR_COLOR "Function domain interval: [{}, {})\n", a, b);
        fmt::print(stderr, STDERR_COLOR "Population size: {}\n", population_size);
        fmt::print(stderr, STDERR_COLOR "Coefficients: {}\n", fmt::join(coefficients, ", "));
        fmt::print(stderr, STDERR_COLOR "Precision: {}\n", precision);
        fmt::print(stderr, STDERR_COLOR "Crossover probability: {}\n", crossover_prob);
        fmt::print(stderr, STDERR_COLOR "Mutation probability: {}\n", mutation_prob);
        fmt::print(stderr, STDERR_COLOR "Number of generations: {}\n", num_generations);
#endif

        assert(fptr != nullptr);

        generation_t g = random_generation();

        print_normal(g, "Initial population:", fptr);

        fmt::print(fptr, "\nSelections probabilities:\n");
        for(std::size_t i = 0; i < population_size; ++i)
        {
                fmt::print(fptr, "{:<3}: p = {:>18.16f}\n", i + 1, g.selection_probs[i]);
        }
        fmt::print(fptr, "\n");

        fmt::print(fptr, "\nSelection intervals:\n{}\n", g.interval_points);

        for(u64 i = 0; i < num_generations; ++i)
        {
                clear_screen();
                print_normal(g, fmt::format("P({})", i));
                menu_pause();

                generation_t prev = g;

                /* select */
                cprint(fptr, fmt::runtime("\n"));
                cprint(fptr,
                       fmt::runtime("Selected elements (chromosome idx, probability)\n"));

                g = op_selection(std::move(g));
                if(phase == 0)
                {
                        print_normal(g, "\nAfter selection:", fptr);
                }

                cprint(fptr, fmt::runtime("\n"));

                print_selected(prev, g);
                menu_pause();

                /* after select */
                print_normal(g, "AFTER SELECTION");
                menu_pause();

                /* crossover */
                prev = g;
                std::unordered_set<chromosome_t> picked;

                cprint(fptr, fmt::runtime("Crossover probability: {}\n"), crossover_prob);
                g = op_crossover(std::move(g), picked);

                if(phase == 0)
                {
                        print_normal(g, "\nAfter crossover:", fptr);
                }

                print_crossover(prev, picked);
                menu_pause();

                /* after crossover */
                print_normal(g, "AFTER CROSSOVER");
                menu_pause();

                /* mutate */
                prev = g;
                g = op_mutations(std::move(g));

                print_mutations(prev, g);
                menu_pause();

                max_avgs.push_back(std::make_pair(g.max_fitness(), g.avg_fitness()));

                ++phase;
        }

        /* print last generation */
        clear_screen();
        print_normal(g, fmt::format("P({})", num_generations));

        fmt::print(fptr, "\nEvolution of (maximum fitness, average fitness):\n{}\n",
                   fmt::join(max_avgs, "\n"));

        fclose(fptr);
}
