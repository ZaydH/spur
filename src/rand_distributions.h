/**
 *  rand_distributions.h
 *
 *  Purpose: Static class for getting randomly created objects
 *  from different distributions and of different types.
 *
 *  @author Zayd Hammoudeh <zayd@ucsc.edu>
 *  @version 0.00.00
 *
 *  Copyright (C) 2018 Zayd Hammoudeh.
 *  All rights reserved.
 *
 * This software may be modified and distributed under the terms of the MIT license.  See the
 * LICENSE file for details.
 */

#ifndef PROJECT_RANDOM_H
#define PROJECT_RANDOM_H

#include <gmpxx.h>
#include <limits.h>

#include <chrono> // NOLINT (build/c++11)
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include "solver_config.h"
#include "primitive_types.h"
#include "sampler_tools.h"


class Random{
 public:
  /**
   * Random State Initializer
   *
   * Creates the random state of the program.  It also seeds the random
   * number generator.
   *
   * @param config Solver configuration
   */
  static void init(SolverConfiguration * config) {
    master_config_ = config;
    gmp_randinit_mt(Random::Mpz::rand_state_);

    // create the generator seed for the random engine to reference
    long long int seed;
    if (master_config_->debug_mode) {
      if (!master_config_->quiet)
        std::cout << "WARNING: Debug mode has a fixed seed for the MPZ class.\n";
      seed = 0;
    } else {
      seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
    mpz_class rand_seed;
    mpz_import(rand_seed.get_mpz_t(), 1, -1, sizeof(seed), 0, 0, &seed);
    gmp_randseed(Random::Mpz::rand_state_, rand_seed.get_mpz_t());

    if (master_config_->debug_mode) {
      if (!master_config_->quiet)
        std::cout << "WARNING: Debug mode uses a fixed random number seed." << std::endl;
    } else {
      Random::rng_ = std::mt19937(Random::rd_());
    }
  }
  /**
   * Generates a uniform random integer in the range [min, max].
   *
   * @param min Minimum value that can be generated inclusive.
   * @param max Maximum value that can be generated inclusively.
   * @return Uniform random integer across the full range of possible values.
   */
  //@{
  /**
   * @brief
   */
  template <typename T> inline static T uniform(T min, T max) {
    std::uniform_int_distribution<T> distribution(min, max);
    return distribution(rng_);
  }
  /**
   * @brief Special int uniform random generator that can optionally avoid reinitialize
   * the uniform int distribution generator.
   */
  inline static int uniform(int min = INT_MIN, int max = INT_MAX) {
    if (min == INT_MIN && max == INT_MAX) {
      return uni_(rng_);
    } else {
      std::uniform_int_distribution<int> distribution(min, max);
      return distribution(rng_);
    }
  }
  //@}
  /**
   * Creates a random variable sampled from the distribution Binom(n,p)
   *
   * @param n Number of Bernoulli trials in the binomial distribution
   * @param p Probability of success.
   *
   * @return Integer random variable generated according to the binomial
   * distribution.
   */
  static SampleSize binom(SampleSize n, double p) {
    std::binomial_distribution<> d(n, p);
    return static_cast<SampleSize>(d(rd_));
  }
  /**
   * Using a modified version of the Fisher-Yates shuffle, the first \p target_size elements of the
   * vector \p oversampled_vec will be selected uniformly at random.  All modifications are done
   * inplace.  The running time is O(\p target_size).
   *
   * @tparam ListType Type of the object in the oversampled vector.
   * @tparam CountType Type of object used to store the counts.
   *
   * @param oversampled_vec Vector with more elements than desired.  \p target_size elements will
   * be selected in the first [0, \p target_size) elements.
   * @param target_size Number of objects to select in \p oversampled_vec
   * @param resize True to resize \p oversampled_vec to the size \p target_count.
   */
  template<typename ListType, typename CountType>
  static void DownsampleList(CountType target_size, std::vector<ListType> &oversampled_vec,
                             bool resize = true) {
    assert(oversampled_vec.size() >= target_size);

    CountType end_point = oversampled_vec.size() - 1;
    while (end_point > target_size) {
      auto id_loc = Random::uniform<CountType>(0, end_point);
      oversampled_vec[id_loc] = oversampled_vec[end_point];  // Copy back & overwrite the used value
      end_point--;
    }
    if (resize)
      oversampled_vec.resize(target_size);
  }
  /**
   * Shuffles the specified vector in place uniformly at random.
   *
   * @tparam T Vector element type
   *
   * @param vec Vector to be shuffled.
   */
  template<typename T>
  static void shuffle(std::vector<T> &vec) {
    std::shuffle(vec.begin(), vec.end(), rd_);
  }

  class Mpz {
    // Allow the Random class to access this class' private methods/fields
    friend class Random;
   public:
    /**
     * Generates and returns a random multiprecision integer.  This may be used
     * by other classes to consolidate the random MP generation.
     *
     * @param max_z Maximum integer value.  All generated random values will be
     * in the range [0, max_z).
     *
     * @param rand_val Output value where the generated random number will
     * be stored.
     */
    inline static void uniform(mpz_class max_z, mpz_class &rand_val) {
      mpz_urandomm(rand_val.get_mpz_t(), Random::Mpz::rand_state_, max_z.get_mpz_t());
    }
    /**
     * Selects a binomially distributed random variable from Binom(n, a/t).
     *
     * It may an approximate or exact version of the distribution depending
     * on the state of the variable Random::Mpz::use_approx_binom_.
     *
     * @param n Number of Bernoulli trials in the binomial distribution
     * @param t Total weight of all samples
     * @param a Weight of success.
     *
     * @return Integer random variable generated according to the binomial
     * distribution.
     */
    inline static SampleSize binom(SampleSize n, const mpz_class &t,
                                   const mpz_class &a) {
      assert(t > 0 && t >= a);
      // Handle the edge cases
      if (a == 0)
        return 0;
      else if (a == t)
        return n;

      if (use_approx_binom_)
        return Random::Mpz::binom_approx(n, t, a);
      else
        return Random::Mpz::binom_exact(n, t, a);
    }

   private:
    /**
     * Use the approximate binomial distribution.  This will be faster but
     * less accurate than the exact method.
     */
    static bool use_approx_binom_;
    /**
     * Used for generating MPZ random numbers
     */
    static gmp_randstate_t rand_state_;
    /**
     * Creates a random variable sampled from the distribution
     * Binom(\p n, \p a/ \p t).  It is labeled as "approximate" since it
     * converts the MPZ objects to doubles and may have floating point
     * errors or underflow.
     *
     * It is expected to be FASTER than the Random::Mpz::binom_exact.
     *
     * @param n Number of Bernoulli trials in the binomial distribution
     * @param t Total weight of all samples
     * @param a Weight of success.
     *
     * @return Integer random variable generated according to the binomial
     * distribution.
     */
    inline static SampleSize binom_approx(SampleSize n, const mpz_class &t,
                                          const mpz_class &a) {
      mpf_class a_mpf = a, t_mpf = t;
      mpf_class p_mpf = a_mpf / t_mpf;
      return Random::binom(n, p_mpf.get_d());
//      double p = mpz_get_d(a.get_mpz_t()) / mpz_get_d(t.get_mpz_t());
//      return Random::binom(n,p);
    }
    /**
     * Creates a random variable sampled from the distribution Binom(n, a/t).
     * It is labeled as "exact" since it does not rely on floating point
     * conversion.  It uses MPZ methods only.
     *
     * It is expected to be SLOWER than the Random::Mpz::binom_approx.
     *
     * @param n Number of Bernoulli trials in the binomial distribution
     * @param t Total weight of all samples
     * @param a Weight of success.
     *
     * @return Integer random variable generated according to the binomial
     * distribution.
     */
    inline static SampleSize binom_exact(SampleSize n, const mpz_class &t,
                                         const mpz_class &a) {
      mpz_class rand_mpz = 0;
      SampleSize num_success = 0;
      for (SampleSize i = 0; i < n; i++) {
        Random::Mpz::uniform(t, rand_mpz);
        if (rand_mpz < a)
          num_success++;
      }
      return num_success;
    }
  };
  /**
   * Builds a randomly selected list of integers in the range [0, \p max_id).
   *
   * This is performed WITHOUT replacement.  Its running time is Theta(\p max_id)
   *
   * @param max_id Maximum value (exclusive) that any integer in the range can be.
   *
   * @param num_elements Number of elements to be selected
   * @param samples_to_replace List of p num_elements randomly selected integers in the range
   * [0, \p max_id).
   */
  static void SelectRangeInts(SampleSize max_id, SampleSize num_elements,
                              std::vector<SampleSize> &samples_to_replace) {
    assert(max_id >= num_elements);
    samples_to_replace.clear();
    samples_to_replace.reserve(max_id);
    for (SampleSize i = 0; i < max_id; i++)
      samples_to_replace.emplace_back(i);
    //ToDo Modify this function to only ever do a maximum of max_id/2 random numbers
    Random::DownsampleList<SampleSize, SampleSize>(num_elements, samples_to_replace);
  }
  /**
   * Shuffles a vector using the standad random device.
   *
   * @tparam T Type contained in the vector
   * @param begin Beginning of the vector to shuffle  (inclusive)
   * @param end End of the vector to shuffle (exclusive)
   */
  template<typename T>
  static void shuffle(typename std::vector<T>::iterator begin, // Dependent types so add the
                      typename std::vector<T>::iterator end) {  // "typename" keyword.
    std::shuffle(begin, end, rd_);
  }

 private:
  /**
   * Disallows creating an instance.  All methods will be static.
   */
  Random() = default;
  /**
   * Produces random integer values i, uniformly distributed on the
   * closed interval [a, b].
   */
  static std::uniform_int_distribution<int> uni_;
  /**
   * This is the configuration for the entire outermost (i.e., master) solver.
   */
  static SolverConfiguration * master_config_;
  /**
   * Source of randomness for the sampler.
   */
  static std::random_device rd_;
  /**
   * Used for generating random integers and doubles.
   */
  static std::mt19937 rng_;
};

#endif //PROJECT_RANDOM_H
