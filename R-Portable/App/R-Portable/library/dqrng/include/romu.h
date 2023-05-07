// Romu Pseudorandom Number Generators
//
// Copyright 2020 Mark A. Overton
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// ------------------------------------------------------------------------------------------------
//
// Website: romu-random.org
// Paper:   http://arxiv.org/abs/2002.11331

#ifndef ROMU_H
#define ROMU_H 1

#include <array>
#include <mystdint.h>
#include <functional>
#include <algorithm>

#define ROTL(d,lrot) ((d<<(lrot)) | (d>>(8*sizeof(d)-(lrot))))

namespace dqrng {
class romuTrio {
public:
  using result_type = uint64_t;

private:
  std::array<result_type, 3> state;

  struct SplitMix {
    SplitMix(const uint64_t& k) : state(k) {}

    uint64_t operator() () {
      uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      return z ^ (z >> 31);
    }

  private:
    uint64_t state;
  };

  result_type next() {
      uint64_t xp = state[0], yp = state[1], zp = state[2];
      state[0] = 15241094284759029579u * zp;
      state[1] = yp - xp;
      state[1] = ROTL(state[1],12);
      state[2] = zp - yp;
      state[2] = ROTL(state[2],44);
      return xp;
  }

public:
  inline static constexpr result_type min() {return 0.0;};
  inline static constexpr result_type max() {return UINT64_MAX;};

  romuTrio(result_type _seed = 0x85c6ea9eb065ebeeULL) {
    seed(_seed);
  }

  void seed(std::function<result_type(void)> rng) {
    std::generate(state.begin(), state.end(), rng);
  }

  void seed(result_type _seed) {
    seed(SplitMix(_seed));
  }

  result_type operator() () {
    return next();
  }

};

class romuDuo {
public:
  using result_type = uint64_t;

private:
  // uint64_t xState, yState;  // set to nonzero seed
  uint64_t xState, yState, zState;  // set to nonzero seed

  struct SplitMix {
    SplitMix(const uint64_t& k) : state(k) {}

    uint64_t operator() () {
      uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      return z ^ (z >> 31);
    }

  private:
    uint64_t state;
  };

  result_type next() {
      /* uint64_t xp = xState; */
      /* xState = 15241094284759029579u * yState; */
      /* yState = ROTL(yState,36) + ROTL(yState,15) - xp; */
      /* return xp; */
      uint64_t xp = xState, yp = yState, zp = zState;
      xState = 15241094284759029579u * zp;
      yState = yp - xp;  yState = ROTL(yState,12);
      zState = zp - yp;  zState = ROTL(zState,44);
      return xp;
  }

public:
  inline static constexpr result_type min() {return 0.0;};
  inline static constexpr result_type max() {return UINT64_MAX;};

  romuDuo(result_type _seed = 0x85c6ea9eb065ebeeULL) {
    seed(_seed);
  }

  void seed(std::function<result_type(void)> rng) {
      xState = rng();
      yState = rng();
      zState = rng();
  }

  void seed(result_type _seed) {
    seed(SplitMix(_seed));
  }

  result_type operator() () {
    return next();
  }
};
}

#endif // ROMU_H
