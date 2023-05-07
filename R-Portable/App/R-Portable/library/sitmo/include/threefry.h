/* stdfin random/threefry.hpp header file
 *
 * Copyright Thijs van den Berg 2014-2015
 * 
 * Distributed under the MIT Software License.
 * See the accompanying file LICENSE or copy at http://opensource.org/licenses/MIT
 *
 */

#ifndef SITMO_RANDOM_THREEFRY_HPP
#define SITMO_RANDOM_THREEFRY_HPP

#include <istream>
#include <algorithm>
#include <cstdint>
#include <limits>

//debug
#include <iostream>
#include <iomanip>

namespace sitmo {

/**
* @brief The threefry random engine is a counter based random engine that uses a 
* stripped-down Threefish cryptographic function that is optimised for speed.
*
* The template parameter @p UIntType the return type of the random engine, must be an 
* unsigned integral type, e.g. uint64_t.
* The template parameter @p w is the width of the return type, the number of bits. 
* Valid values are 8,16,32,64. 
* The template parameter @p r the number of mix-rounds. 13 round or higher give very good 
* quality random numbers.
*
* @blockquote
* title:     Parallel random numbers: as easy as 1, 2, 3
* authors:   Salmon, John K. and Moraes, Mark A. and Dror, Ron O. and Shaw, David E.
* booktitle: Proceedings of 2011 International Conference for High Performance Computing, Networking, Storage and Analysis
* publicher: ACM
* year:      2011
* isbn:      978-1-4503-0771-0
* @endblockquote
*
* @xmlnote
* The stdfin variant has been implemented from scratch and does not
* derive from or use the Random123 library provided at http://www.thesalmons.org/john/random123/
* However, it  was verified that both produce identical output.
* Output was verified against the threefry4x64 unit test cases from https://github.com/girving/random123/blob/master/examples/kat_vectors
* @endxmlnote
*/
namespace detail {

static const uint_least64_t threefry4x64_tweak = 0x1BD11BDAA9FC1A22;

// primary template
template< typename UIntType, std::size_t w>
struct extract4x64_impl {
  inline static UIntType zth(const uint_least64_t (&_output)[4]);
  inline static UIntType nth(const uint_least64_t (&_output)[4], std::size_t n);
  inline static constexpr UIntType w_max();
};

// specialisation
template< typename UIntType>
struct extract4x64_impl<UIntType,64> {
  inline static UIntType zth(const uint_least64_t (&_output)[4])
  { return _output[0]; }
  inline static UIntType nth(const uint_least64_t (&_output)[4], std::size_t n)
  { return _output[n]; }
  inline static constexpr UIntType w_max()
  { return 0xFFFFFFFFFFFFFFFF; }
};
template< typename UIntType>
struct extract4x64_impl<UIntType,32> {
  inline static UIntType zth(const uint_least64_t (&_output)[4])
  { return _output[0] & 0xFFFFFFFF; }
  inline static UIntType nth(const uint_least64_t (&_output)[4], std::size_t n)
  { return (_output[n>>1] >> ((n&1)<<5)) & 0xFFFFFFFF; }
  inline static constexpr UIntType w_max()
  { return 0xFFFFFFFF; }
};
template< typename UIntType>
struct extract4x64_impl<UIntType,16> {
  inline static UIntType zth(const uint_least64_t (&_output)[4])
  { return _output[0] & 0xFFFF; }
  inline static UIntType nth(const uint_least64_t (&_output)[4], std::size_t n)
  { return (_output[n>>2] >> ((n&3)<<4)) & 0xFFFF; }
  inline static constexpr UIntType w_max()
  { return 0xFFFF; }
};
template< typename UIntType>
struct extract4x64_impl<UIntType,8> {
  inline static UIntType zth(const uint_least64_t (&_output)[4])
  { return _output[0] & 0xFF; }
  inline static UIntType nth(const uint_least64_t (&_output)[4], std::size_t n)
  { return (_output[n>>3] >> ((n&7)<<3)) & 0xFF; }
  inline static constexpr UIntType w_max()
  { return 0xFF; }
};

template <class SeedSeq, class Eng>
struct is_seed_sequence {
  static constexpr bool value =
    !std::is_convertible<SeedSeq, typename Eng::result_type>::value &&
    !std::is_same<typename std::remove_cv<SeedSeq>::type, Eng>::value;
};    
}

template <  typename UIntType=uint32_t,  // the return type
            std::size_t w=32,            // number of bits in the return type
            std::size_t r=13             // number of rounds
>
class threefry_engine
{
public:
  // check for supported nr of bits in return type
  static_assert( w==8 || w==16 || w==32 || w==64, "threefry_engine: invalid template argument. Number of bits must be 8,16,32, or 64");
  
  // types
  typedef UIntType result_type;
  
  // engine characteristics
  static constexpr std::size_t word_size = w;
  static constexpr std::size_t rounds = r;
  static constexpr std::size_t samples_per_block = 256/w;
  
  /**
  * @brief Constructs the defafult %threefry_engine.
  */
  threefry_engine() 
  { seed(0); }
  
  /**
  * @brief Constructs a %threefry_engine random number
  *        generator engine with seed @p value.  The default seed value
  *        is 0.
  *
  * @param value The initial seed value.
  */
  explicit
    threefry_engine(UIntType value)
    { seed(value); }
  
  
  /**
  * @brief Constructs a %threefry_engine random number
  *        generator engine seeded from the seed sequence @p seq.
  *
  * @param seq the seed sequence.
  */
  template <class SeedSeq>
  explicit
    threefry_engine(SeedSeq& seq, typename std::enable_if<detail::is_seed_sequence<SeedSeq, threefry_engine>::value>::type* = 0)
    {
      seed(seq);
    }
  
  /**
  * @brief Re-seed the %threefry_engine to it's default seed.
  */
  void seed()
  {
    seed(0);
  }
  
  /**
  * @brief Re-seed the  %threefry_engine random number
  *        generator engine with the seed @p value.
  *        The default seed value is 0.
  *
  * @param value the new seed.
  */
  void seed(result_type value)
  {
    _key[0] = value;
    _key[1] = 0;
    _key[2] = 0;
    _key[3] = 0;
    reset_after_key_change();
  }
  
  /**
  * @brief Seeds the %threefry_engine random number
  *        generator engine with the seed sequence @p seq.
  *
  * @param seq the seed sequence.
  */
  template <class SeedSeq>
  void seed(SeedSeq& seq)
  {
    seed(seq, typename std::is_fundamental<SeedSeq>::type());
    reset_after_key_change();
  }
  
  
  /**
  * @brief Gets the smallest possible value in the output range.
  */
  static constexpr result_type (min)()
  { return 0; }
  
  /**
  * @brief Gets the largest possible value in the output range.
  */
  static constexpr result_type (max)()
  { return detail::extract4x64_impl<UIntType,w>::w_max(); }
  
  
  /**
  * @brief Generate a random sample.
  */
  result_type operator()()
  {
    // can we return a value from the current block?
    if (_o_counter < samples_per_block)
      return detail::extract4x64_impl<UIntType,w>::nth(_output, _o_counter++);
    
    // generate a new block and return the first result_type 
    inc_counter();
    encrypt_counter();
    _o_counter = 1; // the next call
    return detail::extract4x64_impl<UIntType,w>::zth(_output);
  }
  
  /**
  * @brief Discard a number of elements from the random numbers sequence.
  *
  * @param z the number of elements to discard.
  */
  void discard(unsigned long long z)
  {
    // check if we stay in the current block
    if (z < samples_per_block - _o_counter) {
      _o_counter += static_cast<unsigned short>(z);
      return;
    }
    
    // we will have to generate a new block...
    z -= (samples_per_block - _o_counter);  // discard the remainder of the current blok
    _o_counter = z % samples_per_block;     // set the pointer in the correct element in the new block
    z -= _o_counter;                        // update z
    z /= samples_per_block;                 // the number of 256 bit bocks to skip is z/samples_per_block
    ++z;                                    // and one more because we crossed the buffer line
    inc_counter(z);
    encrypt_counter();
  }
  
  /**
  * @brief Writes the textual representation of the state x(i) of x to
  *        @p _os.
  *
  * @param os  The output stream.
  * @param eng A %threefry_engine random number generator.
  * @returns os.
  */
  template<class CharT, class Traits>
  friend std::basic_ostream<CharT, Traits>&
  operator<<(std::basic_ostream<CharT, Traits>& os, const threefry_engine& eng)
  {
    for (unsigned short i=0; i<4; ++i)
      os << eng._key[i] << ' ';
    
    for (unsigned short i=0; i<4; ++i)
      os << eng._counter[i] << ' ';
    
    os << eng._o_counter;
    return os;
  }
  
  /**
  * @brief reads the textual representation of the state x(i) from
  *        @p is.
  *
  * @param is  The input stream.
  * @param eng A %threefry_engine random number generator.
  * @returns is.
  */
  template<class CharT, class Traits>
  friend std::basic_istream<CharT, Traits>&
  operator >> (std::basic_istream<CharT, Traits>& is, threefry_engine& eng)
  {
    for (unsigned short i=0; i<4; ++i) 
      is >> eng._key[i] >> std::ws;
    
    for (unsigned short i=0; i<4; ++i) 
      is >> eng._counter[i] >> std::ws;
    
    is >> eng._o_counter;
    
    eng._key[4] = detail::threefry4x64_tweak ^ eng._key[0] ^ eng._key[1] ^ eng._key[2] ^ eng._key[3];
    eng.encrypt_counter();
    return is;
  } 
  
  /**
  * @brief Compares two %threefry_engine
  * objects of the same type for equality.
  *
  *
  * @param _lhs A threefry engine.
  * @param _rhs Another threefry engine.
  *
  * @returns true if the infinite sequences of generated values
  *          would be equal, false otherwise.
  */
  friend bool 
    operator==(const threefry_engine& _lhs, const threefry_engine& _rhs) 
    {
      if (_lhs._o_counter != _rhs._o_counter) return false;
      
      for (unsigned short i=0; i<4; ++i) {
        if (_lhs._counter[i] != _rhs._counter[i]) return false;
        if (_lhs._key[i]     != _rhs._key[i])     return false;
        if (_lhs._output[i]  != _rhs._output[i])  return false;
      }
      return true;
    }
  
  /**
  * @brief Compares two threefry engines for inequality.
  *
  *
  * @param _lhs A threefry engine.
  * @param _rhs Another threefry engine.
  *
  * @returns true if the infinite sequences of generated values
  *          would not be equal, false otherwise.
  */
  friend bool 
    operator!=(const threefry_engine& _lhs, const threefry_engine& _rhs) 
    { return !(_lhs == _rhs); }
  
  // Extra function to set the key
  void set_key(uint64_t k0=0, uint64_t k1=0, uint64_t k2=0, uint64_t k3=0)
  {
    _key[0] = k0;
    _key[1] = k1;
    _key[2] = k2;
    _key[3] = k3;
    encrypt_counter();
  }
  
  // set the counter
  void set_counter(uint64_t s0=0, uint64_t s1=0, uint64_t s2=0, uint64_t s3=0, unsigned short o_counter=0)
  {
    _counter[0] = s0; 
    _counter[1] = s1; 
    _counter[2] = s2; 
    _counter[3] = s3;
    _o_counter = o_counter % 8;
    encrypt_counter();
  }
  
private:
  template<class Gen>
  void seed(Gen& gen, std::true_type)
  {
    return seed(static_cast<result_type>(gen));
  }
  
  template<class Gen>
  void seed(Gen& gen, std::false_type)
  {
    typename Gen::result_type data[8];
    gen.generate(&data[0], &data[8]);
    
    for (unsigned short i=0; i<4; ++i) {
      _key[i] = ( static_cast<uint64_t>(data[2*i]) << 32) | data[2*i+1];
      _counter[i] = 0;
    }
    _o_counter = 0;
  }
  
  
  inline void rotl64(uint_least64_t& v, const uint8_t bits) const
  { 
    v = (v << bits) | (v >> (sizeof(uint_least64_t) * 8 - bits)); 
  }
  
  inline void mix64(uint_least64_t& x0, uint64_t& x1, const uint8_t bits) const
  {
    x0 += x1;
    rotl64(x1, bits);
    x1 ^= x0;
  }
  
  inline void double_mix64(   uint_least64_t& x0, uint_least64_t& x1, const uint8_t rx,
                              uint_least64_t& z0, uint_least64_t& z1, const uint8_t rz) const
  {
    mix64(x0,x1,rx);
    mix64(z0,z1,rz);
  }
  
  template <std::size_t offset>
  inline void add_key64_t( uint_least64_t (&output)[4], uint_least64_t (&key)[5], const std::size_t c) const
  {
    output[0] += key[ offset   %5];
    output[1] += key[(offset+1)%5];
    output[2] += key[(offset+2)%5];
    output[3] += key[(offset+3)%5];
    output[3] += c;
  }
  
  template <std::size_t Rounds>
  inline void encrypt_counter_t(std::size_t& four_cycles)
  {
    if (Rounds>=1) double_mix64( _output[0], _output[1], 14, _output[2], _output[3], 16);
    if (Rounds>=2) double_mix64( _output[0], _output[3], 52, _output[2], _output[1], 57);
    if (Rounds>=3) double_mix64( _output[0], _output[1], 23, _output[2], _output[3], 40);
    if (Rounds>=4) {
      double_mix64( _output[0], _output[3],  5, _output[2], _output[1], 37);
      add_key64_t<1>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=5) double_mix64( _output[0], _output[1], 25, _output[2], _output[3], 33);
    if (Rounds>=6) double_mix64( _output[0], _output[3], 46, _output[2], _output[1], 12);
    if (Rounds>=7) double_mix64( _output[0], _output[1], 58, _output[2], _output[3], 22);
    if (Rounds>=8) {
      double_mix64( _output[0], _output[3], 32, _output[2], _output[1], 32);
      add_key64_t<2>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=9)  double_mix64( _output[0], _output[1], 14, _output[2], _output[3], 16);
    if (Rounds>=10) double_mix64( _output[0], _output[3], 52, _output[2], _output[1], 57);
    if (Rounds>=11) double_mix64( _output[0], _output[1], 23, _output[2], _output[3], 40);
    if (Rounds>=12) {
      double_mix64( _output[0], _output[3],  5, _output[2], _output[1], 37);
      add_key64_t<3>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=13) double_mix64( _output[0], _output[1], 25, _output[2], _output[3], 33);
    if (Rounds>=14) double_mix64( _output[0], _output[3], 46, _output[2], _output[1], 12);
    if (Rounds>=15) double_mix64( _output[0], _output[1], 58, _output[2], _output[3], 22);
    if (Rounds>=16) {
      double_mix64( _output[0], _output[3], 32, _output[2], _output[1], 32);
      add_key64_t<4>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=17) double_mix64( _output[0], _output[1], 14, _output[2], _output[3], 16);
    if (Rounds>=18) double_mix64( _output[0], _output[3], 52, _output[2], _output[1], 57);
    if (Rounds>=19) double_mix64( _output[0], _output[1], 23, _output[2], _output[3], 40);
    if (Rounds>=20) {
      double_mix64( _output[0], _output[3],  5, _output[2], _output[1], 37);
      add_key64_t<0>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=21) double_mix64( _output[0], _output[1], 25, _output[2], _output[3], 33);
    if (Rounds>=22) double_mix64( _output[0], _output[3], 46, _output[2], _output[1], 12);
    if (Rounds>=23) double_mix64( _output[0], _output[1], 58, _output[2], _output[3], 22);
    if (Rounds>=24) {
      double_mix64( _output[0], _output[3], 32, _output[2], _output[1], 32);
      add_key64_t<1>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=25) double_mix64( _output[0], _output[1], 14, _output[2], _output[3], 16);
    if (Rounds>=26) double_mix64( _output[0], _output[3], 52, _output[2], _output[1], 57);
    if (Rounds>=27) double_mix64( _output[0], _output[1], 23, _output[2], _output[3], 40);
    if (Rounds>=28) {
      double_mix64( _output[0], _output[3],  5, _output[2], _output[1], 37);
      add_key64_t<2>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=29) double_mix64( _output[0], _output[1], 25, _output[2], _output[3], 33);
    if (Rounds>=30) double_mix64( _output[0], _output[3], 46, _output[2], _output[1], 12);
    if (Rounds>=31) double_mix64( _output[0], _output[1], 58, _output[2], _output[3], 22);
    if (Rounds>=32) {
      double_mix64( _output[0], _output[3], 32, _output[2], _output[1], 32);
      add_key64_t<3>(_output, _key, ++four_cycles);
    }
    
    if (Rounds>=33) double_mix64( _output[0], _output[1], 14, _output[2], _output[3], 16);
    if (Rounds>=34) double_mix64( _output[0], _output[3], 52, _output[2], _output[1], 57);
    if (Rounds>=35) double_mix64( _output[0], _output[1], 23, _output[2], _output[3], 40);
    if (Rounds>=36) {
      double_mix64( _output[0], _output[3],  5, _output[2], _output[1], 37);
      add_key64_t<4>(_output, _key, ++four_cycles);
    }
    
    
    if (Rounds>=37) double_mix64( _output[0], _output[1], 25, _output[2], _output[3], 33);
    if (Rounds>=38) double_mix64( _output[0], _output[3], 46, _output[2], _output[1], 12);
    if (Rounds>=39) double_mix64( _output[0], _output[1], 58, _output[2], _output[3], 22);
    if (Rounds>=40) {
      double_mix64( _output[0], _output[3], 32, _output[2], _output[1], 32);
      add_key64_t<0>(_output, _key, ++four_cycles);
    }
  }
  
  void encrypt_counter()
  {
    for (std::size_t i=0; i<4; ++i) _output[i] = _counter[i];
    for (std::size_t i=0; i<4; ++i) _output[i] += _key[i];
    
    std::size_t four_cycles = 0;
    
    // do chunks of 40 rounds
    for (std::size_t big_rounds=0; big_rounds < r/40; ++big_rounds)
      encrypt_counter_t<40>(four_cycles);
    
    // the remaining rounds
    encrypt_counter_t<r - 40*(r/40)>(four_cycles);
  }
  
  // increment the counter with 1
  void inc_counter()
  {
    ++_counter[0]; if (_counter[0] != 0) return;
    ++_counter[1]; if (_counter[1] != 0) return;
    ++_counter[2]; if (_counter[2] != 0) return;
    ++_counter[3];
  }
  
  // increment the counter with z
  void inc_counter(uint_least64_t z)
  {
    if (z > 0xFFFFFFFFFFFFFFFF - _counter[0]) {   // check if we will overflow the first 64 bit of the counter
      ++_counter[1];
      if (_counter[1] == 0) {
        ++_counter[2];
        if (_counter[2] == 0) {
          ++_counter[3];
        }
      }
    }
    _counter[0] += z;
  }
  
  void reset_counter()
  {
    _counter[0] = 0;
    _counter[1] = 0;
    _counter[2] = 0;
    _counter[3] = 0;
    _o_counter  = 0;
  }
  
  // reset the counter to zero, and reset the keyx
  void reset_after_key_change()
  {
    _key[4] = detail::threefry4x64_tweak ^ _key[0] ^ _key[1] ^ _key[2] ^ _key[3];
    reset_counter();
    encrypt_counter();
  }
  
private:
  
  uint_least64_t _counter[4];   // the 256 bit counter (message) that gets encrypted
  uint_least64_t _output[4];    // the 256 bit cipher output 4 * 64 bit = 256 bit output
  uint_least64_t _key[5];       // the 256 bit encryption key
  uint_least16_t _o_counter;     // output chunk counter, e.g. for a 64 bit random engine
  // the 256 bit output buffer gets split in 4x64bit chunks or 8x32bit chunks chunks.
  /// \endcond
};

/**
* 32 bit version of the 13 rounds threefry engine
*/
typedef threefry_engine<uint32_t, 32, 13> threefry_13_32;

/**
 * 64 bit version of the 13 rounds threefry engine
 */
typedef threefry_engine<uint64_t, 64, 13> threefry_13_64;


/**
 * 32 bit version of the 20 rounds threefry engine
 */
typedef threefry_engine<uint32_t, 32, 20> threefry_20_32;

/**
 * 64 bit version of the 20 rounds threefry engine
 */
typedef threefry_engine<uint64_t, 64, 20> threefry_20_64;

typedef threefry_engine<uint32_t, 32, 13> threefry; // default

} // namespace sitmo

#endif // STDFIN_RANDOM_THREEFRY_HPP
