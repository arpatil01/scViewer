//  Copyright (c) 2014 M.A. (Thijs) van den Berg, http://sitmo.com/
//
//  Use, modification and distribution are subject to the MIT Software License. 
// (See accompanying file LICENSE.txt)


#ifndef SITMO_RANDOM_VANDERCORPUT_HHP
#define SITMO_RANDOM_VANDERCORPUT_HHP 1

#include <limits>

namespace sitmo {

    namespace detail {

        template <typename UIntType>  
        constexpr std::size_t 
        max_base_digits_imp(UIntType base, UIntType maxval)
        { return (maxval >= base - 1) ? 1 + max_base_digits_imp<UIntType>(base, (maxval - base + 1)/base) : 0; };

        /**
         * @brief computed the largest n such that base^n - 1 fits in a UIntType
         */
        template <typename UIntType>  
        constexpr std::size_t 
        max_base_digits(std::size_t base)
        { return max_base_digits_imp<UIntType>(base, std::numeric_limits<UIntType>::max()); };
        
        
        template <typename UIntType>
        constexpr UIntType 
        max_base_value_imp(UIntType B, UIntType U)
        { return (U==0) ? 0 : B - 1 + B * max_base_value_imp<UIntType>(B, (U-1)/B ); };

        /**
         * @brief computed the largest value of the form base^n - 1 that fits in a UIntType
         */
        template <typename UIntType>
        constexpr UIntType 
        max_base_value(UIntType B)
        { return max_base_value_imp<UIntType>(B, std::numeric_limits<UIntType>::max() / (B-1) ); };


    } // namespace detail

/**
* @brief The van der Corput low discrepancy random number generator.
*
* A random number generator that produces pseudorandom numbers using 
* van der Corput sequences. It is constructed by reversing the base n 
* representation of the sequence of natural numbers 1,2,3,..
*
* The template parameter @p UIntType must be an unsigned integral type.
* The template parameter @p base must be an integral number of at last 2
* and at most std::numeric_limits<UIntType>::max().
*
* The size of the state is UIntType.
*/
template <typename UIntType = uint_fast64_t, size_t base_ = 2 >
class vandercorput_engine
{
    static_assert(std::is_unsigned<UIntType>::value || std::is_integral<UIntType>::value, 
        "template argument UIntType not an unsigned integral type");
    static_assert(base_ >= 2u || (base_ <= std::numeric_limits<UIntType>::max() ),
        "template argument base out of bounds");
public:
    /** The type of the generated random value. */
    typedef UIntType result_type;
    
    /** The default seed. */
    static constexpr result_type default_seed = 1;
    
    /** The base used for the sequence.. */
    static constexpr size_t base = base_;  
    
    
    /**
     * @brief Gets the smallest possible value in the output range.
     */
    static constexpr result_type 
    (min)() 
    { return 0; }
    
    
    /**
     * @brief Gets the largest possible value in the output range.
     */
    static constexpr result_type 
    (max)() 
    { return detail::max_base_value<UIntType>(base); }
    
    
    /**
     * @brief Constructs a %vandercorput_engine random number
     *        generator engine with seed @p __s.  The default seed value
     *        is 1.
     *
     * @param __s The initial seed value.
     */
    explicit 
    vandercorput_engine(result_type __s = default_seed) 
    { seed(__s); }
    
    
    /**
     * @brief Constructs a %vandercorput_engine random number
     *        generator engine seeded from the seed sequence @p __q.
     *
     * @param __q the seed sequence.
     */
    template<typename _Sseq, typename = typename
        std::enable_if<
            !std::is_same<_Sseq, vandercorput_engine>::value
        >::type>
    vandercorput_engine(_Sseq& __q) { seed(__q); }
    
    
    /**
     * @brief Reseeds the  %vandercorput_engine random number
     *        generator engine with the seed @p __s.
     *        The default seed value is 1.
     *
     * @param __s the new seed.
     */
    void 
    seed(result_type __s = default_seed)
    { _M_counter = __s;}
    
    
    /**
     * @brief Seeds the  %vandercorput_engine random number
     *        generator engine with the seed sequence @p __q.
     *
     * @param __p the seed sequence.
     */
    template<typename _Sseq> 
        typename std::enable_if< std::is_class< _Sseq >::value >::type 
    seed(_Sseq& __q)
        {
        const size_t samples = 1 + ( (sizeof(UIntType) - 1) / sizeof(typename _Sseq::result_type) );
        typename _Sseq::result_type w[samples];
        __q.generate(&w[samples], &w[samples]);
        _M_counter = static_cast<UIntType>(w[0]);
    }
    
    /**
     * @brief Generate a random sample.
     */
    result_type 
    operator()()
    {
        result_type vdc = 0;
        result_type n = _M_counter;
        for (unsigned int i=0; i<detail::max_base_digits<UIntType>(base); ++i)
        {
            vdc *= base;
            vdc += n % base;
            n /= base;
        }
        ++_M_counter;
        return vdc;
    }

    /**
     * @brief Discard a number of elements from the random numbers sequence.
     *
     * @param __z the number of elements to discard.
     */
    void 
    discard(result_type __z) 
    { _M_counter += __z; }

    /**
     * @brief Writes the textual representation of the state x(i) of x to
     *        @p __os.
     *
     * @param __os  The output stream.
     * @param __eng A % vandercorput_engine random number generator.
     * @returns __os.
     */
   template<typename _CharT, typename _Traits>
   friend std::basic_ostream<_CharT, _Traits>&
   operator<<(std::basic_ostream<_CharT, _Traits>& __os, const vandercorput_engine& __eng) {
       __os << __eng._M_counter;
       return __os;
   }

    /**
     * @brief reads the textual representation of the state x(i) of x from
     *        @p __is.
     *
     * @param __is  The input stream.
     * @param __eng A % vandercorput_engine random number generator.
     * @returns __is.
     */
   template<typename _CharT, typename _Traits>
   friend std::basic_istream<_CharT, _Traits>&
   operator>>(std::basic_istream<_CharT, _Traits>& __is, vandercorput_engine& __eng) {
       __is >> __eng._M_counter;
       return __is;
   }

    /**
     * @brief Compares two van der Corput random number generator
     * objects of the same type for equality.
     *
     *
     * @param __lhs A van der Corput random number generator object.
     * @param __rhs Another van der Corput random number generator
     *              object.
     *
     * @returns true if the infinite sequences of generated values
     *          would be equal, false otherwise.
     */
    friend bool 
    operator==(const vandercorput_engine& __lhs, const vandercorput_engine& __rhs) 
    { return (__lhs._M_counter == __rhs._M_counter); }

    /**
     * @brief Compares two van der Corput random number generator
     * objects of the same type for inequality.
     *
     *
     * @param __lhs A van der Corput random number generator object.
     * @param __rhs Another van der Corput random number generator
     *              object.
     *
     * @returns true if the infinite sequences of generated values
     *          would not be equal, false otherwise.
     */
    friend bool 
    operator!=(const vandercorput_engine& __lhs, const vandercorput_engine& __rhs) 
    { return (__lhs._M_counter != __rhs._M_counter); }
    
private:
    UIntType _M_counter;
};

/**
 * The base 2, 32 bit van der Corput engine.
 */
typedef vandercorput_engine<uint_fast32_t, 2> vandercorput;

/**
 * The base 2, 64 bit van der Corput engine.
 */
 typedef vandercorput_engine<uint_fast64_t, 2> vandercorput_64;

} // namespace sitmo
#endif
