#ifndef BYTEME_READER_HPP
#define BYTEME_READER_HPP

/**
 * @file Reader.hpp
 *
 * @brief Read an input source.
 */

namespace byteme {

/**
 * @brief Virtual class for reading bytes from a source.
 */
class Reader {
public:
    virtual ~Reader() = default;

    /**
     * Read the next chunk of bytes from the input source. 
     *
     * To read the entire source, this function should be called repeatedly until `false` is returned.
     * Note that `buffer()` and `available()` may still be valid on the last invocation (i.e., the one that returns `false`),
     * as some bytes may have been read before reaching the end of the file.
     *
     * @return Boolean indicating whether there are still bytes remaining in the source.
     */
    virtual bool operator()() = 0;

    /**
     * @return Pointer to the start of an array containing the decompressed bytes.
     * The number of available bytes is provided in `available()`.
     */
    virtual const unsigned char* buffer() const = 0;

    /**
     * @return Number of decompressed bytes available in the `buffer()`.
     */
    virtual size_t available() const = 0;
};

}

#endif
