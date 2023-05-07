#ifndef BUFFIN_RAW_BUFFER_READER_HPP
#define BUFFIN_RAW_BUFFER_READER_HPP

#include <algorithm>
#include "Reader.hpp"

/**
 * @file RawBufferReader.hpp
 *
 * @brief Read bytes from a raw buffer without any extra transformations.
 */

namespace byteme {

/**
 * @brief Read bytes from a raw buffer, usually text.
 *
 * This is a wrapper around an input buffer, provided for consistency with the other `*Reader` classes.
 * We assume that the lifetime of the data in the `buffer` pointer exceeds the lifetime of the instance.
 */
class RawBufferReader : public Reader {
public:
    /**
     * @param buffer Buffer containing text.
     * @param length Length of the buffer.
     */
    RawBufferReader(const unsigned char* buffer, size_t length) : buffer_(buffer), len_(length) {}

    bool operator()() {
        return false;
    }

    const unsigned char* buffer() const {
        return buffer_;
    }

    size_t available() const {
        return len_;
    }

private:
    const unsigned char* buffer_;
    size_t len_;
};

}

#endif
