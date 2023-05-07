#ifndef BYTEME_ISTREAM_READER_HPP
#define BYTEME_ISTREAM_READER_HPP

#include <istream>
#include <vector>
#include <stdexcept>
#include "Reader.hpp"

/**
 * @file IstreamReader.hpp
 *
 * @brief Read bytes from an input stream.
 */

namespace byteme {

/**
 * @brief Read bytes from a `std::istream`.
 *
 * This is just a wrapper around `std::istream::read`,
 * mostly to avoid having to remember the correct way to check for end of file.
 */
class IstreamReader : public Reader {
public:
    /**
     * @param input An input stream.
     * @param buffer_size Size of the buffer to use for reading.
     */
    IstreamReader(std::istream& input, size_t buffer_size = 65536) : ptr(&input), buffer_(buffer_size) {}

    bool operator()() {
        ptr->read(reinterpret_cast<char*>(buffer_.data()), buffer_.size());
        read = ptr->gcount();

        if (read < buffer_.size()) {
            if (ptr->eof()) {
                return false;
            } else {
                throw std::runtime_error("failed to finish reading the input stream");
            }
        } 

        return true;
    }

    const unsigned char* buffer() const {
        return buffer_.data();
    }

    size_t available() const {
        return read;
    }

private:
    std::istream* ptr;
    std::vector<unsigned char> buffer_;
    size_t read = 0;
};

}

#endif
