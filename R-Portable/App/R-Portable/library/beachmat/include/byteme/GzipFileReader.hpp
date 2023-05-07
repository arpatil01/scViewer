#ifndef BYTEME_GZIP_FILE_READER_HPP
#define BYTEME_GZIP_FILE_READER_HPP

#include "zlib.h"
#include <stdexcept>
#include <vector>
#include <string>
#include "Reader.hpp"

/**
 * @file GzipFileReader.hpp
 *
 * @brief Read a Gzip-compressed file.
 */

namespace byteme {

/**
 * @brief Read uncompressed bytes from a Gzip-compressed file.
 *
 * This is basically a wrapper around Zlib's `gzFile` with correct closing and error checking.
 */
class GzipFileReader : public Reader {
private:
    /**
     * @cond
     */
    struct GZFile {
        GZFile(const char* path) : handle(gzopen(path, "rb")) {
            if (!handle) {
                throw std::runtime_error("failed to open file at '" + std::string(path) + "'");
            }
            return;
        }

        ~GZFile() {
            gzclose(handle);
            return;
        }

        // Delete the remaining constructors.
        GZFile(const GZFile&) = delete;
        GZFile(GZFile&&) = delete;
        GZFile& operator=(const GZFile&) = delete;
        GZFile& operator=(GZFile&&) = delete;

        gzFile handle;
    };
    /**
     * @endcond
     */

public:
    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    GzipFileReader(const char* path, size_t buffer_size = 65536) : gz(path), buffer_(buffer_size), read(0) {}

    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    GzipFileReader(const std::string& path, size_t buffer_size = 65536) : GzipFileReader(path.c_str(), buffer_size) {}

    bool operator()() {
        read = gzread(gz.handle, buffer_.data(), buffer_.size());
        if (read == 0) {
            if (!gzeof(gz.handle)) { 
                int dummy;
                throw std::runtime_error(gzerror(gz.handle, &dummy));
            }
            return false;
        } else {
            return true;
        }
    }

    const unsigned char* buffer() const {
        return buffer_.data();
    }

    size_t available() const {
        return read;
    }

private:
    GZFile gz;
    std::vector<unsigned char> buffer_;
    size_t read;
};

}

#endif
