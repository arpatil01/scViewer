#ifndef BUFFIN_TEMP_FILE_PATH_HPP
#define BUFFIN_TEMP_FILE_PATH_HPP

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

#include <random>
#include <chrono>
#include <string>
#include <cstdint>
#include <fstream>

/**
 * @file temp_file_path.hpp
 *
 * @brief Define a temporary file path, usually for testing.
 */

namespace buffin {

/**
 * @param prefix Prefix of the basename of the file.
 *
 * @return Path to a new file in an appropriate temporary directory with a unique name.
 * The file itself is created at that location, though this may not be thread-safe.
 *
 * This function is wholly intended for unit testing in **buffin** and downstream libraries.
 * Production use should use OS-specific thread-safe alternatives such as `mkstemp()`. 
 */
inline std::string temp_file_path(std::string prefix) {
    auto path = fs::temp_directory_path();
    path.append(prefix);
    
    // Hopefully, we create a new seed when the function re-runs.
    uint64_t seed;
    try {
        std::random_device rd;
        seed = rd();
    } catch (...) {
        auto now = std::chrono::system_clock::now();
        seed = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();
    }

    std::mt19937_64 rng(seed);
    while (1) {
        auto path2 = path;
        path2 += std::to_string(rng());

        if (!fs::exists(path2)) {
            // Force the creation of the file. This is not entirely thread-safe
            // as the existence check is done separately, but there's no way to
            // guarantee safety without OS-level mechanics.
            std::ofstream dummy(path2, std::ofstream::out);
            return path2;
        }
    }

    return path;
}

}

#endif
