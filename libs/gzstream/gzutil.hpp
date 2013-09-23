#include <iostream>
#include <string>
#include <shared_ptr.hpp>
#include <gzstream.hpp>

/**
 * Returns true if the given string ends with given suffix.
 *
 * @param str A string.
 * @param end The desired end of the string.
 *
 * @return True if the given string ends with enb, false otherwise.
 */
bool ends_with(const std::string &str, const std::string &end);

/**
 * Opens a possibly gzipped file and returns a stream to it.
 *
 * @param path path to the file.
 *
 * @return An opened stream to the file.
 */
shared_ptr<std::istream> open_possible_gz(const std::string &path);

/**
 * Creates a possibly gzipped file and returns a stream to it.
 *
 * @param path path to the file.
 *
 * @return An opened stream to the file.
 */
shared_ptr<std::ostream> create_possible_gz(const std::string &path);
