#ifndef STRING_HASH
#define STRING_HASH

#include <cstdint> // for uint32_t
#include <vector> // for std::vector

// A hash function that takes a char * with the provided length in bytes and returns a 32-bit unsigned.
using string_hash = uint32_t (*) (const char *, int);
uint32_t vector_hash(const std::vector<int> &vec, string_hash h);

uint32_t SuperFastHash(const char *data, int len);

#endif // STRING_HASH
