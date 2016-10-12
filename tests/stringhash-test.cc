#include <fstream> // std::ifstream
#include <iostream> // std::cout, endl, etc.
#include <iomanip> // std::setw, std::left, etc.
#include <string>
#include <cstring> // std::strlen

#include "../stringhash.h"

static const int k_max_len = 10000;
static const int k_n_tohash = 100;

/* Output results of hash function on items in lexicon "words.txt" */
static void lexicon_hash_test(string_hash h, const std::string & lexicon_path) {
  char s[k_max_len];

  std::ifstream x(lexicon_path);
  if(x.is_open()) {
    int n = 0;
    while(x.peek() != EOF && n < k_n_tohash) {
      x.getline(s, k_max_len);
      uint32_t hash_result = h(s, std::strlen(s));
      std::cout << ++n << ": " << std::setw(20) << std::left
                << s << "\t->\t" << hash_result << std::endl;
    }
    x.close();
    std::cout << "Success.\n" << std::endl;
  } else {
    std::cout << "Error: could not open " << lexicon_path << std::endl;
  }
}

/* Output results of applying the given hash function to the data of the  
 * vectors [0], [0, 1], [0, 1, 2], ..., [1, 2, 3, ..., 99]. */
static void vec_hash_test(string_hash h) {
  std::vector<int> v;
  for(int i = 0; i < k_n_tohash; i++) {
    v.push_back(i);
    std::cout << std::setw(4) << std::left << i
              << std::setw(12) << std::left << vector_hash(v, h) 
              << std::endl;
  }
  std::cout << "Success.\n" << std::endl;  
}

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "usage: " << argv[0] << " [lexicon path]" << std::endl;
    return 1;
  }

  std::string lexicon_path = argv[1];
  string_hash h = SuperFastHash; // by Paul Hsieh
  std::cout << "Testing hash function as string hash on the first " << k_n_tohash
            << " words of " << lexicon_path << "..." << std::endl;
  lexicon_hash_test(h, lexicon_path);

  std::cout << "Testing hash function as vector hash on [0], [0,1], ..., [0,1,..."
            << k_n_tohash << "]..." << std::endl;
  vec_hash_test(h);

  return 0;
}
