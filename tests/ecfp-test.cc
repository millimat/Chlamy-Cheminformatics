#include <iostream>

#include "../molecule.h"
#include "../ecfp_utils.h"

/* Generate a Molecule based on the sdf or molfile at the target file, report its 
 * basic information, generate an ECFP for it, and print out the identifiers of that
 * ECFP. */
static void basic_ecfp_test(int n_iterations, const std::string & filename, 
                            const element_map *elements, string_hash h) {
  Molecule m(filename, filename, elements, h);
  std::cout << m.atoms_report() << std::endl;

  m.generate_ecfp(n_iterations);
  std::vector<uint32_t> results;
  m.fetch_ecfp(results);

  std::cout << "ECFP_" << 2*n_iterations << " for " 
            << m.get_name() << ":\n";
  for(uint32_t i: results) std::cout << i << "\n";
  std::cout << std::endl;
}

/* Output an ECFP fingerprint for a single molecule at the provided path. */
int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cout << "usage: " << argv[0] << " [niterations] [mol/sdf path]" << std::endl;
    return 1;
  }
  
  const element_map *elements = symbols_to_weights();
  int n_iterations = std::stoi(argv[1]);
  std::string filename = argv[2];
  basic_ecfp_test(n_iterations, filename, elements, SuperFastHash);

  delete elements;
  return 0;
}

