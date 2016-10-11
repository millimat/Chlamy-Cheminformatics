#include <string>
#include <iostream>
#include <fstream>

#include "molecule.h" // for Molecule
#include "ecfp_utils.h"
#include "dirent.h" // For directory traversal

static const std::string kWeakDirectory = "lt45";
static const std::string kMedDirectory = "lvpt45";
static const std::string kStrDirectory = "lspt45";

static void get_subject_ecfp(int niterations, std::string subject_path, std::vector<uint32_t> & to_fill) {
  // std::cout << "Getting ECFP_" << 2 * niterations << " fingerprint for subject molecule at " 
  //           << subject_path << "..." << std::endl;
  const element_map *elements = symbols_to_weights();
  Molecule m(subject_path, subject_path, elements, SuperFastHash);
  m.generate_ecfp(niterations);
  m.fetch_ecfp(to_fill);
  delete elements;
}

static void get_comp_ecfps(int niterations, std::string compdir, const char * prg,
                           std::vector<std::vector<uint32_t> *> & to_fill) {
  // std::cout << "Getting ECFP_" << 2 * niterations << " fingerprints for all molecules in "
  //           << "comparison directory " << compdir << " ... " << std::endl;
  std::vector<Molecule *> comp_molecules;
  fetch_directory_molecules(prg, compdir, comp_molecules);

  for(Molecule *m: comp_molecules) {
    m->generate_ecfp(niterations);
    std::vector<uint32_t> * next_ecfp = new std::vector<uint32_t>;
    m->fetch_ecfp(*next_ecfp);
    to_fill.push_back(next_ecfp);
    delete m;
  }
}

static double get_representation(const std::vector<uint32_t> & search, const std::vector<uint32_t> & compare) {
  if(search.empty()) return 0;

  int intersection_size = 0;
  auto it_s = search.cbegin(), it_c = compare.cbegin();
  while(it_s < search.cend() && it_c < compare.cend()) {
    if(*it_s < *it_c) it_s++;
    else if(*it_c < *it_s) it_c++;
    else { intersection_size++; it_s++; it_c++; }
  }

  return (double)intersection_size / search.size();  
}

int main(int argc, char *argv[]) {
  if(argc != 5) {
    std::cout << "usage: " << argv[0] << " [niterations] [match threshold] [subject molecule]"
              << "[comparison directory]" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  int niterations = get_niterations(argv[0], argv[1]);
  double match_threshold = get_display_threshold(argv[0], argv[2]);

  std::vector<uint32_t> subject_ecfp;
  get_subject_ecfp(niterations, argv[3], subject_ecfp);
  
  std::vector<std::vector<uint32_t> *> comp_ecfps;
  get_comp_ecfps(niterations, argv[4], argv[0], comp_ecfps);

  int n_matches = 0, n_mismatches = 0;
  for(std::vector<uint32_t> * to_compare: comp_ecfps) {
    double match_quality = get_representation(subject_ecfp, *to_compare);
    if(match_quality >= match_threshold) n_matches++;
    else n_mismatches++;
    
    delete to_compare;
  }
  
  std::cout << "Searched for " << argv[3] << " in " << argv[4] << " with match threshold " 
            << match_threshold << ":\n" << n_matches << " matches, " << n_mismatches
            << " mismatches" << std::endl;
  return 0;
}
