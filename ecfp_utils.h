#ifndef ECFP_UTILS_H
#define ECFP_UTILS_H

#include "molecule.h"
#include "stringhash.h"
#include "dirent.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

/* Construct a Map taking element symbols to their atomic numbers and the
 * masses of their most abundant isotopes. */
const element_map *symbols_to_weights() {
  element_map *result = new element_map;
  std::ifstream weight_table("data/weights.txt"); // file mapping symbols to weights
  std::string elem_temp, number_temp, weight_temp;

  while(weight_table.peek() != EOF) {
    getline(weight_table, elem_temp, ' '); // fetch element symbol from line
    getline(weight_table, number_temp, ' '); // fetch atomic number from line
    getline(weight_table, weight_temp); // fetch element weight from line
    (*result)[elem_temp] = { .atomic_number = stoi(number_temp), .atomic_mass = stoi(weight_temp) };
  }

  return result;
}

/* Print the elements of an element_map with their atomic numbers and masses. */
static void print_contents(const element_map *m) {
  for(std::pair<std::string, element> elem: *m) {
    std::cout << elem.first << ": " << elem.second.atomic_number
              << ", " << elem.second.atomic_mass << std::endl;
  }
}

/* Extract a valid number of ECFP iterations (i.e. a nonnegative integer)
 * from command-line argument arg. If arg is out of range or cannot be read
 * as an integer, exit the parent program (with name given by prg). */
int get_niterations(char *prg, char *arg) {
  int niterations;
  try { 
    niterations = std::stoi(arg); 
    if(niterations < 0) throw(std::invalid_argument(arg));
  } catch(std::invalid_argument ia) { 
    std::cout << prg << ": niterations error for argument " << ia.what() 
              << "; provided value must be a nonnegative integer" << std::endl;
    exit(EXIT_FAILURE);
  } 

  return niterations;
}

/* Extract a valid Tanimoto coefficient (i.e. a double between 0 and 1, inclusive)
 * from command-line argument arg. If arg is out of range or cannot be read as a
 * double, exit the parent program (whose name is given by prg). */
double get_display_threshold(char *prg, char *arg) {
  double cytoscape_tc_display_threshold;
  try { 
    cytoscape_tc_display_threshold = std::stod(arg); 
    if(cytoscape_tc_display_threshold < 0 || cytoscape_tc_display_threshold > 1) throw(std::invalid_argument(arg));
  } catch(std::invalid_argument ia) { 
    std::cout << prg << ": threshold error for argument " << ia.what() 
              << "; provided value must be a real number between 0 and 1" << std::endl;
    exit(EXIT_FAILURE);
  }

  return cytoscape_tc_display_threshold;
}

/* Check whether the provided path points to a .mol or .sdf file. */
bool is_valid_molsdf(const std::string & filename) {
  if(filename[0] == '.' || filename[0] == '#' || filename[filename.length() - 1] == '~') return false;
  std::string extension = filename.substr(filename.rfind("."));
  if(extension != ".mol" && extension != ".sdf") return false;
  return true;
}

/* Search the target directory for all .mol and .sdf files, allocate
 * Molecules out of them, and add the resulting Molecule *'s to the 
 * provided vector. The client is responsible for freeing the Molecule *'s
 * when they are finished using them. If a failure occurs, exit with an error
 * message attached to the program name prg. */
void fetch_directory_molecules(const char *prg, std::string dir, std::vector<Molecule *> & to_fill) {
DIR * target_dir = opendir(dir.c_str());

  if(target_dir == NULL) {
    std::cout << prg << ": unable to open target directory " << dir << std::endl;
    exit(EXIT_FAILURE);
  }

  const element_map *elements = symbols_to_weights();
  dirent * file_in_dir;
  while((file_in_dir = readdir(target_dir)) != NULL) {
    std::string filename = file_in_dir->d_name;
    if(is_valid_molsdf(filename)) {
      std::string chemname = filename.substr(0, filename.rfind(".")); // kill extension
      to_fill.push_back(new Molecule(dir + "/" + filename, chemname, elements, SuperFastHash));
    }
  }

  delete elements;
  closedir(target_dir);
}

#endif // ECFP_UTILS_H
