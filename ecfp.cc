#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdlib>

#include "molecule.h"
#include "stringhash.h"
#include "ecfp_utils.h"

/* Given a vector of molecules with their ecfp's generated, create a
 * matrix of Tanimoto similarity coefficients between pairs of molecules
 * and return it as a double **. The (i,j) entry of the matrix lists the
 * Tanimoto coefficient between molecules i and j. This function allocates
 * memory, so the client is responsible for freeing it after use. */
static double ** tanimoto_matrix(const std::vector<Molecule *> & molecules) {
  double ** result = new double *[molecules.size()];
  for(size_t i = 0; i < molecules.size(); i++) result[i] = new double [molecules.size()];

  for(size_t i = 0; i < molecules.size(); i++) {
    for(size_t j = i; j < molecules.size(); j++) {
      result[i][j] = result[j][i] = Molecule::tanimoto_coefficient(molecules[i], molecules[j]);
    }
  }

  return result;
}

/* Write to the file <dirname>_tanimoto.txt a matrix containing Tanimoto
 * coefficients between all pairs of molecules, as stored in the results
 * matrix provided. The top row and leftmost column are reserved for names
 * of the input molecules. The ith row and jth column after the header row
 * and column indicate the Tanimoto coefficient between molecules i and j. */
static void write_tanimoto_matrix(const std::string & dirname, double ** results, 
                                  const std::vector<Molecule *> & molecules) {
  std::ofstream output("output/" + dirname + "_tanimoto.txt");
  output << " ";
  for(const Molecule * m: molecules) output << "\t" << m->get_name();
  output << std::endl;

  for(size_t i = 0; i < molecules.size(); i++) {
    output << molecules[i]->get_name();
    for(size_t j = 0; j < molecules.size(); j++) {
      output << "\t" << (j >= i ? results[i][j] : results[j][i]);
    }
    output << std::endl;
  }

  output.close();
}

static void write_cytoscape_table(const std::string & dirname, const std::string & tag, double ** results, 
                                  const std::vector<Molecule *> & molecules, double threshold) {  
  std::ofstream output("output/" + dirname + "_cytoscape_g" + tag + ".txt");
  output << "source,target,tanimoto coefficient" << std::endl;
  
  for(size_t i = 0; i < molecules.size(); i++) {
    for(size_t j = i+1; j < molecules.size(); j++) {
      double tc = results[i][j];
     
      if(tc >= threshold) {        
        output << molecules[i]->get_name() << "," << molecules[j]->get_name() << ","
               << tc << std::endl;
      }     
    }
  }
  for(Molecule *m: molecules) output << m->get_name() << "," << m->get_name() << ","
                                     << "1.000" << std::endl;
  output.close();
}




int main(int argc, char *argv[]) { 
  if(argc != 4) {
    std::cout << "usage: " << argv[0] << " [niterations] [interaction threshold] [molecule directory]" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  int niterations = get_niterations(argv[0], argv[1]);
  double cytoscape_tc_display_threshold = get_display_threshold(argv[0], argv[2]);

  std::cout << "Making Molecules out of all .mol and .sdf files in directory " 
            << argv[3] << " ..." << std::endl;
  std::vector<Molecule *> all_molecules;
  fetch_directory_molecules(argv[0], argv[3], all_molecules); 
  std::cout << "Generated " << all_molecules.size() << " Molecules." << std::endl;

  std::cout << "Computing ECFP_" << 2 * niterations
            << " fingerprints for all structures..." << std::endl;
  for(Molecule * m: all_molecules) m->generate_ecfp(niterations);
  
  std::cout << "Computing Tanimoto coefficients for all pairs of structures..." << std::endl;
  double ** tanimotos = tanimoto_matrix(all_molecules);
  
  std::string dir_shortname = argv[3]; 
  while(dir_shortname.back() == '/') dir_shortname.pop_back(); // get rid of any slashes at end of dirname
  dir_shortname = dir_shortname.substr(dir_shortname.rfind("/") + 1); // name after last nonterminal slash
  std::cout << "Writing Tanimoto matrix to output/" << dir_shortname << "_tanimoto.txt ... " << std::endl;
  write_tanimoto_matrix(dir_shortname, tanimotos, all_molecules);
  std::cout << "Done." << std::endl;

  std::ostringstream decimal_repr;
  decimal_repr << std::fixed << std::setprecision(2) << cytoscape_tc_display_threshold;
  std::string tag = decimal_repr.str();
  tag = (cytoscape_tc_display_threshold >= 1.0 ? "100" : "0" + tag.substr(tag.find(".") + 1));
  std::cout << "Writing Cytoscape interaction table to output/" << dir_shortname
            << "_cytoscape_g" << tag << ".txt ..." << std::endl;
  write_cytoscape_table(dir_shortname, tag, tanimotos, all_molecules, cytoscape_tc_display_threshold);
  std::cout << "Done." << std::endl;

  for(Molecule * m: all_molecules) delete m;
  for(size_t i = 0; i < all_molecules.size(); i++) delete[] tanimotos[i];
  delete[] tanimotos;
  return 0;
}
