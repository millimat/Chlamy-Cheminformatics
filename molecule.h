#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector> // for std::vector
#include <map> // for molecule weight map
#include <set> // for a feature's set of bonds
#include <unordered_set> // for std::unordered_set
#include <cstring> // std::strlen

#include "stringhash.h"

// Forward declarations
class Atom; // for use by Bond
class Feature; // for use by Atom

/********************************** CLASS: Bond **********************************/

enum BOND_TYPE { SINGLE = 1, DOUBLE = 2, TRIPLE = 3,
                 AROMATIC = 4, SINGLE_OR_DOUBLE = 5, 
                 SINGLE_OR_AROMATIC = 6, DOUBLE_OR_AROMATIC = 7,
                 ANY = 8 }; // matches formatting of bond line type field in CTAB
enum BOND_STEREO { NONE, UP = 1, DOWN = 6 }; // UP/DOWN match CTAB bond line stereo field

/**
 * Class: Bond 
 * -----------
 * A bond (edge) to another atom (vertex) in a Molecule (graph).
 */
class Bond {
 public:
  /* Initialize bond's private fields. */
  Bond(BOND_TYPE type, BOND_STEREO stereo, Atom * source, Atom * target);

  /* Destroy this Bond. */
  ~Bond();

  /* Getter methods */
  BOND_TYPE get_type() const { return this->type; }
  BOND_STEREO get_stereo() const { return this->stereo; }
  Atom * get_source() const { return this->source; }
  Atom * get_target() const { return this->target; }

  /* UnaryPredicate function type telling whether the argument Bond
   * targets a hydrogen Atom. */
  static bool targets_hydrogen(const Bond * b);

  /* Compare two bonds based on their bond order and their
   * target identifiers. Return true if Bond one is "less"
   * than Bond two. Used during the ECFP algorithm at the start
   * of each itearation for each Atom. */
  static bool pre_iteration_cmp(const Bond * one, const Bond * two);

  /* A functor allowing for the comparison of two Bond *'s by 
   * the identity of the Atom *'s they connect. The functor treats 
   * two Bond *'s as equivalent if they go between the same pair
   * of Atom *'s in any direction; hence, this functor can be used 
   * to store a set of Bond *'s as bidirectional edges. */
  struct connection_cmp { bool operator() (const Bond * one, const Bond * two); };

  /* A function that uses connection_cmp to determine whether
   * two const Bond *'s map between the same pair of Atom *'s. */
  static bool connection_eq(const Bond * one, const Bond * two); 

  /* Return the "opposite" of a bond's stereochemistry if it is UP/DOWN. */
  static BOND_STEREO reverse_stereo(BOND_STEREO bs);

 private:
  BOND_TYPE type;
  BOND_STEREO stereo;
  Atom * source;
  Atom * target;
};

/********************************** CLASS: Atom **********************************/

/**
 * struct: element
 * ---------------
 * A structure storing information about the
 * most common nuclide of an element on the 
 * periodic table. 
 */
struct element {
  int atomic_number;
  int atomic_mass;
};

using element_map = std::map<std::string, element>; // map symbols to nuclide info 
enum COLOR { WHITE, GRAY, BLACK }; // DFS visit statuses for Atoms

/**
 * class: Atom
 * -----------
 * A structure representing an atom (vertex) in a Molecule (undirected graph).
 * The atom is equipped with all of its basic chemical information (bonds to
 * other Atoms, atomic number, mass, charge, symbol, etc.) as well as some
 * additional fields to support the execution of the ECFP fingerprinting
 * algorithm.
 */
class Atom {
 public: 

  /* Construct a new atom using an atom line from a CTFile's
   * atom block and a pre-stored table of atomic weights. */
  Atom(const std::string & atomline, const element_map * a_weights);

  /* Destroy this atom. */
  ~Atom();

  /*** Getter methods ***/
  int get_position() const { return this->position; }
  const std::string & get_symbol() const { return this->symbol; }
  int get_atomic_number() const { return this->atomic_number; }
  int get_atomic_mass() const { return this->atomic_mass; }
  int get_charge() const { return this->charge; }
  int get_n_heavy_neighbors() const { return this->n_non_h_neighbors; }
  int get_heavy_valence() const { return this->heavy_valence; }
  int get_n_hydrogens() const { return this->n_h_neighbors; }
  bool get_in_ring() const { return this->in_ring; }
  const std::vector<Bond *> & get_bonds() const { return this->bonds; }
  bool get_color() const { return this->dfs_visited; }
  const Feature * get_feature() const { return this->feature; }

  /*** Field modification methods ***/
  
  /* Add a bond to a new neighboring atom during Molecule initialization. */
  void add_bond(BOND_TYPE type, BOND_STEREO stereo, Atom * target);

  /* Visit this Atom as part of a DFS to locate all rings in the parent Molecule. */
  void ring_dfs_visit(std::vector<Atom *> & grays, std::unordered_set<Atom *> & quicksearch_grays);

  /* Update this atom's number of hydrogens at the end of Molecule initialization. */
  void resolve_n_hydrogens();

  /* Delete all of this Atom's bonds to hydrogens, or delete every bond in 
   * the Atom if the Atom itself is a hydrogen. */
  void delete_bonds_to_hydrogens();

  /* Record this Atom's 1-indexed position in the parent Molecule's adjacency list. */
  void set_position(int pos) { this->position = pos; }

  /* Perform the "zeroth" iteration of the ECFP algorithm by initializing
   * this atom's first feature from its Daylight atomic invariants. */
  void iteration_zero();

  /* Perform a single iteration of the ECFP algorithm's feature-updating
   * step on this atom using the given hash function. */
  void perform_ecfp_iteration();

  /* Initialize a new ECFP feature from this atom's current list of iteration
   * information and return a const pointer to it. This should be done at the
   * end of an ECFP iteration, and the fields in the Atom's feature should 
   * not be used before then except to initialize iteration_info at the start
   * of a new ECFP round. */
  const Feature * generate_feature(string_hash h);

  /* Delete this atom's ECFP iteration and feature information. 
   * The Feature * it contains is not freed, since it is assumed
   * that a parent Molecule is still tracking that pointer and will 
   * free it upon its own destruction. */
  void clear_ecfp_info();

 private:
  int position; // atom's 1-indexed position in the CTAB (if hydrogens are ignored)

  std::string symbol; // atom's periodic table symbol
  int atomic_number;
  int atomic_mass; 
  int charge;
  int n_non_h_neighbors = 0; // number "heavy" neighboring atoms
  int heavy_valence = 0; // number of bonds to heavy atoms
  int n_h_neighbors = 0; // number of neighboring hydrogens
  bool in_ring = false; // is this atom part of a ring?
  std::vector<Bond *> bonds; // bonds to neighboring atoms
  
  COLOR dfs_visited = WHITE;
  Atom * ring_start = NULL; // if Atom in a ring, this is first ring node the DFS found
  
  std::vector<int> iteration_info; // used to update fingerprint during iteration

  /* The set of bonds and atoms represented by this atom's ECFP
   * feature during the current iteration. */
  std::set<Bond *, Bond::connection_cmp> curr_feature_bonds;
  std::set<Atom *> curr_feature_atoms;

  /* The atom's ECFP feature. Until generate_feature is called at 
   * the end of an iteration, the data from this feature represents
   * information from the previous ECFP iteration. */
  Feature * feature = NULL;

  /* Interpret the charge specifier field from the "ccc" field of an atom line
   * from a CTAB. */
  int interpret_charge_specifier(int spec);

  /* Add the sets of bonds and atoms represented by a neighboring 
   * Atom's feature (from the previous iteration) to this Atom's 
   * current sets of feature bonds and atoms. */
  void add_feature_subgraph_from(const Feature * other);
};

/********************************* CLASS: Feature ***********************************/

/**
 * Class: Feature
 * --------------------
 * A structure representing the contents of an ECFP feature,
 * the core piece of information used by the algorithm of 
 * Rogers and Hahn for molecular fingerprinting. 
 * Features store hashed identifiers based on information 
 * about the molecular substructures they represent, 
 * explicit representations of those substructures, and
 * records of the iteration numbers they came from during
 * the ECFP algorithm.
 */
class Feature {
 public:
  uint32_t identifier = 0; // this ECFP feature's hash-generated identifier
  int iteration = 0; // which iteration of the ECFP algorithm this feature was created in
  std::set<Bond *, Bond::connection_cmp> induced_bonds; // pointers to all Bonds this feature encompasses
  std::set<Atom *> induced_atoms;
  
  /* Check whether two Features represent the same subgraph. 
   * It is possible that two features represent the same edges but
   * have different combinations of known and unknown border atoms;
   * this method distinguishes between such features. Bonds are 
   * treated as equal if their source and target form the same set. */
  static bool same_subgraph(const Feature * one, const Feature * two);

  /* Check whether two Features have the same hashed identifier. */
  static bool same_id(const Feature * one, const Feature * two);

  /* Functor allowing the comparison of two Feature *'s first by their source iteration
   * number and then by their hashed identifiers. */
  struct it_id_cmp { bool operator() (const Feature * one, const Feature * two); };

};

/********************************** CLASS: Molecule **********************************/

/**
 * Class: Molecule
 * ---------------
 * A graph-type data structure representing a chemical compound, 
 * with atoms treated as the nodes and bonds treated as the edgess of a graph. 
 * Structures are generated by reading SDFs or Molfiles, and the information 
 * in the molecule is used to generate extended connectivity fingerprints
 * (ECFP's) using the algorithm of Rogers and Hahn.
 */
class Molecule {
 public:
  
  /* Construct a Molecule from a .sdf or .mol structural table file whose path is given
   * by sdf_path using the provided atomic weight table and hash function. */
  Molecule(const std::string & sdf_path, const std::string & name, 
           const element_map * a_weights, string_hash hash);

  /* Destroy this Molecule. */
  ~Molecule();

  /* Get this Molecule's name. */
  std::string get_name() const { return this->name; }

  /* List basic information about every Atom in this Molecule. */
  std::string atoms_report() const;

  /* Generate an ECFP_N using the algorithm of Rogers and Hahn. 
   * N is twice the number of iterations performed by the Molecule. */
  void generate_ecfp(int n_iterations);

  /* Fill the argument vector with a list of the integer identifiers
   * associated with this Molecule's ECFP in sorted order. */
  void fetch_ecfp(std::vector<uint32_t> & out) const;

  /* Delete this Molecule's current fingerprint and clear the feature
   * information for all of its underlying Atoms. */
  void clear_ecfp();

  /* Change this Molecule's fingerprinting hash function. Its existing
   * fingerprint information is cleared in the process. */
  void switch_hash(string_hash alt_hash);

  /* Generate the Tanimoto similarity coefficient between two Molecules
   * by comparing their ECFP's. */
  static double tanimoto_coefficient(const Molecule * one, const Molecule * two);

 private:
  std::string name; // A user-defined name for this molecule.
  std::vector<Atom *> atoms;
  std::set<const Feature *, Feature::it_id_cmp> ecfp; // This Molecule's fingerprint, with associated info.
  std::vector<const Feature *> discarded; // redundant ECFP features
  std::vector<uint32_t> fingerprint; // ECFP with identifiers only, sorted by ID value
  string_hash h; // The hash function this molecule uses.
  bool ecfp_ready = false; // Has this molecule generated its ECFP?

  /* Construct a pair of Bonds between two Atoms using a bond line 
   * from the bond block of a CTAB. */
  void init_bond(const std::string & line);

  /* Do a depth-first search over the entire graph this Molecule represents
   * to mark all Atoms that are part of rings (cycles). */
  void ring_dfs();

  /* Destroy all Atoms representing explicit hydrogens, as well as
   * all bonds they are connected to. */
  void remove_all_hydrogens();
 
  /* Insert the parameter Feature * into this Molecule's EFCP,
   * with special accommodations to avoid fingerprints representing the
   * same chemical feature (as outlined by Rogers and Hahn) and to avoid
   * storage of two fingerprints with the same hashed identifier. */
  void add_to_fingerprint(const Feature * candidate);
 
};

#endif // MOLECULE_H
