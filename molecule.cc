#include "molecule.h" 
#include <fstream> // std::ifstream
#include <sstream> // std::ostringstream
#include <algorithm> // std::min, max, sort, remove_if, erase, equal, etc.
#include <utility> // std::pair, std::swap

#include <cstdio> // sscanf
#include <cassert> // assert

/************************************ Bond definitions ************************************/

/* Initialize the private fields of this Bond using the parameter information. */
Bond::Bond(BOND_TYPE type, BOND_STEREO stereo, Atom * source, Atom * target) {
  this->type = type;
  this->stereo = stereo;
  this->source = source;
  this->target = target;
}

/* Destructor for this Bond. Since a Bond allocates no memory, it does nothing. */
Bond::~Bond() { }

/* UnaryPredicate function type telling whether this Bond targets a 
 * hydrogen Atom. Used to std::remove_if a heavy atom's bonds to hydrogens. */
bool Bond::targets_hydrogen(const Bond * b) {
  return b->get_target()->get_symbol() == "H";
}

/* Compare two bonds based on their bond order (first) and their target
 * identifiers (second). As above, this is done using std::pairs.
 * This function is used to sort a std::vector of Bond *'s,
 * so the arguments are Bond &'s. */
bool Bond::pre_iteration_cmp(const Bond * one, const Bond * two) {
  std::pair<uint32_t, uint32_t> one_pair(one->type, one->target->get_feature()->identifier);
  std::pair<uint32_t, uint32_t> two_pair(two->type, two->target->get_feature()->identifier);
  return one_pair < two_pair;
}

/* Comparison of two Bond *'s by the identify of the Atom *'s they connect, for use
 * in the connection_cmp functor. The Atom *'s are placed in sorted std::pairs
 * as numbers, and these pairs are compared. Hence, the () operator on this
 * functor will reflexively return false for two Bond *'s (i.e. two Bond *'s
 * are considered equal) if and only if they represent the same chemical bond
 * in the containing Molecule. */
bool Bond::connection_cmp::operator() (const Bond * one, const Bond * two) {
  std::pair<intptr_t, intptr_t> one_pair((intptr_t)one->get_source(), (intptr_t)one->get_target());
  std::pair<intptr_t, intptr_t> two_pair((intptr_t)two->get_source(), (intptr_t)two->get_target());
  if(one_pair.first > one_pair.second) std::swap(one_pair.first, one_pair.second);
  if(two_pair.first > two_pair.second) std::swap(two_pair.first, two_pair.second);
  return one_pair < two_pair;
}

/* A function that uses connection_cmp to determine whether two const Bond *'s map
 * between the same pair of Atom *'s. Since connection_cmp::operator() reflexively
 * returns false for two Bond *'s if and only if they map between the same two Atom *'s,
 * it suffices to use two calls to that operator to determine equality. */
bool Bond::connection_eq(const Bond * one, const Bond * two) {
  connection_cmp cc;
  return !cc(one, two) && !cc(two, one);
}

/* Return the "opposite" of a bond's stereochemistry if it is UP/DOWN. */
BOND_STEREO Bond::reverse_stereo(BOND_STEREO bs) {
  if(bs == UP) return DOWN;
  if(bs == DOWN) return UP;
  return NONE;
}

/************************************ Atom definitions ************************************/

static const std::string kAtomLineRelevantInfoFormat = "xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddccc";
static const int kAtomLineSymbolPos = 31; // symbol appears at pos 31 of atom line
static const int kAtomLineSymbolLen = 3;
static const int kAtomLineMassDiffPos = 34;
static const int kAtomLineMassDiffLen = 2;
static const int kAtomLineChargePos = 36;
static const int kAtomLineChargeLen = 3;

/* Construct a new Atom by reading an atom line from a CTAB in a .mol
 * or .sdf. The line information is necessary to set the Atom's symbol,
 * atomic number, atomic mass, and charge. */
Atom::Atom(const std::string & atomline, const element_map * a_weights) {
  char buff[kAtomLineSymbolLen + 1] = { 0 };
  sscanf(atomline.c_str() + kAtomLineSymbolPos, " %3[^ ] ", buff); 
  symbol = buff;

  int mass_diff = std::stoi(atomline.substr(kAtomLineMassDiffPos, kAtomLineMassDiffLen));
  int charge_specifier = std::stoi(atomline.substr(kAtomLineChargePos, kAtomLineChargeLen));
  
  atomic_number = a_weights->at(symbol).atomic_number;
  atomic_mass = a_weights->at(symbol).atomic_mass + mass_diff;
  charge = interpret_charge_specifier(charge_specifier);
}
 
/* Interpret the charge specifier field from the "ccc" field of an atom line
 * from a CTAB. The charge specifier must fall between 0 and 7 (inclusive),
 * with the meanings as indicated in the comments below. */
int Atom::interpret_charge_specifier(int spec) {
  assert(spec >= 0 && spec <= 7);
  int result;

  switch(spec) {
  case 0: result = 0; break; // charge +0
  case 1: result = 3; break; // charge +3
  case 2: result = 2; break; // charge +2
  case 3: result = 1; break; // charge +1
  case 4: result = 0; break; // doublet radical; interpret as charge +0
  case 5: result = -1; break; // charge -1
  case 6: result = -2; break; // charge -2
  case 7: result = -3; break; // charge -3
  }
  
  return result;
}

/* Destructor for the Atom class. 
 * An Atom allocates memory both for its Feature fingerprints
 * and for its Bonds. However, memory management of fingerprints is 
 * passed on to the parent Molecule, so an Atom only needs to deallocate
 * its Bonds as it is destroyed. */
Atom::~Atom() {
  for(Bond * b: bonds) delete b;
}

/* Add a Bond * to another atom to this Atom's bonds vector during
 * Molecule construction. A call to this method on one Atom should
 * always be accompanied by a symmetric call on the other atom,
 * with the bond stereochemistry reversed if it is an up/down
 * single bond. If the new neighbor is hydrogen, increment
 * n_h_neighbors; otherwise, increment n_non_h_neighbors and
 * increase heavy_valence by the order of the bond to the new
 * neighbor. */
void Atom::add_bond(BOND_TYPE type, BOND_STEREO stereo, Atom * target) {
  bonds.push_back(new Bond(type, stereo, this, target));

  // TODO: What if bond is programmed as aromatic? Use assert for now
  assert(type != AROMATIC && type != SINGLE_OR_DOUBLE 
         && type != SINGLE_OR_AROMATIC && type != DOUBLE_OR_AROMATIC
         && type != ANY);
  
  if(target->get_symbol() == "H") n_h_neighbors++;
  else {
    n_non_h_neighbors++;
    heavy_valence += (int)type;
  }
}

/* Subroutine of a molecular depth-first search to locate all cycles (rings) in a molecule.
 * If x is the active node in a depth-first search, and y is a neighbor of x, then the edge
 * x->y must be part of a cycle if y is gray or if a gray node is reachable from y. This
 * is because x was reached from some node z that is currently gray, meaning the path
 * z -> (...) x -> y -> (...) z is a cycle. Conversely, every cycle in a given connected
 * component of an undirected graph must be detectable in this manner: since a depth-first
 * search traverses every edge in the component, every cycle in the component must experience
 * some point where the current node is either adjacent to the start of the cycle or adjacent
 * to a node that can reach the start of the cycle but was visited on an earlier branch.
 * In either case, every previous node in the cycle must be gray.
 * 
 * Hence, we can mark all nodes in a cycle as follows: Let x be the current node, and y be
 * a neighbor of x. If y is gray, mark x and all its predecessors going back until y as belonging
 * to a cycle, and store x as the beginning of that cycle. If y is black, check whether it 
 * belonged to a cycle, and if so, whether the start of that cycle, z, is currently gray; if so,
 * mark x and all its predecessors back to z as belonging to a cycle starting at z. 
 * This method replicates the procedure, with one exception: cycles of length two (e.g. x->y->x)
 * are not considered to be rings in chemistry, so they are not marked as cycles. */
void Atom::ring_dfs_visit(std::vector<Atom *> & grays, std::unordered_set<Atom *> & quicksearch_grays) {
  Atom * predecessor = (grays.empty() ? NULL : grays.back());
  dfs_visited = GRAY; // currently being considered
  grays.push_back(this);
  quicksearch_grays.insert(this);

  for(Bond * b: bonds) {
    Atom * neighbor = b->get_target();
    if(neighbor != predecessor) { // Ignore the edge that immediately goes back to where we came from
      if(neighbor->dfs_visited == WHITE) neighbor->ring_dfs_visit(grays, quicksearch_grays);
      else {
        bool cycle_found = false;
        if(neighbor->dfs_visited == GRAY) { 
          // The edge to the neighbor closes a ring, so we have found a cycle
          cycle_found = neighbor->in_ring = true;
          neighbor->ring_start = neighbor;
        } else if(neighbor->dfs_visited == BLACK && neighbor->in_ring
                  && quicksearch_grays.count(neighbor->ring_start)) {
          // The neighbor is part of a ring and a gray node is reachable from it, so we have found cycle
          cycle_found = true;
        }
        
        if(cycle_found) { // Move back marking nodes as in cycles till we hit the cycle start
          for(std::vector<Atom *>::reverse_iterator it = grays.rbegin(); // current Atom
              *it != neighbor->ring_start; it++) {
            (*it)->in_ring = true;
            (*it)->ring_start = neighbor->ring_start;
          }
        }
      }
    }
  }

  dfs_visited = BLACK; // visit finished
  grays.pop_back();
  quicksearch_grays.erase(this);
}

/* Edit the number of hydrogens said to be attached to certain 
 * p-block nonmetals in order to properly initialize their number
 * of H's to what would be expected in an organic context.
 * The elements subject to editing are B, C, N, O, F, 
 * non-hypervalent Si, non-hypervalent P, non-hypervalent S, Cl, 
 * Br, and I. The actual number of hydrogen attachments is 
 * assumed to be the atom's "default" coordination number
 * (e.g. 4 for carbon) minus the number of heavy atoms attached 
 * to it during explicit bond tabulation plus some function of the
 * atom's charge. 
 *
 * - For boron, we subtract the charge, since boron is rarely positive
 * in an organic context and is usually negative only when tetracoordinate.
 * - For carbon, we subtract the absolute value of the charge, since carbon
 * is generally positive only as a carbocation and negative as a carbanion,
 * and in both situations has fewer bonds than its default coordination.
 * - For silicon, we assume the same rules as for carbon.
 * - For N, O, F, Cl, Br, I, and non-hypervalent P/S, we add the value of
 * the charge, since these atoms tend to be positive when they have
 * more bonds than their default coordination and negative when they have
 * fewer. 
 * 
 * Since the behavior of P/S is difficult to predict when they exist in
 * contexts with >4 bonds--and these atoms are rarely hypervalent
 * when attached to hydrogens in an organic context--we do not attempt
 * to modify the number of hydrogens in those cases. */
void Atom::resolve_n_hydrogens() {
  int final_hydrogens;
  int current_coord = heavy_valence + n_h_neighbors;

  if(symbol == "B") { // default coord 3
    final_hydrogens = 3 - heavy_valence - charge;
  } else if(symbol == "C" || (symbol == "Si" && current_coord <= 4)) { // 4
    final_hydrogens = 4 - heavy_valence - abs(charge);
  } else if(symbol == "N" || (symbol == "P" && current_coord <= 4)) { // 3
    final_hydrogens = 3 - heavy_valence + charge;
  } else if(symbol == "O" || (symbol == "S" && current_coord <= 4)) { // 2
    final_hydrogens = 2 - heavy_valence + charge;
  } else if (symbol == "F" || symbol == "Cl" || symbol == "Br" || symbol == "I") { // 1
    final_hydrogens = 1 - heavy_valence + charge;
  } else return; // don't do anything for any other atoms

  // If for some reason the calculation gave us < 0 hydrogens, just set 0
  n_h_neighbors = std::max(0, final_hydrogens);
}

/* Delete all of this Atom's bonds to hydrogens. */
void Atom::delete_bonds_to_hydrogens() {
  std::vector<Bond *> to_delete;
  for(Bond * b: bonds) { if(Bond::targets_hydrogen(b)) to_delete.push_back(b); }
  bonds.erase(std::remove_if(bonds.begin(), bonds.end(), Bond::targets_hydrogen), bonds.end());
  for(Bond * b: to_delete) delete b;
}

/* Perform the "zeroth" iteration of the ECFP algorithm by initializing
 * this atom with its Daylight atomic invariants.
 * The invariants -- number of bonds to heavy neighbors, number of hydrogen neighbors,
 * atomic number, atomic mass, charge, and ring status -- are put in 
 * a vector. This vector can be hashed to produce a feature using
 * generate_feature. Finally, to set up subsequent iterations of the ECFP
 * algorithm, curr_feature_bonds is initialized with this Atom's connections
 * to other Atoms and curr_feature_atoms is initialized with a pointer to this one. */
void Atom::iteration_zero() {
  iteration_info.push_back(heavy_valence);
  iteration_info.push_back(n_h_neighbors);
  iteration_info.push_back(atomic_number);
  iteration_info.push_back(atomic_mass);
  iteration_info.push_back(charge);
  iteration_info.push_back(in_ring); // as integer

  for(Bond * b: bonds) curr_feature_bonds.insert(b);
  curr_feature_atoms.insert(this);
}


/* Perform a single iteration of Roger and Hahn's ECFP algorithm:
 * 1. Append the current iteration number and the current Atom's
 *    previous-round identifier to its iteration_info array.
 * 2. Sort the Atom's bond vector by the Bonds' orders and their
 *    targets' previous-round identifiers (ensures algorithm is
 *    deterministic)
 * 3. For each Bond, append its order and its target's previous-round
 *    identifier to this Atom's iteration_info array. 
 * 4. Expand the sets of Bonds and Atoms this Atom's feature currently 
 *    represents to include the analogous sets from all its neighbors' 
 *    previous iterations.
 *
 * After this method has been executed by all Atoms, they will all hash
 * their iteration_info vectors into identifiers that will be appended
 * to the Molecule's ECFP if their representative bond sets are not 
 * redundant. */
void Atom::perform_ecfp_iteration() {
  assert(iteration_info.empty());

  iteration_info.push_back(feature->iteration + 1);
  iteration_info.push_back(feature->identifier);
 
  std::sort(bonds.begin(), bonds.end(), Bond::pre_iteration_cmp);

  for(Bond * b: bonds) {
    iteration_info.push_back(b->get_type());
    const Feature * neighbor_prev_feature = b->get_target()->feature;
    iteration_info.push_back(neighbor_prev_feature->identifier);
    add_feature_subgraph_from(neighbor_prev_feature);
  }
}

/* Add the sets of bonds and atoms represented by a neighboring Atom's
 * feature (from the previous iteration) to this Atom's current
 * sets of feature bonds/atoms. */
void Atom::add_feature_subgraph_from(const Feature * other) { 
  curr_feature_bonds.insert(other->induced_bonds.cbegin(), other->induced_bonds.cend());
  curr_feature_atoms.insert(other->induced_atoms.cbegin(), other->induced_atoms.cend());
}
 
/* Initialize a new ECFP feature from this atom's current list of iteration
 * information and return a const pointer to it. This should be done at the
 * end of an ECFP iteration, after iteration_info has been completely finished
 * for all atoms in the Molecule. Thus an Atom's feature for a given iteration
 * is created only at the end of that iteration, so that neighboring atoms
 * expanding their iteration_info during the main run still have access
 * to this Atom's feature from the previous iteration during the main run. */
const Feature * Atom::generate_feature(string_hash h) {
  Feature * next_feature = new Feature;
  next_feature->iteration = (this->feature == NULL ? 0 : this->feature->iteration + 1);
  next_feature->identifier = vector_hash(iteration_info, h);
  next_feature->induced_bonds = curr_feature_bonds;
  next_feature->induced_atoms = curr_feature_atoms;
  this->feature = next_feature; // old feature still reachable via parent Molecule
  iteration_info.clear();

  return next_feature;
}

/* Delete this atom's ECFP iteration and feature information. 
 * The Feature * it contains is not freed, since it is assumed
 * that a parent Molecule is still tracking that pointer and will 
 * free it upon its own destruction. */
void Atom::clear_ecfp_info() {
  iteration_info.clear();
  curr_feature_bonds.clear();
  curr_feature_atoms.clear();
  feature = NULL;
}

/*********************************** Feature definitions ***********************************/

/* Check whether two Features represent the same subgraph. 
 * It is possible that two features represent the same edges but
 * have different combinations of known and unknown border atoms;
 * this method distinguishes between such features by comparing
 * both their Atom sets and Bond sets (again, Bonds are compared
 * using the Atom *'s they connect, without respect for direction). */
bool Feature::same_subgraph(const Feature * one, const Feature * two) {
  return(one->induced_bonds.size() == two->induced_bonds.size() && 
         one->induced_atoms.size() == two->induced_atoms.size() &&
         std::equal(one->induced_bonds.begin(), one->induced_bonds.end(),
                    two->induced_bonds.begin(), Bond::connection_eq) &&
         std::equal(one->induced_atoms.begin(), one->induced_atoms.end(),
                    two->induced_atoms.begin()));  
}

/* Check whether two Features have the same hashed identifier. */
bool Feature::same_id(const Feature * one, const Feature * two) {
  return one->identifier == two->identifier;
}

/* Compare two Feature *'s first by their source iteration number and 
 * then by their hashed identifiers. */
bool Feature::it_id_cmp::operator() (const Feature * one, const Feature * two) {
  if(one->iteration != two->iteration) return one->iteration < two->iteration;
  return one->identifier < two->identifier;
}


/********************************** Molecule definitions **********************************/

static const int kNHeaderLines = 3;

static const std::string kCountsLineRelevantInfoFormat = "aaabbb";
static const int kCountsLineFieldLen = 3;
static const int kCountsLineNAtomsPos = 0;
static const int kCountsLineNBondsPos = 3;

static const std::string kBondLineRelevantInfoFormat = "111222tttsss";
static const int kBondLineFieldLen = 3;
static const int kBondLineAtom1Pos = 0;
static const int kBondLineAtom2Pos = 3;
static const int kBondLineTypePos = 6;
static const int kBondLineStereoPos = 9;

/* Use a CTAB in an SDF or Molfile to create a Molecule and its underlying set of
 * Atoms and Bonds. A CTAB starts with three header lines that can be ignored,
 * followed by a "counts line" listing the number of explicit atoms and bonds 
 * in the molecule. There is then a block of "atom lines" that can be used to 
 * construct Atoms, followed by a block of "bond lines" that allow for the construction
 * of bonds between Atoms. Once the graph has been constructed, a few scans over the
 * molecule are required in order to finalize the number of hydrogen attachments each
 * atom actually has, remove all explicit hydrogens (which do not contribute to an ECFP), 
 * and mark all atoms that are found in rings. The Molecule also stores the hash
 * function passed into its constructor for use in the actual ECFP algorithm. */
Molecule::Molecule(const std::string & sdf_path, const std::string & name,
                   const element_map * a_weights, string_hash hash) {
  std::ifstream sdf(sdf_path);
  assert(sdf.good());
  this->name = name;

  std::string line;

  // Skip header lines
  for(int i = 0; i < kNHeaderLines; i++) std::getline(sdf, line);

  // Read counts line for number of atoms and number of bonds
  std::getline(sdf, line);
  int n_atoms = stoi(line.substr(kCountsLineNAtomsPos, kCountsLineFieldLen)); 
  int n_bonds = stoi(line.substr(kCountsLineNBondsPos, kCountsLineFieldLen));

  // Add Atoms to this Molecule's adjacency list using the atom block
  for(int i = 0; i < n_atoms; i++) {
    std::getline(sdf, line);
    atoms.push_back(new Atom(line, a_weights));
  }
  
  // Use the bond block to create bonds between atoms
  for(int i = 0; i < n_bonds; i++) {
    std::getline(sdf, line);
    init_bond(line);
  }

  // With bonding information completed, finalize the number of hydrogen
  // neighbors each atom has, then remove all explicit hydrogens
  // and finalize the position of each atom in the adjacency list.
  for(int i = 0; i < n_atoms; i++) atoms[i]->resolve_n_hydrogens();
  remove_all_hydrogens();
  for(size_t i = 0; i < atoms.size(); i++) atoms[i]->set_position(i+1);

  // Search the atoms of the molecule to mark all atoms belonging to rings.
  ring_dfs();

  h = hash;
  sdf.close();
}

/* Construct a pair of Bonds between two Atoms using a bond line 
 * from the bond block of a CTAB. The (1-indexed) positions of the
 * two atoms are read from the bond line, followed by the bond's
 * order and its stereochemistry. This information is used to initialize
 * a one-directional bond from the first atom listed to the second,
 * then another one-directional bond from the second to the first, 
 * so that both atoms are now directly connected. */
void Molecule::init_bond(const std::string & line) {
  int atom1 = stoi(line.substr(kBondLineAtom1Pos, kBondLineFieldLen));
  int atom2 = stoi(line.substr(kBondLineAtom2Pos, kBondLineFieldLen));
  BOND_TYPE bt = (BOND_TYPE)stoi(line.substr(kBondLineTypePos, kBondLineFieldLen));

  BOND_STEREO one_to_two;
  if(bt == SINGLE) {
    one_to_two = (BOND_STEREO)stoi(line.substr(kBondLineStereoPos, kBondLineFieldLen));
  } else one_to_two = NONE;
  
  atoms[atom1 - 1]->add_bond(bt, one_to_two, atoms[atom2 - 1]);
  atoms[atom2 - 1]->add_bond(bt, Bond::reverse_stereo(one_to_two), atoms[atom1 - 1]);    
}


void Molecule::ring_dfs() {
  std::vector<Atom *> grays;
  std::unordered_set<Atom *> quicksearch_grays;

  for(Atom * a: atoms) {
    if(a->get_color() == WHITE) a->ring_dfs_visit(grays, quicksearch_grays);
  }
}

/* UnaryPredicate function telling whether this Atom is a hydrogen. */
static bool is_hydrogen(const Atom * a) {
  return a->get_symbol() == "H";
}

/* Delete all hydrogen atoms (and all bonds connected to them) from
 * this Molecule. Atoms that are not hydrogen will need to manually 
 * deallocate their bonds to hydrogens via delete_bonds_to_hydrogens,
 * while hydrogen atoms will invoke their destructors and thus
 * delete their bonds automatically. */
void Molecule::remove_all_hydrogens() {
  std::vector<Atom *> hydrogens_to_delete;
  for(Atom * a: atoms) {
    if(is_hydrogen(a)) hydrogens_to_delete.push_back(a);
    else a->delete_bonds_to_hydrogens();
  }

  atoms.erase(std::remove_if(atoms.begin(), atoms.end(), is_hydrogen), atoms.end());
  for(Atom * h: hydrogens_to_delete) delete h;
}

/* Destroy this Molecule by deleting all of the dynamically allocated
 * fingerprint features stored in its ECFP and all the Atoms in its
 * adjacency list. */
Molecule::~Molecule() {
  for(const Feature * ft: ecfp) delete ft;
  for(Atom * a: atoms) delete a;
}

/* Iterate over the Atoms of this molecule. For each one, report its
 * position in the adjacency list, atomic number, symbol, mass,
 * charge, and basic information about its neighbors and bonding. 
 * Store all the results in a string. */
std::string Molecule::atoms_report() const {
  std::ostringstream oss;
  oss << atoms.size() << " ATOMS: \n\n";

  for(size_t i = 0; i < atoms.size(); i++) {
    const Atom * a = atoms[i];
    oss << "Atom " << i+1 << ": " << a->get_atomic_number() << " " << a->get_symbol()
        << ", mass = " << a->get_atomic_mass() << ", charge = " << a->get_charge()
        << ".\n" << a->get_n_heavy_neighbors() << " heavy neighbor(s), " 
        << a->get_heavy_valence() << " bond(s) to heavy neighbors, "
        << a->get_n_hydrogens() << " H neighbor(s). In ring: " << (a->get_in_ring() ? "Y" : "N")
        <<  ".\nList of heavy neighbors: \n";

    for(const Bond * b: a->get_bonds()) {
      oss << "  - Atom " << b->get_target()->get_position()
          << " (" << b->get_target()->get_symbol() << "),"
          << " bond order = " << (int)b->get_type() << "\n";
    }
    oss << "\n";
  }

  return oss.str();
}

/* Generate an ECFP_N using the algorithm of Rogers and Hahn. 
 * N is twice the number of iterations performed by the Molecule. 
 * The algorithm generates an initial set of fingerprint features
 * from the Daylight atomic invariants of the Molecule's atom set
 * (iteration_zero) and then enters an iterative phase in which
 * each atom's existing fingerprint information and the previous
 * fingerprints of its neighbors are fed forward to generate
 * a new fingerprint encompassing a larger portion of the Molecule
 * centered at that atom. */
void Molecule::generate_ecfp(int n_iterations) {
  for(Atom * a: atoms) {
    a->iteration_zero();
    add_to_fingerprint(a->generate_feature(h));
  }

  // Info collection and fp generation must be two separate loops 
  // to avoid overwriting old fingerprint data
  for(int i = 1; i <= n_iterations; i++) { 
    for(Atom * a: atoms) a->perform_ecfp_iteration();
    for(Atom * a: atoms) add_to_fingerprint(a->generate_feature(h));
  }

  // generate numeric fingerprint sorted by ID values
  for(const Feature * f: ecfp) fingerprint.push_back(f->identifier);
  std::sort(fingerprint.begin(), fingerprint.end());

  for(const Feature * f: discarded) delete f;

  discarded.clear();
  ecfp_ready = true;
}

/* Attempt to insert the parameter Feature * into this Molecule's 
 * ECFP. Precautions are taken to ensure that the fingerprint never
 * contains two features with the same hashed identifier or the same
 * representative bond set:
 * - If a feature is found to represent the same bond set as one
 * already in the fingerprint, then the feature with the smaller iteration
 * number is discarded or, if they came from the same iteration, the
 * feature with the larger hashed identifier is discarded.
 * - If the feature to insert survives the check above, but its hashed
 * identifier duplicates that of some feature already in the set, the old
 * feature is discarded. 
 *
 * Since an ECFP feature that is rejected still needs to be used for the next
 * iteration of the algorithm, we do not delete rejected features right away,
 * but instead queue them to be discarded at the end of the ECFP algorithm.
 */
void Molecule::add_to_fingerprint(const Feature * candidate) {
  const Feature * graph_duplicate = NULL, * id_duplicate = NULL;
  for(const Feature * f: ecfp) {
    if(Feature::same_subgraph(candidate, f)) graph_duplicate = f;
    if(Feature::same_id(candidate, f)) id_duplicate = f;
  }

  if(graph_duplicate != NULL) {
    if(graph_duplicate->iteration < candidate->iteration || 
       graph_duplicate->identifier > candidate->identifier) {
      ecfp.erase(graph_duplicate);
      discarded.push_back(graph_duplicate);
    } else {
      discarded.push_back(candidate);
      return;
    }
  }

  if(id_duplicate != NULL && id_duplicate != graph_duplicate) {
    ecfp.erase(id_duplicate);
    discarded.push_back(id_duplicate);
  }
  
  // Since candidate came from high iteration #, insert with hint to end
  ecfp.insert(ecfp.cend(), candidate);
}

/* Fill the argument vector with a copy of this Molecule's 
 * numeric fingerprint. */
void Molecule::fetch_ecfp(std::vector<uint32_t> & out) const {
  assert(ecfp_ready);
  for(const uint32_t id: fingerprint) out.push_back(id);
}

/* Delete this Molecule's current fingerprint and clear the feature
 * information for all of its underlying Atoms. */
void Molecule::clear_ecfp() {
  for(const Feature * ft: ecfp) delete ft;
  for(Atom * a: atoms) a->clear_ecfp_info();
  fingerprint.clear();
  ecfp_ready = false;
}

/* Change this Molecule's fingerprinting hash function. Its existing
 * fingerprint information is cleared in the process. */
void Molecule::switch_hash(string_hash alt_hash) {
  clear_ecfp();
  h = alt_hash;
}

/* Generate the Tanimoto similarity coefficient between two Molecules
 * by comparing their ECFP's. Let X and Y be molecules with ECFP 
 * feature sets A and B, respectively. The Tanimoto coefficient
 * between X and Y is defined to be the size of the intersection of
 * A and B divided by the size of the union of A and B; i.e. the number 
 * of features the molecules share divided by the total number of distinct
 * features between the two molecules. */
double Molecule::tanimoto_coefficient(const Molecule * one, const Molecule * two) {
  assert(!(one->fingerprint.empty()) && !(two->fingerprint.empty()));
  int intersection_size = 0; 

  std::vector<uint32_t>::const_iterator it1 = one->fingerprint.cbegin(), it2 = two->fingerprint.cbegin();
  while(it1 < one->fingerprint.cend() && it2 < two->fingerprint.cend()) {
    if(*it1 < *it2) it1++;
    else if(*it2 < *it1) it2++;
    else { intersection_size++; it1++; it2++; }
  }

  int union_size = one->fingerprint.size() + two->fingerprint.size() - intersection_size;
  return (double)intersection_size / union_size;
}

