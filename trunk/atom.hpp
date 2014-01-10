#ifndef ATOM_H
#define ATOM_H

#include <point.hpp>
#include <matrix3.hpp>
#include <fxn.hpp>
#include <vector>
#include <string>
#include <bond.hpp>

class Parameters;
class Bond;

class Atom
{
   private:
      int atomIndex; //Index of atom in atomList
      Point pos; //position relative to cellInfo.dim()
      std::vector<Atom*> neighbours; //Pointers to all atoms close to this atom

      void copy(const Atom&);
      
      static void bondOff(Bond); //Turns off a bond to a neighbour atom
      static void bondOn(Bond); //Turns on a bond to a neighbour atom
      void atomOff(); //Turns off all bonds to this atom

      static const double EPSILON = 1e-8;
   public:
      static std::vector<Atom> atomList; //List of all atoms in the cell
      static Matrix adjMatrix; //Adjacency Matrix of all atoms in the cell, as defined by their bonding
      static Parameters cellInfo; //vectors defining cell Atom is found in

      static void readAtoms(std::string); //Reads atoms in from an xcrysden file
      static void readMinAtoms(std::string); //Reads atoms in from a lammps data file
      static void multiplyCell(Point); //Increases a cell size of a factor in every direction
      static void connectAtoms(int); //Defines the neighbour pointers for all atoms in the cell
      static void outputAtoms(std::string); //Outputs the atoms, bonds, and angles to an lammps input file
      static void genBondList(); //Generates a list of all bonds from the atoms
      static void genAngleList(); //Generates a list of all angles from the bonds
      static void genAtomDelList(); //Generates a list of all atoms allowed to be turned off
      static void genBondDelList(unsigned int, double); //Generates a list of all bonds allowed to be deleted
      static void genBondDelList(unsigned int); //Verison without radial deletion condition
      static bool delRandBond(unsigned int); //Deletes an eligible random bond from all atoms
      static bool delRandBond(unsigned int, double); //Version with radial deletion condition
      static bool delRandAtom(); //Deletes an eligible atom
      static int delPercentBond(double, unsigned int); //Deletes a percentage of all bonds randomly
      static int delPercentBond(double, unsigned int, double); //Version with radial deletion condition
      static int delPercentAtom(double); //Deletes a percentage of all atoms randomly
      static int getNumCoordX(int); //Returns the number of atoms with a specific coordination number
      static double getCoordXNeighBond(int); //Returns the average coordination of bonded neighbours to atoms with a specific coordination
      static int getNumAtoms(); //Returns the number bonded within the adjacency matrix

      Atom(); //Intitalizes a null atom
      Atom(const Point&, int); //Intitalizes an atom with a position and index
      Atom(const Atom&); //Initializes a copy of an atom

      void setPos(const Point&); //Sets the position of an atom
      void setRelPos(const Point& other) {pos = other;} //Sets the position relative to the cell vectors
      void setNeighbour(Atom*); //Adds an atom to the neighbours vector
      void clearNeighbours(); //Removes all neighbours from an atom
      void setIndex(int); //Changes the index of an atom
      void setParam(Parameters); //Links to the cell parameters for all atoms
      void delNeighbour(int); //Removes a neighbour from an atom. For changing crystal structure, not bonding structure

      int getNumNeigh() const {return neighbours.size();} //Returns the number of neighbours of an atom
      double getRelPos(int i) const {return pos.coord(i);} //Returns the relative position in vector 0,1,2 of an atom's cell
      const Point& getRelPos() const {return pos;} //Returns the position of an atom relative to its cell vectors

      double getPos(char) const; //Returns the position in x,y,z of an atom
      double getPos(int) const; //Same as above, 0,1,2 input
      Point getPos() const; //Returns the position of an atom as a point
      Atom* getNeighbour(int) const; //Returns a pointer to an atom's neighbour
      int getNumBonds() const; //Returns the coordination of an atom
      int getIndex() const; //returns the index of the atom
      int getNeighbourIndex(int) const; //returns teh index of an atom's neighbour
      bool delCandidate() const; //Determines is an atom is a deletion candidate
      
      Atom& operator=(const Atom& other) {copy(other); return *this;}
      bool operator==(const Atom& other) const {return atomIndex == other.atomIndex;}
      bool operator<(const Atom& other) const {return atomIndex < other.atomIndex;}
      bool operator>(const Atom& other) const {return !(*this < other || *this == other);}
};

#endif
