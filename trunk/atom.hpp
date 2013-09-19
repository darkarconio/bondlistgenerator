#ifndef ATOM_H
#define ATOM_H

#include <point.hpp>
#include <matrix3.hpp>
#include <fxn.hpp>
#include <vector>
#include <string>

class Parameters;

class Atom
{
   private:
      int atomIndex;
      Point pos; //position relative to cellInfo.dim()
      std::vector<Atom*> neighbours;
      void copy(const Atom&);
      static const unsigned int MIN_BOND = 2;
   public:
      static std::vector<Atom> atomList;
      static Matrix adjMatrix;
      static Parameters cellInfo; //vectors defining cell Atom is found in

      static void readAtoms(std::string);
      static void multiplyCell(Point);
      static void connectAtoms(int);
      static void outputAtoms(std::string);
      static void genBondList();
      static void genAngleList();
      static void genDelList();

      Atom();
      Atom(const Point&, int);
      Atom(const Atom&);

      void setPos(const Point&);
      void setRelPos(const Point& other) {pos = other;}
      void setNeighbour(Atom*);
      void delNeighbour(int);
      void clearNeighbours();
      bool delRandNeighbour();
      void setIndex(int);
      void setParam(Parameters);
      
      int getNumNeigh() const {return neighbours.size();}
      double getPos(char) const;
      double getPos(int) const;
      Point getPos() const;
      double getRelPos(int i) const {return pos.coord(i);}
      const Point& getRelPos() const {return pos;}
      Atom* getNeighbour(int) const;
      int getIndex() const;
      int getNeighbourIndex(int) const;
      Atom& operator=(const Atom& other) {copy(other); return *this;}
      bool operator==(const Atom& other) const {return atomIndex == other.atomIndex;}
      bool operator<(const Atom& other) const {return atomIndex < other.atomIndex;}
      bool operator>(const Atom& other) const {return !(*this < other || *this == other);}
};

#endif
