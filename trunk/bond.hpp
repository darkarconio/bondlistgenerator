#ifndef BOND_H
#define BOND_H

#include <atom.hpp>

class Atom;
class Angle;

class Bond
{
   private:
      const Atom * a;
      const Atom * b;
      static const int MIN_BOND = 2;
      static const int MIN_NEIGH_BOND = 3;
      static const int FULL_BOND = 4;
      void copy (const Bond&);
      bool guideCriteria(unsigned int) const;
   protected:
      friend class Angle;
      const Atom * otherAtom (const Atom * c) const {if (c == a) return b; else return a;}
   public:
      Bond (const Atom*, const Atom*);
      Bond (const Bond& other) {copy(other);}
      bool operator== (const Bond&) const;
      bool operator< (const Bond&) const;
      int operator[] (unsigned int) const;
      Bond& operator= (const Bond& other) { copy(other); return *this; }
      bool offCandidate(unsigned int) const;
      bool offCandidate(unsigned int, double) const;
      Point location() const;
};

#endif
