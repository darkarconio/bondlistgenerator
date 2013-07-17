#ifndef ANGLE_H
#define ANGLE_H

#include <set>
#include <bond.hpp>

class Angle
{
   private:
      std::set<Bond> bondList;
      Atom * middle;
   public:
      Angle (Atom* a, Atom* b, Atom* c) {set(a,b,c);}
      void set (Atom*, Atom*, Atom*);
      bool operator== (const Angle&) const;
      bool operator< (const Angle&) const;
      int operator[] (unsigned int) const;
      const Bond first() const {return *bondList.begin();}
      const Bond second() const {return *(++bondList.begin());}
};

#endif
