#include <atom.hpp>
#include <bond.hpp>

Bond::Bond (Atom* c, Atom* d)
{
   if (c > d)
   {
      a = d;
      b = c;
   }
   else
   {
      a = c;
      b = d;
   }
}

bool Bond::operator== (const Bond & other) const 
{
   return (a == other.a && b == other.b) || (a == other.b && b == other.a);
}

bool Bond::operator< (const Bond & other) const
{
   return a < other.a || ( b < other.b && a == other.a );
}

