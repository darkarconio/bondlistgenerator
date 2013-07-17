#include <angle.hpp>
#include <minexcept.hpp>

void Angle::set (Atom* a, Atom* b, Atom* c)
{
   bondList.clear();
   Bond one(a,b);
   Bond two(b,c);
   middle = b;
   if (one == two) throw BadStructureException();
   bondList.insert(one);
   bondList.insert(two);
}

bool Angle::operator== (const Angle & other) const
{
   return ( first() == other.first() && second() == other.second() ) || 
          ( first() == other.second() && second() == other.first() );
}

bool Angle::operator< (const Angle & other) const
{
   return first() < other.first() || ( second() < other.second() && first() == other.first() );
}

int Angle::operator[] (unsigned int i) const
{
   switch(i)
   {
      case 0: 
         return first().otherAtom(middle)->getIndex();
      case 1:
         return middle->getIndex();
      default:
         return second().otherAtom(middle)->getIndex();
   }
}
