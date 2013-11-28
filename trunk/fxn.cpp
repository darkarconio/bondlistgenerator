#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>

#include <ctime>
#include <cstdlib>

#include <point.hpp>
#include <fxn.hpp>
#include <atom.hpp>
#include <bond.hpp>
#include <angle.hpp>
#include <minexcept.hpp>

using namespace std;

#define debug(x) cout << __LINE__ << ' ' << x << endl; cout.flush();

//Checks if two Points are bonded across cell walls, returns the corrected vector difference between them if they are
//Otherwise, simply returns the difference between them
Point Parameters::getRealDiff (Point & a, Point & b)
{
   Point c;
   c = (a-b);
   double sign;
   double xDel = c.scalarProj(dim(0));
   double yDel = c.scalarProj(dim(1));
   double zDel = c.scalarProj(dim(2));

   if (xDel < -mlen.x()/2 || xDel > mlen.x()/2)
   {
      sign = xDel < 0 ? 1 : -1;
      c += dim(0)*sign;
   }
   if (yDel < -mlen.y()/2 || yDel > mlen.y()/2)
   {
      sign = yDel < 0 ? 1 : -1;
      c += dim(1)*sign;
   }
   if (zDel < -mlen.z()/2 || zDel > mlen.z()/2)
   {
      sign = zDel < 0 ? 1 : -1;
      c += dim(2)*sign;
   }

   return c;
}

Point Parameters::checkPeriodBound (Point & c)
{
   double sign;
   double xDel = c.scalarProj(dim(0));
   double yDel = c.scalarProj(dim(1));
   double zDel = c.scalarProj(dim(2));

   if (xDel < 0 || xDel >= mlen.x())
   {
      sign = xDel < 0 ? 1 : -1;
      c = (c+(dim(0)*sign));
   }
   if (yDel < 0 || yDel >= mlen.y())
   {
      sign = yDel < 0 ? 1 : -1;
      c = (c+(dim(1)*sign));
   }
   if (zDel < 0 || zDel >= mlen.z())
   {
      sign = zDel < 0 ? 1 : -1;
      c = (c+(dim(2)*sign));
   }
   return c;
}

Point Parameters::getRealAvg (Point & a, Point & b)
{
   Point avg = getRealDiff(a,b)/2 + b;
   return checkPeriodBound(avg);
}
   
/****** Depricated ******
void Parameters::delRandBond(vector<Atom> & atoms, int c)
{
   set<int>::iterator it = mdelCandidates.begin();
   it += rand() % nCandidates();
   int num = *it;
   if (atoms[num].getNumNeigh() > 2) 
   {
      if (atoms[num].delRandNeighbour())
      {
         
      }
      else 
      {
         mdelCandidates.erase(num);
         delRandBond(atoms, ++c);
      }
      genBondList(atoms);
      if (mangles.size() != 0) genAngleList(atoms);
   }
   else 
   {
   //if (c > 5) debug(atoms[num].getNumNeigh());
   if (c > 5) debug(num);
      mdelCandidates.erase(num);
   debug(c)
      delRandBond(atoms, ++c);
   }
}*/

Parameters::Parameters()
{
   mpnt = 0;
   mvar = 0;
   mdist = 1;
   mcxn = 0;
   mdim = Matrix3();
   mlen = Point();
   mbonds = set<Bond>();
   mangles = set<Angle>();
   moffCandidates = set<Bond>();
}

void Parameters::copy(const Parameters & other)
{
   pnt(other.mpnt);
   mcxn = other.mcxn;
   mdist = other.mdist;
   mbonds = other.mbonds;
   mangles = other.mangles;
   moffCandidates = other.moffCandidates;
   mdim = other.mdim;
   mlen = other.mlen;
}

void Parameters::writeBondLoc() const
{
   int i=0;
   set<Bond>::const_iterator it;
   fstream file ( "bondloc.out", fstream::in | fstream::out | fstream::trunc);
   
   for (it=mbonds.begin();it!=mbonds.end();++it)
   {
      file << ++i << ' ' << it->location().x() << ' ' << it->location().y() << ' ' << it->location().z() << endl;
   }
}

void Parameters::strain(const Point & strain)
{
   Matrix3 strainTensor;
   strainTensor.setDiag(strain);
   mdim *= strainTensor; 
   for (int i=0;i<3;i++)
      mlen.setCoord(dim(i).distance(), i);
}

void Parameters::strain(const Point & strain, const Point& shear)
{
   Matrix3 strainTensor;
   strainTensor.setDiag(strain);
   strainTensor.setSym(1,2,shear.x());
   strainTensor.setSym(0,2,shear.y());
   strainTensor.setSym(0,1,shear.z());
   mdim *= strainTensor; 
   for (int i=0;i<3;i++)
      mlen.setCoord(dim(i).distance(), i);
}

void Parameters::strain(const Matrix3& strainTensor)
{
   mdim *= strainTensor;
   for (int i=0;i<3;i++)
      mlen.setCoord(dim(i).distance(), i);
}

void Parameters::setDim (const Matrix3& other) 
{
   mdim = other;
   for (int i=0; i<3; i++)
      mlen.setCoord( other.getPoint(i).distance(), i );
}

void Parameters::printCellDim () const
{
   cout << "0 " << mdim.get(0,0);
   cout << " 0 " << mdim.get(1,1);
   cout << " 0 " << mdim.get(2,2);
}
