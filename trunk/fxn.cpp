#include <iostream>
#include <vector>
#include <string>
#include <point.hpp>
#include <fxn.hpp>
#include <atom.hpp>
#include <bond.hpp>
#include <set>
#include <ctime>
#include <cstdlib>
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
   
void Parameters::genBondList(vector<Atom> & atoms)
{
   unsigned int i;
   int j;

   for (i=0;i<atoms.size();i++)
   {
      for (j=0;j<atoms[i].getNumNeigh();j++)
      {
         mbonds.insert(Bond(&atoms[i], atoms[i].getNeighbour(j) ) );
      }
   }
}

void Parameters::genAngleList(vector<Atom> & atoms)
{
   unsigned int i;
   int j,k;

   for (i=0;i<atoms.size();i++)
   {
      for (j=0;j<atoms[i].getNumNeigh();j++)
      {
         for (k=0;k<atoms[i].getNeighbour(j)->getNumNeigh();k++)
         {
            try
            {
               mangles.insert( Angle( &atoms[i], atoms[i].getNeighbour(j), 
                                      atoms[i].getNeighbour(j)->getNeighbour(k) ) );
            }
            catch (BadStructureException e) { continue; }
         }
      }
   }
}

void Parameters::delRandBond(vector<Atom> & atoms)
{
   srand ( time(NULL) );
   int num = rand() % atoms.size();
   if (atoms[num].getNumNeigh() > 2) 
   {
      atoms[num].delRandNeighbour();
      genBondList(atoms);
      if (mangles.size() != 0) genAngleList(atoms);
   }
   else delRandBond(atoms);
}

Parameters::Parameters()
{
   mpnt = 0;
   mvar = 0;
   mdist = 1;
   mcxn = 0;
   mk = 1;
}

void Parameters::copy(const Parameters & other)
{
   pnt(other.mpnt);
   mcxn = other.mcxn;
   mdist = other.mdist;
   mbonds = other.mbonds;
   mdim = other.mdim;
   mlen = other.mlen;
   mk = other.mk;
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

double Parameters::averageCoord (std::vector<Atom>& atoms)
{
   double coord = 0;
   for (unsigned int i=0; i<atoms.size(); i++) coord += atoms[i].getNumNeigh();
   return coord/(double)atoms.size();
}
