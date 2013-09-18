#ifndef FXN_H
#define FXN_H

#include <matrix.hpp>
#include <matrix3.hpp>
#include <vector>
#include <string>
#include <point.hpp>
#include <atom.hpp>
#include <cmath>
#include <bond.hpp>
#include <angle.hpp>
#include <set>

class Atom;
class Bond;
class Angle;

//All info needed about the function
class Parameters
{
   private:
      int mpnt; //number of points
      int mvar; //number of variables, must always be 3*pnt
      int mcxn; //number of maximum connections per atom
      Matrix3 mdim; //matrix of the cell vectors
      Point mlen; //The length of each cell vector
      double mdist; //Optimal bond distance
      std::set<Bond> mbonds; //List of bonds
      std::set<Angle> mangles; //List of bonds
      std::set<int> mdelCandidates; //List of atom indexes that can be deleted from
      void copy (const Parameters&);
   public:
      Parameters ();
      Parameters (const Parameters& other) {copy(other);}
      int pnt () const {return mpnt;}
      int var () const {return mvar;}
      int cxn () const {return mcxn;}
      double dist () const {return mdist;}
      Point len () const {return mlen;}
      const std::set<Bond>& bonds () const {return mbonds;}
      const std::set<Angle>& angles () const {return mangles;}
      Point dim (int n) const {return mdim.getPoint(n);}
      const Matrix3& dim () const {return mdim;}
      double volume() const {return fabs( mdim.tripleProduct() );}
      int nBonds () const {return mbonds.size();}
      int nAngles () const {return mangles.size();}
      int nCandidates () const {return mdelCandidates.size();}

      void pnt (int n) {mpnt = n; mvar = n*3;}
      void cxn (int n) {mcxn = n;}
      void dist (double n) {mdist = n;}
      void setDim (Point pt, int n) {mdim.setPoint(pt,n); mlen.setCoord(pt.distance(), n);}
      void setDim (const Matrix3&);
      void genBondList (std::vector<Atom>&);
      void genAngleList (std::vector<Atom>&);
      void delRandBond (std::vector<Atom>&, int);
      void genDelList (std::vector<Atom>&);

      Parameters& operator= (const Parameters& other) {copy(other); return *this;}
      void strain (const Point&);
      void strain (const Point&, const Point&);
      void strain (const Matrix3&);
      
      Point getRealDiff (Point&, Point&); //Returns the distance between atoms bonded across cell boundaries
      Point checkPeriodBound (Point&); //Shifts a point back into the cell if the optimization moves it outside
};

double optimizer (std::vector<Point>&, Parameters&, double, double);

#endif
