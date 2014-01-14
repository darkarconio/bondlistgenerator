#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <atom.hpp>
#include <point.hpp>
#include <fxn.hpp>
#include <bond.hpp>

#define debug(x) cout << __LINE__ << ' ' << x << endl; cout.flush();

using namespace std;

vector<Atom> Atom::atomList;
Matrix Atom::adjMatrix = Matrix();
Parameters Atom::cellInfo = Parameters();

int main (int argc, char* argv[])
{
   string inputFile = argv[1];
   if (inputFile.compare("data.min") == 0)
   {
      bool bondLoc = !(bool)(strcmp(argv[2], "bond"));
      Atom::readMinAtoms(inputFile);
      Atom::genAtomDelList();
      if (bondLoc)
         Atom::cellInfo.writeBondLoc();
      else
         Atom::cellInfo.writeAtomLoc();
      Atom::cellInfo.printCellDim();
      return 0;
   }

   int exBond = atoi(argv[2]);
   int cellFactor = atoi(argv[3]);
   double percentDel = abs(atof(argv[4]));
   string atomOrBond = argv[5];
   unsigned int modDelete = abs(atoi(argv[6]));
   double delDist = abs(atof(argv[7]));
   srand ( time(NULL) );
   time_t start = time(NULL);

   Atom::readAtoms(inputFile);
   Atom::multiplyCell(Point(cellFactor));
   
   cout << "Connecting atoms..." << endl;
   Atom::connectAtoms(exBond);
   time_t connected = time(NULL);
   cout << "Connecting time: " << difftime(time(NULL), start) << endl;
  
   if (!(bool)atomOrBond.compare("atom"))
   {
      Atom::genAtomDelList();
      Atom::delPercentAtom(percentDel);
   }
   else if (!(bool)atomOrBond.compare("both"))
   {
      double percentDel2 = abs(atof(argv[8]));
      
      Atom::genAtomDelList();
      Atom::delPercentAtom(percentDel);
      
      Atom::genBondDelList(modDelete, delDist);
      Atom::delPercentBond(percentDel2, modDelete, delDist);
   }
   else
   {
      Atom::genBondDelList(modDelete, delDist);
      Atom::delPercentBond(percentDel, modDelete, delDist);
   }
   
   Atom::genBondList();
   Atom::genBondDelList(modDelete, delDist);
   Atom::genAtomDelList();
   
   time_t deletions = time(NULL);
   cout << "Deletion time: " << difftime(time(NULL), connected) << endl;
   
   cout << "Generating angle list..." << endl;
   Atom::genAngleList();
   cout << "Angle list time: " << difftime(time(NULL), deletions) << endl;
   
   Atom::outputAtoms(inputFile);
   Atom::cellInfo.writeBondLoc();
   
   cout << "Average Coordination: " << 2*(double)Atom::cellInfo.nBonds()/(double)Atom::cellInfo.nAtoms() << endl;
   cout << "Density: " << (double)Atom::cellInfo.nAtoms()/Atom::cellInfo.volume() << endl;
   cout << "Bond Density: " << 2*(double)Atom::cellInfo.nBonds()/Atom::cellInfo.volume() << endl;
   cout << "Num Coord 2: " << Atom::getNumCoordX(2) << endl;
   cout << "Avg Coord 2 Neighbour Coord: " << Atom::getCoordXNeighBond(2) << endl;
   cout << "Num Coord 3: " << Atom::getNumCoordX(3) << endl;
   cout << "Avg Coord 3 Neighbour Coord: " << Atom::getCoordXNeighBond(3) << endl;
   cout << "Num Coord 4: " << Atom::getNumCoordX(4) << endl;
   cout << "Avg Coord 4 Neighbour Coord: " << Atom::getCoordXNeighBond(4) << endl;
   cout << "Num Bonds: " << Atom::cellInfo.nBonds() << endl;
   cout << "Deletable Bonds: " << Atom::cellInfo.nBondCandidates() << endl;
   cout << "Num Atoms: " << Atom::getNumAtoms() << endl;
   cout << "Deletable Atoms: " << Atom::cellInfo.nAtomCandidates() << endl;
   cout << "Total time: " << difftime(time(NULL), start) << endl;

   return 0;
}
