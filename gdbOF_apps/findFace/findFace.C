/*---------------------------------------------------------------------------*\
# ----------------------------- CIMEC-INTEC ----------------------------------- 
# ------- Centro Internacional de Métodos Computacionales en Ingeniería -------
# -------- (International Center for Numerical Methods in Engineering) --------  
# ----------------------- Santa Fe - Argentina --------------------------------

License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
	findFace

Description
   

Created By:
    J.M. Gimenez(1), S.Márquez Damián(2), N. Nigro(2)
		(1)FICH - UNL - Santa Fe, Argentina
		(2)CIMEC (INTEC-UNL-CONICET), Santa Fe, Argentina
	contact: jmarcelogimenez@gmail.com
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    //argList::validArgs.append("patch");
    argList::validArgs.append("x");
    argList::validArgs.append("y");
    argList::validArgs.append("z");

	Foam::argList args(argc, argv, false,false);

#   include "createTime.H"
#   include "createMesh.H"

	//char* patch = argv[1];
	double x = atof(argv[1]);
	double y = atof(argv[2]);
	double z = atof(argv[3]);

	//Info << endl << "x:" << x << " - y: " << y << " - z:" << z << endl;

	// http://foam.sourceforge.net/doc/Doxygen/html/meshSearch_8H_source.html

	point p(x, y, z);
	//label patchID = mesh.boundaryMesh().findPatchID(patch); 	//encuentra el ID del patch buscado
	meshSearch MS(mesh);
	label faceID = MS.findNearestBoundaryFace(p); 				// "face" mas cercana al punto
	label patchID = mesh.boundaryMesh().whichPatch(faceID);		// a que patch se corresponde esa "face"
	const polyPatch& cPatch = mesh.boundaryMesh()[patchID];

	Info << endl << "RESULTS: " << endl
		<< "Nearest patch ID: " << patchID << endl
		<< "Nearest face number: " << cPatch.whichFace(faceID) << endl;

    return 0;
}


// ************************************************************************* //
