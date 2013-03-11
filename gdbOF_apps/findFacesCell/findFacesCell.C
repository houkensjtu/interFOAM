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
   Find the faces ID's for the cellID specified

Created By:
    J.M. Gimenez(1), S.Márquez Damián(2), N. Nigro(2)
		(1)FICH - UNL - Santa Fe, Argentina
		(2)CIMEC (INTEC-UNL-CONICET), Santa Fe, Argentina
	contact: jmarcelogimenez@gmail.com

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellFeatures.H"
#include "meshTools.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::validArgs.append("cellID");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

	label cellID = atoi(argv[1]);
	cellFeatures CF(mesh,1,cellID);
	const faceList& faces = CF.faces();
	
	Info << endl << "RESULTS: " << endl
		 << "Face List: " << faces << endl;


	const List<int> & owner = mesh.faceOwner();
	const List<int> & neighbour = mesh.faceNeighbour();
	const faceList & ff = mesh.faces(); 
	forAll (ff, faceID){
		if (mesh.isInternalFace(faceID)){
			if((owner[faceID] == cellID) || (neighbour[faceID] == cellID))
				//Info<<"es cara Interna: "<<faceID<<endl;
				Info<<0<<" "<<faceID<<endl;
		}else{
			if (owner[faceID] == cellID){
				label patchID = mesh.boundaryMesh().whichPatch(faceID);		// a que patch se corresponde esa "face"
				const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
				//Info<<"es cara Externa: "<<faceID<<", del patchID: "<<patchID<<", cara local: "<<cPatch.whichFace(faceID)<<endl;
				Info<<1<<" "<<patchID<<" "<<cPatch.whichFace(faceID)<<endl;
			}
		}
	}
	/*
	La caras internas las tomo tal cual (campo[faceID])
	las caras externas, capturo el patchID y el faceID local.
	*/
	
    return 0;
}


// ************************************************************************* //
