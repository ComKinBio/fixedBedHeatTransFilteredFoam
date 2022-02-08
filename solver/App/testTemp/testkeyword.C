/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    manulDecomposeFile

Description
    This utility also read gridCoarseDict to generate decomposition file
    One processer size is required by the time of coarseGridsize

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "entry.H"
#include "labelList.H"
#include "vectorList.H"
#include "DynamicList.H"
#include "IOobjectList.H"
#include "labelIOField.H"
#include "labelFieldIOField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //read input settings
    IOdictionary couplingDict
    (
        IOobject
        (
            "couplingDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    const dictionary TFMDict = couplingDict.subDict("TFM");
    
    Info<< "const dictionary TFMDict.toc():"<<TFMDict.toc()<<"\n" << endl;
    
    Info<< "const dictionary TFMDict.toc()[0]:"<<TFMDict.toc()[0]<<"\n" << endl;
    
    forAllConstIter(IDLList<entry>, TFMDict, iter)
    {
        Info<< "iter().keyword():"<<iter().keyword()<<"\n" << endl;
        
        if (iter().isDict())
        {
            Info<< "iter().isDict():"<<"\n" << endl;
        }
    }
        
    
    

    
                    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
