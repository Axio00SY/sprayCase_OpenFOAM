/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    sprayFoam

Description
    Transient solver for compressible, turbulent flow with a spray particle
    cloud.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "basicSprayCloud.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//int PRflag = 0;
//int Evapflag = 1;
int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    scalarField dist1(mesh.nCells());
    //scalar x =0;

    Info<< "\nStarting time loop\n" << endl;

    //output parcel size, definition, PingYI, 20200920
    const volScalarField & Dl = parcels.ParcelSize();
    scalarField Dsize = Dl.internalField();

    std::ofstream file;
    file.open ("D2.txt", std::ofstream::out | std::ofstream::app);
    //

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        if (!pimple.frozenFlow())
        {
            #include "rhoEqn.H"

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YEqn.H"
                #include "EEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();

            if (runTime.write())
            {
                combustion->Qdot()().write();
            }
        }
        else
        {
            if (runTime.writeTime())
            {
                parcels.write();
            }
        }

        //output parcel size, definition, PingYI, 20200920
        scalar Dlmax  = 1.0e+06*max(parcels.ParcelSize()).value();
        scalar Dlmin  = 1.0e+06*parcels.Dmin();
        scalar Dlmean = 1.0e+06*parcels.Dmean();

   	//Info<< "Dlmax  ====== " << Dlmax  <<endl;
    	//Info<< "Dlmin  ====== " << Dlmin  <<endl;
    	//Info<< "Dlmean ====== " << Dlmean <<endl;
        file<< runTime.timeName()<<", "<<Dlmax<<", "<<Dlmin<<", "<<Dlmean<<std::endl;
        //
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }
    file.close();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
