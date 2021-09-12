tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    combustion->correct();
    Qdot = combustion->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                parcels.SYi(i, Yi)
              + combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            YEqn.relax();

            fvOptions.constrain(YEqn);

            YEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
    
    Info<< "Yt gas min/max   " << min(Yt).value() << ", "
        << max(Yt).value() << endl;
}

const dictionary& subModelDict = parcels.subModelProperties();
vector position = subModelDict.subDict
    (
        "injectionModels"
    ).subDict("model1").lookup("position");

vector direction = subModelDict.subDict
    (
        "injectionModels"
    ).subDict("model1").lookup("direction");

 scalar vaporPenetration(0);
 //scalar mp(0);

//    const word liquids(thermo.lookup("liquids"));
  //  const label fuelIndex(composition.species()[liquids]);

// scalarField dist("dist",);
//
scalarField dist(mesh.nCells());

labelField index(5750);

std::ofstream file;
file.open ("PL.txt", std::ofstream::out | std::ofstream::app);

std::ofstream file2;
file2.open ("2D.txt", std::ofstream::out | std::ofstream::app);


forAll(Y[0],celli)
{
    if(Y[0][celli]  >= 0.001)
    {
        vector raw = position - mesh.C()[celli];
    
        if (mag(raw&direction) > vaporPenetration)
        {
           vaporPenetration = mag(raw&direction);
        }   
    }
} 


const volScalarField & al = parcels.alpha();
scalarField liquidmassf = al.internalField();


const volScalarField & ML = parcels.PARCELMASS();
scalarField liquidmass = ML.internalField();

int thisPntFieldSize = mesh.nCells();

labelList pntFieldSizePerProc( Pstream::nProcs(), 0 );
label thisProcNb = Pstream::myProcNo();
pntFieldSizePerProc[thisProcNb] = thisPntFieldSize;

reduce( pntFieldSizePerProc, sumOp<labelList>() );

int totalSize = sum( pntFieldSizePerProc );

pointField global_pntField( totalSize, vector::zero );
scalarField vmf( totalSize, 0.0);
scalarField lmf( totalSize, 0.0);

scalarField lm( totalSize, 0.0);

scalar lpl = parcels.penetration(0.92);
int pntsBefore;


if( thisPntFieldSize > 0 )
{

//int pntsBefore;

     if( thisProcNb ==0)
     {
         pntsBefore =0;
     }
     else
     {
        pntsBefore =0;
        for( int i =0; i< thisProcNb; i++)
        {
            pntsBefore  += pntFieldSizePerProc[i];
        }
     }

    for( int i = 0; i < thisPntFieldSize; i++)
    {
        global_pntField[ pntsBefore + i ] = mesh.C()[ i ];
        vmf[ pntsBefore + i ] = Y[0][i];
        lmf[ pntsBefore + i ] = liquidmassf[i];
        lm[ pntsBefore + i ] = liquidmass[i]; 
    }
}

reduce( global_pntField, sumOp<pointField>() );
reduce( vmf,sumOp<scalarField>() );
reduce( lmf,sumOp<scalarField>() );
reduce( lm,sumOp<scalarField>() );
reduce(vaporPenetration, maxOp<scalar>());
reduce(lpl, maxOp<scalar>());


if (Pstream::master() == true)
{
file <<  runTime.timeName() << "," << vaporPenetration  << "," << lpl<<  std::endl;

scalar xx = x/500;

label yy = round(xx);


if ( fabs(yy-xx) < 1e-10)
{
    forAll(global_pntField,celli)
    {
        if(global_pntField[celli][1]  <= 0.000001 && global_pntField[celli][1] >= -0.000001 &&  global_pntField[celli][0] <= 0.02  && global_pntField[celli][0] >= -0.02 &&  global_pntField[celli][2] <= 0.100)
        {
        file2 << runTime.timeName() << "," << global_pntField[celli][0] << "," << global_pntField[celli][1] << "," << global_pntField[celli][2] << "," << vmf[celli] << "," << lmf[celli] <<"," << lm[celli]<<  std::endl;

        }
    }

}
//Info << mesh.nCells()<< endl;

x = x+1;

}

file2.close();
file.close();










