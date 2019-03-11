# Implementation of the Euler method in OpenFOAM

The forward Euler integration method is the simplest `ddt` scheme in OpenFOAM. We'll learn how it is implemented and then try to adapt it to implement a PiT scheme.

In this analysis I will break with some of the coding standards normally employed in OpenFOAM sources, mainly because they break with the literate programming exposition. Here I deleted some of the very obvious comments. Comments should make a code *more* readable.

The OpenFOAM implementation consists of three files:

* `EulerDdtScheme.H` templates
* `EulerDdtScheme.C` implementation, template specialisations
* `EulerDdtSchemes.C` global announcement to all code linking, that this integrator exists.

The last file is also the shortest

``` {.cpp file=ddtEuler/EulerDdtSchemes.C}
#include "EulerDdtScheme.H"
#include "fvMesh.H"

makeFvDdtScheme(JHEulerDdtScheme)
```

Here, `makeFvDdtScheme` is a macro that adds the necessary boiler plate.

## Header file

The header file starts with usual top level comments and copyright notices, then include guards.

``` {.cpp file=ddtEuler/EulerDdtScheme.H}
#ifndef EulerDdtScheme_H
#define EulerDdtScheme_H

#include "ddtScheme.H"

<<euler-namespace>>
<<euler-include-c-file>>

#endif
```

### Namespaces

The namespace is `Foam::fv`.

``` {.cpp #euler-namespace}
namespace Foam {
    namespace fv {
        <<euler-ddt-class>>
    }
}
```

### Include C file

Optionally, I'm guessing for debug purposes, the actual 'C' file is included in the header. In general this is a bad idea.

``` {.cpp #euler-include-c-file}
#ifdef NoRepository
    #include "EulerDdtScheme.C"
#endif
```

### Class definition


``` {.cpp #euler-ddt-class}
template<class Type>
class JHEulerDdtScheme
:
    public ddtScheme<Type>
{
    <<euler-private-methods>>

public:
    <<euler-run-time-type-information>>

    <<euler-constructors>>
    <<euler-methods>>
};

<<euler-scalar-specializations>>
```

### Private methods

Copying and assignment are made private.

``` {.cpp #euler-private-methods-old}
//- Disallow default bitwise copy construct
EulerDdtScheme(const EulerDdtScheme&);

//- Disallow default bitwise assignment
void operator=(const EulerDdtScheme&);
```

In modern C++ this code would look like this

``` {.cpp #euler-private-methods}
JHEulerDdtScheme(const JHEulerDdtScheme&) = delete;
void operator=(const JHEulerDdtScheme&) = delete;
```

Comments are then completely superfluous.

### Run-time type information

This adds run-time information on the object.

``` {.cpp #euler-run-time-type-information}
TypeName("JH-Euler");
```

This is a macro that expands to create a `virtual` function `type()` returning a `word` (class ultimately deriving from `std::string`). This virtual method is rooted at the `ddtScheme` template class.

### Constructors

``` {.cpp #euler-constructors}
//- Construct from mesh
JHEulerDdtScheme(const fvMesh& mesh)
:
    ddtScheme<Type>(mesh)
{}

//- Construct from mesh and Istream
JHEulerDdtScheme(const fvMesh& mesh, Istream& is)
:
    ddtScheme<Type>(mesh, is)
{}
```

### Methods

``` {.cpp #euler-methods}
const fvMesh& mesh() const
{
    return fv::ddtScheme<Type>::mesh();
}

tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
(
    const dimensionedScalar&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
(
    const volScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& psi
);

tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fvcDdt
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>&
);

tmp<fvMatrix<Type>> fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>&
);

tmp<fvMatrix<Type>> fvmDdt
(
    const dimensionedScalar&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

tmp<fvMatrix<Type>> fvmDdt
(
    const volScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>&
);

tmp<fvMatrix<Type>> fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& psi
);

typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

tmp<fluxFieldType> fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
);

tmp<fluxFieldType> fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
);

tmp<fluxFieldType> fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
);

tmp<fluxFieldType> fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
);

tmp<surfaceScalarField> meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
);
```

### Scalar specializations

``` {.cpp #euler-scalar-specializations}
template<>
tmp<surfaceScalarField> JHEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf
);

template<>
tmp<surfaceScalarField> JHEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi
);

template<>
tmp<surfaceScalarField> JHEulerDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf
);

template<>
tmp<surfaceScalarField> JHEulerDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi
);
```

## Implementation

``` {.cpp file=ddtEuler/EulerDdtScheme.C}
#include "EulerDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

namespace Foam
{
    namespace fv
    {
        <<euler-fvc>>
        <<euler-fvm>>
        <<euler-corr>>
    }
}
```

### `fvc` **c**urrent derivative

Returns current values for derivative. This stands for *finite volume calculus*. There are several overloaded implementations. 

#### `dimensioned`

I'm guessing the `dimensioned` type encapsulates a scalar value.

``` {.cpp #euler-methods}
tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
(
    const dimensioned<Type>&
);
```

``` {.cpp #euler-fvc}
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
JHEulerDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    <<euler-compute-rDeltaT>>
    <<euler-ddtIO>>

    if (mesh().moving())
    {
        <<euler-moving-mesh>>
    }
    else
    {
        <<euler-static-mesh>>
    }
}
```

* Compute $1 / \Delta T$.
  
``` {.cpp #euler-compute-rDeltaT}
dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();
```

* Initialize IO
  
``` {.cpp #euler-ddtIO}
IOobject ddtIOobject
(
    "ddt("+dt.name()+')',
    mesh().time().timeName(),
    mesh()
);
```

* For moving mesh, create new `GeometricField` `tdtdt`, set its value.

``` {.cpp #euler-moving-mesh}
tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
(
    new GeometricField<Type, fvPatchField, volMesh>
    (
        ddtIOobject,
        mesh(),
        dimensioned<Type>
        (
            "0",
            dt.dimensions()/dimTime,
            Zero
        )
    )
);

tdtdt.ref().primitiveFieldRef() =
    rDeltaT.value()*dt.value()*(1.0 - mesh().Vsc0()/mesh().Vsc());

return tdtdt;
```

* For static mesh, return a new `GeometricField` with value `Zero`.

``` {.cpp #euler-static-mesh}
return tmp<GeometricField<Type, fvPatchField, volMesh>>
(
    new GeometricField<Type, fvPatchField, volMesh>
    (
        ddtIOobject,
        mesh(),
        dimensioned<Type>
        (
            "0",
            dt.dimensions()/dimTime,
            Zero
        ),
        calculatedFvPatchField<Type>::typeName
    )
);
```

#### `GeometricField` for volume mesh

``` {.cpp #euler-methods}
tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>&
);
```

``` {.cpp #euler-fvc}
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
JHEulerDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    vf()
                  - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
JHEulerDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    vf()
                  - vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.value()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
JHEulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    rho()*vf()
                  - rho.oldTime()()
                   *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.value()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
JHEulerDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    alpha()
                   *rho()
                   *vf()

                  - alpha.oldTime()()
                   *rho.oldTime()()
                   *vf.oldTime()()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.value()*
                (
                    alpha.boundaryField()
                   *rho.boundaryField()
                   *vf.boundaryField()

                  - alpha.oldTime().boundaryField()
                   *rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT
               *(
                   alpha*rho*vf
                 - alpha.oldTime()*rho.oldTime()*vf.oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
JHEulerDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+sf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            ddtIOobject,
            rDeltaT*(sf - sf.oldTime())
        )
    );
}
```

### `fvm` **m**atrix

Returns linear system to solve.

``` {.cpp #euler-fvm}
template<class Type>
tmp<fvMatrix<Type>>
JHEulerDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
JHEulerDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*rho.value()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
JHEulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
JHEulerDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() =
        rDeltaT*alpha.primitiveField()*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}
```

### Corrections

``` {.cpp #euler-corr}
template<class Type>
tmp<typename JHEulerDdtScheme<Type>::fluxFieldType>
JHEulerDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
    fluxFieldType phiCorr
    (
        phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            this->fvcDdtPhiCoeff(U.oldTime(), phiUf0, phiCorr)
           *rDeltaT*phiCorr
        )
    );
}
```

### `fvc` `PhiCorr`

``` {.cpp #euler-corr}
template<class Type>
tmp<typename JHEulerDdtScheme<Type>::fluxFieldType>
JHEulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    fluxFieldType phiCorr
    (
        phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
           *rDeltaT*phiCorr
        )
    );
}


template<class Type>
tmp<typename JHEulerDdtScheme<Type>::fluxFieldType>
JHEulerDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr(phiUf0 - fvc::dotInterpolate(mesh().Sf(), rhoU0));

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
                    mesh().time().timeName(),
                    mesh()
                ),
                this->fvcDdtPhiCoeff(rhoU0, phiUf0, phiCorr, rho.oldTime())
               *rDeltaT*phiCorr
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && Uf.dimensions() == rho.dimensions()*dimVelocity
    )
    {
        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr
        (
            phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
                    mesh().time().timeName(),
                    mesh()
                ),
                this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phiUf0,
                    phiCorr,
                    rho.oldTime()
                )*rDeltaT*phiCorr
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename JHEulerDdtScheme<Type>::fluxFieldType>
JHEulerDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiCorr
        (
            phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), rhoU0)
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
                    mesh().time().timeName(),
                    mesh()
                ),
                this->fvcDdtPhiCoeff
                (
                    rhoU0,
                    phi.oldTime(),
                    phiCorr,
                    rho.oldTime()
                )*rDeltaT*phiCorr
            )
        );
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        fluxFieldType phiCorr
        (
            phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
                    mesh().time().timeName(),
                    mesh()
                ),
                this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phi.oldTime(),
                    phiCorr,
                    rho.oldTime()
                )*rDeltaT*phiCorr
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<surfaceScalarField> JHEulerDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return mesh().phi();
}
```
