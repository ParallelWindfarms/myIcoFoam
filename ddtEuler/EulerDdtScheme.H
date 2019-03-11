// ------ language="C++" file="ddtEuler/EulerDdtScheme.H"
#ifndef EulerDdtScheme_H
#define EulerDdtScheme_H

#include "ddtScheme.H"

// ------ begin <<euler-namespace>>[0]
namespace Foam {
    namespace fv {
        // ------ begin <<euler-ddt-class>>[0]
        template<class Type>
        class JHEulerDdtScheme
        :
            public ddtScheme<Type>
        {
            // ------ begin <<euler-private-methods>>[0]
            JHEulerDdtScheme(const JHEulerDdtScheme&) = delete;
            void operator=(const JHEulerDdtScheme&) = delete;
            // ------ end
        
        public:
            // ------ begin <<euler-run-time-type-information>>[0]
            TypeName("JH-Euler");
            // ------ end
        
            // ------ begin <<euler-constructors>>[0]
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
            // ------ end
            // ------ begin <<euler-methods>>[0]
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
            // ------ end
            // ------ begin <<euler-methods>>[1]
            tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
            (
                const dimensioned<Type>&
            );
            // ------ end
            // ------ begin <<euler-methods>>[2]
            tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
            (
                const GeometricField<Type, fvPatchField, volMesh>&
            );
            // ------ end
        };
        
        // ------ begin <<euler-scalar-specializations>>[0]
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
        // ------ end
        // ------ end
    }
}
// ------ end
// ------ begin <<euler-include-c-file>>[0]
#ifdef NoRepository
    #include "EulerDdtScheme.C"
#endif
// ------ end

#endif
// ------ end
