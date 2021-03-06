#include "LDLTSolver.h"

#include <sofa/core/ObjectFactory.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
using std::cerr;
using std::endl;

namespace sofa {
namespace component {
namespace linearsolver {

SOFA_DECL_CLASS(LDLTSolver);
int LDLTSolverClass = core::RegisterObject("Direct LDLT solver").add< LDLTSolver >();

typedef AssembledSystem::vec vec;






LDLTSolver::LDLTSolver() 
    : KKTSolver()
    , regularize( initData(&regularize, std::numeric_limits<real>::epsilon(), "regularize", "add identity*regularize to matrix H to make it definite."))
    , pimpl()
{

}

LDLTSolver::~LDLTSolver() {

}


void LDLTSolver::factor(const AssembledSystem& sys) {

    typedef AssembledSystem::dmat dmat;

    if( !sys.isPIdentity ) // replace H with P^T.H.P to account for projective constraints
    {
        if( regularize.getValue() != (SReal)0.0 ) // add a tiny diagonal matrix to make H psd.
        {
            system_type::cmat identity(sys.m,sys.m);
            identity.setIdentity();
            pimpl->Hinv.compute( sys.P.transpose() * sys.H * sys.P + identity * regularize.getValue() );
        }
        else
            pimpl->Hinv.compute( sys.P.transpose() * sys.H * sys.P );
    }
    else
    {
        if( regularize.getValue() != (SReal)0.0 ) // add a tiny diagonal matrix to make H psd.
        {
            system_type::rmat identity(sys.m,sys.m);
            identity.setIdentity();
            pimpl->Hinv.compute( sys.H + identity * regularize.getValue() );
        }
        else
            pimpl->Hinv.compute( sys.H );
    }

    if( pimpl->Hinv.info() == Eigen::NumericalIssue ) {
        std::cerr << "LDLTSolver::factor: H is not psd. System solution will be wrong. P is identity=" << sys.isPIdentity << " regularize=" << regularize.getValue() << std::endl;

        std::cerr << pimpl->H << std::endl;
    }

    pimpl->dt = sys.dt;
    pimpl->m = sys.m;
    pimpl->n = sys.n;

    //    if( debug.getValue() ){
    //        cerr<< "LDLTSolver::factor, H = " << sys.H << endl;
    //    }

    if( sys.n ) {
        pimpl_type::cmat schur(sys.n, sys.n);
        pimpl_type::cmat PJT = sys.P.transpose() * sys.J.transpose(); //yes, we have to filter J, even if H is filtered already. Otherwise Hinv*JT has large (if not infinite) values on filtered DOFs

        pimpl->HinvPJT.resize(sys.m, sys.n);
        pimpl->HinvPJT = pimpl->Hinv.solve( PJT );

        schur = (sys.C.transpose() + (PJT.transpose() * pimpl->HinvPJT )).selfadjointView<Eigen::Upper>();
        if( debug.getValue() ){
            cerr<< "LDLTSolver::factor, PJT = " << endl << dmat(PJT) << endl;
            cerr<< "LDLTSolver::factor, HinvPJT = " << endl << dmat(pimpl->HinvPJT) << endl;
            cerr<< "LDLTSolver::factor, PJT.transpose() = " << endl << dmat(PJT.transpose()) << endl;
            cerr<< "LDLTSolver::factor, C = " << endl << dmat(sys.C) << endl;
            cerr<< "LDLTSolver::factor, schur = " << endl << dmat(schur) << endl;
        }

        pimpl->schur.compute( schur );

        if( pimpl->schur.info() == Eigen::NumericalIssue ) {
            std::cerr << "LDLTSolver::factor: schur is not psd. System solution will be wrong." << std::endl;
            std::cerr << schur << std::endl;
        }
    } else {
        // nothing lol
    }


}


void LDLTSolver::solve(AssembledSystem::vec& res,
                       const AssembledSystem& sys,
                       const AssembledSystem::vec& rhs) const {

    assert( res.size() == sys.size() );
    assert( rhs.size() == sys.size() );


    vec Pv = (sys.P * rhs.head(sys.m));

    typedef AssembledSystem::dmat dmat;

    if( debug.getValue() ){
        cerr<<"LDLTSolver::solve, rhs = " << rhs.transpose() << endl;
        cerr<<"LDLTSolver::solve, Pv = " << Pv.transpose() << endl;
        cerr<<"LDLTSolver::solve, H = " << endl << dmat(sys.H) << endl;
    }

    // in place solve
    Pv = pimpl->Hinv.solve( Pv );
    if( debug.getValue() ){
        cerr<<"LDLTSolver::solve, free motion solution = " << Pv.transpose() << endl;
        cerr<<"LDLTSolver::solve, verification = " << (sys.H * Pv).transpose() << endl;
        cerr<<"LDLTSolver::solve, sys.m = " << sys.m << ", sys.n = " << sys.n << ", rhs.size = " << rhs.size() << endl;

    }
    res.head( sys.m ) = sys.P * Pv;

    if( sys.n ) {
        vec tmp = rhs.tail( sys.n ) - pimpl->HinvPJT.transpose() * rhs.head( sys.m );


        // lambdas
        res.tail( sys.n ) = pimpl->schur.solve( tmp );

        // constraint forces
        res.head( sys.m ) += sys.P * (pimpl->HinvPJT * res.tail( sys.n));
        if( debug.getValue() ){
            cerr<<"LDLTSolver::solve, free motion constraint error= " << -tmp.transpose() << endl;
            cerr<<"LDLTSolver::solve, lambda = " << res.tail(sys.n).transpose() << endl;
            cerr<<"LDLTSolver::solve, constraint forces = " << (sys.P * (pimpl->HinvPJT * res.tail( sys.n))).transpose() << endl;
        }
    }

} 


}
}
}

