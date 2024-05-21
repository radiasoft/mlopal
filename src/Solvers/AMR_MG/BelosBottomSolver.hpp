//
// Class BelosBottomSolver
//   Interface to Belos solvers of the Trilinos package.
//
// Copyright (c) 2017 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Ippl.h"
#include "Utilities/OpalException.h"

extern Inform* gmsg;

template <class Level>
BelosBottomSolver<Level>::BelosBottomSolver(std::string solvertype,
                                            const std::shared_ptr<prec_t>& prec_p)
    : BottomSolver<Teuchos::RCP<amr::matrix_t>,
                                Teuchos::RCP<amr::multivector_t>,
                                Level>()
    , problem_mp( Teuchos::rcp( new problem_t() ) )
    , prec_mp(prec_p)
    , A_mp(Teuchos::null)
    , reltol_m(1.0e-9)
    , maxiter_m(100)
{
    this->initSolver_m(solvertype);
}


template <class Level>
void BelosBottomSolver<Level>::solve(const Teuchos::RCP<mv_t>& x,
                                     const Teuchos::RCP<mv_t>& b)
{
    /*
     * solve linear system Ax = b
     */
    
    // change sign of rhs due to change of A in BelosBottomSolver::setOperator
    b->scale(-1.0);
    
    problem_mp->setProblem(x, b);
    
    solver_mp->setProblem(problem_mp);
    
    Belos::ReturnType ret = solver_mp->solve();
    
    if ( ret != Belos::Converged ) {
        *gmsg << "Warning: Bottom solver not converged. Achieved tolerance"
              << " after " << solver_mp->getNumIters() << " iterations is "
              << solver_mp->achievedTol() << "." << endl;
    }
    
    // undo sign change
    b->scale(-1.0);
}


template <class Level>
void BelosBottomSolver<Level>::setOperator(const Teuchos::RCP<matrix_t>& A,
                                           Level* level_p)
{
    // make positive definite --> rhs has to change sign as well
    A_mp = Teuchos::rcp(new matrix_t(*A, Teuchos::Copy));
    A_mp->resumeFill();
    A_mp->scale(-1.0);
    A_mp->fillComplete();
    
    if ( problem_mp == Teuchos::null )
        throw OpalException("BelosBottomSolver::setOperator()",
                            "No problem defined.");
    
    problem_mp->setOperator(A_mp);
    
    static IpplTimings::TimerRef precTimer = IpplTimings::getTimer("AMR MG prec setup");

    if ( prec_mp != nullptr ) {
        IpplTimings::startTimer(precTimer);
        prec_mp->create(A_mp, level_p);
        IpplTimings::stopTimer(precTimer);
        problem_mp->setLeftPrec(prec_mp->get());
    }

    this->isInitialized_m = true;
}


template <class Level>
std::size_t BelosBottomSolver<Level>::getNumIters() {
    if ( solver_mp == Teuchos::null )
        throw OpalException("BelosBottomSolver::getNumIters()",
                            "No solver initialized.");
    
    return solver_mp->getNumIters();
}


template <class Level>
void BelosBottomSolver<Level>::initSolver_m(std::string solvertype) {
    
    Belos::SolverFactory<scalar_t, mv_t, op_t> factory;
    
    params_mp = Teuchos::parameterList();
    
    params_mp->set("Block Size", 1);
    params_mp->set("Convergence Tolerance", reltol_m);
    params_mp->set("Adaptive Block Size", false);
    params_mp->set("Use Single Reduction", true);
    params_mp->set("Explicit Residual Scaling", "Norm of RHS");
    params_mp->set("Maximum Iterations", maxiter_m);
    params_mp->set("Verbosity", Belos::Errors + Belos::Warnings);
    params_mp->set("Output Frequency", 1);
    
    solver_mp = factory.create(solvertype, params_mp);
}
