#ifndef COSMOINTERFACE_EVOLVERS_KERNELS_HSCALARKERNELS_H
#define COSMOINTERFACE_EVOLVERS_KERNELS_HSCALARKERNELS_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "TempLat/util/tdd/tdd.h"
#include "CosmoInterface/definitions/potential.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"
#include "CosmoInterface/definitions/fieldfunctionals.h"

namespace TempLat {


    /** \brief A class that computes the kernel for the scalar singlets.
     *
     * 
     * Unit test: make test-scalarkernels
     **/


    class HScalarSingletKernels {
    public:
        /* Put public methods here. These should change very little over time. */
        HScalarSingletKernels() = delete;
        
        template <class Model>
        static auto get(Model& model){
            
            auto ScalarContribution = FieldFunctionals::grad2S(model, 0_c) + FieldFunctionals::pi2S(model, 0_c) - 4 * Potential::potential(model);
            auto GhostContribution = FieldFunctionals::grad2GS(model, 0_c) + FieldFunctionals::pi2GS(model, 0_c) - 4 * Potential::ghostmass(model);
            auto Interactionpotential = - 4 * Potential::ghostpotential(model);

            return (pow(model.aI, 1 + model.alpha) *
                    LatLapl<Model::NDim>(model.fldH(0_c))
                    - pow(model.aI, 3 + model.alpha) * pow<2>(model.fStar / Constants::reducedMPlanck<double>) * (ScalarContribution - GhostContribution + Interactionpotential));
        }

		
    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */



    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };



} /* TempLat */

#ifdef TEMPLATTEST
#include "CosmoInterface/evolvers/kernels/scalarsingletkernels_test.h"
#endif


#endif
