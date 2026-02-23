#ifndef COSMOINTERFACE_EVOLVERS_KERNELS_GHOSTSCALARKERNELS_H
#define COSMOINTERFACE_EVOLVERS_KERNELS_GHOSTSCALARKERNELS_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "TempLat/util/tdd/tdd.h"
#include "CosmoInterface/definitions/potential.h"
#include "TempLat/lattice/algebra/spatialderivatives/latticelaplacian.h"

namespace TempLat {


    /** \brief A class that computes the kernel for the scalar singlets.
     *
     * 
     * Unit test: make test-scalarkernels
     **/


    class GhostScalarSingletKernels {
    public:
        /* Put public methods here. These should change very little over time. */
        GhostScalarSingletKernels() = delete;

        template <class Model, int N>
        static auto getON(Model& model, Tag<N> n){
        
        	// Returns kernel for scalar singlets (formed by laplacian and potential derivative):
            return (pow(model.aI, 1 + model.alpha) *
                    LatLapl<Model::NDim>(model.fldGS(n))
                    -  ((2.0)/(1.0 + 2.0 * model.fldH(0_c))) * (model.piGS(n) * model.piH(0_c) - GaugeDerivatives::GradientGhostScalarGradientH(model,n))
                    - pow(model.aI, 3 + model.alpha) * ((1.0 + 4.0 * model.fldH(0_c))/(1.0 + 2.0 * model.fldH(0_c))) * (- Potential::derivGS(model,n) + Potential::derivMassGS(model,n)));
        }

        template <class Model, int N>
        static auto get(Model& model, Tag<N> n){
        
        	// Returns kernel for scalar singlets (formed by laplacian and potential derivative):
            return (pow(model.aI, 1 + model.alpha) *
                    LatLapl<Model::NDim>(model.fldGS(n))
                    - pow(model.aI, 3 + model.alpha) * (- Potential::derivGS(model,n) + Potential::derivMassGS(model,n)));
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
