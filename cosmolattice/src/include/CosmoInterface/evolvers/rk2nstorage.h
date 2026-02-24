#ifndef COSMOINTERFACE_EVOLVERS_RK2NSTORAGE_H
#define COSMOINTERFACE_EVOLVERS_RK2NSTORAGE_H
/* File created by Adrien Florio, 2021 and modified for axion-U(1) by Ander Urio, 2025 */
/* Released under the MIT license, see LICENSE.md. */

#include "TempLat/util/tdd/tdd.h"
#include "TempLat/lattice/field/collections/vectorfieldcollection.h"
#include "TempLat/util/staticif.h"
#include "CosmoInterface/runparameters.h"
#include "CosmoInterface/evolvers/rk2nstorageparameters.h"
#include "CosmoInterface/evolvers/kernels/kernels.h"
#include "CosmoInterface/definitions/averages.h"



namespace TempLat {


    /** \brief A class which implements low storage ("2N-storage") explicit RK methods. This version is checked only for axion-U(1) model.
     *
     *
     * Unit test: make test-rk2nstorage
     **/


    template<typename Model>
    class RK2NStorage {
    public:

        using T = typename Model::FloatType;

        /* Put public methods here. These should change very little over time. */
        RK2NStorage(Model& model, RunParameters<T> runParams) :
        As(RK2NStorageParameters<T>::getAs(runParams.eType)),
        Bs(RK2NStorageParameters<T>::getBs(runParams.eType)),
        DeltaFldS("DeltaFldS",model.getToolBox(), runParams.getLatParams()),
        DeltaPiS("DeltaPiS",model.getToolBox(), runParams.getLatParams()),
        DeltaFldGS("DeltaFldGS",model.getToolBox(), runParams.getLatParams()),
        DeltaPiGS("DeltaPiGS",model.getToolBox(), runParams.getLatParams()),
        DeltaFldH("DeltaFldH",model.getToolBox(), runParams.getLatParams()),
        DeltaPiH("DeltaPiH",model.getToolBox(), runParams.getLatParams()),
        //oldFldS( RK2NStorageParameters<T>::isAdaptative(runParams.eType) ? std::make_unique<FieldCollection<Field, T, Model::Ns,true>>("oldFldS",model.getToolBox(), runParams.getLatParams()) : nullptr),
        //oldPiS( RK2NStorageParameters<T>::isAdaptative(runParams.eType) ? std::make_unique<FieldCollection<Field, T, Model::Ns,true>>("oldPiS",model.getToolBox(), runParams.getLatParams()) : nullptr),
        type(runParams.eType),
        Grav_ON(model.grav_ON),
        expansion(runParams.expansion),
        tolerance(runParams.tolerance),
        nextDt(model.dt)
        {

        }

        void evolve(Model &model, T tMinust0) {
            /*
             * 2N storage RK ....
             *
             *
             * */


            bool notPreciseEnough = false;
            if(RK2NStorageParameters<T>::isAdaptative(type)) {
                //(*oldFldS) = model.fldS;
                //(*oldPiS) = model.piS;
                //(*oldFldGWs) = (*model.fldGWs);
                //(*oldPiGWs) = (*model.piGWs);
                //oldA = model.aI;
                //oldADot = model.aDotI;
            }
            do {
                if(RK2NStorageParameters<T>::isAdaptative(type)) model.dt = nextDt;


                for (size_t i = 0; i < As.size(); ++i) {  // loop over operations...

                    if (Model::Ns > 0) deltaScalar(model, i, tMinust0);
                 
                    if (Model::NGs > 0) deltaGhostScalar(model, i, tMinust0);

                    if (Grav_ON == true) deltaHScalar(model, i, tMinust0);
                
                    if(expansion) deltaScaleFactor(model, i, tMinust0);

                    if (Model::Ns > 0) advanceScalar(model, i);
              
                    if (Model::NGs > 0) advanceGhostScalar(model, i);

                    if (Grav_ON == true) advanceHScalar(model, i);
               
                    if (expansion) advanceScaleFactor(model, i);
					
                    if(expansion){
                      Averages::setAllAverages(model);
                    }

                }


                if( RK2NStorageParameters<T>::isAdaptative(type)) {

                    T delta = sqrt(average(total(pow<2>(Bs.back() * DeltaFldS)))) / tolerance;

                     nextDt = model.dt *  0.95 * pow(1.0 / delta, 0.333333333);

                    if (delta > 1) {
                        //notPreciseEnough = true;
                        //model.fldS = (*oldFldS);
                        //model.piS = (*oldPiS);
                        //model.aDotI = oldADot;
                        //model.aI = oldA;
                        //Averages::setAllAverages(model);
                    } else notPreciseEnough = false;
                }
            } while(notPreciseEnough);
        }


        //The auxiliary vector fields used are the same as the ones in the gauge/electric field initializations

        void deltaScalar(Model& model, size_t i, T tMinust0){
            ForLoop(n, 0, Model::Ns - 1,
                    if (i == 0) {
                        if(Grav_ON == true){                
                            DeltaPiS(n) = model.dt * ScalarSingletKernels::getON(model, n); 
                        }
                        else{
                            DeltaPiS(n) = model.dt * ScalarSingletKernels::get(model, n);

                        }
                        DeltaFldS(n) = model.dt * model.piS(n);
                    
                    } else {
                        if(Grav_ON){                       
                            DeltaPiS(n) = As[i] * DeltaPiS(n) + model.dt * ScalarSingletKernels::getON(model, n);                          
                        }
                        else{
                            DeltaPiS(n) = As[i] * DeltaPiS(n) + model.dt * ScalarSingletKernels::get(model, n);

                        }
                        DeltaFldS(n) = As[i] * DeltaFldS(n) +  model.dt * model.piS(n);
                    }
            );
        }
             

        void deltaGhostScalar(Model& model, size_t i, T tMinust0){
            ForLoop(n, 0, Model::NGs - 1,
                    if (i == 0) {
                        if(Grav_ON == true){                
                            DeltaPiGS(n) = model.dt * GhostScalarSingletKernels::getON(model, n); 
                        }
                        else{
                            DeltaPiGS(n) = model.dt * GhostScalarSingletKernels::get(model, n);

                        }
                        DeltaFldGS(n) = model.dt * model.piGS(n);
                    
                    } else {
                        if (Grav_ON){                        
                            DeltaPiGS(n) = As[i] * DeltaPiGS(n) + model.dt * GhostScalarSingletKernels::getON(model, n);                          
                        }
                        else{
                            DeltaPiGS(n) = As[i] * DeltaPiGS(n) + model.dt * GhostScalarSingletKernels::get(model, n);

                        }
                        DeltaFldGS(n) = As[i] * DeltaFldGS(n) +  model.dt * model.piGS(n);
                    }
            );
        }

        void deltaHScalar(Model& model, size_t i, T tMinust0){
            
                    if (i == 0) {
                                        
                        DeltaPiH(0_c) = model.dt * HScalarSingletKernels::get(model); 
    
                        DeltaFldH(0_c) = model.dt * model.piH(0_c);
                    
                    } else {
                                                
                        DeltaPiH(0_c) = As[i] * DeltaPiH(0_c) + model.dt * HScalarSingletKernels::get(model);                          
                        
                        DeltaFldH(0_c) = As[i] * DeltaFldH(0_c) +  model.dt * model.piH(0_c);
                    }
            
        }

        void deltaScaleFactor(Model& model, size_t i, T tMinust0){
            if (i == 0) {

               deltaADot = model.dt * ScaleFactorKernels::get(model);

               deltaA = model.dt * model.aDotI;
            } else {
                
                deltaADot = As[i] * deltaADot + model.dt * ScaleFactorKernels::get(model);
                
                deltaA = As[i] * deltaA + model.dt * model.aDotI;
            }
            
        }


        void advanceScalar(Model& model, size_t i){
            ForLoop(n, 0, Model::Ns - 1,
                    model.piS(n) += Bs[i] * DeltaPiS(n);
                    model.fldS(n) += Bs[i] * DeltaFldS(n);
            );

        }

        void advanceGhostScalar(Model& model, size_t i){
            ForLoop(n, 0, Model::NGs - 1,
                    model.piGS(n) += Bs[i] * DeltaPiGS(n);
                    model.fldGS(n) += Bs[i] * DeltaFldGS(n);
            );

        }

        void advanceHScalar(Model& model, size_t i){
           
                    model.piH(0_c) += Bs[i] * DeltaPiH(0_c);
                    model.fldH(0_c) += Bs[i] * DeltaFldH(0_c);

        }

        void advanceScaleFactor(Model& model, size_t i){
            model.aDotI += Bs[i] * deltaADot;
            model.aI += Bs[i] * deltaA;
        }



        // This function is called before doing the measurements. It is used only to set aDotI to its correct value in case
        // of a background expansion.
        void sync(Model& model, T tMinust0) {
          //  if(fixedBackground) model.aDotI = aBackground.dot(tMinust0);
        }

        void storeMomentaAverages(Model& model) {
            if (Model::Ns > 0)  model.pi2AvI = Averages::pi2S(model);
            if (Model::NCs > 0) model.CSpi2AvI = Averages::pi2CS(model);
            if (Model::NSU2Doublet > 0) model.SU2DblPi2AvI = Averages::pi2SU2Doublet(model);
            if (Model::NU1 > 0) model.U1pi2AvI = Averages::pi2U1(model);
            if (Model::NSU2 > 0) model.SU2pi2AvI = Averages::pi2SU2(model);

        }

        void storeFieldsAverages(Model& model) {
            if (Model::Ns > 0) model.grad2AvI = Averages::grad2S(model);
            if (model.NCs > 0) model.CSgrad2AvI = Averages::grad2CS(model);
            if (model.NSU2Doublet > 0) model.SU2DblGrad2AvI = Averages::grad2SU2Doublet(model);
            if (Model::NU1 > 0)  model.U1Mag2AvI = Averages::B2U1(model);
            if (Model::NSU2 > 0)  model.SU2Mag2AvI = Averages::B2SU2(model);
            model.potAvI = average(Potential::potential(model));
        }

    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */

        const std::vector<T> As;
        const std::vector<T> Bs;

        //auxiliary scalar fields for RK
        // --> Scalar singlets
        FieldCollection<Field, T, Model::Ns,true> DeltaFldS; 
        FieldCollection<Field, T, Model::Ns,true> DeltaPiS;

        FieldCollection<Field, T, Model::NGs,true> DeltaFldGS; 
        FieldCollection<Field, T, Model::NGs,true> DeltaPiGS;
        
        FieldCollection<Field, T, 1,true> DeltaFldH; 
        FieldCollection<Field, T, 1,true> DeltaPiH;

        //std::unique_ptr<FieldCollection<Field, T, Model::Ns,true>> oldFldS;
        //std::unique_ptr<FieldCollection<Field, T, Model::Ns,true>> oldPiS;
        //std::unique_ptr<FieldCollection<Field, T, Model::Ns,true>> errorFld;


        T oldA, oldADot, deltaA, deltaADot;//

        EvolverType type;
        bool Grav_ON;
        bool expansion;

        T lowOrderControl;
        T largeOrderControl;
        T tolerance;
        T nextDt;


    };

    struct RK2NStorageTester{
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };



} /* TempLat */
#ifdef TEMPLATTEST
#include "CosmoInterface/evolvers/rk2nstorage_test.h"
#endif

#endif
