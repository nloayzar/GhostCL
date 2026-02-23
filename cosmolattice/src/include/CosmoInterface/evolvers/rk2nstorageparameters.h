#ifndef COSMOINTERFACE_EVOLVERS_RK2NSTORAGEPARAMETERS_H
#define COSMOINTERFACE_EVOLVERS_RK2NSTORAGEPARAMETERS_H
/* File created by Adrien Florio, 2021 */
/* Released under the MIT license, see LICENSE.md. */

#include "TempLat/util/tdd/tdd.h"
#include "CosmoInterface/evolvers/evolvertype.h"

namespace TempLat {


    /** \brief A class which stores coefficients for 2N-storage RK.
     *
     *
     * Unit test: make test-rk2nstorageparameters
     **/

    template<typename T>
    class RK2NStorageParameters {
    public:
        /* Put public methods here. These should change very little over time. */
        RK2NStorageParameters() {

        }


        static std::vector<T> getAs(EvolverType eType){
            std::vector<T> res;

            if(eType == RK3_4 or eType == RK3_4_A) res = {0, -0.7825460361923583, -2.042914325731225, -1.799337253940777};
            else if (eType == RK2) res = {0, -0.5}; 

            return res;
        }

        static std::vector<T> getBs(EvolverType eType){
            std::vector<T> res;

            if(eType == RK3_4 or eType == RK3_4_A) res = {0.06688758201974097, 2.876554598956719, 0.5534657361343982, 0.3912730180961791};
            else if (eType == RK2) res = {0.5, 1.0};

            return res;
        }

        static size_t getEmbeddedStep(EvolverType eType){
            if(eType == RK3_4_A){
                return 2;
            }else return 0;
        }

        static bool isAdaptative(EvolverType eType){
            if(eType == RK3_4_A){
                return true;
            }else return false;
        }


        static bool isRK2n(EvolverType ev){
            bool res = false;
            if(ev == RK2 or ev == RK3_4 or ev == RK3_4_A) res = true;
            return res;
        }


    private:
        /* Put all member variables and private methods here. These may change arbitrarily. */
    };

    struct RK2NStorageParametersTester{
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif
    };

} /* TempLat */



#ifdef TEMPLATTEST
#include "CosmoInterface/evolvers/rk2nstorage_test.h"
#endif

#endif
