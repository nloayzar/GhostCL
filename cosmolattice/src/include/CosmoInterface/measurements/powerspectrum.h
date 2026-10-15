#ifndef TEMPLAT_COSMOINTERFACE_POWERSPECTRUM_H
#define TEMPLAT_COSMOINTERFACE_POWERSPECTRUM_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2019

#include "TempLat/util/function.h"
#include "TempLat/lattice/algebra/helpers/getngrid.h"
#include "TempLat/lattice/algebra/helpers/getkir.h"
#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/measuringtools/radialprojector.h"
#include "CosmoInterface/runparameters.h"


namespace TempLat {


    /** \brief A class which computes the power spectrum, with the appropriate rescaling to make it volume independent.
     *
     *
     **/

     MakeException(WrongPSType);

     class PowerSpectrumMeasurer{
      public:
        template<typename T>
        PowerSpectrumMeasurer(const RunParameters<T>& par) :
        deltaKBin(par.deltaKBin),
        PSType(par.powerSpectrumType),
        PSVersion(par.powerSpectrumVersion)
        {}
        template< typename R>
        auto powerSpectrum(R f)  {
            return this->powerSpectrum(f, GetNGrid::get(f), GetKIR::getKIR(f));
        }



        // This function computes the power spectrum.
        // --> The normalization factor ensures that it recovers the appropriate expression in the continuum limit.
        //     This is discussed in Sect. 3 of arXiv:2006.15122.

        template<typename T>
        auto powerSpectrum(Field<T> f, ptrdiff_t N, T kIR)  {

            ptrdiff_t nDimensions = f.getToolBox()->mNDimensions;
            T nSites = std::pow(T(N), nDimensions);
            T dx = 2 * Constants::pi<T> / kIR / N;  // lattice spacing

            T maxRadius = std::sqrt(T(nDimensions)) / 2 * N;
            T kMaxBins = std::floor(maxRadius) + 1;
            ptrdiff_t nbins = std::floor(maxRadius / deltaKBin);
            T angularNormalization = 2 / (std::pow(T(2), nDimensions)
                * std::pow(Constants::pi<T>, T(nDimensions) / 2)
                * std::tgamma(T(nDimensions) / 2));
            T averageNormalization = angularNormalization * std::pow(kIR * dx, nDimensions) / nSites;
            T sumNormalization = kIR * dx / std::pow(T(N), 2 * nDimensions - 1) / 2 / Constants::pi<T>;

            if (PSVersion != 3){
                auto fk2 = projectRadiallyFourier(pow<2>(abs(f.inFourierSpace())), PSVersion == 1).measure(nbins, kMaxBins); //PSversion == true is a boolean. It tells the spectrum to use the centralValues, and not the bins, when rescaling.

                if (PSType == 2){
                  return Function(ntilde, averageNormalization * std::pow(ntilde, nDimensions)) * fk2;
                }
                else if (PSType == 1){
                  fk2.sumInsteadOfAverage();
                  return Function(ntilde, sumNormalization * ntilde) * fk2;
                }
                else{
                  throw(WrongPSType("You tried to call an undefined PSType " +std::to_string(PSType) + ", abort."));
                  return fk2; //To remove moot warning.
                }
            }
            else{
              WaveNumber ntilde(f.getToolBox()); //WaveNumber is the same as FourierSite, except it does not require to know the dimension at compile time.
              if (PSType == 2){

                  auto fk2 = projectRadiallyFourier(pow(ntilde.norm(), nDimensions) * pow<2>(abs(f.inFourierSpace())), false).measure(nbins, kMaxBins);
                  return averageNormalization * fk2;
              }
              else if (PSType == 1){
                auto fk2 = projectRadiallyFourier(ntilde.norm() * pow<2>(abs(f.inFourierSpace())), false).measure(nbins, kMaxBins);
                  fk2.sumInsteadOfAverage();
                  return sumNormalization * fk2;
              }
              else{
                throw(WrongPSType("You tried to call an undefined PSType " +std::to_string(PSType) + ", abort."));
                return projectRadiallyFourier(abs(f.inFourierSpace()), false).measure(nbins, kMaxBins); //To remove moot warning.
              }
            }
        }


        template<typename R, typename T>
        auto powerSpectrum(R f, ptrdiff_t N, T kIR) { //This function is for expression/composite operator, which need their own memory to perform the fourier transform.
            Field<T> tmp("tmp", GetToolBox::get(f));
            tmp = f;

            return this->powerSpectrum(tmp, N , kIR);
        }

      private:
        double deltaKBin;
        int PSType;
        int PSVersion;
     };


    class PowerSpectrumTester {
    public:
#ifdef TEMPLATTEST
        static inline void Test(TDDAssertion& tdd);
#endif

    };




} /* TempLat */

#ifdef TEMPLATTEST
#include "CosmoInterface/measurements/powerspectrum_test.h"
#endif


#endif
