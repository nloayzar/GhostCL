#ifndef COSMOINTERFACE_INITIALIZERS_THERMALFLUCTUATIONSGENERATOR_H
#define COSMOINTERFACE_INITIALIZERS_THERMALFLUCTUATIONSGENERATOR_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include "TempLat/lattice/field/field.h"
#include "TempLat/lattice/algebra/random/randomgaussianfield.h"
#include "TempLat/lattice/algebra/coordinates/wavenumber.h"
#include "CosmoInterface/initializers/fluctuationsgenerator.h"

namespace TempLat {

    /** \brief Thermal counterpart of FluctuationsGenerator for scalar singlet initial conditions.
     *
     *  This generator uses the same left/right moving wave construction as vacuum fluctuations,
     *  but rescales the mode amplitude according to a Bose-Einstein occupation number n_k,
     *  with factor sqrt(1 + 2 n_k).
     */
    template<typename T>
    class ThermalFluctuationsGenerator {
    public:
        template <class Model>
        ThermalFluctuationsGenerator(Model& model, T pLSide, std::string pSeed) :
            lSide(pLSide),
            baseSeed(pSeed)
        {
        }

        template<class Model>
        auto getThermalFluctuationsNorm(Model& model, Field<T> f, T mass2, T kCutOff, T temperature) const
        {
            FourierSite<Model::NDim> ntilde(f.getToolBox());
            auto k = ntilde.norm() * f.getKIR();
            auto omega = FluctuationsGenerator<T>::omega_k(k,mass2, f.toString());
            auto Hcut = heaviside(kCutOff - k);

            auto vacuumNorm = Hcut * (model.omegaStar / model.fStar * pow(lSide / pow<2>(f.getDx()), 1.5)) * pow(2 * omega, -0.5)/sqrt(2);

            Field<T> thermalNorm("thermalNorm_" + f.toString(), f.getToolBox(), f.getLatParams());
            thermalNorm = vacuumNorm;

            if(temperature > 0) {
                auto nk = 1 / (exp(omega / temperature) - 1);
                thermalNorm = vacuumNorm * sqrt(1 + 2 * nk);
            }

            return thermalNorm;
        }

        template<class Model>
        auto getNormedThermalFluctuations(Model& model, Field<T> f, T mass2, std::string mySeed, T kCutOff, T temperature) const
        {
            auto fFluctuationNorm = getThermalFluctuationsNorm(model, f, mass2, kCutOff, temperature);
            return fFluctuationNorm * RandomGaussianField<T>(baseSeed + mySeed + f.toString(), f.getToolBox());
        }

        template<class Model>
        void conjugateThermalFluctuations(Model& model, Field<T> f, Field<T> p, T mass2, T aDot, T kCutOff, T temperature) const
        {
            auto fLeft = getNormedThermalFluctuations(model, f, mass2, "Random left", kCutOff, temperature);
            auto fRight = getNormedThermalFluctuations(model, f, mass2, "Random right", kCutOff, temperature);

            f.inFourierSpace() = (fLeft + fRight) / sqrt(2);
            f.inFourierSpace().setZeroMode(0);

            FourierSite<Model::NDim> ntilde(f.getToolBox());
            auto k = ntilde.norm() * f.getKIR();
            auto omega = FluctuationsGenerator<T>::omega_k(k,mass2, f.toString());

            p.inFourierSpace() = Constants::I<T> * omega * (fLeft - fRight) / sqrt(2) - aDot * f.inFourierSpace();
            p.inFourierSpace().setZeroMode(0);
        }

    private:
        T lSide;
        std::string baseSeed;
    };

}

#endif
