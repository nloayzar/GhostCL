#ifndef COSMOINTERFACE_MEASUREMENTS_STABILITYESTIMATORMEASURER_H
#define COSMOINTERFACE_MEASUREMENTS_STABILITYESTIMATORMEASURER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <ios>
#include <limits>
#include <type_traits>
#include <utility>

#include "CosmoInterface/definitions/averages.h"
#include "CosmoInterface/definitions/energies.h"
#include "CosmoInterface/definitions/fieldfunctionals.h"
#include "CosmoInterface/definitions/hubbleconstraint.h"
#include "CosmoInterface/measurements/measurementsIO/filesmanager.h"
#include "CosmoInterface/runparameters.h"
#include "TempLat/lattice/measuringtools/maximum.h"
#include "TempLat/util/conditionaloutput/outputstream.h"
#include "TempLat/util/rangeiteration/tagliteral.h"

namespace TempLat {

    template <typename Model>
    class HasStabilityEstimatorParameters {
    private:
        template <typename T>
        static auto test(int) -> decltype(
                std::declval<T>().stabilityEstimatorQ(),
                std::declval<T>().stabilityEstimatorMuG2(),
                std::declval<T>().stabilityEstimatorLambdaPhi(),
                std::declval<T>().stabilityEstimatorLambdaG(),
                std::true_type()
        );

        template <typename>
        static std::false_type test(...);

    public:
        static constexpr bool value = decltype(test<Model>(0))::value;
    };

    template <typename T>
    class StabilityEstimatorMeasurer {
    public:
        template <typename Model>
        StabilityEstimatorMeasurer(Model& model, FilesManager& filesManager, const RunParameters<T>& par, bool append) :
                amIRoot(model.getToolBox()->amIRoot()),
                enabled(par.stabilityEstimator && HasStabilityEstimatorParameters<Model>::value),
                expansion(par.expansion),
                fixedBackground(par.fixedBackground),
                stopThreshold(par.stabilityEstimatorStopThreshold),
                kMax2(computeKMax2<Model>(par)),
                rhoTot0(0),
                haveRhoTot0(false),
                stopRequested(false),
                output(filesManager.getWorkingDir() + "stability_estimator.txt",
                       amIRoot && enabled,
                       append ? std::ios_base::app : std::ios_base::out)
        {
            if (amIRoot && enabled && !append) {
                output << "#t max_phi2 max_chi2 max_pi_phi2 max_pi_chi2 max_grad_phi2 max_grad_chi2 "
                       << "C_phi C_chi C_max C_phi_avg C_chi_avg tachyonic_indicator "
                       << "energy_conservation rho_tot\n";
                output.flush();
            }
        }

        template <class Model>
        typename std::enable_if<HasStabilityEstimatorParameters<Model>::value, void>::type
        measure(Model& model, T t, bool isInitialTime)
        {
            if (!enabled) return;

            const T maxPhi2 = max(pow<2>(model.fldS(0_c)));
            const T maxChi2 = max(pow<2>(model.fldGS(0_c)));
            const T avgPhi2 = average(pow<2>(model.fldS(0_c)));
            const T avgChi2 = average(pow<2>(model.fldGS(0_c)));
            const T maxPiPhi2 = max(pow<2>(model.piS(0_c) * pow(model.aI, model.alpha - 3)));
            const T maxPiChi2 = max(pow<2>(model.piGS(0_c) * pow(model.aI, model.alpha - 3)));
            const T maxGradPhi2 = max(FieldFunctionals::grad2S(model, 0_c));
            const T maxGradChi2 = max(FieldFunctionals::grad2GS(model, 0_c));

            const T q = model.stabilityEstimatorQ();
            const T muG2 = model.stabilityEstimatorMuG2();
            const T lambdaPhi = model.stabilityEstimatorLambdaPhi();
            const T lambdaG = model.stabilityEstimatorLambdaG();

            const T cPhi2 = kMax2 + 1 + 3 * lambdaPhi * maxPhi2 + q * maxChi2;
            const T cChi2 = kMax2 + muG2 + 3 * lambdaG * maxChi2 - q * maxPhi2;
            const T cPhiAvg2 = kMax2 + 1 + 3 * lambdaPhi * avgPhi2 + q * avgChi2;
            const T cChiAvg2 = kMax2 + muG2 + 3 * lambdaG * avgChi2 - q * avgPhi2;

            const T cPhi = model.dt * std::sqrt(std::abs(cPhi2));
            const T cChi = model.dt * std::sqrt(std::abs(cChi2));
            const T cMax = std::max(cPhi, cChi);
            const T cPhiAvg = model.dt * std::sqrt(std::abs(cPhiAvg2));
            const T cChiAvg = model.dt * std::sqrt(std::abs(cChiAvg2));
            const T tachyonicIndicator = q * maxPhi2 - muG2;

            const T rhoTot = Energies::rho(model);
            const T energyConservation = computeEnergyConservation(model, rhoTot, isInitialTime);

            output << std::setprecision(15)
                   << t << " "
                   << maxPhi2 << " "
                   << maxChi2 << " "
                   << maxPiPhi2 << " "
                   << maxPiChi2 << " "
                   << maxGradPhi2 << " "
                   << maxGradChi2 << " "
                   << cPhi << " "
                   << cChi << " "
                   << cMax << " "
                   << cPhiAvg << " "
                   << cChiAvg << " "
                   << tachyonicIndicator << " "
                   << energyConservation << " "
                   << rhoTot << "\n";
            output.flush();

            if (!areEstimatorValuesFinite(maxPhi2, maxChi2, cPhi, cChi, cMax, rhoTot)) {
                stopRequested = true;
                return;
            }

            if (stopThreshold > 0 && cMax > stopThreshold) {
                stopRequested = true;
            }
        }

        template <class Model>
        typename std::enable_if<!HasStabilityEstimatorParameters<Model>::value, void>::type
        measure(Model&, T, bool)
        {
        }

        bool shouldStop() const
        {
            return stopRequested;
        }

        T getStopThreshold() const
        {
            return stopThreshold;
        }

    private:
        template <class Model>
        static T computeKMax2(const RunParameters<T>& par)
        {
            if (par.stabilityEstimatorUseKCutoff && std::isfinite(par.kCutoff)) {
                return par.kCutoff * par.kCutoff;
            }

            return 4 * static_cast<T>(Model::NDim) / (par.dx * par.dx);
        }

        template <class Model>
        T computeEnergyConservation(Model& model, T rhoTot, bool isInitialTime)
        {
            if (isInitialTime || !haveRhoTot0) {
                rhoTot0 = rhoTot;
                haveRhoTot0 = true;
            }

            if (fixedBackground) return 0;

            if (expansion) {
                auto hubbleLaw = HubbleConstraint::get(model);
                return hubbleLaw[0];
            }

            return std::abs(1.0 - rhoTot / rhoTot0);
        }

        bool areEstimatorValuesFinite(T maxPhi2, T maxChi2, T cPhi, T cChi, T cMax, T rhoTot) const
        {
            return finiteIEEE(maxPhi2)
                   && finiteIEEE(maxChi2)
                   && finiteIEEE(cPhi)
                   && finiteIEEE(cChi)
                   && finiteIEEE(cMax)
                   && finiteIEEE(rhoTot);
        }

        bool finiteIEEE(T x) const
        {
            std::uint64_t u;
            static_assert(sizeof(T) == sizeof(u), "StabilityEstimatorMeasurer expects IEEE-754 doubles.");
            std::memcpy(&u, &x, sizeof(T));
            const std::uint64_t exp = u & 0x7ff0000000000000ULL;
            return exp != 0x7ff0000000000000ULL;
        }

        bool amIRoot;
        bool enabled;
        bool expansion;
        bool fixedBackground;
        T stopThreshold;
        T kMax2;
        T rhoTot0;
        bool haveRhoTot0;
        bool stopRequested;
        OutputStream<T> output;
    };

}

#endif
