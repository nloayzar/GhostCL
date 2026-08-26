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
#include <vector>

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
                stopGChiThreshold(par.stabilityEstimatorStopGChiThreshold),
                stopKPhiIntThreshold(par.stabilityEstimatorStopKPhiIntThreshold),
                kCutoff(par.kCutoff),
                rhoTot0(0),
                haveRhoTot0(false),
                haveInitialGhostSpectrum(false),
                stopRequested(false),
                output(filesManager.getWorkingDir() + "stability_estimator.txt",
                       amIRoot && enabled,
                       append ? std::ios_base::app : std::ios_base::out),
                spectralFrontOutput(filesManager.getWorkingDir() + "stability_spectral_front.txt",
                       amIRoot && enabled,
                       append ? std::ios_base::app : std::ios_base::out)
        {
            if (amIRoot && enabled && !append) {
                output << "#t max_phi2 max_chi2 max_pi_phi2 max_pi_chi2 max_grad_phi2 max_grad_chi2 "
                       << "G_chi K_phi_int G_chi_avg K_phi_int_avg tachyonic_indicator "
                       << "tachyonic_indicator_lambda tachyonic_indicator_lambda_avg "
                       << "G_chi_lambda G_chi_lambda_avg "
                       << "energy_conservation rho_tot\n";
                spectralFrontOutput << "#t k_peak_chi k_grow_10 k_grow_100 k_grow_1000 "
                                    << "k_grow_10_over_kCutOff k_grow_100_over_kCutOff "
                                    << "k_grow_1000_over_kCutOff\n";
                output.flush();
                spectralFrontOutput.flush();
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
            const T lambdaG = model.stabilityEstimatorLambdaG();
            const T tachyonicIndicator = q * maxPhi2 - muG2;
            const T tachyonicIndicatorAvg = q * avgPhi2 - muG2;
            const T tachyonicIndicatorLambda = max(
                    q * pow<2>(model.fldS(0_c))
                    - muG2
                    - static_cast<T>(3) * lambdaG * pow<2>(model.fldGS(0_c)));
            const T tachyonicIndicatorLambdaAvg =
                    q * avgPhi2 - muG2 - static_cast<T>(3) * lambdaG * avgChi2;
            const T gChi = model.dt * std::sqrt(std::max(tachyonicIndicator, static_cast<T>(0)));
            const T gChiAvg = model.dt * std::sqrt(std::max(tachyonicIndicatorAvg, static_cast<T>(0)));
            const T gChiLambda = model.dt * std::sqrt(std::max(tachyonicIndicatorLambda, static_cast<T>(0)));
            const T gChiLambdaAvg = model.dt * std::sqrt(std::max(tachyonicIndicatorLambdaAvg, static_cast<T>(0)));
            const T kPhiInt = model.dt * model.dt * q * maxChi2;
            const T kPhiIntAvg = model.dt * model.dt * q * avgChi2;

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
                   << gChi << " "
                   << kPhiInt << " "
                   << gChiAvg << " "
                   << kPhiIntAvg << " "
                   << tachyonicIndicator << " "
                   << tachyonicIndicatorLambda << " "
                   << tachyonicIndicatorLambdaAvg << " "
                   << gChiLambda << " "
                   << gChiLambdaAvg << " "
                   << energyConservation << " "
                   << rhoTot << "\n";
            output.flush();

            if (!areEstimatorValuesFinite(maxPhi2, maxChi2, gChi, kPhiInt,
                                          tachyonicIndicatorLambda, gChiLambda, rhoTot)) {
                stopRequested = true;
                return;
            }

            if ((stopGChiThreshold > 0 && gChi > stopGChiThreshold)
                || (stopKPhiIntThreshold > 0 && kPhiInt > stopKPhiIntThreshold)) {
                stopRequested = true;
            }
        }

        template <class Model>
        typename std::enable_if<!HasStabilityEstimatorParameters<Model>::value, void>::type
        measure(Model&, T, bool)
        {
        }

        template <class Model, class SpectrumMeasurer>
        typename std::enable_if<HasStabilityEstimatorParameters<Model>::value, void>::type
        measureSpectralFront(Model& model, T t, SpectrumMeasurer& PSMeasurer)
        {
            if (!enabled) return;

            auto spectrum = PSMeasurer.powerSpectrum(model.fldGS(0_c));
            if (!haveInitialGhostSpectrum || initialGhostSpectrum.size() != spectrum.size()) {
                initialGhostSpectrum.clear();
                initialGhostSpectrum.reserve(spectrum.size());
                for (auto&& bin : spectrum) {
                    initialGhostSpectrum.emplace_back(bin.getValue().average);
                }
                haveInitialGhostSpectrum = true;
            }

            T kPeak = 0;
            T pPeak = -std::numeric_limits<T>::max();
            T kGrow10 = 0;
            T kGrow100 = 0;
            T kGrow1000 = 0;

            for (size_t i = 0; i < spectrum.size(); ++i) {
                const T k = spectrum[i].getBin().average;
                const T p = spectrum[i].getValue().average;

                if (p > pPeak) {
                    pPeak = p;
                    kPeak = k;
                }

                const T pInitial = initialGhostSpectrum[i];
                if (pInitial > 0 && finiteIEEE(p) && finiteIEEE(pInitial)) {
                    const T ratio = p / pInitial;
                    if (ratio > 10) kGrow10 = k;
                    if (ratio > 100) kGrow100 = k;
                    if (ratio > 1000) kGrow1000 = k;
                }
            }

            spectralFrontOutput << std::setprecision(15)
                                << t << " "
                                << kPeak << " "
                                << kGrow10 << " "
                                << kGrow100 << " "
                                << kGrow1000 << " "
                                << ratioToKCutoff(kGrow10) << " "
                                << ratioToKCutoff(kGrow100) << " "
                                << ratioToKCutoff(kGrow1000) << "\n";
            spectralFrontOutput.flush();
        }

        template <class Model, class SpectrumMeasurer>
        typename std::enable_if<!HasStabilityEstimatorParameters<Model>::value, void>::type
        measureSpectralFront(Model&, T, SpectrumMeasurer&)
        {
        }

        bool shouldStop() const
        {
            return stopRequested;
        }

        T getStopGChiThreshold() const
        {
            return stopGChiThreshold;
        }

        T getStopKPhiIntThreshold() const
        {
            return stopKPhiIntThreshold;
        }

    private:
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

        bool areEstimatorValuesFinite(T maxPhi2, T maxChi2, T gChi, T kPhiInt,
                                      T tachyonicIndicatorLambda, T gChiLambda, T rhoTot) const
        {
            return finiteIEEE(maxPhi2)
                   && finiteIEEE(maxChi2)
                   && finiteIEEE(gChi)
                   && finiteIEEE(kPhiInt)
                   && finiteIEEE(tachyonicIndicatorLambda)
                   && finiteIEEE(gChiLambda)
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

        T ratioToKCutoff(T k) const
        {
            if (!std::isfinite(kCutoff) || kCutoff <= 0) return 0;
            return k / kCutoff;
        }

        bool amIRoot;
        bool enabled;
        bool expansion;
        bool fixedBackground;
        T stopGChiThreshold;
        T stopKPhiIntThreshold;
        T kCutoff;
        T rhoTot0;
        bool haveRhoTot0;
        bool haveInitialGhostSpectrum;
        std::vector<T> initialGhostSpectrum;
        bool stopRequested;
        OutputStream<T> output;
        OutputStream<T> spectralFrontOutput;
    };

}

#endif
