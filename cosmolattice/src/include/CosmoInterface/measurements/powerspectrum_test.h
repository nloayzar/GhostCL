#ifndef TEMPLAT_COSMOINTERFACE_POWERSPECTRUM_TEST_H
#define TEMPLAT_COSMOINTERFACE_POWERSPECTRUM_TEST_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2019

inline void TempLat::PowerSpectrumTester::Test(TempLat::TDDAssertion& tdd) {

    for (ptrdiff_t nDimensions : {2, 3}) {
        for (ptrdiff_t nGrid : {16, 24}) {
            auto toolBox = MemoryToolBox::makeShared(nDimensions, nGrid, 1);
            Field<double> field("spectrumTest", toolBox);
            field.inFourierSpace() = std::complex<double>(1, 0);

            double nSites = std::pow(double(nGrid), nDimensions);
            double angularMeasure = nDimensions == 2 ? 2 * Constants::pi<double> : 4 * Constants::pi<double>;
            double radialMoment = 0;
            for (ptrdiff_t site = 0; site < ptrdiff_t(nSites); ++site) {
                ptrdiff_t remaining = site;
                double radiusSquared = 0;
                for (ptrdiff_t dimension = 0; dimension < nDimensions; ++dimension) {
                    ptrdiff_t coordinate = remaining % nGrid;
                    remaining /= nGrid;
                    if (coordinate > nGrid / 2) coordinate -= nGrid;
                    radiusSquared += coordinate * coordinate;
                }
                radialMoment += std::pow(radiusSquared, double(nDimensions) / 2);
            }

            for (int spectrumType : {1, 2}) {
                for (int spectrumVersion : {1, 2, 3}) {
                    for (double binWidth : {0.5, 1.0}) {
                        ParameterParser parser(0, nullptr);
                        parser.addFromVector({"N=" + std::to_string(nGrid), "kIR=0.7", "dt=0.1", "baseSeed=1",
                            "PS_type=" + std::to_string(spectrumType), "PS_version=" + std::to_string(spectrumVersion),
                            "deltaKBin=" + std::to_string(binWidth)});
                        RunParameters<double> parameters(parser);
                        PowerSpectrumMeasurer measurer(parameters);
                        auto spectrum = measurer.powerSpectrum(field, nGrid, parameters.kIR);

                        ptrdiff_t expectedBins = std::floor(std::sqrt(double(nDimensions)) * nGrid / 2 / binWidth);
                        tdd.verify(spectrum.getCentralBinBounds().size() == size_t(expectedBins));
                        tdd.verify(AlmostEqual(spectrum.back().getBin().maxVal, std::sqrt(double(nDimensions)) * nGrid / 2));

                        double integral = 0;
                        bool normalized = true;
                        for (size_t bin = 0; bin < spectrum.size(); ++bin) {
                            auto& datum = spectrum[bin];
                            double radius = spectrumVersion == 1 ? spectrum.getCentralBinBounds()[bin] : datum.getBin().average;
                            if (spectrumType == 1) {
                                integral += datum.getValue().average / radius;
                            } else if (spectrumVersion != 3) {
                                double expected = angularMeasure * std::pow(radius, nDimensions) / (nSites * nSites);
                                normalized = normalized && AlmostEqual(datum.getValue().average / expected, 1.0);
                            } else {
                                integral += datum.getValue().average * 2 * datum.getBin().multiplicity;
                            }
                        }
                        if (spectrumType == 1) {
                            tdd.verify(AlmostEqual(integral * nSites * nSites / (nSites - 1), 1.0));
                        } else if (spectrumVersion != 3) {
                            tdd.verify(normalized);
                        } else {
                            tdd.verify(AlmostEqual(integral * nSites * nSites / (angularMeasure * radialMoment), 1.0));
                        }

                        auto compositeSpectrum = measurer.powerSpectrum(2 * field, nGrid, parameters.kIR);
                        bool compositeNormalized = compositeSpectrum.size() == spectrum.size();
                        for (size_t bin = 0; bin < compositeSpectrum.size() && compositeNormalized; ++bin) {
                            compositeNormalized = AlmostEqual(compositeSpectrum[bin].getValue().average
                                / spectrum[bin].getValue().average, 4.0);
                        }
                        tdd.verify(compositeNormalized);
                    }
                }
            }
        }
    }

}

#endif
