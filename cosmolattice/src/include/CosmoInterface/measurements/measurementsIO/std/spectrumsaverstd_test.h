#ifndef COSMOINTERFACE_MEASUREMENTS_STD_SPECTRUMSAVERSTD_TEST_H
#define COSMOINTERFACE_MEASUREMENTS_STD_SPECTRUMSAVERSTD_TEST_H
 
/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */ 
   
// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include "CosmoInterface/measurements/powerspectrum.h"

inline void TempLat::SpectrumSaverStdTester::Test(TempLat::TDDAssertion& tdd) {

    char directoryTemplate[] = "/tmp/cosmolattice-spectrum-saver-XXXXXX";
    char* directory = mkdtemp(directoryTemplate);
    tdd.verify(directory != nullptr);
    if (directory == nullptr) return;
    std::string outputDirectory = std::string(directory) + "/";
    std::string filename = outputDirectory + "spectra_test.txt";

    for (ptrdiff_t nDimensions : {2, 3}) {
        for (ptrdiff_t nGrid : {16, nDimensions == 2 ? 512 : 32}) {
            auto toolBox = MemoryToolBox::makeShared(nDimensions, nGrid, 1);
            Field<double> field("saverTest", toolBox);
            field.inFourierSpace() = std::complex<double>(1, 0);

            for (double binWidth : {0.5, 1.0}) {
                for (int verbosity : {0, 1, 2}) {
                    ParameterParser parser(0, nullptr);
                    parser.addFromVector({"N=" + std::to_string(nGrid), "kIR=0.1", "dt=0.1", "baseSeed=1",
                        "deltaKBin=" + std::to_string(binWidth), "spectraVerbosity=" + std::to_string(verbosity)});
                    RunParameters<double> parameters(parser);
                    PowerSpectrumMeasurer measurer(parameters);
                    auto spectrum = measurer.powerSpectrum(field, nGrid, parameters.kIR);
                    FilesManager files(outputDirectory, toolBox, false, false, false);
                    {
                        SpectrumSaverStd<double> saver(files, "test", true, false, parameters);
                        for (int snapshot = 0; snapshot < 2; ++snapshot) {
                            auto fieldSpectrum = std::make_shared<RadialProjectionResult<double>>(spectrum);
                            auto momentumSpectrum = std::make_shared<RadialProjectionResult<double>>(spectrum);
                            saver.save({fieldSpectrum, momentumSpectrum}, double(snapshot));
                        }
                    }

                    std::ifstream output(filename);
                    tdd.verify(output.is_open());
                    std::string line;
                    size_t rows = 0;
                    bool validOutput = true;
                    size_t expectedColumns = verbosity == 2 ? 14 : 4;
                    double maxMomentum = parameters.kIR * (std::floor(std::sqrt(double(nDimensions)) * nGrid / 2) + 1);
                    while (std::getline(output, line)) {
                        if (line.empty()) continue;
                        std::istringstream values(line);
                        double value;
                        size_t columns = 0;
                        while (values >> value) {
                            validOutput = validOutput && std::isfinite(value);
                            if (columns == 0 || (verbosity == 2 && (columns == 1 || columns == 3 || columns == 4))) {
                                validOutput = validOutput && value > 0 && value <= maxMomentum;
                            }
                            ++columns;
                        }
                        validOutput = validOutput && columns == expectedColumns;
                        ++rows;
                    }
                    tdd.verify(rows == 2 * spectrum.size());
                    tdd.verify(validOutput);
                }
            }
        }
    }

    std::remove(filename.c_str());
    rmdir(directory);

}

#endif
