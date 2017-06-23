/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014-2017 Colin MacLean <cmaclean@illinois.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "hapbin.hpp"
#include "hapmap.hpp"
#include "ihsfinder.hpp"
#include "calcselect.hpp"
#include "config.h"
#include "argparse.hpp"
#include "ehhpairfinder.hpp"
#include "ihspairfinder.hpp"

#if MPI_FOUND
#include <mpi.h>
#endif
#include <functional>
#include <cstdlib>

int main(int argc, char** argv)
{
#if MPI_FOUND
    MPI_Init(&argc, &argv);
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<const char*> hap('d', "hap", "Hap file", false, true, "");
    Argument<const char*> map('m', "map", "Map file", false, true, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('b', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<unsigned long> scale('s', "scale", "Gap scale parameter in pb, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<unsigned long long> window('w', "window", "Maximum gap in bp between two core SNP pairs.", false, false, 1000000);
    Argument<unsigned long> bufferSize('u', "buffer", "Maximum starting buffer size in SNPs. Will grow if necessary.", false, false, 2000);
    Argument<unsigned long> jobSize('j', "job-size", "Number of first SNPs to take at a time for each MPI process.", false, false, 100);
    Argument<unsigned int> numBins('i', "bins", "Number of frequency bins to use for iHS standardization", false, false, 50);
    Argument<unsigned int> windowBins('x', "window-bins", "Number of window bins to use for iHS standardization", false, false, 50);
    Argument<std::string> outfile('o', "out", "Output file", false, false, "out");
    Argument<bool> filternonpoly('f', "filter", "Filter out non-polymorphic loci", false, false);
    Argument<std::string> pops('p', "pop", "Filter by population", true, false,"");
    Argument<std::string> key('k', "key", "Population key", false, false, "");
    Argument<const char*> ranges('r', "ranges", "A flle with two space or tab separated columns specifying the start and end IDs of ranges to calculate", false, true, "");
    Argument<int> random_range(ArgumentBase::NO_SHORT_OPT, "random", "Randomly choose 1/x pairs to calculate", false, false, 1);
    ArgParse argparse({&help, &hap, &map, &outfile, &cutoff, &minMAF, &scale, &window, &bufferSize, &jobSize, &numBins, &windowBins, &filternonpoly, &pops, &key, &ranges, &random_range},
                      "Usage: ihs2range --map input.map --hap input.hap --ranges ranges.txt [--out outfile]");
    if (!argparse.parseArguments(argc, argv))
        return 1;

    std::size_t numSnps = HapMap::querySnpLength(hap.value());
    std::cout << "Chromosomes per SNP: " << numSnps << std::endl;

    PopKey *popkey = NULL;
    if ((pops.wasFound() && !key.wasFound()) || (!pops.wasFound() && key.wasFound()))
    {
        std::cerr << "--pop and --key must both be specified to filter by population." << std::endl;
        return 3;
    }
    if (pops.wasFound() && key.wasFound())
        popkey = new PopKey(key.value(), pops.values());

    
    std::ifstream list_file(ranges.value());
    std::string line;
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> jobs;
    while (std::getline(list_file, line)) {
        std::vector<std::string> split = splitString(line, ' ');
        if (split.size() == 1)
            split = splitString(line, '\t');
        if (split.size() == 2) {
            std::size_t start = hm.idToIndex(split[0].c_str());
            std::size_t range_end = hm.idToIndex(split[1].c_str());
            std::size_t job_size = jobSize.value();
            std::size_t job_start = start;
            for (; job_start+job_size < range_end; job_start+=job_size)
            {
                std::size_t job_end = job_start + job_size;
                r.emplace_back(std::make_tuple(job_start,job_end,range_end);
            }
            r.emplace_back(std::make_tuple(job_start,range_end,range_end);
        }
    }

    calcIhs2RangeMpi(hap.value(), map.value(), filternonpoly.value(), popkey, key.value(), pops.values(), cutoff.value(), minMAF.value(), scale.value(), window.value(),
                bufferSize.value(), jobSize.value(), numBins.value(), windowBins.value(), outfile.value(), jobs, random_range.value());

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
#else
    std::cout << "ERROR: ihs2bin must be built with MPI and MPIRPC!" << std::endl;
#endif
    return 0;
}


