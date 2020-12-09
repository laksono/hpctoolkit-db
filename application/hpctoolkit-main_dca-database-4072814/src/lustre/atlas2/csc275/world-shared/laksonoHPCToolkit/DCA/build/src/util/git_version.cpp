// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements the GitVersion.

#include "dca/util/git_version.hpp"
#include <iostream>
#include <string>

namespace dca {
namespace util {
// dca::util::

const std::string GitVersion::git_log = "commit 5ad6518c21b118eae4ab04651ae05ae99383dafe\nMerge: 8e3aef82 1bea44c7\nAuthor: Raffaele Solcà <rasolca@cscs.ch>\nDate:   Wed May 23 16:50:23 2018 +0200\n\n    Merge pull request #149 from eth-cscs/fix_mpi\n    \n    Fix MPI initialization with threading.\n";
const std::string GitVersion::git_status = " M include/dca/parallel/stdthread/stdthread.hpp\n M src/parallel/stdthread/stdthread.cpp\n";

void GitVersion::print() {
  std::cout << "\n"
            << "********************************************************************************\n"
            << "**********                        Git Version                         **********\n"
            << "********************************************************************************\n"
            << "\n"
            << "Last commit:\n"
            << GitVersion::git_log << "\n"
            << "Working tree:\n"
            << GitVersion::git_status << std::endl;
}

std::string GitVersion::string() {
  return std::string("Last commit:\n" + GitVersion::git_log + "\nWorking tree:\n" +
                     GitVersion::git_status);
}

}  // util
}  // dca
