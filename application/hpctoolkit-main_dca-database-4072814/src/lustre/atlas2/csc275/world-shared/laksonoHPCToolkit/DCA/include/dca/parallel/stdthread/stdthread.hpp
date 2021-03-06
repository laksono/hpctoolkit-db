// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class provides an interface for parallelizing with Pthreads.
//
// TODO: Finish sum methods.

#ifndef DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP
#define DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP

#include <iostream>
#include <vector>

#include "dca/parallel/util/threading_data.hpp"

namespace dca {
namespace parallel {

class stdthread {
public:
  stdthread() {}

  void execute(int num_threads, void* (*start_routine)(void*), void* arg) {
    fork(num_threads, start_routine, arg);
    join();
  }

  friend std::ostream& operator<<(std::ostream& some_ostream, const stdthread& this_concurrency);

private:
  static constexpr char parallel_type_str_[] = "stdthread";
  void fork(int num_threads, void* (*start_routine)(void*), void* arg) {

    m_num_threads = num_threads;

    #pragma omp parallel for num_threads(num_threads)
    for (int id = 0; id < num_threads; id++) {
      ThreadingData data;

      data.id = id;
      data.num_threads = num_threads;
      data.arg = arg;
      //
      start_routine( (void*) &data );
    }
  }

  void join() {
  }

  int m_num_threads;
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_STDTHREAD_STDTHREAD_HPP
