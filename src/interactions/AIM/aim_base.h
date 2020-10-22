#ifndef AIM_BASE_H
#define AIM_BASE_H

#include "expansion.h"
#include "grid.h"
#include "interactions/history_interaction.h"
#include "normalization.h"
#include "spacetime.h"

namespace AIM {
  class AimBase : public HistoryInteraction {
   public:
    AimBase(const std::shared_ptr<const DotVector> dots,
            const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>
                history,
            const int interp_order,
            const double c0,
            const double dt,
            const double h,
            std::shared_ptr<const Grid> grid,
            std::shared_ptr<const Expansions::ExpansionTable> expansion_table,
            Expansions::ExpansionFunction expansion_function,
            Expansions::ExpansionFunction expansion_function_fdtd,
            Normalization::SpatialNorm normalization)
        : HistoryInteraction(dots, history, interp_order, c0, dt),
          grid(std::move(grid)),
          expansion_table(std::move(expansion_table)),
          expansion_function(std::move(expansion_function)),
          expansion_function_fdtd(std::move(expansion_function_fdtd)),
          normalization(std::move(normalization)),
          h_(h),
          fdtd( h_? 1 : 0 )
    {
    }

   protected:
    std::shared_ptr<const Grid> grid;
    std::shared_ptr<const Expansions::ExpansionTable> expansion_table;
    Expansions::ExpansionFunction expansion_function;
    Expansions::ExpansionFunction expansion_function_fdtd;
    Normalization::SpatialNorm normalization;
 
    const Eigen::Vector3cd FDTD_Del_Del( // const std::vector<Eigen::Vector3cd> );
                                         const boost::multi_array<Eigen::Vector3cd, 1> );
		// const double omega_;
    const double h_;
    const bool fdtd;
  };
}

#endif
