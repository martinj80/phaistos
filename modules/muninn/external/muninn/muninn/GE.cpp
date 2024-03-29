// GE.cpp
// Copyright (c) 2010-2012 Jes Frellsen
//
// This file is part of Muninn.
//
// Muninn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3 as
// published by the Free Software Foundation.
//
// Muninn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Muninn.  If not, see <http://www.gnu.org/licenses/>.
//
// The following additional terms apply to the Muninn software:
// Neither the names of its contributors nor the names of the
// organizations they are, or have been, associated with may be used
// to endorse or promote products derived from this software without
// specific prior written permission.

#include <deque>
#include <vector>

#include "muninn/GE.h"

namespace Muninn {

void GE::estimate_new_weights(const Binner *binner) {
    // Inform that new weights will be estimated
    MessageLogger::get().info("Estimating new weights.");
    MessageLogger::get().debug("Histogram shape: " + to_string<std::vector<unsigned int> >(current.get_shape()));

    // Bookkeeping
    total_iterations += current.get_n();

    // Tell the update scheme that the history is to be updated
    // TODO: Find a more elegant way of doing this.
    updatescheme->updating_history(current, *history);

    // Put the current histogram into the history
    history->add_histogram(current);

    try {
        // Estimate lnG from the data
        estimator->estimate(*history, *estimate, binner);

        // Log the current statistics
        force_statistics_log();

        // Make a new empty current histogram, with the newly estimated weights
        DArray new_weights = weightscheme->get_weights(*estimate, *history, binner);
        current = Histogram(new_weights); // TODO: Is this the best way of making an empty histogram?

        // TODO: Find a more elegant way of doing this.
        updatescheme->reset_prolonging();
    }
    catch (EstimatorException &exception){
        // Write warnings
        MessageLogger::get().warning(exception.what());
        MessageLogger::get().warning("Keeping old weights.");

        // Clean up the history and prolong the simulation time
        // TODO: Find a more elegant way of doing this.
        history->remove_newest();
        updatescheme->prolong();
    }
}

void GE::extend(const std::vector<unsigned int> &add_under, const std::vector<unsigned int> &add_over, const Binner *binner) {
    // TODO: Check that a maximum number of bins has not been exceeded

    // Extend the the current histogram and the history
    current.extend(add_under, add_over);
    history->extend(add_under, add_over);

    // Extend the estimate
    estimator->extend_estimate(*history, *estimate, add_under, add_over);

    // Set the new weights based on the new estimate
    current.set_lnw(weightscheme->get_weights(*estimate, *history, binner));
}

} // namespace Muninn
