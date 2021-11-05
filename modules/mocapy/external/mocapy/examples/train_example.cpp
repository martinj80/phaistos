/*
 * train_example.cpp
 *
 *  Copyright (C) 2008, Jes frellsen, The Bioinformatics Centre, University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
 *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>

#include <unistd.h>
#include <ctime>

#include "mocapy.h"

using namespace mocapy;

int main(int argc, char *argv[]) {
	// Parse command line arguments
	std::cout << "Usage: train_dbn [dataset] [mismask] [hidden node size] [output dbn]" << std::endl;

	if (argc!=5) {
		std::cout << "You must supply 4 arguments" << std::endl;
		exit(1);
	}

	char* data_filename = argv[1];
	char* mism_filename = argv[2];
	unsigned int h_size = static_cast<unsigned int>(atoi(argv[3]));
	char* out_filename = argv[4];
	unsigned int seed = static_cast<unsigned int>(time(NULL)) + static_cast<unsigned int>(getpid());

	std::cout << "Data filename: " << data_filename << std::endl;
	std::cout << "Mism filename: " << mism_filename << std::endl;
	std::cout << "Hidden node size: " << h_size << std::endl;
	std::cout << "Output filename: " << out_filename << std::endl;
	std::cout << "Seed: " << seed << std::endl;

	// Set the Mocapy random seed
	mocapy_seed(seed);

	// Make the DBN and nodes
	DBN dbn;

	Node* h0 = NodeFactory::new_discrete_node(h_size, "H0");
	Node* h1 = NodeFactory::new_discrete_node(h_size, "H1");

	// Setup dbn architecture
	dbn.set_slices(vec(h0), vec(h1));

	dbn.add_inter("H0", "H1");

	dbn.construct();

	// Setup sampler and EM-Engine - and load training data
	GibbsRandom mcmc = GibbsRandom(&dbn);

	EMEngine em = EMEngine(&dbn, &mcmc);
	em.load_mismask(mism_filename);
	em.load_sequences(data_filename);

	// Do the training
	std::cout << std::endl << "Starting EM loop" << std::endl;
	for (unsigned int i = 0; i<100; i++) {
		em.do_E_step(5, 10);
		double ll = em.get_loglik();
		std::cout << "ll = " << ll << std::endl;
		em.do_M_step();
	}

	// Save the network
	std::cout << std::endl << "Saving dbn: " << out_filename << std::endl;
	dbn.save(out_filename);

	return EXIT_SUCCESS;
}
