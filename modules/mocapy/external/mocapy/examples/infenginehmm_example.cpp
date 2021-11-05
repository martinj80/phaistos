/*
 * infenginehmm_example.cpp
 *
 *  Copyright (C) 2008, Jes Frellsen, The Bioinformatics Centre, University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
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
#include "mocapy.h"
#include <time.h>

using namespace mocapy;
using namespace std;

void copy_flat(const vector<double*> & flat_iterator, double array[]) {
	for(uint i=0; i<flat_iterator.size(); i++) {
		*flat_iterator[i] = array[i];
	}
}

void works_not() {
}

int main(void) {
	time_t t = time(NULL);
	mocapy_seed(t);
	cout << "Mocapy seed = " << t << endl;

	cout << "******************" << endl;
	cout << "Setting up network" << endl;
	cout << "******************" << endl;


	// Setup nodes
	Node* in = NodeFactory::new_discrete_node(2, "in");
	Node* in2 = NodeFactory::new_discrete_node(2, "in2");
	Node* hd0 = NodeFactory::new_discrete_node(4, "hd0");
	Node* hd1 = NodeFactory::new_discrete_node(4, "hd1");
	Node* out = NodeFactory::new_discrete_node(4, "out");

	// Setup dbn
	vector<Node*> start_nodes = vec(in, in2, hd0, out);
	vector<Node*> end_nodes = vec(in, in2, hd1, out);

	DBN dbn = DBN(start_nodes, end_nodes);

	dbn.add_intra(0, 2);
	dbn.add_intra(1, 2);
	dbn.add_intra(2, 3);
	dbn.add_inter(2, 2);
	dbn.construct();

	// Output the CPDs
	cout << *in << endl;
	cout << *in2 << endl;
	cout << *hd0 << endl;
	cout << *hd1 << endl;
	cout << *out << endl;
	cout << endl;

	cout << "******************" << endl;
	cout << "Setting up dataset" << endl;
	cout << "******************" << endl;

	Sequence data;
	data.set_shape(5, 4);
	double seq_array[] = {1,1,0,1,
	                      0,0,0,1,
	                      0,0,0,1,
	                      0,0,0,1,
	                      0,0,0,1};
	copy_flat(data.flat_iterator(), seq_array);


	MDArray<eMISMASK> mism;
	mism.set_shape(5, 4);
	mism.set_wildcard(vec(-1, 0), MOCAPY_OBSERVED);
	mism.set_wildcard(vec(-1, 1), MOCAPY_OBSERVED);
	mism.set_wildcard(vec(-1, 2), MOCAPY_HIDDEN);
	mism.set_wildcard(vec(-1, 3), MOCAPY_OBSERVED);

	MDArray<eMISMASK> mism_sample;
	mism_sample.set_shape(5, 4);
	mism_sample.set_wildcard(vec(-1, 0), MOCAPY_OBSERVED);
	mism_sample.set_wildcard(vec(-1, 1), MOCAPY_OBSERVED);
	mism_sample.set_wildcard(vec(-1, 2), MOCAPY_HIDDEN);
	mism_sample.set_wildcard(vec(-1, 3), MOCAPY_HIDDEN);

	mism_sample.get(3,3) = MOCAPY_OBSERVED;

	cout << "Data: \n" << data << endl;
	cout << "Mism: \n" << mism << endl;
	cout << "Mism_sample: \n" << mism_sample << endl;

	cout << "*****************************" << endl;
	cout << "Example of SampleInfEngineHMM" << endl;
	cout << "*****************************" << endl;

	// Setup the sampler
	SampleInfEngineHMM sampler(&dbn, data, mism_sample, 2);

	MDArray<double> sample = sampler.sample_next();
	cout << "Sample:\n" << sample;
	cout << "ll of sample = " << sampler.calc_ll(mism) << endl << endl;

	cout << "undo()" << endl;
	sampler.undo();
	cout << "ll of initial values =" << sampler.calc_ll(mism) << endl << endl;

	cout << "Setting start=0 and end=1" << endl << endl;
	sampler.set_start_end(0, 1);

	sample = sampler.sample_next();
	cout << "Sample:\n" << sample;
	cout << "ll of sample = " << sampler.calc_ll(mism) << endl << endl;

	cout << "*********************************" << endl;
	cout << "Example of LikelihoodInfEngineHMM" << endl;
	cout << "*********************************" << endl;

	double ll;
	LikelihoodInfEngineHMM infengine(&dbn, 2);

	ll = infengine.calc_ll(data, mism);

	cout << "ll of initail values = " << ll << endl;

	// Clean up!
	delete in;
	delete hd0;
	delete hd1;
	delete out;

	return EXIT_SUCCESS;
}
