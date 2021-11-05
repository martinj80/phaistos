#include "mocapy.h"

using namespace mocapy;

int main(void) {


	// Setup array
	MDArray<int> array1(vec<uint>(4,2,4));

	std::vector<int*> flat1 = array1.flat_iterator();
	for(uint i=0; i<flat1.size(); i++) *(flat1[i]) = i;

	// Print array
	std::cout << "Shape: " << array1.get_shape() << std::endl;
	std::cout << array1 << std::endl;

	// Try to permutate
	MDArray<int> new_array1 = array1.moveaxis(0, 2);

	std::cout << "New shape: " << new_array1.get_shape() << std::endl;
	std::cout << new_array1 << std::endl;

	// Setup array
	MDArray<int> array2(vec<uint>(4,2,2,4));

	std::vector<int*> flat2 = array2.flat_iterator();
	for(uint i=0; i<flat2.size(); i++) *(flat2[i]) = i;

	// Print array
	std::cout << "Shape: " << array2.get_shape() << std::endl;
	std::cout << array2 << std::endl;

	// Try to permutate
	MDArray<int> new_array2 = array2.moveaxis(0, 2);

	std::cout << "New shape: " << new_array2.get_shape() << std::endl;
	std::cout << new_array2 << std::endl;
}
