#include "mocapy.h"

using namespace mocapy;

class Test {
public:
	MDArray<double> array;

	Test(MDArray<double> &array) : array(vec<uint>(4)) {
		this->array = array.moveaxis(0, 1);
	}
};

int main(void) {
	MDArray<double> array1(vec<uint>(4,5,6));
	Test t(array1);
}
