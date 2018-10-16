#include <math.h>
#include <array>
#include <map>

double limit_decimals(double number, int decimals){
    double N = pow(10, decimals);
    return double(int(number * N)) / N;
}

typedef map< array<int, 2>, int> Map;
typedef map< array<int, 2>, int>::const_iterator MapIter;

void show_map(map< array<int, 2>, int> &m){
	for (MapIter iter = m.begin(); iter != m.end(); iter++){
		cout << "Key: " << iter->first[0] << ", " << iter->first[1] << endl << "Value:" << endl;
		cout << iter->second << endl;
	}
}

