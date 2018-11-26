#include "../PDE/include/Randomizer.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "../C++ libs/eigen/Eigen/Dense"
#include "../C++ libs/eigen/Eigen/Sparse"
#include "../C++ libs/eigen/Eigen/Core"

#include "../C++ libs/eigen/Eigen/IterativeLinearSolvers"
using namespace Eigen;
// constructing atomics
#include <iostream>       // std::cout
#include <atomic>         // std::atomic, std::atomic_flag, ATOMIC_FLAG_INIT
#include <thread>         // std::thread, std::this_thread::yield
#include <vector>         // std::vector

std::atomic<bool> ready(false);
std::atomic_flag winner = ATOMIC_FLAG_INIT;

void count1m(int id) {
	while (!ready) { std::this_thread::yield(); }      // wait for the ready signal
	for (volatile int i = 0; i < pow(10,8); ++i) {}          // go!, count to 1 million
	//if (!winner.test_and_set()) { std::cout << "thread #" << id << " won!\n"; }
};

int main()

{
	Timer timer;
	int t0 = timer.get_milliseconds();
	std::vector<std::thread> threads;
	std::cout << "spawning 10 threads that count to 1 million...\n";
	for (int i = 1; i <= 10; ++i) threads.push_back(std::thread(count1m, i));
	ready = true;
	for (auto& th : threads) th.join();
	cout << "duration: " << timer.get_milliseconds() - t0 << endl;
	t0 = timer.get_milliseconds();
	for (int n = 0; n < 10; n++) {
		for (int i = 0; i < pow(10, 8); ++i) {}
	}
	
	cout << "duration: " << timer.get_milliseconds() - t0 << endl;

	return 0;
}



