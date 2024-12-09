#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Dense>

#include "partitionCounting.h"
#include "counting.h"
#include "utils.h"


using namespace Eigen;
using namespace std;

#include <iomanip>
#include <sstream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

int main() {
    // Read the initial and final descriptors from a file
    string descriptor1;
    string descriptor2;
    getline(cin,descriptor1);
    getline(cin,descriptor2);
    partitionDescriptor d1(descriptor1,4);
    partitionDescriptor d2(descriptor2,4);
    if (!d2.descendent(d1)) {
        cout << "[ERR] The two partitions are not compatible" << endl;
        return 1;
    }
    std::cout << d1 << endl;
    std::cout << d2 << endl;
    // This object will be called for counting the partitions: n should be high enough for the binomial coefficients to be computed
    Counter ct(200,true);


    auto t1 = Clock::now();

    // Precompute all the Cbr or read them from file
    std::cout << "[INF] Computing counts" << std::endl;
    uint64_t nn = ct.recursiveCount_DescBreak(d1,d2);
    cout << "[INF] Total count: " << nn << endl;
    
    auto t2 = Clock::now();
    std::cout << std::endl;
    std::cout << "[INF] Took: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

    Counter::printCalls();

    return 0;
}
