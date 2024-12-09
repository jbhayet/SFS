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

// define the format you want, you only need one instance of this...
const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << " <option(s)>"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-n, NUM\tSpecify the number n from which the partitions are generated. Default: 20."
              << std::endl;
}


// To write the results in a CSV file
void writeToCSVfile(const string &name, const MatrixXUL &matrix) {
    std::ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
}

int main(int argc, char *argv[]) {
  // n is the maximum sum of the elements of the compositions
  int n=25;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
        show_usage(argv[0]);
        return 0;
    } else if ((arg == "-n")) {
        if (i + 1 < argc) { // Make sure we aren't at the end of argv!
            n = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
        } else { // Uh-oh, there was no argument to the destination option.
            std::cerr << "The -n option requires one argument." << std::endl;
            return 1;
        }
    } else {
      show_usage(argv[0]);
      return 1;
    }
  }

  // Definition of the alpha parameter, choose a value in (0,2)
  // double alpha = 20.0;
  // This object will be called for counting the partitions
  Counter ct(n);
  // Partition generation
  cout << "[INF] Generating partitions of n=" << n << std::endl;
  // Enumerate all the partitions [a_1,...,a_n] from n, such that sum_i i a_i = n. 
  // They will come in ascending lexicographical order
  std::list<std::vector<unsigned int> >partitionsOfN;
  ascPartition(n,partitionsOfN);

  // Add the trivial one because the algorithm does not give it as an output
  std::vector<unsigned int> trivialPartition; trivialPartition.push_back(n);
  partitionsOfN.push_back(trivialPartition);
  cout << "[INF] First partition" << std::endl;
  printPartition(partitionsOfN.front());

  // Number of partitions
  int dim = partitionsOfN.size();
  cout << "[INF] Number of partitions: " << dim << std::endl;

  // Esta es la parte que nos interesa, vamos a estudiar una cadena de Markov con valores en las composiciones
  // (lo que llamo composiciones es nuestra manera de representar las particiones con el número
  // de bloques de cada tamaño, en plan (2,1,1...) )

  cout << "[INF] Filling Rmatrix" << endl;
  // Definition of the state matrix, the rows Rmatrix[i] are compositions
  MatrixXi Rmatrix(dim,n);
  int i=0;
  // Enumerate all the generated partitions
  for (auto partition: partitionsOfN) {
    // Take the partition one by one
    // Count how many elements in the partition have the value j+1.
    for (int j=0; j<n; j++) {
      Rmatrix(i,j) = countElements(partition,j+1);
    }
    i++;
  }

  // Precompute all the sum(R[i])
  cout << "[INF] Computing rowwise sums" << endl;
  MatrixXi S = Rmatrix.rowwise().sum();

  // Constructs all the compositions
  std::cout << "[INF] Building the list of partition descriptors (compositions)" << std::endl;
  std::vector<partitionDescriptor> P;
  for (auto partition: partitionsOfN) {
    P.push_back(partitionDescriptor(partition,n));
  }
  auto t1 = Clock::now();

  // Precompute all the Cbr or read them from file
  std::cout << "[INF] Computing counts" << std::endl;
  MatrixXUL Combin(dim,dim);
  for (int j=0; j<dim; j++) {
    printBar((float)j/dim);
    ct.resetValues();
    for (int i=0; i<j; i++) {
      Combin(i,j)=ct.recursiveCount_DescBreak(P[i],P[j]);
    }
  }
  auto t2 = Clock::now();
  std::cout << std::endl;
  std::cout << "[INF] Took: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

  Counter::printCalls();

  std::string fileName;
  std::stringstream ss(fileName);
  ss << "Combin-" << setw(3) << setfill('0') << n << ".csv";
  writeToCSVfile(ss.str(),Combin);

  return 0;

//      df = pd.DataFrame(data=Combin.astype(float))
//      df.to_csv('outfile' + str(n) + '.csv', sep=' ', header=False, float_format='%.10f', index=False)
//  else
//  //    try:
//          df  =pd.read_csv('outfile' + str(n) + '.csv', delim_whitespace=True, header=None)
//          Combin=df.astype(float).values
//          computeCombin = False
//      except IOError as e:
//          print("Error in reading file")
//          print(e)
//  printCalls()
//  print(Combin)
//  print(Combin.min())
//  print(Combin.max())

}
