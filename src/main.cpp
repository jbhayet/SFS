#include <iostream>
#include "partitionCounting.h"
#include "counting.h"
#include "utils.h"

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;


int main() {

  partitionDescriptor d_init("4 1 1 1 0 0 0",7);
  partitionDescriptor d_end("0 1 1 2 0 0 0",7);

  // n is the maximum sum of the elements of the compositions
  int n=20;
  // Definition of the alpha parameter, choose a value in (0,2)
  double alpha = 20.0;
  // This object will be called for counting the partitions
  Counter ct(n);

  // Display initial/final configurations
  cout << "--> Description of the initial configuration" << endl;
  cout << d_init << endl;
  cout << "--> Description of the final configuration" << endl;
  cout << d_end << endl;


  // Pre-compute combinations table
  //countSplittingTable  = initCountSplittingTable(n+1,combinationsTable)
  // initCountSplittingTable(n,unsigned int **combinationsTable);


  // Partition generation
  cout << "[INF] Generating partitions of n=" << n << std::endl;
  // Enumerate all the partitions of n. They will come in ascending lexicographical order
  std::list<std::vector<unsigned int> >partitionsOfN;
  ascPartition(n,partitionsOfN);

  // Add the trivial one because the algorithm does not give it as an output
  std::vector<unsigned int> trivialPartition; trivialPartition.push_back(n);
  partitionsOfN.push_back(trivialPartition);

  // Number of partitions
  int dim = partitionsOfN.size();
  cout << "[INF] Number of partitions:" << dim << std::endl;

  // Esta es la parte que nos interesa, vamos a estudiar una cadena de Markov con valores en las composiciones
  // (lo que llamo composiciones es nuestra manera de representar las particiones con el número
  // de bloques de cada tamaño, en plan (2,1,1...) )

  cout << "[INF] Filling Rmatrix" << endl;
  // Definition of the state matrix, the rows Rmatrix[i] are compositions
  // Rmatrix= np.zeros((dim,n), dtype=int)
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
  cout << Rmatrix << endl;

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
  MatrixXi Combin(dim,dim);
  bool computeCombin = true;
  if (computeCombin)
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
