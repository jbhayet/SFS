// @author: jbhayet
#include "partitionDescriptor.h"
#include <Eigen/Dense>

typedef Eigen::Matrix< unsigned int, Eigen::Dynamic, Eigen::Dynamic > 	MatrixXUL;

#define HASH_TABLE_SIZE 3000000
#define HASH_USE 1
uint64_t hash_table[HASH_TABLE_SIZE];
class Counter {
  bool debug;
  // values= -np.ones((100000000,1), dtype=int)
  // calls      = 0
  // shortened  = 0
  MatrixXUL countSplittingTable;
  MatrixXUL combinationsTable;
  static unsigned int calls;
  static unsigned int shortened;

public:
  // Constructor
  Counter(unsigned int n, bool dbg=false) : debug(dbg) {
    countSplittingTable = MatrixXUL::Zero(n+3,n+3);
    combinationsTable   = MatrixXUL::Zero(n+3,n+3);
    initCombinationsTable(n+3);
    initCountSplittingTable(n+3);
  }

  // Count the number of ways to group k*p elements into k sets of p elements
  inline uint64_t countSplitting(unsigned int k,unsigned int p) {
    uint64_t count = 1;
    //  \frac{1}{k!}\prod_i C(ip,i)
    for (unsigned int i =1;i<=k;i++) {
      count*= combinationsTable(i*p,p)/i;
    }
    return count;
  }

  // Pascal triangle for binomial coefficients
  void initCombinationsTable(int n) {
      for (int i=0;i<n;i++) {
          combinationsTable(i,0)=1;
          combinationsTable(i,i)=1;
          for (int j=1;j<i;j++) {
            combinationsTable(i,j)=combinationsTable(i-1,j-1)+combinationsTable(i-1,j);
          }
      }
    }

  // Pre-comute all the splitting counts
  // countSplittingTable[i][j] is the number of ways to group i*j elements into i sets of j elements
  void initCountSplittingTable(unsigned int n) {

    for (unsigned int i=0;i<n;i++)
      for (unsigned int j=0;j<n;j++)
        if (i*j<n)
          countSplittingTable(i,j) = countSplitting(i,j);
  }

  inline void resetValues() {
    memset(hash_table,0,HASH_TABLE_SIZE*sizeof(hash_table[0]));
  }

  static void printCalls() {
    std::cout << "[INF] Calls " << calls << " vs. " << shortened << std::endl;
  }
  // Descend-and-Break recursive algorithm
  unsigned int recursiveCount_DescBreak(const partitionDescriptor&d_init,
                                        const partitionDescriptor&d_end) {
    // global calls
    calls++;

    // If the init and end configurations are not compatible, this is a dead-end
    // (compatibility is evaluated by checking the sum of elements)
    if (d_end.compatible(d_init)==false) {
      if (debug)
        std::cout << "[DBG] Nope!" << std::endl;
      return 0;
    }
    // If the init and end configurations are the same, we found a path!
    if (d_end==d_init) {
      if (debug)
        std::cout << "[DBG] OK" << std::endl;
      return 1;
    }

    // If the computation has already been done, do not repeat it!
#if HASH_USE
    unsigned int key = d_init.key();
    if (key>0 && key<HASH_TABLE_SIZE && hash_table[key]>0) {
      shortened++;
      return hash_table[key];
    }
#endif
    // Find the highest k such that d_end[k]>d_init[k]
    // This supposes that the reverse of d_end is superior (in lexicographical order) than the reverse of d_init
    unsigned int k     = d_end.highestDifferent(d_init);
    unsigned int delta = d_end[k]-d_init[k];
    if (debug) {
      std::cout << "[DBG] Position of farthest difference " << k << std::endl;
    }
    // Generate differential partitions of (d_end[k]-d_init[k])*(k+1) with max. group size k
    std::list<std::vector<unsigned int> > diff_partitions;
    ascPartitionVariant((k+1)*delta,k,d_init,diff_partitions);

    // Count the ways to break d_end into d_init
    uint64_t count = 0;
    // Determine which of the partitions are compatible with the initial distribution
    for (auto diff_partition : diff_partitions) {
        // Convert the partition into a descriptor
        partitionDescriptor d_diff(diff_partition,d_init.size());
        if (debug) {
            std::cout << "[DBG] Describing partitions of " << (k+1)*delta << std::endl;
            std::cout << d_diff << std::endl;
            std::cout << "[DBG] with delta = " << delta << std::endl;
        }
        // TODO: should be more efficient at that point (most of the time is spent here)
        // If it is possible to assign d_diff
        if (d_diff<d_init) {
            // New version: Seems to be 20-30% faster
            partitionDescriptor d_init_remain = d_init.differenceAndShorten(d_diff,k);
            // TODO: VOID HARD COPY
            partitionDescriptor d_end_remain  = d_end.shorten(k);
            if (debug) {
                std::cout << "[DBG] Counting ways in forming: " << std::endl;
                std::cout << d_diff << std::endl;
                std::cout << "[DBG] from: " << std::endl;
                std::cout << d_init << std::endl;
            }
            // Note: this counts the possible ways in forming the
            // partition d_diff *from the elements of d_init*
            //uint64_t ns = countSplittingTable(delta,k+1);
            uint64_t ns = countSplitting(delta,k+1);
            if (debug) {
                std::cout << "[DBG] Splitting options " << ns << std::endl;
            }

            uint64_t nc = d_diff.countPossibleAssignations(d_init,combinationsTable);
            if (debug) {
              std::cout << "[DBG] Possible assignations " << nc << std::endl;
              std::cout << "[DBG] Valid partition" << std::endl;
              std::cout << "[DBG] New init:" << std::endl; 
              std::cout << "[DBG] " << d_init_remain << std::endl;
              std::cout << "[DBG] New end:" << std::endl;
              std::cout << "[DBG] " << d_end_remain << std::endl;
              std::cout << "[DBG] Count: " << ns*nc << std::endl;
            }
            uint64_t nsub = recursiveCount_DescBreak(d_init_remain,d_end_remain);
            if (debug) {
              std::cout << "[DBG] count from recursive call: " << nsub << std::endl;
            }
            count += ns*nc*nsub;
          }
        else
          if (debug)
            std::cout << "[DBG] Invalid partition" << std::endl;
    }
#if HASH_USE
    if (key>0 && key<HASH_TABLE_SIZE)
        hash_table[key]=count;
#endif
    return count;
  }
};

unsigned int Counter::calls = 0;
unsigned int Counter::shortened = 0;
