// @author: jbhayet
#include <list>
#include <vector>

// TODO: this value should be updated
unsigned int l[180];

unsigned int countElements(const std::vector<unsigned int> &partition,unsigned int k) {
  unsigned int c = 0;
  for (auto l: partition)
    if (l==k) c++;
  return c;
}

// Recursive function for the partition of n
void recAsc(unsigned int n,unsigned int m,unsigned int k,std::list<std::vector<unsigned int> > &listOfPartitions) {
    if (m<1)
      return;
    if (m>n)
      return;
    // Start at m
    unsigned int x     = m;
    // Until x=n/2
    while (2*x<=n) {
        // Get x into the description
        l[k]=x;
        // Get the partition of n-x, and will fill element k+1
        recAsc(n-x,x,k+1,listOfPartitions);
        // Increment x
        x = x+1;
    }
    // We cannot add anymore x's, so we add the remainder to
    // complete the partition
    l[k]=n;
    // Just adds it
    if (k>1) {
      std::vector<unsigned int> newPartition(l+1,l+k+1);
      listOfPartitions.push_back(newPartition);
    }
}

// Main function for generating all the partitions of n in ascendent order
void ascPartition(unsigned int n, std::list<std::vector<unsigned int> > &listOfPartitions) {
    return recAsc(n,1,1,listOfPartitions);
}


// Recursive function for the partition of a number n, with elements of max size p
// n: the number we want to decompose into a partition
// m: minimal value for the elements of the partition
// p: maximal value for the elements of the partition
// k: current number
void recAscVariant(unsigned int n,unsigned int p,unsigned int m,unsigned int k,std::list<std::vector<unsigned int> >&listOfPartitions) {
  if (m<1) return;
  if (m>p) return;
  if (m>n) return;
  // Start at m
  unsigned int x     = m;
  // Until x=p, tries all possible values of x (in increasing values)
  while (x<=p) {
    // Get x into the description
    l[k]=x;
    // Get the partition of n-x, and will fill element k+1
    // TODO:
    if (x<=n)
      recAscVariant(n-x,p,x,k+1,listOfPartitions);
    // Increment x
    x++;
  }
  // Here we cannot add more x's so we just add the remainder
  l[k]=n;
  // Just adds the obtained partition to the list
  if (k>1 && n<=p) {
    std::vector<unsigned int> newPartition(l+1,l+k+1);
    listOfPartitions.push_back(newPartition);
  }
}

// Main function for acendent partition
void ascPartitionVariant(unsigned int n,unsigned int p, std::list<std::vector<unsigned int> > &listOfPartitions) {
  recAscVariant(n,p,1,1,listOfPartitions);
  return;
}
