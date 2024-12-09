// @author: jbhayet
#include <list>
#include <iostream>
#include <vector>
#include "partitionDescriptor.h"

// TODO: this value should be updated
unsigned int l[180];

// Counts the number of elements in a partition with value k
unsigned int countElements(const std::vector<unsigned int> &partition,unsigned int k) {
  unsigned int c = 0;
  for (auto l: partition)
    if (l==k) c++;
  return c;
}

// Prints a partition
void printPartition(const std::vector<unsigned int> &partition) {
  for (auto n : partition)
    std::cout << n << " ";
  std::cout << std::endl;
}


// Recursive function for the partition of n
// n: the number we want to decompose into a partition
// m: minimal value for the elements of the partition; this allows for the partition to be described in ascending order
// k: index for the current element of the partition
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
    // When the recursive call is not the first one, we add the remainder to complete the currently built partition
    if (k>1) {
      // We cannot add anymore x's, so we add the remainder to
      // complete the partition
      l[k]=n;
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
// k: current higher index for the elements of the partition
void recAscVariant(unsigned int n,unsigned int p,unsigned int m,unsigned int k,const partitionDescriptor &d,std::list<std::vector<unsigned int> >&listOfPartitions) {
  if (m<1) return;
  if (m>p) return;
  if (m>n) return;
  // Start at m
  // Until x=p, tries all possible values of x (in increasing values)
  for (unsigned int x= m; x<=p; x++) if (d[x-1]>0) {
    // Get x into the description
    l[k]=x;
    // Get the partition of n-x, and will fill element k+1
    if (x<=n)
      recAscVariant(n-x,p,x,k+1,d,listOfPartitions);
  }
  // Here we cannot add more x's so we just add the remainder
  // but we first check if this remainder can be used in d
  if (d[n-1]<1) return;
  l[k]=n;
  // Adds the obtained partition to the list
  if (k>1 && n<=p) {
    std::vector<unsigned int> newPartition(l+1,l+k+1);
    listOfPartitions.push_back(newPartition);
  }
}

// Main function for acendent partition
void ascPartitionVariant(unsigned int n,unsigned int p, const partitionDescriptor &d, std::list<std::vector<unsigned int> > &listOfPartitions) {
  recAscVariant(n,p,1,1,d,listOfPartitions);
  return;
}
