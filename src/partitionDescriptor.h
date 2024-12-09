// @author: jbhayet
#ifndef __PARTITION_DESCRIPTOR__
#define __PARTITION_DESCRIPTOR__
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <memory>
#include <Eigen/Dense>

#define SMAX 30
#define HASH_CHARS 4

typedef Eigen::Matrix< unsigned int, Eigen::Dynamic, Eigen::Dynamic > 	MatrixXUL;

// This is the class for partition descriptions
class partitionDescriptor {
    uint64_t sum;        // Total sum of the histogram
    uint64_t length;     // Length of the histogram
    uint64_t data[SMAX];
  public:
    // Constructor from vector of data
    partitionDescriptor(const uint64_t *data_, uint64_t sz): sum(0),length(sz) {
        // Copy the data
        memcpy(this->data,data_,sz*sizeof(data_[0]));
        // Sum of the histogram elements
        for (unsigned int k=0;k<sz;k++)
            this->sum+=  data_[k]*(k+1);
    }

    // Constructor from a string describing the partition
    partitionDescriptor(const std::string &description_string, uint64_t sz): sum(0),length(sz) {
        // Get the values from the string
        std::stringstream ss(description_string);
        for (unsigned int k=0;k<sz;k++) {
          ss >> this->data[k];
          this->sum+=  this->data[k]*(k+1);
        }
    }

    // Constructor from an integer sz (size) and a vector of integers being a partition
    partitionDescriptor(const std::vector<unsigned int> &partition_list, uint64_t sz): sum(0),length(sz) {
      // Fill to zero
      memset(this->data,0,sz*sizeof(this->data[0]));
      for (auto n : partition_list)
        if (n>0 && n<sz+1) {
          this->data[n-1]++;
        }
      for (unsigned int k=0;k<sz;k++) sum+= this->data[k]*(k+1);
    }

    // Size of the descriptor
    inline uint64_t size() const {
        return this->length;
    }

    // Sum of the descriptor
    inline uint64_t get_sum() const {
        return this->sum;
    }

    // Operator []
    inline const uint64_t& operator[](int idx) const {
      return this->data[idx];
    }

    // Operator ==
    inline bool operator==(const partitionDescriptor &other) const {
        if (other.length!=this->length)
            return false;
        if (other.sum!=this->sum)
            return false;
        for (unsigned int k=0;k<this->length;k++)
            if (other.data[k]!= this->data[k])
                return false;
        return true;
    }

    // Check if the caller is inferior (lexicographically) to the callee
    inline bool operator<(const partitionDescriptor &other) const {
        if (other.length!=this->length)
            return false;
        for (unsigned int k=0;k<this->length;k++)
            if (this->data[k]>other[k])
                return false;
        return true;
    }

    // Count possible assignations to get this descriptor from picking elements in d_init
    inline uint64_t countPossibleAssignations(const partitionDescriptor &d_init,const MatrixXUL &combinationsTable) {
        uint64_t count = 1;
        for (int k=this->length-1;k>=0;k--)
          if (this->data[k]>0) {
            //
            count *= combinationsTable(d_init[k],this->data[k]);
        }
        return count;
      }

    // Removes some of the elements at one position
    inline partitionDescriptor remove(int position,int quantity) const {
        partitionDescriptor p = partitionDescriptor(this->data,this->length);
        p.data[position] -= quantity;
        return p;
    }

    // Shortens at one position
    inline partitionDescriptor shorten(int position) const {
      return partitionDescriptor(this->data,position);
    }

    // Difference between the caller and the callee (elementwise)
    inline partitionDescriptor difference(const partitionDescriptor &other) const {
        partitionDescriptor p = partitionDescriptor(this->data,this->length);
        for (unsigned int k=0;k<this->length;k++)
            p.data[k] -= other.data[k];
        return p;
    }

    // Difference between the caller and the callee (elementwise); then shorten
    inline partitionDescriptor differenceAndShorten(const partitionDescriptor &other,int k) const {
        partitionDescriptor p = partitionDescriptor(this->data,k);
        for (int i=0;i<k;i++) {
            p.data[i] -= other.data[i];
            p.sum     -= other.data[i]*(i+1);
        }
        return p;
    }

    // Check if the two descriptors correspond to the same total number
    inline bool compatible(const partitionDescriptor&other) const {
        if (other.length!=this->length)
            return false;
        return (this->sum==other.sum);
    }

    // Check if the calling descriptor can be a descendent of the first correspond to the same total number
    inline bool descendent(const partitionDescriptor&other) const {
      if (!compatible(other))
        return false;
      for (int k=0;k<(int)this->length;k++)
        if (this->data[k]!=other.data[k]) {
          if (this->data[k]>other.data[k])
            return false;
          else
            break;  
        }
      for (int k=(int)this->length-1;k>=0;k--)
        if (this->data[k]!=other.data[k]) {
          if (this->data[k]<other.data[k])
            return false;
          else
            break;  
        }
      return true;
    }

    // Determine the max key (for using hashing)
    static uint64_t maxKey() {
        uint64_t k = SMAX;
        for (unsigned int i=1;i<=HASH_CHARS;i++)
            k = k*SMAX+SMAX;
        return k;
    }

    // Determine a key (for using hashing)
    inline uint64_t key() const {
        if (this->length>HASH_CHARS)
            return -1;
        uint64_t k = this->data[0];
        for (unsigned int i=1;i<this->length;i++)
            k = k*SMAX+this->data[i];
        return k;
    }

    // Determine the highest index such that d[k]>initial.d[k] and d[j]==initial.d[j] for j>k
    // -1 means they are equal
    inline int highestDifferent(const partitionDescriptor&other) const {
      for (int k=this->length-1;k>=0;k--)
            if (this->data[k]>other.data[k])
                return k;
      return -1;
    }

    // Check whether a vector of integers being a partition can be held within this descriptor
    inline bool holds(const std::vector<unsigned int> &partition_list) const {
      for (auto n : partition_list)
        if (n>0 && n<this->length && this->data[n-1]<1)
          return false;
      return true;
    }

    friend std::ostream& operator<<(std::ostream& os, const partitionDescriptor& dt);
};

std::ostream& operator<<(std::ostream& os, const partitionDescriptor& pDsct) {
  for (unsigned int i=0;i<pDsct.length;i++)
    os << pDsct.data[i] << " | ";
  return os;
}
#endif
