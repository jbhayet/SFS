// @author: jbhayet
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <cstring>
#include <iostream>
#include <memory>
#include <Eigen/Dense>


// This is the class for partition descriptions
class partitionDescriptor {
    unsigned int sum;
    unsigned int length;
    // Keep the pointer in a smart pointer structure
    std::shared_ptr<unsigned int> ptr;
    // This is where the data will be
    unsigned int *data;

  public:
    // Constructor from vector of data
    partitionDescriptor(const unsigned int *data_, unsigned int sz): sum(0),length(sz) {
        this->ptr   = std::shared_ptr<unsigned int>(new unsigned int[sz],std::default_delete<unsigned int[]>());
        this->data  = this->ptr.get();
        // Copy the data
        memcpy(this->data,data_,sz*sizeof(data_[0]));
        // Sum of the histogram elements
        for (unsigned int k=0;k<sz;k++)
            this->sum+=  data_[k]*(k+1);
    }

    // Constructor from a string describing the partition
    partitionDescriptor(const std::string &description_string, unsigned int sz): sum(0),length(sz) {
        // Create array to store the values
        this->ptr   = std::shared_ptr<unsigned int>(new unsigned int[sz],std::default_delete<unsigned int[]>());
        this->data  = this->ptr.get();
        // Get the values from the string
        std::stringstream ss(description_string);
        for (unsigned int k=0;k<sz;k++) {
          ss >> this->data[k];
          this->sum+=  this->data[k]*(k+1);
        }
    }

    // Constructor from an integer n (size) and a vector of integers being a partition
    partitionDescriptor(const std::vector<unsigned int> &partition_list, unsigned int sz): sum(0),length(sz) {
      // Create array to store the values
      this->ptr   = std::shared_ptr<unsigned int>(new unsigned int[sz],std::default_delete<unsigned int[]>());
      this->data  = this->ptr.get();
      // Fill to zero
      memset(this->data,0,sz*sizeof(this->data[0]));
      for (auto n : partition_list)
        if (n>0 && n<sz+1) {
          this->data[n-1]++;
        }
    }

    // Copy constructor
    partitionDescriptor(const partitionDescriptor& other) {
      this->ptr   = other.ptr;
      this->data  = this->ptr.get();
      this->sum   = other.sum;
      this->length= other.length;
    }

    // Assignment operator
    partitionDescriptor& operator=(const partitionDescriptor& other) {
      this->ptr   = other.ptr;
      this->data  = this->ptr.get();
      this->sum   = other.sum;
      this->length= other.length;
      return *this;
    }

    // Size of the descriptor
    inline int size() const {
        return this->length;
    }

    // Operator []
    inline const unsigned int& operator[](int idx) const {
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

    // Check if the caller is superior (elementwise) to the callee
    inline bool operator<(const partitionDescriptor &other) const {
        if (other.length!=this->length)
            return false;
        for (unsigned int k=0;k<this->length;k++)
            if (this->data[k]>other[k])
                return false;
        return true;
    }

    // Count possible assignations to get this descriptor from picking elements in d_init
    inline unsigned int countPossibleAssignations(const partitionDescriptor &d_init,const Eigen::MatrixXi &combinationsTable) {
      unsigned int count = 1;
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
        partitionDescriptor p = partitionDescriptor(this->data,this->length);
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

    // Determine a key (for using hashing)
    inline unsigned int key() const {
        if (this->length>5)
            return -1;
        unsigned int k = this->data[0];
        for (unsigned int i=1;i<this->length;i++)
            k = k*40+this->data[i];
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

    friend std::ostream& operator<<(std::ostream& os, const partitionDescriptor& dt);
};

std::ostream& operator<<(std::ostream& os, const partitionDescriptor& pDsct) {
  for (unsigned int i=0;i<pDsct.length;i++)
    os << pDsct.data[i] << " | ";
  return os;
}
