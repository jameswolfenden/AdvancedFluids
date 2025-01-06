#ifndef HDF5_WRITER_HPP
#define HDF5_WRITER_HPP

#include <string>
#include <vector>
#include "H5Cpp.h"
#include "domain/Domain.hpp"

namespace fluid
{

    class HDF5Writer
    {
    public:
        explicit HDF5Writer(const std::string &filename);

        void writeCoordinates(std::vector<Domain> &domains);
        void writeTimestep(std::vector<Domain> &domains, const double &time, const int &iteration);

    private:
        std::string filename;
        H5::H5File file;

        void writeDomainCoordinates(Domain &domain);
        void writeDomain(Domain &domain, H5::Group &domainGroup);
    };

} // namespace fluid

#endif // HDF5_WRITER_HPP