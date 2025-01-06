#include "io/HDF5Writer.hpp"
#include <iostream>

namespace fluid
{

    HDF5Writer::HDF5Writer(const std::string &filename) : filename(filename)
    {
        file = H5::H5File(filename + ".h5", H5F_ACC_TRUNC);
    }

    void HDF5Writer::writeCoordinates(std::vector<Domain> &domains)
    {
        for (auto &domain : domains)
        {
            writeDomainCoordinates(domain);
        }
    }

    void HDF5Writer::writeDomainCoordinates(Domain &domain)
    {
        std::vector<double> node_coords;
        node_coords.reserve((domain.nx + 1) * (domain.ny + 1) * (domain.nz + 1) * 3);

        for (int i = 0; i < domain.nx + 1; i++)
        {
            for (int j = 0; j < domain.ny + 1; j++)
            {
                for (int k = 0; k < domain.nz + 1; k++)
                {
                    node_coords.push_back(domain.xOrigin + i * domain.boxDims);
                    node_coords.push_back(domain.yOrigin + j * domain.boxDims);
                    node_coords.push_back(domain.zOrigin + k * domain.boxDims);
                }
            }
        }

        std::cout << "writing domain " << domain.id << " with "
                  << node_coords.size() << " nodes" << std::endl;

        hsize_t dims[1] = {node_coords.size()};
        H5::DataSpace dataspace(1, dims);

        std::string dataset_name = "domain_" + std::to_string(domain.id) + "_coordinates";
        H5::DataSet dataset = file.createDataSet(dataset_name,
                                                 H5::PredType::NATIVE_DOUBLE,
                                                 dataspace);
        dataset.write(node_coords.data(), H5::PredType::NATIVE_DOUBLE);
    }

    void HDF5Writer::writeTimestep(std::vector<Domain> &domains,
                                   const double &time,
                                   const int &iteration)
    {
        std::string groupname = "timestep_" + std::to_string(iteration);
        H5::Group timestepGroup = file.createGroup(groupname);

        // Write the time value
        hsize_t dims[1] = {1};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet timeDataset = timestepGroup.createDataSet("time",
                                                              H5::PredType::NATIVE_DOUBLE,
                                                              dataspace);
        timeDataset.write(&time, H5::PredType::NATIVE_DOUBLE);

        // Write each domain's data
        for (auto &domain : domains)
        {
            std::string domainGroupname = "domain_" + std::to_string(domain.id);
            H5::Group domainGroup = timestepGroup.createGroup(domainGroupname);
            writeDomain(domain, domainGroup);
        }
    }

    void HDF5Writer::writeDomain(Domain &domain, H5::Group &domainGroup)
    {
        hsize_t dims[3] = {
            static_cast<hsize_t>(domain.nx),
            static_cast<hsize_t>(domain.ny),
            static_cast<hsize_t>(domain.nz)};
        H5::DataSpace dataspace(3, dims);

        // Write each field
        auto writeField = [&](const char *name, const std::vector<double> &data)
        {
            H5::DataSet dataset = domainGroup.createDataSet(name,
                                                            H5::PredType::NATIVE_DOUBLE,
                                                            dataspace);
            dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
        };

        writeField("rho", domain.rho_);
        writeField("u", domain.u_);
        writeField("v", domain.v_);
        writeField("w", domain.w_);
        writeField("p", domain.p_);

        // Write ghost cell mask with different type
        H5::DataSet ghostDataset = domainGroup.createDataSet("ghostCellMask",
                                                             H5::PredType::NATIVE_UINT8,
                                                             dataspace);
        ghostDataset.write(domain.ghostCellMask.data(), H5::PredType::NATIVE_UINT8);
    }

} // namespace fluid