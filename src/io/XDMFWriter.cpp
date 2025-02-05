#include "io/XDMFWriter.hpp"

namespace fluid
{

    XDMFWriter::XDMFWriter(std::vector<Domain> &domains,
                           const std::string &filename,
                           std::vector<double> &time)
        : filename(filename)
    {

        file.open(filename + ".xmf");
        writeHeader();
        writeGrids(domains, time);
        writeFooter();
        file.close();
    }

    void XDMFWriter::writeHeader()
    {
        file << "<?xml version=\"1.0\" ?>\n"
             << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
             << "<Xdmf Version=\"3.0\">\n"
             << "<Domain>\n";
    }

    void XDMFWriter::writeGrids(std::vector<Domain> &domains, std::vector<double> &time)
    {
        file << "<Grid Name=\"TimeSeries\" GridType=\"Collection\" "
             << "CollectionType=\"Temporal\">\n";

        for (size_t i = 0; i < time.size(); i++)
        {
            file << "<Grid Name=\"Timestep_" << i
                 << "\" GridType=\"Collection\" CollectionType=\"Spatial\">\n"
                 << "<Time Value=\"" << time[i] << "\"/>\n";

            writeTimestep(domains, i);

            file << "</Grid>\n";
        }

        file << "</Grid>\n";
    }

    void XDMFWriter::writeFooter()
    {
        file << "</Domain>\n"
             << "</Xdmf>\n";
    }

    void XDMFWriter::writeTimestep(std::vector<Domain> &domains, const int &iteration)
    {
        for (size_t i = 0; i < domains.size(); i++)
        {
            const auto &domain = domains[i];

            file << "<Grid Name=\"Domain_" << i
                 << "\" GridType=\"Uniform\">\n";

            // Write topology. Number of elements is +1 because the nodes are the vertices of the elements (e.g. 1x1x1 cube has 2x2x2 nodes)
            file << "<Topology TopologyType=\"3DSMesh\" "
                 << "NumberOfElements=\"" << domain.nx + 1 << " "
                 << domain.ny + 1 << " " << domain.nz + 1 << "\"/>\n";

            // Write geometry
            file << "<Geometry GeometryType=\"XYZ\">\n"
                 << "<DataItem Format=\"HDF\" "
                 << "Dimensions=\"" << (domain.nx + 1) * (domain.ny + 1) * (domain.nz + 1)
                 << " 3\" NumberType=\"Float\" Precision=\"8\">"
                 << filename << ".h5:/domain_" << i << "_coordinates</DataItem>\n"
                 << "</Geometry>\n";

            // Write attributes (rho, u, v, w, p)
            // The dimensions are the number of cells rather than nodes as the attributes are cell-centered
            auto writeAttribute = [&](std::string name)
            {
                file << "<Attribute Name=\"" << name
                     << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
                     << "<DataItem Format=\"HDF\" "
                     << "Dimensions=\"" << domain.nx << " "
                     << domain.ny << " " << domain.nz
                     << "\" NumberType=\"Float\" Precision=\"8\">"
                     << filename << ".h5:/timestep_" << iteration
                     << "/domain_" << i << "/" << name << "</DataItem>\n"
                     << "</Attribute>\n";
            };

            writeAttribute("rho");
            writeAttribute("u");
            writeAttribute("v");
            writeAttribute("w");
            writeAttribute("p");

            // Write ghost cell mask (UINT8 rather than float like the other attributes)
            file << "<Attribute Name=\"ghostCellMask\" "
                 << "AttributeType=\"Scalar\" Center=\"Cell\">\n"
                 << "<DataItem Format=\"HDF\" "
                 << "Dimensions=\"" << domain.nx << " "
                 << domain.ny << " " << domain.nz
                 << "\" NumberType=\"UInt8\" Precision=\"8\">"
                 << filename << ".h5:/timestep_" << iteration
                 << "/domain_" << i << "/ghostCellMask</DataItem>\n"
                 << "</Attribute>\n";

            file << "</Grid>\n";
        }
    }

} // namespace fluid