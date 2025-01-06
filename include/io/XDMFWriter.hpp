#ifndef XDMF_WRITER_HPP
#define XDMF_WRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include "domain/Domain.hpp"

namespace fluid
{

    class XDMFWriter
    {
    public:
        XDMFWriter(std::vector<Domain> &domains,
                   const std::string &filename,
                   std::vector<double> &time);

    private:
        std::string filename;
        std::ofstream file;

        void writeHeader();
        void writeGrids(std::vector<Domain> &domains, std::vector<double> &time);
        void writeFooter();
        void writeTimestep(std::vector<Domain> &domains, const int &iteration);
    };

} // namespace fluid

#endif // XDMF_WRITER_HPP