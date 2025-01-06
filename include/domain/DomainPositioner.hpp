#ifndef DOMAIN_POSITIONER_HPP
#define DOMAIN_POSITIONER_HPP

#include <unordered_set>
#include "domain/Domain.hpp"

namespace fluid
{

    class DomainPositioner
    {
    public:
        // Constructor takes the root domain to start positioning from
        explicit DomainPositioner(Domain *root_domain);

    private:
        Domain *domain0;
        std::unordered_set<Domain *> passed;

        void findOffsets(Domain *domain);
    };

} // namespace fluid

#endif // DOMAIN_POSITIONER_HPP