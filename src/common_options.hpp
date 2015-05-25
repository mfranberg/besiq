#ifndef __COMMON_OPTION_H__
#define __COMMON_OPTION_H__

#include <vector>
#include <string>
#include <fstream>

#include <plink/plink_file.hpp>
#include <besiq/io/covariates.hpp>
#include <besiq/io/pairfile.hpp>
#include <besiq/io/resultfile.hpp>
#include <besiq/method/method.hpp>
#include <shared_ptr/shared_ptr.hpp>

#include <cpp-argparse/OptionParser.h>

struct common_options
{
    common_options(plink_file_ptr gf, genotype_matrix_ptr g, method_data_ptr d, pairfile *pf, resultfile *rf)
        : genotype_file( gf ),
          genotypes( g ),
          data( d ),
          pairs( pf ),
          result_file( rf )

    {
    }

    plink_file_ptr genotype_file;
    genotype_matrix_ptr genotypes;
    method_data_ptr data;
    shared_ptr<pairfile> pairs;
    shared_ptr<resultfile> result_file;
};

optparse::OptionParser create_common_options(const std::string &usage, const std::string &description, bool support_cov);

shared_ptr<common_options> parse_common_options(optparse::Values &options, const std::vector<std::string> &args);

#endif /* End of __COMMON_OPTION_H__ */
