#ifndef __CORRECT_H__
#define __CORRECT_H__

#include <vector>

#include <plink/plink_file.hpp>
#include <besiq/method/method.hpp>

class metaresultfile;

/**
 * Holds different options used in the multiple testing
 * correction. Different corrections may use different subsets
 * of options.
 */
struct correction_options
{
    /**
     * Significance threshold.
     */
    float alpha;

    /**
     * Number of tests.
     */
    std::vector<uint64_t> num_tests;

    /**
     * Weight.
     */
    std::vector<float> weight;

    /**
     * Normal, binomial?
     */
    std::string model;
};

void run_bonferroni(metaresultfile *result, float alpha, uint64_t num_tests, size_t column, const std::string &output_path);
void run_top(metaresultfile *result, float alpha, uint64_t num_top, size_t column, const std::string &output_path);
void run_static(metaresultfile *result, genotype_matrix_ptr genotypes, method_data_ptr data, const correction_options &options, const std::string &output_path);
void run_adaptive(metaresultfile *result, genotype_matrix_ptr genotypes, method_data_ptr data, const correction_options &options, const std::string &output_path);

#endif /* End of __CORRECT_H__ */
