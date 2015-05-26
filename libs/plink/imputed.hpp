#ifndef __IMPUTED_H__
#define __IMPUTED_H__

#include <vector>
#include <string>

#include <plink/snp_row.hpp>
#include <shared_ptr/shared_ptr.hpp>

typedef shared_ptr<std::vector<snp_row>> imputed_matrix_ptr;

struct imputed_info
{
    /**
     * Name of the variant.
     */
    std::string name;

    /**
     * Position of the variant.
     */
    unsigned long long pos;

    /**
     * Minor allele frequency.
     */
    float maf;

    /**
     * Info measure.
     */
    float info;
};

struct imputed_data
{
    /**
     * Genotypes.
     */
    imputed_matrix_ptr genotypes;

    /**
     * Samples.
     */
    std::vector<std::string> samples;

    /**
     * Variant info.
     */
    std::vector<imputed_info> info;
};

/**
 * Parses and hard calls the genotypes from an impute2 .gen file.
 *
 * @param path Path to the genotype file.
 * @param call_rate Threshold for calling a genotype.
 *
 * @return The genotype matrix.
 */
imputed_matrix_ptr parse_genotypes(const std::string &path, float call_rate, const std::vector<std::string> &samples);

/**
 * Parses the sample file from an impute2 .gen_samples file.
 *
 * @param path Path to the sample file.
 *
 * @return A list of parsed samples.
 */
std::vector<std::string> parse_sample(const std::string &path);

/**
 * Parses the info file from an impute2 .gen_info file.
 *
 * @param path Path to the info file.
 *
 * @return A list of parsed info for each variant.
 */
std::vector<imputed_info> parse_info(const std::string &path);

/**
 * Prases all of the imputed data associated with an impute2 
 * imputation.
 *
 * @param prefix The path prefix, it is assumed that prefix,
 *               prefix + '_sample' and prefix + '_info' exists.
 * @param call_rate Threshold for calling a genotype.
 *
 * @return Parsed data.
 */
imputed_data parse_imputed_data(const std::string &prefix, float call_rate);

#endif /* __IMPUTED_H__ */
