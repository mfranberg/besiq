#ifndef __BINARY_PAIR_H__
#define __BINARY_PAIR_H__

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>

#define PAIR_CUR_VERSION 0x5cf2d3f2
/**
 * Defines the header.
 */
#pragma pack(push, 1)
struct bpair_header
{
    /**
     * Version of the file format (in case it changes).
     */
    uint32_t version;

    /**
     * Indicates the file format.
     */
    uint32_t format;

    /**
     * The number of pairs stored in the file.
     */
    uint64_t num_pairs;

    /**
     * Length of the header string.
     */
    uint32_t header_length;
};
#pragma pack(pop)

/**
 * Pair file interface.
 */
class pairfile
{
public:
    virtual bool open(size_t split = 1, size_t num_splits = 1) = 0;
    virtual void close() = 0;
    virtual bool read(std::pair<std::string, std::string> &pair) = 0;
    virtual bool write(size_t snp1_id1, size_t snp2_id2) = 0;
    virtual size_t num_pairs() = 0;
    virtual ~pairfile(){ };
};

class tpairfile : public pairfile
{
public:
    tpairfile(const std::string &path, std::vector<std::string> m_snp_names, const char *mode);
    ~tpairfile();
    
    bool open(size_t split = 1, size_t num_splits = 1);
    void close();
    bool read(std::pair<std::string, std::string> &pair);
    bool write(size_t snp1_id1, size_t snp2_id2);
    size_t num_pairs();
private:
    std::string m_path;
    std::string m_mode;
    uint64_t m_num_pairs;
    std::istream *m_input;
    std::ostream *m_output;
    std::vector<std::string> m_snp_names;
    std::map<std::string,size_t> m_snp_to_index;
    uint64_t m_pairs_left;
};

/**
 * A class for reading and writing the pair file.
 */
class bpairfile : public pairfile
{
public:
    /**
     * This constructor is used when reading files.
     *
     * @param path Path to the pair file.
     */
    bpairfile(const std::string &path);

    /**
     * This constructor is used when writing files.
     */
    bpairfile(const std::string &path, const std::vector<std::string> &snp_names);

    ~bpairfile();

    bool open(size_t split = 1, size_t num_splits = 1);
    void close();

    const std::vector<std::string> & get_snp_names();

    bool read(std::pair<std::string, std::string> &pair);
    bool write(size_t snp_id1, size_t snp_id2);
    size_t num_pairs();

private:
    /* Path to the file */
    std::string m_path;

    /* Reading or writing */
    std::string m_mode;

    /* File pointer */
    FILE *m_fp;

    /* Header */
    bpair_header m_header;

    /* Names of the SNPs */
    std::vector<std::string> m_snp_names;

    /* 
     * Number of pairs left to read.
     */
    uint64_t m_pairs_left;
};

pairfile * open_pair_file(const std::string &path, const std::vector<std::string> &snp_names);
bool split_pair_file(const std::string &all_pairs, size_t num_splits, const std::string &output_path);

#endif /* End of __BINARY_PAIR_H__ */
