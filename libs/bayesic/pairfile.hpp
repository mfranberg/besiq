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
    uint32_t num_pairs;

    /**
     * Length of the header string.
     */
    uint32_t header_length;
};

/**
 * Pair file interface.
 */
class pairfile
{
public:
    virtual bool open() = 0;
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
    
    bool open();
    void close();
    bool read(std::pair<std::string, std::string> &pair);
    bool write(size_t snp1_id1, size_t snp2_id2);
    size_t num_pairs();
private:
    std::string m_path;
    std::string m_mode;
    std::istream *m_input;
    std::ostream *m_output;
    std::vector<std::string> m_snp_names;
    std::map<std::string,size_t> m_snp_to_index;
};

/**
 * A class for reading and writing the pair file.
 */
class bpairfile : public pairfile
{
public:
    /**
     * This constructor is used when reading files.
     */
    bpairfile(const std::string &path);

    /**
     * This constructor is used when writing files.
     */
    bpairfile(const std::string &path, const std::vector<std::string> &snp_names);

    ~bpairfile();

    bool open();
    void close();

    const std::vector<std::string> & get_snp_names();

    bool read(std::pair<std::string, std::string> &pair);
    bool write(size_t snp_id1, size_t snp_id2);
    size_t num_pairs();

private:
    std::vector<std::string> parse_snp_names(const char *snp_name_str);
    std::string make_snp_names(const std::vector<std::string> &snp_names);

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
};

pairfile * open_pair_file(const std::string &path, const std::vector<std::string> &snp_names);

#endif /* End of __BINARY_PAIR_H__ */
