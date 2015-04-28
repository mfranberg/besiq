#ifndef __RESULTFILE_H__
#define __RESULTFILE_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

#include <stdlib.h>
#include <stdio.h>

#define RESULT_CUR_VERSION 0x61248fc2

#pragma pack(push, 1)
struct result_header
{
    /**
     * Version number / magic number.
     */
    uint32_t version;

    /**
     * Format number (not used atm).
     */
    uint32_t format;

    /**
     * Length of the snp names part.
     */
    uint32_t snp_names_length;

    /**
     * Length of the col names part.
     */
    uint32_t col_names_length;

    /**
     * Number of pairs in the file.
     */
    uint64_t num_pairs;

    /**
     * Number of columns.
     */
    uint32_t num_float_cols;
};
#pragma pack(pop)

/**
 * Pure virtual class for
 */
class resultfile
{
    public:
        virtual ~resultfile(){ };
        /**
         * Open the result file, return true if successful or
         * false otherwise. Must be called prior to other calls.
         *
         * @return True if successfully opened, false otherwise.
         */
        virtual bool open() = 0;

        /**
         * Returns the number of pairs currently stored in the result
         * file.
         *
         * @return the number of pairs in the result file.
         */
        virtual uint64_t num_pairs() = 0;

        /**
         * Returns a list of column names stored in the result file.
         *
         * @return a list of column names stored in the result file.
         */
        virtual const std::vector<std::string> &get_header() = 0;

        /**
         * Returns a list of variant names stored in the result file.
         *
         * @return a list of variant names stored in the result file.
         */
        virtual const std::vector<std::string> &get_snp_names() = 0;

        /**
         * Set the column names.
         *
         * Important: This call overwrites any previously written
         * data to the file.
         *
         * @param header The new column names.
         *
         * @return True if successful, false otherwise.
         */
        virtual bool set_header(const std::vector<std::string> &header) = 0;

        /**
         * Read a pair from the file.
         *
         * @param pair The names of the variants will be written here.
         * @param values The values in the columns will be stored here.
         *
         * @return True if successful, false otherwise.
         */
        virtual bool read(std::pair<std::string, std::string> *pair, float *values) = 0;

        /**
         * Writes a pair to the file.
         *
         * @param snp1 Index of the first snp.
         * @param snp2 Index of the second snp.
         * @param values List of values to write.
         *
         * @return True if successful, false otherwise.
         */
        virtual bool write(const std::pair<std::string, std::string> &pair, float *values) = 0;

        /**
         * Closes the file.
         */
        virtual void close() = 0;
};

/**
 * A binary file format for results.
 */
class bresultfile : public resultfile
{
    public:
        /**
         * Constructor for reading.
         *
         * @param path Path to the input file.
         */
        bresultfile(const std::string &path);

        /**
         * Constructor for writing.
         *
         * @param path Path to the output file.
         * @param snp_names A list of names for each snp.
         */
        bresultfile(const std::string &path, const std::vector<std::string> &snp_names);

        /**
         * Destructor.
         */
        ~bresultfile();

        /**
         * @see resultfile::open.
         */
        bool open();
        
        /**
         * @see resultfile::close.
         */
        void close();

        /**
         * @see resultfile::read.
         */
        bool read(std::pair<std::string, std::string> *pair, float *values);

        /**
         * @see resultfile::write.
         */
        virtual bool write(const std::pair<std::string, std::string> &pair, float *values);

        /**
         * @see resultfile::num_pairs.
         */
        uint64_t num_pairs();

        /**
         * @see resultfile::get_header.
         */
        const std::vector<std::string> &get_header();
        
        /**
         * @see resultfile::get_snp_names.
         */
        const std::vector<std::string> &get_snp_names();

        /**
         * @see resultfile::set_header.
         */
        bool set_header(const std::vector<std::string> &header);

        /**
         * Returns true if the file seems corrupted.
         *
         * @return True if the file seems corrupted.
         */
        bool is_corrupted();

    private:
        /**
         * Read or writing mode.
         */
        std::string m_mode;

        /**
         * Path to the output file.
         */
        std::string m_path;

        /**
         * Underlying file pointer.
         */
        FILE *m_fp;
        
        /**
         * File header (that has been parsed or will be written).
         */
        result_header m_header;

        /**
         * List of names of the columns in the result file.
         */
        std::vector<std::string> m_col_names;

        /**
         * List of names of the variants in the result file.
         */
        std::vector<std::string> m_snp_names;

        /**
         * Maps snp names to indices in m_snp_names.
         */
        std::map<std::string, size_t> m_snp_to_index;
};

/**
 * Text results.
 */
class tresultfile : public resultfile
{
    public:
        /**
         * Constructor for reading.
         *
         * @param path Path to the input file.
         * @param mode Reading or writing, "r" or "w".
         * @param snp_names List of variants (required for interface).
         */
        tresultfile(const std::string &path, const std::string &mode, const std::vector<std::string> &snp_names);

        /**
         * Destructor.
         */
        ~tresultfile();

        /**
         * @see resultfile::open.
         */
        bool open();
        
        /**
         * @see resultfile::close.
         */
        void close();

        /**
         * @see resultfile::read.
         */
        bool read(std::pair<std::string, std::string> *pair, float *values);

        /**
         * @see resultfile::write.
         */
        virtual bool write(const std::pair<std::string, std::string> &pair, float *values);

        /**
         * @see resultfile::num_pairs.
         */
        uint64_t num_pairs();

        /**
         * @see resultfile::get_header.
         */
        const std::vector<std::string> &get_header();
        
        /**
         * @see resultfile::get_snp_names.
         */
        const std::vector<std::string> &get_snp_names();

        /**
         * @see resultfile::set_header.
         */
        bool set_header(const std::vector<std::string> &header);

    private:
        /**
         * Read or writing mode.
         */
        std::string m_mode;

        /**
         * Path to the output file.
         */
        std::string m_path;

        /**
         * Underlying file pointer.
         */
        std::istream *m_input;
        std::ostream *m_output;

        /**
         * Number of pairs written
         */
        uint64_t m_num_pairs;
        bool m_written;

        /**
         * List of names of the columns in the result file.
         */
        std::vector<std::string> m_col_names;

        /**
         * List of names of the variants in the result file.
         */
        std::vector<std::string> m_snp_names;
};

/**
 * Opens a result file (binary or text) and returns a pointer to it.
 *
 * @param path Path to the file.
 * @param snp_names Name of each variant to be included.
 *
 * @return A result file, or a null pointer if it could not be opened.
 */
resultfile * open_result_file(const std::string &path, const std::vector<std::string> &snp_names);

/**
 * Returns a missing value.
 *
 * @return a missing value.
 */
float result_get_missing();

#endif /* End of __RESULTFILE_H__ */
