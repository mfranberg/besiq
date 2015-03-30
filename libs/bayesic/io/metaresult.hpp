#ifndef __METARESULT_H__
#define __METARESULT_H__

#include <vector>
#include <string>

#include <stdexcept>

class resultfile;

class result_open_error: public std::exception
{
public:
    /**
     * Constructor.
     *
     * @param value The value that caused an error when computing a p-value.
     */
    result_open_error(const std::string &message)
        : m_message( message )
    {
    }

    /**
     * Destructor.
     */
    virtual ~result_open_error() throw()
    {
    }

    /**
     * Returns the error message.
     *
     * @return the error message.
     */
    virtual const char* what() const throw()
    {
        return m_message.c_str( );
    }

private:
    /**
     * The error message.
     */
    std::string m_message;
};

class metaresultfile
{
public:
    metaresultfile(const std::vector<resultfile *> &result_files);
    bool read(std::pair<std::string, std::string> *pair, float *value);
    uint64_t num_pairs();
    std::vector<std::string> get_header();
    std::vector<std::string> get_snp_names();

private:
    std::vector<resultfile *> m_results;
    size_t m_cur_file;
};

std::vector<resultfile *> open_result_files(const std::vector<std::string> &paths);

metaresultfile *open_meta_result_file(const std::vector<std::string> &paths);

#endif
