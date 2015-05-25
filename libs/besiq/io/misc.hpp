#ifndef __MISC_H__
#define __MISC_H__

#include <string>
#include <vector>

std::vector<std::string> unpack_string(const char *snp_name_str);

std::string pack_string(const std::vector<std::string> &snp_names);

#endif /* End of __MISC_H__ */
