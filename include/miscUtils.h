#ifndef MISC_UTILS_H
#define MISC_UTILS_H

namespace miscUtils {
    bool fileExists(const char *fn);
    bool fileCopy(const char *fn_dst, const char *fn_src);
    void mkdir(const char *fn);
};

#endif
