#include "miscUtils.h"
#include <sys/stat.h>
#include <fstream>

/* miscUtils::fileExists 
 * miscUtils::fileCopy
   miscUtils::mkdir */

bool miscUtils::fileExists(const char *fn)
{
    struct stat file_stat;
    return (stat(fn, &file_stat) == 0);
}

bool miscUtils::fileCopy(const char *fn_dst, const char *fn_src)
{
    using namespace std;

    ofstream fdst(fn_dst, fstream::trunc|fstream::binary);
    ifstream fsrc(fn_src, fstream::binary);
    fdst <<  fsrc.rdbuf();
    return true;
}


void miscUtils::mkdir(const char *dir)
{
    ::mkdir(dir, 0775);
}
