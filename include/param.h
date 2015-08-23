#ifndef PARAM_H
#define PARAM_H

namespace param 
{
    typedef map<string, vector<string> > ParamMap;
    extern ParamMap _map;

    void readFile(const char *);

    bool exist(const string &);

    void add(const string &);

    int getIntValue(const string &, int pos=0);
    double getDoubleValue(const string &, int pos=0);
    string getStringValue(const string &);
};

void Tokenize(const string &, vector<string> &, const string &); 

template<class T>
T string_to_val(const string &);

#endif
