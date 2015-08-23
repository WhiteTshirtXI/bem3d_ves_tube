#include "cxxheaders.h"
#include "param.h"
#include "mympi.h"

/* void param::readFile 
 * bool param::exist 
 * void param::add 
 * int param::getIntValue 
 * double param::getDoubleValue 
 * string param::getStringValue 
 * void Tokenize 
 * tempalte<class T> T string_to_val */

namespace param {
    ParamMap _map;
};

// Read parameters from a file
void param::readFile(const char *fn)
{
    ifstream file;

    file.open(fn);

    if (file.fail()) {
	if (mympi::comm_rank() == 0) {
	    printf("param::readFile: file %s does not exist\n", fn);
	}
	throw(-1);
    }

    while (!file.eof()) {
        string line;
        getline(file, line);
	add(line);
    }

    file.close();
}


// Inquire whether a variable exists
bool param::exist(const string &name)
{
    ParamMap::iterator it = _map.find(name);
    return (it != _map.end());
}


// Add parameters from an input line
void param::add(const string &input_line)
{
    string line = input_line;
    vector<string> tokens;
    string name;

    // Delete the anything behind "#" if it exists
    size_t pos = line.find_first_of("#");
    if (pos != string::npos) line.erase(pos);

    Tokenize(line, tokens, " \t\n=,;");
    if (tokens.size() == 0) return;	// no values

    name = tokens[0];
    tokens.erase(tokens.begin());
    _map[name] = tokens;
}


// Get an integer parameter
int param::getIntValue(const string &name, int pos)
{
    ParamMap::iterator it = _map.find(name);

    assert(it != _map.end());
    assert(pos >= 0 && pos < it->second.size());

    return string_to_val<int>(it->second[pos]);
}


// Get a double parameter
double param::getDoubleValue(const string &name, int pos)
{
    ParamMap::iterator it = _map.find(name);

    assert(it != _map.end());
    assert(pos >= 0 && pos < it->second.size());

    return string_to_val<double>(it->second[pos]);
}


// Get a string parameter
string param::getStringValue(const string &name)
{
    ParamMap::iterator it = _map.find(name);

    return it->second[0];
}


// Tokenize an input string
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters)
{
    string::size_type p0 = str.find_first_not_of(delimiters, 0);
    string::size_type p1 = str.find_first_of(delimiters, p0);

    while (string::npos != p0 || string::npos != p1) {
        tokens.push_back(str.substr(p0, p1-p0));

        p0 = str.find_first_not_of(delimiters, p1);
        p1 = str.find_first_of(delimiters, p0);
    }
}


// String to value
template<class T>
T string_to_val(const string & str)
{
    T val;

    istringstream iss(str);
    iss >> val;
    return val;
}
