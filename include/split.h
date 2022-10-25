#ifndef __SPLIT_H
#define __SPLIT_H

#include <vector>
#include <string>
using namespace std;

vector<string> split(string data, string token) {
    vector<string> output;
    size_t pos = string::npos; // size_t to avoid improbable overflow
    do {
        pos = data.find(token);
        output.push_back(data.substr(0, pos));
        if (string::npos != pos)
            data = data.substr(pos + token.size());
    } while (string::npos != pos);
    return output;
}

#endif