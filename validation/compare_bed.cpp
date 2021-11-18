#include <iostream>
#include <map>
#include <stdio.h>
#include <fstream>
#include<regex> 

using namespace std;
const string WHITESPACE = " \n\r\t\f\v";

std::string lefttrim(const std::string &s) {
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == std::string::npos) ? "" : s.substr(start);
}
 
std::string righttrim(const std::string &s) {
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}
 
std::string trim(const std::string &s) {
    return righttrim(lefttrim(s));
}

int main (int argc, char **argv)
{
    
    ifstream INFILE1(argv[1]);
    ifstream INFILE2(argv[2]); 
    string line, chrom, strand;
    int start, end, count; 
    pair<map<tuple<string, int, int, string>, int>::iterator,bool> ret;
    map<tuple<string, int, int, string>, int>::iterator it;

    // use chr-start-end as the key, and +/- as the value
    map<tuple<string, int, int, string>, int> mymap;

    // read the first file and insert the elements to the map
    if (INFILE1.is_open()) {    
        while (getline(INFILE1, line)) {
            trim(line); 
            if (line.empty()) continue; 

            vector<string> fields; 
            stringstream ss(line);
            string token; 
            while(ss >> token) { 
                fields.push_back(token); 
            }

            chrom = fields[0];
            start = std::stoi(fields[1]);	
            end = std::stoi(fields[2]);
            strand = fields[5];

            count = 1;
            it = mymap.find(make_tuple(chrom, start,end, strand));

            // value found, increase the count
            if (it != mymap.end()) {
                count = it->second + 1;
            }

            mymap.insert(pair<tuple<string, int, int, string>, int>(make_tuple(chrom, start,end,strand), count));

        }
    }

    // read the second file and compare its elements with the previous one
    
    if (INFILE2.is_open()) {    
        while (getline(INFILE2, line)) {
            trim(line); 
            if (line.empty()) continue; 

            vector<string> fields; 
            stringstream ss(line);
            string token; 
            while(ss >> token) { 
                fields.push_back(token); 
            }

            chrom = fields[0];
            start = std::stoi(fields[1]);	
            end = std::stoi(fields[2]);
            strand = fields[5];

            it = mymap.find(make_tuple(chrom, start,end, strand));

            // value found, no difference
            if (it != mymap.end()) {
                // decrement the count if it is greater than 0
                if(it->second > 0){
                    it->second = it->second - 1;
                    cout << "Elements matching " << line << endl;
                }
                else{ // elements are not matching any more
                    cout << "Element only exists in file 2 " << line << endl;
                }
            }
            // value does not found, report as only exists in file 2
            else{
                cout << "Element only exists in file 2 " << line << endl;
            }
        }
    }

    // print the elements that exist only in file 1
    for (it=mymap.begin(); it!=mymap.end(); ++it){
        if(it->second > 0)
            cout << "only exists in file1: " << get<0>(it->first) << "\t" << get<1>(it->first) << "\t"  << get<2>(it->first) << "\t" << get<3>(it->first) << '\n';
    }
    return 0;
}