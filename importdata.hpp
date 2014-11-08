#ifndef IMPORTDATA_HPP
#define IMPORTDATA_HPP
#include "h5obj.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include<algorithm>

using namespace std;




vector<string> readnames(const char* inputfilename){
  vector<string> rownames;
  
  ifstream file(inputfilename);
  string temp;
  
  while(file){
    getline(file,temp);
    rownames.push_back(temp);
  }
  file.close();
  return(rownames);
}

size_t rowcount(const char* filename)
{
  ifstream file(filename);
  size_t rownum = std::count(std::istreambuf_iterator<char>(file),
			     std::istreambuf_iterator<char>(),'\n');
  file.close();
  return(rownum);
}

size_t colcount(const char* filename)
{
  stringstream line1;
  string line1s;
  ifstream file(filename);
  getline(file,line1s);
  line1<<line1s;
  size_t colnum=count(istreambuf_iterator<char>(line1),
		      istreambuf_iterator<char>(),'\t');
  file.close();
  return(colnum);
}

  
  
#endif
