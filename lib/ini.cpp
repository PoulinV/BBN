#include "bbn/ini.h"

using namespace std;

string Ini::trim(string const &str) {
    if(str.empty())
        return str;

    size_t firstScan = str.find_first_not_of(' ');
    size_t first     = firstScan == string::npos ? str.length() : firstScan;
    size_t last      = str.find_last_not_of(' ');

    return str.substr(first, last-first+1);
}

bool Ini::parseLine(const string &line, string &name, string &value) {
  string::size_type equalPos, commentPos, endPos;

  // look at for =
  equalPos = line.find_first_of("=");
  if (equalPos == string::npos) {
    return false;
  }

  // find end of line and ignore comment
  endPos = line.length();
  commentPos = line.find("#");
  if (commentPos != string::npos) {
    endPos = commentPos-1;
  }

  // comment should not before equal
  if (commentPos < equalPos) {
    return false;
  }

  name = Ini::trim(line.substr(0, equalPos));
  value = Ini::trim(line.substr(equalPos+1, endPos-equalPos));

  if (name == "") {
    return false;
  }

  return true;
}

void Ini::loadFromIStream(istream &is) {
  string line, name, value;

  while(getline(is, line)) {
    bool valid;

    valid = Ini::parseLine(line, name, value);
    if (valid != false) {
      data[name] = value;
    }
  }
}

void Ini::loadFromString(const string &str) {
  stringstream ss(str);
  Ini::loadFromIStream(ss);
}

void Ini::loadFromFile(const string &filePath) {

}

string Ini::getParameter(const string &name) {
  return data[name];
}
