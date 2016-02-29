#ifndef __BBN_INI_H__
#define __BBN_INI_H__

#include <istream>
#include <sstream>
#include <map>
#include <string>

class Ini {
protected:
  std::map<std::string, std::string> data;
  static bool parseLine(const std::string &line, std::string &name, std::string &value);
  static std::string trim(const std::string &str);
public:
  void loadFromIStream(std::istream &is);
  void loadFromFile(const std::string &filePath);
  void loadFromString(const std::string &str);

  std::string getParameter(const std::string &name);
};

#endif // __BBN_INI_H__
