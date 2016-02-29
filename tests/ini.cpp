#include <iostream>

#include "ini.h"

using namespace std;

bool TestIni::testTrim() {
  if (Ini::trim("") != "") {
    cerr << "Ini::trim(\"\"): should empty string" << endl;
    return false;
  }

  if (Ini::trim("  name") != "name") {
    cerr << "Ini::trim(\"  name\"): should remove left spaces (\"";
    cerr << Ini::trim("  name") << "\"" << endl;
    return false;
  }

  if (Ini::trim("name  ") != "name") {
    cerr << "Ini::trim(\"name  \"): should remove right spaces (\"";
    cerr << Ini::trim("name  ") << "\"" << endl;
    return false;
  }

  if (Ini::trim("  name  ") != "name") {
    cerr << "Ini::trim(\"  name  \"): should remove both right and left spaces (\"";
    cerr << Ini::trim("  name  ") << "\"" << endl;
    return false;
  }

  return true;
}

bool TestIni::testParseLine() {
  bool valid;
  string name, value;

  // valid syntax
  valid = Ini::parseLine("#This is a comment", name, value);
  if (valid == true) {
    cerr << "Ini::parseLine: comment line should not be considered as valid" << endl;
    return false;
  }

  valid = Ini::parseLine("name1=value1", name, value);
  if (valid == false) {
    cerr << "Ini::parseLine(\"name1=value1\"): should be valid." << endl;
    return false;
  } else if (name != "name1") {
    cerr << "Ini::parseLine(\"name1=value1\"): name == " << "\"" << name << "\"" << endl;
    return false;
  } else if (value != "value1") {
    cerr << "Ini::parseLine(\"name1=value1\"): value == " << "\"" << value << "\"" << endl;
    return false;
  }

  valid = Ini::parseLine("name2 = value2", name, value);
  if (valid == false) {
    cerr << "Ini::parseLine(\"name2 = value2\"): should be valid." << endl;
    return false;
  } else if (name != "name2") {
    cerr << "Ini::parseLine(\"name2 = value2\"): name == " << "\"" << name << "\"" << endl;
    return false;
  } else if (value != "value2") {
    cerr << "Ini::parseLine(\"name2 = value2\"): value == " << "\"" << value << "\"" << endl;
    return false;
  }

  valid = Ini::parseLine("name3 = value3 #comment", name, value);
  if (valid == false) {
    cerr << "Ini::parseLine(\"name3 = value3 #comment\"): should be valid." << endl;
    return false;
  } else if (name != "name3") {
    cerr << "Ini::parseLine(\"name3 = value3 #comment\"): name == " << "\"" << name << "\"" << endl;
    return false;
  } else if (value != "value3") {
    cerr << "Ini::parseLine(\"name3 = value3 #comment\"): value == " << "\"" << value << "\"" << endl;
    return false;
  }

  valid = Ini::parseLine("name4 =", name, value);
  if (valid == false) {
    cerr << "Ini::parseLine(\"name4 =\"): should be valid." << endl;
    return false;
  } else if (name != "name4") {
    cerr << "Ini::parseLine(\"name4 =\"): name == " << "\"" << name << "\"" << endl;
    return false;
  } else if (value != "") {
    cerr << "Ini::parseLine(\"name4 =\"): value == " << "\"" << value << "\"" << endl;
    return false;
  }

  valid = Ini::parseLine("name5 = =3=", name, value);
  if (valid == false) {
    cerr << "Ini::parseLine(\"name5 = =3=\"): should be valid." << endl;
    return false;
  } else if (name != "name5") {
    cerr << "Ini::parseLine(\"name5 = =3=\"): name == " << "\"" << name << "\"" << endl;
    return false;
  } else if (value != "=3=") {
    cerr << "Ini::parseLine(\"name5 = =3=\"): value == " << "\"" << value << "\"" << endl;
    return false;
  }


  // not valid syntax
  valid = Ini::parseLine("name6 # comment = value6", name, value);
  if (valid == true) {
    cerr << "Ini::parseLine(\"name6 # comment = value6\"): should not be valid." << endl;
    return false;
  }

  valid = Ini::parseLine(" =", name, value);
  if (valid == true) {
    cerr << "Ini::parseLine(\" =\"): should not be valid." << endl;
    return false;
  }

  return true;
}

int main(int argc, char** argv) {
  if (TestIni::testTrim() == false) {
    return 1;
  }

  if (TestIni::testParseLine() == false) {
    return 1;
  }

  TestIni ini;
  ini.loadFromString("test1=33\ntest2=44\n");
  cout << ini.getParameter("test1") << endl;
  cout << ini.getParameter("test2") << endl;

  return 0;
}
