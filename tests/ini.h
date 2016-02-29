#ifndef __BBN_TEST_INI_H__
#define __BBN_TEST_INI_H__

#include "bbn/ini.h"

// to access protected methods
class TestIni : public Ini {
public:
  static bool testParseLine();
  static bool testTrim();
};

#endif // __BBN_TEST_INI_H__
