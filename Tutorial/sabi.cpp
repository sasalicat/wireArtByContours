#progma once
#include "stdafx.h"
#include <ogldev_util.h>

bool ReadFile(const char* fileName, string& outFile){
	ifstream f(fileName);
	bool ret = false;
	if (f.is_open()) {
		string line;
		while (getline(f, line)) {
			outFile.append(line);
			outFile.append("\n");
		}
		f.close();
		ret = true;
	}

	return ret;
}