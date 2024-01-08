#ifndef JDFTX_ELECTRONIC_MIXGRADIENT_H
#define JDFTX_ELECTRONIC_MIXGRADIENT_H

#include <electronic/ElecInfo.h>

class MixGradient : public std::vector<std::vector< std::vector<double>>> {
	
public:
	void print(Everything& e, FILE* fp) const;

};

#endif
