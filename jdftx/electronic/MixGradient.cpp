#include <electronic/MixGradient.h>
#include <electronic/Everything.h>

void MixGradient::print(Everything& e, FILE* fp) const
{	fprintf(fp, "mix name   | atom | parameters\n");
	for(unsigned sp=0; sp<e.iInfo.species.size(); sp++)
	{	
		auto species = e.iInfo.species[sp];
		if (species->isMixed) {
			for (int i = 0; i < species->atpos.size(); i++) {
				fprintf(fp, "%-10s %3d ", species->name.c_str(), i);
				for (int m = 0; m < species->mixSpecies.size(); m++) {
					fprintf(fp, "%19.12f ", e.iInfo.mixGradient[sp][i][m]);
				}
				fprintf(fp, "\n");
			}
		}
	}
	fprintf(fp, "\n");
}
