/*-------------------------------------------------------------------
Copyright 2023 Brandon Li

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <perturb/SpringConstant.h>
#include <electronic/IonInfo.h>
#include <electronic/Energies.h>
#include <electronic/ElecVars.h>
#include <perturb/PerturbationSolver.h>
#include <electronic/Everything.h>

SpringConstant::SpringConstant(Everything& e) : e(e), eVars(e.eVars), eInfo(e.eInfo), iInfo(e.iInfo), pInfo(e.pertInfo), ps(e) {}

double SpringConstant::computeMatrixElement(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB)
{
	//TODO disallow fillings hsub
	assert(pInfo.commensurate);
	pInfo.datom = modeB;
	ps.updateNonlocalDerivs();
	ps.calcdGradTau();	
	return dot(pInfo.dGradTau.C, pInfo.dC, &pInfo, &eInfo) + dsqEpair(modeA, modeB) + dsqEnl(modeA, modeB) + dsqEloc(modeA, modeB);
}

double SpringConstant::dsqEpair(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB)
{
	static StopWatch watch("dsqEpair"); watch.start();
	double h = 1e-3;
	double k = 1e-3;
	
	assert(h > 0 && k > 0);
	
	Energies Eminusminus, Eminusplus, Eplusminus, Eplusplus;
	getPerturbedEnergy(Eminusminus, modeA, modeB, -h/2.0, -k/2.0);
	getPerturbedEnergy(Eminusplus, modeA, modeB, -h/2.0, k/2.0);
	getPerturbedEnergy(Eplusminus, modeA, modeB, h/2.0, -k/2.0);
	getPerturbedEnergy(Eplusplus, modeA, modeB, h/2.0, k/2.0);
	
	double dsqEpair = ((Eplusplus.E["Eewald"]+Eplusplus.E["EvdW"])
		- (Eplusminus.E["Eewald"]+Eplusminus.E["EvdW"])
		- (Eminusplus.E["Eewald"]+Eminusplus.E["EvdW"])
		+ (Eminusminus.E["Eewald"]+Eminusminus.E["EvdW"])) / (h*k);

	watch.stop();
	return dsqEpair;
}

double SpringConstant::dsqEnl(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB)
{
	static StopWatch watch("dsqEnl"); watch.start();
	
	pInfo.datom = modeA;
	ps.updateNonlocalDerivs();
	
	pInfo.datom = modeB;
	ps.updateNonlocalDerivs();
	
	double Enl = 0;
	if (modeA->sameAtom(modeB)) {
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			matrix dSqVdagC = -D(modeA->dVatom[q], modeB->mode.dirCartesian)^eVars.C[q];
			Enl += real(trace(dSqVdagC*eVars.F[q]*dagger(modeA->Vatom[q]^eVars.C[q])*e.iInfo.species[modeA->mode.sp]->MnlAll))*eInfo.qnums[q].weight;
			Enl += real(trace((modeA->Vatom[q]^eVars.C[q])*eVars.F[q]*dagger(dSqVdagC)*e.iInfo.species[modeA->mode.sp]->MnlAll))*eInfo.qnums[q].weight;
			Enl += real(trace(modeA->dVdagCatom[q]*eVars.F[q]*dagger(modeB->dVdagCatom[q])*e.iInfo.species[modeA->mode.sp]->MnlAll))*eInfo.qnums[q].weight;
			Enl += real(trace(modeB->dVdagCatom[q]*eVars.F[q]*dagger(modeA->dVdagCatom[q])*e.iInfo.species[modeA->mode.sp]->MnlAll))*eInfo.qnums[q].weight;
			//Is expression hermitian
		}
		mpiWorld->allReduce(Enl, MPIUtil::ReduceSum);
	}
	watch.stop();
	return Enl;
}

double SpringConstant::dsqEloc(std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB)
{
	static StopWatch watch("dsqEloc"); watch.start();
	double dsqE = 0;
	if (modeA->Vlocps && modeA->sameAtom(modeB)) {
		ScalarFieldTilde nTilde = J(eVars.get_nTot());
		ScalarFieldTilde dsqVlocps = D(D(modeA->Vlocps, modeA->mode.dirCartesian), modeB->mode.dirCartesian);
		dsqE += dot(nTilde, O(dsqVlocps));
	}
	watch.stop();
	return dsqE;
}

void SpringConstant::getPerturbedEnergy(Energies& ener, std::shared_ptr<AtomPerturbation> modeA, std::shared_ptr<AtomPerturbation> modeB, double deltaA, double deltaB)
{
	static StopWatch watch("getPerturbedEnergy"); watch.start();
	auto spA = iInfo.species[modeA->mode.sp];
	auto spB = iInfo.species[modeB->mode.sp];
	
	vector3<> posA_unperturbed = spA->atpos[modeA->mode.at];
	spA->atpos[modeA->mode.at] = posA_unperturbed + deltaA*modeA->mode.dirLattice;
	mpiWorld->bcastData(spA->atpos);
	
	vector3<> posB_unperturbed = spB->atpos[modeB->mode.at];
	spB->atpos[modeB->mode.at] = posB_unperturbed + deltaB*modeB->mode.dirLattice;
	mpiWorld->bcastData(spB->atpos);
	
	spA->sync_atpos();
	spB->sync_atpos();
	
	e.iInfo.pairPotentialsAndGrad(&ener);
	
	spB->atpos[modeB->mode.at] = posB_unperturbed;
	mpiWorld->bcastData(spB->atpos);
	spA->atpos[modeA->mode.at] = posA_unperturbed;
	mpiWorld->bcastData(spA->atpos);
	
	spB->sync_atpos();
	spA->sync_atpos();
	watch.stop();
}

IonicGradient SpringConstant::getPhononMatrixColumn(std::shared_ptr<AtomPerturbation> modeA, double h)
{
	pInfo.datom = modeA;
	ps.solvePerturbation();
	
	if (h > 0)
	{
		IonicGradient Fminus, Fplus, dF;
		init(Ctmp, eInfo.nStates, eInfo.nBands, &e.basis[0], &eInfo);
		
		auto spA = iInfo.species[modeA->mode.sp];
		vector3<> posA_unperturbed = spA->atpos[modeA->mode.at];
		
		spA->atpos[modeA->mode.at] = posA_unperturbed + h*modeA->mode.dirLattice;
		mpiWorld->bcastData (spA->atpos);
		spA->sync_atpos();
		
		iInfo.update(e.ener);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++)
		{
			Ctmp[q] = eVars.C[q];
			eVars.C[q] += pInfo.dC[q]*h;
			if (modeA->isUltrasoft(iInfo))
				eVars.C[q] += modeA->dCatom[q]*h;
			eVars.VdagC[q].clear();
			e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
		}
		eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
		
		e.iInfo.ionicEnergyAndGrad();
		Fplus = -e.gInfo.invRT * e.iInfo.forces;
		
		spA->atpos[modeA->mode.at] = posA_unperturbed - h*modeA->mode.dirLattice;
		mpiWorld->bcastData (spA->atpos);
		spA->sync_atpos();
		
		iInfo.update(e.ener);
		for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
			eVars.C[q] = Ctmp[q];
			eVars.C[q] -= pInfo.dC[q]*h;
			if (modeA->isUltrasoft(iInfo))
				eVars.C[q] -= modeA->dCatom[q]*h;
			eVars.VdagC[q].clear();
			e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
		}
		eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
		
		e.iInfo.ionicEnergyAndGrad();
		Fminus = -e.gInfo.invRT * e.iInfo.forces;
		
		spA->atpos[modeA->mode.at] = posA_unperturbed;
		mpiWorld->bcastData ( spA->atpos );
		spA->sync_atpos();
		
		iInfo.update(e.ener);
        for (int q=eInfo.qStart; q<eInfo.qStop; q++) {
			eVars.C[q] = Ctmp[q];
			eVars.VdagC[q].clear();
			e.iInfo.project(eVars.C[q], eVars.VdagC[q]);
		}
		eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
		
		dF = (Fplus - Fminus)*(1/(2*h));
		return dF;
	}
	else
	{
		if (pInfo.densityAugRequired(e))
			die("Analytic spring constant calculation not available for ultrasoft potentials.\n");
		
		IonicGradient dF;
		dF.init (e.iInfo);
		for (unsigned s = 0; s < iInfo.species.size(); s++) {
			auto sp = iInfo.species[s];
			for (unsigned at = 0; at < sp->atpos.size(); at++) {
				for (int iDir = 0; iDir < 3; iDir++) {
					std::shared_ptr<AtomPerturbation> modeB = std::make_shared<AtomPerturbation>(s, at, iDir, e);
					pInfo.datom = modeB;
					dF[s][at][iDir] = computeMatrixElement (modeA, modeB);
				}
			}
		}
		pInfo.datom = modeA;
		ps.calcdGradTau();
		return dF;
	}
}
