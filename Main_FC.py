#############################################
# This code is rewritten from Tom's Fortran #
# code and only works for Sz=0 cases.	    #
#############################################
import pickle
import time
import os
import pyscf
from scipy.linalg import block_diag
from pyscf import gto, scf, mcscf, lib, dmrgscf, fci
dmrgscf.settings.BLOCKEXE = os.popen("which block2main").read().strip()
dmrgscf.settings.MPIPREFIX = ''

dista = 1.2
mol = gto.Mole()
#mol.build(
#verbose = 4,
#atom = [['Cr',(  0.000000,  0.000000, 0.000000)],
#	    ['Cr',(  0.000000,  0.000000, dista)], ],
#		basis = 'cc-pvdz-dk',
#		symmetry = 1,
#		)
mol.atom = [['O',(0.0, 0.0, 0.0)], ['O',(0.0, 0.0, dista)],]
mol.spin = 2
mol.basis = 'cc-pvdz'
mol.symmetry = 'C2v'
mol.build()
NFC = 0
	
## RHF
RHF = scf.RHF(mol)#.sfx2c1e()
RHF.diis_space = 20
if (os.path.exists("RHF.p")):
	dm = pickle.load(open( "RHF.p", "rb" ))
	print("load RHF DM")
	RHF.kernel(dm)
else:
	RHF.kernel()
# check RHF stab
#RHFstab = RHF.stability(external=True)
# save and print results
print("E(RHF)=",RHF.e_tot)

#fci
#mci = fci.FCI(mol, RHF.mo_coeff)
#mci.wfnsym = 'B2'
#mci = fci.addons.fix_spin_(mci, ss=0)
#e, civec = mci.kernel()
#print("E(FCI)=", e)

# CASSCF-fci
#cas = mcscf.CASSCF(RHF,12,12)
#cas.frozen = NFC
#cas.fcisolver.threads = 1
#cas.fcisolver.conv_tol = 1e-9
#cas.kernel()

# CASSCF-dmrg
mc = dmrgscf.DMRGSCF(RHF, mol.nao, mol.nelec, tol=1E-6)
mc.frozen = NFC
mc.fcisolver.runtimeDir = lib.param.TMPDIR
mc.fcisolver.scratchDirectory = lib.param.TMPDIR
mc.fcisolver.threads = int(os.environ.get("OMP_NUM_THREADS", 4))
mol.max_memory=30000
mc.fcisolver.memory = int(mol.max_memory / 1000) # mem in GB
mc.kernel()
print(mc.e_tot)

# CCSD(T)
#RCC = cc.CCSD(RHF, frozen = NFC)
#RCC.kernel()
#print("E(CCSD)=", RCC.e_tot)
#et = RCC.ccsd_t()
#print("E(CCSD-T)=", et + RCC.e_tot)

