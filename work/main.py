import numpy, scipy
from scipy.linalg import block_diag

import pyscf
from pyscf import gto, scf, mcscf

from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes

tmpdir = pyscf.lib.param.TMPDIR

bond_dims = [250] * 4 + [500] * 4
noises = [1e-4] * 4 + [1e-5] * 4 + [0]
thrds = [1e-10] * 8

def main(dista):
	assert dista > 0.0

	mol = gto.Mole()
	mol.atom = [['O',(0.0, 0.0, 0.0)], ['O',(0.0, 0.0, dista)],]
	mol.spin = 2
	mol.basis = 'cc-pvdz'
	mol.symmetry = 'C2v'
	mol.build()

	# RHF case (for spin-adapted / non-spin-adapted DMRG)
	mf = scf.RHF(mol).run()
	ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_rhf_integrals(
		mf, ncore=0, ncas=None, g2e_symm=8
		)
	print("ncas = %s, n_elec = %s, spin = %s, ecore = %s" % (ncas, n_elec, spin, ecore))
	print("orb_sym = %s" % orb_sym)

	driver = DMRGDriver(scratch=tmpdir, symm_type=SymmetryTypes.SU2, n_threads=16)
	driver.initialize_system(n_sites=ncas, n_elec=n_elec, spin=spin, orb_sym=orb_sym)

	mpo = driver.get_qc_mpo(h1e=h1e, g2e=g2e, ecore=ecore, iprint=1)
	ket = driver.get_random_mps(tag="GS", bond_dim=250, nroots=1)
	energy = driver.dmrg(
		mpo, ket, n_sweeps=20, bond_dims=bond_dims, 
		noises=noises, thrds=thrds, iprint=1
		)
	
	print('Final result: dista = %12.8f, energy = %12.8f' % (dista, energy))

if __name__ == "__main__":
	main(1.2)
	main(3.0)