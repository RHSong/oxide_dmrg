import os, numpy, scipy
from scipy.linalg import block_diag

import pyscf
from pyscf import gto, scf, mcscf

from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes

tmpdir = pyscf.lib.param.TMPDIR

# bond_dims = [250] * 4 + [500] * 4 + [1000] * 4
# noises = [1e-4] * 4 + [1e-5] * 4 + [1e-6] * 4
# thrds = [1e-10] * 12
bond_dims = [250] * 4 + [500] * 4 + [750] * 4 + [1000] * 4
noises = [1e-4] * 4 + [1e-5] * 12 + [0]
thrds = [1e-10] * 30
n_threads = int(os.environ.get("OMP_NUM_THREADS", "20"))

def main(dista):
    print("\n\n=====================================================")
    print("Running calculation for dista = %12.8f" % dista)
    print("bond_dims = %s" % bond_dims)
    print("noises = %s" % noises)
    print("thrds = %s" % thrds)
    print("n_threads = %s" % n_threads)

    assert dista > 0.0

    mol = gto.Mole()
    mol.atom = [['O',(0.0, 0.0, 0.0)], ['O',(0.0, 0.0, dista)],]
    mol.spin = 2
    mol.basis = 'ccpvdz'
    mol.symmetry = None # 'd2h'
    mol.build()

    # RHF case (for spin-adapted / non-spin-adapted DMRG)
    mf = scf.RHF(mol)
    mf.verbose = 4
    mf.max_cycle = -1
    mf.kernel()

    from pyscf import lo
    mf.mo_coeff = lo.orth_ao(mol, 'meta_lowdin')
    
    print("mf.mo_energy =\n %s" % mf.mo_energy)
    from pyscf.tools.dump_mat import dump_mo
    dump_mo(mol, mf.mo_coeff, digits=6)

    ncas, n_elec, spin, ecore, h1e, g2e, orb_sym = itg.get_rhf_integrals(
        mf, ncore=0, ncas=None, g2e_symm=8
        )
    print("ncas = %s, n_elec = %s, spin = %s, ecore = %s" % (ncas, n_elec, spin, ecore))
    print("orb_sym = %s" % orb_sym)

    driver = DMRGDriver(scratch=tmpdir, symm_type=SymmetryTypes.SU2, n_threads=4, stack_mem=4 << 30)
    driver.initialize_system(n_sites=ncas, n_elec=n_elec, spin=spin, orb_sym=orb_sym)

    mpo = driver.get_qc_mpo(h1e=h1e, g2e=g2e, ecore=ecore, iprint=1)
    ket = driver.get_random_mps(tag="GS", bond_dim=bond_dims[0], nroots=1)
    energy = driver.dmrg(
        mpo, ket, n_sweeps=30, bond_dims=bond_dims,
        noises=noises, thrds=thrds, iprint=1,
        twosite_to_onesite=20
    )

    print('Final result: dista = %12.8f, energy = %12.8f' % (dista, energy))

if __name__ == "__main__":
    main(1.2)
    main(3.0)
