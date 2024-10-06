# https://qiita.com/oyster-tempura/items/2102996ab6b87e54e5a7

import matplotlib.pyplot as plt

from ase.visualize.plot import plot_atoms
from gpaw import GPAW, PW
from ase.build import molecule
from ase.visualize import view

gpaw = GPAW(mode=PW(300), xc='PBE')

# H2
atoms_h2 = molecule('H2')
atoms_h2.center(vacuum=3.)
# view(atoms_h2)
atoms_h2.set_calculator(gpaw)
e_h2 = atoms_h2.get_potential_energy()  # takes several seconds

# H
atoms_h = molecule('H')
atoms_h.center(vacuum=3.)
atoms_h.set_calculator(GPAW(mode=PW(300), xc='PBE', hund=True))
e_h = atoms_h.get_potential_energy()

# 結果の表示
delta_e = 2 * e_h - e_h2
print(f'=== Atomization Energy: {delta_e: 5.2f} eV ===')

fig, ax = plt.subplots(1, 2, figsize=(8, 6))

title = f'Calculation Result for {atoms_h2.get_chemical_formula()}\n' + \
        f'Total Energy : {atoms_h2.get_potential_energy(): 5.2f} eV'
ax[0].set_title(title)
ax[0].axis('off')
plot_atoms(atoms_h2, ax=ax[0], rotation='90x')

title = f'Calculation Result for {atoms_h.get_chemical_formula()}\n' + \
        f'Total Energy : {atoms_h.get_potential_energy(): 5.2f} eV'
ax[1].set_title(title)
ax[1].axis('off')
plot_atoms(atoms_h, ax=ax[1], rotation='90x')

plt.show()
