from pathlib import Path
from sidock_ubuntu.intake import choose_reference_ligand


def test_choose_largest_nonwater_ligand(tmp_path: Path):
    p=tmp_path/'x.pdb'
    lines=['HETATM    1  O   HOH A   1      0.000   0.000   0.000  1.00 20.00           O']
    for i in range(1,9):
        lines.append(f'HETATM{i+1:5d}  C{i%10:<2} LIG A 100      {i:6.3f}   0.000   0.000  1.00 20.00           C')
    p.write_text('\n'.join(lines)+'\n')
    g=choose_reference_ligand(p)
    assert g is not None and g.resname=='LIG' and g.atom_count==8
