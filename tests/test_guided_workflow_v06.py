"""Regression tests for SIDock's manually approved docking workflow.

No test below invokes AutoGrid4 or AutoDock4.
"""
from pathlib import Path

import pytest

from sidock_ubuntu.guided_workflow import (
    create_guided_project, protein_profile, select_protein,
    install_reviewed_receptor, approve, approved, ready_to_dock,
)
from sidock_ubuntu.manual_docking import root as manual_root, settings, save_settings


def atom(serial, name, residue, chain, number, element, x=29.0, y=71.0, z=74.0):
    return (f"HETATM{serial:5d} {name:>4s} {residue:>3s} {chain}{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}\n")


def sample_4ubp(path: Path) -> None:
    lines = [
        "HEADER    HYDROLASE                               25-FEB-99   4UBP\n",
        "MODRES 4UBP KCX C  220  LYS  LYSINE NZ-CARBOXYLIC ACID\n",
    ]
    for i, name in enumerate(("N","CA","CB","CG","CD","CE","NZ","C","O","CX","OQ1","OQ2"),1):
        element = "N" if name in ("N","NZ") else "O" if name.startswith("O") else "C"
        lines.append(atom(i, name, "KCX", "C", 220, element))
    for i, (name, el) in enumerate((("C1","C"),("C2","C"),("O2","O"),("N","N"),("O","O")),20):
        lines.append(atom(i, name, "HAE", "C", 800, el, 28+i*.02, 72.0, 75.0))
    lines.extend([
        atom(30,"NI","NI","C",798,"NI",29.,72.,74.),
        atom(31,"NI","NI","C",799,"NI",30.,72.,74.),
        atom(32,"O","HOH","C",901,"O",30.,69.,74.),
        atom(33,"C","ACE","A",0,"C",20.,61.,72.),
        atom(34,"O","ACE","A",0,"O",20.,61.,73.),
        atom(35,"CH3","ACE","A",0,"C",20.,62.,72.),
    ])
    path.write_text("".join(lines) + "END\n",encoding="utf-8")


def new_project(tmp_path):
    p = create_guided_project(tmp_path, "urease")
    assert (p/"project_state.json").is_file()
    assert (p/"project.json").is_file()
    assert not (p/"input/compounds.sdf").is_file()
    protein = p/"input/protein.pdb"
    sample_4ubp(protein)
    import json
    config = json.loads((p/"project.json").read_text())
    config["receptor_source"] = str(protein)
    (p/"project.json").write_text(json.dumps(config))
    return p


def test_step1_only_creates_project(tmp_path):
    p = create_guided_project(tmp_path, "test-project")
    assert p.name == "test-project"
    assert (p/"manual_docking/workflow.json").is_file()
    assert not (p/"input/protein.pdb").exists()
    assert not (p/"input/compounds.sdf").exists()


def test_step2_4ubp_distinguishes_hae_kcx_ni_and_cap(tmp_path):
    p = new_project(tmp_path)
    profile = protein_profile(p)
    roles = {item["key"]: item["role"] for item in profile["groups"]}
    assert "HAE:C:800" in {g["key"] for g in profile["ligands"]}
    assert "KCX:C:220" not in {g["key"] for g in profile["ligands"]}
    assert roles["KCX:C:220"] == "MODIFIED_PROTEIN"
    assert roles["NI:C:798"] == roles["NI:C:799"] == "ION"
    assert roles["ACE:A:0"] == "POTENTIAL_PROTEIN_CAP"
    assert "ACE:A:0" in {g["key"] for g in profile["ligands"]}


def test_step2_cleanup_retains_metal_and_modified_protein(tmp_path):
    p = new_project(tmp_path)
    result = select_protein(p, ligand_key="HAE:C:800",
                            remove_groups=[], remove_waters=True, padding=8.)
    assert result["mode"] == "CRYSTAL_LIGAND"
    assert result["npts"][0] >= 20
    clean = (manual_root(p)/"inputs/protein_clean.pdb").read_text()
    assert "HAE C 800" not in clean
    assert "HOH C 901" not in clean
    assert "KCX C 220" in clean
    assert " NI C 798" in clean
    assert " NI C 799" in clean
    assert "ACE A   0" in clean
    crystal = (manual_root(p)/"inputs/crystal_ligand.pdb").read_text()
    assert "HAE C 800" in crystal
    with pytest.raises(RuntimeError, match="Prepare or import"):
        approve(p, 2)


def test_step2_manual_grid_without_reference(tmp_path):
    p = new_project(tmp_path)
    result = select_protein(p, ligand_key=None, remove_groups=["HAE:C:800"],
                            center=(10.,20.,30.), remove_waters=True)
    assert result["mode"] == "MANUAL_CENTER"
    assert result["center"] == [10.,20.,30.]
    assert not (manual_root(p)/"inputs/crystal_ligand.pdb").exists()
    values = settings(p)
    assert (values["center_x"], values["center_y"], values["center_z"]) == (10.,20.,30.)


def test_receptor_approval_invalidated_by_changed_cleanup(tmp_path):
    p = new_project(tmp_path)
    select_protein(p,ligand_key="HAE:C:800",remove_groups=[],remove_waters=True)
    fake = p/"input/valid-metal-receptor.pdbqt"
    fake.write_text(
        "ATOM      1  C   GLY C  40      30.000  70.000  74.000  1.00  0.00    +0.100 C\n"
        "HETATM    2 NI    NI C 798      29.000  72.000  74.000  1.00  0.00    +2.000 NI\n"
        "HETATM    3 NI    NI C 799      30.000  72.000  74.000  1.00  0.00    +2.000 NI\n"
    )
    install_reviewed_receptor(p,fake)
    approve(p,2)
    assert approved(p,2)
    with pytest.raises(RuntimeError, match="Step 3"):
        ready_to_dock(p)

    select_protein(p,ligand_key=None,remove_groups=["HAE:C:800"],
                   center=(1.,2.,3.),remove_waters=False)
    assert not approved(p,2)
    with pytest.raises(RuntimeError, match="not associated"):
        approve(p,2)


def test_cannot_remove_modified_protein_even_if_requested(tmp_path):
    p = new_project(tmp_path)
    with pytest.raises(ValueError, match="Modified protein"):
        select_protein(p, ligand_key="HAE:C:800",
                       remove_groups=["KCX:C:220"])



def test_ace_can_be_explicit_grid_reference_with_warning_but_remains_in_receptor(tmp_path):
    p = new_project(tmp_path)
    with pytest.raises(ValueError, match="scientific acknowledgment"):
        select_protein(p, ligand_key="ACE:A:0", remove_groups=[])

    result = select_protein(
        p, ligand_key="ACE:A:0", remove_groups=[], remove_waters=True,
        acknowledge_non_ligand_component=True
    )
    assert result["mode"] == "USER_SELECTED_COMPONENT"
    assert result["component_role"] == "POTENTIAL_PROTEIN_CAP"
    assert result["user_selected_without_native_validation"] is True
    assert result["scientific_acknowledgment"] is True
    assert "ACE:A:0" not in result["removed_residues"]
    assert "ACE A   0" in (manual_root(p) / "inputs/protein_clean.pdb").read_text()
    selected = (manual_root(p) / "inputs/crystal_ligand.pdb").read_text()
    assert "ACE A   0" in selected
    assert result["site_warning"]


def test_ace_removal_requires_acknowledgment_and_is_audited(tmp_path):
    p = new_project(tmp_path)
    with pytest.raises(ValueError, match="scientific acknowledgment"):
        select_protein(p, ligand_key=None, remove_groups=["ACE:A:0"],
                       center=(1.,2.,3.))
    choice = select_protein(
        p, ligand_key="ACE:A:0", remove_groups=["ACE:A:0"],
        acknowledge_non_ligand_component=True
    )
    assert "ACE:A:0" in choice["removed_residues"]
    assert "ACE A   0" not in (manual_root(p) / "inputs/protein_clean.pdb").read_text()



def test_guided_step2_html_offers_ace_as_manual_grid_reference(tmp_path):
    from sidock_ubuntu.guided_ui import step2
    p = new_project(tmp_path)
    html = step2(p)
    assert 'value="ACE:A:0"' in html or "value='ACE:A:0'" in html
    assert "acknowledge_non_ligand_component" in html
    assert "POTENTIAL_PROTEIN_CAP" in html
    assert "Grid reference component" in html



def test_failed_receptor_preparation_displays_log_on_step2(tmp_path):
    from sidock_ubuntu.guided_ui import step2
    from sidock_ubuntu.manual_docking import root as manual_root
    import json

    p = new_project(tmp_path)
    select_protein(p, ligand_key="ACE:A:0", remove_groups=[],
                   acknowledge_non_ligand_component=True)
    folder = manual_root(p)
    (folder / "logs").mkdir(parents=True, exist_ok=True)
    (folder / "logs/receptor_meeko.log").write_text(
        "Template generation failed for unknown residues: {'NI'}\n"
        "No template matched for residue_key='A:0'\n"
        "Error: Creation of data structure for receptor failed.\n"
    )
    (folder / "job.json").write_text(json.dumps({
        "status": "FAILED",
        "action": "receptor",
        "error": "Receptor preparation failed"
    }))
    html = step2(p)
    assert "Why Meeko failed" in html
    assert "Unsupported component templates:" in html
    assert "Unmatched residues:" in html
    assert "NI" in html
    assert "A:0" in html
    assert "Download complete Meeko log" in html
    assert (folder / "inputs/protein_clean.pdb").is_file()
    assert not (folder / "receptor.pdbqt").exists()
