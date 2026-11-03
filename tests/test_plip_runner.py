from pathlib import Path
from types import SimpleNamespace
from sidock_ubuntu.plip_runner import prepare_protonated_complex, run_plip, normalize_plip_interactions
from sidock_ubuntu.errors import SIDockError

def test_obabel_command_is_recovered(tmp_path):
    calls=[]
    def run(cmd,**kwargs): calls.append(cmd); Path(cmd[3]).write_text('PDB'); return SimpleNamespace(returncode=0)
    out=prepare_protonated_complex(tmp_path/'c.pdb',tmp_path/'p.pdb',runner=run)
    assert calls[0]==['obabel',str(tmp_path/'c.pdb'),'-O',str(tmp_path/'p.pdb'),'-h']

def test_plip_command_includes_nohydro(tmp_path):
    calls=[]; out=tmp_path/'o'; out.mkdir()
    (out/'x_report.xml').write_text('<report/>')
    def run(cmd,**kwargs): calls.append(cmd); return SimpleNamespace(returncode=0)
    run_plip(tmp_path/'p.pdb',out,runner=run)
    assert '--nohydro' in calls[0] and calls[0][0].endswith('/plip')

def test_plip_failure_maps_to_e701(tmp_path):
    def run(cmd,**kwargs): return SimpleNamespace(returncode=1)
    try: prepare_protonated_complex(tmp_path/'c.pdb',tmp_path/'p.pdb',runner=run)
    except SIDockError as e: assert e.code_name=='PLIP_FAILED'
    else: raise AssertionError('expected PLIP failure')

def test_normalized_schema_source_is_plip(tmp_path):
    xml='<report><bindingsite><hydrogen_bonds><hydrogen_bond><reschain>D</reschain><restype>GLN</restype><resnr>210</resnr><dist>3.0</dist><angle>120</angle></hydrogen_bond></hydrogen_bonds></bindingsite></report>'
    report=tmp_path/'x_report.xml'; report.write_text(xml); out=normalize_plip_interactions([report],tmp_path/'i.csv',compound_id='c',pose_id='p')
    text=out.read_text(); assert 'PLIP' in text and 'GLN' in text and 'hydrogen_bond' in text

def test_plip_xml_recovers_type_specific_distances_and_real_atom_indices(tmp_path):
    import csv
    xml = """<report><bindingsite>
      <hydrogen_bonds><hydrogen_bond>
      <reschain>A</reschain><restype>ASP</restype><resnr>215</resnr>
      <dist_d-a>2.84</dist_d-a><dist_h-a>1.95</dist_h-a>
      <don_angle>162.0</don_angle><protisdon>False</protisdon>
      <donoridx>110</donoridx><acceptoridx>220</acceptoridx>
      </hydrogen_bond></hydrogen_bonds>
      <hydrophobic_interactions><hydrophobic_interaction>
      <reschain>A</reschain><restype>TRP</restype><resnr>150</resnr>
      <dist>3.72</dist><protcarbonidx>12</protcarbonidx><ligcarbonidx>24</ligcarbonidx>
      </hydrophobic_interaction></hydrophobic_interactions>
      <pistacking><pi_stack><reschain>A</reschain><restype>PHE</restype>
      <resnr>220</resnr><centdist>4.67</centdist><angle>19.0</angle>
      </pi_stack></pistacking></bindingsite></report>"""
    path = tmp_path / "plip_report.xml"
    path.write_text(xml)
    dest = normalize_plip_interactions([path], tmp_path / "contact.csv", compound_id="7i")
    with dest.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    hbond, hydrophob, pi = rows
    assert hbond["distance_A"] == "2.84"
    assert hbond["angle_deg"] == "162.0"
    assert hbond["protein_atom"] == "220"
    assert hbond["ligand_atom"] == "110"
    assert hydrophob["distance_A"] == "3.72"
    assert hydrophob["protein_atom"] == "12"
    assert pi["distance_A"] == "4.67"
    assert pi["angle_deg"] == "19.0"
