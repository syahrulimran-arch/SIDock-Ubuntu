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
    assert '--nohydro' in calls[0] and calls[0][0]=='plip'

def test_plip_failure_maps_to_e701(tmp_path):
    def run(cmd,**kwargs): return SimpleNamespace(returncode=1)
    try: prepare_protonated_complex(tmp_path/'c.pdb',tmp_path/'p.pdb',runner=run)
    except SIDockError as e: assert e.code_name=='PLIP_FAILED'
    else: raise AssertionError('expected PLIP failure')

def test_normalized_schema_source_is_plip(tmp_path):
    xml='<report><bindingsite><hydrogen_bonds><hydrogen_bond><reschain>D</reschain><restype>GLN</restype><resnr>210</resnr><dist>3.0</dist><angle>120</angle></hydrogen_bond></hydrogen_bonds></bindingsite></report>'
    report=tmp_path/'x_report.xml'; report.write_text(xml); out=normalize_plip_interactions([report],tmp_path/'i.csv',compound_id='c',pose_id='p')
    text=out.read_text(); assert 'PLIP' in text and 'GLN' in text and 'hydrogen_bond' in text
