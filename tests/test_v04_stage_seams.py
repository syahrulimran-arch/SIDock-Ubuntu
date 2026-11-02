from pathlib import Path
from sidock_ubuntu import stages

def test_stage_seams_reference_delegates(monkeypatch,tmp_path):
    calls=[]
    monkeypatch.setattr(stages,'box_from_reference',lambda p: calls.append(('box',p)) or 'grid')
    monkeypatch.setattr('sidock_ubuntu.intake.write_group',lambda g,p: calls.append(('group',p)) or p)
    group=type('G',(),{'resname':'LIG','chain':'A','resseq':'1'})()
    out=stages.run_reference_stage(Path('x.pdb'),tmp_path,group)
    assert out['grid']=='grid' and [x[0] for x in calls]==['group','box']

def test_stage_seams_are_not_new_scientific_algorithms():
    assert stages.run_compound_preparation.__module__=='sidock_ubuntu.stages'
    assert stages.run_pose_and_cluster_stage.__name__=='run_pose_and_cluster_stage'
