reinitialize
load benchmark_v031_diagnosis/run_native/analysis/validation/native_ligand.pdbqt, dasatinib
hide everything
show sticks, dasatinib
color gray80, dasatinib
set label_color, white
set dash_width, 3
distance T1, dasatinib and id 1, dasatinib and id 10
color cyan, T1
label dasatinib and id 10, "T1"
distance T2, dasatinib and id 10, dasatinib and id 12
color cyan, T2
label dasatinib and id 12, "T2"
distance T3, dasatinib and id 14, dasatinib and id 18
color cyan, T3
label dasatinib and id 18, "T3"
distance T4, dasatinib and id 18, dasatinib and id 19
color cyan, T4
label dasatinib and id 19, "T4"
distance T5, dasatinib and id 19, dasatinib and id 21
color cyan, T5
label dasatinib and id 21, "T5"
distance T6, dasatinib and id 18, dasatinib and id 29
color red, T6
label dasatinib and id 29, "T6"
distance T7, dasatinib and id 3, dasatinib and id 31
color cyan, T7
label dasatinib and id 31, "T7"
distance T8, dasatinib and id 36, dasatinib and id 37
color cyan, T8
label dasatinib and id 37, "T8"
distance T9, dasatinib and id 37, dasatinib and id 38
color cyan, T9
label dasatinib and id 38, "T9"
distance T10, dasatinib and id 38, dasatinib and id 39
color yellow, T10
label dasatinib and id 39, "T10"
orient dasatinib
ray 1200,900
png benchmark_v031_diagnosis/analysis_validation/torsion_tree_audit.png,1200,900,1
quit