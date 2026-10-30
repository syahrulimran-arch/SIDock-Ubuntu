from sidock_ubuntu.cli import parser


def test_cli_run_parser():
    args=parser().parse_args(['run','protein.pdb','compounds.sdf'])
    assert args.command=='run'
