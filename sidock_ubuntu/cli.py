from __future__ import annotations
import argparse
from pathlib import Path
from .pipeline import run_all


def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(prog='sidock-ubuntu',description='Headless SIDock AutoDock4 automation for Ubuntu')
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('run',help='Run receptor + SDF through the complete docking workflow')
    r.add_argument('receptor',type=Path)
    r.add_argument('ligands',type=Path)
    r.add_argument('--project',type=Path,default=Path('sidock_project'))
    r.add_argument('--padding',type=float,default=8.0)
    r.add_argument('--ga-runs',type=int,default=10)
    r.add_argument('--jobs',type=int,default=1)
    r.add_argument('--top-n',type=int,default=20)
    return p


def main() -> None:
    args=parser().parse_args()
    if args.command=='run':
        run_all(args.project,args.receptor,args.ligands,padding=args.padding,ga_runs=args.ga_runs,parallel_jobs=args.jobs,top_n=args.top_n)


if __name__=='__main__':
    main()
