"""Scientific plots from recorded coordinates/scores/PLIP, using no GUI backend."""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import shutil
import subprocess

from .openbabel_receptor import _atoms

LABEL = "Exploratory docking; native validation not performed"


def generate_figures(rows, receptor, destination, *, activity_rows=(), structures=None, limit=5):
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import numpy as np

    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    figures, warnings = [], []
    ranked = sorted([row for row in rows if row.get("lowest_docking_energy") is not None],
                    key=lambda row: (row["lowest_docking_energy"], row["compound_id"]))

    def save(fig, name, caption):
        fig.text(.5, .012, LABEL, ha="center", fontsize=8)
        fig.tight_layout(rect=(0, .04, 1, 1))
        fig.savefig(destination / name, dpi=220)
        plt.close(fig)
        figures.append({"path": "figures/" + name, "caption": caption + ". " + LABEL + "."})

    if ranked:
        shown = ranked[:40]
        fig, ax = plt.subplots(figsize=(8, max(3.5, len(shown) * .3)))
        ax.barh([r["compound_id"] for r in shown], [r["lowest_docking_energy"] for r in shown])
        ax.invert_yaxis()
        ax.set_xlabel("AutoDock4 score (kcal/mol)")
        ax.set_title("Docking ranking" + (" — first 40 compounds" if len(ranked) > 40 else ""))
        save(fig, "docking_ranking.png", "Scores read from completed DLG models")

    contacts = Counter()
    for row in ranked:
        for contact in row.get("interactions", []):
            residue = (contact.get("protein_residue_name", "") + contact.get("protein_residue_number", "")
                       + ":" + contact.get("protein_chain", ""))
            contacts[(row["compound_id"], residue)] += 1
    if contacts:
        # Rank by frequency, keeping a readable selection and recording the limit.
        residues = [r for r, _ in Counter({r: sum(n for (_, key), n in contacts.items() if key == r)
                                           for _, r in contacts}).most_common(30)]
        cids = [r["compound_id"] for r in ranked[:40]]
        matrix = np.array([[contacts[(cid, residue)] if next(r for r in ranked if r["compound_id"] == cid)["plip_status"] == "COMPLETE"
                            else np.nan for residue in residues] for cid in cids])
        fig, ax = plt.subplots(figsize=(max(7, len(residues) * .35), max(3.5, len(cids) * .32)))
        cmap = plt.get_cmap("Blues").with_extremes(bad="lightgray")
        chart = ax.imshow(np.ma.masked_invalid(matrix), aspect="auto", cmap=cmap, vmin=0)
        ax.set_xticks(range(len(residues)), residues, rotation=65, ha="right", fontsize=8)
        ax.set_yticks(range(len(cids)), cids, fontsize=8)
        ax.set_title("PLIP residue contacts — up to 40 compounds / 30 residues")
        fig.colorbar(chart, ax=ax, label="Recorded contacts; gray = unavailable")
        save(fig, "residue_contacts.png", "PLIP contact counts for the selected docked ligand")

    paired = [row for row in activity_rows if row.get("pIC50") is not None
              and row.get("energy_kcal_mol") is not None]
    if len(paired) >= 2:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        for row in paired:
            ax.scatter(row["energy_kcal_mol"], row["pIC50"])
            ax.annotate(row["ligand"], (row["energy_kcal_mol"], row["pIC50"]),
                        xytext=(3, 3), textcoords="offset points", fontsize=7)
        ax.set_xlabel("AutoDock4 score (kcal/mol)")
        ax.set_ylabel("Recorded pIC50")
        ax.set_title("Docking score and experimental activity")
        save(fig, "activity_comparison.png", "Only exact positive experimental IC50 values are plotted")

    receptor_atoms = _atoms(Path(receptor), input_pdb=False)
    protein = np.array([a["xyz"] for a in receptor_atoms if a["atom_type"] not in {"H", "HD", "HS"}])
    colors = {"C": "#303030", "A": "#303030", "N": "#3068cf", "NA": "#3068cf",
              "OA": "#d33a38", "O": "#d33a38", "S": "#c9ad00", "SA": "#c9ad00"}
    for row in ranked[:limit]:
        cid = row["compound_id"]
        if structures and cid in structures:
            try:
                from rdkit import Chem
                from rdkit.Chem import Draw, rdDepictor
                mol = Chem.RemoveHs(Chem.Mol(structures[cid]))
                rdDepictor.Compute2DCoords(mol)
                fig, ax = plt.subplots(figsize=(6, 4.5))
                ax.imshow(Draw.MolToImage(mol, size=(900, 600)))
                ax.set_axis_off()
                ax.set_title(cid + " — supplied ligand structure")
                save(fig, cid + "_structure.png", cid + " structure from the hash-linked input SDF")
            except Exception as exc:
                plt.close("all")
                warnings.append(cid + " structure figure unavailable: " + str(exc))
        try:
            atoms = [a for a in _atoms(Path(row["selected_pose"]), input_pdb=False)
                     if a["atom_type"] not in {"H", "HD", "HS"}]
            ligand = np.array([a["xyz"] for a in atoms])
            if not len(ligand):
                raise ValueError("Selected pose has no heavy atoms")
            # Plot exact recorded atom positions. No bond graph or interacting
            # residue is inferred from proximity alone.
            distances = np.min(np.sum((protein[:, None, :] - ligand[None, :, :]) ** 2, axis=2), axis=1)
            pocket = protein[distances <= 25]
            fig = plt.figure(figsize=(7, 5.5))
            ax = fig.add_subplot(projection="3d")
            if len(pocket):
                ax.scatter(*pocket.T, c="#b4bcc7", s=12, alpha=.35, label="Receptor atoms within 5 Å")
            metals = np.array([a["xyz"] for a in receptor_atoms
                               if a["atom_type"] in {"Ni", "NI", "Zn", "Mg", "Ca", "Fe", "Mn", "Cu", "Co"}
                               and np.min(np.sum((ligand - np.array(a["xyz"])) ** 2, axis=1)) <= 25])
            if len(metals):
                ax.scatter(*metals.T, c="#268b57", s=95, label="Metal atoms (saved coordinates)")
            ax.scatter(*ligand.T, c=[colors.get(a["atom_type"], "#aa61b5") for a in atoms],
                       s=38, depthshade=False, label="Docked ligand heavy atoms")
            ax.set(xlabel="X (Å)", ylabel="Y (Å)", zlabel="Z (Å)", title=cid + " — selected pose")
            ax.set_box_aspect((1, 1, 1))
            all_xyz = np.vstack((ligand, pocket)) if len(pocket) else ligand
            center = (all_xyz.max(axis=0) + all_xyz.min(axis=0)) / 2
            radius = max(float(np.max(all_xyz.max(axis=0) - all_xyz.min(axis=0))) * .55, 1.0)
            ax.set_xlim(center[0] - radius, center[0] + radius)
            ax.set_ylim(center[1] - radius, center[1] + radius)
            ax.set_zlim(center[2] - radius, center[2] + radius)
            ax.legend(fontsize=8)
            save(fig, cid + "_pose.png", cid + " recorded coordinates; atom-position view, not a bond or contact assignment")
        except (OSError, ValueError, RuntimeError) as exc:
            plt.close("all")
            warnings.append(cid + " pose figure unavailable: " + str(exc))
        # Export a reproducible molecular rendering script as well as the
        # exact-coordinate plot. PyMOL is optional; absence is not faked as PNG.
        script = destination / (cid + "_binding_pose.py")
        raw = destination / (cid + "_pymol_raw.png")
        script.write_text("from pathlib import Path\nfrom pymol import cmd\n" +
                          "here = Path(globals().get('__script__', __file__)).resolve().parent\n" +
                          "receptor = Path(" + repr(str(Path(receptor).resolve())) + ")\n" +
                          "if not receptor.is_file(): receptor = here.parent / 'sources/receptor.pdbqt'\n" +
                          "cmd.load(str(receptor), 'receptor')\n" +
                          "cmd.load(str(here.parent / 'poses' / " + repr(Path(row["selected_pose"]).name) + "), 'ligand')\n" +
                          "cmd.hide('everything', 'all')\ncmd.show('cartoon', 'receptor')\n" +
                          "cmd.show('sticks', 'ligand')\ncmd.select('pocket', 'byres (receptor within 5 of ligand)')\n" +
                          "cmd.show('sticks', 'pocket')\ncmd.color('gray70', 'receptor')\n" +
                          "cmd.color('cyan', 'ligand and elem C')\ncmd.bg_color('white')\n" +
                          "cmd.orient('ligand')\ncmd.zoom('ligand', 8)\ncmd.ray(1400, 1000)\n" +
                          "cmd.png(str(here / " + repr(raw.name) + "), dpi=220)\ncmd.quit()\n", encoding="utf-8")
        pymol = shutil.which("pymol")
        if pymol:
            try:
                result = subprocess.run([pymol, "-cq", str(script)], capture_output=True, text=True, timeout=120)
                script.with_suffix(".log").write_text(result.stdout + result.stderr, encoding="utf-8")
                if result.returncode != 0 or not raw.is_file():
                    raise RuntimeError("PyMOL did not produce the requested image; see rendering log")
                fig, ax = plt.subplots(figsize=(7, 5.5))
                ax.imshow(plt.imread(raw))
                ax.set_axis_off()
                ax.set_title(cid + " — molecular pose rendering")
                save(fig, cid + "_binding_pose.png", cid + " saved pose rendered by PyMOL")
            except Exception as exc:
                plt.close("all")
                warnings.append(cid + " PyMOL image unavailable: " + str(exc))
        elif not any("PyMOL is unavailable" in w for w in warnings):
            warnings.append("PyMOL is unavailable: exact-coordinate pose plots and molecular rendering scripts were exported instead.")
    return figures, warnings
