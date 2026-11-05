"""Server-rendered manually vetted five-step SIDock interface."""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

from .guided_workflow import (
    approved, active_run, compound_profile, load_workflow, protein_profile,
    receptor_approval_readiness, _json, _source,
)
from .manual_docking import (
    root as manual_root, settings, job_status, prepared_ligands,
    _input_metals, _pdbqt_metals,
)


def e(value) -> str:
    return html.escape(str(value), quote=True)


def _compound_file_id(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name).strip())


def flow_nav(project: Path, selected: int) -> str:
    labels = {
        1: "1. Project",
        2: "2. Protein & binding site",
        3: "3. Compounds",
        4: "4. Docking setup",
        5: "5. Docking & analysis",
    }
    elements = []
    for index, label in labels.items():
        if index == 1:
            state = "CREATED"
        elif index == 5:
            state = "READY" if approved(project, 4) else "WAITING"
        else:
            state = "APPROVED" if approved(project, index) else "REVIEW"
        tag = "<strong>" if index == selected else ""
        close = "</strong>" if index == selected else ""
        elements.append(
            f"<a style='padding:.6rem;border:1px solid #bcc7d4;border-radius:7px;"
            f"display:inline-block;margin:0 .4rem .4rem 0;"
            f"{'background:#e8eef5;' if index == selected else ''}' "
            f"href='/projects/{e(project.name)}/guided/{index}'>"
            f"{tag}{e(label)}{close}<small style='display:block'>{e(state)}</small></a>"
        )
    return "<div style='margin:1rem 0 1.5rem'>" + "".join(elements) + "</div>"


def action(project: Path, task: str, label: str, *, disabled=False, run_id=None) -> str:
    return (
        f"<form style='display:inline-block;margin:.5rem 1rem .5rem 0' "
        f"action='/projects/{e(project.name)}/guided/action' method='post'>"
        f"<input type='hidden' name='action' value='{e(task)}'>"
        + (f"<input type='hidden' name='run_id' value='{e(run_id)}'>" if run_id else "")
        + f"<button type='submit' {'disabled' if disabled else ''}>{e(label)}</button></form>"
    )


def approve_form(project: Path, step: int, label="Approve and continue") -> str:
    if approved(project, step):
        return f"<p><b>Step {step} approved for current inputs.</b> If you change anything, it requires reapproval.</p>"
    disabled = False
    reason = ""
    if step == 2:
        readiness = receptor_approval_readiness(project)
        disabled = not readiness["ready"]
        reason = " ".join(readiness["reasons"])
    return (
        f"<form method='post' action='/projects/{e(project.name)}/guided/approve'>"
        f"<input type='hidden' name='step' value='{step}'>"
        f"<button type='submit' {'disabled' if disabled else ''}>{e(label)}</button>"
        + (f"<p class='error'><b>Approval unavailable:</b> {e(reason)}</p>" if reason else "")
        + "</form>"
    )


def job_view(project: Path) -> str:
    job = job_status(project)
    if not job:
        return "<p>No manual job has been launched.</p>"
    status = job.get("status", "UNKNOWN")
    err = job.get("error")
    result = job.get("result")
    return (
        f"<p><b>Last action:</b> {e(job.get('action'))} — {e(status)}</p>"
        + (f"<pre class='error'>{e(err)}</pre>" if err else "")
        + (f"<details><summary>Output</summary><pre>{e(json.dumps(result,indent=2))}</pre></details>" if result else "")
        + f"<p><a href='/projects/{e(project.name)}/manual/log'>See worker log</a></p>"
    )


def step1(project: Path) -> str:
    protein = _source(project, "protein")
    compounds = _source(project, "compounds")
    return (
        "<h2>Step 1 — Create project</h2>"
        "<p>Project is created. Protein and compounds are deliberately collected in separate stages. "
        "No automatic docking starts after project creation.</p>"
        f"<p><b>Project:</b> {e(project.name)}</p>"
        f"<p><b>Protein:</b> {e(protein.name if protein else 'Not selected')}</p>"
        f"<p><b>Compounds:</b> {e(compounds.name if compounds else 'Not selected')}</p>"
        f"<p><a href='/projects/{e(project.name)}/guided/2'>Continue to Step 2: protein and binding site</a></p>"
    )


def step2(project: Path) -> str:
    profile = protein_profile(project)
    choice = _json(manual_root(project) / "inputs/protein_choice.json")
    job = job_status(project)
    receptor = manual_root(project) / "receptor.pdbqt"
    content = [
        "<h2>Step 2 — Protein, crystal ligand and site</h2>",
        "<p>Select a deposited component to use as the grid reference, or choose no component "
        "and enter grid coordinates. Components with potential protein/cap roles are selectable "
        "only with explicit acknowledgment; choosing them is not native-ligand validation. "
        "Review water, ions, solvents and modifications before preparing the receptor.</p>",
    ]
    if not profile["available"]:
        content.append(
            f"<form method='post' enctype='multipart/form-data' "
            f"action='/projects/{e(project.name)}/guided/protein-upload'>"
            "<label>Upload protein PDB <input type='file' name='protein' accept='.pdb,.ent' required></label>"
            "<button>Upload protein</button></form>"
        )
        return "".join(content)
    content.append(f"<p><b>Protein source:</b> {e(profile['protein'])}</p>")
    content.extend(f"<small style='display:block'>{e(line)}</small>" for line in profile["headers"])
    content.append(
        f"<form method='post' action='/projects/{e(project.name)}/guided/protein-review'>"
        "<h3>Detected components</h3>"
    )
    if not profile["groups"]:
        content.append("<p>No nonstandard components detected. Enter a manual grid center below.</p>")
    else:
        waters = [g for g in profile["groups"] if g["role"] == "WATER"]
        if waters:
            content.append(f"<p>{len(waters)} crystallographic water residues; use the water checkbox below.</p>")
        content.append(
            "<table><tr><th>Component</th><th>Role</th><th>Atoms</th><th>Remove from receptor?</th></tr>"
        )
        old_remove = set(choice.get("removed_residues", []))
        for group in profile["groups"]:
            if group["role"] == "WATER":
                continue
            key = group["key"]
            protected = group["protected"]
            selected = key in old_remove if choice else group["default_remove"]
            remove_input = "Protected protein modification (retained)"
            if group["role"] == "WATER":
                remove_input = "Controlled by the water checkbox"
            elif not protected:
                remove_input = (
                    f"<input type='checkbox' name='remove_groups' value='{e(key)}' "
                    f"{'checked' if selected else ''}> Remove"
                )
            content.append(
                f"<tr><td><b>{e(key)}</b></td><td>{e(group['role'])}</td>"
                f"<td>{e(group['atoms'])}</td><td>{remove_input}</td></tr>"
            )
        content.append("</table>")

    options = ["<option value=''>No reference component — enter grid coordinates manually</option>"]
    for group in profile["ligands"]:
        key = group["key"]
        options.append(
            f"<option value='{e(key)}' {'selected' if key == choice.get('ligand_key') else ''}>"
            f"{e(key)} ({group['atoms']} atoms; {e(group['role'])})</option>"
        )
    current = settings(project)
    padding = choice.get("padding", 8.0)
    checked_waters = choice.get("remove_waters", True)
    content.append(
        f"<label><b>Grid reference component (manual selection):</b> "
        f"<select name='ligand_key'>{''.join(options)}</select></label>"
        "<label><input type='checkbox' name='acknowledge_non_ligand_component' value='on' "
        + ("checked" if choice.get("scientific_acknowledgment") else "") + "> "
        "I understand that selecting or removing a potential protein cap/buffer as a reference "
        "does not establish a native ligand binding site, and may change protein chemistry.</label>"
        f"<label><input type='checkbox' name='remove_waters' value='on' "
        f"{'checked' if checked_waters else ''}> Remove crystallographic waters</label>"
        "<p><small>A verified ligand is extracted and removed from the receptor. "
        "Potential protein caps such as ACE are retained unless you expressly tick Remove "
        "in the component table. All changes are recorded. Metal ions remain unless "
        "you choose otherwise.</small></p>"
        f"<label>Binding-site padding (Å) <input type='number' step='0.5' min='1' max='20' "
        f"name='padding' value='{e(padding)}'></label>"
        f"<label><input type='checkbox' name='exploratory_metal_free' value='on' "
        f"{'checked' if current['allow_missing_metal'] else ''}> "
        "I explicitly choose an exploratory model if catalytic metal ions are absent "
        "from the prepared receptor.</label>"
        "<h3>Manual center (used only when no ligand is selected)</h3>"
        + "".join(
            f"<label>{e(axis)} <input name='center_{axis.lower()}' type='number' step='0.001' "
            f"value='{e(current['center_'+axis.lower()] if choice else '')}'></label>"
            for axis in ("X","Y","Z")
        )
        + "<button>Save protein/site review and generate cleaned PDB preview</button></form>"
    )
    if choice:
        nickel_remains = any(
            g["code"] == "NI" and g["key"] not in choice.get("removed_residues", [])
            for g in profile.get("groups", [])
        )
        content.append(
            f"<p><b>Binding-site mode:</b> {e(choice.get('mode'))}; "
            f"<b>center:</b> {e(choice.get('center'))}; "
            f"<b>grid points:</b> {e(choice.get('npts'))}</p>"
            f"<p>Removed components: {e(', '.join(choice.get('removed_residues', [])) or 'None')}</p>"
            + (f"<p><b>Scientific warning:</b> {e(choice['site_warning'])}</p>"
               if choice.get('site_warning') else "")
            + f"<p><a href='/projects/{e(project.name)}/artifacts/manual_docking/inputs/protein_clean.pdb'>"
            "Inspect cleaned receptor PDB</a>"
            + (f" | <a href='/projects/{e(project.name)}/artifacts/manual_docking/inputs/crystal_ligand.pdb'>"
               "Selected reference component PDB</a>" if choice.get("ligand_key") else "")
            + "</p>"
        )
        selected_component = next(
            (g for g in profile.get("groups", []) if g["key"] == choice.get("ligand_key")), None
        )
        if selected_component and selected_component["role"] != "LIGAND_CANDIDATE":
            content.append(
                "<p class='error'><b>Grid-reference caution:</b> The selected "
                f"component {e(selected_component['key'])} is classified as "
                f"{e(selected_component['role'])}, not a verified bound ligand. "
                "The grid may not cover the biologically relevant site. Verify the "
                "reference against experimental ligand/site evidence.</p>"
            )
        # Structural-only export stays available for inspection, separate from
        # the complete preparation action below.
        ob_candidate = manual_root(project) / "inputs/receptor_openbabel_candidate.pdbqt"
        ob_audit = _json(manual_root(project) / "inputs/receptor_openbabel_audit.json")
        content.append(
            "<details><summary>Structural export and audit (advanced)</summary>"
            "<p>Use Open Babel to convert the reviewed PDB directly to a "
            "<b>rigid receptor PDBQT</b>. It does not require residue templates "
            "for modified residues. This creates an inspectable file without "
            "approving metal parameters, protonation or previously unresolved atoms.</p>"
            + action(project, "receptor_openbabel",
                     "Generate rigid receptor PDBQT with Open Babel",
                     disabled=job.get("status") == "RUNNING")
        )
        if ob_candidate.is_file():
            content.append(
                f"<p><b>Open Babel PDBQT generated:</b> "
                f"<a href='/projects/{e(project.name)}/artifacts/manual_docking/inputs/receptor_openbabel_candidate.pdbqt'>"
                "Open/download receptor PDBQT</a>"
                + (f" | <a href='/projects/{e(project.name)}/artifacts/manual_docking/inputs/receptor_openbabel_audit.json'>"
                   "Review export audit</a>" if ob_audit else "")
                + "</p>"
                + (f"<p><b>Atom retention:</b> "
                   f"{e(ob_audit.get('input_heavy_atom_count', 'unknown'))} original heavy atoms; "
                   f"{e(len(ob_audit.get('missing_input_heavy_atoms', [])))} unmatched on export. "
                   f"<b>Metal chemistry:</b> {e(ob_audit.get('metal_parameterization', 'NOT_VALIDATED'))}.</p>"
                   if ob_audit else
                   "<p>Output exists but its provenance audit is absent; review before use.</p>")
            )
        content.append(
            "<p><small>Step 2 approval remains separate. A PDBQT file alone "
            "does not validate protein charge, protonation or coordination chemistry.</small></p></details>"
        )
        from .manual_docking import parameter_file_path
        custom = parameter_file_path(project)
        if not nickel_remains:
            content.append(
                "<h3>Prepare nonmetal receptor</h3>"
                + action(project, "receptor", "Prepare nonmetal receptor with Open Babel",
                         disabled=job.get("status") == "RUNNING")
            )
        content.append(
            "<h3>Prepare receptor with Open Babel</h3>"
            "<p>Use this project's parameter library, then prepare the receptor. "
            "This completes structural export, hydrogen/charge preparation and atom typing, "
            "while preserving the reviewed heavy atoms and metal coordinates. "
            "Existing matching preparations are reused.</p>"
            f"<p><b>AutoDock4 parameter library:</b> {e(custom.name) if custom else 'Select a library below'}</p>"
            f"<form method='post' action='/projects/{e(project.name)}/guided/parameters' "
            "enctype='multipart/form-data'>"
            "<input type='hidden' name='next_step' value='2'>"
            "<label>Complete AutoDock4 parameter library (.dat) "
            "<input type='file' name='upload' accept='.dat,.txt'></label>"
            "<label>Or existing project-local .dat file "
            "<input name='existing_path' placeholder='input/AD4_parameters.dat'></label>"
            "<button>Use parameter library</button></form>"
            + action(project, "receptor_openbabel_pipeline", "Prepare receptor with Open Babel",
                     disabled=not custom or job.get("status") == "RUNNING")
            + "<p>Preparation does not start docking or approve the receptor. "
            "For exploratory approval, run the separate one-ligand technical check below. "
            "Metal charges and coordination remain unvalidated.</p>"
            + action(project, "exploratory_smoke", "Run one-ligand technical check",
                     disabled=not custom or not (manual_root(project) / "inputs/receptor_openbabel_exploratory.pdbqt").is_file()
                     or not (manual_root(project) / "inputs/crystal_ligand.pdb").is_file()
                     or job.get("status") == "RUNNING")
        )
        if not (manual_root(project) / "inputs/crystal_ligand.pdb").is_file():
            content.append("<p>The technical check needs a saved reference component. "
                           "Select and save one in the binding-site review above.</p>")
        content.append(
            "<p><b>Alternatively:</b> import a previously prepared and chemically verified "
            "PDBQT. This is useful for metal-dependent proteins.</p>"
            f"<form method='post' action='/projects/{e(project.name)}/guided/receptor-import' "
            "enctype='multipart/form-data'>"
            "<label>Upload receptor PDBQT <input name='upload' type='file' accept='.pdbqt'></label>"
            "<label>Or select an existing project-local PDBQT "
            "<input name='existing_path' placeholder='provenance/4ubp_protein_diagnostic.pdbqt'></label>"
            "<button>Import receptor PDBQT</button></form>"
        )
    if receptor.is_file():
        content.append(
            f"<p><b>Prepared receptor:</b> {e(receptor)} "
            f"({receptor.stat().st_size:,} bytes)</p>"
        )
        metal_audit_file = manual_root(project) / "inputs/Ni_reconstruction_audit.json"
        receptor_review = _json(manual_root(project) / "inputs/receptor_review.json")
        if receptor_review.get("metal_model") and metal_audit_file.is_file():
            metal_audit = _json(metal_audit_file)
            content.append(
                "<p><b>Nickel reconstructed:</b> "
                f"{e(metal_audit.get('NI_count'))} ions. "
                f"Ni–Ni distance: {e(metal_audit.get('Ni_Ni_distance_A'))} Å. "
                "<b>Fixed-centre nonbonded model — NOT validated by native redocking.</b> "
                f"<a href='/projects/{e(project.name)}/artifacts/manual_docking/inputs/Ni_reconstruction_audit.json'>"
                "Review metal parameters, charges and sources</a></p>"
            )
        original_metals = _input_metals(project)
        prepared_metals = _pdbqt_metals(project)
        if original_metals and prepared_metals != original_metals:
            content.append(
                "<p class='error'><b>Metal mismatch:</b> original PDB has "
                f"{original_metals} nickel ions; prepared receptor has {prepared_metals}. "
                "A metal-free result must be explicitly designated exploratory "
                "under Step 4 and cannot be reported as native metal-site docking.</p>"
            )
    preparation_log = {
        "receptor": "receptor_openbabel_standard.log",
        "receptor_metal_free": "receptor_openbabel_metal_free.log",
        "receptor_ni_intermediate": "receptor_openbabel_Ni_intermediate.log",
    }.get(job.get("action"))
    if preparation_log and not (manual_root(project) / "logs" / preparation_log).exists() and job.get("status") != "RUNNING":
        content.append("<details><summary>Previous receptor preparation attempt (history)</summary>"
                       + job_view(project) + "</details>")
    else:
        content.append(job_view(project))
    if job.get("status") == "RUNNING":
        content.append(
            "<p>Step 2 receptor export is running. Refreshing automatically.</p>"
            "<script>setTimeout(function(){location.reload();}, 3500);</script>"
        )
    from .exploratory_approval import MODE, readiness as exploratory_readiness
    from .guided_workflow import approval_mode
    ex = exploratory_readiness(project)
    if approval_mode(project) == MODE:
        content.append(
            "<p><b>Step 2 technically approved for exploratory docking.</b> "
            "This is NOT scientifically validated. The Ni partial charges and "
            "metal-coordination scoring remain provisional.</p>"
            "<p><a href='/projects/" + e(project.name) + "/guided/3'>"
            "Continue to Step 3 — prepared compounds</a></p>"
        )
    else:
        content.append(
            "<h3>Approve Open Babel receptor for exploratory docking</h3>"
            "<p>Use the existing Open Babel candidate and successful "
            "AutoGrid4/AutoDock4 smoke-test evidence. This does not overwrite "
            "the standard receptor or imply scientifically validated metal "
            "binding energies.</p>"
            + (f"<p><b>Technical evidence:</b> "
               f"{e('Ready — successful smoke test and unchanged hashes' if ex['ready'] else 'Not ready')}</p>")
        )
        if not ex["ready"]:
            content.append(
                "<p>Complete Open Babel preparation and the one-ligand technical check above "
                "before approving. If either fails, inspect the latest worker error.</p>"
                "<details><summary>Outstanding exploratory prerequisites:</summary><p>"
                + e("; ".join(ex["reasons"])) + "</p></details>"
            )
        content.append(
            f"<form action='/projects/{e(project.name)}/guided/approve-exploratory' method='post'>"
            "<label><input type='checkbox' name='acknowledgement' "
            "value='I_ACKNOWLEDGE_UNVALIDATED_EXPLORATORY_DOCKING' required> "
            "I understand this is a <b>technical, exploratory AutoDock4 model</b>. "
            "The protein charge/protonation and metal coordination are unvalidated; "
            "its docking energies must not be reported as validated affinities.</label>"
            f"<button type='submit' {'disabled' if not ex['ready'] or job.get('status') == 'RUNNING' else ''}>"
            "Approve for exploratory docking and continue</button></form>"
        )
    if receptor.is_file() and approval_mode(project) != MODE:
        content.append("<h3>Approve independently reviewed imported receptor</h3>" + approve_form(project, 2))
    return "".join(content)


def step3(project: Path) -> str:
    profile = compound_profile(project)
    job = job_status(project)
    html_bits = [
        "<h2>Step 3 — Compounds and minimization</h2>",
        "<p>Upload an SDF library and inspect the compound count and IDs. "
        "Minimization and PDBQT generation are explicitly started here, only once. "
        "The workflow rejects missing or duplicate identifiers instead of silently omitting compounds.</p>"
    ]
    if not profile["available"]:
        html_bits.append(
            f"<form enctype='multipart/form-data' method='post' "
            f"action='/projects/{e(project.name)}/guided/compound-upload'>"
            "<label>Compound SDF <input type='file' name='compounds' accept='.sdf' required></label>"
            "<button>Upload compounds</button></form>"
        )
        return "".join(html_bits)
    html_bits.append(
        f"<p><b>Source:</b> {e(profile['source'])}</p>"
        f"<p><b>Molecules:</b> {profile['total']} | "
        f"<b>Prepared:</b> {profile['prepared_count']} | "
        f"<b>Verified preparations:</b> {'Yes' if profile['prepared'] else 'No'}</p>"
    )
    if profile["errors"]:
        html_bits.append(
            "<p class='error'>Fix the source SDF before preparation: "
            + e("; ".join(profile["errors"])) + "</p>"
        )
    html_bits.append(
        "<p><b>Preparation is independent of Step 2.</b> You may minimize and convert "
        "ligands to PDBQT now, even while receptor preparation is unresolved. "
        "Step 3 approval, grid generation and docking remain blocked until Step 2 "
        "has a verified, approved receptor.</p>"
    )
    existing_prepared = project / "ligands/prepared"
    busy = job.get("status") == "RUNNING"
    if existing_prepared.is_dir() and any(existing_prepared.glob("*.pdbqt")):
        html_bits.append(action(
            project, "reuse_compounds",
            "Reuse existing project PDBQT files (no minimization)",
            disabled=bool(profile["errors"]) or busy or profile["prepared"]
        ))
    html_bits.append(action(
        project, "prepare",
        "Prepared ligands verified (no rerun needed)" if profile["prepared"]
        else "Minimize and prepare compounds",
        disabled=bool(profile["errors"]) or busy or profile["prepared"]
    ))
    html_bits.append(job_view(project))
    if busy:
        html_bits.append(
            "<p><b>Ligand preparation is running.</b> This page will refresh "
            "automatically to show the latest verified outputs.</p>"
            "<script>setTimeout(function(){location.reload();}, 3500);</script>"
        )
    html_bits.append(
        "<details><summary>Compound IDs and minimized/prepared files</summary><table>"
        "<tr><th>Compound ID</th><th>Minimized SDF</th><th>Docking PDBQT</th></tr>"
        + "".join(
            "<tr><td>" + e(name) + "</td><td>"
            + (f"<a href='/projects/{e(project.name)}/artifacts/manual_docking/ligands/source/{e(_compound_file_id(name))}.sdf'>View</a>"
               if (manual_root(project)/"ligands/source"/(_compound_file_id(name)+".sdf")).is_file() else "Not prepared")
            + "</td><td>"
            + (f"<a href='/projects/{e(project.name)}/artifacts/manual_docking/ligands/prepared/{e(_compound_file_id(name))}.pdbqt'>View</a>"
               if (manual_root(project)/"ligands/prepared"/(_compound_file_id(name)+".pdbqt")).is_file() else "Not prepared")
            + "</td></tr>"
            for name in profile["names"]
        ) + "</table></details>"
    )
    if not profile["prepared"]:
        html_bits.append(
            "<p>Step 3 approval unavailable: complete and verify every ligand preparation first.</p>"
        )
    elif not approved(project, 2):
        html_bits.append(
            "<p>All ligands prepared and verified. <b>Step 3 approval is pending "
            "a current Step 2 scientific or technical-exploratory approval.</b> "
            "No ligand re-preparation is necessary.</p>"
        )
    else:
        from .guided_workflow import approval_mode
        from .exploratory_approval import MODE
        if approval_mode(project) == MODE:
            html_bits.append(
                "<p><b>EXPLORATORY_UNVALIDATED_NI_MODEL:</b> approving these "
                "prepared ligands permits technical grid/docking workflow only. "
                "Metal-site binding predictions are not scientifically validated.</p>"
            )
        html_bits.append(approve_form(project, 3))
    return "".join(html_bits)


def step4(project: Path) -> str:
    t = settings(project)
    run = active_run(project)
    job = job_status(project)
    html_bits = [
        "<h2>Step 4 — Docking setup and grid</h2>",
        "<p>Confirm the docking center, AutoGrid grid points and genetic algorithm settings. "
        "Then generate maps and review the GPF. This stage does not dock compounds.</p>"
    ]
    for k in ("center_x", "center_y", "center_z", "npts_x", "npts_y", "npts_z"):
        if k not in t:
            return "<p class='error'>Step 2 must define the docking site first.</p>"
    fields = []
    for name, label, step in (
        ("center_x","Center X (Å)", "0.001"),
        ("center_y","Center Y (Å)", "0.001"),
        ("center_z","Center Z (Å)", "0.001"),
        ("npts_x","Grid points X (even)", "2"),
        ("npts_y","Grid points Y (even)", "2"),
        ("npts_z","Grid points Z (even)", "2"),
        ("spacing","Spacing (Å)", "0.025"),
        ("ga_runs","GA runs per compound", "1"),
        ("parallel_jobs","Parallel compound jobs", "1"),
        ("ga_num_evals","Maximum evaluations", "1"),
        ("ga_pop_size","GA population", "1"),
        ("ga_num_generations","Maximum generations", "1"),
    ):
        fields.append(
            f"<label>{e(label)} <input name='{name}' type='number' step='{step}' "
            f"value='{e(t[name])}' required></label>"
        )
    exploratory = (
        f"<label><input type='checkbox' value='on' name='allow_missing_metal' "
        f"{'checked' if t['allow_missing_metal'] else ''}>"
        " I explicitly acknowledge that if nickel is absent, the run is an "
        "exploratory metal-free model, not metal-site docking.</label>"
    )
    dims = [t[k]*t["spacing"] for k in ("npts_x","npts_y","npts_z")]
    html_bits.append(f"<p>Current box dimensions: {dims[0]:.2f} × {dims[1]:.2f} × {dims[2]:.2f} Å.</p>")
    html_bits.append(
        f"<form method='post' action='/projects/{e(project.name)}/guided/settings'>"
        + "".join(fields) + exploratory
        + "<button>Save docking setup</button></form>"
    )
    from .guided_workflow import approval_mode
    from .exploratory_approval import MODE
    if approval_mode(project) == MODE:
        html_bits.append(
            "<p><b>EXPLORATORY_UNVALIDATED_NI_MODEL.</b> AutoGrid and AutoDock "
            "will use the hash-approved Open Babel exploratory receptor, NOT the "
            "standard receptor. Review the complete binding-site coverage, "
            "grid dimensions and search parameters; the small smoke grid is "
            "not appropriate for a compound library. Results remain exploratory.</p>"
        )
    from .manual_docking import parameter_file_path
    parameter = parameter_file_path(project)
    html_bits.append(
        "<h3>Optional custom AutoDock4 force-field library</h3>"
        "<p>For nonstandard atoms or metals such as Ni, supply a fully specified, scientifically "
        "justified AD4 parameter library. SIDock applies the same file to AutoGrid and AutoDock. "
        "It does not invent charges or metal parameters.</p>"
        + (f"<p>Current custom library: {e(parameter.name)} ({parameter.stat().st_size:,} bytes)</p>"
           if parameter else "<p>Using standard built-in AutoDock4 parameters.</p>")
        + f"<form method='post' action='/projects/{e(project.name)}/guided/parameters' "
        "enctype='multipart/form-data'>"
        "<label>Upload complete AD4 parameter file <input type='file' name='upload' accept='.dat,.txt'></label>"
        "<label>Or existing project-local path <input name='existing_path' placeholder='input/custom_ad4_parameters.dat'></label>"
        "<button>Use custom parameter library</button></form>"
    )
    html_bits.append(action(
        project, "grid", "Generate AutoGrid4 maps (reuse if unchanged)",
        disabled=not (approved(project, 2) and approved(project, 3))
        or job.get("status") == "RUNNING",
    ))
    if run and (run/"grid/production/grid_provenance.json").exists():
        rel = run.relative_to(project).as_posix()
        html_bits.append(
            f"<p><b>Maps ready for run:</b> {e(run.name)}</p>"
            f"<a href='/projects/{e(project.name)}/artifacts/{e(rel)}/grid/production/receptor.gpf'>"
            "Inspect AutoGrid parameters (GPF)</a>"
        )
    html_bits.append(job_view(project))
    html_bits.append(approve_form(project, 4))
    return "".join(html_bits)


def step5(project: Path, *, selected_run=None) -> str:
    job = job_status(project)
    current_run = active_run(project)
    from .guided_analysis import resolve_run, analysis_status
    run = resolve_run(project, selected_run) if selected_run else current_run
    content = [
        "<h2>Step 5 — Docking and analysis</h2>",
        "<p>The run is unlocked only after Steps 2–4 have been approved for their "
        "current files and settings. Existing completed compound logs are not rerun. "
        "After successful docking, pose/cluster analysis, PLIP, figures and the draft write-up start automatically.</p>",
    ]
    from .guided_workflow import approval_mode
    from .exploratory_approval import MODE
    if approval_mode(project) == MODE:
        content.append(
            "<p><b>EXPLORATORY_UNVALIDATED_NI_MODEL</b> — All AutoDock results "
            "from this branch are technical/exploratory. Charges, protonation, "
            "metal coordination and derived binding energies are not validated. "
            "Review grid and genetic-algorithm settings before starting a full run.</p>"
        )
    if not approved(project, 4):
        content.append("<p class='error'>Docking is locked until the previous steps are approved.</p>")
    if not selected_run or (current_run and current_run.name == selected_run):
        content.append(action(
            project, "dock", "Start / resume AutoDock4 jobs",
            disabled=not approved(project, 4) or job.get("status") == "RUNNING",
        ))
    else:
        content.append("<p><b>Viewing saved run:</b> " + e(run.name) +
                       ". Its analysis uses the saved receptor and DLG files.</p>")
    content.append(job_view(project))
    if job.get("status") == "RUNNING":
        content.append("<p>The worker is running. This page refreshes automatically.</p>"
                       "<script>setTimeout(function(){location.reload();}, 3500);</script>")
    if run:
        rel = run.relative_to(project).as_posix()
        result = run / "docking_result.json"
        if result.is_file():
            report = _json(result)
            content.append(
                f"<p><b>Current run:</b> {e(run.name)} | "
                f"Completed: {report.get('completed',0)}/{report.get('total',0)}</p>"
            )
        csvfile = run / "results.csv"
        if csvfile.exists():
            content.append(
                f"<p><a href='/projects/{e(project.name)}/artifacts/{e(rel)}/results.csv'>"
                "Download docking results CSV</a></p>"
            )
        from .manual_docking import complete_dlg_energy
        completed = any(complete_dlg_energy(p) is not None for p in (run / "docking/dlg").glob("*.dlg"))
        content.append("<h3>Analyse saved docking results</h3>"
                       "<p>Run or retry post-processing of this run's completed DLG files. "
                       "This reuses the prepared receptor/ligand files and does not generate grids or rerun docking.</p>")
        content.append(action(project, "analyze", "Analyse existing results — no redocking",
                              run_id=run.name, disabled=not completed or job.get("status") == "RUNNING"))
        result = analysis_status(project, run)
        if result:
            content.append(f"<p><b>Analysis:</b> {e(result.get('status'))} | "
                           f"Poses analysed: {e(result.get('analysed', 0))}/{e(result.get('total', 0))} | "
                           f"PLIP complete: {e(result.get('plip_complete', 0))}</p>")
            if result.get("status") == "STALE":
                content.append("<p class='error'>" + e(result.get("stale_reason")) +
                               ". Rerun analysis before using the report.</p>")
            elif result.get("analysis_directory"):
                folder = result["analysis_directory"]
                for name, label in (("scientific_report.html", "View analysis, figures and draft write-up"),
                                    ("scientific_report.md", "Download editable write-up"),
                                    ("results.csv", "Download scores, clusters and interactions"),
                                    ("evidence.json", "Download evidence JSON")):
                    content.append(f"<p><a href='/projects/{e(project.name)}/artifacts/{e(folder)}/{name}'>"
                                   + label + "</a></p>")
                content.append(f"<p><a href='/projects/{e(project.name)}/artifacts/{e(result['bundle'])}'>"
                               "Download all analysis files (ZIP)</a></p>")
                diagrams = sorted(name for name in result.get("artifacts", {})
                                  if "_interactions_2d" in Path(name).name and name.endswith(".png"))
                if diagrams:
                    content.append("<details><summary>2D protein–ligand interaction figures</summary>" +
                                   "".join(f"<p><a href='/projects/{e(project.name)}/artifacts/{e(name)}'>"
                                           f"{e(Path(name).name)}</a></p>" for name in diagrams) + "</details>")
            if result.get("warnings"):
                content.append("<details><summary>Analysis warnings and missing outputs</summary><ul>" +
                               "".join("<li>" + e(w) + "</li>" for w in result["warnings"]) + "</ul></details>")
    previous_runs = sorted((p for p in (manual_root(project) / "runs").glob("*")
                            if p.is_dir() and (p / "protocol.json").is_file() and p != run),
                           key=lambda p: p.stat().st_mtime, reverse=True)
    if previous_runs:
        content.append("<details><summary>Other saved runs</summary>" + "".join(
            f"<p><a href='/projects/{e(project.name)}/guided/5?analysis_run={e(p.name)}'>"
            f"View or analyse run {e(p.name)}</a></p>" for p in previous_runs) + "</details>")
    return "".join(content)


def render(project: Path, step: int, *, analysis_run=None) -> str:
    if step not in (1,2,3,4,5):
        raise ValueError("Invalid guided step")
    views = {1: step1, 2: step2, 3: step3, 4: step4, 5: step5}
    content = (
        f"<h1>SIDock — Reviewed docking workflow</h1><p><b>Project:</b> {e(project.name)} "
        f"| <a href='/projects/{e(project.name)}'>Project overview</a></p>"
        + flow_nav(project, step) + (step5(project, selected_run=analysis_run) if step == 5 else views[step](project))
    )
    return content + f"<hr><p><a href='/projects/{e(project.name)}/guided/{min(step+1,5)}'>Next step</a></p>"
