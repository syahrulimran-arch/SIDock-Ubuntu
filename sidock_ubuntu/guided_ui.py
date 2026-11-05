"""Server-rendered manually vetted five-step SIDock interface."""
from __future__ import annotations

import html
import json
import re
from pathlib import Path

from .guided_workflow import (
    approved, active_run, compound_profile, load_workflow, protein_profile,
    _json, _source,
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
        5: "5. Docking run",
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


def action(project: Path, task: str, label: str, *, disabled=False) -> str:
    return (
        f"<form style='display:inline-block;margin:.5rem 1rem .5rem 0' "
        f"action='/projects/{e(project.name)}/guided/action' method='post'>"
        f"<input type='hidden' name='action' value='{e(task)}'>"
        f"<button type='submit' {'disabled' if disabled else ''}>{e(label)}</button></form>"
    )


def approve_form(project: Path, step: int, label="Approve and continue") -> str:
    if approved(project, step):
        return f"<p><b>Step {step} approved for current inputs.</b> If you change anything, it requires reapproval.</p>"
    return (
        f"<form method='post' action='/projects/{e(project.name)}/guided/approve'>"
        f"<input type='hidden' name='step' value='{step}'>"
        f"<button type='submit'>{e(label)}</button></form>"
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
        # Meeko-free receptor export: available for ANY reviewed protein.
        # The file is an inspectable rigid candidate, not a metal-site approval.
        ob_candidate = manual_root(project) / "inputs/receptor_openbabel_candidate.pdbqt"
        ob_audit = _json(manual_root(project) / "inputs/receptor_openbabel_audit.json")
        content.append(
            "<h3>Generate receptor PDBQT without Meeko</h3>"
            "<p>Use Open Babel to convert the reviewed PDB directly to a "
            "<b>rigid receptor PDBQT</b>. It does not require Meeko templates "
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
            "does not validate protein charge, protonation or coordination chemistry.</small></p>"
        )
        content.append(
            "<h3>Prepare receptor PDBQT</h3>"
            "<p>Preparation runs as an independent step. Meeko will stop and show its actual error "
            "if it cannot parameterize nonstandard components; it will not silently delete them.</p>"
            + ("<p class='error'><b>Standard Meeko is disabled:</b> this reviewed PDB contains "
               "catalytic nickel, which Meeko cannot parameterize automatically. "
               "For metal-dependent docking, import a properly parameterized "
               "metal-containing receptor PDBQT. Alternatively, use the separate "
               "exploratory metal-free option below. The original receptor is untouched.</p>"
               if nickel_remains else "")
            + action(project, "receptor", "Prepare receptor using Meeko",
                     disabled=nickel_remains or job.get("status") == "RUNNING")
        )
        if nickel_remains:
            content.append(
                "<h3>Optional: prepare a metal-free exploratory receptor</h3>"
                "<p>Meeko does not support NI in this receptor. This explicit option generates "
                "a separate temporary PDB without the crystallographic nickel ions; it preserves your "
                "original PDB and reviewed cleaned PDB, records a removal audit, and never "
                "runs docking. <b>This is not suitable for interpreting nickel-dependent "
                "urease inhibition.</b></p>"
                f"<form method='post' action='/projects/{e(project.name)}/guided/action'>"
                "<input type='hidden' name='action' value='receptor_metal_free'>"
                "<label><input type='checkbox' required name='confirm_metal_free' "
                "value='I_ACKNOWLEDGE_METAL_FREE'> "
                "I explicitly approve generating a metal-free receptor for exploratory calculations; "
                "I understand nickel-dependent interactions cannot be evaluated.</label>"
                f"<button type='submit' {'disabled' if job.get('status') == 'RUNNING' else ''}>"
                "Prepare exploratory metal-free receptor (no docking)</button></form>"
            )
        if nickel_remains:
            from .nickel_reconstruction import sites as nickel_sites
            from .manual_docking import parameter_file_path
            site_list = nickel_sites(project)
            intermediate = manual_root(project) / "inputs/protein_meeko_Ni_intermediate.pdbqt"
            metadata = _json(manual_root(project) / "inputs/Ni_intermediate.json")
            validated_intermediate = (
                intermediate.is_file()
                and metadata.get("reviewed_clean_sha256") == choice.get("clean_sha256")
            )
            custom = parameter_file_path(project)
            content.append(
                "<h3>Reconstruct crystallographic nickel centre (Ni retained)</h3>"
                "<p>Step A prepares the protein with nickel omitted from a temporary "
                "Meeko input; this is <b>not</b> the final metal-free receptor. "
                "Step B restores the Ni ions at their exact crystallographic coordinates. "
                "It requires a complete, scientifically sourced AutoDock4 parameter library "
                "and individually specified Ni partial charges; neither is guessed.</p>"
                + action(project, "receptor_ni_intermediate",
                         "A. Prepare Ni-free Meeko intermediate (does not install receptor)",
                         disabled=job.get("status") == "RUNNING")
            )
            content.append(
                f"<p><b>Intermediate:</b> "
                f"{'Ready for metal reconstruction' if validated_intermediate else 'Not prepared / stale'}</p>"
                "<p><b>Ni coordinates preserved from the reviewed PDB:</b></p>"
                "<table><tr><th>Site</th><th>X (Å)</th><th>Y (Å)</th><th>Z (Å)</th></tr>"
                + "".join(
                    f"<tr><td>{e(x['key'])}</td>"
                    + "".join(f"<td>{v:.3f}</td>" for v in x["xyz"])
                    + "</tr>" for x in site_list
                ) + "</table>"
            )
            content.append(
                f"<p>Custom AD4 parameter library: "
                f"{e(custom.name) if custom else 'Missing — upload a fully specified library below'}</p>"
                f"<form method='post' action='/projects/{e(project.name)}/guided/parameters' "
                "enctype='multipart/form-data'>"
                "<input type='hidden' name='next_step' value='2'>"
                "<label>Upload complete AutoDock4 parameter library (.dat) "
                "<input type='file' name='upload' accept='.dat,.txt'></label>"
                "<label>Or select existing project-local .dat file "
                "<input name='existing_path' placeholder='input/AD4_nickel_parameters.dat'></label>"
                "<button>Use parameter library</button></form>"
            )
            content.append(
                f"<form method='post' action='/projects/{e(project.name)}/guided/ni-reconstruct'>"
                "<p><b>B. Reconstruct nickel from crystallographic coordinates</b></p>"
                "<label>Nickel AutoDock4 type (must match atom_par in your library) "
                "<input type='text' name='atom_type' maxlength='2' placeholder='Ni' required></label>"
                "<label>Ni partial charges (e), individually and explicitly assigned "
                "<input name='site_charges' required "
                f"placeholder='{e(','.join(x['key'] + '=...' for x in site_list))}'></label>"
                "<label>Parameter source / reference "
                "<input name='parameter_source' required placeholder='Publication/DOI and parameterization method'></label>"
                "<label>Charge source / reference "
                "<input name='charge_source' required placeholder='Publication/DOI or documented charge model'></label>"
                "<label><input type='checkbox' required name='confirm_nickel_model' "
                "value='I_ACKNOWLEDGE_NI_MODEL'> "
                "I acknowledge that a fixed, nonbonded Ni model requires scientific review "
                "and native-ligand redocking before biological interpretation.</label>"
                f"<button type='submit' {'disabled' if not validated_intermediate or not custom or job.get('status') == 'RUNNING' else ''}>"
                "B. Restore Ni and install reviewed receptor (no docking)</button></form>"
            )
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
    content.append(job_view(project))
    from .receptor_diagnostics import read_receptor_diagnostics
    diagnostic_action = job.get("action") in {
        "receptor", "receptor_metal_free", "receptor_ni_intermediate"
    }
    diagnostics = read_receptor_diagnostics(
        project,
        exploratory_metal_free=job.get("action") == "receptor_metal_free",
        nickel_intermediate=job.get("action") == "receptor_ni_intermediate",
    )
    if diagnostics.get("available") and (
        diagnostic_action and job.get("status") in {"FAILED", "EXITED"}
    ):
        content.append(
            "<h3>Why Meeko failed</h3>"
            "<p>The error below is read from this project's preparation log. "
            "The original PDB and reviewed cleaned PDB have not been modified.</p>"
        )
        if diagnostics["unknown_components"]:
            content.append(
                "<p><b>Unsupported component templates:</b> "
                + e(", ".join(diagnostics["unknown_components"])) + "</p>"
            )
        if diagnostics["template_mismatches"]:
            content.append(
                "<p><b>Unmatched residues:</b> "
                + e(", ".join(diagnostics["template_mismatches"])) + "</p>"
            )
        for note in diagnostics["warnings"]:
            content.append(f"<p class='error'>{e(note)}</p>")
        content.append(
            "<details open><summary>Actual Meeko output (last 70 nonblank lines)</summary>"
            f"<pre style='max-height:28rem;overflow:auto;white-space:pre-wrap'>"
            f"{e(diagnostics['last_log_lines'])}</pre></details>"
            f"<p><a href='/projects/{e(project.name)}/artifacts/manual_docking/logs/{e('receptor_meeko_Ni_intermediate.log' if job.get('action') == 'receptor_ni_intermediate' else 'receptor_meeko_metal_free.log' if job.get('action') == 'receptor_metal_free' else 'receptor_meeko.log')}'>"
            "Download complete Meeko log</a></p>"
        )
    content.append(
        "<p>Review the PDBQT and retained cofactors. Approval never runs docking.</p>"
        + approve_form(project, 2)
    )
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
            "the approved receptor from Step 2.</b> No ligand re-preparation is necessary.</p>"
        )
    else:
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


def step5(project: Path) -> str:
    job = job_status(project)
    run = active_run(project)
    content = [
        "<h2>Step 5 — Run docking</h2>",
        "<p>The run is unlocked only after Steps 2–4 have been approved for their "
        "current files and settings. Existing completed compound logs are not rerun.</p>",
    ]
    if not approved(project, 4):
        content.append("<p class='error'>Docking is locked until the previous steps are approved.</p>")
    content.append(action(
        project, "dock", "Start / resume AutoDock4 jobs",
        disabled=not approved(project, 4) or job.get("status") == "RUNNING",
    ))
    content.append(job_view(project))
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
    return "".join(content)


def render(project: Path, step: int) -> str:
    if step not in (1,2,3,4,5):
        raise ValueError("Invalid guided step")
    views = {1: step1, 2: step2, 3: step3, 4: step4, 5: step5}
    content = (
        f"<h1>SIDock — Reviewed docking workflow</h1><p><b>Project:</b> {e(project.name)} "
        f"| <a href='/projects/{e(project.name)}'>Project overview</a></p>"
        + flow_nav(project, step) + views[step](project)
    )
    return content + f"<hr><p><a href='/projects/{e(project.name)}/guided/{min(step+1,5)}'>Next step</a></p>"
