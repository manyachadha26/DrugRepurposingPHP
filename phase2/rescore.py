#!/usr/bin/env python3
"""
Re-score poses we already docked, using Vina's Vinardo scoring function.

WHY THIS SCRIPT EXISTS
----------------------
Gate 2 told us Vina's default scoring function cannot rank PDE4D binders
(AUC 0.612, barely better than a coin flip). But Gate 1 told us the POSES
are right (0.77 A redocking). So we do not need to dock anything again.
We keep the poses and swap the scoring function.

Vinardo is a different scoring function built into AutoDock Vina 1.2.x.
It weights hydrophobic contact and hydrogen bonding differently and drops
Vina's gauss1 term. On many targets it ranks better than Vina. It costs us
nothing to test because scoring a pose takes ~0.1 s, versus ~45 s to dock it.

This script also pulls the ORIGINAL Vina score straight out of each pose
file's REMARK line, so one CSV holds both numbers and we can compare them
on exactly the same molecules.

USAGE
-----
    python rescore.py <receptor.pdbqt> <docked_folder> <actives_folder> [n_jobs] [config.txt]

Example:
    python rescore.py 1XOQ.pdbqt docked2 actives2 8 config_1XOQ.txt

OUTPUT
------
    rescore_results.csv   name, label, vina_score, vinardo_score

"label" is 1 for a known PDE4 active, 0 for a decoy. It is read from the
actives folder: any ligand whose name appears there is an active.
"""

import os
import re
import sys
import glob
import shutil
import subprocess
import tempfile
from multiprocessing import Pool

VINA = shutil.which("vina")

# matches: REMARK VINA RESULT:    -11.229      0.000      0.000
RE_VINA_REMARK = re.compile(
    r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+\.\d+)", re.MULTILINE
)
# matches Vina's --score_only stdout line
RE_SCORE_ONLY = re.compile(
    r"Estimated Free Energy of Binding\s*:\s*(-?\d+\.\d+)"
)


def first_model(path):
    """
    Return the text of MODEL 1 only.

    A docked output file holds up to 9 poses stacked as MODEL/ENDMDL blocks.
    We only ever care about the top-ranked pose, because that is the one a
    real screen would act on. Handing vina a multi-model file is ambiguous,
    so we slice out the first block ourselves.
    """
    out, inside = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                if inside:          # second MODEL reached - stop
                    break
                inside = True
                continue
            if line.startswith("ENDMDL"):
                break
            out.append(line)
    if not out:                      # no MODEL records: single-pose file
        with open(path) as fh:
            return fh.read()
    return "".join(out)


def original_vina_score(path):
    """Read the score Vina already wrote into the pose file. No recompute."""
    with open(path) as fh:
        m = RE_VINA_REMARK.search(fh.read())
    return float(m.group(1)) if m else None


_shown = {"err": False}


def score_one(args):
    receptor, pose_path, config = args
    name = os.path.basename(pose_path).replace(".pdbqt", "")
    vina_score = original_vina_score(pose_path)

    tmp = None
    try:
        # write MODEL 1 to its own file for vina --score_only
        fd, tmp = tempfile.mkstemp(suffix=".pdbqt")
        with os.fdopen(fd, "w") as fh:
            fh.write(first_model(pose_path))

        cmd = [VINA,
               "--receptor", receptor,
               "--ligand", tmp,
               "--scoring", "vinardo",
               "--score_only"]
        # Vina 1.2.x still wants a search space in --score_only mode.
        # Passing the SAME config used for docking keeps the box identical,
        # so the two scoring functions are compared on equal terms.
        if config:
            cmd += ["--config", config]
        else:
            cmd += ["--autobox"]

        res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        m = RE_SCORE_ONLY.search(res.stdout)
        vinardo = float(m.group(1)) if m else None
        if vinardo is None and not _shown["err"]:
            _shown["err"] = True
            sys.stderr.write(
                f"\n[FIRST FAILURE] {name}\n"
                f"  command : {' '.join(cmd)}\n"
                f"  retcode : {res.returncode}\n"
                f"  stdout  : {res.stdout.strip()[-600:]}\n"
                f"  stderr  : {res.stderr.strip()[-600:]}\n\n")
    except subprocess.TimeoutExpired:
        sys.stderr.write(f"[warn] {name}: timed out\n")
        vinardo = None
    except Exception as exc:
        sys.stderr.write(f"[warn] {name}: {exc}\n")
        vinardo = None
    finally:
        if tmp and os.path.exists(tmp):
            os.remove(tmp)

    return name, vina_score, vinardo


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    receptor = sys.argv[1]
    docked_dir = sys.argv[2]
    actives_dir = sys.argv[3]
    n_jobs = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    config = sys.argv[5] if len(sys.argv) > 5 else None

    if VINA is None:
        sys.exit("ERROR: 'vina' not found on PATH. Activate your conda env.")
    for p in (receptor, docked_dir, actives_dir):
        if not os.path.exists(p):
            sys.exit(f"ERROR: not found: {p}")

    poses = sorted(glob.glob(os.path.join(docked_dir, "*.pdbqt")))
    if not poses:
        sys.exit(f"ERROR: no .pdbqt files in {docked_dir}")

    # Which names are actives? Read the folder, never a hardcoded list.
    active_names = {
        os.path.basename(f).replace(".pdbqt", "").replace(".smi", "")
        for f in glob.glob(os.path.join(actives_dir, "*"))
    }
    print(f"Receptor      : {receptor}")
    print(f"Box           : {config if config else 'auto (--autobox)'}")
    print(f"Poses to score: {len(poses)}")
    print(f"Actives known : {len(active_names)}")
    print(f"Parallel jobs : {n_jobs}")
    print("Scoring with Vinardo. This is fast - a few minutes, not hours.\n")

    rows = []
    with Pool(n_jobs) as pool:
        for i, (name, vina_s, vinardo_s) in enumerate(
            pool.imap_unordered(
                score_one, [(receptor, p, config) for p in poses]), 1
        ):
            rows.append((name, 1 if name in active_names else 0,
                         vina_s, vinardo_s))
            if i % 100 == 0 or i == len(poses):
                print(f"  {i} / {len(poses)}")

    import pandas as pd
    df = pd.DataFrame(rows, columns=["name", "label",
                                     "vina_score", "vinardo_score"])
    df = df.sort_values("vinardo_score")
    df.to_csv("rescore_results.csv", index=False)

    n_act = int(df.label.sum())
    missing = int(df.vinardo_score.isna().sum())
    print(f"\nWrote rescore_results.csv")
    print(f"  rows {len(df)}   actives {n_act}   decoys {len(df) - n_act}")
    if missing:
        print(f"  WARNING: {missing} ligands failed to score "
              f"(they stay in the file with a blank score)")
    print("\nNext:  python enrichment2.py rescore_results.csv")


if __name__ == "__main__":
    main()
