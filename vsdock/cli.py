"""
vsdock CLI — entry point principal
"""

import argparse
from vsdock import __version__


def cmd_prepare(args):
    from vsdock.prepare import prepare_receptor
    prepare_receptor(
        ligand_code=args.ligand_code,
        pdbid=args.pdbid,
        pdb_file=args.pdb_file,
        outdir=".",
    )


def cmd_init(args):
    import yaml
    from pathlib import Path

    ligand = args.ligand
    dirs = ["ligand", "zinc", "hits", "docking", "plip", "admet", "report"]
    for d in dirs:
        Path(d).mkdir(exist_ok=True)

    state = {
        "receptor": args.target,
        "config": args.config,
        "box_center": [0.0, 0.0, 0.0],
        "box_size":   [20.0, 20.0, 20.0],
        "has_reference_ligand": ligand is not None,
    }

    if ligand:
        from vsdock.fetch import fetch_ligand_smiles
        print(f"[vsdock] Inicializando projeto com ligante '{ligand}' ...")
        smiles = fetch_ligand_smiles(ligand)
        state["ligand_name"]   = ligand
        state["ligand_smiles"] = smiles
    else:
        print(f"[vsdock] Inicializando projeto sem ligante de referência (receptor-only mode) ...")
        state["ligand_name"]   = None
        state["ligand_smiles"] = None

    Path("vsdock_state.yaml").write_text(yaml.dump(state))
    print(f"[vsdock] Projeto inicializado. Estado salvo em vsdock_state.yaml")

    if not ligand:
        print(f"[vsdock] Modo sem ligante: etapa 'screen' será pulada automaticamente.")
        print(f"[vsdock] Use 'vsdock dock --blind' ou 'vsdock dock --autobox-residues RES1 RES2'")


def cmd_fetch(args):
    from vsdock.fetch import fetch_ligand, fetch_database

    if args.ligand:
        fetch_ligand(args.ligand, outdir="ligand", fmt=args.fmt)
    if args.database:
        fetch_database(
            outdir="zinc", source=args.source,
            mw_min=args.mw_min, mw_max=args.mw_max,
            logp_min=args.logp_min, logp_max=args.logp_max,
            max_mols=args.max_mols,
            available_only=args.available,
        )


def cmd_screen(args):
    from pathlib import Path
    import yaml, glob

    state = {}
    if Path("vsdock_state.yaml").exists():
        state = yaml.safe_load(Path("vsdock_state.yaml").read_text())

    if not state.get("has_reference_ligand", True) or not state.get("ligand_smiles"):
        print("[screen] Modo sem ligante de referência detectado.")
        print("[screen] Etapa de screening pulada — o banco será usado diretamente no docking.")
        return

    from vsdock.screen import screen, load_query_from_state
    query_smiles = load_query_from_state()
    db_file = args.database
    if not db_file:
        candidates = glob.glob("zinc/*.smi") + glob.glob("zinc/*.csv")
        if not candidates:
            print("[screen] ERRO: nenhum banco encontrado em zinc/.")
            raise SystemExit(1)
        db_file = candidates[0]
        print(f"[screen] Banco detectado: {db_file}")

    screen(
        query_smiles=query_smiles, database_file=db_file,
        outdir="hits", threshold=args.similarity,
        fp_type=args.fingerprint, radius=args.radius, max_hits=args.max_hits,
    )


def cmd_analyze(args):
    from vsdock.pains import filter_pains
    import os

    hits_file = args.hits_file
    if not hits_file:
        for c in ["hits/hits.csv", "hits/hits_clean.csv"]:
            if os.path.exists(c):
                hits_file = c
                break
    if not hits_file:
        print("[analyze] ERRO: nenhum arquivo de hits encontrado.")
        raise SystemExit(1)
    print(f"[analyze] Arquivo detectado: {hits_file}")
    filter_pains(hits_file=hits_file, outdir="hits", also_filter_lipinski=not args.no_lipinski)


def _bank_to_hits_csv() -> str:
    import glob
    import pandas as pd
    from pathlib import Path

    candidates = glob.glob("zinc/*.smi") + glob.glob("zinc/*.csv")
    if not candidates:
        return None

    src = candidates[0]
    out = "hits/hits_from_bank.csv"
    Path("hits").mkdir(exist_ok=True)

    if src.endswith(".csv"):
        return src

    lines = [l.strip() for l in open(src) if l.strip()]
    rows = []
    for line in lines:
        parts = line.split()
        smiles = parts[0]
        mol_id = parts[1] if len(parts) > 1 else f"mol_{len(rows)}"
        rows.append({"smiles": smiles, "id": mol_id, "tanimoto": None})

    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"[dock] Banco convertido: {len(rows)} moléculas → {out}")
    return out


def cmd_dock(args):
    from vsdock.dock import dock_all, get_box_from_autobox
    import yaml, os
    from pathlib import Path

    state = {}
    if Path("vsdock_state.yaml").exists():
        state = yaml.safe_load(Path("vsdock_state.yaml").read_text())

    receptor = args.receptor or state.get("receptor")
    if not receptor:
        print("[dock] ERRO: receptor não definido.")
        raise SystemExit(1)

    if args.autobox_ligand or args.autobox_residues or args.blind:
        pdb = args.pdb or state.get("pdb")
        if not pdb and not args.blind:
            print("[dock] ERRO: --autobox requer --pdb com o arquivo PDB original.")
            raise SystemExit(1)
        box = get_box_from_autobox(
            pdb_file=pdb,
            ligand=args.autobox_ligand,
            residues=args.autobox_residues,
            padding=args.padding,
            blind=args.blind,
        )
        center, size = box["center"], box["size"]
        state["box_center"] = list(center)
        state["box_size"]   = list(size)
        Path("vsdock_state.yaml").write_text(yaml.dump(state))
    elif args.center:
        center = tuple(args.center)
        size   = tuple(args.size) if args.size else tuple(state.get("box_size", [20, 20, 20]))
    elif state.get("box_center") and any(v != 0 for v in state["box_center"]):
        center = tuple(state["box_center"])
        size   = tuple(state.get("box_size", [20, 20, 20]))
    else:
        print(
            "[dock] ERRO: defina o sítio de docking com uma das opções:\n"
            "  --autobox-ligand MRV1101A\n"
            "  --autobox-residues VAL2A LYS4A\n"
            "  --blind\n"
            "  --center X Y Z"
        )
        raise SystemExit(1)

    has_ref   = state.get("has_reference_ligand", True)
    hits_file = args.hits_file

    if not hits_file:
        if has_ref:
            for c in ["hits/hits_clean.csv", "hits/hits.csv"]:
                if os.path.exists(c):
                    hits_file = c
                    break
        else:
            hits_file = _bank_to_hits_csv()

    if not hits_file:
        print("[dock] ERRO: nenhum arquivo de moléculas encontrado.")
        raise SystemExit(1)
    print(f"[dock] Moléculas detectadas: {hits_file}")

    dock_all(
        hits_file=hits_file, receptor=receptor, center=center, size=size,
        outdir="docking", exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes, top_n=args.top,
        reference_smiles=state.get("ligand_smiles"),
        reference_name=state.get("ligand_name", "reference"),
    )


def cmd_plip(args):
    from vsdock.plip import analyze_plip
    import os, yaml
    from pathlib import Path

    state = {}
    if Path("vsdock_state.yaml").exists():
        state = yaml.safe_load(Path("vsdock_state.yaml").read_text())

    receptor = args.receptor or state.get("receptor")
    if not receptor:
        print("[plip] ERRO: receptor não definido.")
        raise SystemExit(1)

    docking_file = args.docking_file
    if not docking_file:
        if os.path.exists("docking/docking_results.csv"):
            docking_file = "docking/docking_results.csv"
        else:
            print("[plip] ERRO: nenhum resultado de docking encontrado.")
            raise SystemExit(1)
    print(f"[plip] Resultados de docking: {docking_file}")

    analyze_plip(
        docking_results=docking_file, receptor_pdbqt=receptor,
        poses_dir=args.poses_dir, outdir="plip", top_n=args.top_n,
    )


def cmd_admet(args):
    from vsdock.admet import predict_admet
    import os

    hits_file = args.hits_file
    if not hits_file:
        for c in ["docking/docking_results.csv", "hits/hits_clean.csv",
                  "hits/hits_from_bank.csv", "hits/hits.csv"]:
            if os.path.exists(c):
                hits_file = c
                break
    if not hits_file:
        print("[admet] ERRO: nenhum arquivo de hits encontrado.")
        raise SystemExit(1)
    print(f"[admet] Arquivo detectado: {hits_file}")

    predict_admet(
        hits_file=hits_file, outdir="admet",
        weights_config=args.weights, atc_code=args.atc_code,
    )


def cmd_report(args):
    from vsdock.report import generate_report
    generate_report(outdir="report", project_dir=".", fmt=args.format)


def main():
    parser = argparse.ArgumentParser(prog="vsdock", description="Virtual Screening and Docking Pipeline")
    parser.add_argument("--version", action="version", version=f"vsdock {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    # prepare
    p = subparsers.add_parser("prepare", help="Baixa e prepara receptor a partir de PDB ID ou arquivo")
    p.add_argument("--pdbid", default=None)
    p.add_argument("--pdb-file", default=None, dest="pdb_file")
    p.add_argument("--ligand-code", default=None, dest="ligand_code")
    p.set_defaults(func=cmd_prepare)

    # init
    p = subparsers.add_parser("init", help="Inicializa projeto")
    p.add_argument("--target", required=True)
    p.add_argument("--ligand", default=None)
    p.add_argument("--config", default="configs/default.yaml")
    p.set_defaults(func=cmd_init)

    # fetch
    p = subparsers.add_parser("fetch", help="Baixa ligante e/ou banco molecular")
    p.add_argument("--ligand", default=None)
    p.add_argument("--fmt", choices=["sdf", "smiles"], default="sdf")
    p.add_argument("--database", action="store_true")
    p.add_argument("--source", choices=["chembl", "zinc"], default="chembl")
    p.add_argument("--mw-min", type=float, default=150, dest="mw_min")
    p.add_argument("--mw-max", type=float, default=500, dest="mw_max")
    p.add_argument("--logp-min", type=float, default=-1, dest="logp_min")
    p.add_argument("--logp-max", type=float, default=5, dest="logp_max")
    p.add_argument("--max-mols", type=int, default=500, dest="max_mols")
    p.add_argument("--available", action="store_true",
                   help="Filtra apenas compostos comercialmente disponíveis (ChEMBL)")
    p.set_defaults(func=cmd_fetch)

    # screen
    p = subparsers.add_parser("screen", help="Triagem por similaridade estrutural")
    p.add_argument("--similarity", type=float, default=0.4)
    p.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"], default="morgan")
    p.add_argument("--radius", type=int, default=2)
    p.add_argument("--max-hits", type=int, default=500, dest="max_hits")
    p.add_argument("--database", default=None)
    p.set_defaults(func=cmd_screen)

    # analyze
    p = subparsers.add_parser("analyze", help="Filtro PAINS + Lipinski")
    p.add_argument("--hits-file", default=None, dest="hits_file")
    p.add_argument("--no-lipinski", action="store_true", dest="no_lipinski")
    p.set_defaults(func=cmd_analyze)

    # dock
    p = subparsers.add_parser("dock", help="Docking com AutoDock Vina")
    p.add_argument("--receptor", default=None)
    p.add_argument("--pdb", default=None)
    p.add_argument("--autobox-ligand", default=None, dest="autobox_ligand")
    p.add_argument("--autobox-residues", nargs="+", default=None, dest="autobox_residues")
    p.add_argument("--blind", action="store_true")
    p.add_argument("--padding", type=float, default=5.0)
    p.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    p.add_argument("--size", nargs=3, type=float, metavar=("SX", "SY", "SZ"), default=None)
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--num-modes", type=int, default=9, dest="num_modes")
    p.add_argument("--top", type=int, default=None)
    p.add_argument("--hits-file", default=None, dest="hits_file")
    p.set_defaults(func=cmd_dock)

    # plip
    p = subparsers.add_parser("plip", help="Análise de interações proteína-ligante")
    p.add_argument("--receptor", default=None)
    p.add_argument("--docking-file", default=None, dest="docking_file")
    p.add_argument("--poses-dir", default="docking/poses", dest="poses_dir")
    p.add_argument("--top-n", type=int, default=10, dest="top_n")
    p.set_defaults(func=cmd_plip)

    # admet
    p = subparsers.add_parser("admet", help="Predição de propriedades ADMET")
    p.add_argument("--hits-file", default=None, dest="hits_file")
    p.add_argument("--weights", default=None)
    p.add_argument("--atc-code", default=None, dest="atc_code")
    p.set_defaults(func=cmd_admet)

    # report
    p = subparsers.add_parser("report", help="Gera relatorio/manuscrito")
    p.add_argument("--format", choices=["quarto", "markdown", "html"], default="quarto")
    p.set_defaults(func=cmd_report)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
