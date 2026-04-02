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
    dirs = ["ligand", "zinc", "hits", "docking", "plip", "plif", "admet", "report"]
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
        print(f"[vsdock] Modo sem ligante: etapa 'similarity-search' será pulada automaticamente.")


def cmd_fetch(args):
    from vsdock.fetch import fetch_ligand, fetch_database

    if args.ligand:
        fetch_ligand(args.ligand, outdir="ligand", fmt=args.fmt)
    if args.database:
        fetch_database(
            outdir="zinc",
            source=args.source,
            mw_min=args.mw_min,
            mw_max=args.mw_max,
            logp_min=args.logp_min,
            logp_max=args.logp_max,
            max_mols=args.max_mols,
            available_only=args.available,
            fda_only=args.fda,
            natural_only=args.natural_compounds,
            fragments=args.fragments,
            druglike=not args.fragments,
        )


def cmd_similarity_search(args):
    from pathlib import Path
    import yaml, glob

    state = {}
    if Path("vsdock_state.yaml").exists():
        state = yaml.safe_load(Path("vsdock_state.yaml").read_text())

    if not state.get("has_reference_ligand", True) or not state.get("ligand_smiles"):
        print("[similarity-search] Modo sem ligante de referência detectado.")
        print("[similarity-search] Etapa pulada — banco usado diretamente no docking.")
        return

    from vsdock.screen import screen, load_query_from_state
    query_smiles = load_query_from_state()
    db_file = args.database
    if not db_file:
        candidates = glob.glob("zinc/*.smi") + glob.glob("zinc/*.csv")
        if not candidates:
            print("[similarity-search] ERRO: nenhum banco encontrado em zinc/.")
            raise SystemExit(1)
        db_file = candidates[0]
        print(f"[similarity-search] Banco detectado: {db_file}")

    # Salva threshold no state para o report
    state["screen_threshold"] = args.similarity
    Path("vsdock_state.yaml").write_text(yaml.dump(state))

    screen(
        query_smiles=query_smiles,
        database_file=db_file,
        outdir="hits",
        threshold=args.similarity,
        fp_type=args.fingerprint,
        radius=args.radius,
        max_hits=args.max_hits,
    )


def cmd_clear_library(args):
    from vsdock.clear_library import filter_library
    import os

    hits_file = args.hits_file
    if not hits_file:
        for c in ["hits/hits.csv", "hits/hits_clean.csv"]:
            if os.path.exists(c):
                hits_file = c
                break
    if not hits_file:
        print("[clear-library] ERRO: nenhum arquivo de hits encontrado.")
        raise SystemExit(1)
    print(f"[clear-library] Arquivo detectado: {hits_file}")

    # Se nenhum filtro for especificado, aplica PAINS + Lipinski por padrão
    apply_pains    = args.pains    or (not args.pains and not args.pfizer and not args.gsk)
    apply_lipinski = args.lipinski or (not args.lipinski and not args.pfizer and not args.gsk)

    filter_library(
        hits_file=hits_file,
        outdir="hits",
        apply_pains=apply_pains,
        apply_lipinski=apply_lipinski,
        apply_pfizer=args.pfizer,
        apply_gsk=args.gsk,
    )


def cmd_optimize_geometry(args):
    from vsdock.optimize_geometry import optimize_geometry
    import os

    hits_file = args.hits_file
    if not hits_file:
        for c in ["hits/hits_clean.csv", "hits/hits.csv"]:
            if os.path.exists(c):
                hits_file = c
                break
    if not hits_file:
        print("[optimize] ERRO: nenhum arquivo de hits encontrado.")
        raise SystemExit(1)
    print(f"[optimize] Arquivo detectado: {hits_file}")

    optimize_geometry(
        hits_file=hits_file,
        outdir="hits",
        force_field=args.force_field,
        max_iters=args.max_iters,
        remove_salts=not args.keep_salts,
    )


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
        parts  = line.split()
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
            "[dock] ERRO: defina o sítio com --autobox-ligand, "
            "--autobox-residues, --blind ou --center."
        )
        raise SystemExit(1)

    has_ref   = state.get("has_reference_ligand", True)
    hits_file = args.hits_file

    if not hits_file:
        if has_ref:
            # Prefere geometria otimizada se disponível
            for c in ["hits/hits_optimized.csv", "hits/hits_clean.csv", "hits/hits.csv"]:
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
        hits_file=hits_file,
        receptor=receptor,
        center=center,
        size=size,
        outdir="docking",
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        top_n=args.top,
        reference_smiles=state.get("ligand_smiles"),
        reference_name=state.get("ligand_name", "reference"),
    )


def cmd_plif(args):
    from vsdock.plif import compute_plif
    import os, yaml
    from pathlib import Path

    state = {}
    if Path("vsdock_state.yaml").exists():
        state = yaml.safe_load(Path("vsdock_state.yaml").read_text())

    receptor = args.receptor or state.get("receptor")
    if not receptor:
        print("[plif] ERRO: receptor não definido.")
        raise SystemExit(1)

    docking_file = args.docking_file
    if not docking_file:
        if os.path.exists("docking/docking_results.csv"):
            docking_file = "docking/docking_results.csv"
        else:
            print("[plif] ERRO: nenhum resultado de docking encontrado.")
            raise SystemExit(1)

    compute_plif(
        docking_results=docking_file,
        receptor_pdbqt=receptor,
        poses_dir=args.poses_dir,
        outdir="plif",
        top_n=args.top_n,
        heatmap=args.heatmap,
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

    analyze_plip(
        docking_results=docking_file,
        receptor_pdbqt=receptor,
        poses_dir=args.poses_dir,
        outdir="plip",
        top_n=args.top_n,
    )


def cmd_admet(args):
    from vsdock.admet import predict_admet, list_admet_filters, apply_admet_filters
    import os, yaml
    from pathlib import Path

    # --list-filters: apenas lista os endpoints disponíveis e sai
    if args.list_filters:
        list_admet_filters()
        return

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

    state = {}
    if Path("vsdock_state.yaml").exists():
        state = yaml.safe_load(Path("vsdock_state.yaml").read_text())

    df = predict_admet(
        hits_file=hits_file,
        outdir="admet",
        weights_config=args.weights,
        atc_code=args.atc_code,
    )

    # --apply-filters + --filter-by
    if args.apply_filters and args.filter_by:
        apply_admet_filters(
            df=df,
            filter_by=args.filter_by,
            filters=args.apply_filters,
            outdir="admet",
            reference_name=state.get("ligand_name"),
            thresholds_str=args.thresholds,
        )


def cmd_report(args):
    from vsdock.report import generate_report
    generate_report(outdir="report", project_dir=".", fmt=args.format)


def main():
    parser = argparse.ArgumentParser(
        prog="vsdock",
        description="Virtual Screening and Docking Pipeline",
    )
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
    p.add_argument("--source", choices=["chembl", "pubchem", "zinc"], default="chembl")
    p.add_argument("--mw-min", type=float, default=150, dest="mw_min")
    p.add_argument("--mw-max", type=float, default=500, dest="mw_max")
    p.add_argument("--logp-min", type=float, default=-1, dest="logp_min")
    p.add_argument("--logp-max", type=float, default=5, dest="logp_max")
    p.add_argument("--max-mols", type=int, default=500, dest="max_mols")
    p.add_argument("--available", action="store_true",
                   help="Apenas compostos comercialmente disponíveis")
    p.add_argument("--fda", action="store_true",
                   help="Apenas compostos FDA-aprovados (max_phase=4)")
    p.add_argument("--natural-compounds", action="store_true", dest="natural_compounds",
                   help="Apenas produtos naturais")
    p.add_argument("--fragments", action="store_true",
                   help="Modo fragment-like (Regra de 3: MW≤300, logP≤3)")
    p.add_argument("--druglike", action="store_true", default=True,
                   help="Modo drug-like / Lipinski Ro5 (padrão)")
    p.set_defaults(func=cmd_fetch)

    # similarity-search
    p = subparsers.add_parser("similarity-search", help="Triagem por similaridade estrutural")
    p.add_argument("--similarity", type=float, default=0.4)
    p.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"], default="morgan")
    p.add_argument("--radius", type=int, default=2)
    p.add_argument("--max-hits", type=int, default=500, dest="max_hits")
    p.add_argument("--database", default=None)
    p.set_defaults(func=cmd_similarity_search)

    # clear-library
    p = subparsers.add_parser("clear-library", help="Filtragem da biblioteca (PAINS, Lipinski, Pfizer, GSK)")
    p.add_argument("--hits-file", default=None, dest="hits_file")
    p.add_argument("--pains", action="store_true", help="Aplica filtro PAINS")
    p.add_argument("--lipinski", action="store_true", help="Aplica filtro de Lipinski")
    p.add_argument("--pfizer", action="store_true", help="Aplica filtro de Pfizer")
    p.add_argument("--gsk", action="store_true", help="Aplica filtro de GSK")
    p.set_defaults(func=cmd_clear_library)

    # optimize-library-geometry
    p = subparsers.add_parser("optimize-library-geometry",
                               help="Minimização de energia e otimização de geometria dos ligantes")
    p.add_argument("--hits-file", default=None, dest="hits_file")
    p.add_argument("--force-field", choices=["mmff94", "uff"], default="mmff94", dest="force_field")
    p.add_argument("--max-iters", type=int, default=2000, dest="max_iters")
    p.add_argument("--keep-salts", action="store_true", dest="keep_salts",
                   help="Não remove sais/contrafragmentos do SMILES")
    p.set_defaults(func=cmd_optimize_geometry)

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

    # plif
    p = subparsers.add_parser("plif", help="Protein-Ligand Interaction Fingerprint (ProLIF)")
    p.add_argument("--receptor", default=None)
    p.add_argument("--docking-file", default=None, dest="docking_file")
    p.add_argument("--poses-dir", default="docking/poses", dest="poses_dir")
    p.add_argument("--top-n", type=int, default=10, dest="top_n")
    p.add_argument("--heatmap", action="store_true",
                   help="Gera heatmap PNG das interações (requer display/matplotlib)")
    p.set_defaults(func=cmd_plif)

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
    p.add_argument("--list-filters", action="store_true", dest="list_filters",
                   help="Lista todos os endpoints ADMET disponíveis para filtragem")
    p.add_argument("--apply-filters", nargs="+", default=None, dest="apply_filters",
                   metavar="ENDPOINT",
                   help="Endpoints ADMET para filtrar (ex: hERG DILI AMES)")
    p.add_argument("--filter-by", default=None, dest="filter_by",
                   choices=["reference_compound", "thresholds"],
                   help="Critério de filtragem: reference_compound ou thresholds")
    p.add_argument("--thresholds", default=None,
                   help="Thresholds para filtro (ex: 'hERG<0.5,DILI<0.3')")
    p.set_defaults(func=cmd_admet)

    # report
    p = subparsers.add_parser("report", help="Gera relatório/manuscrito")
    p.add_argument("--format", choices=["quarto", "markdown", "html"], default="quarto")
    p.set_defaults(func=cmd_report)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
