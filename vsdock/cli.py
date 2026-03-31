"""
vsdock CLI — entry point principal
"""

import argparse
from vsdock import __version__


def cmd_init(args):
    from vsdock.fetch import fetch_ligand_smiles
    import yaml
    from pathlib import Path

    print(f"[vsdock] Inicializando projeto '{args.ligand}' ...")

    dirs = ["ligand", "zinc", "hits", "docking", "plip", "admet", "report"]
    for d in dirs:
        Path(d).mkdir(exist_ok=True)

    smiles = fetch_ligand_smiles(args.ligand)
    state = {
        "ligand_name": args.ligand,
        "ligand_smiles": smiles,
        "receptor": args.target,
        "config": args.config,
        "box_center": [0.0, 0.0, 0.0],
        "box_size":   [20.0, 20.0, 20.0],
    }
    Path("vsdock_state.yaml").write_text(yaml.dump(state))
    print(f"[vsdock] Projeto inicializado. Estado salvo em vsdock_state.yaml")
    print(f"[vsdock] IMPORTANTE: edite vsdock_state.yaml e defina box_center e box_size antes de rodar o docking.")


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
        )


def cmd_screen(args):
    from vsdock.screen import screen, load_query_from_state
    import glob

    query_smiles = load_query_from_state()

    db_file = args.database
    if not db_file:
        candidates = glob.glob("zinc/*.smi") + glob.glob("zinc/*.csv")
        if not candidates:
            print("[screen] ERRO: nenhum banco encontrado em zinc/. Rode 'vsdock fetch --database' primeiro.")
            raise SystemExit(1)
        db_file = candidates[0]
        print(f"[screen] Banco detectado: {db_file}")

    screen(
        query_smiles=query_smiles,
        database_file=db_file,
        outdir="hits",
        threshold=args.similarity,
        fp_type=args.fingerprint,
        radius=args.radius,
        max_hits=args.max_hits,
    )


def cmd_analyze(args):
    from vsdock.pains import filter_pains
    import os

    hits_file = args.hits_file
    if not hits_file:
        candidates = ["hits/hits.csv", "hits/hits_clean.csv"]
        for c in candidates:
            if os.path.exists(c):
                hits_file = c
                break
        if not hits_file:
            print("[analyze] ERRO: nenhum arquivo de hits encontrado. Rode 'vsdock screen' primeiro.")
            raise SystemExit(1)
        print(f"[analyze] Arquivo detectado: {hits_file}")

    filter_pains(
        hits_file=hits_file,
        outdir="hits",
        also_filter_lipinski=not args.no_lipinski,
    )


def cmd_dock(args):
    from vsdock.dock import dock_all, get_box_from_autobox
    import yaml, os
    from pathlib import Path

    state = {}
    if Path("vsdock_state.yaml").exists():
        state = yaml.safe_load(Path("vsdock_state.yaml").read_text())

    receptor = args.receptor or state.get("receptor")
    if not receptor:
        print("[dock] ERRO: receptor não definido. Use --receptor ou defina em vsdock_state.yaml")
        raise SystemExit(1)

    # Coordenadas da caixa: autobox > argumento > state > erro
    if args.autobox_ligand or args.autobox_residues or args.blind:
        pdb = args.pdb or state.get("pdb")
        if not pdb:
            print("[dock] ERRO: --autobox requer --pdb com o arquivo PDB original (com ligante).")
            raise SystemExit(1)
        box = get_box_from_autobox(
            pdb_file=pdb,
            ligand=args.autobox_ligand,
            residues=args.autobox_residues,
            padding=args.padding,
            blind=args.blind,
        )
        center = box["center"]
        size   = box["size"]

        # Salva no state para reusar
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
        print("[dock] ERRO: defina a caixa com --autobox-ligand, --center, ou edite vsdock_state.yaml")
        raise SystemExit(1)

    hits_file = args.hits_file
    if not hits_file:
        candidates = ["hits/hits_clean.csv", "hits/hits.csv"]
        for c in candidates:
            if os.path.exists(c):
                hits_file = c
                break
        if not hits_file:
            print("[dock] ERRO: nenhum arquivo de hits encontrado.")
            raise SystemExit(1)
        print(f"[dock] Hits detectados: {hits_file}")

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


def cmd_report(args):
    print(f"[vsdock] Gerando relatorio (formato: {args.format})")
    print("  -> modulo report ainda nao implementado")


def main():
    parser = argparse.ArgumentParser(
        prog="vsdock",
        description="Virtual Screening and Docking Pipeline",
    )
    parser.add_argument("--version", action="version", version=f"vsdock {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    # --- init ---
    p_init = subparsers.add_parser("init", help="Inicializa projeto e baixa ligante")
    p_init.add_argument("--target", required=True, help="Arquivo .pdbqt do receptor")
    p_init.add_argument("--ligand", required=True, help="Nome do ligante (PubChem)")
    p_init.add_argument("--config", default="configs/default.yaml")
    p_init.set_defaults(func=cmd_init)

    # --- fetch ---
    p_fetch = subparsers.add_parser("fetch", help="Baixa ligante e/ou banco molecular")
    p_fetch.add_argument("--ligand", default=None)
    p_fetch.add_argument("--fmt", choices=["sdf", "smiles"], default="sdf")
    p_fetch.add_argument("--database", action="store_true")
    p_fetch.add_argument("--source", choices=["chembl", "zinc"], default="chembl")
    p_fetch.add_argument("--mw-min", type=float, default=150, dest="mw_min")
    p_fetch.add_argument("--mw-max", type=float, default=500, dest="mw_max")
    p_fetch.add_argument("--logp-min", type=float, default=-1, dest="logp_min")
    p_fetch.add_argument("--logp-max", type=float, default=5, dest="logp_max")
    p_fetch.add_argument("--max-mols", type=int, default=500, dest="max_mols")
    p_fetch.set_defaults(func=cmd_fetch)

    # --- screen ---
    p_screen = subparsers.add_parser("screen", help="Triagem por similaridade estrutural")
    p_screen.add_argument("--similarity", type=float, default=0.4)
    p_screen.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"], default="morgan")
    p_screen.add_argument("--radius", type=int, default=2)
    p_screen.add_argument("--max-hits", type=int, default=500, dest="max_hits")
    p_screen.add_argument("--database", default=None)
    p_screen.set_defaults(func=cmd_screen)

    # --- analyze ---
    p_analyze = subparsers.add_parser("analyze", help="Filtro PAINS + Lipinski")
    p_analyze.add_argument("--hits-file", default=None, dest="hits_file")
    p_analyze.add_argument("--no-lipinski", action="store_true", dest="no_lipinski")
    p_analyze.set_defaults(func=cmd_analyze)

    # --- dock ---
    p_dock = subparsers.add_parser("dock", help="Docking com AutoDock Vina")
    p_dock.add_argument("--receptor", default=None, help="Arquivo .pdbqt do receptor")
    p_dock.add_argument("--pdb", default=None, help="Arquivo PDB original (para autobox)")
    # Autobox
    p_dock.add_argument("--autobox-ligand", default=None, dest="autobox_ligand",
                        help="ID do ligante cocristalizado para autobox (ex: MRV1101A)")
    p_dock.add_argument("--autobox-residues", nargs="+", default=None, dest="autobox_residues",
                        help="Resíduos para definir a caixa (ex: VAL2A LYS4A)")
    p_dock.add_argument("--blind", action="store_true",
                        help="Docking cego (caixa = proteína toda)")
    p_dock.add_argument("--padding", type=float, default=5.0,
                        help="Padding do autobox em Angstroms (default: 5.0)")
    # Caixa manual (alternativa ao autobox)
    p_dock.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help="Centro da caixa (alternativa ao autobox)")
    p_dock.add_argument("--size", nargs=3, type=float, metavar=("SX", "SY", "SZ"),
                        default=None, help="Tamanho da caixa em Å (default: 20 20 20)")
    # Vina
    p_dock.add_argument("--exhaustiveness", type=int, default=8)
    p_dock.add_argument("--num-modes", type=int, default=9, dest="num_modes")
    p_dock.add_argument("--top", type=int, default=None)
    p_dock.add_argument("--hits-file", default=None, dest="hits_file")
    p_dock.set_defaults(func=cmd_dock)

    # --- report ---
    p_report = subparsers.add_parser("report", help="Gera relatorio/manuscrito")
    p_report.add_argument("--format", choices=["quarto", "markdown", "html"], default="quarto")
    p_report.set_defaults(func=cmd_report)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
