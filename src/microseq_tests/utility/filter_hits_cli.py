# src/microseq_tests/utility/filter_hits_cli.py 

import pandas as pd 
from pathlib import Path
from microseq_tests.utility.id_normaliser import NORMALISERS 

def filter_hits(df: pd.DataFrame, 
                ident: int, qcov: int, 
                with_status: bool = False):
    """Return (filtered_df, status_table_or_None)."""
    mask = (df.pident >= ident) & (df.qcovhsp >= qcov) 
    filtered = df.loc[mask] 

    if with_status:
        need_id = df.pident < ident 
        need_cov = df.qcovhsp < qcov 
        status = df.assign(
            status = mask.map({True: "PASS", False: "FAIL"}),
            need_id = need_id,
            need_cov = need_cov,
        )
        return filtered, status

    return filtered, None 


# --- CLI entry-point ----------- 

def main(args) -> None:
    inp = Path(args.input)
    df  = pd.read_csv(inp, sep="\t")

    out_path = (Path(args.output) if args.output
                else inp.parent / f"hits_{args.identity}_{args.qcov}.tsv")

    status_path = inp.parent / "hits_full_status.tsv"

    filt, status = filter_hits(df,
                               ident=args.identity,
                               qcov=args.qcov,
                               with_status=args.with_status)
    
    filt.to_csv(out_path, sep="\t", index=False)

    if args.unique:
        ids = filt[args.group_col].map(NORMALISERS[args.id_normaliser])
        unique = ids.nunique()
        print(f"{len(filt):>5} PASS rows  "
          f"({unique} unique {args.group_col} "
          f"after {args.id_normaliser}) → {out_path}")
    else:
        print(f"{len(filt):>5} PASS rows → {out_path}")
    
    if status is not None:
        status.to_csv(status_path, sep="\t", index=False)
        print(f"      status table          → {status_path}") 




                

