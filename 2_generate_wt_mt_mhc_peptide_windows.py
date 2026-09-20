#!/usr/bin/env python3
"""
Default alleles:
    H2-Kb -> 8 aa,  class I
    H2-Db -> 9 aa,  class I
    I-Ab  -> 15 aa, class II
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


DEFAULT_ALLELES = ("H2-Kb:8:I", "H2-Db:9:I", "I-Ab:15:II")

CLASS_II_NAME_HINTS = ("I-A", "I-E", "DR", "DQ", "DP")


SOURCE_COLUMN_RENAMES = {
    "WT_peptide": "Original_WT_window",
    "MUT_peptide": "Original_MT_window",
    "Peptide_start": "Original_Peptide_start",
    "Peptide_end": "Original_Peptide_end",
    "Mutation_position_in_peptide": "Original_Mutation_position_in_peptide",
    "Expected_peptide_length": "Original_Expected_peptide_length",
}


def open_text(path: str):
    """Open a plain-text or gzip-compressed file for reading."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def remove_version(identifier: str) -> str:
    """NP_001001144.2 -> NP_001001144."""
    return identifier.rsplit(".", 1)[0]


def guess_mhc_class(allele: str) -> str:
    allele_upper = allele.upper()
    for hint in CLASS_II_NAME_HINTS:
        if hint in allele_upper:
            return "II"
    return "I"


def parse_alleles(values: Iterable[str]) -> list[tuple[str, int, str]]:
    parsed: list[tuple[str, int, str]] = []
    seen: set[str] = set()

    for value in values:
        parts = value.rsplit(":", 2)

        mhc_class = None
        if len(parts) == 3 and parts[-1].upper() in {"I", "II"}:
            allele, length_text, mhc_class = parts
            mhc_class = mhc_class.upper()
        elif len(parts) >= 2:
            allele, length_text = value.rsplit(":", 1)
        else:
            raise ValueError(
                f"Invalid allele specification {value!r}; "
                "expected ALLELE:LENGTH or ALLELE:LENGTH:CLASS"
            )

        try:
            length = int(length_text)
        except ValueError as exc:
            raise ValueError(
                f"Invalid allele specification {value!r}; "
                "LENGTH must be an integer"
            ) from exc

        allele = allele.strip()
        if mhc_class is None:
            mhc_class = guess_mhc_class(allele)

        seen.add(allele)
        parsed.append((allele, length, mhc_class))

    return parsed


def load_proteins(fasta_path: str):
    exact: dict[str, str] = {}
    versionless_candidates: dict[str, list[tuple[str, str]]] = {}

    with open_text(fasta_path) as handle:
        protein_id = None
        sequence_parts: list[str] = []

        def store_current_record():
            if protein_id is None:
                return
            sequence = "".join(sequence_parts).replace(" ", "").rstrip("*").upper()
            exact[protein_id] = sequence
            base_id = remove_version(protein_id)
            versionless_candidates.setdefault(base_id, []).append(
                (protein_id, sequence)
            )

        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store_current_record()
                protein_id = line[1:].split(maxsplit=1)[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

        store_current_record()


    versionless: dict[str, tuple[str, str]] = {}
    ambiguous: set[str] = set()
    for base_id, candidates in versionless_candidates.items():
        unique_sequences = {sequence for _, sequence in candidates}
        if len(unique_sequences) == 1:
            versionless[base_id] = candidates[0]
        else:
            ambiguous.add(base_id)

    return exact, versionless, ambiguous


def resolve_protein(
    protein_id: str,
    exact: dict[str, str],
    versionless: dict[str, tuple[str, str]],
    ambiguous: set[str],
):
    if protein_id in exact:
        return protein_id, exact[protein_id], None

    base_id = remove_version(protein_id)
    if base_id in ambiguous:
        return None, None, "ambiguous_versionless_protein_id"
    if base_id in versionless:
        resolved_id, sequence = versionless[base_id]
        return resolved_id, sequence, None

    return None, None, "protein_id_not_found"


def parse_protein_position(value) -> int | None:
    if pd.isna(value):
        return None

    text = str(value).strip()
    match = re.fullmatch(r"(\d+)(?:/\d+)?", text)
    if match is None:
        return None
    return int(match.group(1))


def parse_single_amino_acid(value) -> str | None:
    if pd.isna(value):
        return None
    amino_acid = str(value).strip().upper()
    if len(amino_acid) != 1 or amino_acid not in STANDARD_AMINO_ACIDS:
        return None
    return amino_acid


def read_table(path: str, separator: str | None) -> pd.DataFrame:
    if separator is None:
        separator = "," if str(path).lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=separator, dtype=str, keep_default_na=False)


def make_skip_record(row: dict, source_row: int, reason: str, detail: str, allele=""):
    record = dict(row)
    record.update(
        {
            "Source_row": source_row,
            "Allele": allele,
            "Skip_reason": reason,
            "Skip_detail": detail,
        }
    )
    return record


def generate_peptides(args) -> tuple[pd.DataFrame, pd.DataFrame]:
    allele_lengths = parse_alleles(args.alleles)
    input_df = read_table(args.input_table, args.separator)

    required = {
        args.protein_column,
        args.position_column,
        args.reference_aa_column,
        args.mutant_aa_column,
    }
    missing = sorted(required.difference(input_df.columns))


    input_df = input_df.rename(
        columns={
            old: new
            for old, new in SOURCE_COLUMN_RENAMES.items()
            if old in input_df.columns
        }
    )

    exact, versionless, ambiguous = load_proteins(args.protein_fasta)
    print(f"Loaded {len(exact):,} protein FASTA records", file=sys.stderr)

    output_records: list[dict] = []
    skipped_records: list[dict] = []

    for row_number, (_, series) in enumerate(input_df.iterrows(), start=1):
        row = series.to_dict()
        protein_id = str(row.get(args.protein_column, "")).strip()
        protein_position = parse_protein_position(row.get(args.position_column))
        reference_aa = parse_single_amino_acid(row.get(args.reference_aa_column))
        mutant_aa = parse_single_amino_acid(row.get(args.mutant_aa_column))

        if not protein_id:
            skipped_records.append(
                make_skip_record(row, row_number, "missing_protein_id", "")
            )
            continue
        if protein_position is None:
            skipped_records.append(
                make_skip_record(
                    row,
                    row_number,
                    "invalid_protein_position",
                    str(row.get(args.position_column, "")),
                )
            )
            continue
        if reference_aa is None:
            skipped_records.append(
                make_skip_record(
                    row,
                    row_number,
                    "invalid_reference_amino_acid",
                    str(row.get(args.reference_aa_column, "")),
                )
            )
            continue
        if mutant_aa is None:
            skipped_records.append(
                make_skip_record(
                    row,
                    row_number,
                    "invalid_mutant_amino_acid",
                    str(row.get(args.mutant_aa_column, "")),
                )
            )
            continue
        if reference_aa == mutant_aa:
            skipped_records.append(
                make_skip_record(
                    row,
                    row_number,
                    "not_an_amino_acid_substitution",
                    f"{reference_aa}>{mutant_aa}",
                )
            )
            continue

        resolved_id, protein_sequence, resolution_error = resolve_protein(
            protein_id, exact, versionless, ambiguous
        )
        if resolution_error:
            skipped_records.append(
                make_skip_record(
                    row, row_number, resolution_error, protein_id
                )
            )
            continue

        position0 = protein_position - 1
        if position0 < 0 or position0 >= len(protein_sequence):
            skipped_records.append(
                make_skip_record(
                    row,
                    row_number,
                    "protein_position_out_of_bounds",
                    f"position={protein_position}; protein_length={len(protein_sequence)}",
                )
            )
            continue

        fasta_reference_aa = protein_sequence[position0]
        if fasta_reference_aa != reference_aa:
            skipped_records.append(
                make_skip_record(
                    row,
                    row_number,
                    "reference_amino_acid_mismatch",
                    f"table={reference_aa}; FASTA={fasta_reference_aa}",
                )
            )
            continue

        for allele, peptide_length, mhc_class in allele_lengths:
            if len(protein_sequence) < peptide_length:
                skipped_records.append(
                    make_skip_record(
                        row,
                        row_number,
                        "protein_shorter_than_target_peptide",
                        f"protein_length={len(protein_sequence)}; target={peptide_length}",
                        allele,
                    )
                )
                continue

            start_min = max(0, position0 - peptide_length + 1)
            start_max = min(position0, len(protein_sequence) - peptide_length)

            valid_windows = []
            for start0 in range(start_min, start_max + 1):
                wt_peptide = protein_sequence[start0 : start0 + peptide_length]
                mutation_index0 = position0 - start0
                mt_peptide = (
                    wt_peptide[:mutation_index0]
                    + mutant_aa
                    + wt_peptide[mutation_index0 + 1 :]
                )

                if not set(wt_peptide).issubset(STANDARD_AMINO_ACIDS):
                    skipped_records.append(
                        make_skip_record(
                            row,
                            row_number,
                            "nonstandard_amino_acid_in_peptide",
                            f"start={start0 + 1}; peptide={wt_peptide}",
                            allele,
                        )
                    )
                    continue

                valid_windows.append(
                    (start0, mutation_index0, wt_peptide, mt_peptide)
                )

            k_total = len(valid_windows)
            if k_total == 0:
                skipped_records.append(
                    make_skip_record(
                        row,
                        row_number,
                        "no_valid_peptide_windows",
                        f"target_length={peptide_length}",
                        allele,
                    )
                )
                continue

            for k_index, (start0, mutation_index0, wt_peptide, mt_peptide) in enumerate(
                valid_windows, start=1
            ):
                record = dict(row)
                record.update(
                    {
                        "Source_row": row_number,
                        "Allele": allele,
                        "MHC_class": mhc_class,
                        "Target_peptide_length": peptide_length,
                        "Resolved_RefSeq_protein": resolved_id,
                        "Protein_length": len(protein_sequence),
                        "Peptide_start": start0 + 1,
                        "Peptide_end": start0 + peptide_length,
                        "Mutation_position_in_peptide": mutation_index0 + 1,
                        "K_index": k_index,
                        "K_total": k_total,
                        "WT_peptide": wt_peptide,
                        "MUT_peptide": mt_peptide,
                        "Peptide_pair_id": (
                            f"row{row_number}|{allele}|"
                            f"{resolved_id}:{start0 + 1}-{start0 + peptide_length}"
                        ),
                    }
                )
                output_records.append(record)

    peptides_df = pd.DataFrame(output_records)
    skipped_df = pd.DataFrame(skipped_records)
    return peptides_df, skipped_df


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate every allele-length WT/MT peptide pair that contains "
            "a missense mutation."
        )
    )
    parser.add_argument("--input-table", required=True, help="Input TSV or CSV")
    parser.add_argument("--protein-fasta", required=True, help="RefSeq protein FASTA[.gz]")
    parser.add_argument("--output-peptides", required=True, help="Output peptide-pair TSV")
    parser.add_argument("--output-skipped", required=True, help="Output skipped-row TSV")
    parser.add_argument(
        "--alleles",
        nargs="+",
        default=list(DEFAULT_ALLELES),
        metavar="ALLELE:LENGTH[:CLASS]",
        help=(
            "One or more alleles, each as ALLELE:LENGTH or "
            "ALLELE:LENGTH:CLASS (CLASS is 'I' or 'II'). Not limited to "
            "H2-Kb/H2-Db/I-Ab -- any allele name and any length works, "
            "e.g. --alleles H2-Kb:8 HLA-A*02:01:9:I HLA-DRB1*01:01:15:II. "
            "If CLASS is omitted it is guessed from the allele name "
            "(default: H2-Kb:8:I H2-Db:9:I I-Ab:15:II)"
        ),
    )
    parser.add_argument("--protein-column", default="RefSeq_protein")
    parser.add_argument("--position-column", default="Protein_position")
    parser.add_argument("--reference-aa-column", default="Reference_AA")
    parser.add_argument("--mutant-aa-column", default="Mutant_AA")
    parser.add_argument(
        "--separator",
        default=None,
        help=r"Input separator; inferred as comma for .csv and tab otherwise",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        peptides_df, skipped_df = generate_peptides(args)
    except Exception as exc:
        parser.exit(1, f"ERROR: {exc}\n")

    output_peptides = Path(args.output_peptides)
    output_skipped = Path(args.output_skipped)
    output_peptides.parent.mkdir(parents=True, exist_ok=True)
    output_skipped.parent.mkdir(parents=True, exist_ok=True)

    peptides_df.to_csv(output_peptides, sep="\t", index=False)
    skipped_df.to_csv(output_skipped, sep="\t", index=False)

    print(f"Created peptide pairs: {len(peptides_df):,}")
    if not peptides_df.empty:
        summary = (
            peptides_df.groupby(["Allele", "Target_peptide_length"])
            .size()
            .rename("rows")
        )
        print(summary.to_string())
    print(f"Skipped records: {len(skipped_df):,}")
    if not skipped_df.empty:
        print(skipped_df["Skip_reason"].value_counts().to_string())
    print(f"Peptides saved to: {output_peptides}")
    print(f"Skipped saved to: {output_skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
