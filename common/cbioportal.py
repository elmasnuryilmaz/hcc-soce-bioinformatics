from __future__ import annotations

import requests
import pandas as pd

CBIO_BASE = "https://www.cbioportal.org/api"


def fetch_entrez_ids(genes: list[str], base: str = CBIO_BASE) -> dict[str, int]:
    """Resolve HUGO symbols to Entrez IDs using the cBioPortal API."""
    response = requests.post(
        f"{base}/genes/fetch",
        params={"geneIdType": "HUGO_GENE_SYMBOL"},
        json=list(dict.fromkeys(genes)),
        timeout=60,
    )
    response.raise_for_status()
    return {row["hugoGeneSymbol"]: int(row["entrezGeneId"]) for row in response.json()}


def fetch_mrna_expression(
    genes: list[str],
    profile: str = "lihc_tcga_rna_seq_v2_mrna",
    sample_list_id: str = "lihc_tcga_all",
    base: str = CBIO_BASE,
) -> pd.DataFrame:
    """Fetch a samples x genes expression matrix from cBioPortal."""
    gene_map = fetch_entrez_ids(genes, base=base)
    if not gene_map:
        raise RuntimeError("No genes were resolved by cBioPortal.")

    response = requests.post(
        f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
        json={
            "entrezGeneIds": list(gene_map.values()),
            "sampleListId": sample_list_id,
        },
        timeout=180,
    )
    response.raise_for_status()

    entrez_to_symbol = {entrez: symbol for symbol, entrez in gene_map.items()}
    data: dict[str, dict[str, float]] = {}
    for row in response.json():
        symbol = entrez_to_symbol.get(int(row["entrezGeneId"]))
        if symbol is None:
            continue
        data.setdefault(row["sampleId"], {})[symbol] = row["value"]

    df = pd.DataFrame.from_dict(data, orient="index")
    df = df.apply(pd.to_numeric, errors="coerce")
    return df.reindex(columns=[g for g in genes if g in df.columns])


def fetch_molecular_profile_data(
    genes: list[str],
    profile: str,
    sample_list_id: str,
    base: str = CBIO_BASE,
) -> pd.DataFrame:
    """Fetch long-form molecular profile rows and add HUGO gene symbols."""
    gene_map = fetch_entrez_ids(genes, base=base)
    if not gene_map:
        raise RuntimeError("No genes were resolved by cBioPortal.")

    response = requests.post(
        f"{base}/molecular-profiles/{profile}/molecular-data/fetch",
        json={
            "entrezGeneIds": list(gene_map.values()),
            "sampleListId": sample_list_id,
        },
        timeout=180,
    )
    response.raise_for_status()

    df = pd.DataFrame(response.json())
    if df.empty:
        return df
    entrez_to_symbol = {entrez: symbol for symbol, entrez in gene_map.items()}
    df["hugoGeneSymbol"] = df["entrezGeneId"].astype(int).map(entrez_to_symbol)
    if "value" in df.columns:
        df["value"] = pd.to_numeric(df["value"], errors="coerce")
    return df


def fetch_mutations(
    genes: list[str],
    profile: str,
    sample_list_id: str,
    base: str = CBIO_BASE,
) -> pd.DataFrame:
    """Fetch mutation rows from a cBioPortal mutation profile."""
    gene_map = fetch_entrez_ids(genes, base=base)
    if not gene_map:
        raise RuntimeError("No genes were resolved by cBioPortal.")

    response = requests.post(
        f"{base}/molecular-profiles/{profile}/mutations/fetch",
        json={
            "entrezGeneIds": list(gene_map.values()),
            "sampleListId": sample_list_id,
        },
        timeout=180,
    )
    response.raise_for_status()

    df = pd.DataFrame(response.json())
    if df.empty:
        return df
    entrez_to_symbol = {entrez: symbol for symbol, entrez in gene_map.items()}
    df["hugoGeneSymbol"] = df["entrezGeneId"].astype(int).map(entrez_to_symbol)
    return df


def fetch_clinical_data(
    study_id: str = "lihc_tcga",
    base: str = CBIO_BASE,
) -> pd.DataFrame:
    """Fetch sample-level clinical rows plus patient-level survival covariates."""
    sample_response = requests.get(
        f"{base}/studies/{study_id}/clinical-data",
        params={"clinicalDataType": "SAMPLE"},
        timeout=120,
    )
    sample_response.raise_for_status()
    sample_rows = sample_response.json()
    sample_df = pd.DataFrame(sample_rows)
    if sample_df.empty:
        raise RuntimeError(f"No sample clinical data returned for {study_id}.")

    sample_attrs = sample_df.pivot(
        index="sampleId", columns="clinicalAttributeId", values="value"
    )
    sample_patients = (
        sample_df[["sampleId", "patientId"]]
        .drop_duplicates()
        .set_index("sampleId")
    )
    clinical = sample_attrs.join(sample_patients, how="left")

    patient_response = requests.get(
        f"{base}/studies/{study_id}/clinical-data",
        params={"clinicalDataType": "PATIENT"},
        timeout=120,
    )
    patient_response.raise_for_status()
    patient_rows = patient_response.json()
    patient_df = pd.DataFrame(patient_rows)
    if not patient_df.empty:
        patient_attrs = patient_df.pivot(
            index="patientId", columns="clinicalAttributeId", values="value"
        )
        clinical = clinical.join(patient_attrs, on="patientId", rsuffix="_PATIENT")

    return clinical
