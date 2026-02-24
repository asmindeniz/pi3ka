#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import io
import time
import requests
import pandas as pd

# =========================
# AYARLAR
# =========================
GENE_ID = 5290  # PIK3CA (p110α)
OUTDIR = "data"
OUTCSV = os.path.join(OUTDIR, "raw_pubchem_pi3ka.csv")

# İlk denemede hız için limit koymak istersen örn: 20
# Tam tarama için None bırak.
LIMIT_AIDS = None # None yaparsan hepsini dener

# PubChem'e nazik davranalım
SLEEP_AID = 0.2
SLEEP_CID_BATCH = 0.15

# CID->SMILES batch boyutu (PubChem endpoint limitlerine takılmamak için)
CID_BATCH_SIZE = 200

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "pi3ka-pubchem-fetch/1.0"})


# =========================
# YARDIMCI FONKSİYONLAR
# =========================
def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def get_aids_for_gene(gene_id: int) -> list[int]:
    """
    PubChem PUG-REST: assay/target/geneid/<geneid>/aids/TXT
    TXT dönmesi JSON kararsızlıklarını azaltıyor.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/geneid/{gene_id}/aids/TXT"
    r = SESSION.get(url, timeout=60)
    if r.status_code != 200:
        raise RuntimeError(f"AID listesi alınamadı. HTTP {r.status_code}. Yanıtın başı: {r.text[:200]}")

    # TXT genelde satır satır AID listeler
    lines = [ln.strip() for ln in r.text.splitlines() if ln.strip().isdigit()]
    aids = [int(x) for x in lines]
    if not aids:
        # Bazen TXT formatı farklı olabilir; debug için başı yazdır
        raise RuntimeError(f"AID listesi boş geldi. Yanıtın başı: {r.text[:300]}")
    return aids


def fetch_assay_csv(aid: int) -> pd.DataFrame | None:
    """
    AID için assay tablosunu CSV olarak indir.
    Bazı AID'ler çok büyük/garip olabilir, hata olursa None döndür.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/CSV"
    try:
        r = SESSION.get(url, timeout=120)
        if r.status_code != 200:
            return None
        # CSV parse
        return pd.read_csv(io.StringIO(r.text))
    except Exception:
        return None


def normalize_outcome(x) -> str | None:
    """
    PubChem activity outcome bazen string ("Active"), bazen numeric code olabilir.
    Kodlar: 1=Inactive, 2=Active, 3=Inconclusive, 4=Unspecified
    """
    if pd.isna(x):
        return None

    if isinstance(x, str):
        s = x.strip().lower()
        if s == "active":
            return "Active"
        if s == "inactive":
            return "Inactive"
        return None

    # numeric
    try:
        v = int(x)
        if v == 2:
            return "Active"
        if v == 1:
            return "Inactive"
    except Exception:
        pass
    return None


def pick_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """
    Aday kolon isimlerinden df'de olanı bul.
    """
    cols = set(df.columns)
    for c in candidates:
        if c in cols:
            return c
    return None


def cid_to_smiles_batch(cids: list[int]) -> dict[int, str]:
    """
    CID listesini batch halinde CanonicalSMILES'a çevirir.
    """
    out = {}
    if not cids:
        return out

    # PubChem endpoint: /compound/cid/<comma-separated>/property/CanonicalSMILES/JSON
    # Çok uzun URL olmasın diye batchliyoruz.
    for i in range(0, len(cids), CID_BATCH_SIZE):
        batch = cids[i:i + CID_BATCH_SIZE]
        cid_str = ",".join(str(x) for x in batch)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/CanonicalSMILES/JSON"
        try:
            r = SESSION.get(url, timeout=120)
            if r.status_code != 200:
                time.sleep(SLEEP_CID_BATCH)
                continue
            data = r.json()
            props = data.get("PropertyTable", {}).get("Properties", [])
            for p in props:
                cid = p.get("CID")
                smi = p.get("CanonicalSMILES")
                if cid is not None and smi:
                    out[int(cid)] = smi
        except Exception:
            pass

        time.sleep(SLEEP_CID_BATCH)

    return out


def aggregate_outcome_per_cid(rows: list[dict]) -> pd.DataFrame:
    """
    Aynı CID farklı assay'lerde hem Active hem Inactive görünebilir.
    Basit ve pratik kural:
      - Herhangi bir assay'de Active ise: Active
      - değilse Inactive
    Ayrıca kaç farklı AID'den geldiğini sayar.
    """
    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # CID başına AID sayısı
    aid_counts = df.groupby("CID")["AID"].nunique().rename("Assay_Count")

    # CID başına outcome: Active öncelikli
    def reduce_outcomes(series: pd.Series) -> str:
        s = set(series.dropna().tolist())
        if "Active" in s:
            return "Active"
        if "Inactive" in s:
            return "Inactive"
        return "Unspecified"

    outcome = df.groupby("CID")["Outcome"].apply(reduce_outcomes).rename("Activity")
    out = pd.concat([outcome, aid_counts], axis=1).reset_index()
    return out


# =========================
# ANA AKIŞ
# =========================
def main():
    ensure_outdir(OUTDIR)

    print(f"[1/4] PIK3CA (GeneID={GENE_ID}) için AID listesi çekiliyor...")
    aids = get_aids_for_gene(GENE_ID)
    print(f"  Bulunan AID sayısı: {len(aids)}")

    if LIMIT_AIDS is not None:
        aids = aids[:LIMIT_AIDS]
        print(f"  Test modu: İlk {LIMIT_AIDS} AID işlenecek.")

    print(f"[2/4] AID'lerden Active/Inactive CID'ler toplanıyor...")
    records = []
    skipped = 0

    for idx, aid in enumerate(aids, start=1):
        df = fetch_assay_csv(aid)
        if df is None or df.empty:
            skipped += 1
            continue

        # Kolon adları assay'e göre değişebiliyor; olası isimleri deniyoruz
        cid_col = pick_column(df, ["PUBCHEM_CID", "CID", "cid"])
        out_col = pick_column(df, ["Activity Outcome", "PUBCHEM_ACTIVITY_OUTCOME", "OUTCOME", "outcome"])

        if cid_col is None or out_col is None:
            skipped += 1
            continue

        # Sadece Active/Inactive topla
        for _, row in df[[cid_col, out_col]].iterrows():
            cid_val = row[cid_col]
            outcome_raw = row[out_col]

            try:
                cid = int(cid_val)
            except Exception:
                continue

            outcome = normalize_outcome(outcome_raw)
            if outcome in ("Active", "Inactive"):
                records.append({"AID": aid, "CID": cid, "Outcome": outcome})

        if idx % 10 == 0:
            print(f"  İşlenen AID: {idx}/{len(aids)} | kayıt: {len(records)} | atlanan AID: {skipped}")

        time.sleep(SLEEP_AID)

    if not records:
        raise RuntimeError("Hiç Active/Inactive kayıt toplanamadı. LIMIT_AIDS'i artırıp tekrar dene veya endpoint bloklanmış olabilir.")

    print(f"  Toplam kayıt (AID-CID): {len(records)}")

    print(f"[3/4] CID bazında birleştiriliyor (Active öncelikli)...")
    df_cid = aggregate_outcome_per_cid(records)
    print("  CID bazlı dağılım:")
    print(df_cid["Activity"].value_counts(dropna=False))

    # SMILES çek
    print(f"[4/4] CID -> CanonicalSMILES çekiliyor (batch)...")
    unique_cids = df_cid["CID"].astype(int).tolist()
    smiles_map = cid_to_smiles_batch(unique_cids)

    df_cid["SMILES"] = df_cid["CID"].map(smiles_map)

    # 1) Ham tablo: CID + Activity + Assay_Count (+SMILES olabilir/olmayabilir)
    OUTCSV_RAW = os.path.join(OUTDIR, "raw_pubchem_pi3ka_cid_activity.csv")
    df_cid.to_csv(OUTCSV_RAW, index=False)
    print(f"✅ Ham tablo kaydedildi (SMILES eksik olabilir): {OUTCSV_RAW}")

    # 2) Sadece SMILES olanlar
    df_with_smiles = df_cid.dropna(subset=["SMILES"]).reset_index(drop=True)
    OUTCSV_SMILES = os.path.join(OUTDIR, "raw_pubchem_pi3ka_with_smiles.csv")
    df_with_smiles.to_csv(OUTCSV_SMILES, index=False)
    print(f"✅ SMILES olanlar kaydedildi: {OUTCSV_SMILES}")

    print("SMILES doluluk oranı:", len(df_with_smiles), "/", len(df_cid))
    print("✅ Final dağılım (SMILES olanlar):")
    print(df_with_smiles["Activity"].value_counts(dropna=False))


if __name__ == "__main__":
    main()