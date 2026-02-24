import pandas as pd
from chembl_webresource_client.new_client import new_client
import os
from tqdm import tqdm

# 1. Klasör yapısını kontrol et
if not os.path.exists('data'):
    os.makedirs('data')

def fetch_pi3ka_data():
    print("ChEMBL'den PIK3CA (P42336) verileri çekiliyor... Lütfen bekleyin.")
    
    # ChEMBL API bağlantısı
    target = new_client.target
    # PIK3CA hedefi (UniProt: P42336) sorgusu
    res = target.filter(target_components__accession='P42336')
    if not res:
        print("Hedef bulunamadı!")
        return
    target_id = res[0]['target_chembl_id']
    print(f"Hedef ID bulundu: {target_id}")
    
    # Aktivite verilerini çekme (Sadece IC50 olanlar)
    activities = new_client.activity.filter(target_chembl_id=target_id, standard_type='IC50')
    
    # Veriyi DataFrame'e aktar
    df = pd.DataFrame.from_dict(activities)
    
    # Gerekli sütunları seç ve boş olanları temizle
    df = df[['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'standard_units', 'relation']]
    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
    df = df.dropna(subset=['standard_value', 'canonical_smiles'])
    
    # nM'den uM'ye çevrim (ChEMBL genelde nM verir, 1 uM = 1000 nM)
    df['Activity_uM'] = df['standard_value'] / 1000
    
    # --- MANUEL KRİTERLERİ ---
    # Aktifler: IC50 <= 10 uM
    actives = df[df['Activity_uM'] <= 10].copy()
    actives['Filter_Status'] = 'PASS'
    
    # İnaktifler (Negatif Kontrol): IC50 > 100 uM
    inactives = df[df['Activity_uM'] > 100].copy()
    inactives['Filter_Status'] = 'INACTIVE_CONTROL'
    
    # Birleştirme
    final_df = pd.concat([actives, inactives])
    
    # Manueldeki sütun isimlerine hazırlık
    final_df = final_df.rename(columns={
        'molecule_chembl_id': 'Original_ID',
        'canonical_smiles': 'SMILES'
    })
    final_df['Source'] = 'ChEMBL'
    final_df['Activity_Type'] = 'IC50'
    
    # CSV olarak kaydet
    output_path = 'data/01_chembl_raw.csv'
    final_df.to_csv(output_path, index=False)
    
    print("-" * 30)
    print(f"İşlem Tamamlandı!")
    print(f"Aktif Bileşik Sayısı (<= 10 uM): {len(actives)}")
    print(f"Negatif Kontrol Sayısı (> 100 uM): {len(inactives)}")
    print(f"Dosya kaydedildi: {output_path}")

if __name__ == "__main__":
    fetch_pi3ka_data()