import pandas as pd


NAMA_FILE = "Book1.xlsx"   # ganti dengan nama file 
NAMA_SHEET = 0                # untuk Excel: 0 = sheet pertama,
                              # atau tulis nama sheet, misalnya "Sheet1"
# =====================================================

path = "data/" + NAMA_FILE

# Baca file sesuai jenisnya
if NAMA_FILE.lower().endswith((".xlsx", ".xls")):
    df = pd.read_excel(path, sheet_name=NAMA_SHEET)
elif NAMA_FILE.lower().endswith(".csv"):
    # sep=None membuat Python menebak pemisah (koma atau titik koma)
    df = pd.read_csv(path, sep=None, engine="python")
else:
    raise ValueError("Format file harus .xlsx, .xls, atau .csv")

# Tampilkan semua kolom saat dicetak
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)

print("=" * 60)
print("JUMLAH BARIS DAN KOLOM:", df.shape)

print("\n" + "=" * 60)
print("NAMA KOLOM DAN JENIS DATA:")
print(df.dtypes)

print("\n" + "=" * 60)
print("5 BARIS PERTAMA:")
print(df.head())

print("\n" + "=" * 60)
print("JUMLAH DATA KOSONG PER KOLOM:")
print(df.isnull().sum())

print("\n" + "=" * 60)
print("RINGKASAN ANGKA:")
print(df.describe().T)