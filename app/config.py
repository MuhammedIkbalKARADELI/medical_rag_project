FAISS_INDEX_PATH = "models/faiss_index/"
EMBEDDING_MODEL = "sentence-transformers/all-MiniLM-L6-v2"
LLM_MODEL_NAME = "microsoft/BioGPT"
DATA_PATH = "data/processed_docs.json"


# BU DOSYADA: Konfigürasyon sabitlerini buraya koyabilirsin.
EN_YAKIN_DOKUMAN_SAYISI = 10
GETIRILECEK_DOKUMAN_SAYISI = 30

# PostgreSQL bağlantı bilgileri
POSTGRES_USER = "medical_rag_user"
POSTGRES_PASSWORD = "medical_rag_password"
POSTGRES_DB = "medical_rag_db"
POSTGRES_HOST = "localhost"
POSTGRES_PORT = "5434"

API_URL = "http://localhost:8000/query" 

POSTGRES_CONN_INFO = (
    f"dbname='{POSTGRES_DB}' "
    f"user='{POSTGRES_USER}' "
    f"password='{POSTGRES_PASSWORD}' "
    f"host='{POSTGRES_HOST}' "
    f"port='{POSTGRES_PORT}'"
)