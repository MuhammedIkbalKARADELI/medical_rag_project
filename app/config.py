import os

FAISS_INDEX_PATH = "models/faiss_index/"
EMBEDDING_MODEL = "sentence-transformers/all-MiniLM-L6-v2"
LLM_MODEL_NAME = "microsoft/BioGPT"
DATA_PATH = "data/processed_docs.json"

API_URL = "http://localhost:8000/query" 

# BU DOSYADA: Konfigürasyon sabitlerini buraya koyabilirsin.
EN_YAKIN_DOKUMAN_SAYISI = 10
GETIRILECEK_DOKUMAN_SAYISI = 30

# # PostgreSQL bağlantı bilgileri
POSTGRES_USER_LOCAL = "medical_rag_user"
POSTGRES_PASSWORD_LOCAL = "medical_rag_password"
POSTGRES_DB_LOCAL = "medical_rag_db"
POSTGRES_HOST_LOCAL = "localhost"
POSTGRES_PORT_LOCAL = "5434"
POSTGRES_CONN_INFO_LOCAL = (
    f"dbname='{POSTGRES_DB_LOCAL}' "
    f"user='{POSTGRES_USER_LOCAL}' "
    f"password='{POSTGRES_PASSWORD_LOCAL}' "
    f"host='{POSTGRES_HOST_LOCAL}' "
    f"port='{POSTGRES_PORT_LOCAL}'"
)


POSTGRES_USER_DOCKER = os.getenv("POSTGRES_USER", "medical_rag_user")
POSTGRES_PASSWORD_DOCKER = os.getenv("POSTGRES_PASSWORD", "medical_rag_password")
POSTGRES_DB_DOCKER = os.getenv("POSTGRES_DB", "medical_rag_db")
POSTGRES_HOST_DOCKER = os.getenv("POSTGRES_HOST", "db")  # docker-compose'da host ismi db
POSTGRES_PORT_DOCKER = os.getenv("POSTGRES_PORT", "5432")

POSTGRES_CONN_INFO = (
    f"dbname='{POSTGRES_DB_DOCKER}' "
    f"user='{POSTGRES_USER_DOCKER}' "
    f"password='{POSTGRES_PASSWORD_DOCKER}' "
    f"host='{POSTGRES_HOST_DOCKER}' "
    f"port='{POSTGRES_PORT_DOCKER}'"
)