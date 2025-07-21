import psycopg2
from config import POSTGRES_CONN_INFO
import time
import psycopg2

def wait_for_postgres(dsn, timeout=60):
    start = time.time()
    while True:
        try:
            conn = psycopg2.connect(dsn)
            conn.close()
            print("PostgreSQL is up!")
            return
        except psycopg2.OperationalError:
            if time.time() - start > timeout:
                raise
            print("Waiting for PostgreSQL...")
            time.sleep(2)


def create_rag_queries_table():
    conn = psycopg2.connect(POSTGRES_CONN_INFO)
    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE IF NOT EXISTS rag_log_queries(
        id SERIAL PRIMARY KEY,
        question TEXT,
        query_length TEXT,
        answer TEXT,
        retrieval_time_ms INTEGER,
        generation_time_ms INTEGER,
        total_time_ms INTEGER,
        bleu_score FLOAT,
        rouge_l_score FLOAT,
        bertscore_f1 FLOAT
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    );
    """)
    conn.commit()
    cur.close()
    conn.close()
    print("rag_log_queries database was checked")


def create_query_documents_table():
    conn = psycopg2.connect(POSTGRES_CONN_INFO)
    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE IF NOT EXISTS query_documents(
        id SERIAL PRIMARY KEY,
        query_id INTEGER REFERENCES rag_log_queries(id) ON DELETE CASCADE,
        pmid VARCHAR(32),
        sentence_index INTEGER,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    );
    """)
    conn.commit()
    cur.close()
    conn.close()
    print("query_documents table was checked")    


if __name__ == "__main__":
    # Bağlantı bilgini gir:
    wait_for_postgres(POSTGRES_CONN_INFO)
    # Bu kod bloğu, sadece ilgili Python dosyası DOĞRUDAN çalıştırıldığında çalışır. Eğer dosya bir modül olarak başka bir dosyada import edilirse, bu blok çalışmaz.
    create_rag_queries_table()
    create_query_documents_table()