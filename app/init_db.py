import psycopg2
from config import POSTGRES_CONN_INFO

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
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    );
    """)
    conn.commit()
    cur.close()
    conn.close()
    print("rag_queries database was checked")


def create_query_documents_table():
    conn = psycopg2.connect(POSTGRES_CONN_INFO)
    cur = conn.cursor()
    cur.execute("""
    CREATE TABLE IF NOT EXISTS query_documents(
        id SERIAL PRIMARY KEY,
        query_id INTEGER REFERENCES rag_queries(id) ON DELETE CASCADE,
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
    # Bu kod bloğu, sadece ilgili Python dosyası DOĞRUDAN çalıştırıldığında çalışır. Eğer dosya bir modül olarak başka bir dosyada import edilirse, bu blok çalışmaz.
    create_rag_queries_table()
    create_query_documents_table()