import psycopg2
from psycopg2.extras import execute_values
from config import POSTGRES_CONN_INFO

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'app')))


def rag_log_query(question, answer, query_length, retrieval_time, generation_time, total_time, bleu_score, rouge_l_score, bertscore_f1):
    conn = None
    query_id = None
    try:
        conn = psycopg2.connect(POSTGRES_CONN_INFO)
        cur = conn.cursor()
        cur.execute("""
            INSERT INTO rag_log_queries
            (question, answer, query_length, retrieval_time, generation_time, total_time, bleu_score, rouge_l_score, bertscore_f1)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
            RETURNING id
        """, (question, answer, query_length, retrieval_time, generation_time, total_time, bleu_score, rouge_l_score, bertscore_f1))
        query_id = cur.fetchone()[0]  # Yeni eklenen kaydın ID'sini alır
        conn.commit()
        cur.close()
    except Exception as e:
        print(f"Log DB error: {e}")
    finally:
        if conn:
            conn.close()
    return query_id



def insert_query_documents(query_id, pubmed_ids):
    conn = psycopg2.connect(POSTGRES_CONN_INFO)
    cur = conn.cursor()
    for idx, pmid in enumerate(pubmed_ids):
        cur.execute("""
            INSERT INTO query_documents (query_id, pmid, sentence_index)
            VALUES (%s, %s, %s)
        """, (query_id, pmid, idx))
    conn.commit()
    cur.close()
    conn.close()