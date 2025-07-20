# BU DOSYADA: Arama, embedding, filtreleme, cümle çıkarımı
from sentence_transformers import SentenceTransformer
from app.document_loader import fetch_pubmed_abstracts, extract_sentences, relevance_filter
import faiss
import numpy as np
from app.config import EN_YAKIN_DOKUMAN_SAYISI, GETIRILECEK_DOKUMAN_SAYISI

embed_model = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")

# def retrieve_documents(question):
#     # Anahtar kelimeler (daha iyi bir tokenizer eklersin)
#     keywords = question.split()
#     abstracts = fetch_pubmed_abstracts(question, max_count=GETIRILECEK_DOKUMAN_SAYISI)
#     filtered_abstracts = relevance_filter(abstracts, keywords)
#     detailed_data, sentence_list = extract_sentences(filtered_abstracts)

#     if len(sentence_list) == 0:
#         return [], []

#     # Embedding ve FAISS ile dense arama
#     doc_embeddings = embed_model.encode(sentence_list, convert_to_numpy=True)
#     if len(doc_embeddings.shape) == 1:
#         doc_embeddings = doc_embeddings.reshape(1, -1)

#     dimension = doc_embeddings.shape[1]
#     index = faiss.IndexFlatL2(dimension)
#     index.add(doc_embeddings)
#     id2doc = {i: doc for i, doc in enumerate(sentence_list)}

#     query_vec = embed_model.encode([question], convert_to_numpy=True)
#     distances, indices = index.search(query_vec, EN_YAKIN_DOKUMAN_SAYISI)
#     retrieved_docs = [id2doc[i] for i in indices[0]]

#     return retrieved_docs



def retrieve_documents(question):
    keywords = question.split()
    abstracts = fetch_pubmed_abstracts(question, max_count=GETIRILECEK_DOKUMAN_SAYISI)
    filtered_abstracts = relevance_filter(abstracts, keywords)
    detailed_data, sentence_list = extract_sentences(filtered_abstracts)

    if len(sentence_list) == 0:
        return [], []

    # Embedding ve FAISS ile dense arama
    doc_embeddings = embed_model.encode(sentence_list, convert_to_numpy=True)
    if len(doc_embeddings.shape) == 1:
        doc_embeddings = doc_embeddings.reshape(1, -1)

    dimension = doc_embeddings.shape[1]
    index = faiss.IndexFlatL2(dimension)
    index.add(doc_embeddings)

    query_vec = embed_model.encode([question], convert_to_numpy=True)
    distances, indices = index.search(query_vec, EN_YAKIN_DOKUMAN_SAYISI)

    retrieved_details = [detailed_data[i] for i in indices[0]]
    retrieved_docs = [detail["sentence"] for detail in retrieved_details]
    retrieved_pmid = [detail["pmid"] for detail in retrieved_details]

    return retrieved_docs, retrieved_pmid
