from fastapi import FastAPI, Request
from pydantic import BaseModel
from sentence_transformers import SentenceTransformer
from transformers import pipeline
import faiss
import numpy as np
import json

app = FastAPI()

# ---- Embedding Model ----
embed_model = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")

# ---- LLM Model ----
llm = pipeline("text-generation", model="gpt2", max_new_tokens=100)

# ---- Veri Yükle ----
with open("data.json", "r", encoding="utf-8") as f:
    docs = json.load(f)

# ---- Vektörleri Hazırla ve FAISS Index Kur ----
doc_embeddings = embed_model.encode(docs, convert_to_numpy=True)
dimension = doc_embeddings.shape[1]
index = faiss.IndexFlatL2(dimension)
index.add(doc_embeddings)
id2doc = {i: doc for i, doc in enumerate(docs)}

# ---- Request Modeli ----
class QueryRequest(BaseModel):
    question: str
    top_k: int = 3


@app.post("/query")
def query_endpoint(request: QueryRequest):
    # 1. Sorguyu embed et
    query_vec = embed_model.encode([request.question], convert_to_numpy=True)

    # 2. FAISS ile en benzer dökümanları bul
    distances, indices = index.search(query_vec, request.top_k)
    retrieved_docs = [id2doc[i] for i in indices[0]]

    # 3. LLM'e prompt hazırla ve cevap oluştur
    context = "\n".join(retrieved_docs)
    prompt = f"Context:\n{context}\n\nQuestion: {request.question}\nAnswer:"
    response = llm(prompt)[0]["generated_text"].split("Answer:")[-1].strip()

    # 4. JSON yanıtla
    return {
        "question": request.question,
        "top_documents": retrieved_docs,
        "answer": response,
        "retrieval_count": request.top_k
    }
