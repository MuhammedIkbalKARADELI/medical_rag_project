from fastapi import FastAPI, Request
from pydantic import BaseModel
from sentence_transformers import SentenceTransformer
from transformers import pipeline
import faiss
import numpy as np
import json
from transformers import AutoTokenizer, AutoModelForCausalLM
from data_loader import fetch_pubmed_abstracts, extract_sentences, save_json

app = FastAPI() # Bu fast api yi çalıştırmak için = uvicorn app4:app --reload    komutunu kullanman lazım

# ---- Embedding Model ----
embed_model = SentenceTransformer("sentence-transformers/all-MiniLM-L6-v2")


# ---- LLM Model: BioGPT Update ----
llm_tokenizer = AutoTokenizer.from_pretrained("microsoft/BioGPT")
llm_model = AutoModelForCausalLM.from_pretrained("microsoft/BioGPT")

def generate_answer(prompt: str) -> str:
    inputs = llm_tokenizer(prompt, return_tensors="pt")
    # outputs = llm_model.generate(**inputs, max_new_tokens=300)
    outputs = llm_model.generate(
    **inputs,
    max_new_tokens=300,
    do_sample=True,
    top_k=50,
    top_p=0.95,
    temperature=0.8,
)

    return llm_tokenizer.decode(outputs[0], skip_special_tokens=True)



# ---- Request Modeli ----
class QueryRequest(BaseModel):
    question: str
    top_k: int = 3


@app.post("/query")
def query_endpoint(request: QueryRequest):


    # 🔍 Sorgu gir
    # query = "covid-19 vaccine immune response"
    max_articles = 30

    # 1. Makaleleri çek
    abstracts = fetch_pubmed_abstracts(request.question, max_count=max_articles)

    # 2. Cümleleri ayır
    detailed_data, sentence_list = extract_sentences(abstracts)
    save_json(sentence_list, "pubmed_sentences_only.json")

    with open("pubmed_sentences_only.json", "r", encoding="utf-8") as f:
        docs = json.load(f)

    if len(docs) == 0:
        return {"There is no any related article give me a different question..."}


    # ---- Vektörleri Hazırla ve FAISS Index Kur ----
    # doc_embeddings = embed_model.encode(sentence_list, convert_to_numpy=True)
    doc_embeddings = embed_model.encode(docs, convert_to_numpy=True)

    # Şekil kontrolü
    if len(doc_embeddings.shape) == 1:
        doc_embeddings = doc_embeddings.reshape(1, -1)  # Tek belge ise 2D'ye çevir

    dimension = doc_embeddings.shape[1]  # Artık güvenli
    index = faiss.IndexFlatL2(dimension)
    index.add(doc_embeddings)
    id2doc = {i: doc for i, doc in enumerate(docs)}

        
    # 1. Sorguyu embed et
    query_vec = embed_model.encode([request.question], convert_to_numpy=True)

    # 2. FAISS ile en benzer dökümanları bul
    distances, indices = index.search(query_vec, request.top_k)
    retrieved_docs = [id2doc[i] for i in indices[0]]

    # 3. BioGPT için prompt oluştur
    context = "\n".join(retrieved_docs)
    prompt = f"{context}\n\nQuestion: {request.question}\nAnswer:"

    # 4. LLM cevabını oluştur
    response = generate_answer(prompt).split("Answer:")[-1].strip()

    # 5. JSON yanıtla
    return {
        "question": request.question,
        "top_documents": retrieved_docs,
        "answer": response,
        "retrieval_count": request.top_k
    }
