from fastapi import FastAPI, Request
from pydantic import BaseModel
from sentence_transformers import SentenceTransformer
from transformers import pipeline
import faiss
import numpy as np
import json
from transformers import AutoTokenizer, AutoModelForCausalLM
from Bio import Entrez
import nltk


# NLTK tokenizer için indirme (ilk çalıştırmada gerekli)
# nltk.download('punkt')
# nltk.download("punkt_tab")

# PubMed için mail adresin
Entrez.email = "karadeli2001@hotmail.com"


def relevance_filter(abstracts, required_keywords):
    filtered = []
    for entry in abstracts:
        abstract_text = entry["abstract"].lower()
        if any(keyword.lower() in abstract_text for keyword in required_keywords):
            filtered.append(entry)
    return filtered


def fetch_pubmed_abstracts(query, max_count=50):
    """PubMed'den özetleri alır"""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_count)
    record = Entrez.read(handle)
    ids = record["IdList"]
    
    abstracts = []
    for pmid in ids:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read().strip()
        if abstract:
            abstracts.append({
                "pmid": pmid,
                "abstract": abstract
            })
    return abstracts

def extract_sentences(abstracts):
    """Hem detaylı hem sadece cümle listesi oluşturur"""
    detailed = []
    only_sentences = []
    for entry in abstracts:
        sentences = nltk.sent_tokenize(entry["abstract"])
        for i, sentence in enumerate(sentences):
            detailed.append({
                "pmid": entry["pmid"],
                "sentence_id": i,
                "sentence": sentence
            })
            only_sentences.append(sentence)
    return detailed, only_sentences

def save_json(data, filename):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)



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
        # top_k=15,
        top_p=0.95,
        temperature=0.8) 

    return llm_tokenizer.decode(outputs[0], skip_special_tokens=True)



# ---- Request Modeli ----
class QueryRequest(BaseModel):
    question: str


@app.post("/query")
def query_endpoint(request: QueryRequest):


    keywords = request.question.split()
    # print(keywords)
    # 1. Makaleleri çek
    abstracts = fetch_pubmed_abstracts(request.question, max_count=50)
    filtered_abstracts = relevance_filter(abstracts, keywords)

    # 2. Cümleleri ayır
    detailed_data, sentence_list = extract_sentences(filtered_abstracts)
    save_json(sentence_list, "pubmed_sentences_only.json")
    save_json(detailed_data, "pubmed_abstract.json")

    # with open("pubmed_sentences_only.json", "r", encoding="utf-8") as f:
    #     docs = json.load(f)

    if len(sentence_list) == 0:
        return {"There is no any related article give me a different question..."}


    # ---- Vektörleri Hazırla ve FAISS Index Kur ----
    doc_embeddings = embed_model.encode(sentence_list, convert_to_numpy=True)

    # Şekil kontrolü
    if len(doc_embeddings.shape) == 1:
        doc_embeddings = doc_embeddings.reshape(1, -1)  # Tek belge ise 2D'ye çevir

    dimension = doc_embeddings.shape[1]  # Artık güvenli
    index = faiss.IndexFlatL2(dimension)
    index.add(doc_embeddings)
    id2doc = {i: doc for i, doc in enumerate(sentence_list)}

    # 1. Sorguyu embed et
    query_vec = embed_model.encode([request.question], convert_to_numpy=True)

    en_yakin_döküman_sayisi = 20
    # 2. FAISS ile en benzer dökümanları bul
    distances, indices = index.search(query_vec, en_yakin_döküman_sayisi)
    retrieved_docs = [id2doc[i] for i in indices[0]]
    # print(retrieved_docs)

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
        "Query_Length": "Kisa-orta-uzun",
        "Retrieval_Time_MS" : "12_sn",
        "Generation_Time_MS" : "13_sn",
        "Total_Time_MS" : "25_sn"
        # "retrieval_count": top_k
    }
