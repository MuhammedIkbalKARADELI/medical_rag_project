from Bio import Entrez
import nltk
import json

# NLTK tokenizer için indirme (ilk çalıştırmada gerekli)
# nltk.download('punkt')
# nltk.download("punkt_tab")

# PubMed için mail adresin
Entrez.email = "karadeli2001@hotmail.com"

def fetch_pubmed_abstracts(query, max_count=5):
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

# 🔍 Sorgu gir
query = "diabet"
max_articles = 10

# 1. Makaleleri çek
abstracts = fetch_pubmed_abstracts(query, max_count=max_articles)

# 2. Cümleleri ayır
detailed_data, sentence_list = extract_sentences(abstracts)

# 3. JSON dosyalarını oluştur
save_json(detailed_data, "pubmed_sentences_detailed.json")
save_json(sentence_list, "pubmed_sentences_only.json")

# print("✅ JSON dosyaları oluşturuldu:")
# print("- pubmed_sentences_detailed.json")
# print("- pubmed_sentences_only.json")
