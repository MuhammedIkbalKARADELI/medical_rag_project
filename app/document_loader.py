# BU DOSYADA: PubMed'den veri çekme, cümle çıkarma, anahtar kelime filtreleme
from Bio import Entrez
import nltk


Entrez.email = "karadeli2001@hotmail.com"

def fetch_pubmed_abstracts(query, max_count=50):
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


def relevance_filter(abstracts, required_keywords):
    filtered = []
    for entry in abstracts:
        abstract_text = entry["abstract"].lower()
        if any(keyword.lower() in abstract_text for keyword in required_keywords):
            filtered.append(entry)
    return filtered


