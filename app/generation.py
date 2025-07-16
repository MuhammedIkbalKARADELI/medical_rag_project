def generate_answer(query, documents):
    # 1. prompt hazırla: soru + dokümanlar
    # 2. LLM'den yanıt al (örneğin BioGPT)
    # 3. kaynakları listele
    return {
        "text": "İlgili medikal cevap...",
        "sources": ["PubMed:1234", "NIH:5678"],
        "latency": 345
    }
