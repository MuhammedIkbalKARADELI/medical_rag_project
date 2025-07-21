import requests
import time
import json

API_URL = "http://localhost:8000/query"   # Kendi API adresini buraya yaz


with open("/Users/ikbalkrdl/Desktop/Baykar Sağlıkta Yapay Zeka/medical_rag_project/scripts/150_question.json", "r", encoding="utf-8") as f:
    questions = json.load(f)

results = []

for i, question in enumerate(questions, 1):
    payload = {"question": question}
    try:
        response = requests.post(API_URL, json=payload, timeout=5000)
        response.raise_for_status()
        result = response.json()
        print(f"{i}. Soru gönderildi: {question[:40]}... [Cevap geldi]")
        results.append(result)
    except Exception as e:
        print(f"{i}. Soru gönderilirken hata oluştu: {e}")
    time.sleep(0.3)  # API’yı boğmamak için biraz bekleyebilirsin

# Sonuçları JSON dosyasına kaydetmek için (isteğe bağlı):
import json
with open("api_results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print("Tüm sorgular tamamlandı ve sonuçlar kaydedildi.")