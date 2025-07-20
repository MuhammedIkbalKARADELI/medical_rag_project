import tkinter as tk
from tkinter import messagebox, scrolledtext
import requests

API_URL = "http://localhost:8000/query"   # Docker-compose ile başlattığımda çalışan adres

def send_query():
    query = entry.get("1.0", tk.END).strip()
    if not query:
        messagebox.showwarning("Uyarı", "Lütfen bir sorgu girin.")
        return

    send_btn.config(state=tk.DISABLED)
    result_text.delete("1.0", tk.END)
    result_text.insert(tk.END, "Yanıt bekleniyor...\n")

    try:
        response = requests.post(API_URL, json={"question": query}, timeout=60)
        response.raise_for_status()
        data = response.json()
        # Dilersen burada PubMed_ID'leri ve daha fazlasını da gösterebilirsin.
        yanit = data.get("answer", "Cevap alınamadı.")
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, yanit)
    except Exception as e:
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, f"Bir hata oluştu:\n{e}")
    finally:
        send_btn.config(state=tk.NORMAL)

root = tk.Tk()
root.title("Medical RAG Demo UI")
root.geometry("600x400")

tk.Label(root, text="Sorgunuzu Girin:", font=("Arial", 12, "bold")).pack(pady=8)

entry = tk.Text(root, height=4, width=70, font=("Arial", 11))
entry.pack(padx=10)

send_btn = tk.Button(root, text="Sorguyu Gönder", command=send_query, font=("Arial", 11, "bold"), bg="#38b6ff")
send_btn.pack(pady=8)

tk.Label(root, text="Model Yanıtı:", font=("Arial", 12, "bold")).pack(pady=4)
result_text = scrolledtext.ScrolledText(root, height=10, width=70, font=("Arial", 11))
result_text.pack(padx=10, pady=4)

root.mainloop()
