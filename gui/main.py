import tkinter as tk
from tkinter import messagebox, scrolledtext
import requests

API_URL = "http://localhost:8000/query"   # Docker-compose ile başlattığımda çalışan adres

def send_query():
    query = entry.get("1.0", tk.END).strip()
    if not query:
        messagebox.showwarning("Warning", "Please, give query.")
        return

    send_btn.config(state=tk.DISABLED)
    result_text.delete("1.0", tk.END)
    result_text.insert(tk.END, "Response is waiting...\n")

    try:
        response = requests.post(API_URL, json={"question": query}, timeout=500)
        response.raise_for_status()
        data = response.json()
        # Dilersen burada PubMed_ID'leri ve daha fazlasını da gösterebilirsin.
        yanit = data.get("answer", "There are not any answer.")
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, yanit)
    except Exception as e:
        result_text.delete("1.0", tk.END)
        result_text.insert(tk.END, f"There is a error:\n{e}")
    finally:
        send_btn.config(state=tk.NORMAL)

root = tk.Tk()
root.title("Medical RAG Demo UI")
root.geometry("600x400")

tk.Label(root, text="Give the Query:", font=("Arial", 12, "bold")).pack(pady=8)

entry = tk.Text(root, height=4, width=70, font=("Arial", 11))
entry.pack(padx=10)

send_btn = tk.Button(root, text="Send Query", command=send_query, font=("Arial", 11, "bold"), bg="#38b6ff")
send_btn.pack(pady=8)

tk.Label(root, text="Response of the model:", font=("Arial", 12, "bold")).pack(pady=4)
result_text = scrolledtext.ScrolledText(root, height=10, width=70, font=("Arial", 11))
result_text.pack(padx=10, pady=4)

root.mainloop()
