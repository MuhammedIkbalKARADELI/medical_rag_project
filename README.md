# Medical RAG Project: PubMed-Supported Medical Question Answering System

This project is a **research and development (R&D) study** aimed at automatically answering medical questions with real-time support from up-to-date scientific literature.  
It combines user queries with PubMed-sourced documents and generates reliable, evidence-based answers using a modern large language model (LLM, e.g., BioGPT).

---

## 🔬 Project Objective

- To develop a **scientifically grounded, automated medical question-answering system**
- To generate accurate and referenceable answers based on recent literature
- To record every user query, its properties, and the model’s response in a **relational database** for full transparency
- To link each query log with the exact PubMed documents used for generating the answer
- To provide full auditability and reproducibility for academic research and analysis

---

## 🗃️ Project Directory Structure

---

## 💾 **Database & Logging Architecture**

The project features a **PostgreSQL-based, fine-grained logging system** designed for research transparency:
You can find the answers to more than 100 questions we queried in the project and saved in the PostgreSQL database in the output file.


- Each user query is recorded in the **rag_log_queries** table with:
  - ID (Match the Query ID of the query_documents table)
  - The query text
  - Query length (short/medium/long)
  - Model-generated answer
  - Document retrieval time (ms)
  - Answer generation time (ms)
  - Total response time (ms)
  - Query timestamp
  - Bleu score
  - Rouge-L score
  - Bert Score F1 
  - Query timestamp
  
Output: [text](Output/rag-log-queries.csv)

- Every PubMed document (or sentence) used for each query is separately stored in the **query_documents** table, **linked to the log entry**:
  - ID
  - Query ID (foreign key to rag_log_queries table)
  - PubMed ID (PMID)
  - Document/sentence index

Output: [text](Output/query-documents.csv)

**Benefits:**
- All steps are stored in a relational database for **auditing and analysis**
- Model performance and answer quality can be quantitatively evaluated
- Complete reproducibility: past queries and corresponding documents are fully tracked


---

## 🏗️ **Setup & Deployment (with Docker)**

1. **Clone the repository**
    ```bash
    git clone https://github.com/MuhammedIkbalKARADELI/medical_rag_project.git
    cd medical_rag_project
    ```

2. **Build Docker Images**
    ```bash
    docker compose build
    ```

3. **Launch Services**
    ```bash
    docker compose up
    ```
    - API: [http://localhost:8000](http://localhost:8000)
    - PostgreSQL: localhost:5434

Models and data files will be downloaded automatically on the first run.

---

## ⚙️ **Technology Stack**

- **FastAPI:** REST API service
- **sentence-transformers/all-MiniLM-L6-v2:** Embeddings for queries and documents
- **FAISS:** Dense vector-based semantic search
- **BioGPT:** Medical language model for answer generation
- **PostgreSQL:** Advanced logging and document tracking system
- **Tkinter GUI:** Local interactive desktop demo
- **Docker & docker-compose:** Fully isolated, reproducible environment

---


## 🖥️ **Usage & Examples**

### 1. **API Status Check**
```bash
curl http://localhost:8000/status
```


### 2. **Ask a Medical Question (POST /query)**
```bash
import requests
response = requests.post(
    "http://localhost:8000/query",
    json={"question": "What are the symptoms of anemia?"}
)
print(response.json())
```
Sample Response:

```bash
{
    "question": "...",
    "top_documents": [...],
    "PubMed_ID": [...],
    "answer": "...",
    "Query_Length": "...",
    "Retrieval_Time_MS": ...,
    "Generation_Time_MS": ...,
    "Total_Time_MS": ...
}
```
Every query and its associated documents are persistently logged in the postgres database.


## 🗨️ ** GUI Demo**
```bash
cd gui
python main.py
```
While the API is running, use this desktop interface to submit questions and view model responses.


## 🔎 **Bulk and Automated Querying**
```bash
cd scripts
python send_bulk_queries.py
```
Sends 50 different medical questions to the API and only saves all responses as a JSON file and insterted the pstgres-db.


## 🧪 **Testing**

###  **To run FastAPI endpoint tests:**
```bash
cd fast_api_test
pytest
```

### **For manual testing of the Mediacl Rag Project:**
```bash
cd test
python test_query.py
```


## 📝 **Academic Transparency & Reproducibility**
Every query, its relevant PubMed documents, and all model outputs are persistently stored in the database, supporting scientific integrity and repeatability.

Full logs allow for post-hoc statistical analysis, error inspection, and reproducibility of results.



## ⚠️ **Disclaimer**
This system is not intended for diagnostic use; it is strictly for research and educational purposes only.
All model-generated answers must be evaluated by qualified healthcare professionals.


## 👤 **Author / Contact**
Developer: Muhammed İkbal KARADELİ

Contact: karadeli2001@hotmail.com


---

## 📚 References

1. Lewis, Patrick, et al. "Retrieval-augmented generation for knowledge-intensive NLP tasks." *Advances in Neural Information Processing Systems* 33 (2020): 9459-9474.  
   [https://arxiv.org/abs/2005.11401](https://arxiv.org/abs/2005.11401)

2. Luo, Ruoqi, et al. "BioGPT: Generative Pre-trained Transformer for Biomedical Text Generation and Mining." *bioRxiv* (2022): 2022-08.  
   [https://arxiv.org/abs/2210.10341](https://arxiv.org/abs/2210.10341)

3. Reimers, Nils, and Iryna Gurevych. "Sentence-BERT: Sentence Embeddings using Siamese BERT-Networks." *EMNLP 2019*.  
   [https://arxiv.org/abs/1908.10084](https://arxiv.org/abs/1908.10084)

4. Johnson, Jeff, Matthijs Douze, and Hervé Jégou. "Billion-scale similarity search with GPUs." *IEEE Transactions on Big Data* 7.3 (2021): 535-547.  
   [https://arxiv.org/abs/1702.08734](https://arxiv.org/abs/1702.08734)

5. PubMed: [https://pubmed.ncbi.nlm.nih.gov/](https://pubmed.ncbi.nlm.nih.gov/)


---

## 🛠️ Technology References

- [FastAPI](https://fastapi.tiangolo.com/)  
  Sebastián Ramírez. "FastAPI framework, high performance, easy to learn, fast to code, ready for production."
- [SentenceTransformers](https://www.sbert.net/)
- [FAISS (Facebook AI Similarity Search)](https://github.com/facebookresearch/faiss)
- [BioGPT (Microsoft)](https://huggingface.co/microsoft/BioGPT)
- [NLTK (Natural Language Toolkit)](https://www.nltk.org/)
- [psycopg2 (PostgreSQL adapter for Python)](https://www.psycopg.org/)
- [Tkinter (Python standard GUI library)](https://docs.python.org/3/library/tkinter.html)
- [Docker](https://www.docker.com/)
