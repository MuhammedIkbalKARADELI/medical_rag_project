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

- Each user query is recorded in the **rag_log_queries** table with:
  - The query text
  - Query length (short/medium/long)
  - Model-generated answer
  - Document retrieval time (ms)
  - Answer generation time (ms)
  - Total response time (ms)
  - Query timestamp

- Every PubMed document (or sentence) used for each query is separately stored in the **query_documents** table, **linked to the log entry**:
  - Query ID (foreign key to log table)
  - PubMed ID (PMID)
  - Document/sentence index

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

