from fastapi import APIRouter, Request
from pydantic import BaseModel
# from app.retrieval import retrieve_documents
# from app.generation import generate_answer
# from app.db_logger import insert_query_documents, rag_log_query
# from app.adjunct import classify_query_length
import time
from retrieval import retrieve_documents
from generation import generate_answer
from db_logger import insert_query_documents, rag_log_query
from adjunct import classify_query_length

from nltk.translate.bleu_score import sentence_bleu, SmoothingFunction
from rouge import Rouge
import bert_score

router = APIRouter()

@router.get("/status")
async def status():
    return {"Message" : "OK"}



# ---- Request Modeli ----
class QueryRequest(BaseModel):
    question: str

@router.post("/query")
def query_endpoint(request: QueryRequest):
    start_time = time.time()
    # 1. Dökümanları getir
    retrieved_docs, retrieved_pmid = retrieve_documents(request.question)
    retrieval_time = int((time.time() - start_time) * 1000)

    if not retrieved_docs:
        return {"answer": "There is no any related article, give me a different question..."}

    # 2. Cevap oluştur
    gen_start = time.time()
    context = "\n".join(retrieved_docs)
    prompt = f"{context}\n\nQuestion: {request.question}\nAnswer:"
    answer = generate_answer(prompt)
    generation_time = int((time.time() - gen_start) * 1000)
    total_time = int((time.time() - start_time) * 1000)
    size_of_query = classify_query_length(request.question)


    try:
        # Burada, referans cevapları soruya göre eşleştir (örnek için ilk referansı alıyoruz)
        ref_answer = retrieved_docs[0]  # (Düzenle: doğru referansı seçmek için)
    except Exception:
        ref_answer = None

    bleu_score = None
    rouge_l_score = None
    bertscore_f1 = None

    if ref_answer:
        # BLEU
        smoother = SmoothingFunction().method2
        bleu_score = sentence_bleu([ref_answer.split()], answer.split(), smoothing_function=smoother)

        # ROUGE-L
        try:
            rouge = Rouge()
            rouge_scores = rouge.get_scores(answer, ref_answer)
            rouge_l_score = rouge_scores[0]["rouge-l"]["f"]
        except Exception:
            rouge_l_score = None

        # BERTScore (batch ile hızlı çalışmazsa, tek örnek için)
        try:
            P, R, F1 = bert_score.score([answer], [ref_answer], lang="en")
            bertscore_f1 = float(F1[0])
        except Exception:
            bertscore_f1 = None



    # Sorgu sonrası loglama:
    rag_query_id = rag_log_query(
            question=request.question,
            answer=answer,
            query_length=size_of_query,
            retrieval_time=retrieval_time,
            generation_time=generation_time,
            total_time=total_time,
            bleu_score=bleu_score,
            rouge_l_score=rouge_l_score,
            bertscore_f1=bertscore_f1
        )
    insert_query_documents(rag_query_id, retrieved_pmid)


    return {
        "question": request.question,
        "top_documents": retrieved_docs,
        "PubMed_ID" : retrieved_pmid,     
        "answer": answer,
        "Query_Length": size_of_query,  
        "Retrieval_Time_MS": retrieval_time,
        "Generation_Time_MS": generation_time,
        "Total_Time_MS": total_time
    }
