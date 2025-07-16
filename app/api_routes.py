from fastapi import APIRouter, Request
from app.retrieval import retrieve_documents
from app.generation import generate_answer

router = APIRouter()

@router.post("/query")
async def query_api(request: Request):
    body = await request.json()
    question = body.get("question")

    retrieved_docs = retrieve_documents(question)
    answer = generate_answer(question, retrieved_docs)

    return {
        "question": question,
        "answer": answer["text"],
        "sources": answer["sources"],
        "latency_ms": answer["latency"]
    }