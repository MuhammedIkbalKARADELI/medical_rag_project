from fastapi import FastAPI
from api_routes import router

# uvicorn app.main:app --reload 
app = FastAPI(
    title="Medical RAG API",
    description="Retrieval-Augmented Generation (RAG) for Medical Q&A.",
    version="1.0.0"
)

app.include_router(router)

