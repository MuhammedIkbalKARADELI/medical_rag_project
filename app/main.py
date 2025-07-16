from fastapi import FastAPI
from app.api_routes import router

app = FastAPI(title="Medical RAG API")
app.include_router(router)

