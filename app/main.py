# from app.api_routes import router
from fastapi import FastAPI

app = FastAPI(title="Medical RAG API")

@app.get("/")
def sayGreeting():
    return "Hello World"


@app.get("/status")
async def status():
    return {"Message" : "OK"}

