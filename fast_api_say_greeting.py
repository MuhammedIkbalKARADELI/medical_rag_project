from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def sayGreeting():
    return "Hello World"



