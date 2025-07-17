# from fastapi import FastAPI
# from app.api_routes import router

# app = FastAPI(title="Medical RAG API")
# app.include_router(router)



from fastapi import FastAPI

app = FastAPI()

@app.get("/")
def sayGreeting():
    return "Hello World"



# from typing import Union

# from fastapi import FastAPI

# app = FastAPI()


# @app.get("/")
# def read_root():
#     return {"Hello": "World"}


# @app.get("/items/{item_id}")
# def read_item(item_id: int, q: Union[str, None] = None):
#     return {"item_id": item_id, "q": q}


