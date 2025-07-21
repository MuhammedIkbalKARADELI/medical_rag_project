import requests

API_URL = "http://localhost:8000/query"
question = "What are the symptoms of anemia?"
r = requests.post(API_URL, json={"question": question})
data = r.json()
answer = data.get("answer", "There are not any answer.")
print(answer)