import pytest
from fastapi.testclient import TestClient
from random import randint
from app.main import app


# Bu kodu test etmek için direkt terminale pytest yazman yeterlidir...

@pytest.fixture(scope="module")
def client():
    with TestClient(app) as c:
        yield c

def test_when_app_runing_status_endpoint_should_return_OK(client):
    response = client.get("/status")
    assert response.status_code == 200
    data = response.json()
    assert data["Message"] == "OK"

