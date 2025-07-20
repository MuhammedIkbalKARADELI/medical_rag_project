FROM python:3.13

WORKDIR /app

# Gerekli sistem kütüphaneleri
RUN apt-get update && apt-get install -y build-essential
# Kodları kopyala
COPY app/ /app

# Gereksinimleri yükle
RUN pip install --no-cache-dir --upgrade pip
RUN pip install -r requirements.txt

# NLTK için verileri indir (eğer gerekiyorsa)
RUN python -c "import nltk; nltk.download('punkt')"

# API'yı başlat (uvicorn ile)
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]



# cd medical_rag_project
# docker compose build
# docker compose up

### or 

# docker build -t my-medical-rag .
# docker run -p 8000:8000 my-medical-rag