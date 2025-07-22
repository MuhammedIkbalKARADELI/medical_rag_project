FROM python:3.13

WORKDIR /app

# Sistem gereksinimleri
RUN apt-get update && apt-get install -y build-essential
# Kodları kopyala
COPY app/ /app


# Gereksinimleri yükle
RUN pip install --no-cache-dir --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

# NLTK data (özellikle "punkt") otomatik indir
RUN python -c "import nltk; nltk.download('punkt'); nltk.download('punkt_tab'); nltk.download('wordnet'); nltk.download('omw-1.4')"

# FastAPI başlatıcı komutu
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]

