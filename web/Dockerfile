FROM python:3.8-slim
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV DEBUG 0

# used by bedToBigBed tool
RUN apt-get update && apt-get install --no-install-recommends -y \
    libkrb5-dev \
    libpq-dev \
    gcc \
 && pip install psycopg2 \
 && rm -rf /var/lib/apt/lists/*


WORKDIR /app
COPY ./requirements.txt .
RUN pip install -r requirements.txt
COPY . .


RUN adduser --disabled-password myuser; \
    chown myuser /usr/local/lib/python3.8/site-packages/rnpfind; \
    chown myuser /usr/local/lib/python3.8/site-packages/hgfind; \
    chown myuser /app/website/migrations


USER myuser
ENTRYPOINT ["./web-entrypoint.sh"]
