FROM python:3.8-slim AS build-python
COPY ./requirements.txt /
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r requirements.txt

FROM python:3.8-slim
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV DEBUG 0

RUN apt-get update && apt-get install --no-install-recommends -y \
    bash \
    gcc \
    libpq-dev \
    postgresql \
    python3-dev \
    python3-psycopg2 \
 && rm -rf /var/lib/apt/lists/*

RUN pip install psycopg2
RUN pip install awscli --force-reinstall --upgrade
# RUN aws --version

COPY --from=build-python /wheels /wheels
COPY --from=build-python requirements.txt .

RUN pip install --no-cache /wheels/*

WORKDIR /app
COPY . .

RUN python manage.py collectstatic --noinput; \
    python manage.py makemigrations; \
    python manage.py migrate; \
    adduser --disabled-password myuser; \
    mkdir -p website/output-data/pickles; \
    chown myuser website/output-data/; \
    chown myuser website/output-data/pickles; \
    chown -R myuser website/data; \
    chown -R myuser website/ucsc-tools/linux; \
    chmod -R 755 website/ucsc-tools/linux; \
    chown myuser db.sqlite3 || true; \
    chown myuser /app

USER myuser
CMD gunicorn rnp_find.wsgi:application --bind 0.0.0.0:$PORT -w 2 --timeout 960
