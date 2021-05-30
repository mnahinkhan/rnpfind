FROM python:3.8-slim AS build-python
COPY ./requirements.txt /
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r requirements.txt

FROM python:3.8-slim
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV DEBUG 0
RUN apt-get update
RUN apt-get install --no-install-recommends -y gcc python3-dev
RUN apt-get install --no-install-recommends -y postgresql
RUN apt-get install --no-install-recommends -y python3-psycopg2
RUN apt-get install --no-install-recommends -y libpq-dev
RUN pip install psycopg2
RUN apt-get install --no-install-recommends -y bash
RUN pip install awscli --force-reinstall --upgrade
RUN aws --version
COPY --from=build-python /wheels /wheels
COPY --from=build-python requirements.txt .
RUN pip install --no-cache /wheels/*
WORKDIR /app
COPY . .
RUN python manage.py collectstatic --noinput
RUN python manage.py makemigrations
RUN python manage.py migrate
RUN adduser --disabled-password myuser
RUN mkdir -p website/output-data/pickles
RUN chown myuser website/output-data/
RUN chown myuser website/output-data/pickles
RUN chown -R myuser website/data
RUN chown -R myuser website/ucsc-tools/linux
RUN chmod -R 755 website/ucsc-tools/linux
RUN chown myuser db.sqlite3 || true
RUN chown myuser /app
USER myuser
CMD gunicorn rnp_find.wsgi:application --bind 0.0.0.0:$PORT -w 2 --timeout 960
