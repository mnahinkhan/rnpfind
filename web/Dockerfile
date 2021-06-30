FROM python:3.8-slim
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1
ENV DEBUG 0

# used by bedToBigBed tool
RUN apt-get update && apt-get install --no-install-recommends -y \
    libkrb5-dev \
 && rm -rf /var/lib/apt/lists/*


# awscli has conflicts when included in requirements.txt, so these flags are
# used, I think ...
RUN pip install awscli --force-reinstall --upgrade

WORKDIR /app
COPY ./requirements.txt .
RUN pip install -r requirements.txt
COPY . .

RUN adduser --disabled-password myuser; \
    chown myuser /app


USER myuser
ENTRYPOINT ["./web-entrypoint.sh"]
