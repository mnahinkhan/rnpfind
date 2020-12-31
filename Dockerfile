FROM python:3.8.6-slim AS build-python
COPY ./requirements.txt /
RUN pip wheel --no-cache-dir --no-deps --wheel-dir /wheels -r requirements.txt



FROM python:3.8.6-slim
# RUN wget -q -O /etc/apk/keys/sgerrand.rsa.pub https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub
# RUN wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.32-r0/glibc-2.32-r0.apk
# RUN apk add glibc-2.32-r0.apk

# RUN wget "https://www.archlinux.org/packages/core/x86_64/zlib/download" -O /tmp/libz.tar.xz \
#     && mkdir -p /tmp/libz \
#     && tar -xf /tmp/libz.tar.xz -C /tmp/libz \
#     && cp /tmp/libz/usr/lib/libz.so.1.2.11 /usr/glibc-compat/lib \
#     && /usr/glibc-compat/sbin/ldconfig \
#     && rm -rf /tmp/libz /tmp/libz.tar.xz

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
# RUN apk add --no-cache libc6-compat
# RUN apk add libbsd
# RUN apk add libbsd-dev
# RUN apk add gcompat
# ENV LANG en_US.UTF-8
# ENV LANGUAGE en_US:en
# ENV LC_ALL en_US.UTF-8

# ENV GLIBC_REPO=https://github.com/sgerrand/alpine-pkg-glibc
# ENV GLIBC_VERSION=2.30-r0

# RUN set -ex && \
#     apk --update add libstdc++ curl ca-certificates && \
#     for pkg in glibc-${GLIBC_VERSION} glibc-bin-${GLIBC_VERSION}; \
#         do curl -sSL ${GLIBC_REPO}/releases/download/${GLIBC_VERSION}/${pkg}.apk -o /tmp/${pkg}.apk; done && \
#     apk add --allow-untrusted /tmp/*.apk && \
#     rm -v /tmp/*.apk && \
#     /usr/glibc-compat/sbin/ldconfig /lib /usr/glibc-compat/lib

COPY --from=build-python /wheels /wheels
COPY --from=build-python requirements.txt .
RUN pip install --no-cache /wheels/*
WORKDIR /app
COPY . .
RUN python manage.py collectstatic --noinput
RUN adduser --disabled-password myuser
RUN mkdir -p website/output-data/pickles
RUN chown myuser website/output-data/
RUN chown myuser website/output-data/pickles
RUN chown -R myuser website/data
RUN chown -R myuser website/ucsc-tools/linux
RUN chmod -R 755 website/ucsc-tools/linux
RUN chown myuser db.sqlite3 || true
RUN chown myuser /app
# RUN ./website/ucsc-tools/linux/bedToBigBed
USER myuser
CMD gunicorn rnp_find.wsgi:application --bind 0.0.0.0:$PORT -w 2 --timeout 960