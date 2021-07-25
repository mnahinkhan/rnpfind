#!/usr/bin/env sh
# Inspired by: https://serverfault.com/a/919212

envsubst '${UPSTREAM_HOST}' \
< /etc/nginx/conf.d/nginx.conf.template \
> /etc/nginx/conf.d/nginx.conf

exec "$@"
