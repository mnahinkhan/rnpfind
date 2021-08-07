#!/usr/bin/env bash

export NGINX_LOCATIONS
NGINX_LOCATIONS="$(cat nginx/nginx-locations.conf)"

envsubst '${NGINX_LOCATIONS}' \
< nginx/nginx-prod-template.conf \
> nginx/nginx-prod.conf

envsubst '${NGINX_LOCATIONS}' \
< nginx/nginx-dev-template.conf \
> nginx/nginx-dev.conf
