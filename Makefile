# Makefile

SHELL := /bin/bash

test:
	source .env.test \
	&& envsubst < docker-compose.yml > docker-compose-test.yml \
	&& docker-compose -f docker-compose-test.yml pull \
	&& docker-compose -f docker-compose-test.yml up

prod:
	source .env.prod \
	&& envsubst < docker-compose.yml > docker-compose-prod.yml \
	&& docker-compose -f docker-compose-prod.yml pull \
	&& docker-compose -f docker-compose-prod.yml up
