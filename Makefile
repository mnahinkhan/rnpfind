# Makefile

SHELL := /bin/bash

test:
	source .env.test \
	&& envsubst < docker-compose.yml > docker-compose-test.yml \
	&& docker-compose -f docker-compose-test.yml pull \
	&& docker-compose -f docker-compose-test.yml up

test-swarm:
	source .env.test \
	&& envsubst < docker-compose.yml > docker-compose-test.yml \
	&& docker stack deploy -c docker-compose-test.yml test_swarm

prod:
	source .env.prod \
	&& envsubst < docker-compose.yml > docker-compose-prod.yml \
	&& docker-compose -f docker-compose-prod.yml pull \
	&& docker-compose -f docker-compose-prod.yml up

prod-swarm:
	source .env.prod \
	&& envsubst < docker-compose.yml > docker-compose-prod.yml \
	&& docker stack deploy -c docker-compose-prod.yml prod_swarm

dev:
	./preprocess_web_files.sh \
	&& docker-compose -f docker-compose-dev.yml build \
	&& docker-compose -f docker-compose-dev.yml up

push:
	./preprocess_web_files.sh \
	&& docker build --tag rnpfind/nginx nginx/ \
	&& docker build -f web/Dockerfile-celery --tag rnpfind/celery web/ \
	&& docker build --tag rnpfind/web web/ \
	&& docker push rnpfind/nginx \
	&& docker push rnpfind/web \
	&& docker push rnpfind/celery
