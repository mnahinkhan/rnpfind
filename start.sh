#!/usr/bin/env bash

sudo docker run -p 80:8000 -e PORT=8000 -e AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" \
-e AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" --name rnpfind rnpfind-test
