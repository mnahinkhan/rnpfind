#!/usr/bin/env bash
# Preprocess web files
cp README.md web/about.md \
&& mkdir -p web/images \
&& cp cli/README.md web/cli_docs.md \
&& cp images/rnpfind-logo-transparent.png \
web/images/rnpfind-logo-transparent.png \
&& sass web/website/static/website/styles.scss \
web/website/static/website/styles.css --no-source-map
