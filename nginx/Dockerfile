FROM nginx:1.19.0-alpine

ARG CONF_FILE=nginx-prod.conf
RUN rm /etc/nginx/conf.d/default.conf
COPY $CONF_FILE /etc/nginx/conf.d/nginx.conf.template
COPY nginx-entry.sh /nginx-entry.sh

ENTRYPOINT ["/nginx-entry.sh"]
