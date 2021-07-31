FROM nginx:1.19.0-alpine

RUN rm /etc/nginx/conf.d/default.conf
COPY nginx.conf /etc/nginx/conf.d/nginx.conf.template
COPY nginx-entry.sh /nginx-entry.sh

ENTRYPOINT ["/nginx-entry.sh"]
