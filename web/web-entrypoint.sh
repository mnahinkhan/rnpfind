#!/usr/bin/env bash

# Collect static files
echo "Collect static files"
python manage.py collectstatic --noinput

# Apply database migrations
echo "Apply database migrations"
python manage.py makemigrations
python manage.py migrate

# Clean up bad database records
echo "Remove bad records"
echo 'import cleanup' | python manage.py shell

# Start server
echo "Starting server"
exec gunicorn rnp_find.wsgi:application --bind 0.0.0.0:"$PORT" -w 3 --timeout 60
