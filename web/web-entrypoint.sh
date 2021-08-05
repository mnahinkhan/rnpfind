#!/usr/bin/env bash

# Collect static files
echo "Collect static files"
python manage.py collectstatic --noinput

# Move images to static files"
echo "Moving images to static files"
cp -r /app/images /app/staticfiles

# Apply database migrations
echo "Apply database migrations"
python manage.py makemigrations
python manage.py migrate

# To allow testing
echo "Remove gene record"
echo 'import remove_record' | python manage.py shell

# Clean up bad database records
echo "Remove bad records"
echo 'import cleanup' | python manage.py shell

# Start server
echo "Starting server"
if [[ "$DEBUG" -eq 0 ]]; then
    exec gunicorn rnp_find.wsgi:application --bind 0.0.0.0:"$PORT" -w 3 \
    --timeout 60
else
    exec python manage.py runserver 0.0.0.0:"$PORT"
fi
