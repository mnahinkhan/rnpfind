# Generated by Django 3.1.1 on 2020-10-12 18:16

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('website', '0003_auto_20201012_1406'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='gene',
            name='reverse_name',
        ),
    ]
