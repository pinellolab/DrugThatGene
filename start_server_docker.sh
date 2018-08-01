#python app.py
gunicorn DrugThatGene:app --bind 0.0.0.0:9999 --timeout 1800 --access-logfile - --workers 2
