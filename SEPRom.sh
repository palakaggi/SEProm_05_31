#pip install virtualenv

#virtualenv --python=python3.8 SEProm_venv

source SEProm_venv/bin/activate

#pip install -r requirements.txt

python readAndCreateDF.py

python dataset_ML.py

deactivate
