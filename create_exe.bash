pyinstaller --onefile main.py Protein.py Operations.py FileMngr.py DomainBuilder.py Domain.py Clusterer.py
pyinstaller main.spec

pyinstaller --onefile --path=C:/Users/User/OneDrive/Documents/GitHub/DynDom-Py/venv/Lib/site-packages main.py Protein.py Operations.py FileMngr.py DomainBuilder.py Domain.py Clusterer.py
