This is the main repository of the [NOMAD](https://www.nomad-coe.eu/) parser for
[AMS-BAND](https://www.scm.com/).

# Setup and run example

We are currently targeting Python 3.6

Best use a virtual environment:
```
virtualenv -p python3 .pyenv
source .pyenv/bin/activate
```

Clone and install the nomad infrastructure and the necessary dependencies (including this parser)
```
git clone https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR nomad
git submodule update --init
pip install -r requirements.txt
./dependencies.sh -e
```

To run the parser:
```
cd nomad/dependencies/parser/band
python bandparser/parser_band.py test/new-18.105/<file.out>
```
