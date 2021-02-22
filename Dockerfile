FROM debian:latest
RUN apt update
RUN apt install -y python3 git python3-pip python3-numpy python3-scipy
RUN ln -s `which python3` /bin/python
RUN ln -s `which pip3` /bin/pip
RUN mkdir src; cd src; git clone https://github.com/yhoogstrate/dr-disco.git ; cd dr-disco
RUN pip install -r requirements.txt
RUN python setup.py install
RUN nosetests tests/*.py

