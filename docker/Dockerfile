FROM debian:latest
RUN apt-get update
RUN apt-get install -y python3 git python3-pip python3-numpy python3-scipy zlib1g-dev
RUN ln -s `which python3` /bin/python
RUN ln -s `which pip3` /bin/pip
RUN cd /usr/local/src ; git clone https://github.com/alexdobin/STAR.git STAR ; cd STAR ; git checkout -b star_2_7_7a tags/2.7.7a ; cd source ; make ; make install ; cp STAR /usr/local/bin
RUN cd /usr/local/src ; git clone https://github.com/yhoogstrate/dr-disco.git ; cd dr-disco ; pip install -r requirements.txt ; python setup.py install ; nosetests tests/*.py
RUN mkdir /data

