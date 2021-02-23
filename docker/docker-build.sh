#!/bin/bash

docker build -t yhoogstrate/drdisco -f Dockerfile .

docker run -ti yhoogstrate/dr-disco
