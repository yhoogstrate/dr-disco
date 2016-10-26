#!/bin/bash

if [[ "$PWD" =~ scripts$ ]] ; then
    cd .. ;
    py.test --cov-report html --cov=drdisco tests
else
    py.test --cov-report html --cov=drdisco tests
fi
