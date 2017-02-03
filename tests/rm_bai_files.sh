#!/bin/bash

bai_files=$(find tests | grep -v samples | grep bai$)
rm -f $bai_files
