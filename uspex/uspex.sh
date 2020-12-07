#!/bin/bash
date >> uspex.log
USPEX -r 2>&1 | tee -a uspex.log
