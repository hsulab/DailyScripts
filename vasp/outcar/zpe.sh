#!/bin/sh
grep "f  =" OUTCAR | awk '{a+=$10}END{printf "%22.16g\n", a/2000}'
