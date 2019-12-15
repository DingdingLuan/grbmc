#!/bin/bash
a=`tail -n 1 const_relog.txt`
b=`tail -n 1 gau_relog.txt`
c=`tail -n 1 pl_relog.txt`
echo -e "\e[1;31mConst\e[0m":$a ;echo  -e "\e[1;32mGauss\e[0m":$b ; echo -e "\e[1;34mBpl\e[0m":$c
