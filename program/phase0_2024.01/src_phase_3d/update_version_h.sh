#!/bin/sh
if [ -e ../.git ];then
echo -n 'character(len=40), parameter :: commit_id = ' > version.h ; git rev-parse --sq HEAD >> version.h ; echo >> version.h
else
echo "character(len=40), parameter :: commit_id = 'unknown'" > version.h
fi
