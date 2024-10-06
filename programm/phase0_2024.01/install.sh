#!/bin/sh
# PHASE installer
# Copyright (C) 2007 Takenori YAMAMOTO

echo " === PHASE installer ==="
while [ 0 -eq 0 ]
do

echo " Do you want to install PHASE? (yes/no) [yes]"
read reply
if [ -z "${reply}" ]; then
reply="yes"
fi

# For SunOS
system=`uname`
if [ ${system} = 'SunOS' ]; then
PATH=/usr/xpg4/bin:${PATH}
cd bin
cat << 'here' > make
#!/bin/sh
dmake -m serial $*
here
chmod +x make
PATH=`pwd`:${PATH}
cd ..
fi

case ${reply} in
"yes")
cd src_phase
/bin/bash configure
if [ ${?} -ne 0 ]; then
exit 1
fi
break
	;;
"no"|"exit")
exit 1
	;;
esac

done

while [ 0 -eq 0 ]
do

echo ""
echo " Do you want to edit the makefile that has been generated? (yes/no/exit) [no]"
read reply
if [ -z "${reply}" ]; then
reply="no"
fi

case ${reply} in
"no")
break
	;;
"exit")
exit 1
	;;
"yes")
while [ 0 -eq 0 ]
do

echo "Enter an editor that will be used to edit the makefile. [vi]"
read editor
if [ -z "${editor}" ]; then
editor="vi"
fi
${editor} Makefile
if [ ${?} -ne 0 ]; then
echo "Can ${editor} be used on the system?"  
continue
else
break
fi

done
	;;
*)
continue
	;;
esac
break

done


while [ 0 -eq 0 ]
do

echo " Do you want to make PHASE now? (yes/no) [yes]"
read reply
if [ -z "${reply}" ]; then
reply="yes"
fi

case ${reply} in
"exit"|"no")
echo " PHASE has not been made, but using the make command in the src directory you can make PHASE later."
exit 0
	;;
"yes")
make all
if [ ${?} -ne 0 ]; then
exit 1
else
break
fi
	;;
esac

done

if [ ! -x phase ]; then
echo "Executable file 'phase' is not found. Probably, compilation failed."
exit 1
elif [ ! -x ekcal ]; then
echo "Executable file 'ekcal' is not found. Probably, compilation failed."
exit 1
else
make install
if [ ${?} -ne 0 ]; then
exit 1
fi
make clean
if [ ${?} -ne 0 ]; then
exit 1
fi
echo "PHASE was successfully installed."
fi

cd ../test
while [ 0 -eq 0 ]
do

echo "Do you want to check the installed programs? (yes/no) [no]"
read reply
if [ -z "${reply}" ]; then
reply="no"
fi

case ${reply} in
"yes")
./Run.sh
break
	;;
"no")
exit 0
	;;
esac

done
exit 0
