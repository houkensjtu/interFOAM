#!/bin/bash
#
# USAGE: ./installgdbOF
#
# Created By:
# J.M. Gimenez(1), S.Marquez Damian(2), N. Nigro(2)
#
# (1)FICH - UNL - Santa Fe, Argentina
# (2)CIMEC (INTEC-UNL-CONICET), Santa Fe, Argentina
#
#
# contact: jmarcelogimenez@gmail.com
#
# February 2010
#
set -e
echo "-----------------------------------------------------"
echo "Installing gdbOF..."
echo "-----------------------------------------------------"
echo "Making sure that you have all needed libraries"
echo "-----------------------------------------------------"
if [ "Debug" != "$WM_COMPILE_OPTION" ]
then
    echo "Error: The parameter \$WM_PROJECT_DIR: " $WM_COMPILE_OPTION " is not \"Debug\""
    echo "    Change this flag in \$WM_PROJECT_DIR/etc/bashrc to perform the instalation"
    exit 1
fi
#sudo apt-get install gawk python gdb
#sudo apt-get install gawk python gdb
#yum install gawk python gdb
#yum install gawk python gdb
pushd . 
echo "# gdbOF sources #" >> ~/.gdbinit
echo "source ~/OpenFOAM/gdbOF/gdbOF.gdb" >> ~/.gdbinit

echo "making structure"
if [ -d $WM_PROJECT_USER_DIR/gdbOF_apps ];
then
	rm -rf $WM_PROJECT_USER_DIR/gdbOF_apps
fi
if [ -d ~/OpenFOAM/gdbOF ];
then
	rm -rf ~/OpenFOAM/gdbOF
fi

mkdir $WM_PROJECT_USER_DIR/gdbOF_apps

cp gdbOF.tar.gz ~/OpenFOAM/.
cp findCell.tar.gz $WM_PROJECT_USER_DIR/gdbOF_apps/.
cp findFace.tar.gz $WM_PROJECT_USER_DIR/gdbOF_apps/.
cp findFacesCell.tar.gz $WM_PROJECT_USER_DIR/gdbOF_apps/.

echo "descompressing ..."
cd ~/OpenFOAM
tar xzvf gdbOF.tar.gz
rm gdbOF.tar.gz
cd $WM_PROJECT_USER_DIR/gdbOF_apps
tar -xvf findCell.tar.gz
tar -xvf findFace.tar.gz
tar -xvf findFacesCell.tar.gz
rm findCell.tar.gz
rm findFace.tar.gz
rm findFacesCell.tar.gz

echo "installing gdbOF_apps..."
cd findCell
wmake
cd ../findFace
wmake
cd ../findFacesCell
wmake

echo "-----------------------------------------------------"
echo "End Installation.."
echo "-----------------------------------------------------"
popd
set +e

