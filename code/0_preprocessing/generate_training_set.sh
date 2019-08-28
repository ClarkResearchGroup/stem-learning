#!/bin/bash
# run this code with ./generate_training_set.sh path_to_data_directory
# In the data directory, there must be a file called "input." + $(ftype) that is the original stem image
# and the rest of the files are label iamge with names "label_" + $(label type) + "." + $(ftype)
# this script creates a folder called augments that resides in the data directory
# that includes folders of various manipulations on the data images

curdir=$PWD
cd $1
augdir=augments
[ -d $augdir ] && rm -rf $augdir
mkdir $augdir

ftype=$2 #tif or png
echo inverting image
# invert (2)
mkdir $augdir/flip
mkdir $augdir/orig
for f in $(ls *.$ftype)
do
    convert       $f $augdir/orig/$f
    convert -flop $f $augdir/flip/$f
done

echo rotating
# do rotation (4)
for dirpath in $(ls -d $augdir/*)
do
    dn=${dirpath##*/}
    mkdir $augdir/rot0_$dn
    mkdir $augdir/rot1_$dn
    mkdir $augdir/rot2_$dn
    mkdir $augdir/rot3_$dn
    for f in $(ls $dirpath/*)
    do
        fn=${f##*/}
        cp                   $f $augdir/rot0_$dn/$fn
        convert -rotate 90   $f $augdir/rot1_$dn/$fn
        convert -rotate 180  $f $augdir/rot2_$dn/$fn
        convert -rotate 270  $f $augdir/rot3_$dn/$fn
    done
    rm -rf $dirpath
done

echo down-up sampling
# do magnification (3)
for dirpath in $(ls -d $augdir/*)
do
    dn=${dirpath##*/}
    mkdir $augdir/mag0_$dn
    mkdir $augdir/mag1_$dn
    mkdir $augdir/mag2_$dn

    for f in $(ls $dirpath/*)
    do
        fn=${f##*/}
        cp $f $augdir/mag0_$dn/$fn
        cp $f $augdir/mag1_$dn/$fn
        cp $f $augdir/mag2_$dn/$fn
    done

    f=$dirpath/input.$ftype
    l=$(identify -format "%w" $f)> /dev/null
    l2=$((l/2))
    l4=$((l/4))
    convert -resize ${l2}x${l2}             $f   tmp.$ftype
    convert -resize ${l}x${l} -filter cubic tmp.$ftype $augdir/mag1_$dn/input.$ftype
    convert -resize ${l4}x${l4}             $f   tmp.$ftype
    convert -resize ${l}x${l} -filter cubic tmp.$ftype $augdir/mag2_$dn/input.$ftype
    rm tmp.$ftype

    rm -rf $dirpath
done



cd $curdir
