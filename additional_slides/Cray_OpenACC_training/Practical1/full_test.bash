#!/bin/bash

#for MYPE in cray pgi
#do

    for TARGET in Fstatic Fdynamic Cstatic Cdynamic
    do

        for VERSION in 00 01 02
        do

            bash build_submit.bash $MYPE $TARGET $VERSION

        done
    done
#done

make clean
