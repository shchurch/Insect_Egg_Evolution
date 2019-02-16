#!/bin/bash
for infile in *.nex; do
    sed -i "s/\[mcmcp append=yes;\]/mcmcp append=yes;/g" $infile
done
