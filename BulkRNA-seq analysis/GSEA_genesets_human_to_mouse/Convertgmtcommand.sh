#!/bin/bash
cat <Geneset file name>.gmt | parallel -k -j 5 './ConvertingHuman_mouse_GMT.sh {}' > filename_mouse.gmt
