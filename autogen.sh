#!/bin/sh
aclocal && \
autoheader && \
automake -a && \
autoconf && \
echo "Now run ./configure --prefix=$HOME ; make install"

