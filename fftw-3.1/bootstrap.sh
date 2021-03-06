#! /bin/sh
############################################################################
#
# NOTE: If you just want to build FFTW, do not use this file.  Just use
# the ordinary ./configure && make commmands as described in the installation
# section of the manual.
#
# This file is only for users that want to generate their own codelets,
# as described in the "generating your own code" section of the manual.
#
############################################################################

touch ChangeLog
touch genfft/.depend
touch genfft-k7/.depend

echo "PLEASE IGNORE WARNINGS AND ERRORS"

# paranoia: sometimes autoconf doesn't get things right the first time
rm -rf autom4te.cache
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force

rm -f config.cache

# --enable-maintainer-mode enables build of genfft and automatic
# rebuild of codelets whenever genfft changes
(
    ./configure --disable-shared --enable-maintainer-mode --enable-threads $*
    (cd genfft; make depend)
    (cd genfft-k7; make depend)
)
