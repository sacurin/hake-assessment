#!/bin/bash

# Pull a repository. Used for a starting script for Docker container
# so that you are starting with the most current commit in the repo
# each time you start the container

git stash # There seem to be changes somehow. Ignore them.
git stash clear

git pull

# `git pull` resets the file permissions so we have to set them back to
# fully open again after pulling
chmod -R 777 .

# Return control back to the main Dockerfile script. If this is missing,
# then the ENTRYPOINT command that called this script will exit the
# container once this script is finished. Don't removeor the container will
# never open up for you.
exec "$@"
