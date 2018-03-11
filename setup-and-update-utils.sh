#! /bin/bash

thisscript=$0
scriptname=$(basename $thisscript)
source_dir=$(dirname $thisscript)
dest_dir="./"
origin="../CppProjectUtils"

if [ -x $origin/$scriptname ]
then
  echo "Running $scriptname script in $origin"
  $origin/$scriptname
  exit 0
fi

# copy this script
rsync -a $thisscript $dest_dir/$scriptname
sed -i "s|../CppProjectUtils|$source_dir|" $dest_dir/$scriptname

for dir in cmake util-scripts
do
rsync -a $source_dir/$dir/ $dest_dir/$dir
done
