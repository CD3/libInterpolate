#! /bin/bash

tag=$1
shift

if [[ -z $tag ]]
then
  echo "ERROR: version number required"
  echo "example: $0 v1.0.1"
  echo "current tags:"
  git tag
  exit 1
fi

function exit_on_error()
{
  echo "There was an error. Commit will not be tagged."
}

trap exit_on_error ERR

set -e
root=$(git rev-parse --show-toplevel)

echo "cd'ing to root directory ($root)"
cd $root

echo "looking for pre-tag-release.sh to run"
$(find ./ -name 'pre-tag-release.sh')

echo "tagging with ${tag}"
git tag ${tag}
git tag | grep ${tag}
echo "Successfully tagged commit."
