#!/usr/bin/env bash
#
# How the source was downloaded.

set -euo pipefail
cd "$(dirname "$0")" || exit 1

VERSION=Version_4.0.1.0

curl -sL https://github.com/graph-algorithms/edge-addition-planarity-suite/archive/refs/tags/$VERSION.zip -o planarity-$VERSION.zip
unzip planarity-$VERSION.zip

[[ -d c ]] && rm -r c
mv edge-addition-planarity-suite-$VERSION/{c,LICENSE.TXT} .

rm -r c/{planarityApp,samples,.gdbinit}

rm -r edge-addition-planarity-suite-$VERSION
rm planarity-$VERSION.zip
