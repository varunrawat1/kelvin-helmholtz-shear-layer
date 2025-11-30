
#!/usr/bin/env bash

set -e
mkdir -p build
cd build
cmake .. >/dev/null

cmake --build . -j >/dev/null

cd ..

./build/kh_sim

