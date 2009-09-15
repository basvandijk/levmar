#!/bin/env bash
cabal clean
find . -name "*.o" -o -name "*.hi" | xargs --no-run-if-empty rm
