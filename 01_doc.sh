#!/bin/bash

sphinx-apidoc -f -o docs/source ./aquabreeding
sphinx-build ./docs/source ./docs/build
