#!/bin/bash
mono $(dirname "$0")/MzidToTsv/MzidToTsvConverter.exe ${@:1}
