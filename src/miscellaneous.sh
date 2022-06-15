#!/bin/bash

# This script contains miscellaneous helper functions for proper presentation and error handling.
# Author: Louis She (2022-04)
# Contact: louisshe@hku.hk

# Helper functions
function timestamp () {
    echo -e "[$(date +%a) $(date +%b-'%-m') $(date +'%-I':%M:%S%P)]"
}

# Alias
alias time="/usr/bin/time -f '$(timestamp) INFO: Elapsed time: %es\tMemory usage: %M KB\tCPU usage: %P'"
