#!/bin/bash
# Build script for the phasing graph Rust module

echo "Building build_phasing_graph Rust module..."

# Clean previous builds
cargo clean

# Build in release mode for maximum performance
maturin develop --release

if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo "The module 'build_phasing_graph_rs' is now available for import in Python"
else
    echo "Build failed!"
    exit 1
fi