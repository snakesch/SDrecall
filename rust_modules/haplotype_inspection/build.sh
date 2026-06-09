#!/bin/bash
# Build script for haplotype_inspection module
# Must be run within SDrecall conda environment
#
# Produces a manylinux_2_28 wheel (compatible with Rocky 8+, Ubuntu 18.04+)
# The key trick: temporarily hide conda's libstdc++.so during linking so the
# linker resolves against the system's libstdc++ (GCC 8, GLIBCXX_3.4.22 max).
# Without this, conda's GCC 15 libstdc++ assigns GLIBCXX_3.4.30 version tags,
# which restricts the wheel to Ubuntu 22.04+ only.

set -e

# Check if we're in the SDrecall environment
if [[ "$CONDA_DEFAULT_ENV" != "SDrecall" ]]; then
    echo "Error: Must be run in SDrecall conda environment"
    echo "Run: conda activate SDrecall"
    exit 1
fi

# Set environment variables for building
export LIBCLANG_PATH=$CONDA_PREFIX/lib
export OPENSSL_NO_VENDOR=1
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH

# Unset VIRTUAL_ENV if set (maturin errors when both VIRTUAL_ENV and CONDA_PREFIX are set)
unset VIRTUAL_ENV

echo "Building haplotype_inspection Rust module in SDrecall environment..."
echo "LIBCLANG_PATH: $LIBCLANG_PATH"
echo "CONDA_PREFIX: $CONDA_PREFIX"

MODE="${1:-wheel}"  # "wheel" (default) or "develop"

if [[ "$MODE" == "develop" ]]; then
    echo "Mode: develop (install directly into current env)"
    maturin develop --release
    echo ""
    echo "Module installed directly into $(python -c 'import sys; print(sys.prefix)')"
else
    echo "Mode: wheel (manylinux_2_28 portable wheel)"

    # Use system linker to avoid conda's libstdc++ version tags
    export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=/usr/bin/cc

    # Temporarily hide conda's libstdc++.so so linker finds system version
    CONDA_STDCPP="$CONDA_PREFIX/lib/libstdc++.so.6"
    CONDA_STDCPP_LINK="$CONDA_PREFIX/lib/libstdc++.so"
    HIDDEN=false

    cleanup() {
        if $HIDDEN; then
            mv "${CONDA_STDCPP}.HIDDEN" "$CONDA_STDCPP" 2>/dev/null
            mv "${CONDA_STDCPP_LINK}.HIDDEN" "$CONDA_STDCPP_LINK" 2>/dev/null
            echo "Restored conda libstdc++.so"
        fi
    }
    trap cleanup EXIT

    mv "$CONDA_STDCPP" "${CONDA_STDCPP}.HIDDEN"
    mv "$CONDA_STDCPP_LINK" "${CONDA_STDCPP_LINK}.HIDDEN" 2>/dev/null
    HIDDEN=true
    echo "Temporarily hidden conda libstdc++ to force system version during linking"

    maturin build --release --manylinux 2_28

    echo ""
    echo "Build complete! Wheel file is in target/wheels/"
    echo "To install: pip install target/wheels/haplotype_inspection-*.whl"
fi
