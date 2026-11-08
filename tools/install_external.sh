#!/usr/bin/env bash
# Download and build the third-party programs into external/ (git-ignored; nothing of them is stored in this
# repository).  Fixed versions, so that every installation runs the same code:
#   DSurfTomo   https://github.com/HongjianFang/DSurfTomo   a651faf (stable)   MIT      inversion, and the 1-D forward fwd1d
#   pDSurfTomo  https://github.com/zshang825/pDSurfTomo   9956bbb + patch     GPL-3    optional parallel solver (Mode 2)
#   CC-FJpy     https://github.com/ColinLii/CC-FJpy       e79a4bc (CPU)       LGPL-3   F-J reference curve
#
#   bash tools/install_external.sh                 all three
#   bash tools/install_external.sh dsurftomo       DSurfTomo + fwd1d only (enough when the reference curve is borrowed)
#   bash tools/install_external.sh pdsurftomo ccfj
# Needs git, make, gfortran (DSurfTomo, pDSurfTomo), a C/C++ compiler, Cython and the single-precision FFTW library
# (CC-FJpy; in the conda environment of environment.yml these are all present).
# The patch tools/patches/pdsurftomo-lowmem.patch only lowers the initial memory reserve per travel time
# (0.25 -> 0.05 of the model size); the matrix written, and so the model, is the same.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
EXT=${LANT_EXTERNAL:-$ROOT/external}
mkdir -p "$EXT/bin"
WHAT=${*:-dsurftomo pdsurftomo ccfj}

fetch() {   # url commit dir
  if [ ! -d "$3/.git" ]; then git clone -q "$1" "$3"; fi
  git -C "$3" fetch -q origin "$2" 2>/dev/null || true
  git -C "$3" checkout -q -f "$2"
}

for w in $WHAT; do
  case $w in
  dsurftomo)
    echo "== DSurfTomo"
    fetch https://github.com/HongjianFang/DSurfTomo.git a651faf92758c1b4c5eb494a4865ef797b9a4b1c "$EXT/DSurfTomo"
    (cd "$EXT/DSurfTomo/src" && make -s clean >/dev/null 2>&1 || true; make -s)
    cp "$EXT/DSurfTomo/src/DSurfTomo" "$EXT/bin/DSurfTomo"
    echo "== fwd1d (1-D forward with DSurfTomo's own routines)"
    B=$(mktemp -d)
    F="-O -ffixed-line-length-none -ffloat-store -m64 -mcmodel=medium"
    S=$EXT/DSurfTomo/src
    (cd "$B" && gfortran $F -c "$S/surfdisp96.f" -o surfdisp96.o && gfortran -fopenmp $F -c "$S/CalSurfG.f90" -o CalSurfG.o \
      && gfortran $F -c "$S/gaussian.f90" -o gaussian.o && gfortran $F -c "$ROOT/fortran/fwd1d.f90" -o fwd1d.o \
      && gfortran -fopenmp fwd1d.o CalSurfG.o surfdisp96.o gaussian.o -o "$EXT/bin/fwd1d")
    rm -rf "$B"
    ;;
  pdsurftomo)
    echo "== pDSurfTomo"
    fetch https://github.com/zshang825/pDSurfTomo.git 9956bbb858a1a7c70d36cc9b4c8d5070acdd6a0a "$EXT/pDSurfTomo"
    git -C "$EXT/pDSurfTomo" apply --check "$ROOT/tools/patches/pdsurftomo-lowmem.patch" 2>/dev/null \
      && git -C "$EXT/pDSurfTomo" apply "$ROOT/tools/patches/pdsurftomo-lowmem.patch" || echo "   (patch already applied)"
    (cd "$EXT/pDSurfTomo/src_pDSurfTomo" && make -s clean >/dev/null 2>&1 || true; make -s)
    cp "$EXT/pDSurfTomo/src_pDSurfTomo/pDSurfTomo" "$EXT/bin/pDSurfTomo"
    ;;
  ccfj)
    echo "== CC-FJpy (CPU version)"
    fetch https://github.com/ColinLii/CC-FJpy.git e79a4bc713e93835043be12771c33af93187054e "$EXT/CC-FJpy"
    PREFIX=${CONDA_PREFIX:-}
    if [ -n "$PREFIX" ] && ls "$PREFIX"/lib/libfftw3f.* >/dev/null 2>&1; then
      export CFLAGS="-I$PREFIX/include ${CFLAGS:-}" LDFLAGS="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib ${LDFLAGS:-}"
    elif ! ldconfig -p 2>/dev/null | grep -q libfftw3f; then
      echo "single-precision FFTW (libfftw3f) not found: conda install -c conda-forge fftw, then rerun" >&2; exit 1
    fi
    (cd "$EXT/CC-FJpy" && python setupCPU.py build_ext --inplace > "$EXT/ccfj_build.log" 2>&1) \
      || { echo "CC-FJpy build failed, see $EXT/ccfj_build.log" >&2; exit 1; }
    python -c "import sys; sys.path.insert(0, '$EXT/CC-FJpy'); import ccfj; print('   ccfj imported from', ccfj.__file__)"
    ;;
  *) echo "unknown item $w (dsurftomo, pdsurftomo, ccfj)" >&2; exit 1 ;;
  esac
done
echo "done: $EXT/bin"; ls "$EXT/bin"
