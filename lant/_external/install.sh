#!/usr/bin/env bash
# Download and build the third-party programs into external/ (git-ignored; nothing of them is stored in this
# repository).  Fixed versions, so that every installation runs the same code:
#   DSurfTomo   https://github.com/HongjianFang/DSurfTomo   a651faf (as is)     MIT      inversion, and the 1-D forward fwd1d
#   pDSurfTomo  https://github.com/zshang825/pDSurfTomo   9956bbb + patch     GPL-3    optional parallel solver (Mode 2)
#
#   bash tools/install_external.sh                 both
#   bash tools/install_external.sh dsurftomo       DSurfTomo + fwd1d only
#   bash tools/install_external.sh pdsurftomo
# Needs git, make and gfortran (in the conda environment of environment.yml these are present).
# The patch lant/_external/pdsurftomo-lowmem.patch only lowers the initial memory reserve per travel time
# (0.25 -> 0.05 of the model size); the matrix written, and so the model, is the same.
# The synthetic noise of the checkerboard tests is drawn by DSurfTomo / pDSurfTomo themselves, without a fixed seed,
# so the recovery numbers of a test may vary slightly from run to run.
set -euo pipefail
RESOURCE_DIR=$(cd "$(dirname "$0")" && pwd)
EXT=${LANT_EXTERNAL:-${XDG_CACHE_HOME:-$HOME/.cache}/lant/external}
mkdir -p "$EXT/bin"
WHAT=${*:-dsurftomo pdsurftomo}

fetch() {   # url commit dir
  if [ ! -d "$3/.git" ]; then git clone -q "$1" "$3"; fi
  if ! git -C "$3" cat-file -e "$2^{commit}" 2>/dev/null; then
    git -C "$3" fetch -q origin "$2"
  fi
  git -C "$3" checkout -q -f "$2"
}

put_bin() {   # built-file name: copy then rename, so that a run still using the old program is not disturbed
  cp "$1" "$EXT/bin/$2.new.$$" && mv -f "$EXT/bin/$2.new.$$" "$EXT/bin/$2"
}

patch_in() {   # dir patch-file
  if git -C "$1" apply --check "$RESOURCE_DIR/$2" 2>/dev/null; then
    git -C "$1" apply "$RESOURCE_DIR/$2"
  elif git -C "$1" apply --reverse --check "$RESOURCE_DIR/$2" 2>/dev/null; then
    echo "   ($2 already applied)"
  else
    echo "$2 does not apply to $1" >&2; exit 1
  fi
}

for w in $WHAT; do
  case $w in
  dsurftomo)
    echo "== DSurfTomo"
    fetch https://github.com/HongjianFang/DSurfTomo.git a651faf92758c1b4c5eb494a4865ef797b9a4b1c "$EXT/DSurfTomo"
    (cd "$EXT/DSurfTomo/src" && make -s clean >/dev/null 2>&1 || true; make -s)
    put_bin "$EXT/DSurfTomo/src/DSurfTomo" DSurfTomo
    echo "== fwd1d (1-D forward with DSurfTomo's own routines)"
    B=$(mktemp -d)
    trap 'rm -rf "$B"' EXIT
    F="-O -ffixed-line-length-none -ffloat-store -m64 -mcmodel=medium"
    S=$EXT/DSurfTomo/src
    (cd "$B" && gfortran $F -c "$S/surfdisp96.f" -o surfdisp96.o && gfortran -fopenmp $F -c "$S/CalSurfG.f90" -o CalSurfG.o \
      && gfortran $F -c "$S/gaussian.f90" -o gaussian.o && gfortran $F -c "$RESOURCE_DIR/fwd1d.f90" -o fwd1d.o \
      && gfortran -fopenmp fwd1d.o CalSurfG.o surfdisp96.o gaussian.o -o fwd1d)
    put_bin "$B/fwd1d" fwd1d
    rm -rf "$B"
    ;;
  pdsurftomo)
    echo "== pDSurfTomo"
    fetch https://github.com/zshang825/pDSurfTomo.git 9956bbb858a1a7c70d36cc9b4c8d5070acdd6a0a "$EXT/pDSurfTomo"
    patch_in "$EXT/pDSurfTomo" pdsurftomo-lowmem.patch
    (cd "$EXT/pDSurfTomo/src_pDSurfTomo" && make -s clean >/dev/null 2>&1 || true; make -s)
    put_bin "$EXT/pDSurfTomo/src_pDSurfTomo/pDSurfTomo" pDSurfTomo
    ;;
  *) echo "unknown item $w (dsurftomo, pdsurftomo)" >&2; exit 1 ;;
  esac
done
echo "done: $EXT/bin"; ls "$EXT/bin"
