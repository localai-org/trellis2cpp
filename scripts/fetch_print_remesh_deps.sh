#!/usr/bin/env bash
# Fetch the exact header-only dependency set used for CGAL print remeshing.
# Keep the four pins below together: the scheduled dependency workflow updates
# versions and upstream-published SHA-256 digests atomically.
set -euo pipefail

CGAL_VERSION=6.2.1
CGAL_SHA256=db49ffa9d2c5a3386e0725983bf3478e3fd2ade43553f8fc3900e18554decd49
BOOST_VERSION=1.92.0
BOOST_SHA256=c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce
FETCH_FORMAT=2

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DEST=${1:-"$ROOT/.deps/print-remesh"}
MARKER="$DEST/.trellis2-print-remesh-deps"
EXPECTED=$(printf 'CGAL_VERSION=%s\nCGAL_SHA256=%s\nBOOST_VERSION=%s\nBOOST_SHA256=%s\nFETCH_FORMAT=%s' \
	"$CGAL_VERSION" "$CGAL_SHA256" "$BOOST_VERSION" "$BOOST_SHA256" "$FETCH_FORMAT")

is_current() {
	[[ -f "$MARKER" && "$(<"$MARKER")" == "$EXPECTED" ]]
}

if is_current; then
	echo "Print-remesh dependencies already current in $DEST"
	exit 0
fi

mkdir -p "$(dirname "$DEST")"
LOCK="${DEST}.lock"
attempt=0
until mkdir "$LOCK" 2>/dev/null; do
	if is_current; then
		echo "Print-remesh dependencies already current in $DEST"
		exit 0
	fi
	((attempt += 1))
	if ((attempt >= 300)); then
		echo "Timed out waiting for dependency lock $LOCK" >&2
		exit 1
	fi
	sleep 0.2
done

TMP=
cleanup() {
	rm -rf "$LOCK"
	if [[ -n "${TMP:-}" ]]; then
		rm -rf "$TMP"
	fi
}
trap cleanup EXIT
TMP=$(mktemp -d "${DEST}.tmp.XXXXXX")

verify_sha256() {
	local file=$1 expected=$2 actual
	if command -v sha256sum >/dev/null 2>&1; then
		actual=$(sha256sum "$file")
	elif command -v shasum >/dev/null 2>&1; then
		actual=$(shasum -a 256 "$file")
	else
		echo "No SHA-256 tool found (need sha256sum or shasum)" >&2
		exit 127
	fi
	actual=${actual%% *}
	if [[ "$actual" != "$expected" ]]; then
		echo "SHA-256 mismatch for $file: expected $expected, got $actual" >&2
		exit 1
	fi
}

CGAL_ARCHIVE="$TMP/cgal.zip"
curl -fL --retry 3 \
	"https://github.com/CGAL/cgal/releases/download/v${CGAL_VERSION}/CGAL-${CGAL_VERSION}-library.zip" \
	-o "$CGAL_ARCHIVE"
verify_sha256 "$CGAL_ARCHIVE" "$CGAL_SHA256"
mkdir "$TMP/cgal-unpack"
unzip -q "$CGAL_ARCHIVE" -d "$TMP/cgal-unpack"
shopt -s nullglob
cgal_roots=("$TMP/cgal-unpack"/*)
shopt -u nullglob
if [[ ${#cgal_roots[@]} -ne 1 || ! -d "${cgal_roots[0]}" ]]; then
	echo "Expected one root directory in the CGAL archive" >&2
	exit 1
fi
mv "${cgal_roots[0]}" "$TMP/cgal"

BOOST_UNDERSCORE=${BOOST_VERSION//./_}
BOOST_ARCHIVE="$TMP/boost.tgz"
curl -fL --retry 3 \
	"https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_UNDERSCORE}.tar.gz" \
	-o "$BOOST_ARCHIVE"
verify_sha256 "$BOOST_ARCHIVE" "$BOOST_SHA256"
mkdir "$TMP/boost"
tar -xzf "$BOOST_ARCHIVE" -C "$TMP/boost" --strip-components=1 \
	"boost_${BOOST_UNDERSCORE}/boost"

# Raw Boost release archives intentionally have no installed CMake package.
# Supply the header-only targets CGAL consumes so config-mode find_package
# works on CMake 3.30+, where the legacy FindBoost module was removed.
mkdir "$TMP/boost-cmake"
cat > "$TMP/boost-cmake/BoostConfig.cmake" <<EOF
set(Boost_FOUND TRUE)
set(Boost_VERSION "$BOOST_VERSION")
set(Boost_VERSION_STRING "$BOOST_VERSION")
set(Boost_INCLUDE_DIR "\${CMAKE_CURRENT_LIST_DIR}/../boost")
set(Boost_INCLUDE_DIRS "\${Boost_INCLUDE_DIR}")
set(Boost_LIBRARIES "")
if(NOT TARGET Boost::headers)
    add_library(Boost::headers INTERFACE IMPORTED GLOBAL)
    set_target_properties(Boost::headers PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "\${Boost_INCLUDE_DIR}")
endif()
if(NOT TARGET Boost::boost)
    add_library(Boost::boost ALIAS Boost::headers)
endif()
EOF
cat > "$TMP/boost-cmake/BoostConfigVersion.cmake" <<EOF
set(PACKAGE_VERSION "$BOOST_VERSION")
if(PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
    set(PACKAGE_VERSION_COMPATIBLE TRUE)
    if(PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)
        set(PACKAGE_VERSION_EXACT TRUE)
    endif()
endif()
EOF

rm -rf "$CGAL_ARCHIVE" "$BOOST_ARCHIVE" "$TMP/cgal-unpack"
printf '%s\n' "$EXPECTED" > "$TMP/.trellis2-print-remesh-deps"

OLD="${DEST}.old.$$"
if [[ -e "$DEST" ]]; then
	mv "$DEST" "$OLD"
fi
if mv "$TMP" "$DEST"; then
	TMP=
	rm -rf "$OLD"
else
	if [[ -e "$OLD" ]]; then
		mv "$OLD" "$DEST"
	fi
	exit 1
fi

echo "Fetched CGAL $CGAL_VERSION and Boost $BOOST_VERSION into $DEST"
