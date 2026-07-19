#!/usr/bin/env bash
# Assert container portability: every ACTIVE process references either a
# digest-pinned public image (@sha256:) or an ampliresist-owned image via
# ${params.container_registry}. Fails on bare mutable tags or local-only names
# (drag1-*, haema-*). The unused alternative callers (medaka/freebayes) are
# skipped — they are inherited from nano-rave and not on the validated path.
#
#   tests/test_containers.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
fail=0

while IFS= read -r line; do
    ref="$(sed -E 's/.*container +"([^"]+)".*/\1/' <<<"$line")"
    case "$ref" in *medaka*|*freebayes*) continue ;; esac      # unused alt callers
    case "$ref" in
        *'${params.container_registry}'*) : ;;                 # ampliresist-owned
        *@sha256:*) : ;;                                       # digest-pinned public
        *) echo "  UNPINNED / not owned: $ref"; fail=1 ;;
    esac
    case "$ref" in
        *drag1-*|*haema-*) echo "  LOCAL/FOREIGN image name: $ref"; fail=1 ;;
    esac
done < <(grep -hE '^[[:space:]]*container[[:space:]]' "$HERE"/modules/local/*.nf)

# No stale references to other pipelines' assets anywhere in the tracked source.
if grep -RInE '02_plasmodium_qpcr|/pipeline/data/geo|haema-figures|drag1-downstream' \
        "$HERE"/modules "$HERE"/subworkflows "$HERE"/conf "$HERE"/main.nf "$HERE"/nextflow.config 2>/dev/null; then
    echo "  STALE cross-pipeline reference above"; fail=1
fi

[ "$fail" -eq 0 ] && echo "containers: all active references pinned-or-owned; no stale cross-pipeline refs ✓" \
                  || echo "containers: issues above"
exit "$fail"
