#!/usr/bin/env bash
# Resolve requested design; -p selects paired and conflicts with independent.
resolve_design() {
  local requested_design="$1"
  local pair_column="$2"
  case "$requested_design" in
    auto|independent) ;;
    *) echo 'ERROR: design must be auto or independent; use -p for paired.' >&2; return 2 ;;
  esac
  if [[ -n "$pair_column" ]]; then
    if [[ "$requested_design" == 'independent' ]]; then
      echo 'ERROR: -p and -d independent cannot be combined.' >&2
      return 2
    fi
    printf '%s\n' 'paired'
    return 0
  fi
  printf '%s\n' "$requested_design"
}
