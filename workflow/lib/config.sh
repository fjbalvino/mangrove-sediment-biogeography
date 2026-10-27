#!/usr/bin/env bash
# Source from a workflow after setting PACKAGE_ROOT and parsing --config.
# Local settings take precedence over environment-preserving defaults.
magenta_load_config() {
  local requested="${1:-${MAGENTA_CONFIG:-$PACKAGE_ROOT/config/project.env}}"
  if [[ ! -f "$requested" ]]; then
    if [[ -n "${1:-}${MAGENTA_CONFIG:-}" ]]; then
      printf 'Configuration file not found: %s\n' "$requested" >&2
      return 2
    fi
    requested="$PACKAGE_ROOT/config/project.env.example"
  fi
  # Read local choices before deriving defaults from them.
  if [[ "$requested" != "$PACKAGE_ROOT/config/project.env.example" ]]; then
    source "$requested" || return
  fi
  # shellcheck disable=SC1091
  source "$PACKAGE_ROOT/config/project.env.example" || return
}
