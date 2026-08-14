#!/usr/bin/env bash
# Self-test battery for scripts/ci/new_domains_e2e.sh (bead rjoq.8).
#
# Proves the runner's own contract with fixture manifests in a scratch
# directory: argument refusals, manifest schema refusals, deterministic
# listing, tighten-only overrides, exit-class mapping for authority /
# refusal / production-failure / budget cases, summary+log agreement, and
# replay tamper detection. Uses only bash+python3; no cargo, no mocks of
# the runner itself (every assertion executes the real runner binary).
set -u -o pipefail

RUNNER="$(cd "$(dirname "$0")" && pwd)/new_domains_e2e.sh"
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SCRATCH="$REPO_ROOT/.e2e-out/selftest-$$"
mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

PASS=0
FAIL=0
check() { # description expected-exit actual-exit
  if [ "$2" -eq "$3" ]; then
    PASS=$((PASS + 1))
  else
    FAIL=$((FAIL + 1))
    echo "SELFTEST FAIL: $1 (expected exit $2, got $3)" >&2
  fi
}

fixture() { # name body...
  local path="$SCRATCH/$1"
  shift
  printf '%s\n' "$@" > "$path"
  echo "$path"
}

CASE_HEADER='schema = "frankensim.new-domains.case-manifest.v1"'

# --- usage refusals -------------------------------------------------------
bash "$RUNNER" --no-such-flag >/dev/null 2>&1
check "unknown flag refuses as usage" 10 $?
bash "$RUNNER" --phase nope --list >/dev/null 2>&1
check "unknown phase refuses as usage" 10 $?
bash "$RUNNER" --phase e1 --manifest x.toml >/dev/null 2>&1
check "phase+manifest are mutually exclusive" 10 $?
bash "$RUNNER" --seed not-a-number --list >/dev/null 2>&1
check "non-integer seed refuses" 10 $?
bash "$RUNNER" --manifest "$SCRATCH/absent.toml" --list >/dev/null 2>&1
check "missing manifest is unqualified infrastructure" 18 $?

# --- manifest schema refusals --------------------------------------------
M="$(fixture e1.toml 'schema = "wrong.schema.v9"' 'phase = "e1"')"
bash "$RUNNER" --manifest "$M" --validate-only >/dev/null 2>&1
check "wrong schema string refuses" 10 $?

M="$(fixture e1.toml "$CASE_HEADER" 'phase = "e2"')"
bash "$RUNNER" --manifest "$M" --validate-only >/dev/null 2>&1
check "phase/file-stem mismatch refuses" 10 $?

M="$(fixture e1.toml "$CASE_HEADER" 'phase = "e1"' 'undeclared_field = 1' '[[case]]' 'id = "a-case-x"')"
bash "$RUNNER" --manifest "$M" --validate-only >/dev/null 2>&1
check "unknown top-level field refuses (no silent semantics)" 10 $?

good_case() { # id extra-lines...
  local id="$1"; shift
  printf '%s\n' '[[case]]' "id = \"$id\"" 'version = 1' 'purpose = "p"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'gauntlet_tier = "G0"' \
    'seed = 7' 'max_wall_seconds = 60' "$@"
}

M="$SCRATCH/e1.toml"
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case case-dup 'entry_command = ["true"]' 'expected = "authority"'
  good_case case-dup 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
bash "$RUNNER" --manifest "$M" --validate-only >/dev/null 2>&1
check "duplicate semantic case id refuses" 10 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case esc-case 'entry_command = ["/bin/echo", "hi"]' 'expected = "authority"'
} > "$M"
bash "$RUNNER" --manifest "$M" --validate-only >/dev/null 2>&1
check "absolute command word refuses as path escape" 10 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case ref-case 'entry_command = ["false"]' 'expected = "refusal"'
} > "$M"
bash "$RUNNER" --manifest "$M" --validate-only >/dev/null 2>&1
check "refusal case without pattern refuses" 10 $?

# --- listing and selection ------------------------------------------------
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case list-b 'entry_command = ["true"]' 'expected = "authority"'
  good_case list-a 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
LISTED="$(bash "$RUNNER" --manifest "$M" --list 2>/dev/null | cut -f2 | paste -sd, -)"
[ "$LISTED" = "list-b,list-a" ]
check "listing preserves deterministic manifest order" 0 $?
bash "$RUNNER" --manifest "$M" --case no-such-case --list >/dev/null 2>&1
check "unknown --case selection refuses" 10 $?

# --- execution exit classes ----------------------------------------------
run_one() { # manifest extra-args...
  local manifest="$1"; shift
  bash "$RUNNER" --manifest "$manifest" --output-dir "$SCRATCH/out-$RANDOM$RANDOM" "$@" >/dev/null 2>&1
}

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case ok-auth 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
run_one "$M"
check "authority case passes with exit 0" 0 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case prod-fail 'entry_command = ["false"]' 'expected = "authority"'
} > "$M"
run_one "$M"
check "failing production maps to class 12" 12 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case ref-ok 'entry_command = ["grep", "--no-such-flag-zz"]' 'expected = "refusal"' \
    'expected_refusal_pattern = "unrecognized option|unknown option|usage"'
} > "$M"
run_one "$M"
check "matched refusal is exit 0" 0 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case ref-mismatch 'entry_command = ["true"]' 'expected = "refusal"' \
    'expected_refusal_pattern = "never-printed"'
} > "$M"
run_one "$M"
check "expected-refusal-but-authority maps to class 11" 11 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case checker-refuses 'entry_command = ["true"]' 'expected = "authority"' \
    'checker_command = ["false"]'
} > "$M"
run_one "$M"
check "independent checker refusal maps to class 13" 13 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  printf '%s\n' '[[case]]' 'id = "budget-case"' 'version = 1' 'purpose = "p"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'gauntlet_tier = "G0"' \
    'seed = 7' 'max_wall_seconds = 1' 'entry_command = ["sleep", "10"]' \
    'expected = "authority"'
} > "$M"
run_one "$M"
check "wall-budget exhaustion maps to class 14" 14 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case seed-locked 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
run_one "$M" --seed 99
check "seed override on a non-overridable case is unqualified (18)" 18 $?

# --- tighten-only budget override ----------------------------------------
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  printf '%s\n' '[[case]]' 'id = "tighten-case"' 'version = 1' 'purpose = "p"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'gauntlet_tier = "G0"' \
    'seed = 7' 'max_wall_seconds = 3600' 'entry_command = ["sleep", "5"]' \
    'expected = "authority"'
} > "$M"
bash "$RUNNER" --manifest "$M" --output-dir "$SCRATCH/tighten-out" --max-wall-seconds 1 >/dev/null 2>&1
check "--max-wall-seconds tightens the manifest budget" 14 $?

# --- cancellation and drain ----------------------------------------------
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  printf '%s\n' '[[case]]' 'id = "cancel-clean"' 'version = 1' 'purpose = "p"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'gauntlet_tier = "G0"' \
    'seed = 7' 'max_wall_seconds = 60' 'entry_command = ["sleep", "30"]' \
    'expected = "authority"'
} > "$M"
run_one "$M" --cancel-after 1
check "clean cancellation drains and maps to class 15" 15 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  printf '%s\n' '[[case]]' 'id = "cancel-stubborn"' 'version = 1' 'purpose = "p"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'gauntlet_tier = "G0"' \
    'seed = 7' 'max_wall_seconds = 60' \
    'entry_command = ["bash", "-c", "trap \"\" TERM; sleep 60"]' \
    'expected = "authority"'
} > "$M"
run_one "$M" --cancel-after 1
check "TERM-ignoring child is KILLed and maps to class 15 (drain failure)" 15 $?

# --- output budget --------------------------------------------------------
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  printf '%s\n' '[[case]]' 'id = "output-hog"' 'version = 1' 'purpose = "p"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'gauntlet_tier = "G0"' \
    'seed = 7' 'max_wall_seconds = 60' 'max_output_bytes = 1024' \
    'entry_command = ["bash", "-c", "yes overflow | head -c 100000"]' \
    'expected = "authority"'
} > "$M"
run_one "$M"
check "output over the declared cap maps to class 14" 14 $?

# --- determinism repeat ---------------------------------------------------
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case det-stable 'entry_command = ["echo", "stable"]' 'expected = "authority"' \
    'determinism_class = "repo-deterministic"'
} > "$M"
run_one "$M" --determinism-repeat
check "byte-identical repeat stays exit 0" 0 $?

{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case det-drifty 'entry_command = ["bash", "-c", "echo $RANDOM$$"]' \
    'expected = "authority"' 'determinism_class = "repo-deterministic"'
} > "$M"
run_one "$M" --determinism-repeat
check "drifting stdout under --determinism-repeat maps to class 16" 16 $?

# --- redaction by construction --------------------------------------------
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case no-secret-env \
    'entry_command = ["bash", "-c", "[ -z \"${SELFTEST_FAKE_SECRET_TOKEN:-}\" ]"]' \
    'expected = "authority"'
} > "$M"
SELFTEST_FAKE_SECRET_TOKEN="must-never-reach-children" run_one "$M"
check "credential-shaped env vars are scrubbed from children" 0 $?

# --- summary, log agreement, replay, and tamper ---------------------------
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case replay-case 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
OUT="$SCRATCH/replay-out"
bash "$RUNNER" --manifest "$M" --output-dir "$OUT" >/dev/null 2>&1
check "replay fixture run passes" 0 $?
bash "$RUNNER" --replay "$OUT/summary.json" >/dev/null 2>&1
check "replay agrees with an untouched run" 0 $?
printf '%s\n' '{"schema":"x","seq":999,"event":"case-terminal","case":"replay-case","status":"passed"}' >> "$OUT/runner-log.jsonl"
bash "$RUNNER" --replay "$OUT/summary.json" >/dev/null 2>&1
check "appended (tampered) log fails replay as class 17" 17 $?

# --- property/metamorphic battery ----------------------------------------
# Manifest key order must not change the projection or the verdict.
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  printf '%s\n' '[[case]]' 'id = "key-order-case"' 'version = 1' 'purpose = "p"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'gauntlet_tier = "G0"' \
    'seed = 7' 'max_wall_seconds = 60' 'entry_command = ["true"]' \
    'expected = "authority"'
} > "$M"
LIST_A="$(bash "$RUNNER" --manifest "$M" --list 2>/dev/null)"
{ printf '%s\n' 'phase = "e1"' "$CASE_HEADER"
  printf '%s\n' '[[case]]' 'expected = "authority"' 'entry_command = ["true"]' \
    'max_wall_seconds = 60' 'seed = 7' 'gauntlet_tier = "G0"' \
    'owning_bead = "frankensim-ext-epic-gov-rjoq.8"' 'purpose = "p"' \
    'version = 1' 'id = "key-order-case"'
} > "$M"
LIST_B="$(bash "$RUNNER" --manifest "$M" --list 2>/dev/null)"
[ -n "$LIST_A" ] && [ "$LIST_A" = "$LIST_B" ]
check "manifest key order does not change the projection" 0 $?

# Case-selection permutations must produce the same summary counts.
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case perm-a 'entry_command = ["true"]' 'expected = "authority"'
  good_case perm-b 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
bash "$RUNNER" --manifest "$M" --output-dir "$SCRATCH/perm-ab" --case perm-a --case perm-b >/dev/null 2>&1
bash "$RUNNER" --manifest "$M" --output-dir "$SCRATCH/perm-ba" --case perm-b --case perm-a >/dev/null 2>&1
python3 - "$SCRATCH/perm-ab/summary.json" "$SCRATCH/perm-ba/summary.json" <<'PYCMP'
import json, sys
ab, ba = (json.load(open(path)) for path in sys.argv[1:3])
assert ab["counts"] == ba["counts"] == {"passed": 2}, (ab["counts"], ba["counts"])
assert sorted(row["case"] for row in ab["cases"]) == sorted(row["case"] for row in ba["cases"])
PYCMP
check "case-selection permutation is semantically invariant" 0 $?

# The declared seed must actually reach the child.
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case seed-visible \
    'entry_command = ["bash", "-c", "[ \"$NEW_DOMAINS_SEED\" = 7 ]"]' \
    'expected = "authority"'
} > "$M"
run_one "$M"
check "the manifest seed reaches the child environment" 0 $?

# A truncated log must fail replay as tamper.
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case trunc-case 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
OUT="$SCRATCH/trunc-out"
bash "$RUNNER" --manifest "$M" --output-dir "$OUT" >/dev/null 2>&1
sed -i '' -e '$d' "$OUT/runner-log.jsonl" 2>/dev/null || sed -i -e '$d' "$OUT/runner-log.jsonl"
bash "$RUNNER" --replay "$OUT/summary.json" >/dev/null 2>&1
check "a truncated log fails replay as class 17" 17 $?

# An in-place event edit (no append) must also fail replay.
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case edit-case 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
OUT="$SCRATCH/edit-out"
bash "$RUNNER" --manifest "$M" --output-dir "$OUT" >/dev/null 2>&1
sed -i '' -e 's/"status": "passed"/"status": "failed"/' "$OUT/runner-log.jsonl" 2>/dev/null \
  || sed -i -e 's/"status": "passed"/"status": "failed"/' "$OUT/runner-log.jsonl"
bash "$RUNNER" --replay "$OUT/summary.json" >/dev/null 2>&1
check "an in-place edited log fails replay as class 17" 17 $?

# Equivalent relative and absolute manifest paths admit identically.
{ printf '%s\n' "$CASE_HEADER" 'phase = "e1"'
  good_case pathform-case 'entry_command = ["true"]' 'expected = "authority"'
} > "$M"
REL_OUT="$(cd "$SCRATCH" && bash "$RUNNER" --manifest ./e1.toml --validate-only 2>/dev/null)"
ABS_OUT="$(bash "$RUNNER" --manifest "$SCRATCH/e1.toml" --validate-only 2>/dev/null)"
[ -n "$REL_OUT" ] && [ "$REL_OUT" = "$ABS_OUT" ]
check "relative and absolute manifest paths validate identically" 0 $?

echo "selftest: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ] || exit 1
