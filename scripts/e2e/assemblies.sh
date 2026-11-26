#!/usr/bin/env bash
# No-mock E2E runner for the tracked assembly gallery (examples/assemblies/).
# Numerical, conservation, determinism, and provenance validation only (Estimate-only).
set -euo pipefail

CALLER_PWD="$PWD"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

ARTIFACT_DIR=""
PREBUILT_BINARY=""

usage() {
  printf '%s\n' \
    'usage: scripts/e2e/assemblies.sh --artifact-dir PATH [--prebuilt-binary PATH]' \
    'Runs examples/assemblies/*.fsim through real frankensim, checks independent' \
    'closed-form oracles, replay-audit determinism, resume, and refusals, and' \
    'writes stages.jsonl + summary-receipt.json into a new --artifact-dir.'
}

die() {
  printf 'assemblies_e2e result=refused elapsed_s=%s reason=%s\n' "$SECONDS" "$*" >&2
  exit 2
}

sha256_file() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 -- "$1" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum -- "$1" | awk '{print $1}'
  else
    die 'shasum or sha256sum is required for SHA-256 provenance'
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --artifact-dir) [[ $# -ge 2 && -n "${2:-}" ]] || die '--artifact-dir requires a path'; ARTIFACT_DIR="$2"; shift 2 ;;
    --prebuilt-binary) [[ $# -ge 2 && -n "${2:-}" ]] || die '--prebuilt-binary requires a path'; PREBUILT_BINARY="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -n "$ARTIFACT_DIR" ]] || die '--artifact-dir is required'
[[ "$ARTIFACT_DIR" = /* ]] || ARTIFACT_DIR="$CALLER_PWD/$ARTIFACT_DIR"
[[ -z "$PREBUILT_BINARY" || "$PREBUILT_BINARY" = /* ]] || PREBUILT_BINARY="$CALLER_PWD/$PREBUILT_BINARY"
[[ ! -e "$ARTIFACT_DIR" && ! -L "$ARTIFACT_DIR" ]] || die "refusing existing artifact dir: $ARTIFACT_DIR"
command -v python3 >/dev/null 2>&1 || die 'python3 is required for oracle verification'

ASM_DIR="examples/assemblies"
STEEL_PACK="$ASM_DIR/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
ALUM_PACK="$ASM_DIR/packs/aluminum-pure-thermal-nbs-c447.fsmcdpk"
ATKINS_PACK="examples/system/atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk"
SENS_PACK="$ASM_DIR/packs/aluminum-pure-nbs-c447-sensible-heat.fsmcdpk"
PG50_PACK="$ASM_DIR/packs/dowfrost-pg50-40-60c.fsmcdpk"
FILM_INPUTS=("$ASM_DIR/shear_film_pair.fsim" "$SENS_PACK" "$PG50_PACK")
for req in "$ASM_DIR/forced_oscillator.fsim" "$ASM_DIR/two_body_conduction.fsim" \
           "$ASM_DIR/thermo_mech_slider.fsim" "$STEEL_PACK" "$ALUM_PACK" "$ATKINS_PACK" "${FILM_INPUTS[@]}"; do
  [[ -f "$req" ]] || die "missing tracked input: $req"
done

if [[ -n "$PREBUILT_BINARY" ]]; then
  [[ -f "$PREBUILT_BINARY" && -x "$PREBUILT_BINARY" ]] || die "not executable: $PREBUILT_BINARY"
  BIN="$PREBUILT_BINARY"
  SOURCE_IDENTITY="prebuilt-binary-sha256:$(sha256_file "$BIN")"
else
  command -v git >/dev/null 2>&1 || die 'git is required unless --prebuilt-binary is supplied'
  git ls-files --error-unmatch Cargo.toml Cargo.lock "$ASM_DIR/forced_oscillator.fsim" \
    "$ASM_DIR/two_body_conduction.fsim" "$ASM_DIR/thermo_mech_slider.fsim" \
    "$STEEL_PACK" "$ALUM_PACK" "$ATKINS_PACK" "${FILM_INPUTS[@]}" >/dev/null 2>&1 \
    || die 'required sources or gallery files are untracked; commit or pass --prebuilt-binary'
  if [[ -n "$(git status --porcelain --untracked-files=all -- Cargo.toml Cargo.lock crates "$ASM_DIR" "$ATKINS_PACK")" ]]; then
    die 'dirty or untracked inputs in Cargo.toml, Cargo.lock, crates/, or gallery; pass --prebuilt-binary'
  fi
  SOURCE_IDENTITY="git-head:$(git rev-parse --verify HEAD);inputs=clean"
  "${CARGO:-cargo}" build --locked -p fs-cli --bin frankensim >/dev/null
  BIN="${CARGO_TARGET_DIR:-$REPO_ROOT/target}/debug/frankensim"
  [[ -x "$BIN" ]] || die "built binary not found at $BIN"
fi

RUNNER_SHA256="$(sha256_file "${BASH_SOURCE[0]}")"
mkdir -p "$ARTIFACT_DIR/ledgers" "$ARTIFACT_DIR/explains"
JSONL="$ARTIFACT_DIR/stages.jsonl"
RECEIPT="$ARTIFACT_DIR/summary-receipt.json"
: > "$JSONL"

run_case() {
  local tag="$1" fsim="$2"; shift 2
  local db="$ARTIFACT_DIR/ledgers/${tag}.db"
  local out="$ARTIFACT_DIR/explains/${tag}.run.json"
  local err="$ARTIFACT_DIR/explains/${tag}.stderr.jsonl"
  "$BIN" --json run "$fsim" "$db" "$@" > "$out" 2> "$err"
  local run_id
  run_id="$(python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); assert d["status"]=="completed" and d["stages_completed"]==6 and d["checker"]=="pass", d; print(d["run"])' "$out")"
  "$BIN" --json explain "$run_id" "$db" > "$ARTIFACT_DIR/explains/${tag}.explain.json"
  python3 -c 'import json,sys
tag, errp, outp = sys.argv[1:4]
with open(outp, "a") as w:
  for ln in open(errp):
    ln = ln.strip()
    if ln.startswith("{"):
      rec = json.loads(ln); rec["case"] = tag; w.write(json.dumps(rec, sort_keys=True) + "\n")' "$tag" "$err" "$JSONL"
}

run_case fo_transient "$ASM_DIR/forced_oscillator.fsim" --materials "$STEEL_PACK"
run_case fo_modal "$ASM_DIR/forced_oscillator.fsim" --materials "$STEEL_PACK" --analysis modal
run_case fo_harmonic "$ASM_DIR/forced_oscillator.fsim" --materials "$STEEL_PACK" --analysis harmonic
run_case fo_stationary "$ASM_DIR/forced_oscillator.fsim" --materials "$STEEL_PACK" --analysis stationary
run_case two_body "$ASM_DIR/two_body_conduction.fsim" --materials "$ALUM_PACK" --interfaces "$ATKINS_PACK"
run_case slider "$ASM_DIR/thermo_mech_slider.fsim" --materials "$STEEL_PACK"
run_case film "$ASM_DIR/shear_film_pair.fsim" --materials "$SENS_PACK" --materials "$PG50_PACK"

# Determinism repeat + replay-audit on slider across two fresh ledgers
sleep 0.05
run_case slider_repeat "$ASM_DIR/thermo_mech_slider.fsim" --materials "$STEEL_PACK"
"$BIN" --json replay-audit "$ARTIFACT_DIR/ledgers/slider.db" "$ARTIFACT_DIR/ledgers/slider_repeat.db" \
  > "$ARTIFACT_DIR/explains/replay_audit.json"

# Resume idempotency drill on completed slider ledger
SLIDER_RUN="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["run"])' "$ARTIFACT_DIR/explains/slider.run.json")"
"$BIN" --json run --resume "$SLIDER_RUN" "$ARTIFACT_DIR/ledgers/slider.db" > "$ARTIFACT_DIR/explains/slider_resume.json"

# Refusal drills (missing material card, invalid analysis mode, missing :area on contact resistance)
expect_refusal() {
  local tag="$1" want_rc="$2" needle="$3"; shift 3
  local rc=0 err="$ARTIFACT_DIR/explains/refusal_${tag}.stderr"
  "$BIN" "$@" >/dev/null 2>"$err" || rc=$?
  [[ "$rc" -eq "$want_rc" ]] || die "refusal $tag expected rc=$want_rc got rc=$rc"
  grep -q "$needle" "$err" || die "refusal $tag stderr missing '$needle'"
}
expect_refusal missing_card 4 "cli-assembly-card-missing" --json run "$ASM_DIR/thermo_mech_slider.fsim" "$ARTIFACT_DIR/ledgers/ref1.db"
expect_refusal bad_analysis 2 "cli-run-analysis-usage" --json run "$ASM_DIR/thermo_mech_slider.fsim" "$ARTIFACT_DIR/ledgers/ref2.db" --analysis bogus
NO_AREA="$ARTIFACT_DIR/two_body_no_area.fsim"
sed 's/:area 0.0025m\^2 //' "$ASM_DIR/two_body_conduction.fsim" > "$NO_AREA"
expect_refusal missing_area 4 "requires positive :area" --json run "$NO_AREA" "$ARTIFACT_DIR/ledgers/ref3.db" --materials "$ALUM_PACK" --interfaces "$ATKINS_PACK"

# Independent closed-form Python oracles + summary receipt emission
python3 - "$ARTIFACT_DIR" "$SOURCE_IDENTITY" "$RUNNER_SHA256" "$JSONL" "$RECEIPT" \
  "$ASM_DIR/forced_oscillator.fsim" "$ASM_DIR/two_body_conduction.fsim" "$ASM_DIR/thermo_mech_slider.fsim" \
  "$STEEL_PACK" "$ALUM_PACK" "$ATKINS_PACK" "${FILM_INPUTS[@]}" <<'PY'
import hashlib, json, math, pathlib, sys
art, src_id, runner_sha, jsonl_p, rcpt_p = sys.argv[1:6]
inputs = sys.argv[6:]
exp = pathlib.Path(art) / "explains"
ld = lambda name: json.loads((exp / name).read_text())
probes = lambda e: {p["name"]: p["values"] for p in e["probe_series"]["probes"]}
wres = lambda e: abs(e["window_ledger"]["totals"]["energy_residual_j"])

# 1. Forced oscillator (m=1.95 kg, c=4 N*s/m, k=200 N/m, F=10 N)
m, c, k, F = 7800.0 * 0.1 * 0.05 * 0.05, 4.0, 200.0, 10.0
wn = math.sqrt(k / m); zeta = c / (2.0 * math.sqrt(k * m)); wd = wn * math.sqrt(1.0 - zeta * zeta)
e_ft, e_fm, e_fh, e_fs = ld("fo_transient.explain.json"), ld("fo_modal.explain.json"), ld("fo_harmonic.explain.json"), ld("fo_stationary.explain.json")
pm, ph, ps, pt = probes(e_fm), probes(e_fh), probes(e_fs), probes(e_ft)
err_wn = abs(pm["modal.natural_rad_s"][0] - wn)
err_fd = abs(pm["modal.frequency_hz"][0] - wd / (2.0 * math.pi))
err_z = abs(pm["modal.damping_ratio"][0] - zeta)
assert max(err_wn, err_fd, err_z) <= 1e-9, (err_wn, err_fd, err_z)
err_harm = max(abs(got - abs(complex(c, m * 2*math.pi*f) / complex(k - m*(2*math.pi*f)**2, c*2*math.pi*f)))
               for f, got in zip((1.0, 2.0, 5.0), ph["harmonic.transfer_magnitude"]))
assert err_harm <= 1e-9, err_harm
err_xeq = abs(abs(ps["stationary.state.joint.spring_damper.spring.state"][0]) - F / k)
assert err_xeq <= 1e-9 and ps["stationary.residual_norm"][0] <= 1e-9, (err_xeq, ps)
xs, vs = pt["mass_pos"], pt["mass_vel"]
fo_ex = max(abs(x - (F/k)*(1.0 - math.exp(-zeta*wn*i*0.001)*(math.cos(wd*i*0.001) + (zeta*wn/wd)*math.sin(wd*i*0.001)))) for i, x in enumerate(xs))
fo_ev = max(abs(v - (F/(m*wd))*math.exp(-zeta*wn*i*0.001)*math.sin(wd*i*0.001)) for i, v in enumerate(vs))
assert fo_ex <= 1e-5 and fo_ev <= 1e-4 and wres(e_ft) <= 1e-12, (fo_ex, fo_ev, wres(e_ft))

# 2. Two-body conduction (NBS C447 Al + Atkins 1964 R''=3.1606e-4 at 313 kPa, A=0.0025 m^2)
e_tb = ld("two_body.explain.json"); ptb = probes(e_tb)
assert "area-specific-thermal-contact-resistance" in json.dumps(e_tb["lowering_receipt"])
Ch = 2700.0 * (0.06 * 0.05 * 0.05) * 945.584; Cc = 2700.0 * (0.04 * 0.05 * 0.05) * 945.584
G = 0.0025 / 3.1606e-4; Q = 500.0; T0 = 293.15; tau = Ch * Cc / (G * (Ch + Cc))
Th, Tc = ptb["hot_temp"], ptb["cold_temp"]
tb_eh = tb_ec = tb_en = 0.0
for i, (th, tc) in enumerate(zip(Th, Tc)):
  t = i * 0.05; t_avg = T0 + (Q / (Ch + Cc)) * t; dlt = (Q * Cc / (G * (Ch + Cc))) * (1.0 - math.exp(-t / tau))
  tb_eh = max(tb_eh, abs(th - (t_avg + (Cc / (Ch + Cc)) * dlt)))
  tb_ec = max(tb_ec, abs(tc - (t_avg - (Ch / (Ch + Cc)) * dlt)))
  tb_en = max(tb_en, abs(Ch * (th - T0) + Cc * (tc - T0) - Q * t))
assert tb_eh <= 1e-4 and tb_ec <= 1e-4 and tb_en <= 1e-8 and wres(e_tb) <= 1e-9, (tb_eh, tb_ec, tb_en)
assert 289.35 <= Th[-1] <= 305.85 and 289.35 <= Tc[-1] <= 305.85, (Th[-1], Tc[-1])

# 3. Thermo-mech slider (m=1.95 kg, C=947.7 J/K, c=12, k=400, F=5+15 sin(6pi t + 0.25))
e_sl = ld("slider.explain.json"); psl = probes(e_sl)
C_s = m * 486.0; cs, ks, F0, Af, om, phi, dt = 12.0, 400.0, 5.0, 15.0, 2*math.pi*3.0, 0.25, 0.001
wns = math.sqrt(ks / m); zs = cs / (2.0 * math.sqrt(ks * m)); wds = wns * math.sqrt(1.0 - zs * zs)
amp_c = Af * complex(math.cos(phi), math.sin(phi)) / complex(ks - m*om*om, cs*om)
xp = lambda t: F0/ks + (amp_c * complex(math.cos(om*t), math.sin(om*t))).imag
vp = lambda t: (1j * om * amp_c * complex(math.cos(om*t), math.sin(om*t))).imag
C1 = -xp(0.0); C2 = (-vp(0.0) + zs * wns * C1) / wds
sxs, svs, traw, tlag = psl["slider_pos"], psl["slider_vel"], psl["slider_temp_raw"], psl["slider_temp"]
sl_ex = sl_ev = 0.0
for i, (x, v) in enumerate(zip(sxs, svs)):
  t = i * dt; et = math.exp(-zs * wns * t)
  xh = et * (C1*math.cos(wds*t) + C2*math.sin(wds*t))
  vh = et * ((-zs*wns*C1 + wds*C2)*math.cos(wds*t) + (-zs*wns*C2 - wds*C1)*math.sin(wds*t))
  sl_ex = max(sl_ex, abs(x - (xp(t) + xh))); sl_ev = max(sl_ev, abs(v - (vp(t) + vh)))
sinc = math.sin(0.5 * om * dt) / (0.5 * om * dt)
W_ext = Q_damp = sl_etot = sl_eth = 0.0
for i in range(len(sxs) - 1):
  vm = 0.5 * (svs[i] + svs[i + 1]); fm = F0 + Af * math.sin(om * (i + 0.5) * dt + phi) * sinc
  W_ext += dt * fm * vm; Q_damp += dt * cs * vm * vm
  Em = 0.5 * m * svs[i + 1]**2 + 0.5 * ks * sxs[i + 1]**2; Et = C_s * (traw[i + 1] - traw[0])
  sl_etot = max(sl_etot, abs(Em + Et - W_ext)); sl_eth = max(sl_eth, abs(Et - Q_damp))
alpha = 1.0 - math.exp(-dt / 0.005); tf = [traw[0]]
for y in traw[1:]: tf.append(tf[-1] + alpha * (y - tf[-1]))
sl_tc = max(abs(a - b) for a, b in zip(tlag, tf))
assert sl_ex <= 1e-5 and sl_ev <= 1e-4 and sl_etot <= 1e-9 and sl_eth <= 1e-9 and sl_tc <= 1e-12 and wres(e_sl) <= 1e-12

# 4. Eight-term QoI uncertainty budgets & verdicts on gallery transient runs
CANON_TERMS = ["roundoff", "solver-algebraic", "discretization", "geometry", "parameters", "boundary-conditions", "model-form", "measurement"]
def check_budget(e, req_id, min_disc=0.0):
  qb = e["qoi_budget"]
  assert qb and len(qb["requirements"]) == 1 and qb["requirements"][0]["id"] == req_id, qb
  r = qb["requirements"][0]
  assert r["outcome"] == "indeterminate" and [t["kind"] for t in r["terms"]] == CANON_TERMS, r
  byk = {t["kind"]: t for t in r["terms"]}
  assert byk["roundoff"]["state"] == "measured" and byk["roundoff"]["value"] >= 0.0, r
  assert byk["discretization"]["state"] == "measured" and byk["discretization"]["value"] >= min_disc, r
  assert all(k in r["unmeasured_terms"] and byk[k]["state"] == "no-data" for k in ("geometry", "parameters", "boundary-conditions", "model-form", "measurement")), r
  assert any(x["label"] == "step_halving" and x["wall_s"] >= 0.0 for x in qb["extra_runs"]), qb["extra_runs"]
  return {"id": req_id, "outcome": r["outcome"], "solver_state": byk["solver-algebraic"]["state"], "unmeasured_terms": r["unmeasured_terms"]}
bv_ft, bv_tb, bv_sl = check_budget(e_ft, "pos_limit", fo_ex), check_budget(e_tb, "hot_limit", tb_eh), check_budget(e_sl, "temp_limit")

# 5. Determinism replay-audit + resume idempotency
ra = ld("replay_audit.json"); sr = ld("slider.run.json"); sres = ld("slider_resume.json")
assert ra["status"] == "clean" and ra["clean"] is True and ra["compared"] == 7 and not ra["deterministic_mismatches"], ra
assert sres["run"] == sr["run"] and sres["status"] == "completed" and sres["stages_completed"] == 6 and e_sl["plan_hash"] == sr["plan_hash"], sres

# 6. Shear-film pair (Al discs r=0.05 m, h=0.01 m, +-10 rad/s about the shared axis; PG50 mu=3.10 mPa*s, gap 0.1 mm)
e_sf = ld("film.explain.json"); wl = e_sf["window_ledger"]; tot = wl["totals"]; Qf = tot["dissipated_j"]
nodes = {n["name"]: n for n in wl["nodes"]}
Iz = 2700.0 * math.pi * 0.05**2 * 0.01 * 0.05**2 / 2.0; ct = 0.0031 * (math.pi * 0.05**2)**2 / (2.0 * math.pi) / 1e-4
sf_q = abs(Qf / (0.5 * (Iz / 2.0) * 20.0**2 * (1.0 - math.exp(-4.0 * ct * wl["duration_s"] / Iz))) - 1.0)
sf_acc = max(abs(nodes["thermal-timeline"]["stored_change_j"] - 0.9 * Qf), abs(tot["sink_heat_j"] - 0.1 * Qf),
             abs(nodes["viscous-network"]["stored_change_j"] + Qf), wres(e_sf), abs(tot["numerical_loss_j"])) / Qf
assert wl["duration_s"] == 1.0 and tot["boundary_energy_in_j"] == 0.0 and sf_q <= 1e-12 and sf_acc <= 1e-12, (sf_q, sf_acc, tot)
# Each disc's spin decays at 2 c_t / I_z (the energy at twice that); each cell's C (T - T0) is its share of the heat,
# within two roundings at T's size, 2 ulp(T) C, from each sample's own T (MB 39197, 39225).
psf = probes(e_sf); Cc = 2700.0 * 0.02**3 * 945.584; Qt = lambda t: 0.5 * (Iz / 2.0) * 20.0**2 * (1.0 - math.exp(-4.0 * ct * t / Iz))
sf_w = max(abs(w - s * 10.0 * math.exp(-2.0 * ct * i * 0.01 / Iz)) / 10.0
           for s, k in ((1.0, "lower_spin"), (-1.0, "upper_spin")) for i, w in enumerate(psf[k]))
sf_c = max(abs(Cc * (T - 293.15) - f * Qt(i * 0.01)) / (2.0 * math.ulp(T) * Cc)
           for f, k in ((0.7, "lower_cell_temp"), (0.2, "upper_cell_temp")) for i, T in enumerate(psf[k]))
assert all(len(v) == 101 for v in psf.values()) and len(psf) == 4 and sf_w <= 1e-12 and sf_c <= 1.0, (sf_w, sf_c)

receipt = {
  "schema": "frankensim.e2e.assemblies-summary.v1",
  "authority": "Estimate-only",
  "no_claim": "Numerical, conservation, determinism, and provenance check against closed-form analytical oracles only; no experimental or physical validation claim.",
  "source_identity": src_id, "runner_sha256": runner_sha,
  "inputs_sha256": {p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest() for p in inputs},
  "assemblies": {
    "forced_oscillator": {"runs": {m: ld(f"fo_{m}.run.json")["run"] for m in ("transient","modal","harmonic","stationary")}, "verdict": bv_ft,
      "errors": {"modal_wn_rad_s": err_wn, "modal_fd_hz": err_fd, "modal_zeta": err_z, "harmonic_max": err_harm, "stationary_x_m": err_xeq, "transient_x_m": fo_ex, "transient_v_m_s": fo_ev, "energy_residual_j": wres(e_ft)}},
    "two_body_conduction": {"run": ld("two_body.run.json")["run"], "verdict": bv_tb,
      "errors": {"hot_temp_k": tb_eh, "cold_temp_k": tb_ec, "enthalpy_j": tb_en, "energy_residual_j": wres(e_tb)}},
    "thermo_mech_slider": {"run": sr["run"], "verdict": bv_sl,
      "errors": {"pos_m": sl_ex, "vel_m_s": sl_ev, "midpoint_total_energy_j": sl_etot, "midpoint_thermal_energy_j": sl_eth, "thermocouple_lag_k": sl_tc, "energy_residual_j": wres(e_sl)}},
    "shear_film_pair": {"run": ld("film.run.json")["run"],
      "errors": {"free_pair_dissipated_rel": sf_q, "ledger_accounts_rel": sf_acc, "energy_residual_j": wres(e_sf),
                 "disc_spin_rel": sf_w, "cell_heat_over_2ulp_floor": sf_c}},
  },
  "drills": {"replay_audit_clean": True, "resume_idempotent": True, "refusals": ["missing_card", "bad_analysis", "missing_area"]},
}
pathlib.Path(rcpt_p).write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
with open(jsonl_p, "a") as w:
  w.write(json.dumps({"stage": "oracles_and_drills", "status": "pass", "receipt": rcpt_p}, sort_keys=True) + "\n")
print(f"assemblies_e2e result=pass receipt={rcpt_p}")
PY
