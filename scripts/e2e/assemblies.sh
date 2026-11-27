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
SRM_PACK="$ASM_DIR/packs/eps-nist-srm1453-2013.fsmcdpk"
N2_PACK="$ASM_DIR/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk"
KATS_PACK="$ASM_DIR/packs/katsivardi-2025-pla-metal-dry-air-interface.fsintpk"
FILM_INPUTS=("$ASM_DIR/shear_film_pair.fsim" "$SENS_PACK" "$PG50_PACK" "$ASM_DIR/srm1453_slab.fsim" "$SRM_PACK"
  "$ASM_DIR/two_mass_chain.fsim" "$ASM_DIR/gas_piston.fsim" "$N2_PACK"
  "$ASM_DIR/block_on_belt.fsim" "$KATS_PACK")
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
run_case srm1453 "$ASM_DIR/srm1453_slab.fsim" --materials "$SRM_PACK" --analysis stationary
# The device added as a file only (no Rust change), and its transient again at half the step
run_case tmc_modal "$ASM_DIR/two_mass_chain.fsim" --materials "$STEEL_PACK" --analysis modal
run_case tmc_stationary "$ASM_DIR/two_mass_chain.fsim" --materials "$STEEL_PACK" --analysis stationary
run_case tmc_transient "$ASM_DIR/two_mass_chain.fsim" --materials "$STEEL_PACK"
TMC_HALF="$ARTIFACT_DIR/two_mass_chain_half_dt.fsim"
sed 's/:dt 0.001s/:dt 0.0005s/' "$ASM_DIR/two_mass_chain.fsim" > "$TMC_HALF"
grep -q ":dt 0.0005s" "$TMC_HALF" || die "half-step copy of two_mass_chain.fsim has no :dt 0.0005s"
run_case tmc_half "$TMC_HALF" --materials "$STEEL_PACK"
# Gas piston (adiabatic Ruchardt modal, equilibrium transient, and +500 N compression transient)
run_case gp_modal "$ASM_DIR/gas_piston.fsim" --materials "$STEEL_PACK" --materials "$N2_PACK" --analysis modal
run_case gp_transient "$ASM_DIR/gas_piston.fsim" --materials "$STEEL_PACK" --materials "$N2_PACK"
GP_COMP="$ARTIFACT_DIR/gas_piston_comp.fsim"
sed 's/-6361.725123519331N/-6861.725123519331N/' "$ASM_DIR/gas_piston.fsim" > "$GP_COMP"
grep -q -- "-6861.725123519331N" "$GP_COMP" || die "compression copy of gas_piston.fsim missing shifted load"
run_case gp_comp "$GP_COMP" --materials "$STEEL_PACK" --materials "$N2_PACK"
# Dry-friction stick-slip block on belt (Katsivardi et al. 2025 pure PLA / dry metallic plate)
run_case bob "$ASM_DIR/block_on_belt.fsim" --materials "$STEEL_PACK" --interfaces "$KATS_PACK"

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
expect_refusal srm_transient 4 "cli-assembly-fidelity-analysis" --json run "$ASM_DIR/srm1453_slab.fsim" "$ARTIFACT_DIR/ledgers/ref4.db" --materials "$SRM_PACK" --analysis transient
SRM_BAD_DT="$ARTIFACT_DIR/srm1453_bad_dt.fsim"
sed 's/307.15K/309.15K/' "$ASM_DIR/srm1453_slab.fsim" > "$SRM_BAD_DT"
expect_refusal srm_bad_dt 4 "face temperature difference must be 20 K" --json run "$SRM_BAD_DT" "$ARTIFACT_DIR/ledgers/ref5.db" --materials "$SRM_PACK" --analysis stationary
SRM_BAD_LX="$ARTIFACT_DIR/srm1453_bad_lx.fsim"
sed 's/0.01343m/0.0254m/' "$ASM_DIR/srm1453_slab.fsim" > "$SRM_BAD_LX"
expect_refusal srm_bad_lx 4 "outside SRM 1453 production-run range" --json run "$SRM_BAD_LX" "$ARTIFACT_DIR/ledgers/ref6.db" --materials "$SRM_PACK" --analysis stationary
SRM_FAIL="$ARTIFACT_DIR/srm1453_fail.fsim"
sed 's/3.35W/2.95W/' "$ASM_DIR/srm1453_slab.fsim" > "$SRM_FAIL"
run_case srm1453_fail "$SRM_FAIL" --materials "$SRM_PACK" --analysis stationary
SRM_IND="$ARTIFACT_DIR/srm1453_ind.fsim"
sed 's/3.35W/3.20W/' "$ASM_DIR/srm1453_slab.fsim" > "$SRM_IND"
run_case srm1453_ind "$SRM_IND" --materials "$SRM_PACK" --analysis stationary
BOB_OOD="$ARTIFACT_DIR/block_on_belt_ood_n.fsim"
sed 's/:normal-force 5N/:normal-force 6N/' "$ASM_DIR/block_on_belt.fsim" > "$BOB_OOD"
expect_refusal bob_ood_n 4 "outside the validity" --json run "$BOB_OOD" "$ARTIFACT_DIR/ledgers/ref7.db" --materials "$STEEL_PACK" --interfaces "$KATS_PACK"

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

# 4. Eight-term QoI uncertainty budgets & verdicts across all gallery runs
CANON_TERMS = ["roundoff", "solver-algebraic", "discretization", "geometry", "parameters", "boundary-conditions", "model-form", "measurement"]
def check_budget(e, req_id, min_disc=0.0, roundoff_measured=True):
  qb = e["qoi_budget"]
  assert qb and len(qb["requirements"]) == 1 and qb["requirements"][0]["id"] == req_id, qb
  r = qb["requirements"][0]
  assert r["outcome"] == "indeterminate" and [t["kind"] for t in r["terms"]] == CANON_TERMS, r
  byk = {t["kind"]: t for t in r["terms"]}
  if roundoff_measured:
    assert byk["roundoff"]["state"] == "measured" and byk["roundoff"]["value"] >= 0.0, r
  else:
    assert byk["roundoff"]["state"] == "no-data" and "roundoff" in r["unmeasured_terms"], r
  assert byk["solver-algebraic"]["state"] == "measured" and byk["solver-algebraic"]["value"] >= 0.0, r
  assert byk["discretization"]["state"] == "measured" and byk["discretization"]["value"] >= min_disc, r
  assert all(k in r["unmeasured_terms"] and byk[k]["state"] == "no-data" for k in ("geometry", "parameters", "boundary-conditions", "model-form", "measurement")), r
  assert any(x["label"] == "step_halving" and x["wall_s"] >= 0.0 for x in qb["extra_runs"]), qb["extra_runs"]
  return {"id": req_id, "outcome": r["outcome"], "solver_state": byk["solver-algebraic"]["state"], "unmeasured_terms": r["unmeasured_terms"]}
bv_ft, bv_tb, bv_sl = check_budget(e_ft, "pos_limit", fo_ex), check_budget(e_tb, "hot_limit", tb_eh), check_budget(e_sl, "temp_limit")
bv_fm = check_budget(e_fm, "mode1_freq", roundoff_measured=False)
bv_fh = check_budget(e_fh, "harm_peak", roundoff_measured=False)
bv_fs = check_budget(e_fs, "eq_extension", roundoff_measured=False)

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
bv_sf = check_budget(e_sf, "cell_temp_limit")

# 7. NIST SRM 1453 EPS slab (stationary 1-D steady conduction + 8-term determinate pass/fail/indeterminate)
e_srm = ld("srm1453.explain.json"); psrm = probes(e_srm)
k_srm = 0.00111 - 0.0000424 * 41.0 + 0.000115 * 297.15
q_srm = k_srm * (0.25 * 0.25 / 0.01343) * 20.0
srm_rel = abs(psrm["slab_heat_rate"][0] - q_srm) / q_srm
assert srm_rel <= 1e-12 and psrm["stationary.residual_norm"][0] == 0.0, (srm_rel, psrm)
r_pass = e_srm["qoi_budget"]["requirements"][0]
byk_srm = {t["kind"]: t for t in r_pass["terms"]}
assert r_pass["outcome"] == "pass" and r_pass["unmeasured_terms"] == [] and [t["kind"] for t in r_pass["terms"]] == CANON_TERMS, r_pass
assert byk_srm["discretization"]["state"] == "negligible" and byk_srm["model-form"]["state"] == "negligible", byk_srm
assert all(byk_srm[k]["state"] == "measured" and byk_srm[k]["value"] > 0.0 for k in ("roundoff", "parameters", "geometry", "boundary-conditions", "measurement")), byk_srm
r_fail = ld("srm1453_fail.explain.json")["qoi_budget"]["requirements"][0]
r_ind = ld("srm1453_ind.explain.json")["qoi_budget"]["requirements"][0]
assert r_fail["outcome"] == "fail" and r_fail["unmeasured_terms"] == [], r_fail
assert r_ind["outcome"] == "indeterminate" and r_ind["unmeasured_terms"] == [], r_ind

# 8. Two-mass chain, added as a file only (steel card rho = 7800 kg/m^3: m1 = 1.95 kg, m2 = 1.17 kg; anchor -k1,c1- m1 -k2,c2- m2,
# k1 = 200, k2 = 150 N/m, c = 0.02 k (stiffness-proportional, so the undamped modes decouple), F = 10 N on m2 from rest)
m1, m2 = 7800.0 * 0.1 * 0.05 * 0.05, 7800.0 * 0.06 * 0.05 * 0.05
k1, k2, beta, Ft = 200.0, 150.0, 0.02, 10.0
qa, qb, qc = m1 * m2, -(m1 * k2 + m2 * (k1 + k2)), k1 * k2
tm_modes = []
for w2 in sorted(((-qb - math.sqrt(qb*qb - 4*qa*qc)) / (2*qa), (-qb + math.sqrt(qb*qb - 4*qa*qc)) / (2*qa))):
  phi = (k2, k1 + k2 - w2 * m1); w = math.sqrt(w2); zt = beta * w / 2.0
  tm_modes.append((w2, phi, phi[1] * Ft / (m1 * phi[0]**2 + m2 * phi[1]**2), w, zt, w * math.sqrt(1.0 - zt * zt)))
def tm_exact(t):
  x = [0.0, 0.0]
  for w2, phi, f, w, zt, wd in tm_modes:
    q = f / w2 * (1.0 - math.exp(-zt*w*t) * (math.cos(wd*t) + zt*w/wd * math.sin(wd*t)))
    x[0] += phi[0] * q; x[1] += phi[1] * q
  return x
e_tm, e_ts, e_tt, e_th = ld("tmc_modal.explain.json"), ld("tmc_stationary.explain.json"), ld("tmc_transient.explain.json"), ld("tmc_half.explain.json")
ptm, pts = probes(e_tm), probes(e_ts)
# Each mode appears once per conjugate eigenvalue: (w1, w1, w2, w2)
tm_err_modal = max(abs(got - want) for key, pick in (("modal.natural_rad_s", lambda m: m[3]), ("modal.frequency_hz", lambda m: m[5] / (2.0 * math.pi)),
                                                    ("modal.damping_ratio", lambda m: m[4]))
                   for got, want in zip(ptm[key], [pick(m) for m in tm_modes for _ in (0, 1)]))
assert all(len(ptm[k]) == 4 for k in ("modal.natural_rad_s", "modal.frequency_hz", "modal.damping_ratio")) and tm_err_modal <= 1e-9, (tm_err_modal, ptm)
tm_err_stat = max(abs(pts["stationary.state.joint.ground_link.spring.state"][0] + Ft / k1),
                  abs(pts["stationary.state.joint.middle_link.spring.state"][0] + Ft / k2),
                  abs(pts["stationary.state.inner.velocity_m_s"][0]), abs(pts["stationary.state.outer.velocity_m_s"][0]))
assert tm_err_stat <= 1e-9 and pts["stationary.residual_norm"][0] <= 1e-9, (tm_err_stat, pts)
def tm_err(e, horizon=1.0):
  pr = probes(e); xi, xo = pr["inner_pos"], pr["outer_pos"]; dt = horizon / (len(xi) - 1)
  return (max(abs(x - tm_exact(i * dt)[0]) for i, x in enumerate(xi)), max(abs(x - tm_exact(i * dt)[1]) for i, x in enumerate(xo)), len(xi))
tm_full, tm_half = tm_err(e_tt), tm_err(e_th)
tm_ratio = (tm_full[0] / tm_half[0], tm_full[1] / tm_half[1])
# Second order: the error against the exact modal superposition falls by 4 when dt halves
assert tm_full[2] == 1001 and tm_half[2] == 2001 and max(tm_full[:2]) <= 5e-6 and all(3.9 <= r <= 4.1 for r in tm_ratio), (tm_full, tm_half, tm_ratio)
assert wres(e_tt) <= 1e-12 and wres(e_th) <= 1e-12, (wres(e_tt), wres(e_th))

# 9. Gas piston (N2 NASA-9 McBride 2002 + Lemmon 2004, r=0.15 m, h=0.20 m, T0=300 K, p0=90 kPa; steel 0.1 m cube m=7.8 kg)
gp_r, gp_h, gp_t0, gp_p0, gp_mp = 0.15, 0.20, 300.0, 90000.0, 7800.0 * 0.1**3
gp_A, gp_V0 = math.pi * gp_r**2, math.pi * gp_r**2 * gp_h
# Exact NASA-9 N2 (200..1000 K, McBride 2002) cp/R and molar mass 0.0280134 kg/mol:
gp_R = 8.31446261815324 / 0.0280134
gp_cp_R = lambda T: (22103.71497 * T**-2 - 381.846182 * T**-1 + 6.08273836
                     - 8.53091441e-3 * T + 1.384646189e-5 * T**2
                     - 9.62579362e-9 * T**3 + 2.519705809e-12 * T**4)
gp_cv0 = (gp_cp_R(gp_t0) - 1.0) * gp_R; gp_gam = gp_cp_R(gp_t0) / (gp_cp_R(gp_t0) - 1.0)
gp_wiso2 = gp_p0 * gp_A**2 / (gp_mp * gp_V0); gp_wad2 = gp_gam * gp_wiso2
gp_fad = math.sqrt(gp_wad2) / (2.0 * math.pi)
e_gpm, e_gpt, e_gpc = ld("gp_modal.explain.json"), ld("gp_transient.explain.json"), ld("gp_comp.explain.json")
pgpm, pgpt, pgpc = probes(e_gpm), probes(e_gpt), probes(e_gpc)
gp_ad_idx = next(i for i, f in enumerate(pgpm["modal.frequency_hz"]) if f > 1.0)
gp_fad_rel = abs(pgpm["modal.frequency_hz"][gp_ad_idx] - gp_fad) / gp_fad
gp_z_ad = abs(pgpm["modal.damping_ratio"][gp_ad_idx])
assert gp_fad_rel <= 1e-6 and gp_z_ad <= 1e-6, (gp_fad_rel, gp_z_ad, pgpm)
bv_gpm = check_budget(e_gpm, "ruchardt_freq", roundoff_measured=False)
bv_gpt = check_budget(e_gpt, "gas_temp_max", roundoff_measured=False)
bv_gpc = check_budget(e_gpc, "gas_temp_max", roundoff_measured=False)
gp_eq_v = max(abs(v) for v in pgpt["piston_vel"]); gp_eq_t = max(abs(t - gp_t0) for t in pgpt["gas_temp"])
assert gp_eq_v <= 1e-11 and gp_eq_t <= 1e-11 and wres(e_gpt) <= 1e-12, (gp_eq_v, gp_eq_t, wres(e_gpt))
# Independent RK4 adiabatic ideal-gas ODE for +500 N compression step (dt_rk4 = 1e-5 s)
gp_st = [0.0, 0.0, gp_V0, gp_t0]; gp_fext = -gp_p0 * gp_A - 500.0
gp_mg = gp_p0 * gp_V0 / (gp_R * gp_t0)
def gp_rhs(y):
  p = gp_mg * gp_R * y[3] / y[2]; cv = (gp_cp_R(y[3]) - 1.0) * gp_R
  return [y[1], (p * gp_A + gp_fext) / gp_mp, gp_A * y[1], -(p * gp_A * y[1]) / (gp_mg * cv)]
gp_rk = [gp_st[:]]
for _ in range(40):
  for _ in range(50):
    h4 = 1e-5; k1 = gp_rhs(gp_st)
    k2 = gp_rhs([a + 0.5*h4*b for a, b in zip(gp_st, k1)])
    k3 = gp_rhs([a + 0.5*h4*b for a, b in zip(gp_st, k2)])
    k4 = gp_rhs([a + h4*b for a, b in zip(gp_st, k3)])
    gp_st = [a + h4*(b1 + 2*b2 + 2*b3 + b4)/6.0 for a, b1, b2, b3, b4 in zip(gp_st, k1, k2, k3, k4)]
  gp_rk.append(gp_st[:])
gp_comp_x = max(abs(a - b[0]) for a, b in zip(pgpc["piston_pos"], gp_rk))
gp_comp_v = max(abs(a - b[1]) for a, b in zip(pgpc["piston_vel"], gp_rk))
gp_comp_t = max(abs(a - b[3]) for a, b in zip(pgpc["gas_temp"], gp_rk))
assert len(pgpc["gas_temp"]) == 41 and gp_comp_x <= 5e-6 and gp_comp_v <= 1e-4 and gp_comp_t <= 2e-3 and wres(e_gpc) <= 1e-9, (gp_comp_x, gp_comp_v, gp_comp_t, wres(e_gpc))

# 10. Dry-friction stick-slip block on belt (Katsivardi et al. 2025 Table 3: mu_s=0.371, mu_k=0.195 at N=5 N, v=0.0025 m/s; 500 g steel sled m=0.5 kg, k=1000 N/m, ambient heat)
e_bob = ld("bob.explain.json"); pbob = probes(e_bob)
bob_m, bob_k, bob_N, bob_v, bob_mus, bob_muk, bob_dt = 0.5, 1000.0, 5.0, 0.0025, 0.371, 0.195, 0.001
bob_w = math.sqrt(bob_k / bob_m); bob_xs, bob_xk = bob_mus * bob_N / bob_k, bob_muk * bob_N / bob_k; bob_d = bob_xs - bob_xk
bob_T_ex = 2.0 * bob_d / bob_v + 2.0 * (math.pi - math.atan(bob_w * bob_d / bob_v)) / bob_w
bob_W_ex = bob_muk * bob_N * bob_v * bob_T_ex; bob_vslip_ex = bob_v + math.hypot(bob_v, bob_d * bob_w)
bpos, bvel = pbob["block_pos"], pbob["block_vel"]
bob_cr = [(i - 1 + (bob_xk - bpos[i - 1]) / (bpos[i] - bpos[i - 1])) * bob_dt
          for i in range(1, len(bpos)) if abs(bvel[i] - bob_v) < 1e-6 and bpos[i - 1] < bob_xk <= bpos[i]]
bob_periods = [bob_cr[j + 1] - bob_cr[j] for j in range(len(bob_cr) - 1)]
bob_T_err = max(abs(t - bob_T_ex) for t in bob_periods)
bob_tot = e_bob["window_ledger"]["totals"]; bob_diss = bob_tot["dissipated_j"]
bob_W_cyc = bob_diss / len(bob_periods); bob_W_rel = abs(bob_W_cyc - bob_W_ex) / bob_W_ex
bob_amb_rel = abs(bob_tot["sink_heat_j"] - bob_diss) / bob_diss
bob_vslip_pk = max(abs(bob_v - v) for v in bvel)
bob_rcpt = json.dumps(e_bob["lowering_receipt"])
print(f"block_on_belt: cycles={len(bob_periods)} T_exact={bob_T_ex:.9f}s max_period_err={bob_T_err:.3e}s (<dt={bob_dt}s) W_cyc={bob_W_cyc*1e3:.6f}mJ W_rel={bob_W_rel:.3e} peak_slip={bob_vslip_pk*1e3:.3f}mm/s (exact {bob_vslip_ex*1e3:.3f}mm/s, {bob_vslip_pk/bob_v:.2f}xV)")
assert len(bob_periods) == 4 and bob_T_err < bob_dt and bob_W_rel <= bob_dt / bob_T_ex + 1e-9, (len(bob_periods), bob_T_err, bob_W_rel)
assert bob_amb_rel <= 1e-12 and wres(e_bob) <= 1e-12 and abs(bob_vslip_pk - bob_vslip_ex) / bob_vslip_ex <= 1e-2, (bob_amb_rel, wres(e_bob), bob_vslip_pk, bob_vslip_ex)
assert "static-friction-coefficient.held" in bob_rcpt and "kinetic-friction-coefficient.held" in bob_rcpt
assert "default: no thermal side; all heat to ambient" in bob_rcpt
bv_bob = check_budget(e_bob, "block_pos_cap", roundoff_measured=True)

receipt = {
  "schema": "frankensim.e2e.assemblies-summary.v1",
  "authority": "Estimate-only",
  "no_claim": "Numerical, conservation, determinism, and provenance check against closed-form analytical oracles only; no experimental or physical validation claim.",
  "source_identity": src_id, "runner_sha256": runner_sha,
  "inputs_sha256": {p: hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest() for p in inputs},
  "assemblies": {
    "forced_oscillator": {"runs": {m: ld(f"fo_{m}.run.json")["run"] for m in ("transient","modal","harmonic","stationary")},
      "verdict": bv_ft, "mode_verdicts": {"modal": bv_fm, "harmonic": bv_fh, "stationary": bv_fs},
      "errors": {"modal_wn_rad_s": err_wn, "modal_fd_hz": err_fd, "modal_zeta": err_z, "harmonic_max": err_harm, "stationary_x_m": err_xeq, "transient_x_m": fo_ex, "transient_v_m_s": fo_ev, "energy_residual_j": wres(e_ft)}},
    "two_body_conduction": {"run": ld("two_body.run.json")["run"], "verdict": bv_tb,
      "errors": {"hot_temp_k": tb_eh, "cold_temp_k": tb_ec, "enthalpy_j": tb_en, "energy_residual_j": wres(e_tb)}},
    "thermo_mech_slider": {"run": sr["run"], "verdict": bv_sl,
      "errors": {"pos_m": sl_ex, "vel_m_s": sl_ev, "midpoint_total_energy_j": sl_etot, "midpoint_thermal_energy_j": sl_eth, "thermocouple_lag_k": sl_tc, "energy_residual_j": wres(e_sl)}},
    "shear_film_pair": {"run": ld("film.run.json")["run"], "verdict": bv_sf,
      "errors": {"free_pair_dissipated_rel": sf_q, "ledger_accounts_rel": sf_acc, "energy_residual_j": wres(e_sf),
                 "disc_spin_rel": sf_w, "cell_heat_over_2ulp_floor": sf_c}},
    "srm1453_slab": {"run": ld("srm1453.run.json")["run"],
      "verdict": {"id": r_pass["id"], "outcome": r_pass["outcome"], "fail_drill": r_fail["outcome"], "indeterminate_drill": r_ind["outcome"], "unmeasured_terms": r_pass["unmeasured_terms"]},
      "errors": {"heat_rate_rel": srm_rel, "stationary_residual_norm": psrm["stationary.residual_norm"][0]}},
    "two_mass_chain": {"runs": {m: ld(f"tmc_{m}.run.json")["run"] for m in ("modal", "stationary", "transient", "half")}, "added_as_file_only": True,
      "errors": {"modal_max": tm_err_modal, "stationary_max": tm_err_stat, "transient_x_m": max(tm_full[:2]), "half_step_x_m": max(tm_half[:2]),
                 "half_step_error_ratio": list(tm_ratio), "energy_residual_j": wres(e_tt)}},
    "gas_piston": {"runs": {m: ld(f"gp_{m}.run.json")["run"] for m in ("modal", "transient", "comp")},
      "mode_verdicts": {"modal": bv_gpm, "transient": bv_gpt, "comp": bv_gpc},
      "errors": {"ruchardt_freq_rel": gp_fad_rel, "adiabatic_zeta": gp_z_ad, "equilibrium_v_m_s": gp_eq_v, "equilibrium_t_k": gp_eq_t,
                 "comp_pos_m": gp_comp_x, "comp_vel_m_s": gp_comp_v, "comp_temp_k": gp_comp_t, "comp_energy_residual_j": wres(e_gpc)}},
    "block_on_belt": {"run": ld("bob.run.json")["run"], "verdict": bv_bob,
      "errors": {"max_period_s": bob_T_err, "work_per_cycle_rel": bob_W_rel, "ambient_heat_rel": bob_amb_rel,
                 "peak_slip_speed_m_s": bob_vslip_pk, "energy_residual_j": wres(e_bob)}},
  },
  "drills": {"replay_audit_clean": True, "resume_idempotent": True, "refusals": ["missing_card", "bad_analysis", "missing_area", "srm_transient", "srm_bad_dt", "srm_bad_lx", "bob_ood_n"]},
}
pathlib.Path(rcpt_p).write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
with open(jsonl_p, "a") as w:
  w.write(json.dumps({"stage": "oracles_and_drills", "status": "pass", "receipt": rcpt_p}, sort_keys=True) + "\n")
print(f"assemblies_e2e result=pass receipt={rcpt_p}")
PY
