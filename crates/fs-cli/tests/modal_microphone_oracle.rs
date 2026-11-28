//! Baffled-pressure microphone oracles (`.11.9.5c-1`):
//! 1x & 4x `BaffledPressure::new_multirate` parity + receipted moist-air `rho`/`c`,
//! `fs-cli` `run_assembly` + `explain`, batter outward sign, arrival delay `r/c` + far-field `1/r`,
//! on-axis `(1,1)` dipole cancellation vs off-axis dipole radiation, and all 9 named refusals.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_cli::{
    AssemblyAnalysisMode, AssemblyCancelPlan, CardPackKind, CardPackSet, RawCardPack, SolveRunId,
    explain_assembly_run, run_assembly,
};
use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium, SurfaceSample};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fsim::{PackCatalog, admit_assembly_source};
use fs_ledger::Ledger;
use fs_material::gas::GasState;
use fs_nlmodal::component::ModalStructure;
use fs_ports::component::Component;
use fs_scenario::Scenario;
use fs_scenario::assembly::draft::{DraftAssembly, StepPolicy};
use fs_system::acoustic::{AcousticMicrophone, RadiatingModalSurface};
use fs_system::{LoweredSystem, LoweringRegistry, lowering::owners::register_owner_defaults};

const MYLAR: &[u8] =
    include_bytes!("../../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const HICKORY: &[u8] =
    include_bytes!("../../fs-vkplate/packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");
type MState = <ModalStructure as Component>::State;

fn fsim_src(dt: f64, sr_plan: u32, extra: &str, probes: &str) -> String {
    format!(
        "(assembly \"drum_mic\" (version 1) (units \"SI\") (seed 0x7) \
         (versions (constellation \"2026-07\") (workspace \"t\")) \
         (budget (wall 60s) (mem 64MiB)) \
         (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
         (reference-state assembled) \
         (environment :ambient-temp 293.15K :ambient-pressure 90000Pa :relative-humidity 0) \
         (time-plan :horizon 0.005s :dt {dt:.17}s :sample-rate {sr_plan}Hz) \
         (part \"stick\" :kind rigid-body :geometry (sphere 0.015m) :material \"carya-ovata\" \
           :query-at (T 293.15K moisture_content 0.12) :physics (mechanics) :fidelity lumped \
           :initial (:position (0m 0m 0m) :orientation (1 0 0 0) \
             :linear-velocity (-2m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s))) \
         (part \"batter\" :kind shell :geometry (tensioned-disk :radius 0.1778m :thickness 0.00025m \
           :tension 2500N/m :radial-intervals 6 :azimuths 24) :material \"polyethylene-terephthalate\" \
           :query-temp 293.15K :physics (mechanics) :fidelity (modal 3)) \
         (interface \"hit\" :kind hertz-contact :a \"stick\" :b \"batter\" \
           :maximum-force 100000N :maximum-penetration 0.01m :gap 0.0002m :station-x 0.0889m :station-y 0m) \
         {extra} {probes})"
    )
}

fn lower(packs: &PackCatalog, src: &str) -> Result<(DraftAssembly, LoweredSystem), String> {
    let adm = admit_assembly_source(src.as_bytes(), packs).map_err(|e| e.what)?;
    let mut reg = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut reg);
    let low = reg
        .lower(&adm.draft, &adm.card_store)
        .map_err(|e| e.to_string())?;
    Ok((adm.draft, low))
}

fn couple_mic(s: &RadiatingModalSurface, mr: u32, gas: &GasState) -> BaffledPressure {
    let surf: Vec<_> = s
        .samples
        .iter()
        .map(|p| SurfaceSample {
            position_m: [p.xy_m[0], p.xy_m[1], 0.0],
            area_m2: p.area_m2,
            mode_shape: p.shapes.clone(),
        })
        .collect();
    BaffledPressure::new_multirate(
        &surf,
        s.mode_count,
        mr,
        48_000,
        [0.08, 0.05, 0.35],
        RayleighMedium::from(gas),
    )
    .unwrap()
}

fn first_arrival(xs: &[f64]) -> (usize, f64, f64) {
    let pk = xs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let (i, &v) = xs
        .iter()
        .enumerate()
        .find(|(_, v)| v.abs() > 1.0e-3 * pk)
        .unwrap();
    (i, v, pk)
}

fn step_with(steps: usize, mut f: impl FnMut(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    let key = StreamKey {
        seed: 7,
        kernel_id: 5,
        tile: 0,
        iteration: 0,
    };
    ArenaPool::new(ArenaConfig::default()).scope(|a| {
        let cx = Cx::new(&gate, a, key, Budget::INFINITE, ExecMode::Deterministic);
        for _ in 0..steps {
            f(&cx);
        }
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn parity_multirate_cli_sign_geometry_dipole_and_refusals() {
    let packs = PackCatalog::from_raw_bytes(&[MYLAR, HICKORY]).unwrap();
    let gas = GasState::try_new_moist_air(293.15, 90_000.0, 0.0).unwrap();
    let (bound_par, bound_inv_r, bound_dip) = (1.0e-14_f64, 0.05_f64, 1.0e-12_f64);
    let mic_p = "(probe \"m\" :kind microphone :target \"batter\" :position (0.08m 0.05m 0.35m) :sample-rate 48000Hz)";
    let mk_mic = |d: &DraftAssembly, l: &LoweredSystem, i: usize, dt_s: f64| {
        AcousticMicrophone::from_draft_and_lowered(&d.input().probes[i], d, l, dt_s).unwrap()
    };
    // 1. 1x & 4x parity vs fs_couple BaffledPressure::new_multirate + receipted rho/c.
    for ratio in [1u32, 4u32] {
        let (mr, dt) = (48_000 * ratio, 1.0 / f64::from(48_000 * ratio));
        let (draft, mut low) = lower(&packs, &fsim_src(dt, mr, "", mic_p)).unwrap();
        low.receipt.verify_completeness().unwrap();
        let ap = &low.receipt.acoustic_members[0].parameter_sources;
        let bits = |n| {
            ap.iter()
                .find(|p| p.parameter == n)
                .unwrap()
                .value
                .to_bits()
        };
        assert_eq!(bits("exterior_air_density"), gas.density.to_bits());
        assert_eq!(bits("exterior_air_sound_speed"), gas.sound_speed.to_bits());
        let surf = &low.radiating_surfaces[&draft.input().parts[1].entity];
        let (mut mic, mut rmic) = (mk_mic(&draft, &low, 0, dt), couple_mic(surf, mr, &gas));
        let (mut err, mut pk, mut tr) = (0.0_f64, 0.0_f64, vec![0.0; ratio as usize * 3]);
        step_with(240, |cx| {
            for s in 0..ratio as usize {
                low.system.step(dt, cx).unwrap();
                let st = low.system.state::<MState>(surf.node).unwrap();
                for i in 0..3 {
                    tr[s * 3 + i] = st[i].velocity_m_sqrt_kg_per_s;
                }
                if let Some(got) = mic.advance(&low.system).unwrap() {
                    let exp = rmic.step_trace(&tr).unwrap();
                    (err, pk) = (err.max((got - exp).abs()), pk.max(exp.abs()));
                }
            }
        });
        let rel = err / pk;
        eprintln!("MIC PARITY {ratio}x: rel={rel:.2e} (bound={bound_par:.1e}) peak={pk:.3}Pa");
        assert!(pk > 1.0 && rel <= bound_par);
    }

    // 2. End-to-end fs-cli run_assembly + explain_assembly_run microphone probe series.
    let dt = 1.0 / 48_000.0;
    let cp = CardPackSet::admit(
        [("m.fsmatpk", MYLAR), ("h.fsmatpk", HICKORY)]
            .into_iter()
            .map(|(s, b)| RawCardPack {
                kind: CardPackKind::Material,
                source: s.into(),
                bytes: b.to_vec(),
                expect: None,
            })
            .collect(),
    )
    .unwrap();
    let (ledger, mut t) = (Ledger::open(":memory:").unwrap(), 0.0);
    let cli_s = fsim_src(dt, 48_000, "", mic_p);
    let (out, _) = run_assembly(
        &ledger,
        cli_s.as_bytes(),
        &cp,
        AssemblyAnalysisMode::Transient,
        &mut || {
            t += 0.01;
            t
        },
        AssemblyCancelPlan::default(),
    )
    .unwrap();
    let exp = explain_assembly_run(&ledger, SolveRunId::parse_hex(&out.run).unwrap()).unwrap();
    let ps_json = exp.probe_series_json.unwrap();
    assert!(ps_json.contains("\"kind\":\"microphone\"") && ps_json.contains("\"unit\":\"Pa\""));

    // 3. Batter positive outward push (+z toward mic) -> strictly positive pressure at t0 = r/c.
    let push = fsim_src(
        dt,
        48_000,
        "(source \"up\" :part \"batter\" :port \"face\" :history (constant 50N))",
        mic_p,
    )
    .replace("-2m/s", "0m/s");
    let (pd, mut plow) = lower(&packs, &push).unwrap();
    let (mut pm, mut ps) = (mk_mic(&pd, &plow, 0, dt), Vec::new());
    step_with(120, |cx| {
        plow.system.step(dt, cx).unwrap();
        ps.push(pm.advance(&plow.system).unwrap().unwrap());
    });
    let (ip, sp, _) = first_arrival(&ps);
    eprintln!("BATTER OUTWARD SIGN: p={sp:+.3}Pa @ {ip}");
    assert!(sp > 0.0);

    // 4. Arrival delay r/c, far-field 1/r falloff, and on-axis (1,1) dipole cancellation.
    let far_p = "(probe \"m1\" :kind microphone :target \"batter\" :position (0m 0m 1.0m) :sample-rate 48000Hz) \
        (probe \"m2\" :kind microphone :target \"batter\" :position (0m 0m 2.0m) :sample-rate 48000Hz) \
        (probe \"on\" :kind microphone :target \"batter\" :position (0m 0m 0.5m) :sample-rate 48000Hz) \
        (probe \"off\" :kind microphone :target \"batter\" :position (0.3m 0m 0.5m) :sample-rate 48000Hz)";
    let (fd, mut flow) = lower(
        &packs,
        &fsim_src(dt, 48_000, "", far_p).replace("0.0889m", "0m"),
    )
    .unwrap();
    let [mut m1, mut m2, mut mon, mut moff] = [0, 1, 2, 3].map(|i| mk_mic(&fd, &flow, i, dt));
    let (mut y1, mut y2) = (Vec::new(), Vec::new());
    step_with(420, |cx| {
        flow.system.step(dt, cx).unwrap();
        y1.push(m1.advance(&flow.system).unwrap().unwrap());
        y2.push(m2.advance(&flow.system).unwrap().unwrap());
    });
    let (k1, k2) = (
        (48_000.0 / gas.sound_speed).floor() as usize,
        (96_000.0 / gas.sound_speed).floor() as usize,
    );
    assert!(y1[..k1].iter().all(|&v| v == 0.0) && y2[..k2].iter().all(|&v| v == 0.0));
    let ((o1, _, pk1), (o2, _, pk2)) = (first_arrival(&y1), first_arrival(&y2));
    let lag_err = (o2 as isize - o1 as isize - (k2 - k1) as isize).abs();
    let inv_r_err = (pk1 / pk2 - 2.0).abs();
    let (mut p_on, mut p_off) = (0.0_f64, 0.0_f64);
    for k in 0..120 {
        let v = if k == 0 { [0.0, 1.0, 0.5] } else { [0.0; 3] };
        p_on = p_on.max(mon.step_velocities(&v).unwrap().abs());
        p_off = p_off.max(moff.step_velocities(&v).unwrap().abs());
    }
    let dip_rel = p_on / p_off;
    eprintln!(
        "GEOM/DIPOLE: lag_err={lag_err} inv_r_err={inv_r_err:.4} (bound={bound_inv_r}) dip_rel={dip_rel:.2e} (bound={bound_dip:.1e}) p_off={p_off:.2}Pa"
    );
    assert!(lag_err <= 2 && inv_r_err <= bound_inv_r && p_off > 1.0 && dip_rel <= bound_dip);

    // 5. Named refusals: target, z<=0, z<0.05, pattern, decimation ratio, Nyquist, adaptive policy.
    let bad = |s: &str, want: &str| {
        let got = lower(&packs, s).err().unwrap();
        assert!(got.contains(want), "got `{got}`, want `{want}`");
    };
    let mpr = |tgt, z, pat, sr| {
        format!(
            "(probe \"m\" :kind microphone :target \"{tgt}\" :position (0m 0m {z}m) {pat} :sample-rate {sr}Hz)"
        )
    };
    for (tgt, z, pat, want) in [
        ("stick", "0.5", "", "not a radiating modal surface"),
        ("batter", "0", "", "z must be > 0"),
        ("batter", "-0.2", "", "z must be > 0"),
        ("batter", "0.02", "", ">= 0.05 m"),
        ("batter", "0.5", ":pattern cardioid", "omnidirectional"),
    ] {
        bad(&fsim_src(dt, 48_000, "", &mpr(tgt, z, pat, 48_000)), want);
    }
    bad(
        &fsim_src(1.0 / 50_000.0, 50_000, "", mic_p),
        "decimation ratio",
    );
    bad(
        &fsim_src(1.0 / 960_000.0, 960_000, "", mic_p),
        "decimation ratio",
    );
    let nyq = fsim_src(1.0 / 8_000.0, 48_000, "", &mpr("batter", "0.5", "", 8_000))
        .replace(":tension 2500N/m", ":tension 500000N/m");
    bad(&nyq, "Nyquist guard");
    let adm = admit_assembly_source(fsim_src(dt, 48_000, "", mic_p).as_bytes(), &packs).unwrap();
    let mut inp = adm.draft.input().clone();
    let StepPolicy::Fixed { dt: dt_q } = inp.time_plan.step_policy else {
        panic!()
    };
    inp.time_plan.step_policy = StepPolicy::Adaptive {
        dt_min: dt_q,
        dt_max: dt_q,
        dt_initial: dt_q,
        relative_tol: 1.0e-4,
        order: 2,
    };
    let sc = Scenario::new(&inp.name, 1, inp.environment.base.clone());
    let ad = DraftAssembly::admit(inp, &sc).unwrap();
    let mut reg = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut reg);
    assert!(
        reg.lower(&ad, &adm.card_store)
            .err()
            .unwrap()
            .to_string()
            .contains("fixed step")
    );
}
