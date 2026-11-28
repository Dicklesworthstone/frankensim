//! Option B core tests for `fs-runtime` loading real `.fsim` assemblies:
//! full 10-example `.fsim` sweep against the shared `fs-cli` lowering registry,
//! live `set_input` work booking, `seed_override`, `render_f32` halt/silence,
//! inline probe target guards, `from_system` per-step latch checks, the weight
//! along a vertical prismatic joint, and a revolute child's swing.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fsim::{PackCatalog, admit_assembly_source};
use fs_runtime::AssemblyRuntime;
use fs_scenario::assembly::draft::{PhysicsDomain, ProbeDecl, ProbeKind};
use fs_system::elements::Motion;
use fs_system::lowering::Observed;
use fs_system::{LoweredSystem, LoweringRegistry};

const FORCED_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/forced_oscillator.fsim");
const TWO_MASS_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/two_mass_chain.fsim");
const CONDUCTION_FSIM: &[u8] =
    include_bytes!("../../../examples/assemblies/two_body_conduction.fsim");
const SRM_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/srm1453_slab.fsim");
const BELT_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/block_on_belt.fsim");
const GAS_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/gas_piston.fsim");
const FILM_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/shear_film_pair.fsim");
const RUB288_FSIM: &[u8] =
    include_bytes!("../../../examples/assemblies/rubber_mount_creep_288k.fsim");
const RUB298_FSIM: &[u8] =
    include_bytes!("../../../examples/assemblies/rubber_mount_creep_298k.fsim");
const SLIDER_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/thermo_mech_slider.fsim");

const STEEL_PACK: &[u8] = include_bytes!(
    "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
);
const AL_PACK: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/aluminum-pure-thermal-nbs-c447.fsmcdpk");
const ATKINS_PACK: &[u8] =
    include_bytes!("../../../examples/system/atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk");
const SENS_PACK: &[u8] = include_bytes!(
    "../../../examples/assemblies/packs/aluminum-pure-nbs-c447-sensible-heat.fsmcdpk"
);
const PG50_PACK: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/dowfrost-pg50-40-60c.fsmcdpk");
const SRM_PACK: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/eps-nist-srm1453-2013.fsmcdpk");
const N2_PACK: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk");
const KATS_PACK: &[u8] = include_bytes!(
    "../../../examples/assemblies/packs/katsivardi-2025-pla-metal-dry-air-interface.fsintpk"
);
const RUBBER_PACK: &[u8] = include_bytes!(
    "../../../examples/assemblies/packs/rouleau-2015-self-adhesive-rubber-shear.fsintpk"
);

const INLINE_THERMAL_FSIM: &str = r#"(assembly "inline_thermal"
  (version 1) (units "SI") (seed 0x11)
  (versions (constellation "2026-07") (workspace "gallery-v1"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state assembled)
  (environment :ambient-temp 293.15K :ambient-pressure 101325Pa)
  (time-plan :horizon 1.0s :dt 0.05s :sample-rate 20Hz)
  (part "hot" :kind solid :geometry (box 0.06m 0.05m 0.05m)
    :material "aluminum-pure-nbs-c447" :query-temp 293.15K
    :physics (thermal) :fidelity lumped)
  (source "heater" :part "hot" :port "thermal" :history (constant 500W))
  (probe "hot_temp" :kind thermocouple :target "hot" :sample-rate 20Hz)
)"#;

fn cli_registry() -> LoweringRegistry {
    let mut reg = LoweringRegistry::with_kernel_defaults();
    fs_system::lowering::owners::register_owner_defaults(&mut reg);
    reg
}

fn direct_probe_value(
    lowered: &LoweredSystem,
    draft: &fs_scenario::assembly::draft::DraftAssembly,
    decl: &ProbeDecl,
) -> f64 {
    let (eid, qty) = match &decl.kind {
        ProbeKind::Thermocouple { location, .. } => (*location, "temperature"),
        ProbeKind::FieldQuantity {
            target, quantity, ..
        } => (*target, quantity.name.as_str()),
        other => panic!("unexpected probe kind {other:?}"),
    };
    let obs = lowered.observers.get(&eid).and_then(|f| f(&lowered.system));
    let nodes = lowered
        .entity_nodes
        .get(&eid)
        .map_or(&[][..], Vec::as_slice);
    let is_th = draft
        .input()
        .parts
        .iter()
        .any(|p| p.entity == eid && p.physics.contains(&PhysicsDomain::Thermal));
    match qty {
        "position" | "displacement" => nodes
            .iter()
            .find_map(|&n| lowered.system.state::<Motion>(n).map(|m| m.position_m))
            .unwrap(),
        "velocity" => nodes
            .iter()
            .find_map(|&n| lowered.system.state::<Motion>(n).map(|m| m.velocity_m_s))
            .unwrap(),
        "temperature" => match obs {
            Some(Observed::Cell(t)) => t,
            _ if is_th => nodes
                .iter()
                .find_map(|&n| lowered.system.state::<f64>(n).copied())
                .unwrap(),
            _ => panic!("no temperature for {}", decl.name),
        },
        "heat_rate" => match obs {
            Some(Observed::HeatRate(q)) => q,
            _ => panic!("no heat_rate for {}", decl.name),
        },
        other => {
            let Some(Observed::Body {
                velocity_m_s: v,
                angular_velocity_rad_s: w,
            }) = obs
            else {
                panic!("no body observer for {other}");
            };
            let vals = [v[0], v[1], v[2], w[0], w[1], w[2]];
            let idx = [
                "velocity-x",
                "velocity-y",
                "velocity-z",
                "angular-velocity-x",
                "angular-velocity-y",
                "angular-velocity-z",
            ]
            .iter()
            .position(|k| *k == other)
            .unwrap();
            vals[idx]
        }
    }
}

fn assert_registry_parity(
    tag: &str,
    fsim: &[u8],
    packs: &[&[u8]],
    seed_override: Option<u64>,
    dt: f64,
    steps: usize,
    expected_inputs: &[&str],
) -> AssemblyRuntime {
    let mut rt = AssemblyRuntime::load(fsim, packs, seed_override).unwrap_or_else(|e| {
        panic!("{tag}: load failed: {e:?}");
    });
    let names: Vec<&str> = rt.inputs().iter().map(|i| i.name.as_str()).collect();
    assert_eq!(names, expected_inputs, "{tag} live inputs");
    let admitted =
        admit_assembly_source(fsim, &PackCatalog::from_raw_bytes(packs).unwrap()).unwrap();
    let mut direct = cli_registry()
        .lower(&admitted.draft, &admitted.card_store)
        .unwrap();
    assert_eq!(rt.block_modes(), direct.system.block_modes(), "{tag} modes");
    let (pool, gate) = (
        ArenaPool::new(ArenaConfig::default()),
        CancelGate::new_clock_free(),
    );
    let key = StreamKey {
        seed: seed_override.unwrap_or(admitted.seed),
        kernel_id: 119,
        tile: 0,
        iteration: 0,
    };
    let mut buf = vec![0.0; rt.snapshot_len()];
    let decls = &admitted.draft.input().probes;
    for k in 0..steps {
        rt.step(dt).unwrap();
        pool.scope(|a| {
            let cx = Cx::new(&gate, a, key, Budget::INFINITE, ExecMode::Deterministic);
            direct.system.step(dt, &cx).unwrap();
        });
        rt.read_probes(&mut buf).unwrap();
        for (pi, decl) in decls.iter().enumerate() {
            let off = rt.probes().probes[pi].offset;
            let got = buf[off];
            let exp = direct_probe_value(&direct, &admitted.draft, decl);
            assert_eq!(
                got.to_bits(),
                exp.to_bits(),
                "{tag} probe {} step {k}",
                decl.name
            );
        }
    }
    assert_eq!(
        rt.ledger().energy_residual_j().to_bits(),
        direct.system.ledger().energy_residual_j().to_bits(),
        "{tag} ledger residual"
    );
    rt
}

#[test]
fn all_ten_example_assemblies_match_cli_registry_or_refuse_explicitly() {
    type SweepCase<'a> = (&'a str, &'a [u8], &'a [&'a [u8]], f64, usize, &'a [&'a str]);
    let cases: &[SweepCase<'_>] = &[
        (
            "forced_oscillator",
            FORCED_FSIM,
            &[STEEL_PACK],
            1.0e-3,
            20,
            &["load"],
        ),
        (
            "two_mass_chain",
            TWO_MASS_FSIM,
            &[STEEL_PACK],
            1.0e-3,
            20,
            &["load"],
        ),
        (
            "gas_piston",
            GAS_FSIM,
            &[STEEL_PACK, N2_PACK],
            5.0e-4,
            15,
            &["ambient_pressure"],
        ),
        (
            "two_body_conduction",
            CONDUCTION_FSIM,
            &[AL_PACK, ATKINS_PACK],
            0.05,
            15,
            &[],
        ),
        ("srm1453_slab", SRM_FSIM, &[SRM_PACK], 0.01, 10, &[]),
        (
            "block_on_belt",
            BELT_FSIM,
            &[STEEL_PACK, KATS_PACK],
            1.0e-3,
            20,
            &[],
        ),
        (
            "shear_film_pair",
            FILM_FSIM,
            &[SENS_PACK, PG50_PACK],
            0.01,
            15,
            &[],
        ),
        (
            "rubber_mount_creep_288k",
            RUB288_FSIM,
            &[STEEL_PACK, SENS_PACK, RUBBER_PACK],
            0.25,
            2,
            &[],
        ),
        (
            "rubber_mount_creep_298k",
            RUB298_FSIM,
            &[STEEL_PACK, SENS_PACK, RUBBER_PACK],
            0.25,
            2,
            &[],
        ),
    ];
    for &(tag, fsim, packs, dt, steps, ins) in cases {
        assert_registry_parity(tag, fsim, packs, None, dt, steps, ins);
    }
    // 10th example: lag-filtered thermocouple refuses explicitly; without lag it matches 0-ULP.
    assert_eq!(
        AssemblyRuntime::load(SLIDER_FSIM, &[STEEL_PACK], None)
            .err()
            .unwrap()
            .code,
        "probe-unsupported"
    );
    let no_lag = std::str::from_utf8(SLIDER_FSIM)
        .unwrap()
        .replace(" :time-constant 0.005s", "");
    assert_registry_parity(
        "thermo_mech_slider_no_lag",
        no_lag.as_bytes(),
        &[STEEL_PACK],
        None,
        1.0e-3,
        15,
        &[],
    );
}

#[test]
fn forced_oscillator_live_work_seed_override_and_no_source_parity() {
    let mut rt = assert_registry_parity(
        "forced_oscillator",
        FORCED_FSIM,
        &[STEEL_PACK],
        None,
        1.0e-3,
        50,
        &["load"],
    );
    assert!(rt.probes_json().contains("\"name\":\"mass_pos\""));
    assert_eq!(rt.inputs()[0].unit, "N");

    let drive_work = |r: &AssemblyRuntime| {
        let a = r
            .ledger()
            .nodes
            .iter()
            .find(|n| n.name == "runtime.runtime_input:load")
            .unwrap();
        (a.port_energy_in_j.value(), a.stored_change_j.value())
    };
    let (p0, s0) = drive_work(&rt);
    rt.step(1.0e-3).unwrap();
    assert!(rt.ledger().energy_residual_j().abs() < 1.0e-12);
    rt.set_input_by_name("load", -5.0).unwrap();
    rt.step(1.0e-3).unwrap();
    assert!(rt.ledger().energy_residual_j().abs() < 1.0e-12);
    for _ in 0..29 {
        rt.step(1.0e-3).unwrap();
    }
    let (p1, s1) = drive_work(&rt);
    assert!((p1 - p0).abs() > 1.0e-3 && ((p1 - p0) - (s1 - s0)).abs() < 1.0e-12);
    assert!(rt.ledger().energy_residual_j().abs() < 1.0e-12);
    assert!(rt.ledger_json().contains("frankensim.window_ledger.v1"));

    // Seed override re-admits cleanly and matches direct stepping under the overridden seed.
    assert_registry_parity(
        "forced_seed_override",
        FORCED_FSIM,
        &[STEEL_PACK],
        Some(0x9999),
        1.0e-3,
        10,
        &["load"],
    );

    // Graph with no constant source keeps identical solve modes & zero live inputs.
    let no_src: String = std::str::from_utf8(FORCED_FSIM)
        .unwrap()
        .lines()
        .filter(|l| !l.contains("(source \"load\""))
        .collect::<Vec<_>>()
        .join("\n");
    assert_registry_parity(
        "forced_no_source",
        no_src.as_bytes(),
        &[STEEL_PACK],
        None,
        1.0e-3,
        10,
        &[],
    );
}

#[test]
fn render_f32_parity_halt_and_probe_target_guards() {
    let (mut stepped, mut rendered) = (
        AssemblyRuntime::load(FORCED_FSIM, &[STEEL_PACK], None).unwrap(),
        AssemblyRuntime::load(FORCED_FSIM, &[STEEL_PACK], None).unwrap(),
    );
    stepped.set_input(0, 8.0).unwrap();
    rendered.set_input(0, 8.0).unwrap();
    let (mut block, mut buf) = ([0.0_f32; 16], [0.0_f64; 2]);
    assert_eq!(rendered.render_f32(0, 1.0e-3, 2, &mut block), 16);
    for s in block {
        stepped.step(1.0e-3).unwrap();
        stepped.step(1.0e-3).unwrap();
        stepped.read_probes(&mut buf).unwrap();
        assert_eq!(s.to_bits(), (buf[0] as f32).to_bits());
    }
    let mut bad = [1.0_f32; 4];
    assert_eq!(rendered.render_f32(0, -1.0, 1, &mut bad), 0);
    assert!(bad.iter().all(|s| s.to_bits() == 0) && rendered.is_halted());
    assert_eq!(rendered.take_refusal().unwrap().code, "invalid-dt");
    assert!(rendered.take_refusal().is_none());

    // Inline thermal fixture: thermal-only part refuses `position`; 1D part refuses `velocity-x`.
    assert!(AssemblyRuntime::load(INLINE_THERMAL_FSIM.as_bytes(), &[AL_PACK], None).is_ok());
    let bad_th = INLINE_THERMAL_FSIM.replace(
        "(probe \"hot_temp\" :kind thermocouple :target \"hot\" :sample-rate 20Hz)",
        "(probe \"p\" :kind field-quantity :target \"hot\" :quantity \"position\" :dims (1 0 0 0 0 0) :sample-rate 20Hz)",
    );
    assert_eq!(
        AssemblyRuntime::load(bad_th.as_bytes(), &[AL_PACK], None)
            .err()
            .unwrap()
            .code,
        "probe-unlowered"
    );
    let bad_vx = std::str::from_utf8(FORCED_FSIM).unwrap().replace(
        ":quantity \"velocity\" :dims (1 0 -1 0 0 0)",
        ":quantity \"velocity-x\" :dims (1 0 -1 0 0 0)",
    );
    assert_eq!(
        AssemblyRuntime::load(bad_vx.as_bytes(), &[STEEL_PACK], None)
            .err()
            .unwrap()
            .code,
        "probe-unlowered"
    );
}

#[test]
fn from_system_latches_inputs_per_step_and_rejects_bad_bindings() {
    let rendered = AssemblyRuntime::load(FORCED_FSIM, &[STEEL_PACK], None).unwrap();
    let (admitted, cell) = (
        admit_assembly_source(
            FORCED_FSIM,
            &PackCatalog::from_raw_bytes(&[STEEL_PACK]).unwrap(),
        )
        .unwrap(),
        std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0)),
    );
    let sys = || {
        cli_registry()
            .lower(&admitted.draft, &admitted.card_store)
            .unwrap()
            .system
    };
    let good_in = rendered.inputs()[0].clone();
    let mut latched_rt = AssemblyRuntime::from_system(
        "latch",
        sys(),
        vec![fs_runtime::InputBinding {
            descriptor: good_in.clone(),
            cell: std::sync::Arc::clone(&cell),
        }],
        vec![],
        1,
    )
    .unwrap();
    assert_eq!(
        cell.load(std::sync::atomic::Ordering::Relaxed),
        good_in.default_value.to_bits()
    );
    latched_rt.set_input(0, 42.0).unwrap();
    assert_eq!(
        cell.load(std::sync::atomic::Ordering::Relaxed),
        good_in.default_value.to_bits()
    );
    latched_rt.step(1.0e-3).unwrap();
    assert_eq!(
        cell.load(std::sync::atomic::Ordering::Relaxed),
        42.0_f64.to_bits()
    );
    let mut bad_in = good_in;
    bad_in.width = 0;
    let ins = vec![fs_runtime::InputBinding {
        descriptor: bad_in,
        cell,
    }];
    assert_eq!(
        AssemblyRuntime::from_system("a", sys(), ins, vec![], 1)
            .err()
            .unwrap()
            .code,
        "invalid-input"
    );
    let mut bad_pr = rendered.probes().probes[0].clone();
    (bad_pr.len, bad_pr.shape) = (2, vec![3]);
    let prs = vec![fs_runtime::ProbeBinding {
        descriptor: bad_pr,
        reader: Box::new(|_, d| d.fill(0.0)),
    }];
    assert_eq!(
        AssemblyRuntime::from_system("a", sys(), vec![], prs, 1)
            .err()
            .unwrap()
            .code,
        "invalid-probe"
    );
}

/// The child of a vertical prismatic joint takes its weight `m (g . n)`, which
/// the joint's receipt records. Free, it falls at `a = g . n` with `x = a t^2 /
/// 2`, exact to rounding because the midpoint rule is exact for a constant
/// force. On a spring of 200 N/m, damped near critically, it settles at the sag
/// `m (g . n) / k`.
#[test]
fn a_vertical_prismatic_child_falls_free_and_sags_on_a_spring() {
    let fsim = |joint: &str| {
        format!(
            r#"(assembly "vertical"
  (version 1) (units "SI") (seed 0x23)
  (versions (constellation "2026-07") (workspace "gallery-v1"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded) (environment earth-lab-air)
  (time-plan :horizon 5s :dt 0.01s :sample-rate 100Hz)
  (part "mass" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m)
    :material "steel-music-wire-ligo-p1400229" :query-temp 293.15K
    :physics (mechanics) :fidelity lumped)
  (joint "slide" :kind prismatic :parent root :child "mass" :axis (0 0 1){joint})
  (probe "mass_pos" :kind field-quantity :target "mass"
    :quantity "position" :dims (1 0 0 0 0 0) :sample-rate 100Hz)
)"#
        )
    };
    let (dt, a, k) = (0.01, -9.806_65, 200.0);
    let packs = PackCatalog::from_raw_bytes(&[STEEL_PACK]).unwrap();
    for (joint, steps) in [("", 100), (" :damping 40N*s/m :stiffness 200N/m", 500)] {
        let text = fsim(joint);
        let admitted = admit_assembly_source(text.as_bytes(), &packs).unwrap();
        let lowered = cli_registry().lower(&admitted.draft, &admitted.card_store);
        let receipt = &lowered.unwrap().receipt.joints[0].parameter_sources;
        let value = |name: &str| receipt.iter().find(|p| p.parameter == name).unwrap().value;
        let m = value("mass");
        let recorded = [value("gravity[2]"), value("axis_unit[2]"), value("weight")];
        let bits = |v: [f64; 3]| v.map(f64::to_bits);
        assert_eq!(bits(recorded), bits([a, 1.0, m * a]), "{recorded:?}");
        let mut rt = AssemblyRuntime::load(text.as_bytes(), &[STEEL_PACK], None).unwrap();
        let (mut buf, at) = (vec![0.0; rt.snapshot_len()], rt.probes().probes[0].offset);
        let near = |x: f64, exact: f64| (x - exact).abs() <= 1.0e-12 * exact.abs();
        for step in 1..=steps {
            rt.step(dt).unwrap();
            rt.read_probes(&mut buf).unwrap();
            let t = f64::from(step) * dt;
            let fall = a * t * t / 2.0;
            if joint.is_empty() {
                assert!(near(buf[at], fall), "step {step}: {} m", buf[at]);
            }
        }
        let sag = m * a / k;
        if !joint.is_empty() {
            assert!(near(buf[at], sag), "sag {sag} m: {} m", buf[at]);
        }
        assert!(rt.ledger().energy_residual_j().abs() < 1.0e-9, "{joint}");
    }
}

/// A revolute box turns about its centroid with the receipt's `rho n^T J n`. On
/// a spring at rest at its turn `theta0` it swings `theta_n = theta0 + (w0 /
/// omega) sin(n phi)`, `tan(phi / 2) = omega dt / 2`, also on a falling platform.
#[test]
fn a_revolute_child_swings_about_its_centroid_on_a_torsional_spring() {
    let ([lx, ly, lz], half, w0, dt) = ([0.5, 0.4, 0.1], 0.3_f64, 10.0, 0.001);
    let fsim = |n: [f64; 3], parent: &str, extra: &str| {
        let v = |s: f64, unit: &str| n.map(|a| format!("{:?}{unit}", s * a)).join(" ");
        let (axis, spin) = (v(1.0, ""), v(w0, "rad/s"));
        let turn = format!("{:?} {}", half.cos(), v(half.sin(), ""));
        format!(
            r#"(assembly "rotor"
  (version 1) (units "SI") (seed 0x23)
  (versions (constellation "2026-07") (workspace "gallery-v1"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded) (environment earth-lab-air)
  (time-plan :horizon 1s :dt 0.001s :sample-rate 1000Hz)
  (part "rotor" :kind rigid-body :geometry (box {lx}m {ly}m {lz}m)
    :material "steel-music-wire-ligo-p1400229" :query-temp 293.15K
    :physics (mechanics) :fidelity lumped
    :initial (:position (0m 0m 0m) :orientation ({turn})
      :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity ({spin})))
  (probe "turn" :kind field-quantity :target "rotor"
    :quantity "angle" :dims (0 0 0 0 0 0) :sample-rate 1000Hz){extra}
  (joint "pin" :kind revolute :parent {parent} :child "rotor" :axis ({axis})
    :stiffness 2000N*m/rad)
)"#
        )
    };
    let platform = r#"
  (part "platform" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m)
    :material "steel-music-wire-ligo-p1400229" :query-temp 293.15K
    :physics (mechanics) :fidelity lumped)
  (joint "slide" :kind prismatic :parent root :child "platform" :axis (0 0 1))
  (probe "fall" :kind field-quantity :target "platform"
    :quantity "position" :dims (1 0 0 0 0 0) :sample-rate 1000Hz)"#;
    let packs = PackCatalog::from_raw_bytes(&[STEEL_PACK]).unwrap();
    let lower = |text: &str| {
        let admitted = admit_assembly_source(text.as_bytes(), &packs).unwrap();
        cli_registry().lower(&admitted.draft, &admitted.card_store)
    };
    let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let sums = [ly * ly + lz * lz, lx * lx + lz * lz, lx * lx + ly * ly];
    let on = "\"platform\"";
    for (i, n, parent, extra) in [(0, x, "root", ""), (1, y, "root", ""), (2, z, on, platform)] {
        let text = fsim(n, parent, extra);
        let receipt = lower(&text).unwrap().receipt;
        let rotor = &receipt.components[0].parameter_sources;
        let value = |name: &str| rotor.iter().find(|p| p.parameter == name).unwrap().value;
        let [rho, v, inertia] = ["density", "volume", "moment_of_inertia"].map(value);
        let exact = rho * (v * sums[i] / 12.0);
        assert_eq!(inertia.to_bits(), exact.to_bits(), "{n:?}");
        let (theta0, omega) = (value("initial_turn_rad"), (2000.0 / inertia).sqrt());
        assert!((theta0 - 2.0 * half).abs() <= 1.0e-15, "{n:?}: {theta0}");
        let phi = 2.0 * (omega * dt / 2.0).atan();
        let mut rt = AssemblyRuntime::load(text.as_bytes(), &[STEEL_PACK], None).unwrap();
        let mut buf = vec![0.0; rt.snapshot_len()];
        let at: Vec<usize> = rt.probes().probes.iter().map(|p| p.offset).collect();
        for step in 1..=500_u32 {
            let refusal = rt.step(dt).err();
            assert!(refusal.is_none(), "{n:?} {step}: {refusal:?}");
            rt.read_probes(&mut buf).unwrap();
            let k = f64::from(step);
            let theta = theta0 + w0 / omega * (k * phi).sin();
            let (got, tol) = (buf[at[0]], 16.0 * f64::EPSILON * (theta0 + w0 / omega) * k);
            assert!((got - theta).abs() <= tol, "{n:?} {step}: {got}");
            if let Some(&j) = at.get(1) {
                let (fell, fall) = (buf[j], -9.806_65 * (k * dt) * (k * dt) / 2.0);
                assert!((fell / fall - 1.0).abs() <= 1.0e-12, "{step}: {fell}");
            }
        }
        assert!(rt.ledger().energy_residual_j().abs() < 1.0e-9, "{n:?}");
    }
    let refused = |text: &str, why: &str| {
        let refusal = format!("{:?}", lower(text).err().unwrap());
        assert!(refusal.contains(why), "{refusal}");
    };
    let slid = platform.replace(":child \"platform\"", ":child \"rotor\"");
    refused(&fsim(z, "root", &slid), "prismatic joint `slide`");
    let turning = platform.replace("prismatic", "revolute");
    let turning = turning.replace("(0 0 1))", "(0 0 1) :stiffness 1N*m/rad)");
    refused(&fsim(z, on, &turning), "its parent turns");
    let off = fsim(x, "root", "").replace(":axis (1.0 0.0 0.0)", ":axis (0.0 1.0 0.0)");
    refused(&off, "it refuses initial.orientation");
    let angled = platform.replace("\"position\" :dims (1", "\"angle\" :dims (0");
    let angle = fsim(z, on, &angled);
    let place = fsim(z, "root", "").replace("\"angle\" :dims (0", "\"position\" :dims (1");
    let spun = place.replace("\"position\" :dims (1 0 0", "\"velocity\" :dims (1 0 -1");
    for text in [angle, place, spun] {
        let load = AssemblyRuntime::load(text.as_bytes(), &[STEEL_PACK], None);
        assert_eq!(load.err().unwrap().code, "probe-unlowered");
    }
}

/// A revolute box's receipt inertia is `rho n^T J n` for the unit vector `n` of
/// its joint's axis: about the oblique `(0 0.6 0.8)` it is `0.36 J_y + 0.64 J_z`,
/// and about `(0 0 2)` and `(0 0 -1)` it is the inertia about `(0 0 1)`.
#[test]
fn a_revolute_childs_inertia_is_the_quadratic_form_of_its_unit_axis() {
    let [lx, ly, lz] = [0.5, 0.4, 0.1];
    let packs = PackCatalog::from_raw_bytes(&[STEEL_PACK]).unwrap();
    let lowered = |axis: &str| {
        let text = format!(
            r#"(assembly "rotor"
  (version 1) (units "SI") (seed 0x23)
  (versions (constellation "2026-07") (workspace "gallery-v1"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded) (environment earth-lab-air)
  (time-plan :horizon 1s :dt 0.001s :sample-rate 1000Hz)
  (part "rotor" :kind rigid-body :geometry (box {lx}m {ly}m {lz}m)
    :material "steel-music-wire-ligo-p1400229" :query-temp 293.15K
    :physics (mechanics) :fidelity lumped)
  (joint "pin" :kind revolute :parent root :child "rotor" :axis ({axis})
    :stiffness 2000N*m/rad)
)"#
        );
        let admitted = admit_assembly_source(text.as_bytes(), &packs).unwrap();
        let lowered = cli_registry().lower(&admitted.draft, &admitted.card_store);
        let receipt = lowered.unwrap().receipt;
        let rotor = &receipt.components[0].parameter_sources;
        let value = |name: &str| rotor.iter().find(|p| p.parameter == name).unwrap().value;
        ["density", "moment_of_inertia"].map(value)
    };
    let [rho, oblique] = lowered("0 0.6 0.8");
    let [j_y, j_z] = [lx * lx + lz * lz, lx * lx + ly * ly].map(|s| lx * ly * lz * s / 12.0);
    let (exact, tol) = (rho * (0.36 * j_y + 0.64 * j_z), 1.0e-15);
    assert!((oblique / exact - 1.0).abs() <= tol, "{oblique} {exact}");
    let unit = lowered("0 0 1")[1];
    for axis in ["0 0 2", "0 0 -1"] {
        assert_eq!(lowered(axis)[1].to_bits(), unit.to_bits(), "{axis}");
    }
}

/// A 15 g pawl (0.04 x 0.004 x 0.011 82 m) on 0.05 N m/rad, whose stage slope
/// I / dt is far below one, swings like the box above at 1 ms and at 0.1 ms.
/// While the monolithic solve floored that slope at one, its audit refused step
/// 33 at 1 ms and step 7 at 0.1 ms.
#[test]
fn a_pawl_sized_revolute_child_swings_on_a_soft_torsional_spring() {
    let ([lx, ly, lz], half, w0, stiffness) = ([0.04, 0.004, 0.011_82], 0.3_f64, 10.0, 0.05);
    let fsim = |dt: f64, rate: &str| {
        let (c, s) = (half.cos(), half.sin());
        format!(
            r#"(assembly "pawl"
  (version 1) (units "SI") (seed 0x23)
  (versions (constellation "2026-07") (workspace "gallery-v1"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded) (environment earth-lab-air)
  (time-plan :horizon 1s :dt {dt:?}s :sample-rate {rate}Hz)
  (part "pawl" :kind rigid-body :geometry (box {lx}m {ly}m {lz}m)
    :material "steel-music-wire-ligo-p1400229" :query-temp 293.15K
    :physics (mechanics) :fidelity lumped
    :initial (:position (0m 0m 0m) :orientation ({c:?} {s:?} 0.0 0.0)
      :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity ({w0:?}rad/s 0.0rad/s 0.0rad/s)))
  (probe "turn" :kind field-quantity :target "pawl"
    :quantity "angle" :dims (0 0 0 0 0 0) :sample-rate {rate}Hz)
  (joint "pin" :kind revolute :parent root :child "pawl" :axis (1.0 0.0 0.0)
    :stiffness {stiffness:?}N*m/rad)
)"#
        )
    };
    let packs = PackCatalog::from_raw_bytes(&[STEEL_PACK]).unwrap();
    for (dt, rate) in [(0.001, "1000"), (0.000_1, "10000")] {
        let text = fsim(dt, rate);
        let admitted = admit_assembly_source(text.as_bytes(), &packs).unwrap();
        let lowered = cli_registry()
            .lower(&admitted.draft, &admitted.card_store)
            .unwrap();
        let sources = &lowered.receipt.components[0].parameter_sources;
        let rho = sources
            .iter()
            .find(|p| p.parameter == "density")
            .unwrap()
            .value;
        let inertia = rho * (lx * ly * lz) * (ly * ly + lz * lz) / 12.0;
        let omega = (stiffness / inertia).sqrt();
        let (theta0, phi) = (2.0 * half, 2.0 * (omega * dt / 2.0).atan());
        let mut rt = AssemblyRuntime::load(text.as_bytes(), &[STEEL_PACK], None).unwrap();
        let mut buf = vec![0.0; rt.snapshot_len()];
        let at = rt.probes().probes[0].offset;
        for step in 1..=500_u32 {
            let refusal = rt.step(dt).err();
            assert!(refusal.is_none(), "{dt} {step}: {refusal:?}");
            rt.read_probes(&mut buf).unwrap();
            let k = f64::from(step);
            let theta = theta0 + w0 / omega * (k * phi).sin();
            let tol = 16.0 * f64::EPSILON * (theta0 + w0 / omega) * k;
            assert!((buf[at] - theta).abs() <= tol, "{dt} {step}: {}", buf[at]);
        }
        assert!(rt.ledger().energy_residual_j().abs() < 1.0e-9, "{dt}");
    }
}
