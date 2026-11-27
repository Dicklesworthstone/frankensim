//! Option B core tests for `fs-runtime` loading real `.fsim` assemblies:
//! 0-ULP CLI parity, solve-mode parity (with & without constant sources),
//! live `set_input` work booking, `render_f32` halt/silence, and probe target guards.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fsim::{PackCatalog, admit_assembly_source};
use fs_runtime::AssemblyRuntime;
use fs_system::LoweringRegistry;
use fs_system::elements::Motion;

const FORCED_FSIM: &[u8] = include_bytes!("../../../examples/assemblies/forced_oscillator.fsim");
const STEEL_PACK: &[u8] = include_bytes!(
    "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
);
const CONDUCTION_FSIM: &[u8] =
    include_bytes!("../../../examples/assemblies/two_body_conduction.fsim");
const AL_PACK: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/aluminum-pure-thermal-nbs-c447.fsmcdpk");

fn assert_cli_parity(fsim: &[u8], packs: &[&[u8]], rt: &mut AssemblyRuntime, steps: usize) {
    let admitted =
        admit_assembly_source(fsim, &PackCatalog::from_raw_bytes(packs).unwrap()).unwrap();
    let mut cli = LoweringRegistry::with_kernel_defaults()
        .lower(&admitted.draft, &admitted.card_store)
        .unwrap();
    assert_eq!(rt.block_modes(), cli.system.block_modes());
    let mass = cli.entity_nodes[&admitted.draft.input().parts[0].entity][0];
    let (pos_off, vel_off) = (
        rt.probes().probe("mass_pos").unwrap().offset,
        rt.probes().probe("mass_vel").unwrap().offset,
    );
    let (pool, gate) = (
        ArenaPool::new(ArenaConfig::default()),
        CancelGate::new_clock_free(),
    );
    let key = StreamKey {
        seed: admitted.seed,
        kernel_id: 119,
        tile: 0,
        iteration: 0,
    };
    let mut buf = vec![0.0; rt.snapshot_len()];
    for k in 0..steps {
        rt.step(1.0e-3).unwrap();
        pool.scope(|a| {
            let cx = Cx::new(&gate, a, key, Budget::INFINITE, ExecMode::Deterministic);
            cli.system.step(1.0e-3, &cx).unwrap();
        });
        rt.read_probes(&mut buf).unwrap();
        let m = cli.system.state::<Motion>(mass).unwrap();
        assert_eq!(buf[pos_off].to_bits(), m.position_m.to_bits(), "pos {k}");
        assert_eq!(buf[vel_off].to_bits(), m.velocity_m_s.to_bits(), "vel {k}");
    }
    assert_eq!(
        rt.ledger().energy_residual_j().to_bits(),
        cli.system.ledger().energy_residual_j().to_bits()
    );
}

#[test]
fn forced_oscillator_matches_cli_and_modes_and_books_live_work() {
    let mut rt = AssemblyRuntime::load(FORCED_FSIM, &[STEEL_PACK], None).unwrap();
    assert!(rt.probes_json().contains("\"name\":\"mass_pos\""));
    assert_eq!(
        (
            rt.inputs().len(),
            rt.inputs()[0].name.as_str(),
            rt.inputs()[0].unit.as_str()
        ),
        (1, "load", "N")
    );
    assert_cli_parity(FORCED_FSIM, &[STEEL_PACK], &mut rt, 50);

    let drive_work = |r: &AssemblyRuntime| {
        let a = r
            .ledger()
            .nodes
            .iter()
            .find(|n| n.name == "runtime.runtime_input:load")
            .unwrap();
        (a.port_energy_in_j.value(), a.stored_change_j.value())
    };
    // Step-by-step live input change: ledger closes on BOTH step 1 and step 2 (latched per step).
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

    // Graph with no constant source keeps identical solve modes & zero live inputs.
    let no_src: String = std::str::from_utf8(FORCED_FSIM)
        .unwrap()
        .lines()
        .filter(|l| !l.contains("(source \"load\""))
        .collect::<Vec<_>>()
        .join("\n");
    let mut free_rt = AssemblyRuntime::load(no_src.as_bytes(), &[STEEL_PACK], None).unwrap();
    assert!(free_rt.inputs().is_empty());
    assert_cli_parity(no_src.as_bytes(), &[STEEL_PACK], &mut free_rt, 10);
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

    // Thermal-only part refuses `position` (never leaks K as m); 1D part refuses `velocity-x`.
    let (mut kept, mut skip) = (Vec::new(), 0usize);
    for l in std::str::from_utf8(CONDUCTION_FSIM).unwrap().lines() {
        if l.contains("(part \"cold\"") {
            skip = 3;
        }
        if skip > 0 {
            skip -= 1;
            continue;
        }
        if !["atkins_joint", "6061-T6", ":area", "cold_temp", "hot_limit"]
            .iter()
            .any(|k| l.contains(k))
        {
            kept.push(l);
        }
    }
    let base_th = kept.join("\n");
    assert!(AssemblyRuntime::load(base_th.as_bytes(), &[AL_PACK], None).is_ok());
    let bad_th = base_th.replace(
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
        LoweringRegistry::with_kernel_defaults()
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
    // Staged only in pending cell until step() latches it.
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
