//! R1d-1 tests for `AssemblyRuntime::load`: 0-ULP CLI parity on `forced_oscillator.fsim`
//! and thermal/3D probe target guards on `two_body_conduction.fsim`.

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

#[test]
fn load_matches_cli_bit_for_bit_and_enforces_probe_target_guards() {
    let mut rt = AssemblyRuntime::load(FORCED_FSIM, &[STEEL_PACK], None).unwrap();
    let cat = PackCatalog::from_raw_bytes(&[STEEL_PACK]).unwrap();
    let admitted = admit_assembly_source(FORCED_FSIM, &cat).unwrap();
    let mut cli = LoweringRegistry::with_kernel_defaults()
        .lower(&admitted.draft, &admitted.card_store)
        .unwrap();
    assert_eq!(rt.block_modes(), cli.system.block_modes());
    let mass = cli.entity_nodes[&admitted.draft.input().parts[0].entity][0];
    let pos_off = rt.probes().probe("mass_pos").unwrap().offset;
    let vel_off = rt.probes().probe("mass_vel").unwrap().offset;
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
    for k in 0..30 {
        rt.step(1.0e-3).unwrap();
        pool.scope(|a| {
            let cx = Cx::new(&gate, a, key, Budget::INFINITE, ExecMode::Deterministic);
            cli.system.step(1.0e-3, &cx).unwrap();
        });
        rt.read_probes(&mut buf).unwrap();
        let m = cli.system.state::<Motion>(mass).unwrap();
        assert_eq!(
            (buf[pos_off].to_bits(), buf[vel_off].to_bits()),
            (m.position_m.to_bits(), m.velocity_m_s.to_bits()),
            "{k}"
        );
    }
    assert_eq!(
        rt.ledger().energy_residual_j().to_bits(),
        cli.system.ledger().energy_residual_j().to_bits()
    );

    // Thermal-only part refuses `position` (never leaks K as m); 1D part refuses `velocity-x`.
    let (mut kept, mut skip) = (Vec::new(), 0usize);
    let drop_keys = ["atkins_joint", "6061-T6", ":area", "cold_temp", "hot_limit"];
    for l in std::str::from_utf8(CONDUCTION_FSIM).unwrap().lines() {
        if l.contains("(part \"cold\"") {
            skip = 3;
        }
        if skip > 0 {
            skip -= 1;
            continue;
        }
        if !drop_keys.iter().any(|k| l.contains(k)) {
            kept.push(l);
        }
    }
    let base_th = kept.join("\n");
    assert!(AssemblyRuntime::load(base_th.as_bytes(), &[AL_PACK], None).is_ok());
    let bad_th = base_th.replace(
        "(probe \"hot_temp\" :kind thermocouple :target \"hot\" :sample-rate 20Hz)",
        "(probe \"p\" :kind field-quantity :target \"hot\" :quantity \"position\" :dims (1 0 0 0 0 0) :sample-rate 20Hz)",
    );
    let err_th = AssemblyRuntime::load(bad_th.as_bytes(), &[AL_PACK], None)
        .err()
        .unwrap();
    assert_eq!(err_th.code, "probe-unlowered");
    let bad_vx = std::str::from_utf8(FORCED_FSIM).unwrap().replace(
        ":quantity \"velocity\" :dims (1 0 -1 0 0 0)",
        ":quantity \"velocity-x\" :dims (1 0 -1 0 0 0)",
    );
    let err_vx = AssemblyRuntime::load(bad_vx.as_bytes(), &[STEEL_PACK], None)
        .err()
        .unwrap();
    assert_eq!(err_vx.code, "probe-unlowered");
}
