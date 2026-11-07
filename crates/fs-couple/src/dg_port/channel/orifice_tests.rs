//! The sharp-edged tone-hole orifice terminal (`DgTerminalChannelKind::Orifice`): `dp = K |Q_prev| Q`.

use super::*;
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::GasSpec;
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::MimoLoadLimits;
use fs_vfit::model::PoleTerm;

const DT: f64 = 1.0 / 96000.0;
const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const SIDE: PortId = PortId(11);

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn duct(gas: &GasState, length: f64, elements: usize, inlet: PortId, outlet: PortId) -> DgDuct {
    DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length }] },
        gas, DgMesh { elements, order: 2 },
        DgBoundary::Port(inlet), DgBoundary::Port(outlet),
    ).unwrap()
}

fn model() -> MimoModel {
    MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Real { pole: -400.0, residue: 0.0 },
            PoleTerm::Pair { pole: C64::new(-900.0, 2000.0), residue: C64::ZERO }],
        residues: vec![
            vec![C64::from_re(2e9), C64::from_re(5e8), C64::from_re(5e8), C64::from_re(1.8e9)],
            vec![C64::new(2e8, 1e8), C64::new(2e7, 1e7), C64::new(2e7, 1e7), C64::new(1.5e8, 5e7)],
        ],
        d: vec![2e7, 5e6, 5e6, 3e7], e: vec![0.0; 4], reciprocal: true,
    }
}

/// A bore with one finite chimney (the orifice's hole) on `SIDE` and the axial outlet, on the clock `dt`.
fn load(dt: f64) -> DgPortLoad {
    let main = duct(&air(), 0.05, 2, INLET, OUTLET)
        .with_side_ports(&[DgSidePort { id: SIDE, downstream_element: 1 }]).unwrap();
    let branches = vec![DgBranch { port: SIDE, duct: duct(&air(), 0.012, 1, PortId(0), PortId(1)) }];
    let terminal = PreparedMimoLoad::new(&model(), dt, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / dt], sampled_passivity_points: 129,
        max_states: 16, max_memory_bytes: 1 << 20,
    }).unwrap();
    DgPortLoad::new_with_terminal_branches(main, dt, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, vec![SIDE, OUTLET], terminal, branches).unwrap()
}

fn orifice(port: PortId, vena_contracta: f64) -> DgTerminalChannelInput {
    DgTerminalChannelInput { port, kind: DgTerminalChannelKind::Orifice { vena_contracta } }
}

/// `C_vc` small enough that `K |Q|` is a visible fraction of the group's impedance at these flows.
const C_VC: f64 = 0.05;
/// A 1 kHz inlet flow of this amplitude [m^3/s].
const AMPLITUDE: f64 = 4e-5;

fn drive(sample: usize, dt: f64) -> f64 {
    AMPLITUDE * (core::f64::consts::TAU * 1000.0 * sample as f64 * dt).sin()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 1e-9 * scale.max(1e-300),
        "actual={actual}, expected={expected}, scale={scale}");
}

#[test]
fn the_first_sample_has_no_resistance_and_is_the_plain_terminals_sample_bit_for_bit() {
    let mut plain = load(DT);
    let mut with = load(DT).with_terminal_channels(vec![orifice(SIDE, C_VC)], 0).unwrap();
    assert_eq!(with.port_response().unwrap(), plain.port_response().unwrap());
    let q = drive(1, DT);
    let a = plain.preview_step(q).unwrap().accept();
    let b = with.preview_step(q).unwrap().accept();
    assert_eq!(a.port_outputs, b.port_outputs);
    assert_eq!(plain.last_terminal_frames(), with.last_terminal_frames());
    // The first sample's drop is zero: Q_prev was zero.
    assert_eq!(with.last_terminal_channel_frames()[0].constitutive_pressure_drop_pa, 0.0);
}

#[test]
fn each_sample_solves_the_drop_of_the_previous_flow_and_every_audit_and_the_ledger_close() {
    let mut load = load(DT).with_terminal_channels(vec![orifice(SIDE, C_VC)], 0).unwrap();
    let area = load.terminal_area_m2(SIDE).unwrap().unwrap();
    let mut previous = 0.0_f64;
    let (mut work, mut dissipated_from_frames) = (0.0_f64, 0.0_f64);
    let mut peak = 0.0_f64;
    for sample in 1..=480 {
        let output = load.preview_step(drive(sample, DT)).unwrap().accept();
        let frame = load.last_terminal_channel_frames()[0];
        let rho = load.terminal_channels()[0].gas().density;
        let k = rho / (2.0 * C_VC * C_VC * area * area);
        // The drop the sample was solved with is K |Q_prev| Q, and it is the observed duct-minus-exterior drop.
        near(frame.constitutive_pressure_drop_pa, k * previous.abs() * frame.outward_flow_m3_s,
            frame.constitutive_pressure_drop_pa.abs() + 1.0);
        // (The engine's own same-step audit holds this to 1e-10 of the whole pressure scale; here a 1e3 Pa floor.)
        near(frame.mean_pressure_drop_pa, frame.constitutive_pressure_drop_pa,
            frame.mean_pressure_drop_pa.abs() + frame.constitutive_pressure_drop_pa.abs() + 1e3);
        // Passive: a nonnegative resistance never returns work.
        assert!(frame.pressure_work_j >= 0.0, "sample {sample}: {}", frame.pressure_work_j);
        near(frame.pressure_work_j, frame.constitutive_pressure_drop_pa * frame.outward_flow_m3_s * DT,
            frame.pressure_work_j.abs() + 1e-30);
        assert!(frame.step.is_none() && frame.flow_residual_m3_s == 0.0);
        work += frame.pressure_work_j;
        dissipated_from_frames += frame.pressure_work_j;
        previous = frame.outward_flow_m3_s;
        peak = peak.max(previous.abs());
        assert!(output.ledger.relative_energy_drift < 1e-9, "sample {sample}: {}", output.ledger.relative_energy_drift);
    }
    assert!(peak > 0.0 && work > 0.0, "the drive must move the hole");
    let report = load.terminal_channels()[0].orifice_report().unwrap();
    assert_eq!(report.steps, 480);
    assert_eq!(report.vena_contracta, C_VC);
    assert_eq!(report.area_m2, area);
    assert_eq!(report.last_flow_m3_s, previous);
    near(report.dissipated_j, dissipated_from_frames, dissipated_from_frames);
    assert!(report.max_lag_drop_pa.is_finite() && report.lag_work_j >= 0.0 && report.law_work_j > 0.0);
    // The terminal's cumulative dissipation (and so the load's) holds exactly the orifice's work.
    near(load.terminal_channels()[0].cumulative_dissipated_j(), report.dissipated_j, report.dissipated_j);
}

#[test]
fn a_steady_flow_has_no_lag_and_the_drop_is_the_law() {
    let mut load = load(DT).with_terminal_channels(vec![orifice(SIDE, C_VC)], 0).unwrap();
    let area = load.terminal_area_m2(SIDE).unwrap().unwrap();
    // A constant inlet flow settles the hole flow; then |Q| = |Q_prev| and K |Q_prev| Q is K |Q| Q.
    for _ in 0..4000 { load.preview_step(AMPLITUDE).unwrap().accept(); }
    let before = load.terminal_channels()[0].orifice_report().unwrap();
    load.preview_step(AMPLITUDE).unwrap().accept();
    let frame = load.last_terminal_channel_frames()[0];
    let rho = load.terminal_channels()[0].gas().density;
    let law = rho / (2.0 * C_VC * C_VC * area * area) * frame.outward_flow_m3_s.abs() * frame.outward_flow_m3_s;
    assert!(law.abs() > 1e-6, "the settled flow must matter: {law}");
    assert!((frame.constitutive_pressure_drop_pa - law).abs() <= 1e-4 * law.abs(),
        "solved {} against law {law}", frame.constitutive_pressure_drop_pa);
    let after = load.terminal_channels()[0].orifice_report().unwrap();
    assert!(after.lag_work_j - before.lag_work_j <= 1e-9 * (after.law_work_j - before.law_work_j));
}

/// The lag shrinks with the step: the same 1 kHz drive over the same physical time at `dt` and `dt / 2`.
#[test]
fn halving_the_step_cuts_the_lags_share_of_the_laws_dissipation() {
    let share = |dt: f64| {
        let mut load = load(dt).with_terminal_channels(vec![orifice(SIDE, C_VC)], 0).unwrap();
        let samples = (0.004 / dt).round() as usize;
        for sample in 1..=samples { load.preview_step(drive(sample, dt)).unwrap().accept(); }
        let report = load.terminal_channels()[0].orifice_report().unwrap();
        assert!(report.law_work_j > 0.0);
        report.lag_work_j / report.law_work_j
    };
    let (coarse, fine) = (share(DT), share(DT / 2.0));
    assert!(coarse > 0.0 && fine > 0.0 && fine < 0.7 * coarse, "coarse {coarse}, fine {fine}");
}

#[test]
fn the_orifice_refuses_what_it_cannot_be() {
    for bad in [0.0, -0.61, f64::NAN, f64::INFINITY, 2.5] {
        assert!(orifice(SIDE, bad).validate().is_err(), "{bad}");
    }
    for good in [0.05, 0.61, 1.0956, 2.0] { orifice(SIDE, good).validate().unwrap(); }
    // The axial outlet is not a tone hole: no chimney, no orifice.
    assert!(load(DT).with_terminal_channels(vec![orifice(OUTLET, 0.61)], 0).is_err());
    // An undeclared port, and a second channel on the same port.
    assert!(load(DT).with_terminal_channels(vec![orifice(PortId(99), 0.61)], 0).is_err());
    assert!(load(DT).with_terminal_channels(vec![orifice(SIDE, 0.61), orifice(SIDE, 0.61)], 0).is_err());
    // A warm load takes no channels.
    let mut warm = load(DT);
    warm.preview_step(1e-9).unwrap().accept();
    assert!(warm.with_terminal_channels(vec![orifice(SIDE, 0.61)], 0).is_err());
}

#[test]
fn a_dropped_trial_leaves_the_orifice_and_its_group_where_they_were() {
    let mut load = load(DT).with_terminal_channels(vec![orifice(SIDE, C_VC)], 0).unwrap();
    for sample in 1..=40 { load.preview_step(drive(sample, DT)).unwrap().accept(); }
    let clean = load.clone();
    drop(load.preview_step(drive(41, DT)).unwrap());
    assert!(load.preview_step(f64::NAN).is_err());
    assert_eq!(load.terminal_channels()[0].orifice_report(), clean.terminal_channels()[0].orifice_report());
    assert_eq!(load.port_response().unwrap(), clean.port_response().unwrap());
    assert_eq!(load.preview_step(drive(41, DT)).unwrap().accept(),
        { let mut c = clean; c.preview_step(drive(41, DT)).unwrap().accept() });
}
