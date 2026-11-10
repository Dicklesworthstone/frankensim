//! A player's airway as the valve's upstream owner: the account, the port law, the series source, the cold binding
//! and the passivity certificate. These fixtures do not qualify an instrument.
use super::*;
use crate::driving_point::LineLadder;
use crate::render::player_tract::{PlayerTract, TractPressure};
use crate::upstream_airway::{AirwayChain, AirwaySection};

const RATE: u32 = 96_000;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

/// One uniform tube, valve side first; its far end is the airway's pressure source.
fn chain(length_m: f64, area_m2: f64) -> AirwayChain {
    AirwayChain::new(
        vec![AirwaySection::from_area(length_m, area_m2, 1, None).unwrap()],
        fs_duct::LossModel::Bessel,
    ).unwrap()
}

fn tract(gas: &GasState, rate: u32, pressure: TractPressure) -> PlayerTract {
    PlayerTract::realize(&chain(0.17, 3.0e-4), gas, pressure, rate, &LineLadder::default()).unwrap()
}

fn airway_system() -> LoadedApertureDg {
    let down = load(0.05, 101325.0, 1.0 / f64::from(RATE));
    let owner: UpstreamOwner = tract(&air(), RATE, TractPressure::Lung).into();
    let (_, z) = owner.inlet_state(&down).unwrap();
    LoadedApertureDg::new_loaded_with_upstream(
        aperture(&down, z, BladePropagator::ExactZoh), owner, down,
    ).unwrap()
}

#[test]
fn an_airway_owner_balances_the_valve_account_and_keeps_its_own_ledger() {
    let mut model = airway_system();
    assert!(model.has_upstream() && model.upstream_load().is_none());
    assert_eq!(model.upstream_owner().unwrap().stored_energy_j(), 0.0);
    let dt = model.load().propagator().dt();
    let (mut absorbed, mut wave, mut reverse_flow) = (0.0, 0.0, false);
    for step in 0..128 {
        let drive = TubeDrive {
            upstream_pressure_pa: if step < 40 { 0.5 } else if step < 80 { -0.5 } else { 0.0 },
            body_flow_m3_s: if step < 16 { 2e-10 } else if step < 32 { -1e-10 } else { 0.0 },
        };
        let frame = model.step_with_forces(drive,
            [if step < 24 { 1e-6 } else { 0.0 }, 0.0]).unwrap();
        assert!(frame.upstream.is_none());
        let airway = frame.airway.expect("the airway's ledger frame");
        let through = frame.aperture.bore_flow_m3_s - drive.body_flow_m3_s;
        assert_eq!(airway.flow_m3_s, through);
        assert_eq!(airway.source_pa, drive.upstream_pressure_pa);
        assert_eq!(airway.series_pa, 0.0);
        assert_eq!(airway.face_pressure_pa,
            airway.open_circuit_pa - airway.impedance_pa_s_m3 * through);
        assert_eq!(airway.delivered_work_j,
            (airway.face_pressure_pa - airway.series_pa) * airway.flow_m3_s * dt);
        assert_eq!(frame.upstream_work_j, 0.0, "no series source: the gesture is the airway's far end");
        assert_eq!(frame.body_work_j, frame.aperture.bore_pressure_pa * (drive.body_flow_m3_s * dt));
        let scale = frame.stored_energy_j.abs() + frame.storage_change_j.abs()
            + frame.inferred_total_dissipation_j.abs() + airway.delivered_work_j.abs()
            + frame.body_work_j.abs() + frame.aperture.mechanical_work_j.abs()
            + sum_terminal_work(&frame.terminals).abs();
        assert!(frame.balance_residual_j().abs() <= 1e-8 * scale.max(1e-24),
            "step {step}: balance {:e} against {scale:e}", frame.balance_residual_j());
        assert!(frame.port_work_residual_j.abs() <= 1e-10 * scale.max(1e-24));
        absorbed += airway.absorbed_work_j;
        assert_eq!(airway.cumulative_absorbed_work_j, absorbed);
        assert!(absorbed >= -1e-6 * airway.wave_energy_j, "step {step}");
        wave = airway.wave_energy_j;
        reverse_flow |= through < 0.0;
    }
    assert!(wave > 0.0 && reverse_flow);
    assert_eq!(model.upstream_airway().unwrap().cumulative_absorbed_work_j(), absorbed);
    // The original DG tract's frame carries no airway term.
    let plain = system(0.035).step(TubeDrive { upstream_pressure_pa: 0.5, body_flow_m3_s: 0.0 }).unwrap();
    assert!(plain.airway.is_none() && plain.upstream.is_some());
}

#[test]
fn the_airway_port_law_reproduces_the_chain_input_impedance_under_a_planted_flow() {
    // A 0.17 m tube of 3e-4 m^2 with its far end the pressure source at rest; a flow of 1 uL/s amplitude planted
    // into the valve side at 150 Hz (a third of the first resonance, c / 4 L = 504 Hz). The port law
    // p = p_oc - Z U gives, with U the flow OUT of the airway, Z_in = -P / U for the phasors of the steady part.
    let gas = air();
    let chain = chain(0.17, 3.0e-4);
    let mut tract = PlayerTract::realize(&chain, &gas, TractPressure::Lung, RATE, &LineLadder::default()).unwrap();
    let (rate, f, u0) = (f64::from(RATE), 150.0, 1e-6);
    let omega = core::f64::consts::TAU * f;
    let per_cycle = 640; // 96000 / 150
    let (warm, measure) = (80, 20);
    let (mut pa, mut pb) = (0.0, 0.0);
    for n in 0..(warm + measure) * per_cycle {
        let t = n as f64 / rate;
        let flow = u0 * (omega * t).sin();
        let port = tract.dg_port(0.0, 0.0).unwrap();
        let frame = tract.dg_commit(port, flow).unwrap();
        if n >= warm * per_cycle {
            pa += frame.face_pressure_pa * (omega * t).sin();
            pb += frame.face_pressure_pa * (omega * t).cos();
        }
    }
    let norm = 2.0 / (measure * per_cycle) as f64;
    // The crate's clock is exp(-iwt) (upstream_airway.rs module doc): x = A sin + B cos has the phasor B + iA, so
    // P = B_p + i A_p, U = i u0 and Z = -P / U = (-A_p + i B_p) / u0.
    let measured = (-pa * norm / u0, pb * norm / u0);
    let analytic = chain.input_impedance(&gas, omega).unwrap();
    let error = ((measured.0 - analytic.re).powi(2) + (measured.1 - analytic.im).powi(2)).sqrt();
    let size = analytic.re.hypot(analytic.im);
    assert!(size > 1e5, "the planted probe must see a real impedance: {size:e}");
    assert!(error <= 0.03 * size, "measured {measured:?} against analytic ({:e}, {:e})", analytic.re, analytic.im);
}

#[test]
fn the_series_source_at_the_valve_is_set_per_step_and_a_dg_tract_refuses_it() {
    let mut model = airway_system();
    let dt = model.load().propagator().dt();
    assert!(model.set_airway_series_pa(f64::NAN).is_err());
    assert!(model.set_airway_series_pa(0.5).is_ok());
    let drive = TubeDrive { upstream_pressure_pa: 1.0, body_flow_m3_s: 0.0 };
    let frame = model.step(drive).unwrap();
    let airway = frame.airway.unwrap();
    let through = frame.aperture.bore_flow_m3_s;
    assert_eq!(airway.series_pa, 0.5);
    assert_eq!(frame.upstream_work_j, 0.5 * (through * dt));
    assert_eq!(airway.delivered_work_j, (airway.face_pressure_pa - 0.5) * through * dt);
    // The step that accepted it took it: the next one has none.
    assert_eq!(model.step(drive).unwrap().airway.unwrap().series_pa, 0.0);
    // A DG tract takes the noise with the gesture pressure at its inlet.
    assert!(system(0.035).set_airway_series_pa(0.5).is_err());
}

#[test]
fn the_airway_binds_cold_in_the_bores_gas_on_its_clock_and_against_the_summed_impedance() {
    let down = load(0.05, 101325.0, 1.0 / f64::from(RATE));
    let owner: UpstreamOwner = tract(&air(), RATE, TractPressure::Lung).into();
    let (gas, z) = owner.inlet_state(&down).unwrap();
    assert_eq!(&gas, down.duct().inlet_gas());
    assert_eq!(z, owner.impedance_pa_s_m3() + down.impedance_pa_s_m3());
    // A valve declared against the bore alone, or in the wrong mode, is refused.
    assert!(LoadedApertureDg::new_loaded_with_upstream(
        aperture(&down, down.impedance_pa_s_m3(), BladePropagator::ExactZoh),
        owner.clone(), down.clone(),
    ).is_err());
    assert!(LoadedApertureDg::new_loaded_with_upstream(
        aperture(&down, z, BladePropagator::Midpoint), owner.clone(), down.clone(),
    ).is_err());
    // Another gas, another clock, an airway that has already advanced.
    let other = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 81060.0).unwrap();
    assert!(UpstreamOwner::from(tract(&other, RATE, TractPressure::Lung)).inlet_state(&down).is_err());
    assert!(UpstreamOwner::from(tract(&air(), 48_000, TractPressure::Lung)).inlet_state(&down).is_err());
    let mut advanced = tract(&air(), RATE, TractPressure::Lung);
    assert!(advanced.is_cold());
    let port = advanced.dg_port(0.0, 0.0).unwrap();
    advanced.dg_commit(port, 1e-7).unwrap();
    assert!(!advanced.is_cold());
    assert!(UpstreamOwner::from(advanced).inlet_state(&down).is_err());
}

#[test]
fn a_held_source_with_no_flow_adds_nothing_to_the_source_free_ledger() {
    // The far-end source reaches the valve face, but with no flow the source-free airway takes and holds nothing:
    // its waves are not the source's. (The draft read the reflectance history, which carries the source's own
    // reflected waves, and charged the gesture's power to the airway.)
    let mut tract = tract(&air(), RATE, TractPressure::Lung);
    let mut peak_open_circuit = 0.0_f64;
    for _ in 0..600 {
        let port = tract.dg_port(1000.0, 0.0).unwrap();
        let frame = tract.dg_commit(port, 0.0).unwrap();
        assert_eq!(frame.absorbed_work_j, 0.0);
        assert_eq!(frame.wave_energy_j, 0.0);
        assert_eq!(frame.cumulative_absorbed_work_j, 0.0);
        peak_open_circuit = peak_open_circuit.max(frame.open_circuit_pa.abs());
    }
    // The face rings against the closed valve between about zero and twice the source, so the peak is read, not one sample.
    assert!(peak_open_circuit > 100.0, "the source must reach the face: {peak_open_circuit}");
}

#[test]
fn an_active_regulation_refuses_by_name_and_leaves_the_airway_as_it_was() {
    let base = tract(&air(), RATE, TractPressure::Lung);
    let zs = base.series_impedance_pa_s_m3();
    // A mouth-pressure regulation whose DC resistance is twenty times the airway's series impedance acts as a source
    // against a held flow; the admission (`admitted_margin`) would have refused it, so it is planted.
    let mut active = base.clone().with_unadmitted_regulation(1e-3, 20.0 * zs);
    let mut passive = base;
    let flow = 1e-6;
    let mut refused_at = None;
    for n in 0..4000 {
        let before = (active.cumulative_absorbed_work_j(), active.flow_m3_s(), active.mouth_pressure_pa());
        let port = active.dg_port(0.0, 0.0).unwrap();
        match active.dg_commit(port, flow) {
            Ok(_) => {}
            Err(error) => {
                let text = error.to_string();
                assert!(text.contains("player tract refused") && text.contains("acting as a source"), "{text}");
                // A refusal moved nothing: ledger, flow memory and face pressure are as before the sample.
                assert_eq!((active.cumulative_absorbed_work_j(), active.flow_m3_s(),
                    active.mouth_pressure_pa()), before);
                refused_at = Some(n);
                break;
            }
        }
        // The lung tract (provably passive on its grid) takes the same flow without a refusal.
        let port = passive.dg_port(0.0, 0.0).unwrap();
        passive.dg_commit(port, flow).unwrap();
    }
    assert!(refused_at.is_some_and(|n| n > 0), "{refused_at:?}");
    assert!(passive.cumulative_absorbed_work_j() > 0.0);
}
