use super::*;
use fs_duct::dg::{DgMesh, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};
use fs_vfit::model::PoleTerm;

const DT: f64 = 1.0 / 96_000.0;
const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);

fn duct(lossy: bool, inlet: DgBoundary, outlet: DgBoundary) -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let source = Duct {
        segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }],
    };
    let source = DgDuct::new(
        &source, &gas, DgMesh { elements: 1, order: 2 }, inlet, outlet,
    ).unwrap();
    if lossy { source.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap() } else { source }
}

fn resistance(ports: usize, value: f64) -> MimoModel {
    let mut d = vec![0.0; ports * ports];
    for i in 0..ports { d[i * ports + i] = value; }
    MimoModel {
        ports, poles: vec![], residues: vec![], d,
        e: vec![0.0; ports * ports], reciprocal: true,
    }
}
fn reactive_model() -> MimoModel {
    let mut model = resistance(1, 2e7);
    model.poles = vec![
        PoleTerm::Real { pole: -400.0, residue: 0.0 },
        PoleTerm::Pair { pole: C64::new(-900.0, 2000.0), residue: C64::ZERO },
    ];
    model.residues = vec![
        vec![C64::from_re(4e10)], vec![C64::new(4e9, 1e9)],
    ];
    model
}
fn prepared(model: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(model, dt, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / dt],
        sampled_passivity_points: 129, max_states: 16, max_memory_bytes: 1 << 20,
    }).unwrap()
}
fn loaded(lossy: bool, model: &MimoModel) -> DgPortLoad {
    DgPortLoad::new_with_terminal(
        duct(lossy, DgBoundary::Port(INLET), DgBoundary::Port(OUTLET)),
        DT, 1e-11, MAX_DG_PROPAGATOR_DOFS, prepared(model, DT),
    ).unwrap()
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2e-10 * scale.max(1e-300),
        "actual={actual}, expected={expected}, scale={scale}");
}
fn frame_bits(frame: Option<&DgTerminalFrame>) -> Option<(PortId, [u64; 5])> {
    frame.map(|frame| (frame.port, [
        frame.outward_flow_m3_s.to_bits(), frame.mean_pressure_pa.to_bits(),
        frame.port_work_j.to_bits(), frame.load_work_j.to_bits(),
        frame.port_work_residual_j.to_bits(),
    ]))
}

#[test]
fn loaded_dg_retains_original_duct_and_real_and_pair_pole_histories() {
    for lossy in [false, true] {
        let model = reactive_model();
        let mut load = loaded(lossy, &model);
        let mut rational = model.zoh(DT).unwrap();
        let propagator = load.propagator().clone();
        let mut reference = load.state().clone();
        let mut total_inlet = 0.0;
        let mut total_terminal = 0.0;
        for sample in 0..96 {
            let flow = if sample < 37 { (f64::from(sample % 7) - 3.0) * 1e-9 } else { 0.0 };
            let previous = load.state().clone();
            let relation = load.port_response().unwrap();
            let mut terminal_history = [0.0];
            let terminal = load.terminal().unwrap();
            terminal.mean_history_output_into(&mut terminal_history).unwrap();
            let terminal_direct = terminal.mean_direct_term()[0];
            let trial = load.preview_step(flow).unwrap();
            assert_eq!(trial.previous_state(), &previous);
            let frame = trial.terminal_frame().unwrap();
            assert_eq!(frame.port, OUTLET);
            let outward = frame.outward_flow_m3_s;
            let expected_pressure = rational.step_mean(&[outward])[0];
            near(frame.mean_pressure_pa, expected_pressure,
                terminal_history[0].abs() + (terminal_direct * outward).abs());
            let expected = propagator.step(&mut reference, &[
                PortInput { id: INLET, effort: flow },
                PortInput { id: OUTLET, effort: -outward },
            ]).unwrap();
            propagator.sync_nodal(&mut reference).unwrap();
            assert_eq!(trial.state(), &reference);
            assert_eq!(trial.output(), &expected);
            let inlet_pressure = expected.port_outputs[0].conjugate_rate;
            let affine_pressure = relation.history_pressure_pa + relation.impedance_pa_s_m3 * flow;
            near(inlet_pressure, affine_pressure,
                relation.history_pressure_pa.abs() + (relation.impedance_pa_s_m3 * flow).abs());
            let inlet_work = trial.inlet_work_j();
            near(inlet_work, DT * inlet_pressure * flow, inlet_work.abs());
            // Both reported works are positive OUT of the duct, INTO the load.
            near(frame.port_work_j, DT * outward * expected.port_outputs[1].conjugate_rate,
                frame.port_work_j.abs());
            near(frame.load_work_j, DT * outward * expected_pressure, frame.load_work_j.abs());
            near(frame.port_work_j - frame.load_work_j, 0.0,
                frame.port_work_j.abs() + frame.load_work_j.abs());
            near(frame.port_work_residual_j, frame.port_work_j - frame.load_work_j,
                frame.port_work_j.abs() + frame.load_work_j.abs());
            near(expected.ledger.port_work_step_j, inlet_work - frame.port_work_j,
                inlet_work.abs() + frame.port_work_j.abs());
            total_inlet += inlet_work;
            total_terminal += frame.load_work_j;
            let accepted_frame = frame_bits(trial.terminal_frame());
            assert_eq!(trial.accept(), expected);
            assert_eq!(load.state(), &reference);
            assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
            assert_eq!(frame_bits(load.last_terminal_frame()), accepted_frame);
            near(expected.ledger.stored_energy_j
                    + expected.ledger.cumulative_dissipated_j + total_terminal,
                total_inlet,
                total_inlet.abs() + total_terminal.abs() + expected.ledger.stored_energy_j);
        }
        assert!(rational.pole_states().iter().any(|state| *state != C64::ZERO));
    }
}

#[test]
fn abandoned_and_failed_terminal_trials_preserve_both_accepted_histories() {
    let mut load = loaded(true, &reactive_model());
    for _ in 0..12 { load.preview_step(2e-9).unwrap().accept(); }
    let previous = load.state().clone();
    let poles = load.terminal().unwrap().pole_states().to_vec();
    let response = load.port_response().unwrap();
    let frame = frame_bits(load.last_terminal_frame());
    let mut clean = load.clone();
    {
        let trial = load.preview_step(-1e-9).unwrap();
        assert_eq!(trial.previous_state(), &previous);
        assert_eq!(trial.state().step_count(), previous.step_count() + 1);
        assert!(trial.terminal_frame().is_some());
    }
    assert_eq!(load.state(), &previous);
    assert_eq!(load.terminal().unwrap().pole_states(), poles.as_slice());
    assert_eq!(load.port_response().unwrap(), response);
    assert_eq!(frame_bits(load.last_terminal_frame()), frame);
    for flow in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX] {
        assert!(load.preview_step(flow).is_err());
        assert_eq!(load.state(), &previous);
        assert_eq!(load.terminal().unwrap().pole_states(), poles.as_slice());
        assert_eq!(load.port_response().unwrap(), response);
        assert_eq!(frame_bits(load.last_terminal_frame()), frame);
    }
    let expected = clean.preview_step(-2e-9).unwrap().accept();
    assert_eq!(load.preview_step(-2e-9).unwrap().accept(), expected);
    assert_eq!(load.state(), clean.state());
    assert_eq!(load.terminal().unwrap().pole_states(), clean.terminal().unwrap().pole_states());
    assert_eq!(frame_bits(load.last_terminal_frame()), frame_bits(clean.last_terminal_frame()));
}

#[test]
fn terminal_admission_refuses_clock_port_count_topology_and_nonquiescent_state() {
    let model = reactive_model();
    let source = duct(false, DgBoundary::Port(INLET), DgBoundary::Port(OUTLET));
    for terminal in [prepared(&model, 2.0 * DT), prepared(&resistance(2, 2e7), DT)] {
        assert!(DgPortLoad::new_with_terminal(
            source.clone(), DT, 1e-11, MAX_DG_PROPAGATOR_DOFS, terminal,
        ).is_err());
    }
    for (inlet, outlet) in [
        (DgBoundary::Closed, DgBoundary::Port(OUTLET)),
        (DgBoundary::Port(INLET), DgBoundary::IdealOpen),
    ] {
        assert!(DgPortLoad::new_with_terminal(
            duct(false, inlet, outlet), DT, 1e-11, MAX_DG_PROPAGATOR_DOFS,
            prepared(&model, DT),
        ).is_err());
    }
    let mut warm = prepared(&model, DT);
    warm.step_mean_into(&[1e-9], &mut [0.0]).unwrap();
    assert!(DgPortLoad::new_with_terminal(
        source, DT, 1e-11, MAX_DG_PROPAGATOR_DOFS, warm,
    ).is_err());
}

#[test]
fn reactive_terminal_returns_signed_work_without_becoming_duct_dissipation() {
    let mut model = resistance(1, 1e6);
    model.poles = vec![PoleTerm::Real { pole: -400.0, residue: 0.0 }];
    model.residues = vec![vec![C64::from_re(4e10)]];
    let mut load = loaded(false, &model);
    for _ in 0..64 { load.preview_step(1e-8).unwrap().accept(); }
    let response = load.propagator().multiport_response(load.state()).unwrap();
    let terminal = load.terminal().unwrap();
    let mut history = [0.0];
    terminal.mean_history_output_into(&mut history).unwrap();
    let direct = terminal.mean_direct_term()[0];
    assert!(history[0] != 0.0 && direct > 0.0);
    // Oppose retained load pressure while leaving its mean at 3/4 of history.
    // This controlled energy return does not assume a waveform happens to reverse.
    let outward = -history[0] / (4.0 * direct);
    let z = &response.impedance_pa_s_m3;
    assert_ne!(z[2], 0.0);
    let inlet_flow = (outward * (z[3] + direct)
        - response.history_pressure_pa[1] + history[0]) / z[2];
    let trial = load.preview_step(inlet_flow).unwrap();
    let frame = trial.terminal_frame().unwrap();
    near(frame.outward_flow_m3_s, outward, outward.abs());
    assert!(frame.load_work_j < 0.0);
    assert!(frame.port_work_j < 0.0, "the terminal returns energy to the duct");
    near(trial.output().ledger.port_work_step_j, trial.inlet_work_j() - frame.port_work_j,
        trial.inlet_work_j().abs() + frame.port_work_j.abs());
    assert_eq!(trial.output().ledger.cumulative_dissipated_j, 0.0);
    trial.accept();
    assert!(load.last_terminal_frame().unwrap().load_work_j < 0.0);
}

#[test]
fn rational_memory_changes_the_bore_response_from_resistive_and_matched_outlets() {
    let mut reactive = loaded(false, &reactive_model());
    let mut resistive = loaded(false, &resistance(1, 2e7));
    let mut matched = DgPortLoad::new(
        duct(false, DgBoundary::Port(INLET), DgBoundary::Matched),
        DT, 1e-11, MAX_DG_PROPAGATOR_DOFS,
    ).unwrap();
    assert_ne!(reactive.port_response().unwrap().impedance_pa_s_m3.to_bits(),
        resistive.port_response().unwrap().impedance_pa_s_m3.to_bits());
    let mut resistive_difference = 0.0_f64;
    let mut matched_difference = 0.0_f64;
    let mut scale = 0.0_f64;
    for sample in 0..96 {
        let flow = if sample < 24 { 2e-9 } else { 0.0 };
        let a = reactive.preview_step(flow).unwrap().accept().port_outputs[0].conjugate_rate;
        let b = resistive.preview_step(flow).unwrap().accept().port_outputs[0].conjugate_rate;
        let c = matched.preview_step(flow).unwrap().accept().port_outputs[0].conjugate_rate;
        resistive_difference = resistive_difference.max((a - b).abs());
        matched_difference = matched_difference.max((a - c).abs());
        scale = scale.max(a.abs()).max(b.abs()).max(c.abs());
    }
    assert!(resistive_difference > 1e-6 * scale);
    assert!(matched_difference > 1e-6 * scale);
}
