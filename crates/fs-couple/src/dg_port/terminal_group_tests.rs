use super::*;
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::MimoLoadLimits;
use fs_vfit::model::PoleTerm;

const DT: f64 = 1.0 / 96000.0;
const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const SIDE: PortId = PortId(11);

mod many_ports;

fn duct(walls: bool, side: bool, outlet: DgBoundary) -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let mut duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }] },
        &gas, DgMesh { elements: if side { 2 } else { 1 }, order: 2 },
        DgBoundary::Port(INLET), outlet,
    ).unwrap();
    if walls { duct = duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap(); }
    if side { duct = duct.with_side_ports(&[DgSidePort { id: SIDE, downstream_element: 1 }]).unwrap(); }
    duct
}

fn model() -> MimoModel {
    // The admitted real/pair fixture, scaled uniformly to acoustic impedance.
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

fn scalar_model() -> MimoModel {
    let model = model();
    MimoModel {
        ports: 1, poles: model.poles,
        residues: model.residues.into_iter().map(|matrix| vec![matrix[0]]).collect(),
        d: vec![model.d[0]], e: vec![0.0], reciprocal: true,
    }
}

fn prepared(model: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(model, dt, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / dt], sampled_passivity_points: 129,
        max_states: 16, max_memory_bytes: 1 << 20,
    }).unwrap()
}

fn loaded(walls: bool) -> DgPortLoad {
    DgPortLoad::new_with_terminal_ports(
        duct(walls, true, DgBoundary::Port(OUTLET)), DT, 1e-11, MAX_DG_PROPAGATOR_DOFS,
        vec![SIDE, OUTLET], prepared(&model(), DT),
    ).unwrap()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2e-10 * scale.max(1e-300),
        "actual={actual}, expected={expected}, scale={scale}");
}

#[test]
fn full_terminal_group_matches_analytic_two_by_two_and_original_state_owners() {
    for walls in [false, true] {
        let mut load = loaded(walls);
        assert_eq!(load.terminal_ports(), &[SIDE, OUTLET]);
        assert!(load.last_terminal_frames().is_empty());
        let mut rational = model().zoh(DT).unwrap();
        let propagator = load.propagator().clone();
        let mut state = load.state().clone();
        let mut inlet_work = 0.0;
        let mut terminal_work = 0.0;
        for sample in 0..96 {
            let flow = if sample < 37 { (f64::from(sample % 7) - 3.0) * 1e-9 } else { 0.0 };
            let previous = load.state().clone();
            let response = propagator.multiport_response(&state).unwrap();
            assert_eq!(response.ports, vec![INLET, OUTLET, SIDE]);
            let relation = load.port_response().unwrap();
            assert_eq!(load.state(), &previous);
            let terminal = load.terminal().unwrap();
            let mut history = [0.0; 2];
            terminal.mean_history_output_into(&mut history).unwrap();
            let direct = terminal.mean_direct_term();
            // Independently assemble and solve the 2x2 system in the supplied
            // load order [SIDE, OUTLET], unlike DG order [INLET, OUTLET, SIDE].
            let index = [2, 1];
            let mut matrix = [0.0; 4];
            let mut rhs = [0.0; 2];
            for i in 0..2 {
                rhs[i] = response.history_pressure_pa[index[i]] - history[i]
                    + response.impedance_pa_s_m3[3 * index[i]] * flow;
                for j in 0..2 {
                    matrix[2 * i + j] = response.impedance_pa_s_m3[3 * index[i] + index[j]] + direct[2 * i + j];
                }
            }
            let determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2];
            let analytic = [
                (rhs[0] * matrix[3] - matrix[1] * rhs[1]) / determinant,
                (matrix[0] * rhs[1] - rhs[0] * matrix[2]) / determinant,
            ];
            let flow_scale = analytic.iter().map(|value| value.abs()).sum::<f64>();
            let trial = load.preview_step(flow).unwrap();
            assert_eq!(trial.previous_state(), &previous);
            let frames = trial.terminal_frames();
            assert_eq!(frames.iter().map(|frame| frame.port).collect::<Vec<_>>(), vec![SIDE, OUTLET]);
            assert_eq!(trial.terminal_frame(), Some(&frames[1]));
            let outward = [frames[0].outward_flow_m3_s, frames[1].outward_flow_m3_s];
            for i in 0..2 { near(outward[i], analytic[i], flow_scale); }
            let pressures = rational.step_mean(&outward);
            let expected = propagator.step(&mut state, &[
                PortInput { id: INLET, effort: flow },
                PortInput { id: SIDE, effort: -outward[0] },
                PortInput { id: OUTLET, effort: -outward[1] },
            ]).unwrap();
            propagator.sync_nodal(&mut state).unwrap();
            assert_eq!(trial.state(), &state);
            assert_eq!(trial.output(), &expected);
            for (i, frame) in frames.iter().enumerate() {
                let scale = frame.mean_pressure_pa.abs() + pressures[i].abs();
                near(frame.mean_pressure_pa, pressures[i], scale);
                near(frame.port_work_j, outward[i] * expected.port_outputs[index[i]].delta_x,
                    frame.port_work_j.abs());
                near(frame.load_work_j, pressures[i] * (outward[i] * DT), frame.load_work_j.abs());
                near(frame.port_work_residual_j, 0.0,
                    frame.port_work_j.abs() + frame.load_work_j.abs());
            }
            near(expected.port_outputs[0].conjugate_rate,
                relation.history_pressure_pa + relation.impedance_pa_s_m3 * flow,
                relation.history_pressure_pa.abs() + (relation.impedance_pa_s_m3 * flow).abs());
            let total = sum_terminal_work(frames);
            let work_scale = trial.inlet_work_j().abs() + frames.iter().map(|frame| frame.port_work_j.abs()).sum::<f64>();
            near(expected.ledger.port_work_step_j, trial.inlet_work_j() - total, work_scale);
            inlet_work += trial.inlet_work_j();
            terminal_work += total;
            let frames = frames.to_vec();
            assert_eq!(trial.accept(), expected);
            assert_eq!(load.last_terminal_frames(), frames);
            assert_eq!(load.last_terminal_frame(), Some(&frames[1]));
            assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
            near(state.stored_energy_j() + state.cumulative_dissipated_j() + terminal_work,
                inlet_work, state.stored_energy_j() + inlet_work.abs() + terminal_work.abs());
        }
        assert!(rational.pole_states().iter().any(|state| *state != C64::ZERO));
        if walls { assert!(load.state().wall_stored_energy_j() > 0.0); }
    }
}

#[test]
fn abandoned_and_failed_group_trials_preserve_every_accepted_history_and_frame() {
    let mut load = loaded(true);
    for flow in [1e-8, -2e-8, 3e-8] { load.preview_step(flow).unwrap().accept(); }
    let mut clean = load.clone();
    let state = load.state().clone();
    let poles = load.terminal().unwrap().pole_states().to_vec();
    let frames = load.last_terminal_frames().to_vec();
    let response = load.port_response().unwrap();
    {
        let abandoned = load.preview_step(7e-9).unwrap();
        assert_ne!(abandoned.state(), &state);
        assert_eq!(abandoned.terminal_frames().len(), 2);
    }
    for bad_flow in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(load.preview_step(bad_flow).is_err());
        assert_eq!(load.state(), &state);
        assert_eq!(load.terminal().unwrap().pole_states(), poles);
        assert_eq!(load.last_terminal_frames(), frames);
        assert_eq!(load.port_response().unwrap(), response);
    }
    assert_eq!(load.preview_step(4e-9).unwrap().accept(), clean.preview_step(4e-9).unwrap().accept());
    assert_eq!(load.state(), clean.state());
    assert_eq!(load.terminal().unwrap().pole_states(), clean.terminal().unwrap().pole_states());
    assert_eq!(load.last_terminal_frames(), clean.last_terminal_frames());
}

#[test]
fn group_admission_owns_all_channels_and_side_only_load_has_no_axial_frame() {
    let source = duct(false, true, DgBoundary::Port(OUTLET));
    for ports in [vec![], vec![SIDE], vec![SIDE, SIDE], vec![INLET, OUTLET], vec![SIDE, PortId(999)]] {
        assert!(DgPortLoad::new_with_terminal_ports(source.clone(), DT, 1e-11,
            MAX_DG_PROPAGATOR_DOFS, ports, prepared(&model(), DT)).is_err());
    }
    assert!(DgPortLoad::new_with_terminal_ports(source.clone(), DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, vec![SIDE, OUTLET], prepared(&model(), 2.0 * DT)).is_err());
    let mut warm = prepared(&model(), DT);
    warm.step_mean_into(&[1e-8, 0.0], &mut [0.0; 2]).unwrap();
    assert!(DgPortLoad::new_with_terminal_ports(source.clone(), DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, vec![SIDE, OUTLET], warm).is_err());
    assert!(DgPortLoad::new_with_terminal(source, DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, prepared(&scalar_model(), DT)).is_err());
    let side_only = duct(false, true, DgBoundary::IdealOpen);
    assert!(DgPortLoad::new(side_only.clone(), DT, 1e-11, MAX_DG_PROPAGATOR_DOFS).is_err());
    let mut load = DgPortLoad::new_with_terminal_ports(side_only, DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, vec![SIDE], prepared(&scalar_model(), DT)).unwrap();
    load.preview_step(1e-8).unwrap().accept();
    assert_eq!(load.last_terminal_frames().len(), 1);
    assert_eq!(load.last_terminal_frames()[0].port, SIDE);
    assert!(load.last_terminal_frame().is_none());
}

#[test]
fn scalar_outlet_preserves_original_divide_schur_and_propagation_bits() {
    for walls in [false, true] {
        let source = duct(walls, false, DgBoundary::Port(OUTLET));
        let terminal = prepared(&scalar_model(), DT);
        let mut load = DgPortLoad::new_with_terminal(source.clone(), DT, 1e-11,
            MAX_DG_PROPAGATOR_DOFS, terminal.clone()).unwrap();
        let mut explicit = DgPortLoad::new_with_terminal_ports(source, DT, 1e-11,
            MAX_DG_PROPAGATOR_DOFS, vec![OUTLET], terminal.clone()).unwrap();
        let propagator = load.propagator().clone();
        let mut state = load.state().clone();
        let mut rational = terminal;
        for flow in [1e-8, -2e-8, 0.0, 3e-8, 0.0, 0.0] {
            let response = propagator.multiport_response(&state).unwrap();
            let h = &response.history_pressure_pa;
            let z = &response.impedance_pa_s_m3;
            let mut history = [0.0];
            rational.mean_history_output_into(&mut history).unwrap();
            // Original scalar algorithm, kept independent of the new cache.
            let denominator = z[3] + rational.mean_direct_term()[0];
            let outward_history = (h[1] - history[0]) / denominator;
            let outward_per_inlet = z[2] / denominator;
            let expected_relation = DgPortResponse {
                history_pressure_pa: h[0] - z[1] * outward_history,
                impedance_pa_s_m3: z[0] - z[1] * outward_per_inlet,
            };
            let actual_relation = load.port_response().unwrap();
            assert_eq!(actual_relation.history_pressure_pa.to_bits(), expected_relation.history_pressure_pa.to_bits());
            assert_eq!(actual_relation.impedance_pa_s_m3.to_bits(), expected_relation.impedance_pa_s_m3.to_bits());
            let outward = outward_history + outward_per_inlet * flow;
            rational.step_mean_into(&[outward], &mut [0.0]).unwrap();
            let expected = propagator.step(&mut state, &[
                PortInput { id: INLET, effort: flow }, PortInput { id: OUTLET, effort: -outward },
            ]).unwrap();
            propagator.sync_nodal(&mut state).unwrap();
            let trial = load.preview_step(flow).unwrap();
            assert_eq!(trial.terminal_frames()[0].outward_flow_m3_s.to_bits(), outward.to_bits());
            assert_eq!(trial.state(), &state);
            assert_eq!(trial.accept(), expected);
            assert_eq!(explicit.preview_step(flow).unwrap().accept(), expected);
            assert_eq!(load.last_terminal_frames(), explicit.last_terminal_frames());
            assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
        }
    }
}
