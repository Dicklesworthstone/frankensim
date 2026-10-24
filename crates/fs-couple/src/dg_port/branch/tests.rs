use super::*;
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_la::prepared_lu::LuWorkspace;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::MimoLoadLimits;
use fs_vfit::model::PoleTerm;

const DT: f64 = 1.0 / 96000.0;
const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const FIRST: PortId = PortId(11);
const SECOND: PortId = PortId(19);

fn gas(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}
fn main_duct(walls: bool, count: usize) -> DgDuct {
    let mut duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.05 }] },
        &gas(101325.0), DgMesh { elements: count + 1, order: 2 },
        DgBoundary::Port(INLET), DgBoundary::Port(OUTLET),
    ).unwrap().with_side_ports(&[FIRST, SECOND][..count].iter().enumerate()
        .map(|(i, &id)| DgSidePort { id, downstream_element: i + 1 }).collect::<Vec<_>>()).unwrap();
    if walls { duct = duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap(); }
    duct
}
fn finite_duct(walls: bool, length: f64, gas: &GasState) -> DgDuct {
    let mut duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length }] },
        gas, DgMesh { elements: 1, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Port(PortId(1)),
    ).unwrap();
    if walls { duct = duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap(); }
    duct
}
fn model(count: usize) -> MimoModel {
    if count == 2 {
        return MimoModel {
            ports: 2,
            poles: vec![PoleTerm::Real { pole: -400.0, residue: 0.0 },
                PoleTerm::Pair { pole: C64::new(-900.0, 2000.0), residue: C64::ZERO }],
            residues: vec![
                vec![C64::from_re(2e9), C64::from_re(5e8), C64::from_re(5e8), C64::from_re(1.8e9)],
                vec![C64::new(2e8, 1e8), C64::new(2e7, 1e7), C64::new(2e7, 1e7), C64::new(1.5e8, 5e7)],
            ],
            d: vec![2e7, 5e6, 5e6, 3e7], e: vec![0.0; 4], reciprocal: true,
        };
    }
    let d = vec![2e7, 5e6, 2e6, 5e6, 3e7, 4e6, 2e6, 4e6, 2.5e7];
    MimoModel {
        ports: 3, poles: vec![PoleTerm::Real { pole: -400.0, residue: 0.0 }],
        residues: vec![d.iter().map(|value| C64::from_re(100.0 * value)).collect()],
        d, e: vec![0.0; 9], reciprocal: true,
    }
}
fn prepared(model: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(model, dt, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / dt], sampled_passivity_points: 129,
        max_states: 32, max_memory_bytes: 1 << 20,
    }).unwrap()
}
fn loaded(walls: bool, count: usize) -> DgPortLoad {
    let ports = if count == 1 { vec![FIRST, OUTLET] } else { vec![SECOND, OUTLET, FIRST] };
    let branches = [FIRST, SECOND][..count].iter().enumerate().map(|(i, &port)|
        DgBranch { port, duct: finite_duct(walls, 0.012 + 0.006 * i as f64, &gas(101325.0)) }
    ).collect();
    DgPortLoad::new_with_terminal_branches(main_duct(walls, count), DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, ports, prepared(&model(count + 1), DT), branches).unwrap()
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2e-9 * scale.max(1e-300),
        "actual={actual}, expected={expected}, scale={scale}");
}
fn index(response: &DgMultiportResponse, port: PortId) -> usize {
    response.ports.iter().position(|id| *id == port).unwrap()
}

// Independent UNCONDENSED system: unknowns are all junction flows followed by
// all outward load flows. Each row is a physical pressure-continuity equation;
// no BranchCoupling response, transfer, inverse or effective inlet law is used.
fn uncondensed(
    main: &DgMultiportResponse, branches: &[DgMultiportResponse],
    branch_ports: &[PortId], terminal_ports: &[PortId], load_history: &[f64],
    load_direct: &[f64], inlet_flow: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let b = branches.len();
    let t = terminal_ports.len();
    let width = b + t;
    let n = main.ports.len();
    let inlet = index(main, INLET);
    let mut matrix = vec![0.0; width * width];
    let mut rhs = vec![0.0; width];
    for i in 0..b {
        let row = index(main, branch_ports[i]);
        rhs[i] = main.history_pressure_pa[row] - branches[i].history_pressure_pa[0]
            + main.impedance_pa_s_m3[row * n + inlet] * inlet_flow;
        for j in 0..b {
            matrix[i * width + j] = main.impedance_pa_s_m3[row * n + index(main, branch_ports[j])];
            if i == j { matrix[i * width + j] += branches[i].impedance_pa_s_m3[0]; }
        }
        for (j, &port) in terminal_ports.iter().enumerate() {
            matrix[i * width + b + j] = match branch_ports.iter().position(|id| *id == port) {
                Some(k) if k == i => -branches[i].impedance_pa_s_m3[1],
                Some(_) => 0.0,
                None => main.impedance_pa_s_m3[row * n + index(main, port)],
            };
        }
    }
    for (i, &port) in terminal_ports.iter().enumerate() {
        let row = b + i;
        for j in 0..t { matrix[row * width + b + j] = load_direct[i * t + j]; }
        if let Some(j) = branch_ports.iter().position(|id| *id == port) {
            rhs[row] = branches[j].history_pressure_pa[1] - load_history[i];
            matrix[row * width + j] = -branches[j].impedance_pa_s_m3[2];
            matrix[row * width + b + i] += branches[j].impedance_pa_s_m3[3];
        } else {
            let local = index(main, port);
            rhs[row] = main.history_pressure_pa[local] - load_history[i]
                + main.impedance_pa_s_m3[local * n + inlet] * inlet_flow;
            for j in 0..b {
                matrix[row * width + j] =
                    main.impedance_pa_s_m3[local * n + index(main, branch_ports[j])];
            }
            for (j, &other) in terminal_ports.iter().enumerate() {
                if !branch_ports.contains(&other) {
                    matrix[row * width + b + j] +=
                        main.impedance_pa_s_m3[local * n + index(main, other)];
                }
            }
        }
    }
    let mut lu = LuWorkspace::new(width).unwrap();
    lu.factor(&matrix).unwrap();
    let mut solution = vec![0.0; width];
    lu.solve_factored_into(&rhs, &mut solution).unwrap();
    (solution, matrix, rhs)
}

#[test]
fn finite_junctions_match_uncondensed_system_and_every_original_owner_over_history() {
    for walls in [false, true] {
        for count in [1, 2] {
            let mut load = loaded(walls, count);
            assert!(load.last_branch_frames().is_empty());
            let ports = load.terminal_ports().to_vec();
            let branch_ports: Vec<_> = load.branches().iter().map(DgBranchLoad::port).collect();
            let main_propagator = load.propagator().clone();
            let mut main_state = load.state().clone();
            let propagators: Vec<_> = load.branches().iter().map(|b| b.propagator().clone()).collect();
            let mut states: Vec<_> = load.branches().iter().map(|b| b.state().clone()).collect();
            let mut rational = model(count + 1).zoh(DT).unwrap();
            let direct = load.terminal().unwrap().mean_direct_term().to_vec();
            let mut endpoint_difference = 0.0_f64;
            for sample in 0..96 {
                let flow = if sample < 24 { (f64::from(sample % 5) - 2.0) * 1e-9 } else { 0.0 };
                let main_response = main_propagator.multiport_response(&main_state).unwrap();
                let responses: Vec<_> = propagators.iter().zip(&states)
                    .map(|(p, s)| p.multiport_response(s).unwrap()).collect();
                let mut free = rational.clone();
                let history = free.step_mean(&vec![0.0; ports.len()]);
                let (analytic, matrix, rhs) = uncondensed(&main_response, &responses,
                    &branch_ports, &ports, &history, &direct, flow);
                let width = analytic.len();
                // Check the independent solve's original rows as well as its
                // agreement with the separately condensed runtime.
                for row in 0..width {
                    let value: f64 = (0..width).map(|j| matrix[row * width + j] * analytic[j]).sum();
                    let scale: f64 = rhs[row].abs() + (0..width)
                        .map(|j| (matrix[row * width + j] * analytic[j]).abs()).sum::<f64>();
                    near(value, rhs[row], scale);
                }
                let previous_energy = load.stored_energy_j();
                let previous_loss = load.cumulative_dissipated_j();
                let trial = load.preview_step(flow).unwrap();
                assert_eq!(trial.previous_stored_energy_j(), previous_energy);
                assert_eq!(trial.previous_cumulative_dissipated_j(), previous_loss);
                let junctions: Vec<_> = trial.branch_frames().iter()
                    .map(|f| f.junction_flow_m3_s).collect();
                let outward: Vec<_> = trial.terminal_frames().iter()
                    .map(|f| f.outward_flow_m3_s).collect();
                let flow_scale: f64 = analytic.iter().map(|v| v.abs()).sum();
                for (actual, expected) in junctions.iter().chain(&outward).zip(&analytic) {
                    near(*actual, *expected, flow_scale);
                }
                let mut main_inputs = vec![PortInput { id: INLET, effort: flow }];
                for (i, &port) in ports.iter().enumerate() {
                    let value = branch_ports.iter().position(|id| *id == port)
                        .map_or(outward[i], |j| junctions[j]);
                    main_inputs.push(PortInput { id: port, effort: -value });
                }
                let main_step = main_propagator.step(&mut main_state, &main_inputs).unwrap();
                main_propagator.sync_nodal(&mut main_state).unwrap();
                assert_eq!(trial.state(), &main_state);
                let pressure = rational.step_mean(&outward);
                let mut energy = main_step.ledger.stored_energy_j;
                let mut work = main_step.ledger.port_work_step_j;
                let mut loss = main_step.ledger.cumulative_dissipated_j;
                for i in 0..count {
                    let external = ports.iter().position(|id| *id == branch_ports[i]).unwrap();
                    let expected = propagators[i].step(&mut states[i], &[
                        PortInput { id: PortId(0), effort: junctions[i] },
                        PortInput { id: PortId(1), effort: -outward[external] },
                    ]).unwrap();
                    propagators[i].sync_nodal(&mut states[i]).unwrap();
                    let frame = &trial.branch_frames()[i];
                    assert_eq!(frame.duct, expected);
                    assert_eq!(frame.bore_port,
                        main_step.port_outputs[index(&main_response, branch_ports[i])]);
                    assert_eq!(trial.branch_states[i], states[i]);
                    near(frame.junction_work_residual_j, 0.0,
                        frame.bore_work_j.abs() + frame.branch_work_j.abs());
                    endpoint_difference += (junctions[i] - outward[external]).abs();
                    let outer = &trial.terminal_frames()[external];
                    assert_eq!(outer.mean_pressure_pa, expected.port_outputs[1].conjugate_rate);
                    energy += expected.ledger.stored_energy_j;
                    work += expected.ledger.port_work_step_j;
                    loss += expected.ledger.cumulative_dissipated_j;
                }
                assert_eq!(trial.output().ledger.stored_energy_j, energy);
                assert_eq!(trial.output().ledger.port_work_step_j, work);
                assert_eq!(trial.output().ledger.cumulative_dissipated_j, loss);
                for (i, frame) in trial.terminal_frames().iter().enumerate() {
                    near(frame.mean_pressure_pa, pressure[i],
                        frame.mean_pressure_pa.abs() + pressure[i].abs());
                }
                let terminal_work = sum_terminal_work(trial.terminal_frames());
                let scale = trial.inlet_work_j().abs() + trial.terminal_frames().iter()
                    .map(|f| f.port_work_j.abs()).sum::<f64>()
                    + trial.branch_frames().iter().map(|f| f.bore_work_j.abs()).sum::<f64>();
                near(work, trial.inlet_work_j() - terminal_work, scale);
                let frames = trial.branch_frames().to_vec();
                trial.accept();
                assert_eq!(load.last_branch_frames(), frames);
                assert_eq!(load.stored_energy_j(), energy);
                assert_eq!(load.cumulative_dissipated_j(), loss);
                assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
            }
            assert!(endpoint_difference > 0.0, "a finite compressible chimney must retain distinct endpoint flows");
            assert!(states.iter().all(|state| state.stored_energy_j() > 0.0));
            if walls { assert!(states.iter().all(|state| state.wall_stored_energy_j() > 0.0)); }
        }
    }
}

#[test]
fn dropped_or_failed_trials_preserve_bore_branch_wall_and_rational_histories() {
    let mut load = loaded(true, 2);
    for flow in [1e-8, -2e-8, 3e-8] { load.preview_step(flow).unwrap().accept(); }
    let mut clean = load.clone();
    let main = load.state().clone();
    let states: Vec<_> = load.branches().iter().map(|b| b.state().clone()).collect();
    let poles = load.terminal().unwrap().pole_states().to_vec();
    let frames = load.last_branch_frames().to_vec();
    let terminals = load.last_terminal_frames().to_vec();
    let relation = load.port_response().unwrap();
    { let trial = load.preview_step(7e-9).unwrap(); assert_ne!(trial.branch_states, states); }
    for flow in [f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(load.preview_step(flow).is_err());
        assert_eq!(load.state(), &main);
        assert_eq!(load.branches().iter().map(|b| b.state().clone()).collect::<Vec<_>>(), states);
        assert_eq!(load.terminal().unwrap().pole_states(), poles);
        assert_eq!(load.last_branch_frames(), frames);
        assert_eq!(load.last_terminal_frames(), terminals);
        assert_eq!(load.port_response().unwrap(), relation);
    }
    assert_eq!(load.preview_step(4e-9).unwrap().accept(), clean.preview_step(4e-9).unwrap().accept());
    assert_eq!(load.last_branch_frames(), clean.last_branch_frames());
    assert_eq!(load.last_terminal_frames(), clean.last_terminal_frames());
}

#[test]
fn branch_admission_checks_real_junction_material_clock_topology_and_combined_state_budget() {
    let main = main_duct(false, 1);
    let valid = DgBranch { port: FIRST, duct: finite_duct(false, 0.012, &gas(101325.0)) };
    let build = |main: DgDuct, branches: Vec<DgBranch>, dt: f64, budget: usize| {
        DgPortLoad::new_with_terminal_branches(main, dt, 1e-11, budget,
            vec![FIRST, OUTLET], prepared(&model(2), dt), branches)
    };
    let sum = main.total_dofs() + valid.duct.total_dofs();
    assert!(build(main.clone(), vec![valid.clone()], DT, sum).is_ok());
    assert!(build(main.clone(), vec![valid.clone()], DT, sum - 1).is_err());
    assert!(build(main.clone(), vec![valid.clone(), valid.clone()], DT, MAX_DG_PROPAGATOR_DOFS).is_err());
    for port in [INLET, OUTLET, PortId(999)] {
        assert!(build(main.clone(), vec![DgBranch { port, ..valid.clone() }],
            DT, MAX_DG_PROPAGATOR_DOFS).is_err());
    }
    let wrong_pressure = DgBranch { port: FIRST, duct: finite_duct(false, 0.012, &gas(90000.0)) };
    assert!(build(main.clone(), vec![wrong_pressure], DT, MAX_DG_PROPAGATOR_DOFS).is_err());
    let wall = DgBranch { port: FIRST, duct: finite_duct(true, 0.012, &gas(101325.0)) };
    assert!(build(main.clone(), vec![wall], 1.0 / 3000.0, MAX_DG_PROPAGATOR_DOFS).is_err());
    let bad_topology = DgBranch { port: FIRST, duct: main_duct(false, 1) };
    assert!(build(main, vec![bad_topology], DT, MAX_DG_PROPAGATOR_DOFS).is_err());
}

#[test]
fn an_empty_branch_list_preserves_existing_group_arithmetic_and_histories() {
    for walls in [false, true] {
        let duct = main_duct(walls, 1);
        let terminal = prepared(&model(2), DT);
        let mut original = DgPortLoad::new_with_terminal_ports(duct.clone(), DT, 1e-11,
            MAX_DG_PROPAGATOR_DOFS, vec![FIRST, OUTLET], terminal.clone()).unwrap();
        let mut empty = DgPortLoad::new_with_terminal_branches(duct, DT, 1e-11,
            MAX_DG_PROPAGATOR_DOFS, vec![FIRST, OUTLET], terminal, vec![]).unwrap();
        for flow in [1e-8, -2e-8, 0.0, 3e-8, 0.0, 0.0] {
            assert_eq!(empty.port_response().unwrap(), original.port_response().unwrap());
            assert_eq!(empty.preview_step(flow).unwrap().accept(), original.preview_step(flow).unwrap().accept());
            assert_eq!(empty.state(), original.state());
            assert_eq!(empty.terminal().unwrap().pole_states(), original.terminal().unwrap().pole_states());
            assert_eq!(empty.last_terminal_frames(), original.last_terminal_frames());
            assert_eq!(empty.stored_energy_j(), original.state().stored_energy_j());
            assert!(empty.last_branch_frames().is_empty());
        }
    }
}

#[test]
fn material_and_length_data_change_actual_external_flow_and_keep_outer_source_identity() {
    let reference = loaded(false, 1);
    let build = |length: f64, medium: GasState| {
        DgPortLoad::new_with_terminal_branches(main_duct(false, 1), DT, 1e-11,
            MAX_DG_PROPAGATOR_DOFS, vec![FIRST, OUTLET], prepared(&model(2), DT),
            vec![DgBranch { port: FIRST, duct: finite_duct(false, length, &medium) }]).unwrap()
    };
    let mut longer = build(0.024, gas(101325.0));
    let mixture = fs_material::gas_species::input::read_gas_record(
        &mut "mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2".split_ascii_whitespace(), None, None,
    ).unwrap();
    let mut changed_gas = build(0.012, *mixture.state());
    assert_eq!(changed_gas.terminal_gas(FIRST).unwrap(), *mixture.state());
    assert_eq!(changed_gas.terminal_gas(OUTLET).unwrap(), gas(101325.0));
    assert_eq!(changed_gas.terminal_area_m2(FIRST).unwrap(), Some(changed_gas.branches()[0].duct().outlet_area_m2()));
    assert!(changed_gas.terminal_gas(INLET).is_err());
    let mut reference = reference;
    let mut differences = [0.0; 2];
    for sample in 0..128 {
        let input = if sample < 8 { 1e-8 } else { 0.0 };
        reference.preview_step(input).unwrap().accept();
        longer.preview_step(input).unwrap().accept();
        changed_gas.preview_step(input).unwrap().accept();
        let expected = reference.last_terminal_frames()[0].outward_flow_m3_s;
        differences[0] += (longer.last_terminal_frames()[0].outward_flow_m3_s - expected).abs();
        differences[1] += (changed_gas.last_terminal_frames()[0].outward_flow_m3_s - expected).abs();
    }
    assert!(differences.iter().all(|value| *value > 1e-15));
}
