use super::*;
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_la::prepared_lu::LuWorkspace;
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

fn duct(gas: &GasState, length: f64, elements: usize, inlet: PortId,
    outlet: PortId, walls: bool) -> DgDuct {
    let source = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length }] },
        gas, DgMesh { elements, order: 2 },
        DgBoundary::Port(inlet), DgBoundary::Port(outlet),
    ).unwrap();
    if walls { source.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap() }
    else { source }
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

fn load(walls: bool, branch_gas: Option<GasState>) -> DgPortLoad {
    // Connected ducts share equilibrium pressure, including pressure-change fixtures.
    let pressure = branch_gas.as_ref().map_or(101325.0, |gas| gas.pressure);
    let main_gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let main = duct(&main_gas, 0.05, 2, INLET, OUTLET, walls)
        .with_side_ports(&[DgSidePort { id: SIDE, downstream_element: 1 }]).unwrap();
    let branches = branch_gas.map_or_else(Vec::new, |gas| vec![DgBranch {
        port: SIDE, duct: duct(&gas, 0.012, 1, PortId(0), PortId(1), walls),
    }]);
    let terminal = PreparedMimoLoad::new(&model(), DT, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI / DT], sampled_passivity_points: 129,
        max_states: 16, max_memory_bytes: 1 << 20,
    }).unwrap();
    DgPortLoad::new_with_terminal_branches(main, DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, vec![SIDE, OUTLET], terminal, branches).unwrap()
}

fn slit(port: PortId, height_m: f64) -> DgTerminalChannelInput {
    DgTerminalChannelInput { port, kind: DgTerminalChannelKind::Slit {
        width_m: 0.002, height_m, length_m: 0.001, transverse_cells: 16,
        maximum_mach: 0.2,
    } }
}

fn circular(port: PortId, radius_m: f64) -> DgTerminalChannelInput {
    DgTerminalChannelInput { port, kind: DgTerminalChannelKind::Circular {
        radius_m, length_m: 0.002, radial_cells: 16, maximum_mach: 0.2,
    } }
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2e-9 * scale.max(1e-300),
        "actual={actual}, expected={expected}, scale={scale}");
}

#[test]
fn full_gap_group_matches_uncondensed_pressure_rows_and_original_fluid_histories() {
    for (walls, round) in [(false, false), (true, false), (false, true), (true, true)] {
        for branched in [false, true] {
            let channels = if round { vec![circular(SIDE, 8e-5), slit(OUTLET, 7e-5)] }
                else { vec![slit(SIDE, 4e-5), slit(OUTLET, 7e-5)] };
            let mut load = load(walls, branched.then(air)).with_terminal_channels(
                channels, 32).unwrap();
            let main = load.propagator().clone();
            let mut main_state = load.state().clone();
            let branch = load.branches().first().map(|branch| branch.propagator().clone());
            let mut branch_state = load.branches().first().map(|branch| branch.state().clone());
            let mut channels: Vec<_> = load.terminal_channels().iter()
                .map(|channel| channel.channel().unwrap().clone()).collect();
            let mut rational = model().zoh(DT).unwrap();
            let direct = load.terminal().unwrap().mean_direct_term().to_vec();
            for sample in 0..48 {
                let q = if sample < 19 { (f64::from(sample % 7) - 3.0) * 1e-9 } else { 0.0 };
                let dg = main.multiport_response(&main_state).unwrap();
                let b = branch.as_ref().map(|p| p.multiport_response(branch_state.as_ref().unwrap()).unwrap());
                let history = rational.clone().step_mean(&[0.0; 2]);
                let gap: Vec<_> = channels.iter().map(|channel| channel.pressure_load().unwrap()).collect();
                // Original physical equations, without using the runtime's
                // condensed response: [junction, gap-side, gap-out] or both gaps.
                let (matrix, rhs) = if let Some(b) = &b {
                    (vec![
                        dg.impedance_pa_s_m3[8] + b.impedance_pa_s_m3[0],
                        -b.impedance_pa_s_m3[1], dg.impedance_pa_s_m3[7],
                        -b.impedance_pa_s_m3[2],
                        b.impedance_pa_s_m3[3] + direct[0] + gap[0].impedance_pa_s_m3, direct[1],
                        dg.impedance_pa_s_m3[5], direct[2],
                        dg.impedance_pa_s_m3[4] + direct[3] + gap[1].impedance_pa_s_m3,
                    ], vec![
                        dg.history_pressure_pa[2] - b.history_pressure_pa[0] + dg.impedance_pa_s_m3[6] * q,
                        b.history_pressure_pa[1] - history[0] + gap[0].history_pressure_pa,
                        dg.history_pressure_pa[1] - history[1] + dg.impedance_pa_s_m3[3] * q + gap[1].history_pressure_pa,
                    ])
                } else {
                    (vec![dg.impedance_pa_s_m3[8] + direct[0] + gap[0].impedance_pa_s_m3,
                        dg.impedance_pa_s_m3[7] + direct[1],
                        dg.impedance_pa_s_m3[5] + direct[2],
                        dg.impedance_pa_s_m3[4] + direct[3] + gap[1].impedance_pa_s_m3],
                    vec![dg.history_pressure_pa[2] - history[0] + dg.impedance_pa_s_m3[6] * q + gap[0].history_pressure_pa,
                        dg.history_pressure_pa[1] - history[1] + dg.impedance_pa_s_m3[3] * q + gap[1].history_pressure_pa])
                };
                let mut lu = LuWorkspace::new(rhs.len()).unwrap();
                lu.factor(&matrix).unwrap();
                let mut expected = vec![0.0; rhs.len()];
                lu.solve_factored_into(&rhs, &mut expected).unwrap();
                let previous_energy = load.stored_energy_j();
                let previous_loss = load.cumulative_dissipated_j();
                let relation = load.port_response().unwrap();
                let trial = load.preview_step(q).unwrap();
                let flows: Vec<_> = trial.terminal_frames().iter().map(|frame| frame.outward_flow_m3_s).collect();
                for (i, &flow) in flows.iter().enumerate() {
                    let expected_flow = expected[i + usize::from(branched)];
                    near(flow, expected_flow, flow.abs() + expected_flow.abs());
                }
                let exterior_pressure = rational.step_mean(&flows);
                let mut inputs = vec![PortInput { id: INLET, effort: q },
                    PortInput { id: OUTLET, effort: -flows[1] }];
                let mut energy = 0.0;
                if branched {
                    let junction = trial.branch_frames()[0].junction_flow_m3_s;
                    near(junction, expected[0], junction.abs() + expected[0].abs());
                    inputs.push(PortInput { id: SIDE, effort: -junction });
                    let step = branch.as_ref().unwrap().step(branch_state.as_mut().unwrap(), &[
                        PortInput { id: PortId(0), effort: junction },
                        PortInput { id: PortId(1), effort: -flows[0] },
                    ]).unwrap();
                    branch.as_ref().unwrap().sync_nodal(branch_state.as_mut().unwrap()).unwrap();
                    assert_eq!(trial.branch_frames()[0].duct, step);
                    energy += step.ledger.stored_energy_j;
                } else { inputs.push(PortInput { id: SIDE, effort: -flows[0] }); }
                let main_step = main.step(&mut main_state, &inputs).unwrap();
                main.sync_nodal(&mut main_state).unwrap();
                assert_eq!(trial.state(), &main_state);
                energy += main_step.ledger.stored_energy_j;
                for i in 0..2 {
                    let drop = channels[i].pressure_load().unwrap().pressure_drop_pa(flows[i]).unwrap();
                    let step = channels[i].step(drop).unwrap();
                    assert_eq!(trial.terminal_channel_frames()[i].step, Some(step));
                    let frame = &trial.terminal_frames()[i];
                    near(frame.mean_pressure_pa - exterior_pressure[i], drop,
                        frame.mean_pressure_pa.abs() + exterior_pressure[i].abs() + drop.abs());
                    assert_eq!(frame.channel_work_j, step.pressure_work_j);
                    near(frame.port_work_residual_j, 0.0,
                        frame.port_work_j.abs() + frame.channel_work_j.abs() + frame.load_work_j.abs());
                    energy += step.kinetic_energy_j;
                }
                near(trial.output().ledger.stored_energy_j, energy, energy);
                let inlet = trial.output().port_outputs.iter().find(|port| port.id == INLET).unwrap();
                near(inlet.conjugate_rate, relation.history_pressure_pa + relation.impedance_pa_s_m3 * q,
                    inlet.conjugate_rate.abs() + relation.history_pressure_pa.abs() + (relation.impedance_pa_s_m3 * q).abs());
                let external_work = sum_terminal_work(trial.terminal_frames());
                let storage_change = trial.output().ledger.stored_energy_j - previous_energy;
                let loss = trial.output().ledger.cumulative_dissipated_j - previous_loss;
                near(storage_change + loss + external_work, trial.inlet_work_j(),
                    storage_change.abs() + loss.abs() + external_work.abs() + trial.inlet_work_j().abs());
                trial.accept();
                for (actual, expected) in load.terminal_channels().iter().zip(&channels) {
                    assert_eq!(actual.channel().unwrap().velocities_m_s(), expected.velocities_m_s());
                    assert_eq!(actual.channel().unwrap().last_step(), expected.last_step());
                }
                assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
            }
            assert!(load.terminal_channels().iter().all(|channel|
                channel.cumulative_viscous_dissipation_j() > 0.0
                && channel.cumulative_numerical_dissipation_j() > 0.0));
        }
    }
}

#[test]
fn exact_closure_keeps_mutual_exterior_history_and_all_closed_has_no_infinite_impedance() {
    for all_closed in [false, true] {
        let mut controls = vec![DgTerminalChannelInput { port: SIDE, kind: DgTerminalChannelKind::Closed }];
        if all_closed { controls.push(DgTerminalChannelInput { port: OUTLET, kind: DgTerminalChannelKind::Closed }); }
        let mut load = load(false, None).with_terminal_channels(controls, 0).unwrap();
        let mut rational = model().zoh(DT).unwrap();
        let mut reaction = 0.0;
        for q in [1e-8, -2e-8, 0.0, 3e-8, 0.0] {
            let response = load.propagator().multiport_response(load.state()).unwrap();
            if all_closed {
                assert_eq!(load.port_response().unwrap().impedance_pa_s_m3.to_bits(), response.impedance_pa_s_m3[0].to_bits());
            }
            let trial = load.preview_step(q).unwrap();
            let frames = trial.terminal_frames();
            assert_eq!(frames[0].outward_flow_m3_s, 0.0);
            if all_closed { assert_eq!(frames[1].outward_flow_m3_s, 0.0); }
            let pressure = rational.step_mean(&[frames[0].outward_flow_m3_s, frames[1].outward_flow_m3_s]);
            reaction += (frames[0].mean_pressure_pa - pressure[0]).abs();
            assert_eq!(trial.terminal_channel_frames()[0].step, None);
            assert_eq!(frames[0].port_work_j, 0.0);
            assert_eq!(frames[0].load_work_j, 0.0);
            assert_eq!(frames[0].channel_work_j, 0.0);
            trial.accept();
            assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
            assert_eq!(load.terminal().unwrap().ports(), 2);
        }
        assert!(reaction > 0.0);
        if !all_closed { assert!(rational.pole_states().iter().any(|state| *state != C64::ZERO)); }
    }
}

#[test]
fn actual_branch_gas_and_gap_geometry_change_source_impedance_and_outer_flow() {
    let mixture = fs_material::gas_species::input::read_gas_record(
        &mut "mixture 303.15 101325 mole 2 N2 0.6 CO2 0.4".split_ascii_whitespace(), None, None,
    ).unwrap();
    let build = |gas, height| load(false, Some(gas)).with_terminal_channels(vec![slit(SIDE, height)], 16).unwrap();
    let mut reference = build(air(), 4e-5);
    let mut gas_change = build(*mixture.state(), 4e-5);
    let mut geometry_change = build(air(), 8e-5);
    assert_eq!(gas_change.terminal_channels()[0].gas(), mixture.state());
    assert_eq!(gas_change.terminal_channels()[0].channel().unwrap().gas(), mixture.state());
    assert_ne!(reference.port_response().unwrap(), gas_change.port_response().unwrap());
    assert_ne!(reference.port_response().unwrap(), geometry_change.port_response().unwrap());
    let mut differences = [0.0; 2];
    for sample in 0..48 {
        let q = if sample < 9 { 1e-8 } else { 0.0 };
        reference.preview_step(q).unwrap().accept();
        gas_change.preview_step(q).unwrap().accept();
        geometry_change.preview_step(q).unwrap().accept();
        let value = reference.last_terminal_frames()[0].outward_flow_m3_s;
        differences[0] += (gas_change.last_terminal_frames()[0].outward_flow_m3_s - value).abs();
        differences[1] += (geometry_change.last_terminal_frames()[0].outward_flow_m3_s - value).abs();
    }
    assert!(differences.iter().all(|difference| *difference > 1e-15));
}

#[test]
fn circular_passage_keeps_its_local_mixture_pressure_and_radial_history_in_the_joint_solve() {
    let mixture = fs_material::gas_species::input::read_gas_record(
        &mut "mixture 303.15 101325 mole 2 N2 0.6 CO2 0.4".split_ascii_whitespace(), None, None,
    ).unwrap();
    let dense = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202650.0).unwrap();
    let build = |gas, radius| load(true, Some(gas))
        .with_terminal_channels(vec![circular(SIDE, radius)], 16).unwrap();
    let mut reference = build(air(), 8e-5);
    let mut changed = [build(*mixture.state(), 8e-5), build(dense, 8e-5), build(air(), 1.2e-4)];
    for (load, gas) in changed.iter().zip([*mixture.state(), dense, air()]) {
        let channel = load.terminal_channels()[0].channel().unwrap().circular().unwrap();
        assert_eq!(*channel.gas(), gas);
        assert_eq!(channel.spec().radial_cells, 16);
        assert_ne!(load.port_response().unwrap(), reference.port_response().unwrap());
    }
    let mut differences = [0.0; 3];
    for sample in 0..48 {
        let q = if sample < 9 { 1e-8 } else { 0.0 };
        reference.preview_step(q).unwrap().accept();
        let expected = reference.last_terminal_frames()[0].outward_flow_m3_s;
        for (load, difference) in changed.iter_mut().zip(&mut differences) {
            load.preview_step(q).unwrap().accept();
            *difference += (load.last_terminal_frames()[0].outward_flow_m3_s - expected).abs();
            let frame = load.last_terminal_channel_frames()[0];
            let step = frame.step.unwrap();
            near(frame.mean_pressure_drop_pa, frame.constitutive_pressure_drop_pa,
                frame.mean_pressure_drop_pa.abs() + frame.constitutive_pressure_drop_pa.abs());
            near(step.energy_residual_j, 0.0, step.pressure_work_j.abs()
                + step.stored_delta_j.abs() + step.viscous_dissipation_j
                + step.numerical_dissipation_j);
        }
    }
    assert!(differences.iter().all(|difference| *difference > 1e-15));
    assert!(reference.terminal_channels()[0].channel().unwrap().velocities_m_s()
        .windows(2).any(|pair| pair[0] != pair[1]), "radial shear must remain resolved");
    let accepted = reference.clone();
    assert!(reference.preview_step(1.0).is_err());
    assert_eq!(reference.state(), accepted.state());
    assert_eq!(reference.last_terminal_channel_frames(), accepted.last_terminal_channel_frames());
    assert_eq!(reference.terminal_channels()[0].channel().unwrap().velocities_m_s(),
        accepted.terminal_channels()[0].channel().unwrap().velocities_m_s());
}

#[test]
fn dropped_and_late_mach_refused_trials_preserve_every_original_field_and_retry() {
    let mut load = load(true, Some(air())).with_terminal_channels(vec![slit(SIDE, 4e-5)], 16).unwrap();
    for q in [1e-8, -2e-8, 3e-8] { load.preview_step(q).unwrap().accept(); }
    let mut clean = load.clone();
    drop(load.preview_step(1e-9).unwrap());
    for q in [f64::NAN, f64::INFINITY, f64::MAX, 1.0] {
        assert!(load.preview_step(q).is_err());
        assert_eq!(load.state(), clean.state());
        assert_eq!(load.branches()[0].state(), clean.branches()[0].state());
        assert_eq!(load.terminal().unwrap().pole_states(), clean.terminal().unwrap().pole_states());
        assert_eq!(load.last_terminal_channel_frames(), clean.last_terminal_channel_frames());
        assert_eq!(load.terminal_channels()[0].channel().unwrap().velocities_m_s(),
            clean.terminal_channels()[0].channel().unwrap().velocities_m_s());
        assert_eq!(load.port_response().unwrap(), clean.port_response().unwrap());
    }
    assert_eq!(load.preview_step(2e-9).unwrap().accept(), clean.preview_step(2e-9).unwrap().accept());
    assert_eq!(load.last_terminal_channel_frames(), clean.last_terminal_channel_frames());
    assert!(load.with_terminal_channels(vec![slit(OUTLET, 4e-5)], 16).is_err());
}

#[test]
fn retained_gap_kinetic_energy_can_return_signed_work_without_becoming_negative_heat() {
    let mut load = load(false, None).with_terminal_channels(vec![
        slit(SIDE, 4e-5),
        DgTerminalChannelInput { port: OUTLET, kind: DgTerminalChannelKind::Closed },
    ], 16).unwrap();
    load.preview_step(1e-8).unwrap().accept();
    let gap = load.terminal_channels()[0].channel().unwrap().pressure_load().unwrap();
    let target = 0.5 * gap.history_pressure_pa / gap.impedance_pa_s_m3;
    assert_ne!(target, 0.0);
    let terminal = load.terminal_relation(load.terminal().unwrap()).unwrap();
    let coefficient = load.terminal_coupling.as_ref().unwrap().outward_per_inlet[0];
    let source_flow = (target - terminal.outward_history_m3_s[0]) / coefficient;
    let previous_energy = load.stored_energy_j();
    let previous_loss = load.cumulative_dissipated_j();
    let trial = load.preview_step(source_flow).unwrap();
    let step = trial.terminal_channel_frames()[0].step.unwrap();
    assert!(step.pressure_work_j < 0.0);
    assert!(step.stored_delta_j < 0.0);
    assert!(step.viscous_dissipation_j >= 0.0);
    assert!(step.numerical_dissipation_j >= 0.0);
    near(step.pressure_work_j, step.stored_delta_j + step.viscous_dissipation_j
        + step.numerical_dissipation_j,
        step.pressure_work_j.abs() + step.stored_delta_j.abs()
            + step.viscous_dissipation_j + step.numerical_dissipation_j);
    let storage_change = trial.output().ledger.stored_energy_j - previous_energy;
    let loss = trial.output().ledger.cumulative_dissipated_j - previous_loss;
    let external = sum_terminal_work(trial.terminal_frames());
    near(storage_change + loss + external, trial.inlet_work_j(),
        storage_change.abs() + loss.abs() + external.abs() + trial.inlet_work_j().abs());
    trial.accept();
}

#[test]
fn source_preflight_checks_thin_geometry_unique_ports_and_combined_field_budget() {
    let source = load(false, None);
    let input = slit(SIDE, 4e-5);
    assert!(source.clone().with_terminal_channels(vec![input], 15).is_err());
    assert!(source.clone().with_terminal_channels(vec![input], 16).is_ok());
    assert!(source.clone().with_terminal_channels(vec![input, input], 32).is_err());
    assert!(source.clone().with_terminal_channels(vec![input], MAX_DG_TERMINAL_TRANSVERSE_CELLS + 1).is_err());
    for port in [INLET, PortId(999)] {
        assert!(source.clone().with_terminal_channels(vec![DgTerminalChannelInput { port, ..input }], 16).is_err());
    }
    for height in [0.0, -1.0, f64::NAN, f64::INFINITY, 2e-4] {
        assert!(slit(SIDE, height).validate().is_err());
    }
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY, 3e-4] {
        assert!(circular(SIDE, radius).validate().is_err());
    }
    assert!(source.clone().with_terminal_channels(vec![circular(SIDE, 8e-5)], 15).is_err());
    assert!(source.clone().with_terminal_channels(vec![circular(SIDE, 8e-5), input], 32).is_err());
    let wide = DgTerminalChannelInput { port: SIDE, kind: DgTerminalChannelKind::Circular {
        radius_m: 0.002, length_m: 0.04, radial_cells: 16, maximum_mach: 0.2,
    } };
    assert!(wide.validate().is_ok());
    assert!(source.clone().with_terminal_channels(vec![wide], 16).is_err(),
        "a passage cannot exceed its actual attached duct section");
    let mut empty = source.clone().with_terminal_channels(Vec::new(), 0).unwrap();
    let mut original = source;
    for q in [1e-8, 0.0, -2e-8, 0.0] {
        assert_eq!(empty.port_response().unwrap(), original.port_response().unwrap());
        assert_eq!(empty.preview_step(q).unwrap().accept(), original.preview_step(q).unwrap().accept());
        assert_eq!(empty.last_terminal_frames(), original.last_terminal_frames());
    }
}
