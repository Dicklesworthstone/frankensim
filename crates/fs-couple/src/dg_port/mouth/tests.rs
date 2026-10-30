use super::*;
use crate::dg_port::{DgBranch, DgBoundary, DgDuct, PortId, PortInput};
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_la::prepared_lu::LuWorkspace;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_vfit::mimo::MimoModel;
use fs_vfit::model::PoleTerm;

const DT: f64 = 1.0 / 96000.0;
const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const SIDE: PortId = PortId(11);

fn model(ports: usize) -> MimoModel {
    // Positive definite D times a scalar positive-real real/pair response.
    // The pair can subtract at most 2*sqrt(10^2+5^2)/900 < 0.025 of D.
    let complete = [4e7, 3e6, 5e6, 3e6, 2e7, 4e6, 5e6, 4e6, 3e7];
    let d: Vec<_> = (0..ports).flat_map(|i| (0..ports).map(move |j| complete[3*i+j])).collect();
    MimoModel { ports, poles: vec![
        PoleTerm::Real { pole: -400.0, residue: 0.0 },
        PoleTerm::Pair { pole: C64::new(-900.0, 2000.0), residue: C64::ZERO },
    ], residues: vec![d.iter().map(|v| C64::from_re(100.0*v)).collect(),
        d.iter().map(|v| C64::new(10.0*v, 5.0*v)).collect()],
        d, e: vec![0.0; ports*ports], reciprocal: true }
}

fn prepared(model: &MimoModel) -> PreparedMimoLoad {
    PreparedMimoLoad::new(model, DT, MimoLoadLimits {
        band_rad_s: [1.0, core::f64::consts::PI/DT], sampled_passivity_points: 129,
        max_states: 32, max_memory_bytes: 1 << 20,
    }).unwrap()
}

fn duct(length: f64, elements: usize, inlet: PortId, outlet: PortId) -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    DgDuct::new(&Duct { segments: vec![Segment::Cylinder { radius: 0.001, length }] },
        &gas, DgMesh { elements, order: 2 }, DgBoundary::Port(inlet), DgBoundary::Port(outlet))
        .unwrap().with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
}

fn load(branched: bool) -> (DgPortLoad, MimoModel) {
    let mut main = duct(0.05, if branched { 2 } else { 1 }, INLET, OUTLET);
    let mut branches = Vec::new();
    let ports = if branched {
        main = main.with_side_ports(&[DgSidePort { id: SIDE, downstream_element: 1 }]).unwrap();
        branches.push(DgBranch { port: SIDE, duct: duct(0.015, 1, PortId(0), PortId(1)) });
        vec![SIDE, OUTLET]
    } else { vec![OUTLET] };
    let mut load = DgPortLoad::new_with_terminal_branches(main, DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, ports.clone(), prepared(&model(ports.len())), branches).unwrap();
    let joint = model(ports.len() + 1);
    load.install_driven_mouth(prepared(&joint)).unwrap();
    (load, joint)
}

fn near(a: f64, b: f64, scale: f64) {
    assert!((a-b).abs() <= 2e-9*scale.max(1e-300), "{a} != {b}, scale={scale}");
}

#[test]
fn joint_mouth_matches_uncondensed_pressure_equations_and_every_original_history() {
    for branched in [false, true] {
        let (mut load, model) = load(branched);
        let mut rational = model.zoh(DT).unwrap();
        let main = load.propagator().clone();
        let mut main_state = load.state().clone();
        let branch = load.branches().first().map(|branch| branch.propagator().clone());
        let mut branch_state = load.branches().first().map(|branch| branch.state().clone());
        for sample in 0..64 {
            let q = if sample < 23 { (f64::from(sample % 7)-3.0)*1e-9 } else { 0.0 };
            let dg = main.multiport_response(&main_state).unwrap();
            let b = branch.as_ref().map(|p| p.multiport_response(branch_state.as_ref().unwrap()).unwrap());
            let history = rational.clone().step_mean(&vec![0.0; model.ports]);
            let direct = load.terminal().unwrap().mean_direct_term().to_vec();
            // Independent UNCONDENSED unknowns are [J, Uside, Uout] or [Uout].
            // Each row is pressure continuity on an actual original owner.
            let (matrix, rhs) = if let Some(b) = &b {
                (vec![
                    dg.impedance_pa_s_m3[8]+b.impedance_pa_s_m3[0],
                    -b.impedance_pa_s_m3[1], dg.impedance_pa_s_m3[7],
                    -b.impedance_pa_s_m3[2],
                    b.impedance_pa_s_m3[3]+direct[4], direct[5],
                    dg.impedance_pa_s_m3[5], direct[7],
                    dg.impedance_pa_s_m3[4]+direct[8],
                ], vec![
                    dg.history_pressure_pa[2]-b.history_pressure_pa[0]+dg.impedance_pa_s_m3[6]*q,
                    b.history_pressure_pa[1]-history[1]+direct[3]*q,
                    dg.history_pressure_pa[1]-history[2]+(dg.impedance_pa_s_m3[3]+direct[6])*q,
                ])
            } else {
                (vec![dg.impedance_pa_s_m3[3]+direct[3]],
                    vec![dg.history_pressure_pa[1]-history[1]+(dg.impedance_pa_s_m3[2]+direct[2])*q])
            };
            let mut lu = LuWorkspace::new(rhs.len()).unwrap();
            lu.factor(&matrix).unwrap();
            let mut expected = vec![0.0; rhs.len()];
            lu.solve_factored_into(&rhs, &mut expected).unwrap();
            let relation = load.port_response().unwrap();
            let source_effort = relation.history_pressure_pa + relation.impedance_pa_s_m3*q;
            let trial = load.preview_step(q).unwrap();
            let mouth = *trial.driven_mouth_frame().unwrap();
            let terminal_flows: Vec<_> = trial.terminal_frames().iter().map(|f| f.outward_flow_m3_s).collect();
            let mut complete_flows = vec![-q];
            complete_flows.extend_from_slice(&terminal_flows);
            let pressure = rational.step_mean(&complete_flows);
            assert_eq!(mouth.outward_flow_m3_s, -q);
            assert_eq!(mouth.mean_pressure_pa, pressure[0]);
            near(source_effort, mouth.source_pressure_pa,
                source_effort.abs()+mouth.mean_pressure_pa.abs()+trial.inlet_work_j().abs()/(DT*q.abs()).max(1e-30));
            let mut inputs = vec![PortInput { id: INLET, effort: q }];
            if branched {
                let junction = trial.branch_frames()[0].junction_flow_m3_s;
                near(junction, expected[0], junction.abs()+expected[0].abs());
                inputs.push(PortInput { id: SIDE, effort: -junction });
                branch.as_ref().unwrap().step(branch_state.as_mut().unwrap(), &[
                    PortInput { id: PortId(0), effort: junction },
                    PortInput { id: PortId(1), effort: -terminal_flows[0] },
                ]).unwrap();
                branch.as_ref().unwrap().sync_nodal(branch_state.as_mut().unwrap()).unwrap();
            }
            inputs.push(PortInput { id: OUTLET, effort: -*terminal_flows.last().unwrap() });
            main.step(&mut main_state, &inputs).unwrap();
            main.sync_nodal(&mut main_state).unwrap();
            for (i, frame) in trial.terminal_frames().iter().enumerate() {
                let expected_flow = expected[i+usize::from(branched)];
                near(frame.outward_flow_m3_s, expected_flow, frame.outward_flow_m3_s.abs()+expected_flow.abs());
                near(frame.mean_pressure_pa, pressure[i+1], frame.mean_pressure_pa.abs()+pressure[i+1].abs());
            }
            let full_load_work: f64 = pressure.iter().zip(&complete_flows).map(|(p,q)| p*(q*DT)).sum();
            let dg_work = trial.output().ledger.port_work_step_j;
            near(trial.source_work_j(), dg_work+full_load_work,
                trial.source_work_j().abs()+dg_work.abs()+full_load_work.abs());
            trial.accept();
            assert_eq!(load.state(), &main_state);
            if let Some(state) = &branch_state { assert_eq!(load.branches()[0].state(), state); }
            assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
        }
    }
}

#[test]
fn abandoned_joint_trials_preserve_bore_branch_rational_and_published_frames() {
    let (mut load, model) = load(true);
    let before = load.clone();
    drop(load.preview_step(1e-9).unwrap());
    assert!(load.preview_step(f64::NAN).is_err());
    assert_eq!(load.state(), before.state());
    assert_eq!(load.branches()[0].state(), before.branches()[0].state());
    assert_eq!(load.terminal().unwrap().pole_states(), before.terminal().unwrap().pole_states());
    assert!(load.last_terminal_frames().is_empty());
    assert!(load.last_branch_frames().is_empty());
    assert!(load.last_driven_mouth_frame().is_none());
    // A joint model is attached once; replacement cannot erase accepted memory.
    assert!(load.install_driven_mouth(prepared(&model)).is_err());
    load.preview_step(1e-9).unwrap().accept();
    assert!(load.last_driven_mouth_frame().is_some());
    assert!(load.install_driven_mouth(prepared(&model)).is_err());
}

#[test]
fn joint_mouth_and_fixed_gap_channels_commute_at_construction_and_share_signed_work() {
    use crate::dg_port::{DgTerminalChannelInput, DgTerminalChannelKind};
    let controls = vec![
        DgTerminalChannelInput { port: SIDE, kind: DgTerminalChannelKind::Slit {
            width_m: 0.002, height_m: 4e-5, length_m: 0.001,
            transverse_cells: 16, maximum_mach: 0.2,
        } },
        DgTerminalChannelInput { port: OUTLET, kind: DgTerminalChannelKind::Closed },
    ];
    let (original, joint_model) = load(true);
    let mut mouth_first = original.with_terminal_channels(controls.clone(), 16).unwrap();
    let main = duct(0.05, 2, INLET, OUTLET)
        .with_side_ports(&[DgSidePort { id: SIDE, downstream_element: 1 }]).unwrap();
    let mut channel_first = DgPortLoad::new_with_terminal_branches(main, DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, vec![SIDE, OUTLET], prepared(&model(2)),
        vec![DgBranch { port: SIDE, duct: duct(0.015, 1, PortId(0), PortId(1)) }])
        .unwrap().with_terminal_channels(controls, 16).unwrap();
    channel_first.install_driven_mouth(prepared(&joint_model)).unwrap();
    let mut rational = joint_model.zoh(DT).unwrap();
    for sample in 0..32 {
        let q = if sample < 13 { (f64::from(sample % 5) - 2.0) * 1e-9 } else { 0.0 };
        assert_eq!(mouth_first.port_response().unwrap(), channel_first.port_response().unwrap());
        let trial = mouth_first.preview_step(q).unwrap();
        let flows = [-q, trial.terminal_frames()[0].outward_flow_m3_s, 0.0];
        assert_eq!(trial.terminal_frames()[1].outward_flow_m3_s, 0.0);
        let pressure = rational.step_mean(&flows);
        let exterior_work: f64 = pressure.iter().zip(flows).map(|(p, q)| p * (q * DT)).sum();
        let aggregate_work = trial.output().ledger.port_work_step_j;
        near(trial.source_work_j(), aggregate_work + exterior_work,
            trial.source_work_j().abs() + aggregate_work.abs() + exterior_work.abs());
        assert_eq!(trial.driven_mouth_frame().unwrap().mean_pressure_pa, pressure[0]);
        let actual = trial.accept();
        assert_eq!(actual, channel_first.preview_step(q).unwrap().accept());
        assert_eq!(mouth_first.last_terminal_channel_frames(), channel_first.last_terminal_channel_frames());
        assert_eq!(mouth_first.terminal().unwrap().pole_states(), rational.pole_states());
    }
}


#[test]
fn independent_mouth_matches_original_uncondensed_owners_and_mutual_history() {
    for branched in [false, true] {
        let (mut load, model) = load(branched);
        let mut rational = model.zoh(DT).unwrap();
        let main = load.propagator().clone();
        let mut main_state = load.state().clone();
        let branch = load.branches().first().map(|branch| branch.propagator().clone());
        let mut branch_state = load.branches().first().map(|branch| branch.state().clone());
        for sample in 0..48 {
            let qb = if sample < 27 { (f64::from(sample % 7)-3.0)*1e-9 } else { 0.0 };
            let qm = if sample < 19 { (f64::from(sample % 5)+1.0)*2e-9 } else { 0.0 };
            let dg = main.multiport_response(&main_state).unwrap();
            let b = branch.as_ref().map(|p| p.multiport_response(branch_state.as_ref().unwrap()).unwrap());
            let history = rational.clone().step_mean(&vec![0.0; model.ports]);
            let direct = load.terminal().unwrap().mean_direct_term().to_vec();
            // Independent uncondensed continuity unknowns [J, Uside, Uout]
            // retain each original main/branch/radiation owner.
            let (matrix, rhs) = if let Some(b) = &b {
                (vec![
                    dg.impedance_pa_s_m3[8]+b.impedance_pa_s_m3[0],
                    -b.impedance_pa_s_m3[1], dg.impedance_pa_s_m3[7],
                    -b.impedance_pa_s_m3[2],
                    b.impedance_pa_s_m3[3]+direct[4], direct[5],
                    dg.impedance_pa_s_m3[5], direct[7],
                    dg.impedance_pa_s_m3[4]+direct[8],
                ], vec![
                    dg.history_pressure_pa[2]-b.history_pressure_pa[0]+dg.impedance_pa_s_m3[6]*qb,
                    b.history_pressure_pa[1]-history[1]-direct[3]*qm,
                    dg.history_pressure_pa[1]-history[2]+dg.impedance_pa_s_m3[3]*qb-direct[6]*qm,
                ])
            } else {
                (vec![dg.impedance_pa_s_m3[3]+direct[3]],
                    vec![dg.history_pressure_pa[1]-history[1]+dg.impedance_pa_s_m3[2]*qb-direct[2]*qm])
            };
            let mut lu = LuWorkspace::new(rhs.len()).unwrap();
            lu.factor(&matrix).unwrap();
            let mut expected = vec![0.0; rhs.len()];
            lu.solve_factored_into(&rhs, &mut expected).unwrap();
            let response = load.independent_mouth_response().unwrap();
            let legacy = load.port_response().unwrap();
            let h = response.history_pressure_pa;
            let z = response.impedance_pa_s_m3;
            near(legacy.history_pressure_pa, h[0]-h[1],
                legacy.history_pressure_pa.abs()+h[0].abs()+h[1].abs());
            near(legacy.impedance_pa_s_m3, z[0][0]-z[0][1]-z[1][0]+z[1][1],
                legacy.impedance_pa_s_m3.abs()+z.iter().flatten().map(|x| x.abs()).sum::<f64>());
            assert!(z[0][1].abs() > 0.0 && z[1][0].abs() > 0.0,
                "self-only mouths must not replace mutual radiation");
            let trial = load.preview_step_with_mouth_flow(qb, qm).unwrap();
            let mouth = *trial.independent_mouth_frame().unwrap();
            assert!(trial.driven_mouth_frame().is_none());
            let flows: Vec<_> = trial.terminal_frames().iter().map(|f| f.outward_flow_m3_s).collect();
            let mut complete_flows = vec![qm];
            complete_flows.extend_from_slice(&flows);
            let pressure = rational.step_mean(&complete_flows);
            assert_eq!(mouth.outward_flow_m3_s, qm);
            assert_eq!(mouth.mean_pressure_pa, pressure[0]);
            for (row, actual) in [mouth.inlet_pressure_pa, mouth.mean_pressure_pa].into_iter().enumerate() {
                let predicted = h[row]+z[row][0]*qb+z[row][1]*qm;
                near(actual, predicted, actual.abs()+h[row].abs()
                    +(z[row][0]*qb).abs()+(z[row][1]*qm).abs());
            }
            let mut inputs = vec![PortInput { id: INLET, effort: qb }];
            if branched {
                let junction = trial.branch_frames()[0].junction_flow_m3_s;
                near(junction, expected[0], junction.abs()+expected[0].abs());
                inputs.push(PortInput { id: SIDE, effort: -junction });
                branch.as_ref().unwrap().step(branch_state.as_mut().unwrap(), &[
                    PortInput { id: PortId(0), effort: junction },
                    PortInput { id: PortId(1), effort: -flows[0] },
                ]).unwrap();
                branch.as_ref().unwrap().sync_nodal(branch_state.as_mut().unwrap()).unwrap();
            }
            inputs.push(PortInput { id: OUTLET, effort: -*flows.last().unwrap() });
            main.step(&mut main_state, &inputs).unwrap();
            main.sync_nodal(&mut main_state).unwrap();
            for (i, frame) in trial.terminal_frames().iter().enumerate() {
                let q = expected[i+usize::from(branched)];
                near(frame.outward_flow_m3_s, q, frame.outward_flow_m3_s.abs()+q.abs());
                near(frame.mean_pressure_pa, pressure[i+1], frame.mean_pressure_pa.abs()+pressure[i+1].abs());
            }
            let exterior_work: f64 = pressure.iter().zip(&complete_flows).map(|(p,q)| p*(q*DT)).sum();
            let dg_work = trial.output().ledger.port_work_step_j;
            near(trial.source_work_j(), dg_work+exterior_work,
                trial.source_work_j().abs()+dg_work.abs()+exterior_work.abs());
            assert_eq!(mouth.source_work_j, trial.source_work_j());
            trial.accept();
            assert_eq!(load.state(), &main_state);
            if let Some(state) = &branch_state { assert_eq!(load.branches()[0].state(), state); }
            assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
            assert_eq!(load.last_independent_mouth_frame(), Some(&mouth));
        }
    }
}

#[test]
fn independent_mouth_trials_rollback_and_legacy_steps_clear_their_record() {
    let (mut load, _) = load(true);
    load.preview_step_with_mouth_flow(1e-9, 3e-9).unwrap().accept();
    let before = load.clone();
    drop(load.preview_step_with_mouth_flow(-2e-9, 4e-9).unwrap());
    for (qb, qm) in [(f64::NAN, 0.0), (0.0, f64::NAN), (0.0, f64::INFINITY)] {
        assert!(load.preview_step_with_mouth_flow(qb, qm).is_err());
    }
    assert_eq!(load.state(), before.state());
    assert_eq!(load.branches()[0].state(), before.branches()[0].state());
    assert_eq!(load.terminal().unwrap().pole_states(), before.terminal().unwrap().pole_states());
    assert_eq!(load.last_terminal_frames(), before.last_terminal_frames());
    assert_eq!(load.last_independent_mouth_frame(), before.last_independent_mouth_frame());
    load.preview_step(1e-9).unwrap().accept();
    assert!(load.last_independent_mouth_frame().is_none());
    assert!(load.last_driven_mouth_frame().is_some());
}

#[test]
fn independent_mouth_retains_finite_gap_work_and_closed_terminal_rows() {
    use crate::dg_port::{DgTerminalChannelInput, DgTerminalChannelKind};
    let (original, model) = load(true);
    let mut load = original.with_terminal_channels(vec![
        DgTerminalChannelInput { port: SIDE, kind: DgTerminalChannelKind::Slit {
            width_m: 0.002, height_m: 4e-5, length_m: 0.001,
            transverse_cells: 16, maximum_mach: 0.2,
        } },
        DgTerminalChannelInput { port: OUTLET, kind: DgTerminalChannelKind::Closed },
    ], 16).unwrap();
    let mut rational = model.zoh(DT).unwrap();
    for sample in 0..24 {
        let qb = (f64::from(sample % 5)-2.0)*1e-9;
        let qm = if sample < 13 { 3e-9 } else { -1e-9 };
        let trial = load.preview_step_with_mouth_flow(qb, qm).unwrap();
        let mouth = *trial.independent_mouth_frame().unwrap();
        let flows = [qm, trial.terminal_frames()[0].outward_flow_m3_s, 0.0];
        assert_eq!(trial.terminal_frames()[1].outward_flow_m3_s, 0.0);
        let pressure = rational.step_mean(&flows);
        assert_eq!(mouth.mean_pressure_pa, pressure[0]);
        let exterior_work: f64 = pressure.iter().zip(flows).map(|(p,q)| p*(q*DT)).sum();
        let original_work = trial.output().ledger.port_work_step_j;
        near(trial.source_work_j(), original_work+exterior_work,
            trial.source_work_j().abs()+original_work.abs()+exterior_work.abs());
        trial.accept();
        assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
    }
}
