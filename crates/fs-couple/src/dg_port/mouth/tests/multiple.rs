//! Original uncondensed DG/branch/MIMO oracle for independently driven mouths.
use super::*;

fn many_model(ports: usize) -> MimoModel {
    let d: Vec<_> = (0..ports).flat_map(|i| (0..ports).map(move |j| {
        if i == j { 4e7 + i as f64 * 1e6 }
        else { 3e6 / (1 + i.abs_diff(j)) as f64 }
    })).collect();
    MimoModel {
        ports,
        poles: vec![PoleTerm::Real { pole: -400.0, residue: 0.0 },
            PoleTerm::Pair { pole: C64::new(-900.0, 2000.0), residue: C64::ZERO }],
        residues: vec![d.iter().map(|v| C64::from_re(100.0 * v)).collect(),
            d.iter().map(|v| C64::new(10.0 * v, 5.0 * v)).collect()],
        d, e: vec![0.0; ports * ports], reciprocal: true,
    }
}

fn many_load(branched: bool) -> (DgPortLoad, MimoModel) {
    let mut main = duct(0.05, if branched { 2 } else { 1 }, INLET, OUTLET);
    let mut branches = Vec::new();
    let ports = if branched {
        main = main.with_side_ports(&[DgSidePort { id: SIDE, downstream_element: 1 }]).unwrap();
        branches.push(DgBranch { port: SIDE, duct: duct(0.015, 1, PortId(0), PortId(1)) });
        vec![SIDE, OUTLET]
    } else { vec![OUTLET] };
    let mut load = DgPortLoad::new_with_terminal_branches(main, DT, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, ports.clone(), prepared(&many_model(ports.len())),
        branches).unwrap();
    let model = many_model(ports.len() + 2);
    load.install_independent_mouths(prepared(&model), 2).unwrap();
    (load, model)
}

#[test]
fn independent_mouths_match_all_uncondensed_histories_pressures_and_work() {
    for branched in [false, true] {
        let (mut load, model) = many_load(branched);
        let mut rational = model.zoh(DT).unwrap();
        let main = load.propagator().clone();
        let mut main_state = load.state().clone();
        let branch = load.branches().first().map(|b| b.propagator().clone());
        let mut branch_state = load.branches().first().map(|b| b.state().clone());
        for sample in 0..32 {
            let qb = if sample < 24 { (f64::from(sample % 7) - 3.0) * 1e-9 } else { 0.0 };
            let qm = [
                if sample < 17 { (f64::from(sample % 5) - 2.0) * 2e-9 } else { 0.0 },
                if (6..23).contains(&sample) { (3.0 - f64::from(sample % 7)) * 3e-9 } else { 0.0 },
            ];
            let dg = main.multiport_response(&main_state).unwrap();
            let b = branch.as_ref().map(|p|
                p.multiport_response(branch_state.as_ref().unwrap()).unwrap());
            let history = rational.clone().step_mean(&vec![0.0; model.ports]);
            let direct = load.terminal().unwrap().mean_direct_term().to_vec();
            let z = |row: usize, col: usize| direct[row * model.ports + col];
            let mouth_effort = |row| z(row, 0) * qm[0] + z(row, 1) * qm[1];
            // Unknowns are original branch junction flow and terminal flows.
            // These equations retain every original owner before condensation.
            let (matrix, rhs) = if let Some(b) = &b {
                (vec![
                    dg.impedance_pa_s_m3[8] + b.impedance_pa_s_m3[0],
                    -b.impedance_pa_s_m3[1], dg.impedance_pa_s_m3[7],
                    -b.impedance_pa_s_m3[2], b.impedance_pa_s_m3[3] + z(2, 2), z(2, 3),
                    dg.impedance_pa_s_m3[5], z(3, 2), dg.impedance_pa_s_m3[4] + z(3, 3),
                ], vec![
                    dg.history_pressure_pa[2] - b.history_pressure_pa[0]
                        + dg.impedance_pa_s_m3[6] * qb,
                    b.history_pressure_pa[1] - history[2] - mouth_effort(2),
                    dg.history_pressure_pa[1] - history[3]
                        + dg.impedance_pa_s_m3[3] * qb - mouth_effort(3),
                ])
            } else {
                (vec![dg.impedance_pa_s_m3[3] + z(2, 2)],
                    vec![dg.history_pressure_pa[1] - history[2]
                        + dg.impedance_pa_s_m3[2] * qb - mouth_effort(2)])
            };
            let mut lu = LuWorkspace::new(rhs.len()).unwrap();
            lu.factor(&matrix).unwrap();
            let mut expected = vec![0.0; rhs.len()];
            lu.solve_factored_into(&rhs, &mut expected).unwrap();
            let response = load.independent_mouths_response().unwrap();
            assert_eq!(response.history_pressure_pa.len(), 3);
            assert_eq!(response.impedance_pa_s_m3.len(), 9);
            assert!(response.impedance_pa_s_m3[5].abs() > 0.0
                && response.impedance_pa_s_m3[7].abs() > 0.0);
            let trial = load.preview_step_with_mouth_flows(qb, &qm).unwrap();
            let account = trial.independent_mouths_frame().unwrap().clone();
            let flows: Vec<_> = trial.terminal_frames().iter()
                .map(|f| f.outward_flow_m3_s).collect();
            let mut complete = qm.to_vec();
            complete.extend_from_slice(&flows);
            let pressure = rational.step_mean(&complete);
            for (i, mouth) in account.mouths.iter().enumerate() {
                assert_eq!(mouth.outward_flow_m3_s, qm[i]);
                assert_eq!(mouth.mean_pressure_pa, pressure[i]);
                near(mouth.load_work_j, pressure[i] * (qm[i] * DT),
                    mouth.load_work_j.abs() + (pressure[i] * qm[i] * DT).abs());
            }
            let controls = [qb, qm[0], qm[1]];
            let efforts = [account.inlet_pressure_pa,
                account.mouths[0].mean_pressure_pa, account.mouths[1].mean_pressure_pa];
            for row in 0..3 {
                let h = response.history_pressure_pa[row];
                let terms: Vec<_> = (0..3).map(|col|
                    response.impedance_pa_s_m3[3 * row + col] * controls[col]).collect();
                let predicted = h + terms.iter().sum::<f64>();
                near(efforts[row], predicted,
                    efforts[row].abs() + h.abs() + terms.iter().map(|v| v.abs()).sum::<f64>());
            }
            let mut inputs = vec![PortInput { id: INLET, effort: qb }];
            if branched {
                let junction = trial.branch_frames()[0].junction_flow_m3_s;
                near(junction, expected[0], junction.abs() + expected[0].abs());
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
            for (i, terminal) in trial.terminal_frames().iter().enumerate() {
                let q = expected[i + usize::from(branched)];
                near(terminal.outward_flow_m3_s, q, terminal.outward_flow_m3_s.abs() + q.abs());
                near(terminal.mean_pressure_pa, pressure[i + 2],
                    terminal.mean_pressure_pa.abs() + pressure[i + 2].abs());
            }
            let exterior_work: f64 = pressure.iter().zip(&complete)
                .map(|(p, q)| p * (q * DT)).sum();
            let duct_work = trial.output().ledger.port_work_step_j;
            near(trial.source_work_j(), duct_work + exterior_work,
                trial.source_work_j().abs() + duct_work.abs() + exterior_work.abs());
            assert_eq!(account.source_work_j, trial.source_work_j());
            trial.accept();
            assert_eq!(load.state(), &main_state);
            if let Some(state) = &branch_state { assert_eq!(load.branches()[0].state(), state); }
            assert_eq!(load.terminal().unwrap().pole_states(), rational.pole_states());
            assert_eq!(load.last_independent_mouths_frame(), Some(&account));
        }
    }
}

#[test]
fn plural_trials_cannot_guess_missing_flows_or_publish_abandoned_histories() {
    let (mut load, _) = many_load(true);
    load.preview_step_with_mouth_flows(1e-9, &[2e-9, -3e-9]).unwrap().accept();
    let before = load.clone();
    {
        let _trial = load.preview_step_with_mouth_flows(-1e-9, &[4e-9, 5e-9]).unwrap();
    }
    assert!(load.preview_step_with_mouth_flows(1e-9, &[2e-9]).is_err());
    assert!(load.preview_step_with_mouth_flows(1e-9, &[2e-9, f64::NAN]).is_err());
    assert!(load.preview_step(1e-9).is_err());
    assert!(load.independent_mouth_response().is_err());
    assert_eq!(load.state(), before.state());
    assert_eq!(load.branches()[0].state(), before.branches()[0].state());
    assert_eq!(load.terminal().unwrap().pole_states(), before.terminal().unwrap().pole_states());
    assert_eq!(load.last_independent_mouths_frame(), before.last_independent_mouths_frame());
}

#[test]
fn one_mouth_slice_preserves_original_scalar_independent_history() {
    let (mut scalar, _) = load(true);
    let mut slice = scalar.clone();
    for (qb, qm) in [(1e-9, 3e-9), (-2e-9, 1e-9), (0.0, 0.0)] {
        let a = scalar.preview_step_with_mouth_flow(qb, qm).unwrap();
        let b = slice.preview_step_with_mouth_flows(qb, &[qm]).unwrap();
        assert_eq!(a.source_work_j(), b.source_work_j());
        assert_eq!(a.inlet_work_j(), b.inlet_work_j());
        assert_eq!(a.output().ledger.port_work_step_j, b.output().ledger.port_work_step_j);
        a.accept();
        b.accept();
        assert_eq!(scalar.state(), slice.state());
        assert_eq!(scalar.branches()[0].state(), slice.branches()[0].state());
        assert_eq!(scalar.terminal().unwrap().pole_states(), slice.terminal().unwrap().pole_states());
    }
}
