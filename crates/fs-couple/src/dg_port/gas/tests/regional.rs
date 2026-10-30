use super::*;

fn input(field: &AxialGasField, refresh: DgGasRefresh) -> DgBranchGasInput<'_> {
    DgBranchGasInput {
        port: SIDE,
        field,
        refresh,
    }
}

fn warm(channel: DgTerminalChannelInput, diffusion: bool) -> DgPortLoad {
    let mut result = load(&field(300.0, 101_325.0, 0.05));
    if diffusion {
        let branch = result.branches()[0]
            .duct()
            .clone()
            .with_viscothermal(8, [50.0, 2000.0], 0.01)
            .unwrap();
        result = DgPortLoad::new_with_terminal_branches(
            result.duct().clone(),
            DT,
            1e-11,
            MAX_DG_PROPAGATOR_DOFS,
            vec![SIDE, OUTLET],
            result.terminal().unwrap().clone(),
            vec![DgBranch {
                port: SIDE,
                duct: branch,
            }],
        )
        .unwrap();
    }
    result = result
        .with_terminal_channels(vec![channel, slit(OUTLET)], 16)
        .unwrap();
    for _ in 0..30 {
        result.preview_step(1e-9).unwrap().accept();
    }
    result
}

fn unchanged(actual: &DgPortLoad, expected: &DgPortLoad) {
    assert_eq!(actual.state(), expected.state());
    assert_eq!(
        actual.duct().operator_identity(),
        expected.duct().operator_identity()
    );
    assert_eq!(actual.duct().gas_field(), expected.duct().gas_field());
    assert_eq!(
        actual.port_response().unwrap(),
        expected.port_response().unwrap()
    );
    assert_eq!(
        actual.cumulative_parameter_work_j(),
        expected.cumulative_parameter_work_j()
    );
    assert_eq!(actual.last_branch_frames(), expected.last_branch_frames());
    assert_eq!(
        actual.last_terminal_frames(),
        expected.last_terminal_frames()
    );
    assert_eq!(
        actual.last_terminal_channel_frames(),
        expected.last_terminal_channel_frames()
    );
    assert_eq!(
        actual.terminal().unwrap().pole_states(),
        expected.terminal().unwrap().pole_states()
    );
    for (a, b) in actual.branches().iter().zip(expected.branches()) {
        assert_eq!(a.state(), b.state());
        assert_eq!(a.duct().operator_identity(), b.duct().operator_identity());
        assert_eq!(a.duct().gas_field(), b.duct().gas_field());
    }
    for (a, b) in actual
        .terminal_channels()
        .iter()
        .zip(expected.terminal_channels())
    {
        assert_eq!(a.gas(), b.gas());
        assert_eq!(
            a.channel().unwrap().velocities_m_s(),
            b.channel().unwrap().velocities_m_s()
        );
        assert_eq!(
            a.channel().unwrap().last_step(),
            b.channel().unwrap().last_step()
        );
        assert_eq!(
            a.cumulative_parameter_work_j(),
            b.cumulative_parameter_work_j()
        );
    }
}

#[test]
fn g1_branch_gas_preserves_original_compression_momentum_memory_and_whole_load_work() {
    let after = field(450.0, 101_325.0, 0.012);
    let circular = DgTerminalChannelInput {
        port: SIDE,
        kind: DgTerminalChannelKind::Circular {
            radius_m: 8e-5,
            length_m: 0.002,
            radial_cells: 8,
            maximum_mach: 0.2,
        },
    };
    for (channel, diffusion) in [(slit(SIDE), false), (circular, false), (slit(SIDE), true)] {
        let mut live = warm(channel, diffusion);
        let accepted = live.clone();
        let old = &accepted.branches()[0];
        let refresh = if diffusion {
            DgGasRefresh::UniformDiffusion {
                maximum_relative_error: 0.01,
            }
        } else {
            DgGasRefresh::Lossless
        };
        let report = live
            .rebind_regional_gas_fields_controlled(
                None,
                &[input(&after, refresh)],
                1.0,
                &mut || false,
            )
            .unwrap();
        let new = &live.branches()[0];
        assert!(report.main.is_none());
        assert_eq!(report.branches.len(), 1);
        assert_eq!(report.branches[0].port, SIDE);
        assert_eq!(live.state(), accepted.state());
        assert_eq!(
            live.duct().operator_identity(),
            accepted.duct().operator_identity()
        );
        assert_eq!(new.duct().gas_field(), Some(&after));
        assert_eq!(new.state().step_count(), old.state().step_count());
        assert_eq!(
            new.state().initial_energy_j(),
            old.state().initial_energy_j()
        );
        assert_eq!(
            new.state().cumulative_dissipated_j(),
            old.state().cumulative_dissipated_j()
        );
        let old_gas = old.duct().inlet_gas();
        let new_gas = new.duct().inlet_gas();
        let old_bulk = old_gas.density * old_gas.sound_speed.powi(2);
        let new_bulk = new_gas.density * new_gas.sound_speed.powi(2);
        let area = old.duct().inlet_area_m2();
        let pressure_scale = old
            .state()
            .pressure()
            .iter()
            .map(|p| p.abs())
            .fold(0.0, f64::max);
        let flow_scale = old
            .state()
            .volume_flow()
            .iter()
            .map(|q| q.abs())
            .fold(0.0, f64::max);
        for (&p, &previous) in new.state().pressure().iter().zip(old.state().pressure()) {
            near(p / new_bulk, previous / old_bulk, pressure_scale / old_bulk);
        }
        for (&q, &previous) in new
            .state()
            .volume_flow()
            .iter()
            .zip(old.state().volume_flow())
        {
            near(
                new_gas.density * q / area,
                old_gas.density * previous / area,
                old_gas.density * flow_scale / area,
            );
        }
        // The public native transfer provides a reference for the complete
        // diffusion history, including arms inaccessible through nodal p/q.
        if diffusion {
            let mut reference_duct = old.duct().clone();
            let mut reference_state = old.state().clone();
            let native = reference_duct
                .rebind_uniform_viscothermal_gas_field(
                    &mut reference_state,
                    &after,
                    0.01,
                    1.0,
                    &mut || false,
                )
                .unwrap();
            assert_eq!(new.state(), &reference_state);
            assert_eq!(report.branches[0].duct, native);
        }
        let old_channel = &accepted.terminal_channels()[0];
        let new_channel = &live.terminal_channels()[0];
        let mut old_energy = 0.0;
        let mut new_energy = 0.0;
        for ((volume, &v0), &v1) in volumes(channel.kind)
            .iter()
            .zip(old_channel.channel().unwrap().velocities_m_s())
            .zip(new_channel.channel().unwrap().velocities_m_s())
        {
            let momentum = old_gas.density * volume * v0;
            near(new_gas.density * volume * v1, momentum, momentum.abs());
            old_energy += 0.5 * old_gas.density * volume * v0 * v0;
            new_energy += 0.5 * new_gas.density * volume * v1 * v1;
        }
        assert!(old_energy > 0.0);
        near(old_channel.stored_energy_j(), old_energy, old_energy);
        near(new_channel.stored_energy_j(), new_energy, new_energy);
        let channel_work = report.terminal_channels[0]
            .transfer
            .unwrap()
            .parameter_work_j;
        near(
            channel_work,
            new_energy - old_energy,
            new_energy + old_energy,
        );
        assert_eq!(
            report.terminal_channels[1]
                .transfer
                .unwrap()
                .parameter_work_j,
            0.0
        );
        assert_eq!(
            live.terminal_channels()[1].gas(),
            accepted.terminal_channels()[1].gas()
        );
        assert_eq!(
            live.terminal().unwrap().pole_states(),
            accepted.terminal().unwrap().pole_states()
        );
        assert_eq!(live.last_branch_frames(), accepted.last_branch_frames());
        assert_eq!(
            live.last_terminal_channel_frames(),
            accepted.last_terminal_channel_frames()
        );
        assert_eq!(
            new_channel.channel().unwrap().last_step(),
            old_channel.channel().unwrap().last_step()
        );
        assert_eq!(
            new_channel.cumulative_pressure_work_j(),
            old_channel.cumulative_pressure_work_j()
        );
        assert_eq!(
            report.parameter_work_j,
            report.branches[0].duct.parameter_work_j + channel_work
        );
        near(
            report.energy_after_j - report.energy_before_j,
            report.parameter_work_j,
            report.energy_after_j + report.energy_before_j,
        );
        assert_eq!(live.cumulative_parameter_work_j(), report.parameter_work_j);
        assert_ne!(live.impedance_pa_s_m3(), accepted.impedance_pa_s_m3());
        for _ in 0..12 {
            let output = live.preview_step(-2e-10).unwrap().accept();
            assert_eq!(
                output.ledger.cumulative_parameter_work_j,
                report.parameter_work_j
            );
            near(
                output.ledger.energy_balance_residual_j,
                0.0,
                output.ledger.stored_energy_j
                    + output.ledger.cumulative_port_work_j.abs()
                    + output.ledger.cumulative_dissipated_j
                    + report.absolute_parameter_work_j,
            );
            let frame = &live.last_branch_frames()[0];
            near(
                frame.bore_port.conjugate_rate,
                frame.branch_inlet.conjugate_rate,
                frame.bore_port.conjugate_rate.abs() + frame.branch_inlet.conjugate_rate.abs(),
            );
        }
    }
}

#[test]
fn g3_coordinated_pressure_update_admits_only_the_complete_regional_endpoint() {
    let main = field(420.0, 120_000.0, 0.05);
    let branch = field(450.0, 120_000.0, 0.012);
    let accepted = warm(slit(SIDE), false);
    let mut live = accepted.clone();
    assert!(
        live.rebind_gas_field_controlled(&main, DgGasRefresh::Lossless, 1.0, &mut || false,)
            .is_err()
    );
    unchanged(&live, &accepted);
    assert!(
        live.rebind_regional_gas_fields_controlled(
            None,
            &[input(&branch, DgGasRefresh::Lossless)],
            1.0,
            &mut || false,
        )
        .is_err()
    );
    unchanged(&live, &accepted);
    let report = live
        .rebind_regional_gas_fields_controlled(
            Some((&main, DgGasRefresh::Lossless)),
            &[input(&branch, DgGasRefresh::Lossless)],
            1.0,
            &mut || false,
        )
        .unwrap();
    assert_eq!(live.duct().gas_field(), Some(&main));
    assert_eq!(live.branches()[0].duct().gas_field(), Some(&branch));
    assert_eq!(live.state().step_count(), accepted.state().step_count());
    let contributions = [
        report.main.unwrap().parameter_work_j,
        report.branches[0].duct.parameter_work_j,
        report.terminal_channels[0]
            .transfer
            .unwrap()
            .parameter_work_j,
        report.terminal_channels[1]
            .transfer
            .unwrap()
            .parameter_work_j,
    ];
    assert_eq!(report.parameter_work_j, contributions.iter().sum::<f64>());
    assert_eq!(
        report.absolute_parameter_work_j,
        contributions.iter().map(|w| w.abs()).sum::<f64>()
    );
    let output = live.preview_step(0.0).unwrap().accept();
    assert_eq!(
        output.ledger.cumulative_parameter_work_j,
        report.parameter_work_j
    );
    near(
        output.ledger.energy_balance_residual_j,
        0.0,
        output.ledger.stored_energy_j
            + output.ledger.cumulative_dissipated_j
            + output.ledger.cumulative_port_work_j.abs()
            + report.absolute_parameter_work_j,
    );
}

#[test]
fn g4_regional_identity_budget_and_final_cancellation_refusals_preserve_exact_retry() {
    let after = field(450.0, 101_325.0, 0.012);
    let selected = input(&after, DgGasRefresh::Lossless);
    let accepted = warm(slit(SIDE), false);
    let mut reference = accepted.clone();
    let mut polls = 0;
    let expected = reference
        .rebind_regional_gas_fields_controlled(None, &[selected], 1.0, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let allowance = expected.branches[0].duct.parameter_work_j.abs().max(
        expected.terminal_channels[0]
            .transfer
            .unwrap()
            .parameter_work_j
            .abs(),
    );
    assert!(expected.absolute_parameter_work_j > allowance);
    for refusal in 0..4 {
        let mut live = accepted.clone();
        let mut requested = vec![selected];
        if refusal == 0 {
            requested.push(selected);
        }
        if refusal == 1 {
            requested[0].port = PortId(999);
        }
        let mut count = 0;
        assert!(
            live.rebind_regional_gas_fields_controlled(
                None,
                &requested,
                if refusal == 2 { allowance } else { 1.0 },
                &mut || {
                    count += 1;
                    refusal == 3 && count == polls
                },
            )
            .is_err()
        );
        if refusal == 3 {
            assert_eq!(count, polls);
        }
        unchanged(&live, &accepted);
        assert_eq!(
            live.rebind_regional_gas_fields_controlled(None, &[selected], 1.0, &mut || false,)
                .unwrap(),
            expected
        );
        let mut control = reference.clone();
        for _ in 0..8 {
            assert_eq!(
                live.preview_step(-1e-10).unwrap().accept(),
                control.preview_step(-1e-10).unwrap().accept()
            );
        }
        unchanged(&live, &control);
    }
}
