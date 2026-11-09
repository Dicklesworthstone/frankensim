//! Genuine T/rate source surfaces must close on the same native heat endpoint.

mod sgs;

use super::*;

fn joint_response(
    material: &MaterialCard,
    law: &MushyPhaseLaw<'_>,
    h: f64,
    rate: f64,
) -> Result<(ResolvedPhaseTransportState, f64), ChtError> {
    let source = resolve_phase_transport_state_with_shear_for_volume_change(
        law,
        h,
        PhaseTransportSourceQuery {
            source: ThermalMaterialPoint {
                material,
                context: &material_point(301.0),
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            viscosity_basis: LiquidViscosityBasis::Kinematic,
            extra_requirements: &[],
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        "rate",
        rate,
    )
    .map_err(|error| refuse(error.to_string()))?;
    // Independent joint response, not inferred from molecular shear viscosity.
    Ok((
        source,
        1.0 / 3.0 + 5.0 * (h - 1.0) + 3.0 * 3.0_f64.sqrt() * rate / 16.0,
    ))
}

#[test]
fn joint_temperature_rate_sources_close_exact_momentum_heat_and_finite_wall() {
    let (material, melting) = thermal_material_with_rate(true);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    for finite_wall in [false, true] {
        let n = if finite_wall { 4 } else { 3 };
        let mut domain = VoxelDomain::new(n, 1, 1, 1.0).unwrap();
        if finite_wall {
            domain.set(3, 0, 0, Voxel::Solid(0));
        }
        let field = FlowField::quiescent(&domain);
        let mut mass = vec![2.0; n];
        let mut phases = vec![Some(caloric.state_at_specific_enthalpy(1.0).unwrap()); n];
        let mut momentum = vec![[0.0; 3]; n];
        momentum[0][0] = 2.0;
        momentum[2][0] = -2.0;
        let mut authored = vec![0.0; n];
        let force = vec![[0.0; 3]; n];
        let mut solids = vec![];
        if finite_wall {
            mass[3] = 0.0;
            phases[3] = None;
            authored[2] = -1.0 / 16.0;
            solids.push(ThermalVolumeSolidCell {
                cell: 3,
                conductivity: &conductivity,
                mass_kg: 2.0,
                phase: caloric.state_at_specific_enthalpy(9.0 / 8.0).unwrap(),
                net_enthalpy_increment_j: 0.0,
            });
        }
        let calls = Cell::new(0);
        let result = field
            .advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths(
                &domain,
                None,
                &conductivity,
                ThermalVolumeBoundaries::default(),
                &solids,
                &[],
                ThermalVolumeMomentumInput {
                    volume: input(&mass, &phases, &authored),
                    momentum_kg_m_s: &momentum,
                    acceleration_m_s2: [0.0; 3],
                    force_density_n_m3: &force,
                    boundary_velocity_m_s: [None; 6],
                },
                [FvBoundary::Symmetry; 6],
                super::super::super::super::momentum::coupled_controls(0.5),
                molecular_controls(),
                thermal_solve_controls(),
                feedback_controls(),
                &mut ImplicitThermalVolumeWorkspace::new(n).unwrap(),
                &CancelGate::new_clock_free(),
                |cell, phase, rate| {
                    assert!(cell < 3);
                    calls.set(calls.get() + 1);
                    joint_response(&material, &law, phase.specific_enthalpy_j_kg(), rate)
                },
            )
            .unwrap();
        // Native axial faces begin at +1/2,-1/2. At the joint endpoint their
        // magnitude b=1/4, so RMS dev-strain is 2b/sqrt(3) at the edges and
        // 4b/sqrt(3) centrally. The original two rows give
        // b=(1/2)/(1+lambda_edge/4+lambda_center/2), lambda=4mu/3+zeta.
        // mu=1/8+4*(h-1)+(sqrt(3)/8)*rate and
        // zeta=1/3+5*(h-1)+(3sqrt(3)/16)*rate give lambda=[1,3/2,1].
        // Hence Q=dt*lambda*(du/dx)^2=[1/32,3/16,1/32], and the
        // original heat rows 2*(h-1)=Q+sum_neighbors(h_neighbor-h) close
        // at [33/32,17/16,33/32]. Freezing T, rate, or both changes this.
        let expected_h = [33.0 / 32.0, 17.0 / 16.0, 33.0 / 32.0];
        let expected_rate = [
            1.0 / (2.0 * 3.0_f64.sqrt()),
            1.0 / 3.0_f64.sqrt(),
            1.0 / (2.0 * 3.0_f64.sqrt()),
        ];
        let expected_mu = [5.0 / 16.0, 0.5, 5.0 / 16.0];
        let expected_bulk = [7.0 / 12.0, 5.0 / 6.0, 7.0 / 12.0];
        let expected_heat = [1.0 / 32.0, 3.0 / 16.0, 1.0 / 32.0];
        let accepted = &result.evolution;
        let viscous = &accepted.evolution.evolution.viscosity;
        let volume = &accepted.evolution.evolution.evolution.thermal.volume;
        let mut signed = 0.0;
        let mut absolute = 0.0;
        for cell in 0..3 {
            near(result.endpoint_shear_rate_per_s[cell], expected_rate[cell]);
            let source = accepted.endpoint_sources[cell].as_ref().unwrap();
            near(
                source.phase_state().specific_enthalpy_j_kg(),
                expected_h[cell],
            );
            near(
                source.transport_state().unwrap().dynamic_viscosity_pa_s(),
                expected_mu[cell],
            );
            near(accepted.bulk_viscosity_pa_s[cell], expected_bulk[cell]);
            // These are original source receipts at BOTH accepted coordinates,
            // never the initial predictor or the later transported velocity.
            for property in source.material().properties() {
                for (axis, expected) in [
                    ("T", 300.0 + expected_h[cell]),
                    ("rate", expected_rate[cell]),
                ] {
                    let coordinate = property
                        .answer()
                        .receipt
                        .query_point
                        .iter()
                        .find(|(name, _)| name == axis)
                        .unwrap()
                        .1;
                    near(coordinate, expected);
                }
            }
            near(viscous.cell_molecular_heat_j[cell], expected_heat[cell]);
            near(accepted.used_molecular_heat_j[cell], expected_heat[cell]);
            near(
                volume.phase_states[cell].unwrap().specific_enthalpy_j_kg(),
                expected_h[cell],
            );
            near(volume.transport.mass_kg[cell], 2.0);
            near(
                accepted.evolution.conductive_enthalpy_increment_j[cell],
                2.0 * (expected_h[cell] - 1.0) - authored[cell] - expected_heat[cell],
            );
            let mismatch =
                accepted.used_molecular_heat_j[cell] - viscous.cell_molecular_heat_j[cell];
            signed += mismatch;
            absolute += mismatch.abs();
        }
        near(viscous.molecular_heat_j, 1.0 / 4.0);
        near(viscous.backward_euler_loss_j, 1.0 / 8.0);
        assert!(viscous.relative_residual <= molecular_controls().relative_tolerance);
        let report = &accepted.evolution.report;
        near(
            report.conduction.solve.energy.external_heat_j,
            1.0 / 4.0 + authored.iter().sum::<f64>(),
        );
        near(
            volume.report.material_enthalpy_change_j + report.solid_enthalpy_change_j,
            1.0 / 4.0 + authored.iter().sum::<f64>(),
        );
        if finite_wall {
            near(
                accepted.evolution.solid_endpoints[0]
                    .phase
                    .specific_enthalpy_j_kg(),
                35.0 / 32.0,
            );
            near(report.solid_to_fluid_heat_j, 1.0 / 16.0);
            near(report.solid_enthalpy_change_j, -1.0 / 16.0);
            assert!(accepted.endpoint_sources[3].is_none());
            assert_eq!(result.endpoint_shear_rate_per_s[3], 0.0);
            assert_eq!(accepted.bulk_viscosity_pa_s[3], 0.0);
            assert_eq!(accepted.used_molecular_heat_j[3], 0.0);
            assert_eq!(solids[0].phase.specific_enthalpy_j_kg(), 9.0 / 8.0);
        }
        let work = &accepted.feedback;
        assert!(work.trials > 1);
        assert_eq!(work.molecular_solves, work.trials + 1);
        assert_eq!(
            result.rheology_source_updates,
            work.molecular_sweeps + work.molecular_solves
        );
        assert_eq!(calls.get(), 3 * result.rheology_source_updates);
        assert_eq!(work.source_evaluations, calls.get());
        assert!(work.molecular_sweeps > work.molecular_solves);
        assert!(work.molecular_sweeps <= molecular_controls().max_sweeps);
        assert!(work.molecular_krylov_iterations >= viscous.krylov_iterations);
        assert!(
            work.heat_residual_evaluations <= thermal_solve_controls().maximum_residual_evaluations
        );
        assert!(work.heat_linear_iterations <= thermal_solve_controls().linear.maximum_iterations);
        assert!((work.heat_coupling_defect_j - signed).abs() < 1e-15);
        assert!((work.absolute_heat_coupling_defect_j - absolute).abs() < 1e-15);
        assert!(work.absolute_heat_coupling_defect_j <= feedback_controls().maximum_heat_defect_j);
        assert!(work.absolute_enthalpy_defect_j <= 1e-9);
        assert_eq!(field, FlowField::quiescent(&domain));
        assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 1.0);
        assert_eq!(momentum[0], [2.0, 0.0, 0.0]);
    }
}

#[test]
fn joint_inner_work_is_whole_call_and_thermal_callback_cancellation_is_atomic() {
    let (material, melting) = thermal_material_with_rate(true);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 3];
    let phases = [Some(caloric.state_at_specific_enthalpy(1.0).unwrap()); 3];
    let momentum = [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0]];
    let authored = [0.0; 3];
    let forces = [[0.0; 3]; 3];
    let calls = Cell::new(0);
    let run = |viscosity: MaterialViscosityControls,
               cancel_at_thermal_endpoint: bool,
               workspace: &mut ImplicitThermalVolumeWorkspace,
               gate: &CancelGate| {
        field.advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths(
            &domain,
            None,
            &conductivity,
            ThermalVolumeBoundaries::default(),
            &[],
            &[],
            ThermalVolumeMomentumInput {
                volume: input(&mass, &phases, &authored),
                momentum_kg_m_s: &momentum,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &forces,
                boundary_velocity_m_s: [None; 6],
            },
            [FvBoundary::Symmetry; 6],
            super::super::super::super::momentum::coupled_controls(0.5),
            viscosity,
            thermal_solve_controls(),
            feedback_controls(),
            workspace,
            gate,
            |_, phase, rate| {
                calls.set(calls.get() + 1);
                if cancel_at_thermal_endpoint && phase.specific_enthalpy_j_kg() > 1.0 {
                    gate.request();
                    return Err(refuse("source refusal after requested cancellation"));
                }
                joint_response(&material, &law, phase.specific_enthalpy_j_kg(), rate)
            },
        )
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(3).unwrap();
    let gate = CancelGate::new_clock_free();
    let expected = run(molecular_controls(), false, &mut workspace, &gate).unwrap();
    let total = expected.evolution.feedback.molecular_sweeps;
    assert!(total > expected.evolution.evolution.evolution.viscosity.sweeps + 1);
    let bounded = MaterialViscosityControls {
        max_sweeps: total - 1,
        ..molecular_controls()
    };
    calls.set(0);
    assert!(run(bounded, false, &mut workspace, &gate).is_err());
    // The ceiling exceeds a final trial's work. Resetting it for each thermal
    // or rheological trial would incorrectly publish a successful endpoint.
    assert!(calls.get() > 3);
    let cancelled = CancelGate::new_clock_free();
    calls.set(0);
    assert!(matches!(
        run(molecular_controls(), true, &mut workspace, &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert!(calls.get() > 3);
    assert_eq!(mass, [2.0; 3]);
    assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 1.0);
    assert_eq!(momentum, [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0]]);
    assert_eq!(field, FlowField::quiescent(&domain));
    calls.set(0);
    let retried = run(molecular_controls(), false, &mut workspace, &gate).unwrap();
    assert_eq!(retried, expected);
    assert_eq!(calls.get(), retried.evolution.feedback.source_evaluations);
}
