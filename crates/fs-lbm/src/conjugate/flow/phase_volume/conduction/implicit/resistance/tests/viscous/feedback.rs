//! Source-density changes affect viscosity without replacing original inertia.

use super::*;

#[test]
fn thermal_mushy_viscosity_keeps_original_inertia_and_shared_heat_history() {
    let (card, curve) = fixture_with_densities([4.0, 2.0]);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 4.0, 2.0, 1.0).unwrap();
    let context = context();
    let conductivity = prepare(&card, &curve, &context, [40.0, 90.0]);
    let wall_conductivity = prepare(&card, &curve, &context, [110.0, 130.0]);
    let mut domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    domain.set(0, 1, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let initial = curve.state_at_specific_enthalpy(50.0).unwrap();
    let phases = [Some(initial), None];
    let mass = [initial.bulk_density_kg_m3(), 0.0];
    let momentum = [[0.1 * mass[0], 0.0, 0.0], [0.0; 3]];
    let heat = [mass[0] * 10.0, 0.0];
    let force = [[0.0; 3]; 2];
    let mut input = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &heat,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: [None; 6],
    };
    let mut faces = [FlowFace::Wall; 6];
    let mut boundaries = [FvBoundary::Symmetry; 6];
    for side in 0..2 {
        faces[side] = FlowFace::Free;
        boundaries[side] = FvBoundary::Outlet;
        boundaries[2 + side] = FvBoundary::wall();
        input.volume.boundary_inflow[side] = Some(initial);
        input.boundary_velocity_m_s[side] = Some([0.0; 3]);
    }
    let px = [Some(0.0), Some(0.0), None, None];
    let py = [None; 3];
    let pz = [None; 4];
    let reservoir = Some((
        faces,
        StaggeredBoundaryPressure {
            x_pa: &px,
            y_pa: &py,
            z_pa: &pz,
        },
    ));
    let solids = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall_conductivity,
        mass_kg: 2.0,
        phase: curve.state_at_specific_enthalpy(116.5).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let contact_resistance = 16.0 / 15.0;
    let contacts = [ThermalVolumeContact::Resistance(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: contact_resistance,
    })];
    let limits = limits(2, 1);
    let solve = ImplicitThermalVolumeSolveControls {
        cell_tolerance_j: 1e-12,
        ..ImplicitThermalVolumeSolveControls::default()
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let molecular = MaterialViscosityControls {
        max_sweeps: 40_000,
        ..molecular_controls()
    };
    let feedback = ThermalViscosityFeedbackControls {
        maximum_trials: 1024,
        relaxation: 1.0,
        maximum_heat_defect_j: 1e-11,
    };
    let run = |workspace: &mut ImplicitThermalVolumeWorkspace,
               solve,
               molecular,
               feedback,
               gate: &CancelGate,
               mode: u8| {
        let mut queries = 0;
        let result = field
            .advance_isotropic_implicit_conjugate_thermoviscous_phase_resistance_thermal_volume_momentum_with_baths(
                &domain, reservoir, &conductivity, ThermalVolumeBoundaries::default(),
                &solids, &contacts, &law, input, boundaries, limits, molecular,
                solve, feedback, workspace, gate,
                |_, phase| {
                    queries += 1;
                    // Native caloric solve has run before the third query.
                    if mode == 1 && queries == 3 { gate.request(); }
                    let source_phase = if mode == 2 && queries == 3 { &initial } else { phase };
                    Ok((source(&card, &law, source_phase)?, 0.25 + phase.liquid_mass_fraction()))
                },
            );
        (result, queries)
    };
    let (result, queries) = run(
        &mut workspace,
        solve,
        molecular,
        feedback,
        &CancelGate::new_clock_free(),
        0,
    );
    let result = result.unwrap();
    let coupled = &result.evolution;
    let combined = &coupled.evolution;
    let native = &combined.evolution.evolution;
    let volume = &native.thermal.volume;
    let source = combined.viscosity_sources[0].as_ref().unwrap();
    let caloric = source.phase_state();
    let mu = source.mobile_transport().unwrap().dynamic_viscosity_pa_s();
    assert!(caloric.bulk_density_kg_m3() < mass[0]);
    assert_eq!(source.source_density_kg_m3(), caloric.bulk_density_kg_m3());
    close(mu, 0.5 * caloric.bulk_density_kg_m3());
    close(
        result.bulk_viscosity_pa_s[0],
        0.25 + caloric.liquid_mass_fraction(),
    );
    assert_eq!(result.feedback.source_evaluations, queries);
    assert!(combined.viscosity.molecular_heat_j > 0.0);
    assert!(native.resistance.resistance_heat_j > 0.0);
    assert!(coupled.conductive_enthalpy_increment_j[0] > 0.0);
    assert!(result.feedback.trials > native.coupling.coupling_iterations);
    assert!(result.feedback.molecular_solves > native.coupling.coupling_iterations);
    assert_eq!(coupled.work.heat_trials, result.feedback.trials);
    assert_eq!(
        coupled.work.residual_evaluations,
        result.feedback.heat_residual_evaluations
    );

    // Independent native molecular kick at the accepted exact source, retaining
    // ORIGINAL M/P. Using its lower hot rho as old inertia changes this result.
    let reference = field
        .prepare_material_viscosity_with_source(
            &domain,
            faces,
            Some(0),
            &mass,
            &momentum,
            boundaries,
            Some(&result.bulk_viscosity_pa_s),
            limits.resistance.thermal.mechanics.reconstruction,
            0.25,
            molecular,
            &CancelGate::new_clock_free(),
            |_| Ok(vec![mu, 0.0]),
        )
        .unwrap();
    assert_eq!(combined.viscosity, reference.report);
    assert_eq!(combined.viscosity_reconstruction, reference.reconstruction);
    let wrong_mass = [caloric.bulk_density_kg_m3(), 0.0];
    let wrong = field
        .prepare_material_viscosity_with_source(
            &domain,
            faces,
            Some(0),
            &wrong_mass,
            &momentum,
            boundaries,
            Some(&result.bulk_viscosity_pa_s),
            limits.resistance.thermal.mechanics.reconstruction,
            0.25,
            molecular,
            &CancelGate::new_clock_free(),
            |_| Ok(vec![mu, 0.0]),
        )
        .unwrap();
    assert!((wrong.report.molecular_heat_j - combined.viscosity.molecular_heat_j).abs() > 1e-7);

    let used_heat = heat[0]
        + result.used_molecular_heat_j[0]
        + native.caloric_resistance_heat_j[0]
        + coupled.conductive_enthalpy_increment_j[0];
    close(caloric.specific_enthalpy_j_kg(), 50.0 + used_heat / mass[0]);
    close(volume.report.requested_enthalpy_increment_j, used_heat);
    close(volume.volume_source_m3_s[0], used_heat / (400.0 * 0.25));
    let endpoint = volume.phase_states[0].unwrap();
    close(
        volume.transport.mass_kg[0] / endpoint.bulk_density_kg_m3(),
        1.0,
    );
    let physical_heat =
        heat[0] + combined.viscosity.molecular_heat_j + native.resistance.resistance_heat_j;
    close(
        volume.transport.mass_kg[0] * endpoint.specific_enthalpy_j_kg() - mass[0] * 50.0
            + coupled.report.solid_enthalpy_change_j,
        physical_heat + volume.transport.report.carried_enthalpy_j,
    );
    let wall = &coupled.solid_endpoints[0];
    close(
        wall.mass_kg * (wall.phase.specific_enthalpy_j_kg() - 116.5),
        -coupled.conductive_enthalpy_increment_j[0],
    );
    close(
        volume.transport.mass_kg[0] - mass[0],
        volume.transport.report.incoming_mass_kg - volume.transport.report.outgoing_mass_kg,
    );
    assert!(result.feedback.absolute_heat_coupling_defect_j <= feedback.maximum_heat_defect_j);
    assert_eq!(
        coupled.report.conduction.absolute_enthalpy_defect_j,
        combined.absolute_enthalpy_defect_j
    );
    assert!(combined.absolute_energy_defect_j <= molecular.maximum_energy_defect_j);
    assert_eq!(native.momentum.momentum_kg_m_s[1], [0.0; 3]);
    assert!(combined.viscosity_sources[1].is_none());

    // The previous frozen-phase entry still solves its declared split, but
    // cannot reproduce the thermally resolved viscosity and physical heat.
    let frozen = field
        .advance_isotropic_implicit_conjugate_viscous_phase_resistance_thermal_volume_momentum_with_baths(
            &domain, reservoir, &conductivity, ThermalVolumeBoundaries::default(),
            &solids, &contacts, &law, input, boundaries,
            Some(&[0.25 + initial.liquid_mass_fraction(), 0.0]), limits,
            molecular, solve, &mut workspace, &CancelGate::new_clock_free(),
            |_, phase| super::source(&card, &law, phase),
        ).unwrap();
    assert!(
        (frozen.evolution.viscosity.molecular_heat_j - combined.viscosity.molecular_heat_j).abs()
            > 1e-7
    );

    // Refuse stale endpoint source evidence, cancellation after a real heat
    // solve, and exhaustion of EACH whole-call work ceiling. A retry must use
    // original finite-wall history and produce the identical accepted result.
    let (cancelled, _) = run(
        &mut workspace,
        solve,
        molecular,
        feedback,
        &CancelGate::new_clock_free(),
        1,
    );
    assert!(matches!(cancelled, Err(ChtError::Cancelled)));
    assert!(
        run(
            &mut workspace,
            solve,
            molecular,
            feedback,
            &CancelGate::new_clock_free(),
            2
        )
        .0
        .is_err()
    );
    let short_heat = ImplicitThermalVolumeSolveControls {
        maximum_residual_evaluations: result.feedback.heat_residual_evaluations - 1,
        ..solve
    };
    assert!(
        run(
            &mut workspace,
            short_heat,
            molecular,
            feedback,
            &CancelGate::new_clock_free(),
            0
        )
        .0
        .is_err()
    );
    let short_molecular = MaterialViscosityControls {
        max_sweeps: result.feedback.molecular_sweeps - 1,
        ..molecular
    };
    assert!(
        run(
            &mut workspace,
            solve,
            short_molecular,
            feedback,
            &CancelGate::new_clock_free(),
            0
        )
        .0
        .is_err()
    );
    let short_trials = ThermalViscosityFeedbackControls {
        maximum_trials: result.feedback.trials - 1,
        ..feedback
    };
    assert!(
        run(
            &mut workspace,
            solve,
            molecular,
            short_trials,
            &CancelGate::new_clock_free(),
            0
        )
        .0
        .is_err()
    );
    let (retry, retry_queries) = run(
        &mut workspace,
        solve,
        molecular,
        feedback,
        &CancelGate::new_clock_free(),
        0,
    );
    assert_eq!(retry.unwrap(), result);
    assert_eq!(retry_queries, queries);
}
