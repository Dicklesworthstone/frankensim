//! Exposed wall baths use the receiving material and its finite storage.

use super::*;

#[test]
fn finite_pcm_wall_baths_close_separate_wall_and_whole_system_energy() {
    let fluid_card = card("T", &[(300.0, 2.0), (400.0, 2.0)]);
    let fluid_curve = curve(
        &fluid_card,
        &[(0.0, 300.0, 1.0, 2.0), (100.0, 400.0, 1.0, 2.0)],
    );
    let fluid = source(&fluid_card, &fluid_curve, [0.0, 100.0]);
    let solid_card = card("f", &[(0.0, 1.0), (1.0, 2.0)]);
    let solid_curve = curve(
        &solid_card,
        &[
            (-20.0, 300.0, 0.0, 2.0),
            (0.0, 320.0, 0.0, 2.0),
            (100.0, 320.0, 1.0, 2.0),
            (180.0, 400.0, 1.0, 2.0),
        ],
    );
    let solid = source(&solid_card, &solid_curve, [-20.0, 180.0]);
    let mut domain = VoxelDomain::graded([vec![1.0, 2.0], vec![1.0], vec![1.0]]).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let mass = [2.0, 0.0];
    let phases = [
        Some(fluid_curve.state_at_specific_enthalpy(9.5).unwrap()),
        None,
    ];
    let zero = [0.0; 2];
    let data = input(&mass, &phases, &zero);
    let walls = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &solid,
        mass_kg: 4.0,
        phase: solid_curve.state_at_specific_enthalpy(48.0).unwrap(),
        net_enthalpy_increment_j: 5.0,
    }];
    let contacts = [ThermalVolumeContact::Resistance(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.25,
    })];
    let baths = [
        bath(0, Face3::XMin, 314.0, 0.75),
        bath(1, Face3::XMax, 334.0, 0.5),
        bath(1, Face3::YMin, 319.0, 1.0 / 3.0),
    ];
    // Independent endpoint: h_f*=13, h_s*=50, T_f*=313, T_s*=320,
    // f_s*=1/2, k_s*=3/2. A_x=1, A_y=2, d_s,x=1 and d_s,y=1/2.
    // G_interface=1/(.5/2+.25+1/1.5)=6/7, giving Q_s->f=6 J at dt=1.
    // Solid baths give +12 J on x+ and -3 J on y-; the fluid bath gives +1 J.
    // M_f*delta_h=7=6+1; M_s*delta_h=8=5+12-3-6. Only physical heat enters.
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let result = field
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            &baths,
            &walls,
            &contacts,
            data,
            controls(1.0),
            solve_controls(),
            &mut workspace,
            &gate,
        )
        .unwrap();
    let returned = &result.solid_endpoints[0];
    let report = &result.report;
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .specific_enthalpy_j_kg(),
        13.0,
    );
    close(returned.phase.specific_enthalpy_j_kg(), 50.0);
    close(returned.phase.temperature_k(), 320.0);
    close(returned.phase.liquid_mass_fraction(), 0.5);
    close(result.conductive_enthalpy_increment_j[0], 7.0);
    assert_eq!(result.conductive_enthalpy_increment_j[1], 0.0);
    close(returned.conductive_enthalpy_increment_j, 3.0);
    close(report.solid_to_fluid_heat_j, 6.0);
    close(report.solid_authored_heat_j, 5.0);
    close(report.solid_enthalpy_change_j, 8.0);
    close(report.solid_boundary_heat_j, 9.0);
    close(report.absolute_solid_boundary_heat_j, 15.0);
    close(report.conduction.absolute_boundary_heat_j, 16.0);
    for (actual, expected) in report
        .conduction
        .boundary_heat_j
        .into_iter()
        .zip([1.0, 12.0, -3.0, 0.0, 0.0, 0.0])
    {
        close(actual, expected);
    }
    close(
        report.solid_enthalpy_change_j,
        report.solid_authored_heat_j + report.solid_boundary_heat_j - report.solid_to_fluid_heat_j,
    );
    close(
        result.thermal.volume.report.material_enthalpy_change_j + report.solid_enthalpy_change_j,
        15.0,
    );
    close(report.conduction.solve.energy.external_heat_j, 15.0);
    assert_eq!(report.conduction.boundary_faces, 3);
    assert_eq!(report.conduction.solve.energy.substeps, 1);
    assert!(report.conduction.absolute_enthalpy_defect_j <= 1e-9);
    assert_eq!(returned.conductivity_source_identity, solid.identity());
    assert_eq!(
        returned.phase.phase_curve_identity(),
        solid_curve.identity()
    );
    assert_eq!(returned.mass_kg, 4.0);
    assert_eq!(result.thermal.volume.transport.mass_kg, mass);
    assert!(result.thermal.volume.phase_states[1].is_none());
    assert_eq!(result.thermal.volume.transport.report.mass_routes, 0);

    let momentum = [[0.0; 3]; 2];
    let coupled = field
        .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            None,
            &fluid,
            &baths,
            &walls,
            &contacts,
            ThermalVolumeMomentumInput {
                volume: data,
                momentum_kg_m_s: &momentum,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &momentum,
                boundary_velocity_m_s: [None; 6],
            },
            super::super::super::momentum::coupled_controls(1.0),
            solve_controls(),
            &mut workspace,
            &gate,
        )
        .unwrap();
    assert_eq!(coupled.evolution.thermal, result.thermal);
    assert_eq!(coupled.solid_endpoints, result.solid_endpoints);
    assert_eq!(coupled.report, result.report);
    assert_eq!(coupled.evolution.momentum.momentum_kg_m_s, momentum);

    // Source bounds remain those of the WALL, even though every fluid state
    // is admitted. A failed wall solve leaves both inventories available to retry.
    let bounded = source(&solid_card, &solid_curve, [47.0, 49.0]);
    assert!(
        field
            .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
                &domain,
                [FlowFace::Wall; 6],
                &fluid,
                &baths,
                &[ThermalVolumeSolidCell {
                    conductivity: &bounded,
                    ..walls[0]
                }],
                &contacts,
                data,
                controls(1.0),
                solve_controls(),
                &mut workspace,
                &gate,
            )
            .is_err()
    );
    let retry = field
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            &baths,
            &walls,
            &contacts,
            data,
            controls(1.0),
            solve_controls(),
            &mut workspace,
            &gate,
        )
        .unwrap();
    assert_eq!(retry, result);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 9.5);
    assert_eq!(walls[0].phase.specific_enthalpy_j_kg(), 48.0);
}

#[test]
fn solid_baths_require_exposed_bound_receivers_and_preserve_fluid_only_admission() {
    let material = card("T", &[(300.0, 1.0), (400.0, 1.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (100.0, 400.0, 1.0, 1.0)],
    );
    let conductivity = source(&material, &caloric, [0.0, 100.0]);
    let mut domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    domain.set(2, 0, 0, Voxel::Solid(0));
    domain.set(3, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let phases = [Some(caloric.state_at_specific_enthalpy(10.0).unwrap()); 2];
    let states = [phases[0], phases[1], None, None];
    let masses = [1.0, 1.0, 0.0, 0.0];
    let zero = [0.0; 4];
    let data = input(&masses, &states, &zero);
    let walls = [2, 3].map(|cell| ThermalVolumeSolidCell {
        cell,
        conductivity: &conductivity,
        mass_kg: 1.0,
        phase: caloric.state_at_specific_enthalpy(20.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    });
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(4).unwrap();
    let run = |baths: &[ThermalVolumeBathFace],
               walls: &[ThermalVolumeSolidCell<'_, '_>],
               workspace: &mut _| {
        field.advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            baths,
            walls,
            &[],
            data,
            controls(1.0),
            solve_controls(),
            workspace,
            &gate,
        )
    };
    // Neither orientation of a fluid/solid or solid/solid occupied face can
    // acquire a second bath channel. Fluid/fluid faces keep the same rule.
    for (cell, side) in [
        (0, Face3::XMax),
        (1, Face3::XMax),
        (2, Face3::XMin),
        (2, Face3::XMax),
        (3, Face3::XMin),
    ] {
        assert!(
            matches!(run(&[bath(cell, side, 350.0, 0.0)], &walls, &mut workspace),
            Err(ChtError::InvalidInput { reason, .. }) if reason.contains("boundary fluid-cell face"))
        );
    }
    let exposed = bath(3, Face3::XMax, 350.0, 0.5);
    assert!(matches!(run(&[exposed, exposed], &walls, &mut workspace),
        Err(ChtError::InvalidInput { reason, .. }) if reason.contains("duplicate thermal bath")));
    assert!(matches!(run(&[exposed], &walls[..1], &mut workspace),
        Err(ChtError::InvalidInput { reason, .. }) if reason.contains("boundary fluid-cell face")));
    // Original fluid-only APIs do not silently activate solid storage.
    assert!(
        matches!(field.advance_isotropic_implicit_conductive_thermal_volume_with_baths(
        &domain, [FlowFace::Wall; 6], &conductivity, &[exposed], data, controls(1.0), solve_controls(),
        &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(), &gate,
    ), Err(ChtError::InvalidInput { reason, .. }) if reason.contains("boundary fluid-cell face"))
    );
    let accepted = run(&[exposed], &walls, &mut workspace).unwrap();
    assert!(accepted.report.solid_boundary_heat_j > 0.0);
    assert_eq!(accepted, run(&[exposed], &walls, &mut workspace).unwrap());
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 10.0);
    assert_eq!(walls[1].phase.specific_enthalpy_j_kg(), 20.0);
}

#[test]
fn equal_content_caloric_curves_keep_their_original_bath_owner() {
    let material = card("T", &[(300.0, 1.0), (400.0, 1.0)]);
    let knots = [(0.0, 300.0, 1.0, 1.0), (100.0, 400.0, 1.0, 1.0)];
    let fluid_curve = curve(&material, &knots);
    let wall_curve = curve(&material, &knots);
    assert_eq!(fluid_curve.identity(), wall_curve.identity());
    assert!(!core::ptr::eq(&fluid_curve, &wall_curve));
    let fluid = source(&material, &fluid_curve, [0.0, 100.0]);
    let wall_source = source(&material, &wall_curve, [0.0, 100.0]);
    assert_eq!(fluid.identity(), wall_source.identity());
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let phases = [
        Some(fluid_curve.state_at_specific_enthalpy(10.0).unwrap()),
        None,
    ];
    let wall = ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall_source,
        mass_kg: 1.0,
        phase: wall_curve.state_at_specific_enthalpy(20.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    };
    let result = field
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            &[bath(1, Face3::XMax, 330.0, 0.5)],
            &[wall],
            &[],
            input(&[1.0, 0.0], &phases, &[0.0; 2]),
            controls(1.0),
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(2).unwrap(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // M=Cp=G_interface=G_bath=dt=1: 2*h_f-h_s=10 and -h_f+3*h_s=50.
    // The wall bath must retain its original caloric object, even though all
    // immutable material/curve/source content hashes equal the fluid's.
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .specific_enthalpy_j_kg(),
        16.0,
    );
    close(
        result.solid_endpoints[0].phase.specific_enthalpy_j_kg(),
        22.0,
    );
    close(result.report.solid_boundary_heat_j, 8.0);
    close(result.report.solid_to_fluid_heat_j, 6.0);
    close(result.report.solid_enthalpy_change_j, 2.0);
    close(result.report.conduction.solve.energy.external_heat_j, 8.0);
}
