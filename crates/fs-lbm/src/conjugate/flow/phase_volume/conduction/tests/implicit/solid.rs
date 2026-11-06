//! Finite-wall calorimetry, native momentum and exclusive contact ownership.

mod baths;

use super::*;
use crate::conjugate::FaceContactResistance;

#[test]
fn graded_fluid_pcm_contact_uses_both_sources_and_conserves_finite_storage() {
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
    // Manufactured Q*=10: wall h*=25+10/4=27.5, f*=.275 and k*=1.275.
    // G*=1/(.5/2 + R + 1/1.275), R=.25. The fluid h*=h0-10/2;
    // Tfluid*=320+10/G*. Thus h0=30+400/51 and the BE pair closes exactly.
    let h0 = 30.0 + 400.0 / 51.0;
    let phases = [
        Some(fluid_curve.state_at_specific_enthalpy(h0).unwrap()),
        None,
    ];
    let zero = [0.0; 2];
    let data = input(&mass, &phases, &zero);
    let walls = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &solid,
        mass_kg: 4.0,
        phase: solid_curve.state_at_specific_enthalpy(25.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let contacts = [ThermalVolumeContact::Resistance(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.25,
    })];
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let result = field
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            ThermalVolumeBoundaries::default(),
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
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .specific_enthalpy_j_kg(),
        h0 - 5.0,
    );
    close(returned.phase.specific_enthalpy_j_kg(), 27.5);
    close(returned.phase.temperature_k(), 320.0);
    close(returned.phase.liquid_mass_fraction(), 0.275);
    close(returned.conductive_enthalpy_increment_j, 10.0);
    close(result.conductive_enthalpy_increment_j[0], -10.0);
    assert_eq!(result.conductive_enthalpy_increment_j[1], 0.0);
    assert_eq!(returned.mass_kg, 4.0);
    assert_eq!(returned.conductivity_source_identity, solid.identity());
    assert_eq!(
        returned.phase.phase_curve_identity(),
        solid_curve.identity()
    );
    assert_eq!(result.thermal.volume.transport.mass_kg, mass);
    assert!(result.thermal.volume.phase_states[1].is_none());
    assert_eq!(result.thermal.volume.transport.report.mass_routes, 0);
    close(result.report.solid_enthalpy_change_j, 10.0);
    close(result.report.solid_to_fluid_heat_j, -10.0);
    close(result.report.solid_authored_heat_j, 0.0);
    close(
        result.thermal.volume.report.material_enthalpy_change_j
            + result.report.solid_enthalpy_change_j,
        0.0,
    );
    close(result.report.conduction.solve.energy.external_heat_j, 0.0);
    assert_eq!(result.report.conduction.solve.energy.substeps, 1);
    assert_eq!(result.report.conduction.conductive_faces, 1);
    assert!(result.report.conduction.absolute_enthalpy_defect_j <= 1e-9);

    // The pressure/momentum adapter consumes the same heat result once. No
    // thermal wall inventory or fictitious momentum enters the fluid arrays.
    let momentum = [[0.0; 3]; 2];
    let coupled = field
        .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            None,
            &fluid,
            ThermalVolumeBoundaries::default(),
            &walls,
            &contacts,
            ThermalVolumeMomentumInput {
                volume: data,
                momentum_kg_m_s: &momentum,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &momentum,
                boundary_velocity_m_s: [None; 6],
            },
            super::super::momentum::coupled_controls(1.0),
            solve_controls(),
            &mut workspace,
            &gate,
        )
        .unwrap();
    assert_eq!(coupled.evolution.thermal, result.thermal);
    assert_eq!(coupled.solid_endpoints, result.solid_endpoints);
    assert_eq!(coupled.report, result.report);
    assert_eq!(coupled.evolution.momentum.momentum_kg_m_s, momentum);
    let perfect = field
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            ThermalVolumeBoundaries::default(),
            &walls,
            &[],
            data,
            controls(1.0),
            solve_controls(),
            &mut workspace,
            &gate,
        )
        .unwrap();
    assert!(perfect.solid_endpoints[0].conductive_enthalpy_increment_j > 10.0);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(walls[0].phase.specific_enthalpy_j_kg(), 25.0);
}

#[test]
fn authored_heat_crosses_solid_solid_and_fluid_faces_without_solid_advection() {
    let material = card("T", &[(300.0, 1.0), (400.0, 1.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (100.0, 400.0, 1.0, 1.0)],
    );
    let conductivity = source(&material, &caloric, [0.0, 100.0]);
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    domain.set(2, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let phases = [
        Some(caloric.state_at_specific_enthalpy(10.0).unwrap()),
        None,
        None,
    ];
    let wall = |cell, h, heat| ThermalVolumeSolidCell {
        cell,
        conductivity: &conductivity,
        mass_kg: 1.0,
        phase: caloric.state_at_specific_enthalpy(h).unwrap(),
        net_enthalpy_increment_j: heat,
    };
    let walls = [wall(2, 30.0, 4.0), wall(1, 20.0, 0.0)];
    let data = input(&[1.0, 0.0, 0.0], &phases, &[0.0; 3]);
    let gate = CancelGate::new_clock_free();
    let result = field
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            ThermalVolumeBoundaries::default(),
            &walls,
            &[],
            data,
            controls(1.0),
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(3).unwrap(),
            &gate,
        )
        .unwrap();
    // M=Cp=G=dt=1. Independent equations are
    // 2h0-h1=10; -h0+3h1-h2=20; -h1+2h2=34.
    // Hence h=[15.5,21,27.5], face transfers to the fluid/left solid=[5.5,6.5].
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .specific_enthalpy_j_kg(),
        15.5,
    );
    assert_eq!(
        result
            .solid_endpoints
            .iter()
            .map(|s| s.cell)
            .collect::<Vec<_>>(),
        [1, 2]
    );
    close(
        result.solid_endpoints[0].phase.specific_enthalpy_j_kg(),
        21.0,
    );
    close(
        result.solid_endpoints[1].phase.specific_enthalpy_j_kg(),
        27.5,
    );
    close(
        result.solid_endpoints[0].conductive_enthalpy_increment_j,
        1.0,
    );
    close(
        result.solid_endpoints[1].conductive_enthalpy_increment_j,
        -6.5,
    );
    close(result.report.solid_to_fluid_heat_j, 5.5);
    close(result.report.solid_enthalpy_change_j, -1.5);
    close(result.report.solid_authored_heat_j, 4.0);
    close(result.report.conduction.solve.energy.external_heat_j, 4.0);
    close(
        result.thermal.volume.report.material_enthalpy_change_j
            + result.report.solid_enthalpy_change_j,
        4.0,
    );
    assert_eq!(result.thermal.volume.transport.mass_kg, [1.0, 0.0, 0.0]);
    assert_eq!(result.report.conduction.conductive_faces, 2);

    // Empty walls preserve the existing insulated-solid path bit for bit.
    let legacy = field
        .advance_isotropic_implicit_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            &[],
            data,
            controls(1.0),
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(1).unwrap(),
            &gate,
        )
        .unwrap();
    let empty = field
        .advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            ThermalVolumeBoundaries::default(),
            &[],
            &[],
            data,
            controls(1.0),
            solve_controls(),
            &mut ImplicitThermalVolumeWorkspace::new(1).unwrap(),
            &gate,
        )
        .unwrap();
    assert_eq!(empty.thermal, legacy.thermal);
    assert_eq!(empty.report.conduction, legacy.report);
    assert_eq!(
        empty.conductive_enthalpy_increment_j,
        legacy.conductive_enthalpy_increment_j
    );
    assert!(empty.solid_endpoints.is_empty());
}

#[test]
fn finite_wall_ownership_source_and_contact_refusals_are_atomic() {
    let material = card("T", &[(300.0, 1.0), (400.0, 1.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (100.0, 400.0, 1.0, 1.0)],
    );
    let conductivity = source(&material, &caloric, [0.0, 100.0]);
    let bounded = source(&material, &caloric, [49.0, 50.0]);
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    domain.set(2, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let phases = [
        Some(caloric.state_at_specific_enthalpy(10.0).unwrap()),
        None,
        None,
    ];
    let wall = ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &conductivity,
        mass_kg: 1.0,
        phase: caloric.state_at_specific_enthalpy(50.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let run = |walls: &[ThermalVolumeSolidCell<'_, '_>],
               contacts: &[ThermalVolumeContact<'_>],
               baths: &[ThermalVolumeBathFace],
               workspace: &mut _,
               gate: &_| {
        field.advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            baths.into(),
            walls,
            contacts,
            input(&[1.0, 0.0, 0.0], &phases, &[0.0; 3]),
            controls(1.0),
            solve_controls(),
            workspace,
            gate,
        )
    };
    let gate = CancelGate::new_clock_free();
    let face = FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.5,
    };
    for contacts in [
        vec![
            face,
            FaceContactResistance {
                cells: [1, 0],
                ..face
            },
        ],
        vec![FaceContactResistance {
            cells: [0, 2],
            ..face
        }],
        vec![FaceContactResistance {
            cells: [1, 2],
            ..face
        }],
        vec![FaceContactResistance {
            resistance_m2_k_w: -1.0,
            ..face
        }],
    ] {
        let contacts: Vec<_> = contacts.into_iter().map(ThermalVolumeContact::from).collect();
        assert!(run(&[wall], &contacts, &[], &mut workspace, &gate).is_err());
    }
    // Both endpoints are active and workspace is sufficient: this refusal
    // specifically prevents a nonadjacent contact from being silently ignored.
    assert!(
        run(
            &[wall, ThermalVolumeSolidCell { cell: 2, ..wall }],
            &[FaceContactResistance {
                cells: [0, 2],
                ..face
            }
            .into()],
            &[],
            &mut ImplicitThermalVolumeWorkspace::new(3).unwrap(),
            &gate,
        )
        .is_err()
    );
    assert!(run(&[wall, wall], &[], &[], &mut workspace, &gate).is_err());
    assert!(
        run(
            &[ThermalVolumeSolidCell { cell: 0, ..wall }],
            &[],
            &[],
            &mut workspace,
            &gate
        )
        .is_err()
    );
    assert!(
        run(
            &[wall],
            &[],
            &[bath(0, Face3::XMax, 400.0, 0.0)],
            &mut workspace,
            &gate
        )
        .is_err()
    );
    assert!(
        run(
            &[ThermalVolumeSolidCell {
                conductivity: &bounded,
                ..wall
            }],
            &[],
            &[],
            &mut workspace,
            &gate
        )
        .is_err()
    );
    gate.request();
    assert!(matches!(
        run(&[wall], &[], &[], &mut workspace, &gate),
        Err(ChtError::Cancelled)
    ));
    let fresh = CancelGate::new_clock_free();
    let accepted = run(&[wall], &[face.into()], &[], &mut workspace, &fresh).unwrap();
    assert_eq!(
        accepted,
        run(&[wall], &[face.into()], &[], &mut workspace, &fresh).unwrap()
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(wall.phase.specific_enthalpy_j_kg(), 50.0);
}
