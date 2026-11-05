//! G1/G3/G4: physical half-cell resistance, latent energy and atomic refusal.

use super::*;
use crate::d3q19::Face3;

fn bath(cell: usize, face: Face3, temperature_k: f64, resistance: f64) -> ThermalVolumeBathFace {
    ThermalVolumeBathFace {
        cell,
        face,
        temperature_k,
        resistance_m2_k_w: resistance,
    }
}

#[test]
fn six_graded_faces_follow_independent_fourier_resistance_and_calorimetry() {
    // rho=1, Cp=2, k(T)=2+0.04*(T-300). At the original 350 K, k=4.
    let material = card("T", &[(300.0, 2.0), (500.0, 10.0)]);
    let domain = VoxelDomain::graded([vec![2.0], vec![3.0], vec![4.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let gate = CancelGate::new_clock_free();
    for datum in [0.0, 1024.0] {
        let caloric = curve(
            &material,
            &[(datum, 300.0, 1.0, 1.0), (datum + 400.0, 500.0, 1.0, 1.0)],
        );
        let source = source(&material, &caloric, [datum, datum + 400.0]);
        let states = [Some(
            caloric.state_at_specific_enthalpy(datum + 100.0).unwrap(),
        )];
        let data = input(&[24.0], &states, &[0.0]);
        let legacy = field
            .advance_isotropic_conductive_thermal_volume(
                &domain,
                [FlowFace::Wall; 6],
                &source,
                data,
                controls(0.125),
                &gate,
            )
            .unwrap();
        let empty = field
            .advance_isotropic_conductive_thermal_volume_with_baths(
                &domain,
                [FlowFace::Wall; 6],
                &source,
                &[],
                data,
                controls(0.125),
                &gate,
            )
            .unwrap();
        assert_eq!(legacy.report, empty.report);
        assert_eq!(
            legacy.conductive_enthalpy_increment_j,
            empty.conductive_enthalpy_increment_j
        );
        assert_eq!(
            legacy.thermal.volume.transport.specific_enthalpy_j_kg,
            empty.thermal.volume.transport.specific_enthalpy_j_kg
        );
        for side in Face3::ALL {
            let axis = side as usize / 2;
            let distance = [1.0, 1.5, 2.0][axis];
            let area = [12.0, 8.0, 6.0][axis];
            for resistance in [0.0, 0.5] {
                for temperature in [340.0, 360.0] {
                    let heat = 0.125 * area * (temperature - 350.0) / (distance / 4.0 + resistance);
                    let result = field
                        .advance_isotropic_conductive_thermal_volume_with_baths(
                            &domain,
                            [FlowFace::Wall; 6],
                            &source,
                            &[bath(0, side, temperature, resistance)],
                            data,
                            controls(0.125),
                            &gate,
                        )
                        .unwrap();
                    close(result.conductive_enthalpy_increment_j[0], heat);
                    close(
                        result.thermal.volume.transport.specific_enthalpy_j_kg[0] - datum,
                        100.0 + heat / 24.0,
                    );
                    close(
                        result.thermal.volume.phase_states[0]
                            .unwrap()
                            .temperature_k(),
                        350.0 + heat / 48.0,
                    );
                    close(result.thermal.volume.transport.mass_kg[0], 24.0);
                    close(
                        result.thermal.volume.report.material_enthalpy_change_j,
                        heat,
                    );
                    close(
                        result.thermal.volume.report.requested_enthalpy_increment_j,
                        heat,
                    );
                    close(result.report.conduction_balance_defect_j, 0.0);
                    close(result.report.absolute_boundary_heat_j, heat.abs());
                    close(
                        result.report.maximum_enthalpy_mixing_weight,
                        0.125 * area / (distance / 4.0 + resistance) / 48.0,
                    );
                    assert_eq!(result.report.boundary_faces, 1);
                    assert_eq!(result.report.conductive_faces, 0);
                    for face in 0..6 {
                        close(
                            result.report.boundary_heat_j[face],
                            if face == side as usize { heat } else { 0.0 },
                        );
                    }
                    assert_eq!(
                        result.report.conductivity_source_identity,
                        source.identity()
                    );
                }
            }
        }
    }
    assert_eq!(field.flux_x(0, 0, 0), 0.0);
    assert_eq!(field.flux_x(1, 0, 0), 0.0);
}

#[test]
fn opposed_baths_retain_gross_heat_and_are_independent_of_input_order() {
    let material = card("T", &[(300.0, 4.0), (500.0, 4.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (400.0, 500.0, 1.0, 1.0)],
    );
    let source = source(&material, &caloric, [0.0, 400.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [Some(caloric.state_at_specific_enthalpy(100.0).unwrap()); 2];
    let boundaries = [
        bath(0, Face3::XMin, 360.0, 0.25),
        bath(1, Face3::XMax, 340.0, 0.25),
    ];
    let run = |baths: &[ThermalVolumeBathFace]| {
        field
            .advance_isotropic_conductive_thermal_volume_with_baths(
                &domain,
                [FlowFace::Wall; 6],
                &source,
                baths,
                input(&[1.0; 2], &states, &[0.0; 2]),
                controls(0.125),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let result = run(&boundaries);
    assert_eq!(result, run(&[boundaries[1], boundaries[0]]));
    let heat = 10.0 / 3.0;
    close(result.report.boundary_heat_j[0], heat);
    close(result.report.boundary_heat_j[1], -heat);
    close(result.report.absolute_boundary_heat_j, 2.0 * heat);
    close(result.report.exchanged_heat_j, 0.0);
    close(result.thermal.volume.report.material_enthalpy_change_j, 0.0);
    for (cell, sign) in [(0, 1.0), (1, -1.0)] {
        close(
            result.thermal.volume.transport.specific_enthalpy_j_kg[cell],
            100.0 + sign * heat,
        );
    }
}

#[test]
fn finite_bath_heat_changes_latent_fraction_without_fictitious_heat_capacity() {
    let material = card("f", &[(0.0, 2.0), (1.0, 6.0)]);
    let caloric = curve(
        &material,
        &[
            (-100.0, 300.0, 0.0, 1.0),
            (0.0, 400.0, 0.0, 1.0),
            (100.0, 400.0, 1.0, 1.0),
            (200.0, 500.0, 1.0, 1.0),
        ],
    );
    // Only the true latent plateau is admitted by this compiled source.
    let source = source(&material, &caloric, [0.0, 100.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [25.0, 75.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let data = input(&[1.0; 2], &states, &[0.0; 2]);
    let baths = [
        bath(0, Face3::XMin, 410.0, 0.25),
        bath(1, Face3::XMax, 390.0, 0.25),
    ];
    let result = field
        .advance_isotropic_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &baths,
            data,
            controls(1.0),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // k=3 and 5; two independent half cells in series with R''=1/4.
    for (cell, h) in [(0, 49.0), (1, 75.0 - 200.0 / 7.0)] {
        let state = result.thermal.volume.phase_states[cell].unwrap();
        close(state.specific_enthalpy_j_kg(), h);
        close(state.liquid_mass_fraction(), h / 100.0);
        assert_eq!(state.temperature_k(), 400.0);
        close(result.thermal.volume.transport.mass_kg[cell], 1.0);
    }
    close(result.report.maximum_enthalpy_mixing_weight, 0.0);
    close(result.report.exchanged_heat_j, 0.0);
    close(
        result.thermal.volume.report.material_enthalpy_change_j,
        24.0 - 200.0 / 7.0,
    );
    invalid(
        field.advance_isotropic_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &baths,
            data,
            controls(100.0),
            &CancelGate::new_clock_free(),
        ),
        "calorically advanced state",
    );
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 25.0);
}

#[test]
fn invalid_baths_and_combined_heat_limits_leave_original_inventories_unchanged() {
    let material = card("T", &[(300.0, 4.0), (500.0, 4.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 1.0), (400.0, 500.0, 1.0, 1.0)],
    );
    let source = source(&material, &caloric, [0.0, 400.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [Some(caloric.state_at_specific_enthalpy(100.0).unwrap()); 2];
    let data = input(&[1.0; 2], &states, &[0.0; 2]);
    let good = bath(0, Face3::XMin, 360.0, 0.0);
    let gate = CancelGate::new_clock_free();
    let run = |baths: &[ThermalVolumeBathFace], controls| {
        field.advance_isotropic_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            baths,
            data,
            controls,
            &gate,
        )
    };
    for bad in [
        bath(0, Face3::XMax, 360.0, 0.0),
        bath(2, Face3::XMax, 360.0, 0.0),
        bath(0, Face3::XMin, f64::NAN, 0.0),
        bath(0, Face3::XMin, 0.0, 0.0),
        bath(0, Face3::XMin, 360.0, -1.0),
    ] {
        assert!(run(&[bad], controls(0.125)).is_err());
    }
    invalid(
        run(&[good, good], controls(0.125)),
        "duplicate thermal bath",
    );
    let mut capped = controls(0.125);
    capped.max_conductive_faces = 1; // one internal face PLUS one bath needs two
    invalid(run(&[good], capped), "face-work budget");
    invalid(run(&[good], controls(1.0)), "mixing weight exceeds one");
    let unequal = [100.0, 200.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    // Internal weight .4 and bath weight .8 are individually admissible,
    // but the same cell must reject their combined explicit weight 1.2.
    invalid(
        field.advance_isotropic_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &[good],
            input(&[1.0; 2], &unequal, &[0.0; 2]),
            controls(0.2),
            &gate,
        ),
        "mixing weight exceeds one",
    );
    let mut solid_domain = domain.clone();
    solid_domain.set(0, 0, 0, Voxel::Solid(0));
    invalid(
        FlowField::quiescent(&solid_domain).advance_isotropic_conductive_thermal_volume_with_baths(
            &solid_domain,
            [FlowFace::Wall; 6],
            &source,
            &[good],
            input(&[0.0, 1.0], &[None, states[1]], &[0.0; 2]),
            controls(0.125),
            &gate,
        ),
        "exterior fluid-cell face",
    );
    let accepted = run(&[good], controls(0.125)).unwrap();
    gate.request();
    assert!(matches!(
        run(&[good], controls(0.125)),
        Err(ChtError::Cancelled)
    ));
    let replay = field
        .advance_isotropic_conductive_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &source,
            &[good],
            data,
            controls(0.125),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(accepted, replay);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(states[0].unwrap().specific_enthalpy_j_kg(), 100.0);
}
