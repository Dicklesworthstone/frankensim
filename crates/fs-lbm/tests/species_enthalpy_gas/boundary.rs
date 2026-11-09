//! Original gas calorics on the native fixed-composition surface map.
use super::*;
use fs_lbm::conjugate::{GasSpeciesEnthalpyBoundarySource, SpeciesDiffusionBoundary};

#[test]
fn explicit_physical_face_temperature_closes_open_energy_with_formation_references() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let accepted = carrier(&domain, &[1.0]);
    let masses = [[0.8, 0.2]];
    let boundaries = [SpeciesDiffusionBoundary {
        cell: 0,
        face: 0,
        reservoir_mass_fractions: [0.25, 0.75],
        film_resistance_m2_s_kg: 2.5,
    }];
    let face_t = [450.0];
    let names = ids();
    let context = source::point(350.0, 100_000.0);
    let diffusivity = [0.2];
    let gate = CancelGate::new_clock_free();
    // Unit area/mass/volume: half-cell resistance 0.5/(rho*D)=2.5,
    // in series with film resistance 2.5, hence G=0.2 kg/s.
    let conductance: f64 = 0.2;
    let dm = (0.25 - 0.8) * (-(-conductance).exp_m1());
    let expected_mass = [0.8 + dm, 0.2 - dm];
    let mut reference_temperature: Option<f64> = None;
    for a7 in [0.0, -4000.0, 10000.0] {
        let bank = source::moist_model(0.5, a7);
        let calorics =
            GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
                .unwrap();
        let bound =
            GasSpeciesEnthalpyBoundarySource::new(&calorics, &boundaries, &face_t, &gate).unwrap();
        let h = [enthalpy(&bank, masses[0], 350.0)];
        let before = (accepted.clone(), masses, h);
        let result = accepted
            .diffuse_thermochemical_species_with_boundaries(
                &domain,
                &masses,
                &h,
                1.0,
                &bound,
                controls(),
                ENERGY_BUDGET,
                &gate,
            )
            .unwrap();
        let face_h = partial_h(face_t[0], a7);
        let expected_heat = dm * (face_h[0] - face_h[1]);
        let expected_t = temperature(h[0] + expected_heat, expected_mass, a7);
        let mass_report = &result.step.boundary.faces[0];
        near(mass_report.conductance_kg_s, conductance, 1e-15);
        for k in 0..2 {
            near(
                result.step.step.species.species_mass_kg[0][k],
                expected_mass[k],
                5e-14,
            );
            near(mass_report.mass_into_cell_kg[k], [dm, -dm][k], 5e-14);
        }
        near(
            result.step.thermal_boundary.enthalpy_into_domain_j,
            expected_heat,
            2e-7,
        );
        near(
            result.step.thermal_boundary.faces[0].enthalpy_into_cell_j,
            expected_heat,
            2e-7,
        );
        near(
            result.step.thermal_boundary.positive_enthalpy_into_domain_j,
            expected_heat.max(0.0),
            2e-7,
        );
        near(
            result.step.thermal_boundary.negative_enthalpy_into_domain_j,
            (-expected_heat).max(0.0),
            2e-7,
        );
        near(
            result.step.step.diffusive_enthalpy_change_j[0],
            expected_heat,
            2e-7,
        );
        near(
            result.step.step.specific_enthalpy_j_kg[0],
            h[0] + expected_heat,
            3e-7,
        );
        near(result.state.temperatures_k()[0], expected_t, 1e-8);
        near(
            result.state.transport(0).unwrap().specific_heat_j_kg_k(),
            capacity(expected_mass) * (3.0 + 0.001 * expected_t),
            2e-8,
        );
        // A cell/donor-temperature substitution would incorrectly preserve 350 K.
        assert!((expected_t - 350.0).abs() > 1.0);
        assert!(expected_heat.abs() > 1000.0);
        if let Some(previous) = reference_temperature {
            near(result.state.temperatures_k()[0], previous, 1e-8);
        } else {
            reference_temperature = Some(result.state.temperatures_k()[0]);
        }
        assert!(result.report.absolute_enthalpy_defect_j <= ENERGY_BUDGET);
        assert_eq!(result.state.species(), &names);
        assert_eq!(result.state.source().source_packs(), bank.source_packs());
        assert_eq!(accepted, before.0);
        assert_eq!(masses, before.1);
        assert_eq!(h, before.2);
    }
}

#[test]
fn absent_species_uses_full_bank_and_zero_diffusivity_blocks_boundary_heat() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let accepted = carrier(&domain, &[1.0]);
    let masses = [[1.0, 0.0]];
    let bank = source::moist_model(0.5, -4000.0);
    let h = [enthalpy(&bank, masses[0], 400.0)];
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let boundaries = [SpeciesDiffusionBoundary {
        cell: 0,
        face: 1,
        reservoir_mass_fractions: [0.0, 1.0],
        film_resistance_m2_s_kg: 0.0,
    }];
    let face_t = [400.0];
    let gate = CancelGate::new_clock_free();
    let before = (accepted.clone(), masses, h);
    for d in [0.0, 0.2] {
        let diffusivity = [d];
        let calorics =
            GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
                .unwrap();
        let bound =
            GasSpeciesEnthalpyBoundarySource::new(&calorics, &boundaries, &face_t, &gate).unwrap();
        let field = bound
            .prepare(&domain, &accepted, &masses, &h, &gate)
            .unwrap();
        let partials = bound.boundary_partials(&field, &domain, 0, &gate).unwrap();
        for k in 0..2 {
            near(partials[k], partial_h(400.0, -4000.0)[k], 2e-7);
        }
        assert!(bound.boundary_partials(&field, &domain, 1, &gate).is_err());
        let result = accepted
            .diffuse_thermochemical_species_with_boundaries(
                &domain,
                &masses,
                &h,
                0.5,
                &bound,
                controls(),
                ENERGY_BUDGET,
                &gate,
            )
            .unwrap();
        if d == 0.0 {
            assert_eq!(result.step.step.species.species_mass_kg, masses);
            assert_eq!(result.step.step.specific_enthalpy_j_kg, h);
            assert_eq!(result.step.thermal_boundary.enthalpy_into_domain_j, 0.0);
        } else {
            assert!(result.step.step.species.species_mass_kg[0][1] > 0.0);
            assert_ne!(result.step.thermal_boundary.enthalpy_into_domain_j, 0.0);
        }
        near(result.state.temperatures_k()[0], 400.0, 1e-8);
        assert_eq!(result.state.species(), &names);
        assert!(
            GasSpeciesEnthalpyBoundarySource::new(&calorics, &boundaries, &[], &gate,).is_err()
        );
        for bad_t in [f64::NAN, 199.0, 701.0] {
            assert!(
                GasSpeciesEnthalpyBoundarySource::new(&calorics, &boundaries, &[bad_t], &gate,)
                    .is_err()
            );
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            accepted.diffuse_thermochemical_species_with_boundaries(
                &domain,
                &masses,
                &h,
                0.5,
                &bound,
                controls(),
                ENERGY_BUDGET,
                &cancelled,
            ),
            Err(ChtError::Cancelled)
        ));
    }
    assert_eq!(accepted, before.0);
    assert_eq!(masses, before.1);
    assert_eq!(h, before.2);
}
