use super::*;
use crate::phase::{HeatCapacityKnot, SolidLiquidPhase};

fn knot(h: f64, t: f64, f: f64, rho: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: rho,
    }
}

fn matrix(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([1; 32]),
        SolidLiquidPhase::Solid,
        vec![
            knot(datum, 300.0, 0.0, 8000.0),
            knot(datum + 200_000.0, 500.0, 0.0, 8000.0),
        ],
    )
    .unwrap()
}

fn pcm(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([2; 32]),
        vec![
            knot(datum, 300.0, 0.0, 1000.0),
            knot(datum + 200_000.0, 400.0, 0.0, 1000.0),
            knot(datum + 400_000.0, 400.0, 1.0, 800.0),
            knot(datum + 650_000.0, 500.0, 1.0, 800.0),
        ],
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn original_affine_capacities_compose_without_linearizing_either_constituent() {
    let matrix = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([0x91; 32]),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 500.0,
                bulk_density_kg_m3: 8000.0,
            },
            HeatCapacityKnot {
                temperature_k: 500.0,
                specific_heat_capacity_j_kg_k: 1500.0,
                bulk_density_kg_m3: 6000.0,
            },
        ],
        1e6,
    )
    .unwrap();
    let mut pcm = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x92; 32]),
        vec![
            knot(0.0, 300.0, 0.0, 1000.0),
            knot(150_000.0, 400.0, 0.0, 900.0),
            knot(350_000.0, 400.0, 1.0, 800.0),
            knot(650_000.0, 500.0, 1.0, 700.0),
        ],
    )
    .unwrap();
    pcm.retain_sensible_density_law();
    pcm.retain_sensible_heat_capacities(vec![Some([1000.0, 2000.0]), None, Some([2000.0, 4000.0])])
        .unwrap();
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    assert_eq!(law.knots().len(), 4);
    close(law.maximum_temperature_slope(), 1.0 / 625.0, 1e-18);
    for (h, t, cm, cp, rho_m, rho_p) in [
        (39_062.5, 350.0, 750.0, 1500.0, 7500.0, 950.0),
        (217_187.5, 450.0, 1250.0, 3000.0, 6500.0, 750.0),
    ] {
        let response = law.response_at_specific_enthalpy(h).unwrap();
        let state = response.state;
        close(state.temperature_k, t, 2e-12);
        close(state.matrix.temperature_k, t, 2e-12);
        close(state.pcm.temperature_k, t, 2e-12);
        close(state.matrix.bulk_density_kg_m3, rho_m, 2e-11);
        close(state.pcm.bulk_density_kg_m3, rho_p, 2e-12);
        close(
            0.75 * state.matrix.specific_enthalpy_j_kg + 0.25 * state.pcm.specific_enthalpy_j_kg,
            h,
            2e-9,
        );
        let capacity = 0.75 * cm + 0.25 * cp;
        close(response.temperature_slope_kg_k_j, 1.0 / capacity, 2e-18);
        close(response.matrix_enthalpy_derivative, cm / capacity, 2e-15);
        close(response.pcm_enthalpy_derivative, cp / capacity, 2e-15);
        let difference = (law
            .values_at_specific_enthalpy(h + 0.01)
            .unwrap()
            .pcm_volume_fraction
            - law
                .values_at_specific_enthalpy(h - 0.01)
                .unwrap()
                .pcm_volume_fraction)
            / 0.02;
        close(response.pcm_volume_fraction_slope_kg_j, difference, 2e-13);
    }
    let melting = law.response_at_specific_enthalpy(118_750.0).unwrap();
    assert_eq!(melting.state.temperature_k, 400.0);
    assert_eq!(melting.state.pcm.liquid_mass_fraction, 0.5);
    assert_eq!(melting.temperature_slope_kg_k_j, 0.0);
    assert_eq!(melting.matrix_enthalpy_derivative, 0.0);
    assert_eq!(melting.pcm_enthalpy_derivative, 4.0);
    assert_eq!(
        law.specific_entropy_change_j_kg_k(93_750.0, 143_750.0)
            .unwrap(),
        125.0
    );
}

#[test]
fn latent_storage_and_density_use_fixed_constituent_masses() {
    let matrix = matrix(0.0);
    let pcm = pcm(0.0);
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    assert!(std::ptr::eq(law.matrix_curve(), &matrix));
    assert!(std::ptr::eq(law.pcm_curve(), &pcm));
    assert_eq!(law.enthalpy_range_j_kg(), [0.0, 312_500.0]);
    assert_eq!(
        law.knots()
            .iter()
            .map(|k| k.specific_enthalpy_j_kg)
            .collect::<Vec<_>>(),
        vec![0.0, 125_000.0, 175_000.0, 312_500.0]
    );
    close(law.maximum_temperature_slope(), 1.0 / 1250.0, 1e-18);
    close(
        law.temperature_at_specific_enthalpy(62_500.0).unwrap(),
        350.0,
        1e-12,
    );
    let half = law.values_at_specific_enthalpy(150_000.0).unwrap();
    assert_eq!(half.temperature_k, 400.0);
    assert_eq!(half.matrix.liquid_mass_fraction, 0.0);
    assert_eq!(half.pcm.liquid_mass_fraction, 0.5);
    assert_eq!(half.liquid_mass_fraction, 0.125);
    let pcm_specific_volume = 0.5 / 1000.0 + 0.5 / 800.0;
    let volume = 0.75 / 8000.0 + 0.25 * pcm_specific_volume;
    close(half.specific_volume_m3_kg, volume, 1e-19);
    close(half.bulk_density_kg_m3, 1.0 / volume, 1e-10);
    close(
        ThermalLaw::bulk_density_kg_m3(&law, 150_000.0)
            .unwrap()
            .unwrap(),
        1.0 / volume,
        1e-10,
    );
    assert_eq!(
        ThermalLaw::bulk_density_kg_m3(&law, -1.0),
        Err(DepositionError::MaterialDomain),
    );
    close(
        half.pcm_volume_fraction,
        0.25 * pcm_specific_volume / volume,
        1e-15,
    );
    assert!((half.pcm_volume_fraction - half.pcm_mass_fraction).abs() > 0.4);
    assert_eq!(
        law.values_at_specific_enthalpy(175_000.0)
            .unwrap()
            .liquid_mass_fraction,
        0.25
    );
    close(
        law.temperature_at_specific_enthalpy(243_750.0).unwrap(),
        450.0,
        1e-12,
    );
    let initial_volume = law.volume_m3(0.0, 2.0).unwrap();
    let melted_volume = law.volume_m3(175_000.0, 2.0).unwrap();
    close(
        melted_volume - initial_volume,
        2.0 * 0.25 * (1.0 / 800.0 - 1.0 / 1000.0),
        1e-18,
    );
    close(
        law.specific_entropy_change_j_kg_k(125_000.0, 175_000.0)
            .unwrap(),
        50_000.0 / 400.0,
        1e-12,
    );
}

/// G1: the native geometry consumer must refresh thickness across melting,
/// with constituent mass and latent storage retained by the original owner.
#[test]
fn density_drives_pcm_layer_resistance_without_changing_mass() {
    use fs_thermal::enthalpy::conduction::{
        ConductionObserver, ConductiveReservoir, ConductivityLaw, PlanarDensityLayer,
    };
    use fs_thermal::enthalpy::{
        Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation,
        Event, MaterialId, Reservoir, StepConfig, TimedEvent,
    };
    struct UnitConductivity;
    impl ConductivityLaw for UnitConductivity {
        fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
            Ok(1.0)
        }
    }
    let matrix = matrix(0.0);
    let pcm = pcm(0.0);
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let initial_h = 150_000.0;
    let volume = |h: f64| {
        0.75 / 8000.0 + 0.25 / 1000.0 + (h - 125_000.0) * (1.0 / 800.0 - 1.0 / 1000.0) / 200_000.0
    };
    let config = StepConfig {
        max_dt_s: 0.01,
        energy_tolerance_j: 1e-8,
    };
    let initial = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(0),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: initial_h,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: volume(initial_h),
                    surface_area_m2: 2.0,
                },
                boundary: Boundary {
                    ambient_k: 400.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: Some(Reservoir {
                        area_m2: 1.0,
                        conductance_w_k: 0.0,
                        temperature_k: 500.0,
                    }),
                },
            }),
        }],
        &[(MaterialId(0), &law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 10,
        },
    )
    .unwrap()
    .advance_to(0.0, config, || Ok(()))
    .unwrap()
    .0;
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(0), &UnitConductivity)],
        &[],
        &[ConductiveReservoir {
            cell_id: CellId(0),
            half_distance_m: 0.5 * volume(initial_h),
            surface_resistance_m2_k_w: 0.0,
        }],
        || Ok(()),
    )
    .unwrap()
    .with_material_density_layers(
        &initial,
        &[PlanarDensityLayer {
            cell_id: CellId(0),
            footprint_area_m2: 1.0,
        }],
        || Ok(()),
    )
    .unwrap();
    let (heated, observer, audit) = initial
        .advance_observed(&observer, 0.02, config, || Ok(()))
        .unwrap();
    let first_h = initial_h + 0.01 * 200.0 / volume(initial_h);
    let expected_h = first_h + 0.01 * 200.0 / volume(first_h);
    let cell = &heated.cells()[&CellId(0)];
    close(cell.specific_h_j_kg, expected_h, 1e-8);
    close(heated.temperature_k(CellId(0)).unwrap(), 400.0, 1e-12);
    assert_eq!(cell.mass_kg, 1.0);
    assert_eq!(cell.geometry, initial.cells()[&CellId(0)].geometry);
    close(audit.external_heat_j, expected_h - initial_h, 1e-8);
    close(audit.residual_j, 0.0, 1e-8);
    close(
        observer.last_planar_half_distance_m(CellId(0)).unwrap(),
        0.5 * volume(first_h),
        1e-16,
    );
}

#[test]
fn independent_constituent_datums_leave_the_physics_unchanged() {
    let matrix0 = matrix(0.0);
    let pcm0 = pcm(0.0);
    let shifted_matrix = matrix(1_073_741_824.0);
    let shifted_pcm = pcm(-268_435_456.0);
    let reference = PcmComposite::try_new(&matrix0, &pcm0, 0.25, -2000.0).unwrap();
    let shifted = PcmComposite::try_new(&shifted_matrix, &shifted_pcm, 0.25, -2000.0).unwrap();
    let shifted_composite =
        PcmComposite::try_new(&shifted_matrix, &shifted_pcm, 0.25, 30_000.0).unwrap();
    assert_ne!(reference.identity(), shifted.identity());
    assert_ne!(shifted.identity(), shifted_composite.identity());
    for h in [
        -2000.0, 60_500.0, 123_000.0, 148_000.0, 173_000.0, 241_750.0, 310_500.0,
    ] {
        let a = reference.values_at_specific_enthalpy(h).unwrap();
        let b = shifted.values_at_specific_enthalpy(h).unwrap();
        let c = shifted_composite
            .values_at_specific_enthalpy(h + 32_000.0)
            .unwrap();
        for state in [b, c] {
            close(state.temperature_k, a.temperature_k, 1e-12);
            close(state.bulk_density_kg_m3, a.bulk_density_kg_m3, 1e-10);
            close(state.liquid_mass_fraction, a.liquid_mass_fraction, 1e-15);
        }
    }
    let entropy = reference
        .specific_entropy_change_j_kg_k(-2000.0, 310_500.0)
        .unwrap();
    close(
        shifted
            .specific_entropy_change_j_kg_k(-2000.0, 310_500.0)
            .unwrap(),
        entropy,
        1e-10,
    );
    close(
        shifted_composite
            .specific_entropy_change_j_kg_k(30_000.0, 342_500.0)
            .unwrap(),
        entropy,
        1e-10,
    );
    close(
        entropy,
        1250.0 * (400.0_f64 / 300.0).ln() + 125.0 + 1375.0 * (500.0_f64 / 400.0).ln(),
        1e-10,
    );
}

#[test]
fn source_density_laws_and_all_caloric_breakpoints_survive_preparation() {
    let density_source = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([3; 32]),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 1000.0,
                bulk_density_kg_m3: 8000.0,
            },
            HeatCapacityKnot {
                temperature_k: 500.0,
                specific_heat_capacity_j_kg_k: 1000.0,
                bulk_density_kg_m3: 6000.0,
            },
        ],
        1e-8,
    )
    .unwrap();
    let pcm = pcm(0.0);
    let law = PcmComposite::try_new(&density_source, &pcm, 0.25, 0.0).unwrap();
    let state = law.values_at_specific_enthalpy(150_000.0).unwrap();
    // The original sourced sensible rho(T) is affine density, while the PCM's
    // authored latent coexistence uses additive specific volume.
    close(state.matrix.bulk_density_kg_m3, 7000.0, 1e-12);
    close(
        state.pcm.bulk_density_kg_m3,
        1.0 / (0.5 / 1000.0 + 0.5 / 800.0),
        1e-10,
    );
    close(
        state.bulk_density_kg_m3,
        1.0 / (0.75 / 7000.0 + 0.25 * (0.5 / 1000.0 + 0.5 / 800.0)),
        1e-10,
    );

    let kinked_matrix = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([4; 32]),
        SolidLiquidPhase::Solid,
        vec![
            knot(0.0, 300.0, 0.0, 8000.0),
            knot(50_000.0, 350.0, 0.0, 8000.0),
            knot(350_000.0, 500.0, 0.0, 8000.0),
        ],
    )
    .unwrap();
    let kinked = PcmComposite::try_new(&kinked_matrix, &pcm, 0.25, 0.0).unwrap();
    assert_eq!(kinked.knots().len(), 5);
    assert_eq!(kinked.knots()[1].temperature_k, 350.0);
    assert_eq!(kinked.knots()[1].specific_enthalpy_j_kg, 62_500.0);
    close(
        kinked.temperature_at_specific_enthalpy(31_250.0).unwrap(),
        325.0,
        1e-12,
    );
    close(
        kinked.temperature_at_specific_enthalpy(112_500.0).unwrap(),
        375.0,
        1e-12,
    );
    let entropy = kinked
        .specific_entropy_change_j_kg_k(0.0, 425_000.0)
        .unwrap();
    close(
        entropy,
        1250.0 * (350.0_f64 / 300.0).ln()
            + 2000.0 * (400.0_f64 / 350.0).ln()
            + 125.0
            + 2125.0 * (500.0_f64 / 400.0).ln(),
        1e-10,
    );
    close(
        kinked
            .specific_entropy_change_j_kg_k(425_000.0, 0.0)
            .unwrap(),
        -entropy,
        1e-10,
    );
    close(
        kinked
            .specific_entropy_change_j_kg_k(0.0, 187_500.0)
            .unwrap()
            + kinked
                .specific_entropy_change_j_kg_k(187_500.0, 425_000.0)
                .unwrap(),
        entropy,
        1e-10,
    );
}

#[test]
fn incomplete_support_and_ambiguous_matrix_plateaus_refuse() {
    let matrix = matrix(0.0);
    let pcm = pcm(0.0);
    for fraction in [0.0, 1.0, -0.1, f64::NAN] {
        assert!(PcmComposite::try_new(&matrix, &pcm, fraction, 0.0).is_err());
    }
    assert!(PcmComposite::try_new(&matrix, &pcm, 0.25, f64::INFINITY).is_err());
    assert!(PcmComposite::try_new(&matrix, &matrix, 0.25, 0.0).is_err());
    let short_matrix = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([5; 32]),
        SolidLiquidPhase::Solid,
        vec![
            knot(0.0, 310.0, 0.0, 8000.0),
            knot(190_000.0, 500.0, 0.0, 8000.0),
        ],
    )
    .unwrap();
    assert!(PcmComposite::try_new(&short_matrix, &pcm, 0.25, 0.0).is_err());
    let ambiguous = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([6; 32]),
        vec![
            knot(0.0, 300.0, 0.0, 8000.0),
            knot(100_000.0, 400.0, 0.0, 8000.0),
            knot(110_000.0, 400.0, 0.0, 8000.0),
            knot(210_000.0, 500.0, 0.0, 8000.0),
            knot(310_000.0, 600.0, 1.0, 7000.0),
        ],
    )
    .unwrap();
    assert!(PcmComposite::try_new(&ambiguous, &pcm, 0.25, 0.0).is_err());
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    assert!(law.temperature_at_specific_enthalpy(-1.0).is_err());
    assert!(law.temperature_at_specific_enthalpy(312_501.0).is_err());
    assert!(law.volume_m3(0.0, 0.0).is_err());
    // A pure isothermal PCM still has a unique allocation because the matrix
    // is sensible and stays at its one uniquely determined enthalpy.
    let isothermal = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([7; 32]),
        vec![
            knot(0.0, 400.0, 0.0, 1000.0),
            knot(200_000.0, 400.0, 1.0, 800.0),
        ],
    )
    .unwrap();
    let plateau = PcmComposite::try_new(&matrix, &isothermal, 0.25, 0.0).unwrap();
    assert_eq!(plateau.maximum_temperature_slope(), 0.0);
    assert_eq!(
        plateau
            .values_at_specific_enthalpy(25_000.0)
            .unwrap()
            .liquid_mass_fraction,
        0.125
    );
}
