use super::*;
use crate::phase::water::freezing::thermal::FreezingThermalControls;
use crate::phase::water::freezing::transport::PhaseConductivityMixture;
use crate::phase::water::freezing::{FreezingWater, FreezingWaterControls};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::ice::IceMeltingSolve;
use fs_thermochem::water::ice::bounds::IceIhHeatCapacityBoundSolve;
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

fn laws() -> (FreezingWaterThermalLaw, FreezingWaterConductivity) {
    let chart = FreezingWater::new(
        Pressure::new(101_325.0),
        FreezingWaterControls {
            minimum_temperature_k: 243.15,
            maximum_temperature_k: 303.15,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            maximum_gibbs_mismatch_j_kg: 30.0,
            melting: IceMeltingSolve {
                pressure_tolerance_pa: 1e-5,
                temperature_tolerance_k: 1e-10,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap();
    let caloric = FreezingWaterThermalLaw::new(
        chart,
        FreezingThermalControls {
            ice_heat_capacity: IceIhHeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
            liquid_heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )
    .unwrap();
    let transport =
        FreezingWaterConductivity::new(chart, PhaseConductivityMixture::Series, || false).unwrap();
    (caloric, transport)
}

fn budget(n: usize) -> DepositionBudget {
    DepositionBudget {
        max_cells: n,
        max_links: n - 1,
        max_events: 2 * n + 8,
        max_substeps: 100_000,
    }
}
fn controls(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 2e-6,
    }
}
fn implicit_controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        total_tolerance_j: 2e-6,
        ..ImplicitConductionControls::default()
    }
}
fn slab<'a>(
    caloric: &'a FreezingWaterThermalLaw,
    transport: &'a FreezingWaterConductivity,
    n: usize,
    cold: bool,
) -> FreezingSlab<'a> {
    FreezingSlab::new(
        caloric,
        transport,
        FreezingSlabGeometry {
            thickness_m: 0.01,
            area_m2: 0.001,
            cells: n,
        },
        caloric
            .chart()
            .state_at_temperature(Temperature::new(273.16), || false)
            .unwrap(),
        [
            cold.then_some(FreezingSlabReservoir {
                temperature_k: 253.15,
                surface_resistance_m2_k_w: 0.0,
            }),
            None,
        ],
        &[],
        budget(n),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn original_water_sources_compute_a_spatial_front_and_close_h_and_u_plus_pv() {
    let (caloric, transport) = laws();
    let initial = slab(&caloric, &transport, 24, true);
    let before = initial.inventory(|| Ok(())).unwrap();
    let (next, audit) = initial.advance_to(60.0, controls(0.2), || Ok(())).unwrap();
    let after = next.inventory(|| Ok(())).unwrap();
    let profile = next.profile(|| Ok(())).unwrap();
    assert_eq!(next.time_s(), 60.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.mass_kg(), next.mass_kg());
    assert_eq!(after.mass_kg, before.mass_kg);
    assert!((after.ice_mass_kg + after.liquid_mass_kg - after.mass_kg).abs() < 2e-17);
    assert!(after.ice_mass_kg > 0.0 && after.liquid_mass_kg > 0.0);
    assert!(profile[0].material.temperature_k() < 260.0);
    assert_eq!(profile[0].material.ice_mass_fraction(), 1.0);
    assert_eq!(profile.last().unwrap().material.liquid_mass_fraction(), 1.0);
    assert!(profile.iter().any(|cell| {
        let x = cell.material.liquid_mass_fraction();
        x > 0.0 && x < 1.0
    }));
    for pair in profile.windows(2) {
        assert!(pair[0].material.temperature_k() <= pair[1].material.temperature_k() + 1e-8);
        assert!(
            pair[0].material.liquid_mass_fraction()
                <= pair[1].material.liquid_mass_fraction() + 1e-8
        );
    }
    assert!((0.002..0.006).contains(&after.equivalent_frozen_reference_length_m));
    assert!(after.free_volume_m3() > before.free_volume_m3());
    assert!(audit.thermal.external_heat_j < 0.0);
    assert!(audit.pressure_work_j > 0.0);
    assert!(audit.source_enthalpy_residual_j.abs() < 2e-6);
    assert!(audit.internal_energy_residual_j.abs() < 2e-6);
    assert!(audit.thermal.residual_j.abs() < 2e-6);
    assert!(!audit.entropy.has_unclassified_sources);
    assert!(audit.entropy.conduction_generation_j_k > 0.0);
    assert!(
        audit.entropy.thermal_storage_change_j_k - audit.entropy.reservoir_entropy_input_j_k > 0.0
    );
    for cell in next.thermal().cells().values() {
        assert_eq!(cell.geometry, initial.thermal().cells()[&cell.id].geometry);
        assert_eq!(cell.boundary.source_w, 0.0);
    }
}

#[test]
fn spatial_and_time_refinement_resolve_source_front_independently_of_energy_tolerance() {
    let (caloric, transport) = laws();
    let run = |n, dt| {
        let (state, audit) = slab(&caloric, &transport, n, true)
            .advance_to(30.0, controls(dt), || Ok(()))
            .unwrap();
        (
            state
                .inventory(|| Ok(()))
                .unwrap()
                .equivalent_frozen_reference_length_m,
            audit.thermal.external_heat_j,
            audit.entropy.temporal_defect_j_k.abs(),
        )
    };
    let spatial: Vec<_> = [12, 24, 48].map(|n| run(n, 0.02)).into();
    eprintln!("source spatial refinement: {spatial:?}");
    assert!((spatial[1].0 - spatial[2].0).abs() < (spatial[0].0 - spatial[2].0).abs());
    assert!((spatial[1].1 - spatial[2].1).abs() < (spatial[0].1 - spatial[2].1).abs());
    let temporal: Vec<_> = [0.2, 0.1, 0.05].map(|dt| run(12, dt)).into();
    eprintln!("source time refinement: {temporal:?}");
    assert!((temporal[1].0 - temporal[2].0).abs() < (temporal[0].0 - temporal[2].0).abs());
    assert!(temporal[2].2 < temporal[1].2 && temporal[1].2 < temporal[0].2);
}

#[test]
fn insulated_inventory_is_stationary_and_bath_temperature_events_reverse_the_front() {
    let (caloric, transport) = laws();
    let insulated = slab(&caloric, &transport, 8, false);
    let (unchanged, audit) = insulated
        .advance_to(10.0, controls(0.2), || Ok(()))
        .unwrap();
    assert_eq!(
        unchanged.inventory(|| Ok(())).unwrap(),
        insulated.inventory(|| Ok(())).unwrap()
    );
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.thermal.stored_energy_change_j, 0.0);
    assert_eq!(audit.entropy, EntropyAudit::default());

    let slab = FreezingSlab::new(
        &caloric,
        &transport,
        FreezingSlabGeometry {
            thickness_m: 0.01,
            area_m2: 0.001,
            cells: 12,
        },
        caloric
            .chart()
            .state_at_temperature(Temperature::new(273.16), || false)
            .unwrap(),
        [
            Some(FreezingSlabReservoir {
                temperature_k: 253.15,
                surface_resistance_m2_k_w: 0.0,
            }),
            None,
        ],
        &[FreezingSlabTemperatureEvent {
            time_s: 10.0,
            face: SlabFace::Left,
            temperature_k: 293.15,
        }],
        budget(12),
        || Ok(()),
    )
    .unwrap();
    let (frozen, _) = slab.advance_to(10.0, controls(0.1), || Ok(())).unwrap();
    let frozen_mass = frozen.inventory(|| Ok(())).unwrap().ice_mass_kg;
    assert!(frozen_mass > 0.0);
    assert_eq!(
        frozen.thermal().cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .temperature_k,
        293.15
    );
    let (thawed, audit) = frozen.advance_to(40.0, controls(0.1), || Ok(())).unwrap();
    assert!(thawed.inventory(|| Ok(())).unwrap().ice_mass_kg < frozen_mass);
    assert!(audit.thermal.external_heat_j > 0.0);
    assert!(!audit.entropy.has_unclassified_sources);
}

#[test]
fn spatial_failure_and_cancellation_never_publish_a_partial_thermal_state() {
    let (caloric, transport) = laws();
    let slab = slab(&caloric, &transport, 8, true);
    let original = slab.profile(|| Ok(())).unwrap();
    let mut calls = 0;
    let error = slab
        .advance_to(5.0, controls(0.1), || {
            calls += 1;
            if calls >= 200 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .err()
        .unwrap();
    assert_eq!(
        error,
        FreezingSlabError::Thermal(DepositionError::Cancelled)
    );
    assert_eq!(slab.time_s(), 0.0);
    assert_eq!(slab.profile(|| Ok(())).unwrap(), original);

    let mut tight = budget(8);
    tight.max_substeps = 1;
    let limited = FreezingSlab::new(
        &caloric,
        &transport,
        slab.geometry,
        original[0].material,
        [
            Some(FreezingSlabReservoir {
                temperature_k: 253.15,
                surface_resistance_m2_k_w: 0.0,
            }),
            None,
        ],
        &[],
        tight,
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        limited.advance_to(1.0, controls(0.01), || Ok(())),
        Err(FreezingSlabError::Thermal(DepositionError::Budget(_)))
    ));
    assert_eq!(limited.time_s(), 0.0);
    assert!(
        FreezingSlab::new(
            &caloric,
            &transport,
            FreezingSlabGeometry {
                cells: 1,
                ..slab.geometry
            },
            original[0].material,
            [None, None],
            &[],
            budget(8),
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn implicit_source_front_crosses_the_explicit_cfl_and_closes_original_calorimetry() {
    let (caloric, transport) = laws();
    let initial = slab(&caloric, &transport, 24, true);
    let before = initial.profile(|| Ok(())).unwrap();
    let k = transport
        .state_at_material_state(before[0].material, || false)
        .unwrap()
        .conductivity_w_m_k();
    let dx = initial.geometry().thickness_m / 24.0;
    // The first cell sees one neighbour plus a bath through its half-cell.
    let outgoing_g = 3.0 * k * initial.geometry().area_m2 / dx;
    let initial_explicit_cfl =
        0.8 * before[0].reference_mass_kg / (caloric.maximum_temperature_slope() * outgoing_g);
    let control = implicit_controls(1.0);
    let mut workspace = ImplicitConductionWorkspace::new(24).unwrap();
    let (next, audit) = initial
        .advance_implicit_to(60.0, control, &mut workspace, || Ok(()))
        .unwrap();
    let after = next.profile(|| Ok(())).unwrap();
    let inventory = next.inventory(|| Ok(())).unwrap();
    assert!(60.0 / audit.solve.energy.substeps as f64 > 2.0 * initial_explicit_cfl);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(next.time_s(), 60.0);
    assert_eq!(initial.mass_kg(), next.mass_kg());
    assert_eq!(after[0].material.ice_mass_fraction(), 1.0);
    assert_eq!(after.last().unwrap().material.liquid_mass_fraction(), 1.0);
    let mixed = after
        .iter()
        .find(|cell| {
            let x = cell.material.liquid_mass_fraction();
            x > 0.0 && x < 1.0
        })
        .expect("the spatial front must include a genuine latent-heat cell");
    assert_eq!(
        mixed.material.temperature_k(),
        caloric.chart().melting().temperature().value()
    );
    assert!((0.002..0.006).contains(&inventory.equivalent_frozen_reference_length_m));
    for pair in after.windows(2) {
        assert!(pair[0].material.temperature_k() <= pair[1].material.temperature_k() + 1e-8);
        assert!(
            pair[0].material.liquid_mass_fraction()
                <= pair[1].material.liquid_mass_fraction() + 1e-8
        );
    }
    assert_eq!(audit.physical.thermal, audit.solve.energy);
    assert_eq!(audit.physical.entropy, audit.solve.entropy);
    assert!(audit.solve.maximum_cell_residual_j <= control.cell_tolerance_j);
    assert!(audit.physical.source_enthalpy_residual_j.abs() < control.total_tolerance_j);
    assert!(audit.physical.internal_energy_residual_j.abs() < control.total_tolerance_j);
    assert!(audit.physical.pressure_work_j > 0.0);
    assert!(audit.physical.thermal.external_heat_j < 0.0);
    let entropy_from_sources: f64 = before
        .iter()
        .zip(&after)
        .map(|(a, b)| {
            a.reference_mass_kg
                * (b.material.thermal_specific_entropy_j_kg_k()
                    - a.material.thermal_specific_entropy_j_kg_k())
        })
        .sum();
    let entropy = audit.physical.entropy;
    assert!((entropy.thermal_storage_change_j_k - entropy_from_sources).abs() < 1e-10);
    assert!(
        (entropy.reservoir_entropy_input_j_k - audit.physical.thermal.external_heat_j / 253.15)
            .abs()
            < 1e-10
    );
    assert!(!entropy.has_unclassified_sources);
    assert!(entropy.conduction_generation_j_k > 0.0);
    assert!(entropy.thermal_storage_change_j_k - entropy.reservoir_entropy_input_j_k > 0.0);
}

#[test]
fn implicit_time_refinement_converges_to_the_fine_explicit_source_profile() {
    let (caloric, transport) = laws();
    let initial = slab(&caloric, &transport, 12, true);
    let (reference, reference_audit) = initial
        .advance_to(30.0, controls(0.005), || Ok(()))
        .unwrap();
    let reference_profile = reference.profile(|| Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(12).unwrap();
    let errors = [1.0, 0.5, 0.25].map(|dt| {
        let (next, audit) = initial
            .advance_implicit_to(30.0, implicit_controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        let profile_error_j: f64 = next
            .profile(|| Ok(()))
            .unwrap()
            .iter()
            .zip(&reference_profile)
            .map(|(a, b)| {
                a.reference_mass_kg
                    * (a.stored_specific_enthalpy_j_kg - b.stored_specific_enthalpy_j_kg).abs()
            })
            .sum();
        (
            profile_error_j,
            (audit.physical.thermal.external_heat_j - reference_audit.thermal.external_heat_j)
                .abs(),
            audit.physical.entropy.temporal_defect_j_k.abs(),
        )
    });
    eprintln!("implicit source time refinement (field J, bath J, defect J/K): {errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[1].0 < pair[0].0);
        assert!(pair[1].1 < pair[0].1);
        assert!(pair[1].2 < pair[0].2);
    }
    // This comparison tests the time-discrete field, independently of the
    // much smaller algebraic energy acceptance tolerance used by every solve.
    assert!(errors[2].0 < 0.5 * errors[0].0);
    assert!(errors[2].0 > implicit_controls(0.25).total_tolerance_j);
}

#[test]
fn implicit_bath_event_is_exact_and_chunked_continuation_preserves_the_front() {
    let (caloric, transport) = laws();
    let event_time = 10.125;
    let initial = FreezingSlab::new(
        &caloric,
        &transport,
        FreezingSlabGeometry {
            thickness_m: 0.01,
            area_m2: 0.001,
            cells: 12,
        },
        caloric
            .chart()
            .state_at_temperature(Temperature::new(273.16), || false)
            .unwrap(),
        [
            Some(FreezingSlabReservoir {
                temperature_k: 253.15,
                surface_resistance_m2_k_w: 0.0,
            }),
            None,
        ],
        &[FreezingSlabTemperatureEvent {
            time_s: event_time,
            face: SlabFace::Left,
            temperature_k: 293.15,
        }],
        budget(12),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(12).unwrap();
    let (whole, whole_audit) = initial
        .advance_implicit_to(40.125, implicit_controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    let (frozen, cold_audit) = initial
        .advance_implicit_to(
            event_time,
            implicit_controls(1.0),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(frozen.time_s(), event_time);
    assert_eq!(
        frozen.thermal().cells()[&CellId(0)]
            .boundary
            .reservoir
            .unwrap()
            .temperature_k,
        293.15
    );
    assert_eq!(cold_audit.physical.thermal.events, 1);
    assert!(cold_audit.physical.thermal.external_heat_j < 0.0);
    let (thawed, warm_audit) = frozen
        .advance_implicit_to(40.125, implicit_controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    assert!(warm_audit.physical.thermal.external_heat_j > 0.0);
    assert!(
        thawed.inventory(|| Ok(())).unwrap().ice_mass_kg
            < frozen.inventory(|| Ok(())).unwrap().ice_mass_kg
    );
    assert_eq!(warm_audit.physical.thermal.events, 0);
    assert_eq!(whole_audit.physical.thermal.events, 1);
    assert_eq!(
        thawed.profile(|| Ok(())).unwrap(),
        whole.profile(|| Ok(())).unwrap()
    );
    assert!(
        (whole_audit.physical.thermal.external_heat_j
            - cold_audit.physical.thermal.external_heat_j
            - warm_audit.physical.thermal.external_heat_j)
            .abs()
            < 1e-10
    );
    assert!(
        (whole_audit.physical.entropy.reservoir_entropy_input_j_k
            - cold_audit.physical.entropy.reservoir_entropy_input_j_k
            - warm_audit.physical.entropy.reservoir_entropy_input_j_k)
            .abs()
            < 1e-10
    );
}

#[test]
fn implicit_refusal_preserves_physical_state_and_reuses_numerical_workspace() {
    let (caloric, transport) = laws();
    let initial = slab(&caloric, &transport, 8, true);
    let original = initial.profile(|| Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(8).unwrap();
    let mut calls = 0;
    let cancelled = initial
        .advance_implicit_to(5.0, implicit_controls(0.5), &mut workspace, || {
            calls += 1;
            if calls >= 300 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .err()
        .unwrap();
    assert_eq!(
        cancelled,
        FreezingSlabError::Implicit(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    );
    let exhausted = ImplicitConductionControls {
        maximum_residual_evaluations: 1,
        ..implicit_controls(0.5)
    };
    assert!(matches!(
        initial.advance_implicit_to(5.0, exhausted, &mut workspace, || Ok(())),
        Err(FreezingSlabError::Implicit(
            ImplicitConductionError::Thermal(DepositionError::Budget(
                "implicit residual evaluations"
            ))
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.profile(|| Ok(())).unwrap(), original);
    let (reused, reused_audit) = initial
        .advance_implicit_to(5.0, implicit_controls(0.5), &mut workspace, || Ok(()))
        .unwrap();
    let (fresh, fresh_audit) = initial
        .advance_implicit_to(
            5.0,
            implicit_controls(0.5),
            &mut ImplicitConductionWorkspace::new(8).unwrap(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(
        reused.profile(|| Ok(())).unwrap(),
        fresh.profile(|| Ok(())).unwrap()
    );
    assert_eq!(reused_audit, fresh_audit);
}
