use super::*;
use crate::phase::water::freezing::{FreezingWaterControls, FreezingWaterState};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLink, TimedEvent,
};
use fs_thermochem::water::ice::IceMeltingSolve;

const MASS: f64 = 0.05;
const TOLERANCE: f64 = 2e-6;

fn law() -> FreezingWaterThermalLaw {
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
    FreezingWaterThermalLaw::new(
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
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance {tolerance:.4e}"
    );
}
fn at_t(law: &FreezingWaterThermalLaw, temperature: f64) -> FreezingWaterState {
    law.chart()
        .state_at_temperature(Temperature::new(temperature), || false)
        .unwrap()
}
fn at_h(law: &FreezingWaterThermalLaw, enthalpy: f64) -> FreezingWaterState {
    law.chart()
        .state_at_specific_enthalpy(enthalpy, || false)
        .unwrap()
        .into_state()
}
fn boundary(reservoir: Option<(f64, f64)>, source: f64) -> Boundary {
    Boundary {
        ambient_k: 293.15,
        htc_w_m2_k: 0.0,
        reservoir: reservoir.map(|(temperature, conductance)| Reservoir {
            area_m2: 0.1,
            conductance_w_k: conductance,
            temperature_k: temperature,
        }),
        source_w: source,
    }
}
fn birth(id: u64, state: FreezingWaterState, boundary: Boundary) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: id,
            material_id: MaterialId(1),
            mass_kg: MASS,
            specific_h_j_kg: state.specific_enthalpy_j_kg(),
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64, 0.1, 0.0],
                volume_m3: MASS * state.specific_volume_m3_kg(),
                surface_area_m2: 0.1,
            },
            boundary,
        }),
    }
}
fn budget() -> DepositionBudget {
    DepositionBudget {
        max_cells: 2,
        max_links: 1,
        max_events: 3,
        max_substeps: 8192,
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: TOLERANCE,
    }
}

#[test]
fn thermal_preparation_retains_both_source_bounds_and_admits_a_zero_slope_latent_chart() {
    let law = law();
    let chart = *law.chart();
    let ice = law.ice_heat_capacity_bound().unwrap();
    let liquid = law.liquid_heat_capacity_bound().unwrap();
    assert!(ice.interval_evaluations() > 0 && liquid.interval_evaluations() > 0);
    assert_eq!(
        law.maximum_temperature_slope(),
        ice.maximum_temperature_slope()
            .max(liquid.maximum_temperature_slope())
    );
    for temperature in [243.15, 253.15, 263.15, 273.0] {
        let source = at_t(&law, temperature).ice_source().unwrap();
        let cp = source.cp().value();
        assert!(cp >= ice.minimum_heat_capacity_j_kg_k());
        assert!(cp <= ice.maximum_heat_capacity_j_kg_k());
    }
    for temperature in [273.16, 283.15, 293.15, 303.15] {
        let source = at_t(&law, temperature).liquid_source().unwrap();
        let cp = source.cp().value();
        assert!(cp >= liquid.minimum_heat_capacity_j_kg_k());
        assert!(cp <= liquid.maximum_heat_capacity_j_kg_k());
    }
    let controls = FreezingThermalControls {
        ice_heat_capacity: IceIhHeatCapacityBoundSolve {
            maximum_intervals: 8191,
        },
        liquid_heat_capacity: If97HeatCapacityBoundSolve {
            maximum_intervals: 8191,
        },
    };
    assert!(matches!(
        FreezingWaterThermalLaw::new(chart, controls, || true),
        Err(FreezingThermalError::Cancelled)
    ));
    let wide = FreezingWater::new(
        Pressure::new(chart.pressure_pa()),
        FreezingWaterControls {
            minimum_temperature_k: 100.0,
            ..*chart.controls()
        },
        || false,
    )
    .unwrap();
    assert!(matches!(
        FreezingWaterThermalLaw::new(
            wide,
            FreezingThermalControls {
                ice_heat_capacity: IceIhHeatCapacityBoundSolve {
                    maximum_intervals: 1,
                },
                ..controls
            },
            || false,
        ),
        Err(FreezingThermalError::IceHeatCapacity(
            IceIhHeatCapacityBoundError::NotProven { .. }
        ))
    ));

    let melting_temperature = chart.melting().temperature().value();
    let plateau = FreezingWater::new(
        Pressure::new(chart.pressure_pa()),
        FreezingWaterControls {
            minimum_temperature_k: melting_temperature,
            maximum_temperature_k: melting_temperature,
            ..*chart.controls()
        },
        || false,
    )
    .unwrap();
    // No sensible branch exists, so neither unused proof budget is consumed.
    let plateau = FreezingWaterThermalLaw::new(
        plateau,
        FreezingThermalControls {
            ice_heat_capacity: IceIhHeatCapacityBoundSolve {
                maximum_intervals: 0,
            },
            liquid_heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 0,
            },
        },
        || false,
    )
    .unwrap();
    assert_eq!(plateau.ice_heat_capacity_bound(), None);
    assert_eq!(plateau.liquid_heat_capacity_bound(), None);
    assert_eq!(plateau.maximum_temperature_slope(), 0.0);
    let [lower_h, upper_h] = plateau.chart().enthalpy_bounds_j_kg();
    let ice_density = at_h(&plateau, lower_h).bulk_density_kg_m3();
    let liquid_density = at_h(&plateau, upper_h).bulk_density_kg_m3();
    assert!(ice_density < liquid_density);
    for fraction in [0.0, 0.25, 0.75, 1.0] {
        let h = match fraction {
            0.0 => lower_h,
            1.0 => upper_h,
            _ => lower_h + fraction * (upper_h - lower_h),
        };
        assert_eq!(
            ThermalLaw::temperature_k(&plateau, h).unwrap(),
            melting_temperature
        );
        near(
            ThermalLaw::bulk_density_kg_m3(&plateau, h)
                .unwrap()
                .unwrap(),
            1.0 / ((1.0 - fraction) / ice_density + fraction / liquid_density),
            2e-10,
        );
    }
    near(
        ThermalLaw::specific_entropy_change_j_kg_k(&plateau, lower_h, upper_h)
            .unwrap()
            .unwrap(),
        plateau.chart().melting().latent_enthalpy_j_kg() / melting_temperature,
        2e-12,
    );
    assert!(matches!(
        ThermalLaw::temperature_k(&plateau, lower_h.next_down()),
        Err(DepositionError::MaterialDomain)
    ));
    assert_eq!(
        ThermalLaw::bulk_density_kg_m3(&plateau, lower_h.next_down()),
        Err(DepositionError::MaterialDomain),
    );
}

#[test]
fn thermal_entropy_retains_single_ulp_ice_and_liquid_heat() {
    let law = law();
    for temperature in [253.15, 293.15] {
        let h = at_t(&law, temperature).specific_enthalpy_j_kg();
        let next_h = h.next_up();
        let heat = next_h - h;
        assert!(heat > 0.0 && heat < law.chart().controls().enthalpy_tolerance_j_kg);
        let entropy = ThermalLaw::specific_entropy_change_j_kg_k(&law, h, next_h)
            .unwrap()
            .unwrap();
        assert!(entropy > 0.0);
        near(entropy / (heat / temperature), 1.0, 1e-10);
        assert_eq!(
            ThermalLaw::specific_entropy_change_j_kg_k(&law, next_h, h)
                .unwrap()
                .unwrap(),
            -entropy
        );
        assert_eq!(
            ThermalLaw::specific_entropy_change_j_kg_k(&law, h, h)
                .unwrap()
                .unwrap(),
            0.0
        );
    }
    let [lo, hi] = law.chart().enthalpy_bounds_j_kg();
    assert_eq!(
        ThermalLaw::specific_entropy_change_j_kg_k(&law, hi.next_up(), hi.next_up()),
        Err(DepositionError::MaterialDomain)
    );
    assert_eq!(
        law.chart().specific_entropy_change_j_kg_k(lo, lo, || true),
        Err(crate::phase::water::freezing::FreezingWaterError::Cancelled)
    );
}

#[test]
fn thermal_entropy_retains_stored_fusion_heat_across_fraction_roundoff() {
    let law = law();
    let melting = law.chart().melting();
    let ice_h = melting.ice().specific_enthalpy().value();
    let liquid_h = melting.liquid().specific_enthalpy().value();
    let middle_h = ice_h + 0.5 * (liquid_h - ice_h);
    let temperature = melting.temperature().value();
    for (initial_h, final_h) in [
        (ice_h, liquid_h),
        (middle_h, middle_h.next_up()),
        (middle_h, liquid_h),
    ] {
        let entropy = ThermalLaw::specific_entropy_change_j_kg_k(&law, initial_h, final_h)
            .unwrap()
            .unwrap();
        // Fusion uses the actual stored enthalpy and the isobar's latent h/T,
        // including its source seam, independent of rounded phase fractions.
        assert_eq!(entropy, (final_h - initial_h) / temperature);
        assert_eq!(
            ThermalLaw::specific_entropy_change_j_kg_k(&law, final_h, initial_h)
                .unwrap()
                .unwrap(),
            -entropy
        );
    }
}

fn reservoir_history(law: &FreezingWaterThermalLaw, dt: f64) -> (FreezingWaterState, [f64; 4]) {
    let warm = at_t(law, 293.15);
    let initial = DepositionSimulation::new(
        vec![
            birth(1, warm, boundary(Some((253.15, 5.0)), 0.0)),
            TimedEvent {
                time_s: 300.0,
                event: Event::Boundary(CellId(1), boundary(Some((293.15, 5.0)), 0.0)),
            },
        ],
        &[(MaterialId(1), law)],
        budget(),
    )
    .unwrap();
    let (mut simulation, _) = initial.advance_to(0.0, config(dt), || Ok(())).unwrap();
    let reference_geometry = simulation.cells()[&CellId(1)].geometry;
    let mut heat = 0.0;
    let mut entropy = [0.0; 4];
    let mut final_state = warm;
    for time in [100.0, 300.0, 400.0, 650.0] {
        let (next, audit, s) = simulation
            .advance_to_with_entropy(time, config(dt), || Ok(()))
            .unwrap();
        heat += audit.external_heat_j;
        for (total, value) in entropy.iter_mut().zip([
            s.thermal_storage_change_j_k,
            s.reservoir_entropy_input_j_k,
            s.conduction_generation_j_k,
            s.temporal_defect_j_k,
        ]) {
            *total += value;
        }
        assert!(!s.has_unclassified_sources);
        assert_eq!(s.unclassified_source_at_cell_j_k, 0.0);
        near(audit.residual_j, 0.0, TOLERANCE);
        let state = at_h(law, next.cells()[&CellId(1)].specific_h_j_kg);
        let dh = MASS * (state.specific_enthalpy_j_kg() - warm.specific_enthalpy_j_kg());
        let du =
            MASS * (state.specific_internal_energy_j_kg() - warm.specific_internal_energy_j_kg());
        let work = law.chart().pressure_pa()
            * MASS
            * (state.specific_volume_m3_kg() - warm.specific_volume_m3_kg());
        near(dh, heat, TOLERANCE);
        near(du + work, heat, TOLERANCE);
        near(
            entropy[0],
            MASS * law
                .chart()
                .specific_entropy_change_j_kg_k(
                    warm.specific_enthalpy_j_kg(),
                    state.specific_enthalpy_j_kg(),
                    || false,
                )
                .unwrap(),
            2e-8,
        );
        assert_eq!(next.cells()[&CellId(1)].geometry, reference_geometry);
        assert_eq!(next.mass_kg(), MASS);
        if time == 100.0 || time == 400.0 {
            near(
                state.temperature_k(),
                law.chart().melting().temperature().value(),
                0.0,
            );
            assert!(state.liquid_mass_fraction() > 0.0 && state.liquid_mass_fraction() < 1.0);
            assert!(state.liquid_volume_fraction() < state.liquid_mass_fraction());
        } else if time == 300.0 {
            assert_eq!(state.ice_mass_fraction(), 1.0);
            assert!(state.temperature_k() < 260.0);
            assert!(state.specific_volume_m3_kg() > warm.specific_volume_m3_kg());
            assert!(heat < 0.0 && work > 0.0);
        } else {
            assert_eq!(state.liquid_mass_fraction(), 1.0);
            assert!(state.temperature_k() > 290.0 && state.temperature_k() < 293.15);
        }
        final_state = state;
        simulation = next;
    }
    assert!(entropy[2] > 0.0);
    assert!(entropy[0] - entropy[1] > 0.0);
    near(entropy[0] - entropy[1] - entropy[2], entropy[3], 2e-12);
    (final_state, entropy)
}

#[test]
fn reservoir_driven_freezing_and_thawing_preserve_energy_volume_and_thermal_entropy() {
    let law = law();
    let (coarse, coarse_entropy) = reservoir_history(&law, 0.5);
    let (fine, fine_entropy) = reservoir_history(&law, 0.25);
    near(coarse.temperature_k(), fine.temperature_k(), 0.05);
    // The timeline reports its first-order entropy quadrature defect instead
    // of identifying that discretization residual as physical generation.
    assert!(fine_entropy[3].abs() < 0.6 * coarse_entropy[3].abs());
    assert!(coarse_entropy[3] < 0.0 && fine_entropy[3] < 0.0);
}

#[test]
fn finite_ice_water_contact_uses_one_internal_heat_transfer_and_real_fusion_storage() {
    let law = law();
    let cold = at_t(&law, 253.15);
    let warm = at_t(&law, 293.15);
    let initial = DepositionSimulation::new(
        vec![
            birth(1, cold, boundary(None, 0.0)),
            birth(2, warm, boundary(None, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.05,
                    conductance_w_k: 2.0,
                }),
            },
        ],
        &[(MaterialId(1), &law)],
        budget(),
    )
    .unwrap();
    let (active, _) = initial.advance_to(0.0, config(0.5), || Ok(())).unwrap();
    let (next, audit, entropy) = active
        .advance_to_with_entropy(1000.0, config(0.5), || Ok(()))
        .unwrap();
    assert_eq!(audit.external_heat_j, 0.0);
    near(audit.residual_j, 0.0, TOLERANCE);
    let cold_final = at_h(&law, next.cells()[&CellId(1)].specific_h_j_kg);
    let warm_final = at_h(&law, next.cells()[&CellId(2)].specific_h_j_kg);
    assert!(cold_final.liquid_mass_fraction() > 0.0 && cold_final.liquid_mass_fraction() < 1.0);
    near(
        cold_final.temperature_k(),
        law.chart().melting().temperature().value(),
        0.0,
    );
    assert!(warm_final.temperature_k() > cold_final.temperature_k());
    assert!(warm_final.temperature_k() < 293.15);
    let mut total_h = 0.0;
    let mut total_u_pv = 0.0;
    for (before, after) in [(cold, cold_final), (warm, warm_final)] {
        total_h += MASS * (after.specific_enthalpy_j_kg() - before.specific_enthalpy_j_kg());
        total_u_pv += MASS
            * (after.specific_internal_energy_j_kg() - before.specific_internal_energy_j_kg()
                + law.chart().pressure_pa()
                    * (after.specific_volume_m3_kg() - before.specific_volume_m3_kg()));
    }
    near(total_h, 0.0, TOLERANCE);
    near(total_u_pv, 0.0, TOLERANCE);
    assert_eq!(next.mass_kg(), 2.0 * MASS);
    assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
    assert!(!entropy.has_unclassified_sources);
    assert!(entropy.thermal_storage_change_j_k > 0.0);
    assert!(entropy.conduction_generation_j_k > 0.0);
    assert_eq!(active.time_s(), 0.0);
    assert_eq!(
        active.cells()[&CellId(1)].specific_h_j_kg,
        cold.specific_enthalpy_j_kg()
    );
}

#[test]
fn cancelled_and_out_of_domain_thermal_advances_leave_freezing_history_unchanged() {
    let law = law();
    let warm = at_t(&law, 293.15);
    let initial = DepositionSimulation::new(
        vec![birth(1, warm, boundary(Some((253.15, 5.0)), 0.0))],
        &[(MaterialId(1), &law)],
        budget(),
    )
    .unwrap();
    let (active, _) = initial.advance_to(0.0, config(0.5), || Ok(())).unwrap();
    let saved = active.cells().clone();
    let mut polls = 0;
    let cancelled = active.advance_to_with_entropy(300.0, config(0.5), || {
        polls += 1;
        if polls == 12 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(cancelled, Err(DepositionError::Cancelled)));
    assert_eq!(active.cells(), &saved);
    assert_eq!(active.time_s(), 0.0);
    let too_hot = DepositionSimulation::new(
        vec![birth(1, warm, boundary(None, 1e6))],
        &[(MaterialId(1), &law)],
        budget(),
    )
    .unwrap();
    let (active, _) = too_hot.advance_to(0.0, config(0.5), || Ok(())).unwrap();
    let saved = active.cells().clone();
    assert!(matches!(
        active.advance_to_with_entropy(1.0, config(0.5), || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert_eq!(active.cells(), &saved);
    assert_eq!(active.time_s(), 0.0);
}
