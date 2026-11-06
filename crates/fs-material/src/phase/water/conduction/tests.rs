use super::*;
use crate::phase::water::{WaterChartControls, WaterError, WaterHeatCapacityResponse};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLink, TimedEvent,
};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use fs_thermochem::water::if97::{If97Region, If97State};

fn chart(pressure: f64) -> IsobaricWater {
    IsobaricWater::new(
        Pressure::new(pressure),
        WaterChartControls {
            minimum_temperature_k: 300.0,
            maximum_temperature_k: 650.0,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 65535,
            },
        },
        || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}

fn enthalpy(chart: &IsobaricWater, temperature: f64) -> f64 {
    chart
        .state_at_temperature(Temperature::new(temperature), || false)
        .unwrap()
        .specific_enthalpy_j_kg()
}

#[test]
fn g1_sensible_cp_and_fused_tangents_retain_original_if97_derivatives() {
    for pressure in [101325.0, 3e6] {
        let chart = chart(pressure);
        let law = WaterConductionLaw::try_new(
            &chart,
            ConductivityCurve::new(vec![[300.0, 0.3], [500.0, 0.5], [650.0, 0.2]]).unwrap(),
        )
        .unwrap();
        assert!(std::ptr::eq(law.caloric(), &chart));
        assert!(std::ptr::addr_eq(law.caloric_law(), &chart as &dyn ThermalLaw));
        assert_eq!(law.enthalpy_bounds_j_kg(), chart.enthalpy_bounds_j_kg());
        let ts = chart.saturation().temperature().value();
        for t in [315.0, 0.5 * (300.0 + ts), 0.5 * (ts + 650.0), 625.0] {
            let material = chart
                .state_at_temperature(Temperature::new(t), || false)
                .unwrap();
            let region = if t < ts {
                If97Region::Region1
            } else {
                If97Region::Region2
            };
            let source = If97State::at_pressure_temperature(
                Pressure::new(pressure),
                Temperature::new(t),
                region,
            )
            .unwrap();
            let cp = source.cp().value();
            assert_eq!(material.specific_heat_cp_j_kg_k().unwrap(), cp);
            assert_eq!(
                material.isobaric_heat_capacity_response(),
                WaterHeatCapacityResponse::Sensible {
                    specific_heat_cp_j_kg_k: cp
                }
            );
            // Independent source enthalpy differences, not the production tangent.
            let dh_dt = (enthalpy(&chart, t + 0.001) - enthalpy(&chart, t - 0.001)) / 0.002;
            near(cp, dh_dt, 2e-7 * cp);
            let h = material.specific_enthalpy_j_kg();
            let state = law.state_at_specific_enthalpy(h).unwrap();
            near(state.temperature_k, t, 1e-9);
            near(state.temperature_slope_kg_k_j, cp.recip(), 1e-12);
            let dk_dt = if t < 500.0 { 0.001 } else { -0.002 };
            near(state.conductivity_slope_kg_w_m_k_j, dk_dt / cp, 1e-14);
            let plus = law.state_at_specific_enthalpy(h + 1.0).unwrap();
            let minus = law.state_at_specific_enthalpy(h - 1.0).unwrap();
            near(
                state.temperature_slope_kg_k_j,
                0.5 * (plus.temperature_k - minus.temperature_k),
                1e-10,
            );
            near(
                state.conductivity_slope_kg_w_m_k_j,
                0.5 * (plus.conductivity_w_m_k - minus.conductivity_w_m_k),
                1e-12,
            );
            assert_eq!(state.conductivity_w_m_k, law.conductivity_w_m_k(h).unwrap());
        }
        // The existing k(T) chart keeps its left derivative at the upper
        // endpoint, on the same exact original caloric endpoint.
        let state = law
            .state_at_specific_enthalpy(enthalpy(&chart, 650.0))
            .unwrap();
        near(
            state.conductivity_slope_kg_w_m_k_j,
            -0.002 * state.temperature_slope_kg_k_j,
            1e-14,
        );
    }
}

#[test]
fn g0_latent_response_retains_inventory_entropy_and_plateau_endpoint_tangents() {
    let full = chart(101325.0);
    let ts = full.saturation().temperature().value();
    let latent_only = IsobaricWater::new(
        Pressure::new(full.pressure_pa()),
        WaterChartControls {
            minimum_temperature_k: ts,
            maximum_temperature_k: ts,
            ..full.controls
        },
        || false,
    )
    .unwrap();
    for chart in [&full, &latent_only] {
        let law = WaterConductionLaw::try_new(
            chart,
            ConductivityCurve::new(vec![[300.0, 0.2], [650.0, 0.7]]).unwrap(),
        )
        .unwrap();
        let hl = chart.saturation().liquid().specific_enthalpy().value();
        let hv = chart.saturation().vapor().specific_enthalpy().value();
        for quality in [0.0, 0.125, 0.5, 0.875, 1.0] {
            let h = if quality == 1.0 {
                hv
            } else {
                hl + quality * (hv - hl)
            };
            let material = chart
                .state_at_specific_enthalpy(h, || false)
                .unwrap()
                .into_state();
            assert_eq!(
                material.specific_heat_cp_j_kg_k(),
                Err(WaterError::AmbiguousSaturation)
            );
            let WaterHeatCapacityResponse::Latent {
                vapor_mass_fraction,
                latent_enthalpy_j_kg,
            } = material.isobaric_heat_capacity_response()
            else {
                panic!("latent response required")
            };
            near(vapor_mass_fraction, quality, 3e-16);
            assert_eq!(latent_enthalpy_j_kg, hv - hl);
            let state = law.state_at_specific_enthalpy(h).unwrap();
            assert_eq!(state.temperature_k, ts);
            assert_eq!(state.temperature_slope_kg_k_j, 0.0);
            assert_eq!(state.conductivity_slope_kg_w_m_k_j, 0.0);
            let observable = law.time_observables(h).unwrap();
            assert_eq!(observable.temperature_k, ts);
            near(observable.phase_mass_fraction.unwrap(), quality, 3e-16);
        }
        let entropy = law
            .caloric_law()
            .specific_entropy_change_j_kg_k(hl, hv)
            .unwrap()
            .unwrap();
        assert_eq!(entropy, (hv - hl) / ts);
    }
    assert_eq!(latent_only.maximum_temperature_slope(), 0.0);
}

#[test]
fn g4_conductivity_support_and_caloric_domain_are_not_extended() {
    let chart = chart(101325.0);
    for knots in [
        vec![[300.01, 0.5], [650.0, 0.5]],
        vec![[300.0, 0.5], [649.99, 0.5]],
    ] {
        assert!(matches!(
            WaterConductionLaw::try_new(&chart, ConductivityCurve::new(knots).unwrap()),
            Err(DepositionError::MaterialDomain)
        ));
    }
    let law = WaterConductionLaw::try_new(
        &chart,
        ConductivityCurve::new(vec![[300.0, 0.5], [650.0, 0.5]]).unwrap(),
    )
    .unwrap();
    let [lower, upper] = chart.enthalpy_bounds_j_kg();
    for h in [lower - 1.0, upper + 1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            law.state_at_specific_enthalpy(h),
            Err(DepositionError::MaterialDomain)
        );
        assert_eq!(
            law.conductivity_w_m_k(h),
            Err(DepositionError::MaterialDomain)
        );
        assert_eq!(
            law.time_observables(h),
            Err(DepositionError::MaterialDomain)
        );
    }
    for h in [lower, upper] {
        assert!(law.state_at_specific_enthalpy(h).is_ok());
    }
}

fn birth(id: u64, mass_kg: f64, h: f64) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: 1,
            material_id: MaterialId(1),
            mass_kg,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [id as f64, 0.0, 0.0],
                end_m: [id as f64 + 0.001, 0.0, 0.0],
                volume_m3: 1e-4,
                surface_area_m2: 0.02,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        }),
    }
}

#[test]
fn g1_implicit_unequal_mass_exchange_crosses_boiling_and_conserves_original_enthalpy() {
    let chart = chart(101325.0);
    let law = WaterConductionLaw::try_new(
        &chart,
        ConductivityCurve::new(vec![[300.0, 0.5], [650.0, 0.5]]).unwrap(),
    )
    .unwrap();
    let ts = chart.saturation().temperature().value();
    let masses = [0.01, 0.02];
    let old_h = [enthalpy(&chart, ts - 20.0), enthalpy(&chart, ts + 100.0)];
    let hot_final = ts + 40.0;
    // Select the final steam state independently from source h(T), then fix
    // the interval from the BE Fourier law. The cold cell's final inventory
    // follows the remaining energy, without solving another thermal system.
    let heat = masses[1] * (old_h[1] - enthalpy(&chart, hot_final));
    let cold_h = old_h[0] + heat / masses[0];
    let hl = chart.saturation().liquid().specific_enthalpy().value();
    let hv = chart.saturation().vapor().specific_enthalpy().value();
    assert!(cold_h > hl && cold_h < hv);
    let conductance = 0.01 / (0.001 / 0.5 + 0.001 / 0.5);
    let dt = heat / (conductance * (hot_final - ts));
    let initial = DepositionSimulation::new(
        vec![
            birth(1, masses[0], old_h[0]),
            birth(2, masses[1], old_h[1]),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.01,
                    conductance_w_k: 0.0,
                }),
            },
        ],
        &[(MaterialId(1), &chart)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 128,
        },
    )
    .unwrap();
    let (initial, _) = initial
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: dt,
                energy_tolerance_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap();
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(1), &law as &dyn ConductivityLaw)],
        &[ConductiveLink {
            link_id: 1,
            half_distances_m: [0.001, 0.001],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap();
    let implicit = ImplicitConductionChain::new(
        &initial,
        &observer,
        &[(MaterialId(1), &law as &dyn ImplicitConductionLaw)],
        || Ok(()),
    )
    .unwrap();
    let before = initial.cells().clone();
    let controls = ImplicitConductionControls {
        maximum_dt_s: dt,
        cell_tolerance_j: 1e-7,
        total_tolerance_j: 1e-6,
        ..ImplicitConductionControls::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, report) = implicit
        .advance_to(&initial, dt, controls, &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    near(next.temperature_k(CellId(1)).unwrap(), ts, 0.0);
    near(next.temperature_k(CellId(2)).unwrap(), hot_final, 1e-8);
    for i in 0..2 {
        let actual = next.cells()[&CellId(i as u64 + 1)].specific_h_j_kg;
        near(
            masses[i] * (actual - old_h[i]),
            if i == 0 { heat } else { -heat },
            1e-6,
        );
    }
    let observed = law
        .time_observables(next.cells()[&CellId(1)].specific_h_j_kg)
        .unwrap();
    near(
        observed.phase_mass_fraction.unwrap(),
        (cold_h - hl) / (hv - hl),
        1e-10,
    );
    assert_eq!(report.energy.external_heat_j, 0.0);
    assert!(report.energy.residual_j.abs() <= controls.total_tolerance_j);
    assert!(report.entropy.thermal_storage_change_j_k > 0.0);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    near(next.energy_j(), initial.energy_j(), 1e-6);
    assert_eq!(next.mass_kg(), initial.mass_kg());
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);

    // Identical copied source data is not the actual registered caloric owner.
    let other_chart = chart;
    let other_law = WaterConductionLaw::try_new(&other_chart, law.conductivity().clone()).unwrap();
    assert!(
        ImplicitConductionChain::new(
            &initial,
            &observer,
            &[(MaterialId(1), &other_law as &dyn ImplicitConductionLaw)],
            || Ok(())
        )
        .is_err()
    );
}
