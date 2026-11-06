//! Independent held-stress heat ODE and original moving-liquid publication.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceConductionControls, FreeSurfaceConductionFace,
};

const DT: f64 = 0.001;
const RHO: f64 = 1000.0;
const CP: f64 = 2000.0;

fn controls(substeps: usize) -> FreeSurfaceConductionControls {
    FreeSurfaceConductionControls {
        maximum_faces: 1,
        maximum_substeps: substeps + 2,
        maximum_dt_s: DT / substeps as f64,
        energy_tolerance_j: 2e-8,
    }
}

fn cells() -> impl Iterator<Item = [usize; 3]> {
    (0..8).flat_map(|z| (0..8).flat_map(move |y| (0..8).map(move |x| [x, y, z])))
}

fn cell_mass(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>, p: [usize; 3]) -> f64 {
    liquid
        .flow()
        .region_mass_kg(|x, y, z| [x, y, z] == p)
        .unwrap()
}

// RK4 integrates the dimensional scalar closure, independently of the native
// thermal integrator and constitutive root. The second component integrates
// molecular heat directly; subtracting a much larger authored heat would lose it.
fn heat_ode(
    numerator: f64,
    dilation_numerator: f64,
    initial_temperature: f64,
    temperature_slope: f64,
    shear_slope: f64,
    authored_w_m3: f64,
    independent_bulk: Option<(f64, f64)>,
) -> (f64, f64) {
    let rhs = |t| {
        let gamma = independent_rate(numerator, t, temperature_slope, shear_slope);
        let mu = 0.1 + temperature_slope * (t - 400.0) + shear_slope * gamma;
        let bulk_power = match independent_bulk {
            Some((zeta, bulk_tau)) => {
                let divergence = 1.5_f64.sqrt() * dilation_numerator / (bulk_tau * DT);
                zeta * divergence * divergence
            }
            None => {
                let dilation = dilation_numerator / ((0.5 + 3.0 * mu) * DT);
                mu * dilation * dilation
            }
        };
        let q = mu * gamma * gamma + bulk_power;
        [(authored_w_m3 + q) / (RHO * CP), q]
    };
    let dt = DT / 1024.0;
    let mut t = initial_temperature;
    let mut heat_density = 0.0;
    for _ in 0..1024 {
        let a = rhs(t);
        let b = rhs(t + 0.5 * dt * a[0]);
        let c = rhs(t + 0.5 * dt * b[0]);
        let d = rhs(t + dt * c[0]);
        t += dt * (a[0] + 2.0 * b[0] + 2.0 * c[0] + d[0]) / 6.0;
        heat_density += dt * (a[1] + 2.0 * b[1] + 2.0 * c[1] + d[1]) / 6.0;
    }
    (t, heat_density)
}

#[test]
fn g1_nonlinear_temperature_and_shear_feedback_converges_to_independent_heat_ode() {
    let temperature_slope = -0.0008;
    let shear_slope = 0.003;
    let source = rheological_card(temperature_slope, shear_slope, 0.0, 1000.0);
    let gravity = [0.03, -0.02, 0.01];
    with_liquid(
        &source,
        true,
        true,
        gravity,
        solve(),
        |_| 200000.0,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            let mut authored = Vec::new();
            let mut exact_heat = 0.0;
            let mut exact_cell_t = 0.0;
            let mut cell_numerator = 0.0;
            let mut cell_dilation_numerator = 0.0;
            for p in cells() {
                let mass = cell_mass(liquid, p);
                let [numerator, dilation_numerator] = stress_numerators(liquid, p, gravity);
                let initial_t = coordinate(liquid.material_at(p).unwrap(), "T");
                let (t, q_integral) = heat_ode(
                    numerator,
                    dilation_numerator,
                    initial_t,
                    temperature_slope,
                    shear_slope,
                    RHO * CP * 60.0 / DT,
                    None,
                );
                exact_heat += q_integral * mass / RHO;
                authored.push(FreeSurfaceCellHeat {
                    coordinates: p,
                    energy_j: mass * CP * 60.0,
                });
                if p == CELL {
                    exact_cell_t = t;
                    cell_numerator = numerator;
                    cell_dilation_numerator = dilation_numerator;
                }
            }
            assert!(exact_heat > 0.0);
            let mut errors = Vec::new();
            for intervals in [8, 32] {
                let mut trial = liquid.clone();
                let report = trial
                    .step_with_conduction_and_viscous_feedback(
                        &[],
                        &authored,
                        controls(intervals),
                        gate,
                    )
                    .unwrap();
                let molecular = report.viscous_feedback.as_ref().unwrap();
                let error = (molecular.heat_j / exact_heat - 1.0).abs();
                eprintln!(
                    "viscous thermal feedback intervals={intervals} relative_heat_error={error:e}"
                );
                errors.push(error);
                assert_eq!(molecular.heat_j, report.flow.viscous_heat_j);
                near(
                    report.flow.external_heat_j,
                    authored.iter().map(|h| h.energy_j).sum(),
                    2e-10,
                );
                assert!(report.flow.balance_defect_j.abs() < 2e-8);
                near(report.flow.mass_after_kg, report.flow.mass_before_kg, 1e-12);
                assert_eq!(molecular.last_applied.len(), 512);
                assert_eq!(molecular.thermal_endpoint.len(), 512);
                let endpoint = molecular
                    .thermal_endpoint
                    .iter()
                    .find(|s| s.coordinates == CELL)
                    .unwrap();
                let temperature = coordinate(&endpoint.material, "T");
                assert!((temperature - exact_cell_t).abs() < 2e-7);
                let expected_rate =
                    independent_rate(cell_numerator, temperature, temperature_slope, shear_slope);
                near(endpoint.shear_rate_per_s, expected_rate, 2e-10);
                near(coordinate(&endpoint.material, "rate"), expected_rate, 2e-10);
                let mu = endpoint.material.dynamic_viscosity_pa_s();
                let dilation = cell_dilation_numerator / ((0.5 + 3.0 * mu) * DT);
                near(
                    endpoint.power_w_m3,
                    mu * (expected_rate * expected_rate + dilation * dilation),
                    2e-9,
                );
                near(
                    endpoint.collision_tau,
                    0.5 + 3.0 * endpoint.material.dynamic_viscosity_pa_s(),
                    2e-12,
                );
                assert_eq!(
                    endpoint.material.material().card_identity(),
                    source.content_hash()
                );
                // The thermal snapshot is before the original VOF advection; the
                // subsequent owner still publishes exactly one physical flow step.
                assert_eq!(trial.accepted_steps(), liquid.accepted_steps() + 1);
                near(
                    trial.physical_time_s() - liquid.physical_time_s(),
                    DT,
                    1e-15,
                );
            }
            assert!(errors[1] < 0.35 * errors[0], "{errors:?}");
            assert!(errors[1] < 0.01, "{errors:?}");
            let mut held = liquid.clone();
            let frozen = held
                .step_with_conduction_and_viscous_heating(&[], &authored, controls(32), gate)
                .unwrap();
            assert!(
                (frozen.conduction.flow.viscous_heat_j / exact_heat - 1.0).abs() > 4.0 * errors[0]
            );
        },
    );
}

#[test]
fn independent_bulk_feedback_keeps_volume_relaxation_when_shear_viscosity_changes() {
    let temperature_slope = -0.0008;
    let shear_slope = 0.003;
    let source = rheological_card(temperature_slope, shear_slope, 0.0, 1000.0);
    let gravity = [0.03, -0.02, 0.01];
    for zeta in [0.02, 0.24] {
        with_liquid(
            &source, true, true, gravity, solve(), |_| 200000.0,
            |liquid, gate| {
                liquid.set_bulk_viscosity(zeta, 1.2).unwrap();
                let model = liquid.flow().lattice().collision_model();
                let original_bits = liquid.flow().lattice().grid().state_bits();
                for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                    assert!(liquid.set_bulk_viscosity(invalid, 1.2).is_err());
                    assert_eq!(liquid.flow().lattice().collision_model(), model);
                    assert_eq!(liquid.flow().lattice().grid().state_bits(), original_bits);
                    assert_eq!(liquid.flow().scaling().bulk_viscosity_pa_s(), Some(zeta));
                }
                liquid.step(gate).unwrap();
                let before = liquid.clone();
                let bulk_tau = 0.5 + 4.5 * zeta;
                let mut authored = Vec::new();
                let mut expected_heat = 0.0;
                let mut initial_bulk_heat = 0.0;
                for p in cells() {
                    let mass = cell_mass(&before, p);
                    let [n, spherical] = stress_numerators(&before, p, gravity);
                    let (_, integral) = heat_ode(
                        n, spherical, coordinate(before.material_at(p).unwrap(), "T"),
                        temperature_slope, shear_slope, RHO * CP * 60.0 / DT,
                        Some((zeta, bulk_tau)),
                    );
                    expected_heat += integral * mass / RHO;
                    let divergence = 1.5_f64.sqrt() * spherical / (bulk_tau * DT);
                    initial_bulk_heat += zeta * divergence * divergence * mass / RHO * DT;
                    authored.push(FreeSurfaceCellHeat { coordinates: p, energy_j: mass * CP * 60.0 });
                }
                assert!(initial_bulk_heat > 1e-12);
                let report = liquid.step_with_conduction_and_viscous_feedback(
                    &[], &authored, controls(32), gate,
                ).unwrap();
                let feedback = report.viscous_feedback.as_ref().unwrap();
                assert!((feedback.heat_j / expected_heat - 1.0).abs() < 0.01);
                let mut endpoint_bulk_heat = 0.0;
                for sample in &feedback.thermal_endpoint {
                    let [n, spherical] = stress_numerators(&before, sample.coordinates, gravity);
                    let temperature = coordinate(&sample.material, "T");
                    let shear = independent_rate(n, temperature, temperature_slope, shear_slope);
                    let divergence = 1.5_f64.sqrt() * spherical / (bulk_tau * DT);
                    let mu = sample.material.dynamic_viscosity_pa_s();
                    near(sample.shear_rate_per_s, shear, 2e-10);
                    near(sample.divergence_magnitude_per_s, divergence, 2e-10);
                    near(sample.shear_power_w_m3, mu * shear * shear, 2e-9);
                    near(sample.bulk_power_w_m3, zeta * divergence * divergence, 2e-9);
                    near(sample.power_w_m3, sample.shear_power_w_m3 + sample.bulk_power_w_m3, 2e-14);
                    endpoint_bulk_heat += sample.bulk_power_w_m3 * sample.liquid_volume_m3 * DT;
                }
                near(endpoint_bulk_heat, initial_bulk_heat, 2e-12);
                assert!(coordinate(&feedback.thermal_endpoint[0].material, "T") > 450.0);
                assert_eq!(liquid.flow().lattice().collision_model(), model);
                assert_eq!(liquid.flow().scaling().bulk_viscosity_pa_s(), Some(zeta));
                assert_eq!(report.flow.viscous_heat_j, feedback.heat_j);
                assert!(report.flow.balance_defect_j.abs() < 2e-8);
                let accepted = liquid.clone();
                let cancelled = CancelGate::new_clock_free();
                cancelled.request();
                assert!(liquid.step_with_conduction_and_viscous_feedback(
                    &[], &[], controls(1), &cancelled,
                ).is_err());
                assert_eq!(liquid.flow().lattice().grid().state_bits(), accepted.flow().lattice().grid().state_bits());
                assert_eq!(liquid.total_enthalpy_j().unwrap(), accepted.total_enthalpy_j().unwrap());
                assert_eq!(liquid.flow().lattice().collision_model(), model);
            },
        );
    }
}

#[test]
fn g1_newtonian_feedback_uses_new_temperature_and_preserves_independent_conduction_balance() {
    let source = card(0.02, 1000.0);
    with_liquid(
        &source,
        false,
        true,
        [0.0; 3],
        solve(),
        |p| 120000.0 + 10000.0 * p[0] as f64,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            let before = liquid.clone();
            let authored = [FreeSurfaceCellHeat {
                coordinates: CELL,
                energy_j: 0.04,
            }];
            let report = liquid
                .step_with_conduction_and_viscous_feedback(
                    &[FreeSurfaceConductionFace {
                        cells: [CELL, [4, 3, 3]],
                        area_m2: 1e-6,
                        distances_m: [0.0005; 2],
                    }],
                    &authored,
                    controls(8),
                    gate,
                )
                .unwrap();
            let feedback = report.viscous_feedback.unwrap();
            assert!(feedback.heat_j > 0.0);
            near(report.flow.external_heat_j, 0.04, 1e-14);
            near(
                report.thermal.external_heat_j,
                0.04 + feedback.heat_j,
                1e-13,
            );
            assert!(report.flow.balance_defect_j.abs() < 2e-8);
            assert_eq!(report.flow.viscous_heat_j, feedback.heat_j);
            for sample in &feedback.thermal_endpoint {
                let [n, volume] = stress_numerators(&before, sample.coordinates, [0.0; 3]);
                let mu = sample.material.dynamic_viscosity_pa_s();
                let rate = n / ((0.5 + 3.0 * mu) * DT);
                let dilation = volume / ((0.5 + 3.0 * mu) * DT);
                near(sample.shear_rate_per_s, rate, 2e-10);
                near(sample.power_w_m3, mu * (rate * rate + dilation * dilation), 2e-9);
                near(
                    sample.liquid_volume_m3,
                    cell_mass(&before, sample.coordinates) / RHO,
                    2e-12,
                );
            }
            let endpoint = feedback
                .thermal_endpoint
                .iter()
                .find(|s| s.coordinates == CELL)
                .unwrap();
            assert!(
                coordinate(&endpoint.material, "T")
                    > coordinate(before.material_at(CELL).unwrap(), "T") + 10.0
            );
            assert!(
                endpoint.material.dynamic_viscosity_pa_s()
                    < before.material_at(CELL).unwrap().dynamic_viscosity_pa_s()
            );
        },
    );
}

#[test]
fn g3_constant_viscosity_recovers_held_heat_without_changing_native_populations() {
    let source = card(0.1, 1000.0);
    with_liquid(
        &source,
        false,
        true,
        [0.0; 3],
        solve(),
        |_| 200000.0,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            let mut frozen = liquid.clone();
            let held = frozen
                .step_with_conduction_and_viscous_heating(&[], &[], controls(4), gate)
                .unwrap();
            let refreshed = liquid
                .step_with_conduction_and_viscous_feedback(&[], &[], controls(4), gate)
                .unwrap();
            near(
                refreshed.flow.viscous_heat_j,
                held.conduction.flow.viscous_heat_j,
                2e-13,
            );
            assert_eq!(
                liquid.flow().lattice().grid().state_bits(),
                frozen.flow().lattice().grid().state_bits()
            );
            near(
                liquid.total_enthalpy_j().unwrap(),
                frozen.total_enthalpy_j().unwrap(),
                1e-11,
            );
            assert!(held.conduction.viscous_feedback.is_none());
            assert!(refreshed.viscous_feedback.is_some());
        },
    );
}

#[test]
fn g4_thermal_source_and_late_budget_refusals_preserve_complete_liquid_and_retry() {
    let source = rheological_card(-0.0008, 0.003, 0.0, 1000.0);
    with_liquid(
        &source,
        true,
        true,
        [0.0; 3],
        solve(),
        |_| 200000.0,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            let original = liquid.clone();
            let bits = liquid.flow().lattice().grid().state_bits();
            let h = liquid.total_enthalpy_j().unwrap();
            let state = liquid.material_at(CELL).unwrap().clone();
            let rheology = *liquid.rheology_at(CELL).unwrap();
            let authored = [FreeSurfaceCellHeat {
                coordinates: CELL,
                energy_j: 0.04,
            }];
            let late_budget = FreeSurfaceConductionControls {
                maximum_substeps: 1,
                ..controls(4)
            };
            assert!(
                liquid
                    .step_with_conduction_and_viscous_feedback(&[], &authored, late_budget, gate)
                    .is_err()
            );
            let unsupported = [FreeSurfaceCellHeat {
                coordinates: CELL,
                energy_j: 0.5,
            }];
            assert!(
                liquid
                    .step_with_conduction_and_viscous_feedback(&[], &unsupported, controls(4), gate)
                    .is_err()
            );
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(
                liquid
                    .step_with_conduction_and_viscous_feedback(
                        &[],
                        &authored,
                        controls(4),
                        &cancelled
                    )
                    .is_err()
            );
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), h);
            assert_eq!(liquid.material_at(CELL), Some(&state));
            assert_eq!(liquid.rheology_at(CELL), Some(&rheology));
            assert_eq!(liquid.physical_time_s(), original.physical_time_s());
            assert_eq!(liquid.accepted_steps(), original.accepted_steps());
            let mut expected = original;
            let a = liquid
                .step_with_conduction_and_viscous_feedback(&[], &authored, controls(4), gate)
                .unwrap();
            let b = expected
                .step_with_conduction_and_viscous_feedback(&[], &authored, controls(4), gate)
                .unwrap();
            assert_eq!(a.flow, b.flow);
            assert_eq!(a.viscous_feedback, b.viscous_feedback);
            assert_eq!(
                liquid.flow().lattice().grid().state_bits(),
                expected.flow().lattice().grid().state_bits()
            );
        },
    );
}
