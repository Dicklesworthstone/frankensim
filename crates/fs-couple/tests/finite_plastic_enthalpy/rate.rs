//! Source-bound finite Perzyna flow and heat use one physical interval per
//! accepted pair, independently of thermal substeps and nonlinear retries.

use super::*;
use fs_material::plastic::finite::thermal::source::THERMAL_PERZYNA_FINITE_J2_LAW_ID;

const MU: f64 = 80.0e6;
const BULK: f64 = 160.0e6;
const YIELD: f64 = 2.0e6;
const HARDENING: f64 = 8.0e6;
const VISCOSITY: f64 = 8.0e6;
const STRESS_SCALE: f64 = 2.0e6;
const ACTIVATION: f64 = 6000.0;

fn fixture(exponent: f64) -> Fixture {
    let mut source = model();
    source.law = LawId(THERMAL_PERZYNA_FINITE_J2_LAW_ID.into());
    // Isolate viscosity's thermal feedback. Fixed-history mechanical entropy
    // is exactly zero; caloric heating must equal physical flow dissipation.
    for name in [
        "shear_modulus_slope_pa_k",
        "bulk_modulus_slope_pa_k",
        "yield_stress_slope_pa_k",
        "hardening_modulus_slope_pa_k",
        "logarithmic_expansion_per_k",
        "logarithmic_expansion_slope_per_k2",
    ] {
        source.parameters.get_mut(name).unwrap().value = 0.0;
    }
    for (name, value, dims) in [
        (
            "overstress_reference_viscosity_pa_s",
            VISCOSITY,
            Dims([-1, 1, -1, 0, 0, 0]),
        ),
        ("overstress_stress_scale_pa", STRESS_SCALE, PA),
        ("overstress_exponent", exponent, Dims::NONE),
        ("overstress_activation_temperature_k", ACTIVATION, K),
    ] {
        source
            .parameters
            .insert(name.into(), LawParameter { value, dims });
    }
    Fixture::with_model(CP, 10, scalar(DENSITY_PROPERTY, Density::DIMS, RHO), source)
}

fn adiabatic(fixture: &Fixture) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                material_id: MaterialId(0),
                segment_id: 1,
                body_id: 7,
                mass_kg: RHO * V,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [L; 3],
                    volume_m3: V,
                    surface_area_m2: 3.0 * L * L,
                },
                boundary: boundary(0.0),
            }),
        }],
        &[(MaterialId(0), fixture.caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

#[derive(Clone, Copy)]
struct Point {
    a: f64,
    plastic_log: f64,
    alpha: f64,
    temperature: f64,
    free: f64,
    piola: f64,
    dissipation: f64,
    potential: f64,
}
impl Point {
    fn virgin() -> Self {
        Self {
            a: 0.0,
            plastic_log: 0.0,
            alpha: 0.0,
            temperature: 300.0,
            free: 0.0,
            piola: 0.0,
            dissipation: 0.0,
            potential: 0.0,
        }
    }
}

// F=diag(exp(a),1,1), Cp^-1=diag(exp(-2p),exp(p),exp(p)). This
// scalar diagonal oracle never calls a native material/thermal return.
fn mechanical(old: Point, a: f64, temperature: f64, dt: f64, exponent: f64) -> Point {
    let eta = VISCOSITY * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp();
    let trial = 2.0 * MU * a - 3.0 * MU * old.plastic_log;
    let excess = (trial.abs() - YIELD - HARDENING * old.alpha).max(0.0);
    let drag = |gamma: f64| STRESS_SCALE * (eta * gamma / (dt * STRESS_SCALE)).powf(1.0 / exponent);
    let (mut lo, mut hi) = (0.0, excess / (3.0 * MU + HARDENING));
    for _ in 0..100 {
        let gamma = 0.5 * (lo + hi);
        if (3.0 * MU + HARDENING) * gamma + drag(gamma) < excess {
            lo = gamma;
        } else {
            hi = gamma;
        }
    }
    let gamma = 0.5 * (lo + hi);
    let plastic_log = old.plastic_log + trial.signum() * gamma;
    let alpha = old.alpha + gamma;
    let deviator = 2.0 * a / 3.0 - plastic_log;
    let free =
        0.5 * BULK * a * a + 1.5 * MU * deviator * deviator + 0.5 * HARDENING * alpha * alpha;
    let viscous = drag(gamma) * gamma;
    Point {
        a,
        plastic_log,
        alpha,
        temperature,
        free,
        piola: (BULK * a + 2.0 * MU * deviator) / a.exp(),
        dissipation: YIELD * gamma + viscous,
        potential: free + YIELD * gamma + exponent / (exponent + 1.0) * viscous,
    }
}

fn adiabatic_endpoint(old: Point, a: f64, dt: f64, exponent: f64) -> Point {
    let (mut lo, mut hi) = (old.temperature, old.temperature + 30.0);
    for _ in 0..100 {
        let temperature = 0.5 * (lo + hi);
        let point = mechanical(old, a, temperature, dt, exponent);
        if RHO * CP * (temperature - old.temperature) < point.dissipation {
            lo = temperature;
        } else {
            hi = temperature;
        }
    }
    mechanical(old, a, 0.5 * (lo + hi), dt, exponent)
}

#[test]
fn source_bound_power_flow_heats_and_relaxes_once_per_interval_across_feedback_retries() {
    context(|cx| {
        for exponent in [1.0, 3.0] {
            let fixture = fixture(exponent);
            let mut thermal = adiabatic(&fixture);
            let mut solid = bind(&fixture, &thermal, false, cx);
            let prescribed: Vec<_> = (0..12)
                .filter(|dof| *dof != 3)
                .map(|dof| (dof, 0.0))
                .collect();
            let mut old = Point::virgin();
            let mut time = 0.0;
            for (index, dt) in [0.02, 0.02, 0.05].into_iter().enumerate() {
                let expected = adiabatic_endpoint(old, 0.16, dt, exponent);
                let load = [(3, V / L * expected.piola)];
                let end = time + dt;
                let old_history = solid.mechanics().element_states().to_vec();
                let refused = FiniteThermoplasticConfig {
                    maximum_iterations: 1,
                    ..config()
                };
                assert!(
                    solid
                        .advance(
                            &thermal,
                            &(),
                            end,
                            &prescribed,
                            &load,
                            [0.0; 3],
                            settings(),
                            controls(),
                            thermal_config(),
                            refused,
                            cx
                        )
                        .is_err()
                );
                close(thermal.time_s(), time, 0.0);
                close(solid.time_s(), time, 0.0);
                assert_eq!(solid.mechanics().element_states(), old_history.as_slice());

                let candidate = solid
                    .advance(
                        &thermal,
                        &(),
                        end,
                        &prescribed,
                        &load,
                        [0.0; 3],
                        settings(),
                        controls(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let actual = candidate.solid().mechanics();
                let report = candidate.report();
                let history = &actual.element_states()[0];
                close(
                    candidate.thermal().temperature_k(CELL).unwrap(),
                    expected.temperature,
                    3.0e-7,
                );
                close(actual.temperatures_k()[0], expected.temperature, 3.0e-7);
                close(
                    actual.displacement_m()[1][0],
                    L * (expected.a.exp() - 1.0),
                    3.0e-11,
                );
                close(history.alpha, expected.alpha, 3.0e-9);
                for (actual, expected) in history.inverse_plastic_cauchy_green.iter().zip([
                    (-2.0 * expected.plastic_log).exp(),
                    expected.plastic_log.exp(),
                    expected.plastic_log.exp(),
                    0.0,
                    0.0,
                    0.0,
                ]) {
                    close(*actual, expected, 3.0e-9);
                }
                close(actual.element_first_piola_pa()[0][0], expected.piola, 0.5);
                close(actual.free_energy_j(), V * expected.free, 3.0e-9);
                close(
                    report.material_energy.plastic_dissipation_j,
                    V * expected.dissipation,
                    3.0e-9,
                );
                close(
                    report.material_energy.caloric_heat_j,
                    V * expected.dissipation,
                    3.0e-9,
                );
                close(report.material_energy.reversible_heat_j, 0.0, 0.0);
                close(
                    report.material_energy.incremental_potential_j,
                    V * expected.potential,
                    3.0e-9,
                );
                close(
                    report.thermal_audit.stored_energy_change_j,
                    RHO * V * CP * (expected.temperature - old.temperature),
                    3.0e-9,
                );
                close(report.external_heat_j, 0.0, 1.0e-12);
                close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
                close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                assert!(report.iterations > 1);
                assert!(expected.temperature > old.temperature);
                assert!(history.alpha > old.alpha);
                assert!(
                    report.material_energy.incremental_potential_j
                        < actual.free_energy_j() + report.material_energy.plastic_dissipation_j
                            - 1.0e-5
                );
                let work = V * expected.piola * (expected.a.exp() - old.a.exp());
                let defect = work - V * (expected.free - old.free + expected.dissipation);
                close(
                    report.material_energy.internal_endpoint_work_j,
                    work,
                    3.0e-9,
                );
                close(
                    report.material_energy.endpoint_work_defect_j,
                    defect,
                    3.0e-9,
                );
                // Numerical endpoint loss is audited but never added to caloric heat.
                assert!(defect.abs() > 1.0e-5);
                close(
                    work,
                    report.material_energy.internal_energy_change_j
                        + report.thermal_audit.stored_energy_change_j
                        + defect,
                    4.0e-9,
                );
                let material = &actual.families()[0];
                assert_eq!(
                    material.query().card.content_hash(),
                    fixture.card.content_hash()
                );
                assert_eq!(material.family().model_card(), &fixture.card.models()[0]);
                assert!(std::ptr::eq(material.caloric_source(), &fixture.caloric));
                close(material.reference_density_kg_m3(), RHO, 0.0);
                assert!(report.inertia.is_none());
                if index == 0 {
                    let frozen_temperature = mechanical(old, 0.16, old.temperature, dt, exponent);
                    assert!(expected.alpha > frozen_temperature.alpha + 1.0e-5);
                    // Thermal integration may subdivide the same interval. It
                    // must not repeatedly integrate or shorten the plastic clock.
                    let fine_thermal = StepConfig {
                        max_dt_s: 0.002,
                        ..thermal_config()
                    };
                    let subdivided = solid
                        .advance(
                            &thermal,
                            &(),
                            end,
                            &prescribed,
                            &load,
                            [0.0; 3],
                            settings(),
                            controls(),
                            fine_thermal,
                            config(),
                            cx,
                        )
                        .unwrap();
                    close(
                        subdivided.solid().mechanics().element_states()[0].alpha,
                        history.alpha,
                        3.0e-10,
                    );
                    close(
                        subdivided.thermal().temperature_k(CELL).unwrap(),
                        expected.temperature,
                        3.0e-7,
                    );
                } else {
                    assert!(expected.piola < old.piola);
                    close(report.material_energy.internal_endpoint_work_j, 0.0, 3.0e-9);
                }
                assert_eq!(solid.mechanics().element_states(), old_history.as_slice());
                let (next_thermal, next_solid, (), _) = candidate.into_parts();
                thermal = next_thermal;
                solid = next_solid;
                old = expected;
                time = end;
            }
        }
    });
}
