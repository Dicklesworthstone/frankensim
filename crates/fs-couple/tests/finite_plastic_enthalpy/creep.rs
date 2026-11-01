//! Source-bound, sub-yield finite Norton creep with free lateral contraction.
//! The mechanical/caloric oracle is independent of the production return map.

use super::*;
use fs_material::plastic::finite::thermal::source::THERMAL_FINITE_NORTON_LAW_ID;

const MU: f64 = 80.0e6;
const BULK: f64 = 160.0e6;
const VISCOSITY: f64 = 30.0e6;
const STRESS_SCALE: f64 = 1.0e6;
const ACTIVATION: f64 = 6000.0;
const DEAD_PIOLA: f64 = 1.0e6;
const COMPARISON_YIELD: f64 = 2.0e6;

fn fixture(exponent: f64) -> Fixture {
    let mut source = model();
    source.law = LawId(THERMAL_FINITE_NORTON_LAW_ID.into());
    for name in [
        "yield_stress_pa",
        "hardening_modulus_pa",
        "yield_stress_slope_pa_k",
        "hardening_modulus_slope_pa_k",
    ] {
        source.parameters.remove(name).unwrap();
    }
    // Viscosity supplies the only thermal feedback in this fixture. Mechanical
    // entropy is zero, and the original source caloric chart receives q*dalpha.
    for name in [
        "shear_modulus_slope_pa_k",
        "bulk_modulus_slope_pa_k",
        "logarithmic_expansion_per_k",
        "logarithmic_expansion_slope_per_k2",
    ] {
        source.parameters.get_mut(name).unwrap().value = 0.0;
    }
    for (name, value, dims) in [
        (
            "creep_reference_viscosity_pa_s",
            VISCOSITY,
            Dims([-1, 1, -1, 0, 0, 0]),
        ),
        ("creep_stress_scale_pa", STRESS_SCALE, PA),
        ("creep_exponent", exponent, Dims::NONE),
        ("creep_activation_temperature_k", ACTIVATION, K),
    ] {
        source
            .parameters
            .insert(name.into(), LawParameter { value, dims });
    }
    Fixture::with_model(CP, 11, scalar(DENSITY_PROPERTY, Density::DIMS, RHO), source)
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
    b: f64,
    creep: f64,
    temperature: f64,
    q: f64,
    piola: [f64; 3],
    free: f64,
    dissipation: f64,
    potential: f64,
}
impl Point {
    fn virgin() -> Self {
        Self {
            a: 0.0,
            b: 0.0,
            creep: 0.0,
            temperature: 300.0,
            q: 0.0,
            piola: [0.0; 3],
            free: 0.0,
            dissipation: 0.0,
            potential: 0.0,
        }
    }
}
fn flow(q: f64, temperature: f64, dt: f64, exponent: f64) -> f64 {
    let eta = VISCOSITY * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp();
    dt * STRESS_SCALE / eta * (q / STRESS_SCALE).powf(exponent)
}
fn point(old: Point, a: f64, b: f64, creep: f64, temperature: f64, exponent: f64) -> Point {
    let g = a + 2.0 * b;
    let e = 2.0 * (a - b) / 3.0 - creep;
    let q = 3.0 * MU * e;
    let free = 0.5 * BULK * g * g + 1.5 * MU * e * e;
    let dissipation = q * (creep - old.creep);
    Point {
        a,
        b,
        creep,
        temperature,
        q,
        piola: [
            (BULK * g + 2.0 * MU * e) / a.exp(),
            (BULK * g - MU * e) / b.exp(),
            (BULK * g - MU * e) / b.exp(),
        ],
        free,
        dissipation,
        potential: free + exponent / (exponent + 1.0) * dissipation,
    }
}
// Actual uniaxial dead-load equilibrium: q=P*exp(a), with two free
// lateral contractions. Solve the scalar equation on its stable tensile branch.
fn mechanical(old: Point, temperature: f64, dt: f64, exponent: f64, hold: bool) -> Point {
    if hold {
        let (mut lo, mut hi) = (0.0, old.q);
        for _ in 0..100 {
            let q = 0.5 * (lo + hi);
            if q + 3.0 * MU * flow(q, temperature, dt, exponent) < old.q {
                lo = q;
            } else {
                hi = q;
            }
        }
        let q = 0.5 * (lo + hi);
        return point(
            old,
            old.a,
            old.b,
            old.creep + flow(q, temperature, dt, exponent),
            temperature,
            exponent,
        );
    }
    let compliance = 1.0 / (9.0 * BULK) + 1.0 / (3.0 * MU);
    let residual = |q: f64| {
        (q / DEAD_PIOLA).ln() - old.creep - flow(q, temperature, dt, exponent) - compliance * q
    };
    let (mut lo, mut hi) = (DEAD_PIOLA, 2.0 * DEAD_PIOLA);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..100 {
        let q = 0.5 * (lo + hi);
        if residual(q) < 0.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let creep = old.creep + flow(q, temperature, dt, exponent);
    let a = creep + compliance * q;
    let b = q / (9.0 * BULK) - q / (6.0 * MU) - 0.5 * creep;
    point(old, a, b, creep, temperature, exponent)
}
fn adiabatic_endpoint(old: Point, dt: f64, exponent: f64, hold: bool) -> Point {
    let (mut lo, mut hi) = (old.temperature, old.temperature + 4.0);
    for _ in 0..100 {
        let temperature = 0.5 * (lo + hi);
        let state = mechanical(old, temperature, dt, exponent, hold);
        if RHO * CP * (temperature - old.temperature) < state.dissipation {
            lo = temperature;
        } else {
            hi = temperature;
        }
    }
    mechanical(old, 0.5 * (lo + hi), dt, exponent, hold)
}
fn supports(old: Point, hold: bool) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|dof| hold || ![3, 7, 11].contains(dof))
        .map(|dof| {
            (
                dof,
                match dof {
                    3 => L * (old.a.exp() - 1.0),
                    7 | 11 => L * (old.b.exp() - 1.0),
                    _ => 0.0,
                },
            )
        })
        .collect()
}

#[test]
fn sourced_norton_body_creeps_below_yield_heats_and_relaxes_with_physical_clock_and_atomic_retries()
{
    context(|cx| {
        for exponent in [1.0, 3.0] {
            let fixture = fixture(exponent);
            let mut thermal = adiabatic(&fixture);
            let mut solid = bind(&fixture, &thermal, false, cx);
            let mut old = Point::virgin();
            let mut time = 0.0;
            for (index, dt) in [1.0, 0.75, 1.0, 1.0].into_iter().enumerate() {
                let hold = index == 3;
                if hold {
                    solid = solid.with_inertia(&[[1.0, 0.0, 0.0]; 4], cx).unwrap();
                }
                let expected = adiabatic_endpoint(old, dt, exponent, hold);
                assert!(expected.q > 0.0 && expected.q < COMPARISON_YIELD);
                assert!(expected.creep > old.creep);
                let prescribed = supports(old, hold);
                let load = if hold {
                    vec![]
                } else {
                    vec![(3, V / L * DEAD_PIOLA)]
                };
                let end = time + dt;
                let saved_history = solid.mechanics().element_states().to_vec();
                let saved_u = solid.mechanics().displacement_m().to_vec();
                let saved_h = thermal.cells()[&CELL].specific_h_j_kg;
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
                            FiniteThermoplasticConfig {
                                maximum_iterations: 1,
                                ..config()
                            },
                            cx
                        )
                        .is_err()
                );
                assert_eq!(solid.mechanics().element_states(), saved_history.as_slice());
                assert_eq!(solid.mechanics().displacement_m(), saved_u.as_slice());
                close(thermal.cells()[&CELL].specific_h_j_kg, saved_h, 0.0);
                close(thermal.time_s(), time, 0.0);
                close(solid.time_s(), time, 0.0);
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
                    2e-7,
                );
                close(actual.temperatures_k()[0], expected.temperature, 2e-7);
                close(history.alpha, expected.creep, 3e-9);
                for (actual, expected) in history.inverse_plastic_cauchy_green.iter().zip([
                    (-2.0 * expected.creep).exp(),
                    expected.creep.exp(),
                    expected.creep.exp(),
                    0.0,
                    0.0,
                    0.0,
                ]) {
                    close(*actual, expected, 3e-9);
                }
                for (axis, logarithm) in
                    [expected.a, expected.b, expected.b].into_iter().enumerate()
                {
                    close(
                        actual.displacement_m()[axis + 1][axis],
                        L * (logarithm.exp() - 1.0),
                        3e-11,
                    );
                    close(
                        actual.element_first_piola_pa()[0][4 * axis],
                        expected.piola[axis],
                        0.5,
                    );
                }
                close(actual.free_energy_j(), V * expected.free, 3e-9);
                close(
                    report.material_energy.plastic_dissipation_j,
                    V * expected.dissipation,
                    3e-9,
                );
                close(
                    report.material_energy.caloric_heat_j,
                    V * expected.dissipation,
                    3e-9,
                );
                close(report.material_energy.reversible_heat_j, 0.0, 0.0);
                close(
                    report.material_energy.incremental_potential_j,
                    V * expected.potential,
                    3e-9,
                );
                close(
                    report.thermal_audit.stored_energy_change_j,
                    RHO * V * CP * (expected.temperature - old.temperature),
                    3e-9,
                );
                close(report.external_heat_j, 0.0, 1e-12);
                close(report.thermal_residual_j, 0.0, config().energy_tolerance_j);
                close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                assert!(report.iterations > 1 && expected.temperature > old.temperature);
                let work = V
                    * (expected.piola[0] * (expected.a.exp() - old.a.exp())
                        + (expected.piola[1] + expected.piola[2])
                            * (expected.b.exp() - old.b.exp()));
                let defect = work - V * (expected.free - old.free + expected.dissipation);
                close(report.material_energy.internal_endpoint_work_j, work, 3e-9);
                close(report.material_energy.endpoint_work_defect_j, defect, 3e-9);
                assert!(defect.abs() > 1e-8);
                close(
                    work,
                    report.material_energy.internal_energy_change_j
                        + report.thermal_audit.stored_energy_change_j
                        + defect,
                    4e-9,
                );
                // Stored elastic energy, physical creep heat and integrated
                // flow potential remain separately observable at every step.
                assert!(expected.potential < expected.free + expected.dissipation);
                let material = &actual.families()[0];
                assert!(material.family().law().is_norton_creep());
                assert_eq!(material.family().model_card(), &fixture.card.models()[0]);
                assert!(std::ptr::eq(material.caloric_source(), &fixture.caloric));
                close(material.reference_density_kg_m3(), RHO, 0.0);
                if index == 0 {
                    let cold = mechanical(old, old.temperature, dt, exponent, false);
                    assert!(expected.creep > cold.creep + 1e-6);
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
                            StepConfig {
                                max_dt_s: 0.05,
                                ..thermal_config()
                            },
                            config(),
                            cx,
                        )
                        .unwrap();
                    close(
                        subdivided.solid().mechanics().element_states()[0].alpha,
                        history.alpha,
                        3e-10,
                    );
                    close(
                        subdivided.thermal().temperature_k(CELL).unwrap(),
                        expected.temperature,
                        2e-7,
                    );
                }
                if hold {
                    assert!(expected.q < old.q);
                    close(work, 0.0, 0.0);
                    let inertia = report.inertia.as_ref().unwrap();
                    close(inertia.numerical_loss_j, 0.5 * RHO * V, 1e-11);
                    close(inertia.kinetic_change_j, -inertia.numerical_loss_j, 1e-11);
                    assert!(inertia.numerical_loss_j > 1e-5);
                } else {
                    assert!(report.inertia.is_none());
                    assert!(expected.a > old.a && expected.b < old.b);
                    close(
                        report.material_energy.dead_load_endpoint_work_j,
                        V * DEAD_PIOLA * (expected.a.exp() - old.a.exp()),
                        3e-9,
                    );
                }
                assert_eq!(solid.mechanics().element_states(), saved_history.as_slice());
                let (next_thermal, next_solid, (), _) = candidate.into_parts();
                thermal = next_thermal;
                solid = next_solid;
                old = expected;
                time = end;
            }
            assert!(solid.mechanics().displacement_m()[1][0] > 0.09 * L);
        }
    });
}
