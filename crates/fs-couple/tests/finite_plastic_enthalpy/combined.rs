//! Simultaneous finite Norton creep and J2 yielding with source-bound self-heating.
//! A separate uniaxial Hencky/complementarity and caloric root supplies the oracle.

use super::*;
use fs_material::plastic::finite::thermal::source::THERMAL_FINITE_NORTON_J2_LAW_ID;

const MU: f64 = 80.0e6;
const BULK: f64 = 160.0e6;
const YIELD: f64 = 2.0e6;
const HARDENING: f64 = 40.0e6;
const VISCOSITY: f64 = 80.0e6;
const STRESS_SCALE: f64 = 4.0e6;
const ACTIVATION: f64 = 6000.0;

fn fixture(exponent: f64) -> Fixture {
    let mut source = model();
    source.law = LawId(THERMAL_FINITE_NORTON_J2_LAW_ID.into());
    source.state_schema_version = 2;
    source
        .parameters
        .get_mut("hardening_modulus_pa")
        .unwrap()
        .value = HARDENING;
    // Isolate Arrhenius feedback: mechanical entropy vanishes, so the original
    // caloric source must receive Y*dp + q*dc, with neither hardening nor defects.
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
            "creep_reference_viscosity_pa_s",
            VISCOSITY,
            Dims([-1, 1, -1, 0, 0, 0]),
        ),
        ("creep_stress_scale_pa", STRESS_SCALE, PA),
        ("creep_exponent", exponent, Dims::NONE),
        ("creep_activation_temperature_k", ACTIVATION, K),
        ("maximum_accumulated_creep_strain", 2.0, Dims::NONE),
    ] {
        source
            .parameters
            .insert(name.into(), LawParameter { value, dims });
    }
    Fixture::with_model(CP, 12, scalar(DENSITY_PROPERTY, Density::DIMS, RHO), source)
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
    plastic: f64,
    creep: f64,
    temperature: f64,
    q: f64,
    piola: [f64; 3],
    free: f64,
    plastic_loss: f64,
    creep_loss: f64,
    potential: f64,
}
impl Point {
    fn virgin() -> Self {
        Self {
            a: 0.0,
            b: 0.0,
            plastic: 0.0,
            creep: 0.0,
            temperature: 300.0,
            q: 0.0,
            piola: [0.0; 3],
            free: 0.0,
            plastic_loss: 0.0,
            creep_loss: 0.0,
            potential: 0.0,
        }
    }
    fn heat(self) -> f64 {
        self.plastic_loss + self.creep_loss
    }
}
fn flow(q: f64, temperature: f64, dt: f64, exponent: f64) -> f64 {
    let eta = VISCOSITY * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp();
    dt * STRESS_SCALE / eta * (q / STRESS_SCALE).powf(exponent)
}
fn point(old: Point, strains: [f64; 4], temperature: f64, exponent: f64) -> Point {
    let [a, b, plastic, creep] = strains;
    let g = a + 2.0 * b;
    let e = 2.0 * (a - b) / 3.0 - plastic - creep;
    let q = 3.0 * MU * e;
    let free = 0.5 * BULK * g * g + 1.5 * MU * e * e + 0.5 * HARDENING * plastic * plastic;
    let plastic_loss = YIELD * (plastic - old.plastic);
    let creep_loss = q * (creep - old.creep);
    Point {
        a,
        b,
        plastic,
        creep,
        temperature,
        q,
        piola: [
            (BULK * g + 2.0 * MU * e) / a.exp(),
            (BULK * g - MU * e) / b.exp(),
            (BULK * g - MU * e) / b.exp(),
        ],
        free,
        plastic_loss,
        creep_loss,
        potential: free + plastic_loss + exponent / (exponent + 1.0) * creep_loss,
    }
}
// Three free coordinates give true uniaxial traction: q=P*exp(a). Combining
// lateral equilibrium, associative coaxial flow and yield complementarity leaves
// one scalar root, independently of the production tensor return/tangent.
fn mechanical(old: Point, temperature: f64, dt: f64, exponent: f64, load: Option<f64>) -> Point {
    let Some(piola) = load else {
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
            [
                old.a,
                old.b,
                old.plastic,
                old.creep + flow(q, temperature, dt, exponent),
            ],
            temperature,
            exponent,
        );
    };
    let compliance = 1.0 / (9.0 * BULK) + 1.0 / (3.0 * MU);
    let plastic = |q: f64| old.plastic.max((q - YIELD) / HARDENING);
    let residual = |q: f64| {
        (q / piola).ln()
            - plastic(q)
            - old.creep
            - flow(q, temperature, dt, exponent)
            - compliance * q
    };
    let (mut lo, mut hi) = (piola, 2.0 * piola);
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
    let plastic = plastic(q);
    let creep = old.creep + flow(q, temperature, dt, exponent);
    let a = plastic + creep + compliance * q;
    let b = q / (9.0 * BULK) - q / (6.0 * MU) - 0.5 * (plastic + creep);
    point(old, [a, b, plastic, creep], temperature, exponent)
}
fn adiabatic_endpoint(old: Point, dt: f64, exponent: f64, load: Option<f64>) -> Point {
    let (mut lo, mut hi) = (old.temperature, old.temperature + 4.0);
    let residual = |temperature| {
        RHO * CP * (temperature - old.temperature)
            - mechanical(old, temperature, dt, exponent, load).heat()
    };
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..100 {
        let temperature = 0.5 * (lo + hi);
        if residual(temperature) < 0.0 {
            lo = temperature;
        } else {
            hi = temperature;
        }
    }
    mechanical(old, 0.5 * (lo + hi), dt, exponent, load)
}
fn supports(old: Point, hold: bool) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|dof| hold || ![3, 7, 11].contains(dof))
        .map(|dof| {
            (
                dof,
                match dof {
                    3 => L * old.a.exp_m1(),
                    7 | 11 => L * old.b.exp_m1(),
                    _ => 0.0,
                },
            )
        })
        .collect()
}

#[test]
fn sourced_simultaneous_creep_yield_heats_then_retains_hardening_during_creep_and_relaxation() {
    context(|cx| {
        for exponent in [1.0, 3.0] {
            let fixture = fixture(exponent);
            let mut thermal = adiabatic(&fixture);
            let mut solid = bind(&fixture, &thermal, false, cx);
            let mut old = Point::virgin();
            let mut time = 0.0;
            for (index, (dt, load)) in [
                (0.5, Some(4.0e6)),
                (0.5, Some(4.0e6)),
                (1.0, Some(1.5e6)),
                (1.0, None),
            ]
            .into_iter()
            .enumerate()
            {
                let hold = load.is_none();
                if hold {
                    solid = solid.with_inertia(&[[1.0, 0.0, 0.0]; 4], cx).unwrap();
                }
                let expected = adiabatic_endpoint(old, dt, exponent, load);
                let prescribed = supports(old, hold);
                let nodal_load: Vec<_> = load.map(|p| (3, V / L * p)).into_iter().collect();
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
                            &nodal_load,
                            [0.0; 3],
                            settings(),
                            controls(),
                            thermal_config(),
                            FiniteThermoplasticConfig {
                                maximum_iterations: 1,
                                ..config()
                            },
                            cx,
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
                        &nodal_load,
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
                close(history.alpha, expected.plastic, 3e-9);
                close(history.creep_alpha, expected.creep, 3e-9);
                let inelastic = expected.plastic + expected.creep;
                for (actual, expected) in history.inverse_plastic_cauchy_green.iter().zip([
                    (-2.0 * inelastic).exp(),
                    inelastic.exp(),
                    inelastic.exp(),
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
                        L * logarithm.exp_m1(),
                        3e-11,
                    );
                    close(
                        actual.element_first_piola_pa()[0][4 * axis],
                        expected.piola[axis],
                        0.5,
                    );
                }
                close(actual.free_energy_j(), V * expected.free, 3e-9);
                let energy = &report.material_energy;
                close(
                    energy.plastic_dissipation_j,
                    V * expected.plastic_loss,
                    3e-9,
                );
                close(energy.creep_dissipation_j, V * expected.creep_loss, 3e-9);
                close(energy.inelastic_dissipation_j(), V * expected.heat(), 3e-9);
                close(energy.caloric_heat_j, V * expected.heat(), 3e-9);
                close(energy.reversible_heat_j, 0.0, 0.0);
                close(energy.incremental_potential_j, V * expected.potential, 3e-9);
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
                let defect = work - V * (expected.free - old.free + expected.heat());
                close(energy.internal_endpoint_work_j, work, 3e-9);
                close(energy.endpoint_work_defect_j, defect, 3e-9);
                assert!(defect.abs() > 1e-8);
                close(
                    work,
                    energy.internal_energy_change_j
                        + report.thermal_audit.stored_energy_change_j
                        + defect,
                    4e-9,
                );
                assert!(expected.potential < expected.free + expected.heat());
                assert!(history.creep_alpha > saved_history[0].creep_alpha);
                if index < 2 {
                    assert!(history.alpha > saved_history[0].alpha);
                    assert!(energy.plastic_dissipation_j > 0.0);
                    close(expected.q, YIELD + HARDENING * expected.plastic, 1e-6);
                } else {
                    close(history.alpha, saved_history[0].alpha, 0.0);
                    close(energy.plastic_dissipation_j, 0.0, 0.0);
                    assert!(expected.q < YIELD + HARDENING * expected.plastic);
                }
                let material = &actual.families()[0];
                assert_eq!(
                    material.family().law().maximum_accumulated_creep_strain(),
                    Some(2.0)
                );
                assert_eq!(material.family().model_card(), &fixture.card.models()[0]);
                assert_eq!(material.family().model_card().state_schema_version, 2);
                assert!(std::ptr::eq(material.caloric_source(), &fixture.caloric));
                close(material.reference_density_kg_m3(), RHO, 0.0);
                if index == 0 {
                    let cold = mechanical(old, old.temperature, dt, exponent, load);
                    assert!(expected.creep > cold.creep + 1e-6);
                    let subdivided = solid
                        .advance(
                            &thermal,
                            &(),
                            end,
                            &prescribed,
                            &nodal_load,
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
                    let other = &subdivided.solid().mechanics().element_states()[0];
                    close(other.alpha, history.alpha, 3e-10);
                    close(other.creep_alpha, history.creep_alpha, 3e-10);
                    close(
                        subdivided.thermal().temperature_k(CELL).unwrap(),
                        expected.temperature,
                        2e-7,
                    );
                }
                if let Some(piola) = load {
                    assert!(report.inertia.is_none());
                    close(
                        energy.dead_load_endpoint_work_j,
                        V * piola * (expected.a.exp() - old.a.exp()),
                        3e-9,
                    );
                } else {
                    assert!(expected.q < old.q);
                    close(work, 0.0, 0.0);
                    let inertia = report.inertia.as_ref().unwrap();
                    close(inertia.numerical_loss_j, 0.5 * RHO * V, 1e-11);
                    close(inertia.kinetic_change_j, -inertia.numerical_loss_j, 1e-11);
                    assert!(inertia.numerical_loss_j > 1e-5);
                }
                assert_eq!(solid.mechanics().element_states(), saved_history.as_slice());
                let (next_thermal, next_solid, (), _) = candidate.into_parts();
                thermal = next_thermal;
                solid = next_solid;
                old = expected;
                time = end;
            }
            assert!(old.plastic > 0.05 && old.creep > old.plastic);
            assert!(solid.mechanics().displacement_m()[1][0] > 0.15 * L);
        }
    });
}
