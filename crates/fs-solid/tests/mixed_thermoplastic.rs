//! Independent thermal-volume, entropy, kinetic and coaxial inelastic checks.
//! The coefficients are manufactured; these bounded Tet10 fixtures do not
//! establish material calibration, general mesh stability or time accuracy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::mixed::plastic::thermal::{
    MixedThermoplasticControls, MixedThermoplasticTetError, MixedThermoplasticTetSolid,
    MixedThermoplasticTetStep,
};
use fs_solid::hyper3::mixed::{MixedUpError, MixedUpSettings};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

const VOLUME: f64 = 1.0 / 6.0;
const REFERENCE_T: f64 = 300.0;
const DENSITY: f64 = 6.0;
const DT: f64 = 0.2;
const MASS_BYTES: usize = 128 * 1024;

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 67,
                kernel_id: 25,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new(), f)
}

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

fn coefficients(temperature: f64) -> [f64; 4] {
    let d = temperature - REFERENCE_T;
    [
        80.0 - 0.04 * d,
        160.0 - 0.08 * d,
        2.0 - 0.004 * d,
        40.0 - 0.04 * d,
    ]
}

fn log_theta(temperature: f64) -> f64 {
    let d = temperature - REFERENCE_T;
    3.0e-4 * d + 0.5e-7 * d * d
}

fn expansion(temperature: f64) -> f64 {
    3.0e-4 + 1.0e-7 * (temperature - REFERENCE_T)
}

fn viscosity(temperature: f64) -> f64 {
    80.0 * (600.0 * (temperature.recip() - REFERENCE_T.recip())).exp()
}

fn family(combined: bool) -> AffineThermalFiniteJ2 {
    let native = AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: REFERENCE_T,
            temperature_range_k: [280.0, 600.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: 80.0,
                bulk_modulus_pa: 160.0,
                yield_stress_pa: 2.0,
                hardening_modulus_pa: 40.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: -0.04,
                bulk_modulus_pa_k: -0.08,
                yield_stress_pa_k: -0.004,
                hardening_modulus_pa_k: -0.04,
            },
            logarithmic_expansion_per_k: 3.0e-4,
            logarithmic_expansion_slope_per_k2: 1.0e-7,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap();
    if combined {
        native
            .with_norton_creep(
                FiniteJ2ThermalOverstress {
                    reference_viscosity_pa_s: 80.0,
                    stress_scale_pa: 4.0,
                    exponent: 1.0,
                    activation_temperature_k: 600.0,
                },
                2.0,
            )
            .unwrap()
    } else {
        native
    }
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 50,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2.0e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: 80.0,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

fn controls() -> MixedThermoplasticControls {
    MixedThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 2.0e-7,
            maximum_endpoint_work_defect_j: 10.0,
            yield_tolerance_pa: 1.0e-7,
        },
        maximum_numerical_kinetic_loss_j: 10.0,
    }
}

fn initial<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    treatment: VolumetricTreatment,
    combined: bool,
    cx: &Cx<'_>,
) -> MixedThermoplasticTetSolid<'a> {
    MixedThermoplasticTetSolid::try_new(
        nodes,
        elements,
        vec![family(combined)],
        &[REFERENCE_T],
        treatment,
        TetAssemblyBudget::standard(),
        4,
        64,
        cx,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.14e}, expected={expected:.14e}, tolerance={tolerance:e}"
    );
}

fn rigid_supports() -> [(usize, f64); 6] {
    [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)]
}

fn affine_supports(nodes: &[[f64; 3]], logarithms: [f64; 3]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |axis| (3 * node + axis, logarithms[axis].exp_m1() * x[axis]))
        })
        .collect()
}

fn balance(step: &MixedThermoplasticTetStep<'_>) {
    let e = &step.energy;
    let kinetic = step.inertia.as_ref().map_or(0.0, |inertia| {
        inertia.kinetic_change_j + inertia.numerical_loss_j
    });
    near(
        e.dead_load_endpoint_work_j + e.support_endpoint_work_j,
        e.internal_energy_change_j
            + e.caloric_heat_j
            + e.endpoint_work_defect_j
            + kinetic
            + e.balance_error_j,
        2.0e-7,
    );
    near(e.balance_error_j, 0.0, 2.0e-7);
    near(e.pressure_constraint_work_j, 0.0, 2.0e-7);
    near(
        e.caloric_heat_j,
        e.inelastic_dissipation_j() + e.reversible_heat_j,
        2.0e-12,
    );
    near(
        step.element_caloric_heat_j.iter().sum(),
        e.caloric_heat_j,
        2.0e-12,
    );
}

#[test]
fn free_heating_equilibrates_quadratic_nodes_to_thermal_volume_in_both_formulations() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let temperature = 400.0;
        let theta = log_theta(temperature).exp();
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            let solid = initial(&nodes, &elements, treatment, false, cx);
            let step = solid
                .advance(
                    &[temperature],
                    &rigid_supports(),
                    &[],
                    [0.0; 3],
                    None,
                    DT,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            for (x, u) in nodes.iter().zip(step.solid.displacement_m()) {
                for axis in 0..3 {
                    near(u[axis], (theta - 1.0) * x[axis], 2.0e-9);
                }
            }
            near(step.equilibrium.min_det_f, theta.powi(3), 2.0e-8);
            for pressure in step.solid.pressure_pa() {
                near(*pressure, 0.0, 2.0e-7);
            }
            for state in step.solid.quadrature_states() {
                assert_eq!(state, &family(false).initial_state());
            }
            near(step.solid.free_energy_j(), 0.0, 2.0e-10);
            near(step.solid.internal_energy_j(), 0.0, 2.0e-6);
            near(step.solid.entropy_j_k(), 0.0, 2.0e-9);
            near(step.energy.caloric_heat_j, 0.0, 2.0e-6);
            assert_eq!(step.energy.inelastic_dissipation_j(), 0.0);
            assert_eq!(step.solid.temperatures_k(), &[temperature]);
            assert_eq!(step.solid.families(), &[family(false)]);
            near(step.solid.reference_volume_m3(), VOLUME, 2.0e-14);
            assert!(step.solid.velocity_mps().is_none());
            assert!(step.solid.reference_mass_kg().is_none());
            balance(&step);
            assert_eq!(solid.temperatures_k(), &[REFERENCE_T]);
        }
    });
}

#[test]
fn clamped_heating_reports_physical_pressure_entropy_and_caloric_transfer() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let fixed = affine_supports(&nodes, [0.0; 3]);
        let temperature = 400.0;
        let solid = initial(
            &nodes,
            &elements,
            VolumetricTreatment::OriginalMaterial,
            false,
            cx,
        );
        let step = solid
            .advance(
                &[temperature],
                &fixed,
                &[],
                [0.0; 3],
                None,
                DT,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let bulk = coefficients(temperature)[1];
        let g = -3.0 * log_theta(temperature);
        // Total J=1, so pressure is K*ln(J/theta^3), in total-volume convention.
        let pressure = bulk * g;
        let free = 0.5 * bulk * g * g;
        let entropy = 3.0 * bulk * expansion(temperature) * g + 0.04 * g * g;
        // The mixed pressure owner uses the actual thermal-volume increment;
        // its signed endpoint-volume defect remains separate from caloric heat.
        let thermal_work = pressure * (3.0 * log_theta(temperature)).exp_m1();
        let reversible = thermal_work - temperature * entropy;
        for p in step.solid.pressure_pa() {
            near(*p, pressure, 2.0e-7);
        }
        for piola in step.solid.first_piola_pa() {
            for (entry, value) in piola.iter().enumerate() {
                near(*value, if entry % 4 == 0 { pressure } else { 0.0 }, 2.0e-7);
            }
        }
        near(step.solid.free_energy_j(), VOLUME * free, 2.0e-8);
        near(step.solid.entropy_j_k(), VOLUME * entropy, 2.0e-9);
        near(
            step.solid.internal_energy_j(),
            VOLUME * (free + temperature * entropy),
            2.0e-7,
        );
        near(step.energy.reversible_heat_j, VOLUME * reversible, 2.0e-7);
        near(step.energy.caloric_heat_j, VOLUME * reversible, 2.0e-7);
        near(step.energy.coefficient_storage_change_j, 0.0, 2.0e-12);
        near(
            step.energy.thermoelastic_storage_change_j,
            VOLUME * temperature * entropy,
            2.0e-7,
        );
        near(
            step.energy.endpoint_work_defect_j,
            VOLUME * (-thermal_work - free),
            2.0e-7,
        );
        assert_eq!(step.energy.inelastic_dissipation_j(), 0.0);
        near(step.energy.internal_endpoint_work_j, 0.0, 0.0);
        near(step.energy.support_endpoint_work_j, 0.0, 0.0);
        balance(&step);

        // Exact volume cannot accept a fully clamped change in thermal volume.
        let exact = initial(
            &nodes,
            &elements,
            VolumetricTreatment::ExactIncompressible,
            false,
            cx,
        );
        assert!(
            exact
                .advance(
                    &[temperature],
                    &fixed,
                    &[],
                    [0.0; 3],
                    None,
                    DT,
                    settings(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(exact.temperatures_k(), &[REFERENCE_T]);
        assert!(exact.displacement_m().iter().flatten().all(|u| *u == 0.0));
    });
}

#[derive(Clone, Copy, Default)]
struct Point {
    temperature: f64,
    axial: f64,
    lateral: f64,
    plastic: f64,
    creep: f64,
    e: f64,
    g: f64,
    pressure: f64,
    piola: [f64; 3],
    free: f64,
    entropy: f64,
    plastic_loss: f64,
    creep_loss: f64,
}

fn loaded(old: Point, temperature: f64, dead_piola: f64, exact: bool) -> Point {
    let [mu, bulk, yield_stress, hardening] = coefficients(temperature);
    let rate = DT / viscosity(temperature);
    let compliance = mu.recip() / 3.0 + if exact { 0.0 } else { bulk.recip() / 9.0 };
    // Free lateral faces give q=P*exp(a). Plastic complementarity and the
    // linear Norton endpoint then close a scalar equation without a return map.
    let residual = |q: f64| {
        (q / dead_piola).ln()
            - log_theta(temperature)
            - old.plastic.max((q - yield_stress) / hardening)
            - old.creep
            - rate * q
            - compliance * q
    };
    let (mut lo, mut hi) = (dead_piola, 4.0 * dead_piola);
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
    let plastic = old.plastic.max((q - yield_stress) / hardening);
    let creep = old.creep + rate * q;
    let axial = log_theta(temperature) + plastic + creep + compliance * q;
    let lateral = log_theta(temperature) - 0.5 * (plastic + creep) - q / (6.0 * mu)
        + if exact { 0.0 } else { q / (9.0 * bulk) };
    let e = q / (3.0 * mu);
    let g = if exact { 0.0 } else { q / (3.0 * bulk) };
    let free = 1.5 * mu * e * e + 0.5 * hardening * plastic * plastic + 0.5 * bulk * g * g;
    // Exact volume still has pressure entropy: 3*a*p*theta^3=a*q here.
    let entropy = expansion(temperature) * q
        + 0.04 * 1.5 * e * e
        + 0.04 * 0.5 * plastic * plastic
        + 0.08 * 0.5 * g * g;
    Point {
        temperature,
        axial,
        lateral,
        plastic,
        creep,
        e,
        g,
        pressure: q / (3.0 * (axial + 2.0 * lateral).exp()),
        piola: [dead_piola, 0.0, 0.0],
        free,
        entropy,
        plastic_loss: yield_stress * (plastic - old.plastic),
        creep_loss: q * (creep - old.creep),
    }
}

fn uniform_traction(elements: &[[usize; 10]], piola: f64) -> Vec<(usize, f64)> {
    // On the unit tetrahedron, integral(grad N_corner)=0 and
    // integral(grad N_ab)=V*(grad lambda_a+grad lambda_b).
    let gradient_x = [-1.0, 1.0, 0.0, 0.0];
    TET10_EDGES
        .iter()
        .enumerate()
        .filter_map(|(edge, &[a, b])| {
            let force = piola * VOLUME * (gradient_x[a] + gradient_x[b]);
            (force != 0.0).then_some((3 * elements[0][4 + edge], force))
        })
        .collect()
}

fn check_loaded(step: &MixedThermoplasticTetStep<'_>, nodes: &[[f64; 3]], old: Point, next: Point) {
    for (x, u) in nodes.iter().zip(step.solid.displacement_m()) {
        for (axis, strain) in [next.axial, next.lateral, next.lateral]
            .into_iter()
            .enumerate()
        {
            near(u[axis], strain.exp_m1() * x[axis], 2.0e-8);
        }
    }
    for state in step.solid.quadrature_states() {
        near(state.alpha, next.plastic, 2.0e-8);
        near(state.creep_alpha, next.creep, 2.0e-8);
        let total = next.plastic + next.creep;
        for (actual, expected) in state.inverse_plastic_cauchy_green.iter().zip([
            (-2.0 * total).exp(),
            total.exp(),
            total.exp(),
            0.0,
            0.0,
            0.0,
        ]) {
            near(*actual, expected, 5.0e-8);
        }
    }
    for p in step.solid.pressure_pa() {
        near(*p, next.pressure, 2.0e-6);
    }
    for piola in step.solid.first_piola_pa() {
        for (entry, value) in piola.iter().enumerate() {
            near(
                *value,
                if entry % 4 == 0 {
                    next.piola[entry / 4]
                } else {
                    0.0
                },
                2.0e-6,
            );
        }
    }
    let old_coefficients = coefficients(old.temperature);
    let new_coefficients = coefficients(next.temperature);
    let coefficient_change = 1.5 * (new_coefficients[0] - old_coefficients[0]) * old.e * old.e
        + 0.5 * (new_coefficients[1] - old_coefficients[1]) * old.g * old.g
        + 0.5 * (new_coefficients[3] - old_coefficients[3]) * old.plastic * old.plastic;
    let delta_ts = next.temperature * next.entropy - old.temperature * old.entropy;
    let hydrostatic = next.pressure * (next.axial + 2.0 * next.lateral).exp();
    let remainder_contraction = (next.piola[0] - hydrostatic * (-next.axial).exp())
        * old.axial.exp()
        + (next.piola[1] + next.piola[2] - 2.0 * hydrostatic * (-next.lateral).exp())
            * old.lateral.exp();
    let remainder_thermal_work =
        remainder_contraction * (log_theta(next.temperature) - log_theta(old.temperature)).exp_m1();
    let bulk_thermal_work = next.pressure
        * ((3.0 * log_theta(next.temperature)).exp() - (3.0 * log_theta(old.temperature)).exp())
        * old.g.exp();
    let thermal_work = remainder_thermal_work + bulk_thermal_work;
    let work = next.piola[0] * (next.axial.exp() - old.axial.exp())
        + (next.piola[1] + next.piola[2]) * (next.lateral.exp() - old.lateral.exp());
    let heat = next.plastic_loss + next.creep_loss;
    let reversible = thermal_work - coefficient_change - delta_ts;
    near(step.solid.free_energy_j(), VOLUME * next.free, 2.0e-7);
    near(step.solid.entropy_j_k(), VOLUME * next.entropy, 2.0e-9);
    near(
        step.solid.internal_energy_j(),
        VOLUME * (next.free + next.temperature * next.entropy),
        2.0e-7,
    );
    near(
        step.energy.plastic_dissipation_j,
        VOLUME * next.plastic_loss,
        2.0e-7,
    );
    near(
        step.energy.creep_dissipation_j,
        VOLUME * next.creep_loss,
        2.0e-7,
    );
    near(step.energy.reversible_heat_j, VOLUME * reversible, 2.0e-7);
    near(
        step.energy.caloric_heat_j,
        VOLUME * (heat + reversible),
        2.0e-7,
    );
    near(
        step.energy.coefficient_storage_change_j,
        VOLUME * coefficient_change,
        2.0e-7,
    );
    near(
        step.energy.thermoelastic_storage_change_j,
        VOLUME * delta_ts,
        2.0e-7,
    );
    near(step.energy.internal_endpoint_work_j, VOLUME * work, 2.0e-7);
    near(
        step.energy.endpoint_work_defect_j,
        VOLUME * (work - thermal_work - (next.free - old.free - coefficient_change) - heat),
        2.0e-7,
    );
    assert!(next.plastic_loss > 0.0 && next.creep_loss > 0.0);
    balance(step);
}

#[test]
fn heating_softens_and_activates_combined_flow_on_free_mixed_coordinates() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let loads = uniform_traction(&elements, 4.0);
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            let exact = treatment == VolumetricTreatment::ExactIncompressible;
            let mut solid = initial(&nodes, &elements, treatment, true, cx);
            let mut point = Point {
                temperature: REFERENCE_T,
                ..Point::default()
            };
            for temperature in [300.0, 350.0, 400.0] {
                let expected = loaded(point, temperature, 4.0, exact);
                let step = solid
                    .advance(
                        &[temperature],
                        &rigid_supports(),
                        &loads,
                        [0.0; 3],
                        None,
                        DT,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                check_loaded(&step, &nodes, point, expected);
                near(
                    step.equilibrium.min_det_f,
                    (expected.axial + 2.0 * expected.lateral).exp(),
                    2.0e-8,
                );
                assert_eq!(step.solid.temperatures_k(), &[temperature]);
                assert_eq!(step.solid.families(), &[family(true)]);
                point = expected;
                solid = step.solid;
            }
            assert!(point.plastic > 0.08 && point.creep > 0.04);
        }
    });
}

#[test]
fn connected_tetrahedra_retain_distinct_temperatures_and_element_heat_receipts() {
    with_cx(|cx| {
        let (nodes, elements) = promote_p1_to_p2(
            &[
                [0.0; 3],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
            ],
            &[[0, 1, 2, 3], [0, 2, 1, 4]],
            None,
            cx,
        )
        .unwrap();
        let temperatures = [310.0, 290.0];
        let logarithms = temperatures.map(log_theta);
        let stresses: [[f64; 3]; 2] = core::array::from_fn(|e| {
            let mu = coefficients(temperatures[e])[0];
            let l = logarithms[e];
            [
                4.0 * mu * l * (-3.0 * l).exp(),
                -2.0 * mu * l,
                -2.0 * mu * l,
            ]
        });
        let gradients = [
            [
                [-1.0, -1.0, -1.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ],
            [
                [1.0, -1.0, -1.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
            ],
        ];
        // The piecewise affine F=diag(theta^3,1,1) agrees at the shared x=0
        // face. Its pressure is zero. These explicit actuator loads include
        // the interface traction jump and use independent basis integrals.
        let mut forces = vec![0.0; 3 * nodes.len()];
        for (e, element) in elements.iter().enumerate() {
            for (edge, &[a, b]) in TET10_EDGES.iter().enumerate() {
                for axis in 0..3 {
                    forces[3 * element[4 + edge] + axis] += VOLUME
                        * stresses[e][axis]
                        * (gradients[e][a][axis] + gradients[e][b][axis]);
                }
            }
        }
        let loads: Vec<_> = forces
            .into_iter()
            .enumerate()
            .filter(|(_, f)| *f != 0.0)
            .collect();
        let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (3, 0.0), (5, 0.0), (6, 0.0)];
        let solid = MixedThermoplasticTetSolid::try_new(
            &nodes,
            &elements,
            vec![family(false), family(false)],
            &[REFERENCE_T; 2],
            VolumetricTreatment::ExactIncompressible,
            TetAssemblyBudget::standard(),
            4,
            128,
            cx,
        )
        .unwrap();
        let step = solid
            .advance(
                &temperatures,
                &supports,
                &loads,
                [0.0; 3],
                None,
                DT,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(step.solid.temperatures_k(), &temperatures);
        assert_eq!(step.solid.pressure_nodes().len(), 5);
        assert_eq!(step.element_caloric_heat_j.len(), 2);
        near(step.solid.reference_volume_m3(), 2.0 * VOLUME, 3.0e-14);
        for (x, u) in nodes.iter().zip(step.solid.displacement_m()) {
            let e = usize::from(x[0] < 0.0);
            near(u[0], (3.0 * logarithms[e]).exp_m1() * x[0], 2.0e-8);
            near(u[1], 0.0, 2.0e-8);
            near(u[2], 0.0, 2.0e-8);
        }
        for pressure in step.solid.pressure_pa() {
            near(*pressure, 0.0, 2.0e-6);
        }
        let points = step.solid.quadrature_points_per_element();
        let (mut free, mut entropy, mut internal) = (0.0, 0.0, 0.0);
        for e in 0..2 {
            let l = logarithms[e];
            let point_free = 6.0 * coefficients(temperatures[e])[0] * l * l;
            let point_entropy = 0.04 * 6.0 * l * l;
            let thermal_work = stresses[e].iter().sum::<f64>() * l.exp_m1();
            near(
                step.element_caloric_heat_j[e],
                VOLUME * (thermal_work - temperatures[e] * point_entropy),
                2.0e-7,
            );
            for q in e * points..(e + 1) * points {
                let f = &step.solid.deformation_gradients()[q];
                near(f[0], (3.0 * l).exp(), 2.0e-8);
                near(f[4], 1.0, 2.0e-8);
                near(f[8], 1.0, 2.0e-8);
                let p = &step.solid.first_piola_pa()[q];
                for axis in 0..3 {
                    near(p[4 * axis], stresses[e][axis], 2.0e-6);
                }
                assert_eq!(
                    step.solid.quadrature_states()[q],
                    family(false).initial_state()
                );
            }
            free += VOLUME * point_free;
            entropy += VOLUME * point_entropy;
            internal += VOLUME * (point_free + temperatures[e] * point_entropy);
        }
        near(step.solid.free_energy_j(), free, 2.0e-7);
        near(step.solid.entropy_j_k(), entropy, 2.0e-9);
        near(step.solid.internal_energy_j(), internal, 2.0e-7);
        assert_eq!(step.energy.inelastic_dissipation_j(), 0.0);
        balance(&step);
    });
}

#[test]
fn heated_motion_keeps_original_mass_and_rejects_complete_candidates_atomically() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let solid = initial(
            &nodes,
            &elements,
            VolumetricTreatment::OriginalMaterial,
            false,
            cx,
        )
        .with_inertia(&[DENSITY], &vec![[0.0; 3]; nodes.len()], MASS_BYTES, cx)
        .unwrap();
        let dt = 0.1;
        let temperature = 350.0;
        let heated = solid
            .advance(
                &[temperature],
                &affine_supports(&nodes, [log_theta(temperature); 3]),
                &[],
                [0.0; 3],
                None,
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let speed = log_theta(temperature).exp_m1() / dt;
        // integral |X|^2 dV=3V/10; the original quadratic mass integrates it.
        let kinetic = 0.5 * DENSITY * speed * speed * 3.0 * VOLUME / 10.0;
        near(
            heated.solid.reference_mass_kg().unwrap(),
            DENSITY * VOLUME,
            2.0e-14,
        );
        for (x, velocity) in nodes.iter().zip(heated.solid.velocity_mps().unwrap()) {
            for axis in 0..3 {
                near(velocity[axis], speed * x[axis], 2.0e-12);
            }
        }
        let inertia = heated.inertia.as_ref().unwrap();
        near(inertia.kinetic_after_j, kinetic, 2.0e-12);
        near(inertia.numerical_loss_j, kinetic, 2.0e-12);
        near(inertia.inertial_work_j, 2.0 * kinetic, 2.0e-12);
        near(
            heated.energy.inertial_endpoint_work_j,
            inertia.inertial_work_j,
            0.0,
        );
        near(heated.energy.support_endpoint_work_j, 2.0 * kinetic, 2.0e-7);
        near(heated.energy.caloric_heat_j, 0.0, 2.0e-7);
        assert_eq!(heated.energy.inelastic_dissipation_j(), 0.0);
        balance(&heated);

        let motion = heated.solid;
        let snapshot = || {
            (
                motion.displacement_m().to_vec(),
                motion.pressure_pa().to_vec(),
                motion.quadrature_states().to_vec(),
                motion.temperatures_k().to_vec(),
                motion.velocity_mps().unwrap().to_vec(),
                [
                    motion.free_energy_j(),
                    motion.internal_energy_j(),
                    motion.entropy_j_k(),
                    motion.reference_mass_kg().unwrap(),
                ],
                motion.deformation_gradients().to_vec(),
                motion.first_piola_pa().to_vec(),
            )
        };
        let before = snapshot();
        let next_temperature = 420.0;
        let moved = affine_supports(&nodes, [log_theta(next_temperature); 3]);
        assert!(
            motion
                .advance(
                    &[700.0],
                    &moved,
                    &[],
                    [0.0; 3],
                    None,
                    dt,
                    settings(),
                    controls(),
                    cx
                )
                .is_err()
        );
        assert!(
            motion
                .advance(
                    &[next_temperature],
                    &moved,
                    &[],
                    [0.0; 3],
                    None,
                    0.0,
                    settings(),
                    controls(),
                    cx
                )
                .is_err()
        );
        let strict = MixedThermoplasticControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            motion.advance(
                &[next_temperature],
                &moved,
                &[],
                [0.0; 3],
                None,
                dt,
                settings(),
                strict,
                cx
            ),
            Err(MixedThermoplasticTetError::EnergyGate {
                quantity: "backward-Euler numerical kinetic loss",
                ..
            })
        ));
        let cancelled = CancelGate::new();
        cancelled.request();
        with_gate(&cancelled, |cancelled_cx| {
            assert!(matches!(
                motion.advance(
                    &[next_temperature],
                    &moved,
                    &[],
                    [0.0; 3],
                    None,
                    dt,
                    settings(),
                    controls(),
                    cancelled_cx
                ),
                Err(MixedThermoplasticTetError::Native(MixedUpError::Cancelled))
            ));
        });
        assert_eq!(snapshot(), before);
        let retry = motion
            .advance(
                &[next_temperature],
                &moved,
                &[],
                [0.0; 3],
                None,
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let repeat = motion
            .advance(
                &[next_temperature],
                &moved,
                &[],
                [0.0; 3],
                None,
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.solid.displacement_m(), repeat.solid.displacement_m());
        assert_eq!(retry.solid.pressure_pa(), repeat.solid.pressure_pa());
        assert_eq!(
            retry.solid.quadrature_states(),
            repeat.solid.quadrature_states()
        );
        assert_eq!(retry.solid.temperatures_k(), repeat.solid.temperatures_k());
        assert_eq!(retry.solid.velocity_mps(), repeat.solid.velocity_mps());
        assert_eq!(retry.energy, repeat.energy);
        assert_eq!(retry.inertia, repeat.inertia);
        near(
            retry.solid.reference_mass_kg().unwrap(),
            DENSITY * VOLUME,
            2.0e-14,
        );
        assert!(retry.inertia.as_ref().unwrap().numerical_loss_j > 0.0);
        balance(&retry);
        assert_eq!(snapshot(), before);
    });
}
