//! Simultaneous creep and yielding through actual free finite-solid coordinates.
//! The uniaxial Hencky/complementarity oracle is independent of the return map;
//! these synthetic coefficients are not a measured material calibration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetSolid, FiniteThermoplasticTetStep,
};
use fs_solid::{HyperTetSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const V: f64 = 1.0 / 6.0;
const MU: f64 = 8000.0;
const BULK: f64 = 16000.0;
const Y: f64 = 200.0;
const H: f64 = 4000.0;
const ETA: f64 = 8000.0;
const ACTIVATION: f64 = 3000.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 176,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn family() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [280.0, 600.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: MU,
                bulk_modulus_pa: BULK,
                yield_stress_pa: Y,
                hardening_modulus_pa: H,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: 0.0,
                bulk_modulus_pa_k: 0.0,
                yield_stress_pa_k: 0.0,
                hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: 0.0,
            logarithmic_expansion_slope_per_k2: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
    .with_norton_creep(
        FiniteJ2ThermalOverstress {
            reference_viscosity_pa_s: ETA,
            stress_scale_pa: 400.0,
            exponent: 1.0,
            activation_temperature_k: ACTIVATION,
        },
        2.0,
    )
    .unwrap()
}
fn initial(cx: &Cx<'_>) -> FiniteThermoplasticTetSolid<'static> {
    FiniteThermoplasticTetSolid::try_new(
        &NODES,
        &TETS,
        vec![family()],
        &[300.0],
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}
fn controls() -> FiniteThermoplasticControls {
    FiniteThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-6,
            maximum_endpoint_work_defect_j: 100.0,
            yield_tolerance_pa: 1e-6,
        },
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        ..HyperTetSettings::default()
    }
}
fn supports() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|dof| ![3, 7, 11].contains(dof))
        .map(|dof| (dof, 0.0))
        .collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:e}"
    );
}
fn viscosity(temperature: f64) -> f64 {
    ETA * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp()
}

#[derive(Clone, Copy, Default)]
struct Point {
    a: f64,
    b: f64,
    plastic: f64,
    creep: f64,
    q: f64,
    piola: [f64; 3],
    free: f64,
    plastic_loss: f64,
    creep_loss: f64,
}
fn point(old: Point, a: f64, b: f64, plastic: f64, creep: f64) -> Point {
    let g = a + 2.0 * b;
    let e = 2.0 * (a - b) / 3.0 - plastic - creep;
    let q = 3.0 * MU * e;
    Point {
        a,
        b,
        plastic,
        creep,
        q,
        piola: [
            (BULK * g + 2.0 * MU * e) / a.exp(),
            (BULK * g - MU * e) / b.exp(),
            (BULK * g - MU * e) / b.exp(),
        ],
        free: 0.5 * BULK * g * g + 1.5 * MU * e * e + 0.5 * H * plastic * plastic,
        plastic_loss: Y * (plastic - old.plastic),
        creep_loss: q * (creep - old.creep),
    }
}
// Equilibrium with free lateral faces gives q=P exp(a). Plastic consistency
// fixes alpha_p=max(alpha_p_old,(q-Y)/H); creep and elastic strains then determine
// a. This independently solves the endpoint complementarity, not two returns.
fn loaded(old: Point, temperature: f64, dt: f64, dead_piola: f64) -> Point {
    let compliance = 1.0 / (9.0 * BULK) + 1.0 / (3.0 * MU);
    let rate = dt / viscosity(temperature);
    let residual = |q: f64| {
        (q / dead_piola).ln() - old.plastic.max((q - Y) / H) - old.creep - rate * q - compliance * q
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
    let plastic = old.plastic.max((q - Y) / H);
    let creep = old.creep + rate * q;
    let total = plastic + creep;
    point(
        old,
        total + compliance * q,
        q / (9.0 * BULK) - q / (6.0 * MU) - 0.5 * total,
        plastic,
        creep,
    )
}
fn held(old: Point, temperature: f64, dt: f64) -> Point {
    let q = old.q / (1.0 + 3.0 * MU * dt / viscosity(temperature));
    point(
        old,
        old.a,
        old.b,
        old.plastic,
        old.creep + dt * q / viscosity(temperature),
    )
}
fn prescribed(point: Point) -> Vec<(usize, f64)> {
    (0..12)
        .map(|dof| {
            (
                dof,
                match dof {
                    3 => point.a.exp() - 1.0,
                    7 | 11 => point.b.exp() - 1.0,
                    _ => 0.0,
                },
            )
        })
        .collect()
}
fn check(step: &FiniteThermoplasticTetStep<'_>, old: Point, expected: Point) {
    let state = &step.solid.element_states()[0];
    close(state.alpha, expected.plastic, 2e-9);
    close(state.creep_alpha, expected.creep, 2e-9);
    let total = expected.plastic + expected.creep;
    for (actual, want) in state.inverse_plastic_cauchy_green.iter().zip([
        (-2.0 * total).exp(),
        total.exp(),
        total.exp(),
        0.0,
        0.0,
        0.0,
    ]) {
        close(*actual, want, 5e-9);
    }
    for (axis, logarithm) in [expected.a, expected.b, expected.b].into_iter().enumerate() {
        close(
            step.solid.displacement_m()[axis + 1][axis],
            logarithm.exp() - 1.0,
            2e-9,
        );
        close(
            step.solid.element_first_piola_pa()[0][4 * axis],
            expected.piola[axis],
            2e-6,
        );
    }
    let loss = expected.plastic_loss + expected.creep_loss;
    close(step.solid.free_energy_j(), V * expected.free, 1e-7);
    close(
        step.energy.plastic_dissipation_j,
        V * expected.plastic_loss,
        1e-7,
    );
    close(
        step.energy.creep_dissipation_j,
        V * expected.creep_loss,
        1e-7,
    );
    close(step.energy.inelastic_dissipation_j(), V * loss, 1e-7);
    close(step.energy.caloric_heat_j, V * loss, 1e-7);
    close(step.energy.reversible_heat_j, 0.0, 0.0);
    close(
        step.energy.incremental_potential_j,
        V * (expected.free + expected.plastic_loss + 0.5 * expected.creep_loss),
        1e-7,
    );
    let work = V
        * (expected.piola[0] * (expected.a.exp() - old.a.exp())
            + (expected.piola[1] + expected.piola[2]) * (expected.b.exp() - old.b.exp()));
    close(step.energy.internal_endpoint_work_j, work, 1e-6);
    close(
        step.energy.endpoint_work_defect_j,
        work - V * (expected.free - old.free + loss),
        1e-6,
    );
    close(step.energy.balance_error_j, 0.0, 1e-6);
    assert!(expected.q <= Y + H * expected.plastic + 1e-7);
}

#[test]
fn free_finite_body_yields_and_creeps_simultaneously_then_creeps_below_its_retained_yield() {
    with_cx(|cx| {
        let mut solid = initial(cx);
        let mut old = Point::default();
        let supports = supports();
        for (index, (dt, temperature, dead_piola)) in [
            (0.5, 300.0, 400.0),
            (0.5, 300.0, 400.0),
            (0.5, 330.0, 150.0),
            (0.5, 350.0, 150.0),
        ]
        .into_iter()
        .enumerate()
        {
            let expected = loaded(old, temperature, dt, dead_piola);
            let saved = solid.element_states().to_vec();
            let saved_u = solid.displacement_m().to_vec();
            let mut refused = controls();
            refused.material.maximum_endpoint_work_defect_j = 0.0;
            let load = [(3, V * dead_piola)];
            assert!(
                solid
                    .advance(
                        &[temperature],
                        &supports,
                        &load,
                        [0.0; 3],
                        dt,
                        settings(),
                        refused,
                        cx
                    )
                    .is_err()
            );
            assert_eq!(solid.element_states(), saved.as_slice());
            assert_eq!(solid.displacement_m(), saved_u.as_slice());
            let step = solid
                .advance(
                    &[temperature],
                    &supports,
                    &load,
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check(&step, old, expected);
            assert!(expected.creep > old.creep);
            if index < 2 {
                assert!(expected.plastic > old.plastic && expected.plastic_loss > 0.0);
                close(expected.q, Y + H * expected.plastic, 1e-7);
            } else {
                assert_eq!(expected.plastic, old.plastic);
                assert_eq!(step.energy.plastic_dissipation_j, 0.0);
                assert!(expected.q < Y + H * expected.plastic);
                assert!(expected.creep > loaded(old, 300.0, dt, dead_piola).creep);
            }
            close(
                step.energy.dead_load_endpoint_work_j,
                V * dead_piola * (expected.a.exp() - old.a.exp()),
                1e-6,
            );
            solid = step.solid;
            old = expected;
        }
        assert!(solid.displacement_m()[1][0] > 0.20);
        assert!(solid.displacement_m()[2][1] < -0.08);
        assert!(old.creep > old.plastic && old.plastic > 0.05);
        assert_eq!(
            solid.families()[0].maximum_accumulated_creep_strain(),
            Some(2.0)
        );
    });
}

#[test]
fn held_combined_history_relaxes_creep_without_hardening_or_heating_inertial_loss() {
    with_cx(|cx| {
        let first = loaded(Point::default(), 300.0, 0.5, 400.0);
        let loaded = initial(cx)
            .advance(
                &[300.0],
                &supports(),
                &[(3, V * 400.0)],
                [0.0; 3],
                0.5,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let mut solid = loaded
            .solid
            .with_inertia(&[6.0], &[[0.01, 0.0, 0.0]; 4], cx)
            .unwrap();
        let fixed = prescribed(first);
        let mut old = first;
        for index in 0..3 {
            let expected = held(old, 350.0, 1.0);
            let step = solid
                .advance(
                    &[350.0],
                    &fixed,
                    &[],
                    [0.0; 3],
                    1.0,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check(&step, old, expected);
            assert_eq!(
                step.solid.element_states()[0].alpha,
                solid.element_states()[0].alpha
            );
            assert_eq!(step.energy.plastic_dissipation_j, 0.0);
            assert!(step.energy.creep_dissipation_j > 0.0 && expected.q < old.q);
            close(step.energy.internal_endpoint_work_j, 0.0, 1e-9);
            if index == 0 {
                let inertia = step.inertia.as_ref().unwrap();
                close(
                    inertia.numerical_loss_j,
                    0.5 * 6.0 * V * 0.01_f64.powi(2),
                    1e-12,
                );
                close(inertia.kinetic_change_j, -inertia.numerical_loss_j, 1e-12);
            }
            solid = step.solid;
            old = expected;
        }
        assert!(old.q < 0.5 && old.creep > first.creep);
        close(old.plastic, first.plastic, 0.0);
        // Yield hardening survives the relaxation; creep carries none of it.
        close(
            old.free,
            0.5 * BULK * (old.a + 2.0 * old.b).powi(2)
                + old.q * old.q / (6.0 * MU)
                + 0.5 * H * first.plastic.powi(2),
            1e-8,
        );
    });
}
