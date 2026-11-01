//! Finite Norton creep below a comparison plastic yield, with actual free
//! axial/transverse motion and an independent uniaxial Hencky reference.
//! These coefficients are synthetic, not a measured material calibration.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::FiniteJ2Limits;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalOverstress, FiniteNortonThermalParameters,
};
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
const ETA: f64 = 4000.0;
const ACTIVATION: f64 = 3000.0;
const DEAD_PIOLA: f64 = 100.0;
const COMPARISON_YIELD: f64 = 200.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 175,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn family() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new_norton(
        FiniteNortonThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [280.0, 600.0],
            shear_modulus_pa: MU,
            bulk_modulus_pa: BULK,
            shear_modulus_slope_pa_k: 0.0,
            bulk_modulus_slope_pa_k: 0.0,
            logarithmic_expansion_per_k: 0.0,
            logarithmic_expansion_slope_per_k2: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
        FiniteJ2ThermalOverstress {
            reference_viscosity_pa_s: ETA,
            stress_scale_pa: 100.0,
            exponent: 1.0,
            activation_temperature_k: ACTIVATION,
        },
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
        "{actual:.16e} != {expected:.16e}; tolerance {tolerance:e}"
    );
}
fn viscosity(temperature: f64) -> f64 {
    ETA * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp()
}

#[derive(Clone, Copy, Default)]
struct Point {
    a: f64,
    b: f64,
    creep: f64,
    q: f64,
    piola: [f64; 3],
    free: f64,
    dissipation: f64,
}
fn endpoint(old: Point, a: f64, b: f64, creep: f64) -> Point {
    let g = a + 2.0 * b;
    let e = 2.0 * (a - b) / 3.0 - creep;
    let q = 3.0 * MU * e;
    Point {
        a,
        b,
        creep,
        q,
        piola: [
            (BULK * g + 2.0 * MU * e) / a.exp(),
            (BULK * g - MU * e) / b.exp(),
            (BULK * g - MU * e) / b.exp(),
        ],
        free: 0.5 * BULK * g * g + 1.5 * MU * e * e,
        dissipation: q * (creep - old.creep),
    }
}
// Uniaxial Kirchhoff stress q gives g=q/(3K), e_x=q/(3mu).
// At held dead Piola P, q=P*exp(a) and a=p_old+dt*q/eta+q/(9K)+q/(3mu).
// Three actual free nodal coordinates recover these two lateral contractions.
fn loaded(old: Point, temperature: f64, dt: f64) -> Point {
    let compliance = 1.0 / (9.0 * BULK) + 1.0 / (3.0 * MU);
    let rate = dt / viscosity(temperature);
    let residual = |a: f64| a - old.creep - (compliance + rate) * DEAD_PIOLA * a.exp();
    let (mut lo, mut hi) = (old.creep, old.creep + 0.8);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let a = 0.5 * (lo + hi);
    let q = DEAD_PIOLA * a.exp();
    let creep = old.creep + rate * q;
    let b = q / (9.0 * BULK) - q / (6.0 * MU) - 0.5 * creep;
    endpoint(old, a, b, creep)
}
fn held(old: Point, temperature: f64, dt: f64) -> Point {
    let q = old.q / (1.0 + 3.0 * MU * dt / viscosity(temperature));
    let creep = old.creep + dt * q / viscosity(temperature);
    endpoint(old, old.a, old.b, creep)
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
    close(state.alpha, expected.creep, 2e-9);
    for (actual, expected) in state.inverse_plastic_cauchy_green.iter().zip([
        (-2.0 * expected.creep).exp(),
        expected.creep.exp(),
        expected.creep.exp(),
        0.0,
        0.0,
        0.0,
    ]) {
        close(*actual, expected, 5e-9);
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
    close(step.solid.free_energy_j(), V * expected.free, 1e-7);
    close(
        step.energy.plastic_dissipation_j,
        V * expected.dissipation,
        1e-7,
    );
    close(step.energy.caloric_heat_j, V * expected.dissipation, 1e-7);
    close(step.energy.reversible_heat_j, 0.0, 0.0);
    close(
        step.energy.incremental_potential_j,
        V * (expected.free + 0.5 * expected.dissipation),
        1e-7,
    );
    let work = V
        * (expected.piola[0] * (expected.a.exp() - old.a.exp())
            + (expected.piola[1] + expected.piola[2]) * (expected.b.exp() - old.b.exp()));
    close(step.energy.internal_endpoint_work_j, work, 1e-6);
    close(
        step.energy.endpoint_work_defect_j,
        work - V * (expected.free - old.free + expected.dissipation),
        1e-6,
    );
    close(step.energy.balance_error_j, 0.0, 1e-6);
}

#[test]
fn sustained_subyield_dead_load_drives_finite_creep_with_free_lateral_contraction() {
    with_cx(|cx| {
        let mut solid = initial(cx);
        let mut old = Point::default();
        let load = [(3, V * DEAD_PIOLA)];
        let supports = supports();
        for (dt, temperature) in [
            (1.0, 300.0),
            (1.0, 300.0),
            (0.5, 330.0),
            (0.5, 330.0),
            (0.5, 350.0),
            (0.5, 350.0),
        ] {
            let expected = loaded(old, temperature, dt);
            assert!(expected.q < COMPARISON_YIELD);
            assert!(expected.creep > old.creep && expected.a > old.a);
            let saved_history = solid.element_states().to_vec();
            let saved_u = solid.displacement_m().to_vec();
            let mut refuse = controls();
            refuse.material.maximum_endpoint_work_defect_j = 0.0;
            assert!(
                solid
                    .advance(
                        &[temperature],
                        &supports,
                        &load,
                        [0.0; 3],
                        dt,
                        settings(),
                        refuse,
                        cx
                    )
                    .is_err()
            );
            assert_eq!(solid.element_states(), saved_history.as_slice());
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
            close(
                step.energy.dead_load_endpoint_work_j,
                V * DEAD_PIOLA * (expected.a.exp() - old.a.exp()),
                1e-6,
            );
            assert!(step.energy.caloric_heat_j > 0.0);
            solid = step.solid;
            old = expected;
        }
        assert!(solid.displacement_m()[1][0] > 0.20);
        assert!(solid.displacement_m()[2][1] < -0.08);
        assert!(solid.families()[0].is_norton_creep());
    });
}

#[test]
fn held_creeping_body_relaxes_below_yield_to_zero_deviatoric_stress_and_arrests_without_extra_heat()
{
    with_cx(|cx| {
        let first = loaded(Point::default(), 300.0, 1.0);
        let step = initial(cx)
            .advance(
                &[300.0],
                &supports(),
                &[(3, V * DEAD_PIOLA)],
                [0.0; 3],
                1.0,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let mut solid = step
            .solid
            .with_inertia(&[6.0], &[[0.01, 0.0, 0.0]; 4], cx)
            .unwrap();
        let fixed = prescribed(first);
        let mut old = first;
        for index in 0..3 {
            let expected = held(old, 300.0, 1.0);
            let step = solid
                .advance(
                    &[300.0],
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
            close(step.energy.internal_endpoint_work_j, 0.0, 1e-9);
            assert!(expected.q < old.q && expected.dissipation > 0.0);
            assert!(step.energy.endpoint_work_defect_j > 0.0);
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
        assert!(old.q < 0.5);
        assert!(old.creep > first.creep);
        close(solid.displacement_m()[1][0], first.a.exp() - 1.0, 1e-12);
    });
}
