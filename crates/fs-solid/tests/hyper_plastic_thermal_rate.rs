//! Physical-duration creep/relaxation through the original finite tet owner.
//! The diagonal linear-Perzyna oracle is independent of the material return.

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
const K: f64 = 16000.0;
const Y: f64 = 200.0;
const H: f64 = 800.0;
const ETA: f64 = 4000.0;
const ACTIVATION: f64 = 3000.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 174,
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
                bulk_modulus_pa: K,
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
    .with_overstress(FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: ETA,
        stress_scale_pa: 100.0,
        exponent: 1.0,
        activation_temperature_k: ACTIVATION,
    })
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
fn prescribe(a: f64, free_x: bool) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|dof| !(free_x && *dof == 3))
        .map(|dof| (dof, if dof == 3 { a.exp() - 1.0 } else { 0.0 }))
        .collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}; tolerance {tolerance:e}"
    );
}

#[derive(Clone, Copy, Default)]
struct Point {
    a: f64,
    p: f64,
    alpha: f64,
    free: f64,
    piola: f64,
    q: f64,
}
struct Return {
    point: Point,
    increment: f64,
    loss: f64,
    potential: f64,
}
fn oracle(old: Point, a: f64, temperature: f64, dt: f64) -> Return {
    let eta = ETA * (ACTIVATION * (1.0 / temperature - 1.0 / 300.0)).exp();
    let trial = 2.0 * MU * a - 3.0 * MU * old.p;
    let gamma = (trial.abs() - Y - H * old.alpha).max(0.0) / (3.0 * MU + H + eta / dt);
    let p = old.p + trial.signum() * gamma;
    let alpha = old.alpha + gamma;
    let elastic_deviator = 2.0 * a / 3.0 - p;
    let free =
        0.5 * K * a * a + 1.5 * MU * elastic_deviator * elastic_deviator + 0.5 * H * alpha * alpha;
    let piola = (K * a + 2.0 * MU * elastic_deviator) / a.exp();
    let viscous = eta * gamma * gamma / dt;
    Return {
        point: Point {
            a,
            p,
            alpha,
            free,
            piola,
            q: 3.0 * MU * elastic_deviator.abs(),
        },
        increment: gamma,
        loss: Y * gamma + viscous,
        potential: free + Y * gamma + 0.5 * viscous,
    }
}
fn loaded_oracle(old: Point, piola: f64, temperature: f64, dt: f64) -> Return {
    let (mut lo, mut hi) = (0.0, 0.5);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if oracle(old, mid, temperature, dt).point.piola < piola {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    oracle(old, 0.5 * (lo + hi), temperature, dt)
}
fn check(step: &FiniteThermoplasticTetStep<'_>, old: Point, expected: &Return) {
    let p = expected.point;
    let solid = &step.solid;
    close(solid.displacement_m()[1][0], p.a.exp() - 1.0, 3e-10);
    close(solid.element_states()[0].alpha, p.alpha, 2e-10);
    for (actual, expected) in solid.element_states()[0]
        .inverse_plastic_cauchy_green
        .iter()
        .zip([(-2.0 * p.p).exp(), p.p.exp(), p.p.exp(), 0.0, 0.0, 0.0])
    {
        close(*actual, expected, 3e-10);
    }
    close(solid.element_first_piola_pa()[0][0], p.piola, 2e-6);
    close(solid.free_energy_j(), V * p.free, 2e-6);
    close(step.energy.plastic_dissipation_j, V * expected.loss, 2e-6);
    close(step.energy.caloric_heat_j, V * expected.loss, 2e-6);
    close(step.energy.reversible_heat_j, 0.0, 0.0);
    close(
        step.energy.incremental_potential_j,
        V * expected.potential,
        2e-6,
    );
    let work = V * p.piola * (p.a.exp() - old.a.exp());
    close(step.energy.internal_endpoint_work_j, work, 2e-6);
    close(
        step.energy.endpoint_work_defect_j,
        work - V * (p.free - old.free + expected.loss),
        3e-6,
    );
    close(step.energy.balance_error_j, 0.0, 1e-6);
}

#[test]
fn held_finite_shape_relaxes_with_actual_duration_and_arrhenius_temperature_without_heating_potential()
 {
    with_cx(|cx| {
        let solid = initial(cx);
        let first = oracle(Point::default(), 0.16, 300.0, 0.04);
        let loaded = solid
            .advance(
                &[300.0],
                &prescribe(0.16, false),
                &[],
                [0.0; 3],
                0.04,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        check(&loaded, Point::default(), &first);
        assert!(first.point.q > Y + H * first.point.alpha + 100.0);
        assert!(
            loaded.energy.incremental_potential_j
                < loaded.solid.free_energy_j() + loaded.energy.plastic_dissipation_j - 1.0
        );
        let solid = loaded.solid;
        let accepted = solid.element_states().to_vec();
        let mut increments = Vec::new();
        for (dt, t) in [(0.02, 300.0), (0.4, 300.0), (0.02, 420.0)] {
            let expected = oracle(first.point, 0.16, t, dt);
            let candidate = solid
                .advance(
                    &[t],
                    &prescribe(0.16, false),
                    &[],
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check(&candidate, first.point, &expected);
            assert!(expected.point.piola < first.point.piola);
            assert!(expected.increment > 0.0);
            assert!(candidate.energy.plastic_dissipation_j > 0.0);
            assert_eq!(solid.element_states(), accepted.as_slice());
            increments.push(expected.increment);
        }
        assert!(increments[1] > increments[0]);
        assert!(increments[2] > increments[0]);
        // Arrest a moving, already plastic body while its stress relaxes. The
        // original kinetic loss cancels the kinetic change, and is never heat.
        let moving = solid
            .with_inertia(&[6.0], &[[0.01, 0.0, 0.0]; 4], cx)
            .unwrap();
        let expected = oracle(first.point, 0.16, 300.0, 0.02);
        let stopped = moving
            .advance(
                &[300.0],
                &prescribe(0.16, false),
                &[],
                [0.0; 3],
                0.02,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        check(&stopped, first.point, &expected);
        let inertia = stopped.inertia.unwrap();
        close(
            inertia.numerical_loss_j,
            0.5 * 6.0 * V * 0.01_f64.powi(2),
            1e-12,
        );
        close(inertia.kinetic_change_j, -inertia.numerical_loss_j, 1e-12);
        assert!(inertia.numerical_loss_j > 0.0);
    });
}

#[test]
fn held_load_creeps_through_real_free_coordinate_and_retries_from_the_same_overstressed_history() {
    with_cx(|cx| {
        let virgin = initial(cx);
        let first = oracle(Point::default(), 0.16, 300.0, 0.04);
        let load = [(3, V * first.point.piola)];
        let supports = prescribe(0.0, true);
        let initial = virgin
            .advance(
                &[300.0],
                &supports,
                &load,
                [0.0; 3],
                0.04,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        check(&initial, Point::default(), &first);
        let mut solid = initial.solid;
        let mut point = first.point;
        for (dt, t) in [(0.01, 300.0), (0.08, 300.0), (0.16, 360.0)] {
            let expected = loaded_oracle(point, first.point.piola, t, dt);
            assert!(expected.point.a > point.a);
            assert!(expected.point.alpha > point.alpha);
            let old_u = solid.displacement_m().to_vec();
            let old_state = solid.element_states().to_vec();
            let mut refused_controls = controls();
            refused_controls.material.maximum_endpoint_work_defect_j = 0.0;
            assert!(
                solid
                    .advance(
                        &[t],
                        &supports,
                        &load,
                        [0.0; 3],
                        dt,
                        settings(),
                        refused_controls,
                        cx
                    )
                    .is_err()
            );
            assert_eq!(solid.displacement_m(), old_u.as_slice());
            assert_eq!(solid.element_states(), old_state.as_slice());
            let candidate = solid
                .advance(
                    &[t],
                    &supports,
                    &load,
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check(&candidate, point, &expected);
            close(
                candidate.energy.dead_load_endpoint_work_j,
                load[0].1 * (expected.point.a.exp() - point.a.exp()),
                1e-6,
            );
            point = expected.point;
            solid = candidate.solid;
        }
        assert!(solid.displacement_m()[1][0] > 0.18);
    });
}
