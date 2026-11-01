//! Finite temperature/history and momentum checks with independent coaxial laws.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetSolid,
};
use fs_solid::{HyperTetSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const V: f64 = 1.0 / 6.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
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
            temperature_range_k: [250.0, 650.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: 8000.0,
                bulk_modulus_pa: 16000.0,
                yield_stress_pa: 200.0,
                hardening_modulus_pa: 800.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: -6.0,
                bulk_modulus_pa_k: -8.0,
                yield_stress_pa_k: -0.4,
                hardening_modulus_pa_k: -1.0,
            },
            logarithmic_expansion_per_k: 3.0e-4,
            logarithmic_expansion_slope_per_k2: 1.0e-7,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.7,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}
fn controls() -> FiniteThermoplasticControls {
    FiniteThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1.0e-6,
            maximum_endpoint_work_defect_j: 800.0,
            yield_tolerance_pa: 1.0e-5,
        },
        maximum_numerical_kinetic_loss_j: 50.0,
    }
}
fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1.0e-11,
        ..HyperTetSettings::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:e}"
    );
}
fn prescribed_x(a: f64) -> Vec<(usize, f64)> {
    (0..12)
        .map(|dof| (dof, if dof == 3 { a.exp() - 1.0 } else { 0.0 }))
        .collect()
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

#[test]
fn free_thermal_expansion_uses_finite_stretch_and_native_free_coordinates() {
    with_cx(|cx| {
        let solid = initial(cx);
        // Remove rigid motions only: six coordinates actually equilibrate.
        let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let step = solid
            .advance(
                &[600.0],
                &supports,
                &[],
                [0.0; 3],
                0.1,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let logarithm: f64 = 3.0e-4 * 300.0 + 0.5 * 1.0e-7 * 300.0 * 300.0;
        let theta = logarithm.exp();
        for (node, x) in NODES.iter().enumerate() {
            for i in 0..3 {
                close(
                    step.solid.displacement_m()[node][i],
                    (theta - 1.0) * x[i],
                    1.0e-9,
                );
            }
        }
        assert!(step.equilibrium.steps[0].linear_iterations > 0);
        close(step.equilibrium.min_det_f, theta.powi(3), 1.0e-8);
        assert_eq!(step.solid.element_states(), &[family().initial_state()]);
        assert_eq!(step.energy.plastic_dissipation_j, 0.0);
        close(step.solid.free_energy_j(), 0.0, 1.0e-12);
        close(step.energy.caloric_heat_j, 0.0, 1.0e-5);
        close(step.energy.balance_error_j, 0.0, 1.0e-6);
        assert_eq!(step.solid.temperatures_k(), &[600.0]);
        assert_eq!(step.solid.families(), &[family()]);
        assert_eq!(step.solid.reference_volume_m3(0), Some(V));
        assert!(step.solid.velocity_mps().is_none());
        // Initial thermal state is not silently guessed from a nonreference T.
        assert!(
            FiniteThermoplasticTetSolid::try_new(
                &NODES,
                &TETS,
                vec![family()],
                &[600.0],
                TetAssemblyBudget::standard(),
                cx
            )
            .is_err()
        );
    });
}

#[derive(Clone, Copy)]
struct Point {
    a: f64,
    t: f64,
    plastic: f64,
    alpha: f64,
    mu: f64,
    bulk: f64,
    hardening: f64,
    log_theta: f64,
    g: f64,
    b: f64,
    free: f64,
    entropy: f64,
    internal: f64,
    piola: [f64; 3],
}
impl Point {
    // Independent diagonal multiplicative return, not a native response query.
    fn at(a: f64, t: f64, old: Option<Self>) -> (Self, f64) {
        let d = t - 300.0;
        let (mu, bulk, yield_stress, hardening) = (
            8000.0 - 6.0 * d,
            16000.0 - 8.0 * d,
            200.0 - 0.4 * d,
            800.0 - d,
        );
        let (mut plastic, mut alpha) = old.map_or((0.0, 0.0), |p| (p.plastic, p.alpha));
        let trial = 2.0 * a / 3.0 - plastic;
        let dg = ((3.0 * mu * trial.abs() - yield_stress - hardening * alpha)
            / (3.0 * mu + hardening))
            .max(0.0);
        plastic += trial.signum() * dg;
        alpha += dg;
        let log_theta = 3.0e-4 * d + 0.5e-7 * d * d;
        let g = a - 3.0 * log_theta;
        let b = 2.0 * a / 3.0 - plastic;
        let free = 1.5 * mu * b * b + 0.5 * bulk * g * g + 0.5 * hardening * alpha * alpha;
        let entropy = 3.0 * bulk * (3.0e-4 + 1.0e-7 * d) * g
            + 6.0 * 1.5 * b * b
            + 4.0 * g * g
            + 0.5 * alpha * alpha;
        let piola = [
            (bulk * g + 2.0 * mu * b) / a.exp(),
            bulk * g - mu * b,
            bulk * g - mu * b,
        ];
        (
            Self {
                a,
                t,
                plastic,
                alpha,
                mu,
                bulk,
                hardening,
                log_theta,
                g,
                b,
                free,
                entropy,
                internal: free + t * entropy,
                piola,
            },
            yield_stress * dg,
        )
    }
    fn account(self, old: Self, dissipation: f64) -> [f64; 5] {
        let coefficients = (self.mu - old.mu) * 1.5 * old.b * old.b
            + 0.5 * (self.bulk - old.bulk) * old.g * old.g
            + 0.5 * (self.hardening - old.hardening) * old.alpha * old.alpha;
        let delta_ts = self.t * self.entropy - old.t * old.entropy;
        let thermal_work = (self.piola[0] * old.a.exp() + self.piola[1] + self.piola[2])
            * (self.log_theta - old.log_theta).exp_m1();
        let work = self.piola[0] * (self.a.exp() - old.a.exp());
        let defect = work - thermal_work - (self.free - old.free - coefficients) - dissipation;
        let reversible = thermal_work - coefficients - delta_ts;
        [coefficients, delta_ts, work, defect, reversible]
    }
}

#[test]
fn softening_held_finite_shape_then_thermal_momentum_reverse_flow_is_atomic() {
    with_cx(|cx| {
        let mut solid = initial(cx);
        let (mut point, _) = Point::at(0.0, 300.0, None);
        // First create plastic history; heating at fixed F then softens its
        // former yield surface. The old accepted state remains valid at old T.
        for (a, t) in [(0.20, 300.0), (0.20, 450.0)] {
            let (next, dissipation) = Point::at(a, t, Some(point));
            let step = solid
                .advance(
                    &[t],
                    &prescribed_x(a),
                    &[],
                    [0.0; 3],
                    0.05,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert!(dissipation > 0.0);
            close(step.solid.element_states()[0].alpha, next.alpha, 1.0e-11);
            close(step.energy.plastic_dissipation_j, V * dissipation, 1.0e-8);
            let [coefficient, delta_ts, work, defect, reversible] =
                next.account(point, dissipation);
            close(step.solid.free_energy_j(), V * next.free, 1.0e-7);
            close(step.solid.internal_energy_j(), V * next.internal, 1.0e-7);
            close(
                step.energy.coefficient_storage_change_j,
                V * coefficient,
                1.0e-7,
            );
            close(
                step.energy.thermoelastic_storage_change_j,
                V * delta_ts,
                1.0e-7,
            );
            close(step.energy.support_endpoint_work_j, V * work, 1.0e-7);
            close(step.energy.endpoint_work_defect_j, V * defect, 1.0e-7);
            close(step.energy.reversible_heat_j, V * reversible, 1.0e-7);
            close(
                step.element_energy[0].caloric_heat_j_m3(),
                dissipation + reversible,
                1.0e-7,
            );
            if t > 300.0 {
                close(step.energy.internal_endpoint_work_j, 0.0, 1.0e-12);
                assert!(next.alpha > point.alpha);
                assert!(step.energy.reversible_heat_j.abs() > step.energy.plastic_dissipation_j);
            }
            solid = step.solid;
            point = next;
        }
        let before_attach = solid.element_states().to_vec();
        solid = solid.with_inertia(&[6.0], &[[0.0; 3]; 4], cx).unwrap();
        assert_eq!(solid.element_states(), before_attach);
        let supports: Vec<_> = (0..12)
            .filter(|dof| *dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let dt = 0.05;
        let mut velocity = 0.0;
        for (index, (a, t)) in [(0.25, 500.0), (0.24, 550.0), (-0.10, 450.0)]
            .into_iter()
            .enumerate()
        {
            let (next, dissipation) = Point::at(a, t, Some(point));
            let du = a.exp() - point.a.exp();
            let next_velocity = du / dt;
            // Reference mass=1, one free diagonal consistent-mass entry=1/10.
            let inertial_force = 0.1 * (next_velocity - velocity) / dt;
            let load = V * next.piola[0] + inertial_force;
            if index == 2 {
                let old_u = solid.displacement_m().to_vec();
                let old_v = solid.velocity_mps().unwrap().to_vec();
                let old_f = solid.element_deformation().to_vec();
                let old_p = solid.element_first_piola_pa().to_vec();
                let old_history = solid.element_states().to_vec();
                let old_t = solid.temperatures_k().to_vec();
                let mut refused_controls = controls();
                refused_controls.maximum_numerical_kinetic_loss_j = 0.0;
                let refusal = solid
                    .advance(
                        &[t],
                        &supports,
                        &[(3, load)],
                        [0.0; 3],
                        dt,
                        settings(),
                        refused_controls,
                        cx,
                    )
                    .unwrap_err();
                assert!(
                    refusal
                        .to_string()
                        .contains("backward-Euler numerical kinetic loss")
                );
                assert_eq!(solid.displacement_m(), old_u);
                assert_eq!(solid.velocity_mps().unwrap(), old_v);
                assert_eq!(solid.element_deformation(), old_f);
                assert_eq!(solid.element_first_piola_pa(), old_p);
                assert_eq!(solid.element_states(), old_history);
                assert_eq!(solid.temperatures_k(), old_t);
            }
            let step = solid
                .advance(
                    &[t],
                    &supports,
                    &[(3, load)],
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            let inertia = step.inertia.as_ref().unwrap();
            let [coefficient, delta_ts, work, defect, reversible] =
                next.account(point, dissipation);
            close(step.solid.displacement_m()[1][0], a.exp() - 1.0, 1.0e-9);
            close(
                step.solid.velocity_mps().unwrap()[1][0],
                next_velocity,
                1.0e-8,
            );
            close(step.solid.element_states()[0].alpha, next.alpha, 1.0e-10);
            for i in 0..3 {
                close(
                    step.solid.element_states()[0].inverse_plastic_cauchy_green[i],
                    (if i == 0 {
                        -2.0 * next.plastic
                    } else {
                        next.plastic
                    })
                    .exp(),
                    1.0e-10,
                );
                close(
                    step.solid.element_first_piola_pa()[0][4 * i],
                    next.piola[i],
                    1.0e-5,
                );
            }
            close(step.solid.free_energy_j(), V * next.free, 1.0e-6);
            close(step.solid.internal_energy_j(), V * next.internal, 1.0e-6);
            close(step.energy.plastic_dissipation_j, V * dissipation, 1.0e-6);
            close(
                step.energy.coefficient_storage_change_j,
                V * coefficient,
                1.0e-6,
            );
            close(
                step.energy.thermoelastic_storage_change_j,
                V * delta_ts,
                1.0e-6,
            );
            close(step.energy.internal_endpoint_work_j, V * work, 1.0e-6);
            close(step.energy.endpoint_work_defect_j, V * defect, 1.0e-6);
            close(step.energy.reversible_heat_j, V * reversible, 1.0e-6);
            close(
                step.energy.caloric_heat_j,
                V * (dissipation + reversible),
                1.0e-6,
            );
            close(
                inertia.kinetic_change_j,
                0.05 * (next_velocity * next_velocity - velocity * velocity),
                1.0e-7,
            );
            close(
                inertia.numerical_loss_j,
                0.05 * (next_velocity - velocity).powi(2),
                1.0e-7,
            );
            close(step.energy.dead_load_endpoint_work_j, load * du, 1.0e-6);
            close(step.energy.support_endpoint_work_j, 0.0, 1.0e-12);
            close(step.energy.balance_error_j, 0.0, 1.0e-6);
            assert_eq!(step.solid.reference_density_kg_m3(), Some([6.0].as_slice()));
            close(step.solid.reference_mass_kg().unwrap(), 1.0, 1.0e-13);
            assert_eq!(step.solid.families(), &[family()]);
            if index == 2 {
                assert!(next.plastic < point.plastic && next.alpha > point.alpha);
            }
            point = next;
            velocity = next_velocity;
            solid = step.solid;
        }
    });
}
