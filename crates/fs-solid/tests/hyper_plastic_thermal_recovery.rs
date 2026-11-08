//! Recovering creep through the original finite thermal tetrahedral solve.
//! Independent Hencky/Arrhenius references exercise actual free coordinates,
//! thermal work, retained memory and original-reference momentum. Synthetic
//! coefficients verify implementation, not a calibrated alloy or service life.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::StrainHardeningCreep;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, AffineThermalRecoveringCreep, FiniteCreepThermalRecovery,
    FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress, FiniteJ2ThermalParameters,
    FiniteNortonThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetSolid, FiniteThermoplasticTetStep,
    SolidImpulseControls,
};
use fs_solid::{HyperTetSettings, TetAssemblyBudget};
use fs_solver::CheckedCgConfig;

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const V: f64 = 1.0 / 6.0;
const MU: f64 = 8000.0;
const BULK: f64 = 16000.0;
const MU_SLOPE: f64 = -4.0;
const BULK_SLOPE: f64 = -6.0;
const ETA: f64 = 4000.0;
const PRIMARY: f64 = 0.003;
const PIOLA: f64 = 600.0;
const BETA: f64 = 2e-5;
const BETA_SLOPE: f64 = 1e-8;

type Solid = FiniteThermoplasticTetSolid<'static, AffineThermalRecoveringCreep>;
type Step = FiniteThermoplasticTetStep<'static, AffineThermalRecoveringCreep>;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 190,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}
fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.25,
        maximum_elastic_log_strain: 0.8,
        maximum_plastic_log_strain: 2.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}
fn kinetics() -> FiniteJ2ThermalOverstress {
    FiniteJ2ThermalOverstress {
        reference_viscosity_pa_s: ETA,
        stress_scale_pa: PIOLA,
        exponent: 1.0,
        activation_temperature_k: 3000.0,
    }
}
fn recovering(base: AffineThermalFiniteJ2) -> AffineThermalRecoveringCreep {
    AffineThermalRecoveringCreep::new(
        base.with_norton_strain_hardening(StrainHardeningCreep {
            reference_strain: PRIMARY,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap(),
        FiniteCreepThermalRecovery {
            reference_recovery_rate_s: 2.0,
            activation_temperature_k: 2000.0,
        },
    )
    .unwrap()
}
fn family() -> AffineThermalRecoveringCreep {
    recovering(
        AffineThermalFiniteJ2::new_norton(
            FiniteNortonThermalParameters {
                reference_temperature_k: 300.0,
                temperature_range_k: [280.0, 600.0],
                shear_modulus_pa: MU,
                bulk_modulus_pa: BULK,
                shear_modulus_slope_pa_k: MU_SLOPE,
                bulk_modulus_slope_pa_k: BULK_SLOPE,
                logarithmic_expansion_per_k: BETA,
                logarithmic_expansion_slope_per_k2: BETA_SLOPE,
            },
            limits(),
            kinetics(),
        )
        .unwrap(),
    )
}
fn initial(family: AffineThermalRecoveringCreep, cx: &Cx<'_>) -> Solid {
    FiniteThermoplasticTetSolid::try_new(
        &NODES,
        &TETS,
        vec![family],
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
            maximum_endpoint_work_defect_j: 50.0,
            yield_tolerance_pa: 1e-6,
        },
        maximum_numerical_kinetic_loss_j: 1e-6,
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
fn moduli(t: f64) -> (f64, f64) {
    (MU + MU_SLOPE * (t - 300.0), BULK + BULK_SLOPE * (t - 300.0))
}
fn theta(t: f64) -> f64 {
    let d = t - 300.0;
    BETA * d + 0.5 * BETA_SLOPE * d * d
}
fn viscosity(t: f64) -> f64 {
    ETA * (3000.0 * (1.0 / t - 1.0 / 300.0)).exp()
}
fn recovery(t: f64) -> f64 {
    2.0 * (2000.0 * (1.0 / 300.0 - 1.0 / t)).exp()
}

#[derive(Clone, Copy)]
struct Point {
    t: f64,
    a: f64,
    b: f64,
    creep: f64,
    resistance: f64,
    q: f64,
    heat: f64,
    potential: f64,
}
impl Point {
    fn virgin() -> Self {
        Self {
            t: 300.0,
            a: 0.0,
            b: 0.0,
            creep: 0.0,
            resistance: 0.0,
            q: 0.0,
            heat: 0.0,
            potential: 0.0,
        }
    }
    fn volume_log(self) -> f64 {
        self.q / (3.0 * moduli(self.t).1)
    }
    fn deviator_squared(self) -> f64 {
        1.5 * (self.q / (3.0 * moduli(self.t).0)).powi(2)
    }
    fn free(self) -> f64 {
        let (mu, k) = moduli(self.t);
        mu * self.deviator_squared() + 0.5 * k * self.volume_log().powi(2)
    }
    fn entropy(self) -> f64 {
        3.0 * moduli(self.t).1 * (BETA + BETA_SLOPE * (self.t - 300.0)) * self.volume_log()
            - MU_SLOPE * self.deviator_squared()
            - 0.5 * BULK_SLOPE * self.volume_log().powi(2)
    }
}

// n=m=1 backward Euler gives dc*(a_eff+r_old+dc)=a_eff*dt*q/eta(T),
// with a_eff=a*(1+dt*k(T)). The separate resistance is (r_old+dc)/(1+dt*k).
fn increment(old_resistance: f64, q: f64, t: f64, dt: f64) -> f64 {
    let a = PRIMARY * (1.0 + dt * recovery(t));
    let shifted = a + old_resistance;
    let drive = a * dt * q / viscosity(t);
    2.0 * drive / (shifted + (shifted * shifted + 4.0 * drive).sqrt())
}
fn flow_potential(old_resistance: f64, dc: f64, t: f64, dt: f64) -> f64 {
    let a = PRIMARY * (1.0 + dt * recovery(t));
    viscosity(t) / dt * (0.5 * (1.0 + old_resistance / a) * dc * dc + dc.powi(3) / (3.0 * a))
}
fn endpoint(old: Point, t: f64, dt: f64, dead_piola: f64) -> Point {
    let (mu, k) = moduli(t);
    let compliance = 1.0 / (9.0 * k) + 1.0 / (3.0 * mu);
    let q = if dead_piola == 0.0 {
        0.0
    } else {
        let residual = |q: f64| {
            (q / dead_piola).ln()
                - theta(t)
                - old.creep
                - increment(old.resistance, q, t, dt)
                - compliance * q
        };
        let (mut lo, mut hi) = (dead_piola, 2.0 * dead_piola);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if residual(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    let dc = increment(old.resistance, q, t, dt);
    let creep = old.creep + dc;
    Point {
        t,
        a: theta(t) + creep + compliance * q,
        b: theta(t) - 0.5 * creep + q / (9.0 * k) - q / (6.0 * mu),
        creep,
        resistance: (old.resistance + dc) / (1.0 + dt * recovery(t)),
        q,
        heat: q * dc,
        potential: flow_potential(old.resistance, dc, t, dt),
    }
}
fn check(step: &Step, old: Point, next: Point) {
    let state = &step.solid.element_states()[0];
    close(state.native.alpha, next.creep, 2e-9);
    close(state.resistance_strain, next.resistance, 2e-9);
    close(state.native.creep_alpha, 0.0, 0.0);
    for (actual, expected) in state.native.inverse_plastic_cauchy_green.iter().zip([
        (-2.0 * next.creep).exp(),
        next.creep.exp(),
        next.creep.exp(),
        0.0,
        0.0,
        0.0,
    ]) {
        close(*actual, expected, 5e-9);
    }
    for (axis, log) in [next.a, next.b, next.b].into_iter().enumerate() {
        close(
            step.solid.displacement_m()[axis + 1][axis],
            log.exp() - 1.0,
            3e-9,
        );
    }
    let piola = next.q / next.a.exp();
    close(step.solid.element_first_piola_pa()[0][0], piola, 3e-6);
    for value in &step.solid.element_first_piola_pa()[0][1..] {
        close(*value, 0.0, 3e-6);
    }
    assert_eq!(step.solid.temperatures_k(), &[next.t]);
    close(step.solid.free_energy_j(), V * next.free(), 1e-7);
    close(
        step.solid.internal_energy_j(),
        V * (next.free() + next.t * next.entropy()),
        1e-7,
    );
    let delta_coeff = (next.t - old.t)
        * (MU_SLOPE * old.deviator_squared() + 0.5 * BULK_SLOPE * old.volume_log().powi(2));
    let delta_ts = next.t * next.entropy() - old.t * old.entropy();
    let thermal_work = piola * old.a.exp() * (theta(next.t) - theta(old.t)).exp_m1();
    let reversible = thermal_work - delta_coeff - delta_ts;
    let work = piola * (next.a.exp() - old.a.exp());
    close(
        step.energy.coefficient_storage_change_j,
        V * delta_coeff,
        1e-7,
    );
    close(
        step.energy.thermoelastic_storage_change_j,
        V * delta_ts,
        1e-7,
    );
    close(step.energy.inelastic_dissipation_j(), V * next.heat, 1e-7);
    close(step.energy.reversible_heat_j, V * reversible, 1e-7);
    close(
        step.energy.caloric_heat_j,
        V * (next.heat + reversible),
        1e-7,
    );
    close(
        step.energy.incremental_potential_j,
        V * (next.free() + next.potential),
        1e-7,
    );
    close(step.energy.internal_endpoint_work_j, V * work, 1e-7);
    close(
        step.energy.endpoint_work_defect_j,
        V * (work - thermal_work - next.free() + old.free() + delta_coeff - next.heat),
        1e-7,
    );
    close(step.energy.balance_error_j, 0.0, 1e-6);
    close(step.energy.equilibrium_work_residual_j, 0.0, 1e-6);
}

#[test]
fn heated_loaded_body_retains_permanent_creep_while_resistance_recovers() {
    with_cx(|_, cx| {
        let supports = supports();
        let loads = [(3, V * PIOLA)];
        let mut solid = initial(family(), cx);
        let mut point = Point::virgin();
        for (t, dt) in [(300.0, 0.2), (350.0, 0.15), (410.0, 0.08), (330.0, 0.2)] {
            let expected = endpoint(point, t, dt, PIOLA);
            let step = solid
                .advance(
                    &[t],
                    &supports,
                    &loads,
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check(&step, point, expected);
            assert!(expected.heat > expected.potential && expected.potential > 0.0);
            solid = step.solid;
            point = expected;
        }
        let before = solid.element_states().to_vec();
        let immediate = solid
            .advance(
                &[330.0],
                &supports,
                &loads,
                [0.0; 3],
                0.2,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        check(&immediate, point, endpoint(point, 330.0, 0.2, PIOLA));
        assert_eq!(solid.element_states(), before);
        for _ in 0..4 {
            let expected = endpoint(point, 410.0, 0.1, 0.0);
            let step = solid
                .advance(
                    &[410.0],
                    &supports,
                    &[],
                    [0.0; 3],
                    0.1,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            check(&step, point, expected);
            solid = step.solid;
            point = expected;
        }
        close(
            solid.element_states()[0].native.alpha,
            before[0].native.alpha,
            2e-9,
        );
        assert!(solid.element_states()[0].resistance_strain < before[0].resistance_strain);
        let reloaded = solid
            .advance(
                &[330.0],
                &supports,
                &loads,
                [0.0; 3],
                0.2,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        check(&reloaded, point, endpoint(point, 330.0, 0.2, PIOLA));
        assert!(
            reloaded.solid.element_states()[0].native.alpha
                > immediate.solid.element_states()[0].native.alpha
        );
        assert!(reloaded.solid.displacement_m()[1][0] > immediate.solid.displacement_m()[1][0]);
    });
}

#[test]
fn heated_recovering_history_survives_impulse_and_native_free_body_motion() {
    with_cx(|gate, cx| {
        let supports = supports();
        let loaded = initial(family(), cx)
            .advance(
                &[380.0],
                &supports,
                &[(3, V * PIOLA)],
                [0.0; 3],
                0.2,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let rested = loaded
            .solid
            .advance(
                &[410.0],
                &supports,
                &[],
                [0.0; 3],
                0.1,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let before = rested.solid.element_states()[0].clone();
        let displacement = rested.solid.displacement_m().to_vec();
        let moving = rested
            .solid
            .with_nodal_inertia(&[3.0, 9.0, 3.0, 9.0], &[[0.0; 3]; 4], cx)
            .unwrap();
        let velocity = [0.02, -0.01, 0.005];
        // Exact integrals of rho*N_i for rho=3+6X+6Z; total original mass is 1 kg.
        let weights = [9.0 / 40.0, 11.0 / 40.0, 9.0 / 40.0, 11.0 / 40.0];
        let impulses: Vec<_> = weights
            .into_iter()
            .enumerate()
            .flat_map(|(node, weight)| {
                velocity
                    .into_iter()
                    .enumerate()
                    .map(move |(axis, v)| (3 * node + axis, weight * v))
            })
            .collect();
        let impulse = moving
            .impulse_candidate(
                &[],
                &impulses,
                [0.0; 3],
                SolidImpulseControls {
                    solve: CheckedCgConfig {
                        tolerance: 1e-13,
                        max_iterations: 64,
                        max_corrections: 2,
                    },
                    maximum_workspace_components: 4096,
                    linear_tolerance_kg_m_s: 1e-11,
                    angular_tolerance_kg_m2_s: 1e-11,
                    work_tolerance_j: 1e-11,
                },
                cx,
            )
            .unwrap();
        assert_eq!(impulse.solid.element_states()[0], before);
        assert_eq!(impulse.solid.temperatures_k(), moving.temperatures_k());
        assert_eq!(impulse.solid.displacement_m(), moving.displacement_m());
        assert_eq!(
            impulse.solid.internal_energy_j(),
            moving.internal_energy_j()
        );
        let dt = 0.05;
        let step = impulse
            .solid
            .advance(&[410.0], &[], &[], [0.0; 3], dt, settings(), controls(), cx)
            .unwrap();
        for node in 0..4 {
            for axis in 0..3 {
                close(
                    step.solid.displacement_m()[node][axis],
                    displacement[node][axis] + dt * velocity[axis],
                    2e-8,
                );
                close(
                    step.solid.velocity_mps().unwrap()[node][axis],
                    velocity[axis],
                    1e-7,
                );
            }
        }
        close(
            step.solid.element_states()[0].native.alpha,
            before.native.alpha,
            2e-9,
        );
        close(
            step.solid.element_states()[0].resistance_strain,
            before.resistance_strain / (1.0 + dt * recovery(410.0)),
            2e-9,
        );
        close(step.solid.reference_mass_kg().unwrap(), 1.0, 1e-14);
        close(
            step.inertia.as_ref().unwrap().kinetic_after_j,
            0.5 * velocity.iter().map(|v| v * v).sum::<f64>(),
            1e-9,
        );
        close(step.inertia.as_ref().unwrap().numerical_loss_j, 0.0, 1e-9);
        close(step.energy.caloric_heat_j, 0.0, 1e-7);
        assert!(
            impulse
                .solid
                .advance(&[700.0], &[], &[], [0.0; 3], dt, settings(), controls(), cx)
                .is_err()
        );
        assert_eq!(impulse.solid.element_states()[0], before);
        gate.request();
        assert!(
            impulse
                .solid
                .advance(&[410.0], &[], &[], [0.0; 3], dt, settings(), controls(), cx)
                .is_err()
        );
        assert_eq!(impulse.solid.element_states()[0], before);
    });
}

#[test]
fn heated_combined_recovery_keeps_yield_and_creep_heat_separate() {
    with_cx(|_, cx| {
        let y = 200.0;
        let h = 4000.0;
        let base = AffineThermalFiniteJ2::new(
            FiniteJ2ThermalParameters {
                reference_temperature_k: 300.0,
                temperature_range_k: [280.0, 600.0],
                reference: FiniteJ2Parameters {
                    shear_modulus_pa: MU,
                    bulk_modulus_pa: BULK,
                    yield_stress_pa: y,
                    hardening_modulus_pa: h,
                },
                slopes: FiniteJ2TemperatureSlopes {
                    shear_modulus_pa_k: MU_SLOPE,
                    bulk_modulus_pa_k: BULK_SLOPE,
                    yield_stress_pa_k: 0.0,
                    hardening_modulus_pa_k: 0.0,
                },
                logarithmic_expansion_per_k: BETA,
                logarithmic_expansion_slope_per_k2: BETA_SLOPE,
            },
            limits(),
        )
        .unwrap()
        .with_norton_creep(kinetics(), 3.0)
        .unwrap();
        let t = 380.0;
        let dt = 0.1;
        let axial: f64 = 0.15;
        let (mu, k) = moduli(t);
        let compliance = 1.0 / (9.0 * k) + 1.0 / (3.0 * mu);
        let residual =
            |q: f64| theta(t) + compliance * q + (q - y) / h + increment(0.0, q, t, dt) - axial;
        let (mut lo, mut hi) = (y, 2000.0);
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
        let plastic = (q - y) / h;
        let creep = increment(0.0, q, t, dt);
        let lateral = theta(t) - 0.5 * (plastic + creep) + q / (9.0 * k) - q / (6.0 * mu);
        let mut supports = supports();
        supports.push((3, axial.exp() - 1.0));
        let step = initial(recovering(base), cx)
            .advance(
                &[t],
                &supports,
                &[],
                [0.0; 3],
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let state = &step.solid.element_states()[0];
        close(state.native.alpha, plastic, 2e-9);
        close(state.native.creep_alpha, creep, 2e-9);
        close(
            state.resistance_strain,
            creep / (1.0 + dt * recovery(t)),
            2e-9,
        );
        close(step.solid.displacement_m()[2][1], lateral.exp() - 1.0, 3e-9);
        close(step.solid.displacement_m()[3][2], lateral.exp() - 1.0, 3e-9);
        close(
            step.solid.element_first_piola_pa()[0][0],
            q / axial.exp(),
            3e-6,
        );
        close(step.energy.plastic_dissipation_j, V * y * plastic, 1e-7);
        close(step.energy.creep_dissipation_j, V * q * creep, 1e-7);
        let free = V * (0.5 * compliance * q * q + 0.5 * h * plastic * plastic);
        close(step.solid.free_energy_j(), free, 1e-7);
        close(
            step.energy.incremental_potential_j,
            free + V * (y * plastic + flow_potential(0.0, creep, t, dt)),
            1e-7,
        );
        close(step.energy.balance_error_j, 0.0, 1e-6);
    });
}
