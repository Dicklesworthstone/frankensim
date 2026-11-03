//! Independent momentum and finite-plastic endpoint oracles for native inertia.
//!
//! Synthetic fixed reference density supplies the original material inventory.
//! A free body exercises the actual dynamic Newton solve; a one-free-coordinate
//! finite stretch/reversal path checks complete history and each work term.

#[path = "hyper_plastic_dynamics/nodal_density.rs"]
mod nodal_density;
#[path = "hyper_plastic_dynamics/impulse.rs"]
mod impulse;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::plastic::dynamics::{FinitePlasticMotionControls, FinitePlasticTetMotion};
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetSolid,
};
use fs_solid::{HyperTetSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const VOLUME: f64 = 1.0 / 6.0;
const DENSITY: [f64; 1] = [6.0];
const MU: f64 = 8.0e3;
const KAPPA: f64 = 16.0e3;
const YIELD: f64 = 200.0;
const HARDENING: f64 = 800.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
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
        );
        f(&cx)
    })
}

fn law() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: KAPPA,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: HARDENING,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn initial_motion<'a>(
    material: &'a FiniteJ2,
    velocity: &[[f64; 3]],
    cx: &Cx<'_>,
) -> FinitePlasticTetMotion<'a> {
    let solid = FinitePlasticTetSolid::try_new(
        &NODES,
        &TETS,
        FinitePlasticLawField::Uniform(material),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap();
    FinitePlasticTetMotion::try_new(solid, &DENSITY, velocity, cx).unwrap()
}

fn controls() -> FinitePlasticMotionControls {
    FinitePlasticMotionControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1.0e-5,
            maximum_endpoint_work_defect_j: 200.0,
            yield_tolerance_pa: 1.0e-4,
        },
        maximum_numerical_kinetic_loss_j: 20.0,
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

fn norm_squared(vector: [f64; 3]) -> f64 {
    vector.iter().map(|value| value * value).sum()
}

#[test]
fn unrestrained_rigid_motion_and_body_force_keep_reference_mass_and_momentum() {
    with_cx(|cx| {
        let material = law();
        let velocity0 = [0.2, -0.1, 0.05];
        let dt = 0.02;
        let mass = DENSITY[0] * VOLUME;
        for acceleration in [[0.0; 3], [1.0, -0.5, 0.25]] {
            let body_force = acceleration.map(|a| DENSITY[0] * a);
            let mut motion = initial_motion(&material, &[velocity0; 4], cx);
            for interval in 1..=3 {
                let count = f64::from(interval);
                let old_v =
                    core::array::from_fn(|i| velocity0[i] + (count - 1.0) * dt * acceleration[i]);
                let expected_v =
                    core::array::from_fn(|i| velocity0[i] + count * dt * acceleration[i]);
                let expected_u: [f64; 3] = core::array::from_fn(|i| {
                    count * dt * velocity0[i]
                        + 0.5 * count * (count + 1.0) * dt * dt * acceleration[i]
                });
                let step = motion
                    .advance(&[], &[], body_force, dt, settings(), controls(), cx)
                    .unwrap();
                assert!(step.equilibrium.reactions_n.is_empty());
                for (u, v) in step
                    .motion
                    .solid()
                    .displacement_m()
                    .iter()
                    .zip(step.motion.velocity_mps())
                {
                    for component in 0..3 {
                        close(u[component], expected_u[component], 1.0e-10);
                        close(v[component], expected_v[component], 1.0e-8);
                    }
                }
                close(step.motion.reference_mass_kg(), mass, 1.0e-13);
                close(step.inertia.reference_mass_kg, mass, 1.0e-13);
                assert_eq!(step.motion.reference_density_kg_m3(), DENSITY);
                assert_eq!(
                    step.motion.solid().element_states(),
                    &[material.initial_state()]
                );
                close(step.material_energy.stored_after_j, 0.0, 1.0e-12);
                assert_eq!(step.material_energy.plastic_dissipation_j, 0.0);
                let before = 0.5 * mass * norm_squared(old_v);
                let after = 0.5 * mass * norm_squared(expected_v);
                let loss = 0.5 * mass * dt * dt * norm_squared(acceleration);
                let source_work =
                    mass * dt * (0..3).map(|i| acceleration[i] * expected_v[i]).sum::<f64>();
                close(step.inertia.kinetic_before_j, before, 1.0e-9);
                close(step.inertia.kinetic_after_j, after, 1.0e-9);
                close(step.inertia.kinetic_change_j, after - before, 1.0e-9);
                close(step.inertia.numerical_loss_j, loss, 1.0e-9);
                close(step.inertia.inertial_work_j, source_work, 1.0e-9);
                close(
                    step.material_energy.dead_load_endpoint_work_j,
                    source_work,
                    1.0e-9,
                );
                close(step.material_energy.support_endpoint_work_j, 0.0, 1.0e-12);
                close(step.material_energy.internal_endpoint_work_j, 0.0, 1.0e-12);
                close(step.balance_error_j, 0.0, 1.0e-9);
                motion = step.motion;
            }
        }
    });
}

#[test]
fn stress_free_plastic_shape_enters_unrestrained_motion_without_history_reset() {
    with_cx(|cx| {
        let material = law();
        let initial = FinitePlasticTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let prescribed = |a: f64| {
            let diagonal = [a.exp(), (-0.5 * a).exp(), (-0.5 * a).exp()];
            NODES
                .iter()
                .enumerate()
                .flat_map(|(node, x)| {
                    (0..3).map(move |i| (3 * node + i, (diagonal[i] - 1.0) * x[i]))
                })
                .collect::<Vec<_>>()
        };
        let a = 0.20;
        let plastic = (3.0 * MU * a - YIELD) / (3.0 * MU + HARDENING);
        let loaded = initial
            .advance(
                &prescribed(a),
                &[],
                [0.0; 3],
                settings(),
                controls().material,
                cx,
            )
            .unwrap();
        let unloaded = loaded
            .solid
            .advance(
                &prescribed(plastic),
                &[],
                [0.0; 3],
                settings(),
                controls().material,
                cx,
            )
            .unwrap();
        let old_states = unloaded.solid.element_states().to_vec();
        let old_f = unloaded.solid.element_deformation().to_vec();
        let old_u = unloaded.solid.displacement_m().to_vec();
        let storage = VOLUME * 0.5 * HARDENING * plastic * plastic;
        assert!(old_states[0].alpha > 0.18);
        close(old_states[0].alpha, plastic, 1.0e-12);
        close(unloaded.solid.stored_energy_j(), storage, 1.0e-9);
        assert!(
            unloaded.solid.element_first_piola_pa()[0]
                .iter()
                .all(|value| value.abs() < 1.0e-8)
        );

        let velocity = [0.2, -0.1, 0.05];
        let dt = 0.02;
        let motion =
            FinitePlasticTetMotion::try_new(unloaded.solid, &DENSITY, &[velocity; 4], cx).unwrap();
        assert_eq!(motion.solid().element_states(), old_states);
        let step = motion
            .advance(&[], &[], [0.0; 3], dt, settings(), controls(), cx)
            .unwrap();
        assert_eq!(step.motion.solid().element_states(), old_states);
        for (node, old_u) in old_u.iter().enumerate() {
            for component in 0..3 {
                close(
                    step.motion.solid().displacement_m()[node][component],
                    old_u[component] + dt * velocity[component],
                    1.0e-10,
                );
                close(
                    step.motion.velocity_mps()[node][component],
                    velocity[component],
                    1.0e-8,
                );
            }
        }
        for component in 0..9 {
            close(
                step.motion.solid().element_deformation()[0][component],
                old_f[0][component],
                1.0e-10,
            );
        }
        close(step.material_energy.stored_before_j, storage, 1.0e-9);
        close(step.material_energy.stored_after_j, storage, 1.0e-9);
        assert_eq!(step.material_energy.plastic_dissipation_j, 0.0);
        close(step.inertia.kinetic_change_j, 0.0, 1.0e-9);
        close(step.inertia.numerical_loss_j, 0.0, 1.0e-9);
        close(step.balance_error_j, 0.0, 1.0e-9);
    });
}

#[test]
fn finite_plastic_pulse_unload_reverse_and_kinetic_refusal_retain_full_motion() {
    with_cx(|cx| {
        let material = law();
        let mut motion = initial_motion(&material, &[[0.0; 3]; 4], cx);
        // Eleven stationary supports leave node 1's x coordinate to Newton.
        let prescribed: Vec<_> = (0..12)
            .filter(|dof| *dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let dt = 0.02;
        let total_mass = DENSITY[0] * VOLUME;
        let free_mass = total_mass / 10.0; // exact consistent P1 diagonal
        let mut old_plastic = 0.0;
        let mut old_alpha = 0.0;
        let mut old_u = 0.0;
        let mut old_v = 0.0;
        let mut old_storage = 0.0;
        for (interval, a) in [0.08_f64, 0.20, 0.19, -0.08, -0.16].into_iter().enumerate() {
            // Independent commuting finite-strain return: signed axial plastic
            // logarithmic stretch p differs from accumulated alpha on reversal.
            let trial = 2.0 * a / 3.0 - old_plastic;
            let q_trial = 3.0 * MU * trial.abs();
            let increment =
                (q_trial - YIELD - HARDENING * old_alpha).max(0.0) / (3.0 * MU + HARDENING);
            let plastic = old_plastic + trial.signum() * increment;
            let alpha = old_alpha + increment;
            let stretch = a.exp();
            let u = stretch - 1.0;
            let v = (u - old_u) / dt;
            let acceleration = (v - old_v) / dt;
            let returned = 2.0 * a / 3.0 - plastic;
            let tau_x = KAPPA * a + 2.0 * MU * returned;
            let piola_x = tau_x / stretch;
            let piola_lateral = KAPPA * a - MU * returned;
            let internal_force = VOLUME * piola_x;
            let applied_force = internal_force + free_mass * acceleration;
            let forces = [(3, applied_force)];

            if interval == 3 {
                // Reject after a genuine history has accumulated and the
                // reversed endpoint requires a nonzero kinetic numerical loss.
                let before_u = motion.solid().displacement_m().to_vec();
                let before_v = motion.velocity_mps().to_vec();
                let before_f = motion.solid().element_deformation().to_vec();
                let before_p = motion.solid().element_first_piola_pa().to_vec();
                let before_state = motion.solid().element_states().to_vec();
                let before_storage = motion.solid().stored_energy_j();
                let mut refused = controls();
                refused.maximum_numerical_kinetic_loss_j = 0.0;
                let error = motion
                    .advance(&prescribed, &forces, [0.0; 3], dt, settings(), refused, cx)
                    .unwrap_err();
                assert!(
                    error.to_string().contains("kinetic"),
                    "unexpected refusal: {error}"
                );
                assert_eq!(motion.solid().displacement_m(), before_u);
                assert_eq!(motion.velocity_mps(), before_v);
                assert_eq!(motion.solid().element_deformation(), before_f);
                assert_eq!(motion.solid().element_first_piola_pa(), before_p);
                assert_eq!(motion.solid().element_states(), before_state);
                assert_eq!(motion.solid().stored_energy_j(), before_storage);
                assert_eq!(motion.reference_density_kg_m3(), DENSITY);
                close(motion.reference_mass_kg(), total_mass, 1.0e-13);
            }

            let step = motion
                .advance(
                    &prescribed,
                    &forces,
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            close(step.motion.solid().displacement_m()[1][0], u, 1.0e-8);
            close(step.motion.velocity_mps()[1][0], v, 1.0e-6);
            for &(dof, value) in &prescribed {
                assert_eq!(
                    step.motion.solid().displacement_m()[dof / 3][dof % 3],
                    value
                );
                assert_eq!(step.motion.velocity_mps()[dof / 3][dof % 3], 0.0);
            }
            let expected_f = [stretch, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
            let expected_p = [
                piola_x,
                0.0,
                0.0,
                0.0,
                piola_lateral,
                0.0,
                0.0,
                0.0,
                piola_lateral,
            ];
            for component in 0..9 {
                close(
                    step.motion.solid().element_deformation()[0][component],
                    expected_f[component],
                    1.0e-8,
                );
                close(
                    step.motion.solid().element_first_piola_pa()[0][component],
                    expected_p[component],
                    1.0e-4,
                );
            }
            let state = &step.motion.solid().element_states()[0];
            close(state.alpha, alpha, 1.0e-8);
            for (actual, expected) in state.inverse_plastic_cauchy_green.iter().zip([
                (-2.0 * plastic).exp(),
                plastic.exp(),
                plastic.exp(),
                0.0,
                0.0,
                0.0,
            ]) {
                close(*actual, expected, 1.0e-7);
            }
            if interval == 2 {
                assert_eq!(increment, 0.0);
                assert_eq!(
                    step.motion.solid().element_states(),
                    motion.solid().element_states()
                );
            }
            if interval >= 3 {
                assert!(plastic < old_plastic);
                assert!(alpha > old_alpha);
            }
            close(step.motion.reference_mass_kg(), total_mass, 1.0e-13);
            assert_eq!(step.motion.reference_density_kg_m3(), DENSITY);

            let storage = VOLUME
                * (1.5 * MU * returned * returned
                    + 0.5 * KAPPA * a * a
                    + 0.5 * HARDENING * alpha * alpha);
            let dissipation = VOLUME * YIELD * increment;
            let material_work = internal_force * (u - old_u);
            let material_defect = material_work - (storage - old_storage) - dissipation;
            let kinetic_before = 0.5 * free_mass * old_v * old_v;
            let kinetic_after = 0.5 * free_mass * v * v;
            let kinetic_loss = 0.5 * free_mass * (v - old_v).powi(2);
            let inertial_work = free_mass * acceleration * (u - old_u);
            let external_work = applied_force * (u - old_u);
            close(step.material_energy.stored_before_j, old_storage, 1.0e-5);
            close(step.material_energy.stored_after_j, storage, 1.0e-5);
            close(
                step.material_energy.plastic_dissipation_j,
                dissipation,
                1.0e-5,
            );
            close(
                step.material_energy.incremental_potential_j,
                storage + dissipation,
                1.0e-5,
            );
            close(
                step.material_energy.internal_endpoint_work_j,
                material_work,
                1.0e-5,
            );
            close(
                step.material_energy.endpoint_work_defect_j,
                material_defect,
                1.0e-5,
            );
            close(
                step.material_energy.dead_load_endpoint_work_j,
                external_work,
                1.0e-5,
            );
            close(step.material_energy.support_endpoint_work_j, 0.0, 1.0e-12);
            close(step.inertia.kinetic_before_j, kinetic_before, 1.0e-5);
            close(step.inertia.kinetic_after_j, kinetic_after, 1.0e-5);
            close(
                step.inertia.kinetic_change_j,
                kinetic_after - kinetic_before,
                1.0e-5,
            );
            close(step.inertia.numerical_loss_j, kinetic_loss, 1.0e-5);
            close(step.inertia.inertial_work_j, inertial_work, 1.0e-5);
            close(step.balance_error_j, 0.0, 1.0e-5);
            close(
                external_work,
                storage - old_storage + dissipation + material_defect + kinetic_after
                    - kinetic_before
                    + kinetic_loss,
                1.0e-10,
            );

            // Consistent-mass off-diagonal entries transfer part of node 1's
            // impulse to the supports. The whole-body momentum uses the mass
            // column sum m/4, not the free diagonal m/10 or a lumped mass.
            for component in 0..3 {
                let reaction_sum: f64 = step
                    .equilibrium
                    .reactions_n
                    .iter()
                    .filter(|(dof, _)| dof % 3 == component)
                    .map(|(_, value)| value)
                    .sum();
                let force_sum = reaction_sum + if component == 0 { applied_force } else { 0.0 };
                let momentum_change = if component == 0 {
                    total_mass * (v - old_v) / 4.0
                } else {
                    0.0
                };
                close(force_sum * dt, momentum_change, 1.0e-5);
            }
            old_plastic = plastic;
            old_alpha = alpha;
            old_u = u;
            old_v = v;
            old_storage = storage;
            motion = step.motion;
        }
        assert!(old_alpha > 0.3);
        assert!(old_plastic < -0.08);
        assert!(motion.solid().element_deformation()[0][0] < 0.86);
    });
}
