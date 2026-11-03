//! Independent P2 momentum, shear-mode and mixed-history motion checks.
//! Backward Euler has measured numerical damping; its discrete acceleration
//! and convergent trajectories are tested separately from continuum kinematics.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::hyper3::mixed::MixedUpSettings;
use fs_solid::hyper3::mixed::visco::dynamics::{
    FiniteViscoMotionControls, MixedViscoMotionStep, MixedViscoTetMotion,
};
use fs_solid::hyper3::mixed::visco::{MixedViscoLawField, MixedViscoTetError, MixedViscoTetSolid};
use fs_solid::hyper3::visco::FiniteViscoStepControls;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

mod impulse;

const VOLUME: f64 = 1.0 / 6.0;
const DENSITY: f64 = 6.0;
const MU: f64 = 12.0;
const MASS_BYTES: usize = 128 * 1024;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let cancel = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &cancel,
            arena,
            StreamKey {
                seed: 61,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn mesh(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    promote_p1_to_p2(
        &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        &[[0, 1, 2, 3]],
        None,
        cx,
    )
    .unwrap()
}

fn law(bulk: f64, viscous: bool) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: MU / 2.0,
                c01: 0.0,
                kappa: bulk,
            },
            4.0,
        )
        .unwrap(),
        if viscous {
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3.0 * MU },
                tau_s: 0.07,
            }]
        } else {
            vec![]
        },
    )
    .unwrap()
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 30,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        length_scale_m: 1.0,
        stress_scale_pa: MU,
        maximum_assembly_bytes: 2 * 1024 * 1024,
        maximum_factor_bytes: 2 * 1024 * 1024,
    }
}

fn controls() -> FiniteViscoMotionControls {
    FiniteViscoMotionControls {
        material: FiniteViscoStepControls {
            equilibrium_work_tolerance_j: 2e-9,
            maximum_endpoint_work_defect_j: 1.0,
        },
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn solid<'a>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a ViscoHyperelastic,
    treatment: VolumetricTreatment,
    cx: &Cx<'_>,
) -> MixedViscoTetSolid<'a> {
    MixedViscoTetSolid::try_new(
        nodes,
        elements,
        MixedViscoLawField::Uniform(material),
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
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn work_identity(step: &MixedViscoMotionStep<'_>) {
    let energy = &step.material_energy;
    near(
        energy.inertial_endpoint_work_j,
        step.inertia.inertial_work_j,
        0.0,
    );
    near(
        energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j,
        energy.stored_after_j - energy.stored_before_j
            + energy.dissipated_j
            + energy.endpoint_work_defect_j
            + step.inertia.kinetic_change_j
            + step.inertia.numerical_loss_j
            + step.balance_error_j,
        2e-9,
    );
    near(step.inertia.balance_error_j, 0.0, 2e-11);
    near(step.balance_error_j, 0.0, 2e-9);
}

#[test]
fn free_exact_incompressible_body_has_correct_p2_acceleration_and_kinetic_balance() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e300, true);
        let velocity = [0.2, -0.1, 0.05];
        let acceleration = [0.3, -0.2, 0.1];
        let dt = 0.02;
        let mut motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::ExactIncompressible,
                cx,
            ),
            &[DENSITY],
            &vec![velocity; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        for count in 1..=3 {
            let n = f64::from(count);
            let expected_v: [f64; 3] =
                core::array::from_fn(|i| velocity[i] + n * dt * acceleration[i]);
            // BE uses the endpoint velocity, so acceleration travel is n(n+1)/2.
            let expected_u: [f64; 3] = core::array::from_fn(|i| {
                n * dt * velocity[i] + 0.5 * n * (n + 1.0) * dt * dt * acceleration[i]
            });
            let step = motion
                .advance(
                    dt,
                    &[],
                    &[],
                    acceleration.map(|a| DENSITY * a),
                    None,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert!(step.equilibrium.reactions_n.is_empty());
            assert_eq!(step.equilibrium.kkt_inertia.positive, 30);
            assert_eq!(step.equilibrium.kkt_inertia.negative, 4);
            for (u, v) in step
                .motion
                .solid()
                .displacement_m()
                .iter()
                .zip(step.motion.velocity_mps())
            {
                for i in 0..3 {
                    near(u[i], expected_u[i], 3e-11);
                    near(v[i], expected_v[i], 3e-9);
                }
            }
            for p in step.motion.solid().pressure_pa() {
                near(*p, 0.0, 2e-8);
            }
            near(step.motion.reference_mass_kg(), DENSITY * VOLUME, 2e-14);
            assert_eq!(step.motion.reference_density_kg_m3(), &[DENSITY]);
            near(step.material_energy.stored_after_j, 0.0, 2e-12);
            near(step.material_energy.dissipated_j, 0.0, 2e-14);
            near(
                step.inertia.kinetic_after_j,
                0.5 * DENSITY * VOLUME * expected_v.iter().map(|v| v * v).sum::<f64>(),
                2e-10,
            );
            near(
                step.inertia.numerical_loss_j,
                0.5 * DENSITY * VOLUME * dt * dt * acceleration.iter().map(|a| a * a).sum::<f64>(),
                2e-11,
            );
            work_identity(&step);
            motion = step.motion;
        }
    });
}

#[test]
fn quadratic_shear_oscillator_matches_scalar_be_and_first_order_time_refinement() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, false);
        let edge = TET10_EDGES
            .iter()
            .position(|e| e.contains(&2) && e.contains(&3))
            .unwrap();
        let node = elements[0][4 + edge];
        let dof = 3 * node;
        let amplitude = 0.02;
        let duration = 0.2;
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != dof)
            .map(|d| (d, 0.0))
            .collect();
        let loaded: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, if d == dof { amplitude } else { 0.0 }))
            .collect();
        // Exact finite shear u_x=4*q*y*z on the unit tet: J=1 and
        // W=8*mu*V*q²/5. The scalar mass follows the degree-four P2 integral.
        // Finite pressure compliance is required with this constrained space.
        let mass = 8.0 * DENSITY * VOLUME / 105.0;
        let stiffness = 16.0 * MU * VOLUME / 5.0;
        let omega = (stiffness / mass).sqrt();
        let mut errors = Vec::new();
        for intervals in [20, 40, 80] {
            let initial = solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            )
            .advance(
                0.0,
                &loaded,
                &[],
                [0.; 3],
                None,
                settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
            let mut motion = MixedViscoTetMotion::try_new(
                initial,
                &[DENSITY],
                &vec![[0.; 3]; nodes.len()],
                MASS_BYTES,
                cx,
            )
            .unwrap();
            let dt = duration / f64::from(intervals);
            let (mut q, mut v) = (amplitude, 0.0);
            let mut worst_oracle_error = 0.0_f64;
            for _ in 0..intervals {
                let old_q = q;
                q = (mass * q + mass * dt * v) / (mass + dt * dt * stiffness);
                v = (q - old_q) / dt;
                let step = motion
                    .advance(dt, &fixed, &[], [0.; 3], None, settings(), controls(), cx)
                    .unwrap();
                let actual = step.motion.solid().displacement_m()[node][0];
                worst_oracle_error = worst_oracle_error.max((actual - q).abs());
                near(actual, q, 8e-10);
                near(step.motion.velocity_mps()[node][0], v, 2e-8);
                near(
                    step.material_energy.stored_after_j,
                    0.5 * stiffness * q * q,
                    2e-11,
                );
                near(step.inertia.kinetic_after_j, 0.5 * mass * v * v, 2e-11);
                near(
                    step.material_energy.endpoint_work_defect_j,
                    0.5 * stiffness * (q - old_q).powi(2),
                    2e-11,
                );
                assert_eq!(step.material_energy.dissipated_j, 0.0);
                assert!(step.inertia.numerical_loss_j >= 0.0);
                work_identity(&step);
                motion = step.motion;
            }
            let actual_q = motion.solid().displacement_m()[node][0];
            let actual_v = motion.velocity_mps()[node][0];
            let phase = omega * duration;
            let error = (actual_q / amplitude - phase.cos())
                .hypot(actual_v / (amplitude * omega) + phase.sin());
            assert!(
                actual_q < 0.0,
                "a deforming quadratic mode must cross its reference"
            );
            errors.push(error);
            println!(
                "P2 shear intervals={intervals} trajectory_error={error:e} scalar_oracle_error={worst_oracle_error:e}"
            );
        }
        for pair in errors.windows(2) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                (0.9..1.1).contains(&order),
                "BE time order={order:e}, errors={errors:?}"
            );
        }
    });
}

#[test]
fn released_viscoelastic_motion_converges_to_exact_pressure_constraint() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let rollers: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(node, p)| {
                (0..3).filter_map(move |i| (p[i] == 0.0).then_some((3 * node + i, 0.0)))
            })
            .collect();
        let held: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(node, p)| {
                (0..3).filter_map(move |i| {
                    if i == 0 {
                        Some((3 * node, 0.04 * p[0]))
                    } else {
                        (p[i] == 0.0).then_some((3 * node + i, 0.0))
                    }
                })
            })
            .collect();
        let mut trajectories = Vec::new();
        for (treatment, bulk) in [
            (VolumetricTreatment::ExactIncompressible, 1e300),
            (VolumetricTreatment::OriginalMaterial, 1e7),
        ] {
            let material = law(bulk, true);
            let initial = solid(&nodes, &elements, &material, treatment, cx)
                .advance(
                    0.0,
                    &held,
                    &[],
                    [0.; 3],
                    None,
                    settings(),
                    controls().material,
                    cx,
                )
                .unwrap()
                .solid;
            let mut motion = MixedViscoTetMotion::try_new(
                initial,
                &[0.1 * DENSITY],
                &vec![[0.; 3]; nodes.len()],
                MASS_BYTES,
                cx,
            )
            .unwrap();
            let mut trajectory = Vec::new();
            let (mut physical_loss, mut numerical_loss, mut max_speed) = (0.0, 0.0, 0.0_f64);
            for _ in 0..24 {
                let step = motion
                    .advance(
                        0.01,
                        &rollers,
                        &[],
                        [0.; 3],
                        None,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                work_identity(&step);
                assert!(
                    step.material_energy.dissipated_j >= 0.0
                        && step.inertia.numerical_loss_j >= 0.0
                );
                physical_loss += step.material_energy.dissipated_j;
                numerical_loss += step.inertia.numerical_loss_j;
                max_speed = step
                    .motion
                    .velocity_mps()
                    .iter()
                    .flatten()
                    .fold(max_speed, |m, v| m.max(v.abs()));
                trajectory.push((
                    step.motion.solid().displacement_m().to_vec(),
                    step.motion.velocity_mps().to_vec(),
                ));
                motion = step.motion;
            }
            assert!(physical_loss > 1e-6 && numerical_loss > 1e-6 && max_speed > 0.01);
            assert!((motion.solid().displacement_m()[1][0] - 0.04).abs() > 0.01);
            near(motion.solid().time_s(), 0.24, 2e-15);
            println!(
                "released mixed {treatment:?}: physical_loss={physical_loss:e}, numerical_loss={numerical_loss:e}, max_speed={max_speed:e}"
            );
            trajectories.push(trajectory);
        }
        let (mut displacement_error, mut velocity_error) = (0.0_f64, 0.0_f64);
        for (a, b) in trajectories[0].iter().zip(&trajectories[1]) {
            for (x, y) in a.0.iter().flatten().zip(b.0.iter().flatten()) {
                displacement_error = displacement_error.max((x - y).abs());
            }
            for (x, y) in a.1.iter().flatten().zip(b.1.iter().flatten()) {
                velocity_error = velocity_error.max((x - y).abs());
            }
        }
        assert!(
            displacement_error < 2e-5 && velocity_error < 2e-4,
            "finite-bulk versus exact motion: u={displacement_error:e}, v={velocity_error:e}"
        );
        println!(
            "released bulk1e7 versus exact: displacement_error={displacement_error:e}, velocity_error={velocity_error:e}"
        );
    });
}

#[test]
fn moving_supports_keep_signed_p2_reactions_and_reject_motion_transactionally() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        // All displacement DOFs are prescribed; finite pressure compliance
        // avoids the additional exact-pressure null modes of this fixture.
        let motion = MixedViscoTetMotion::try_new(
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            ),
            &[DENSITY],
            &vec![[0.; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let dt = 0.02;
        let velocity = [0.3, -0.2, 0.1];
        let supports: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, dt * velocity[d % 3]))
            .collect();
        let before = (
            motion.solid().displacement_m().to_vec(),
            motion.solid().pressure_pa().to_vec(),
            motion.solid().quadrature_states().to_vec(),
            motion.velocity_mps().to_vec(),
            motion.solid().time_s(),
        );
        let strict = FiniteViscoMotionControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..controls()
        };
        let refusal = motion.advance(dt, &supports, &[], [0.; 3], None, settings(), strict, cx);
        assert!(
            matches!(
                refusal,
                Err(MixedViscoTetError::EnergyGate {
                    quantity: "backward-Euler numerical kinetic loss",
                    ..
                })
            ),
            "{refusal:?}"
        );
        assert_eq!(
            before,
            (
                motion.solid().displacement_m().to_vec(),
                motion.solid().pressure_pa().to_vec(),
                motion.solid().quadrature_states().to_vec(),
                motion.velocity_mps().to_vec(),
                motion.solid().time_s()
            )
        );
        let step = motion
            .advance(
                dt,
                &supports,
                &[],
                [0.; 3],
                None,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        for &(dof, reaction) in &step.equilibrium.reactions_n {
            let integral = if elements[0][..4].contains(&(dof / 3)) {
                -VOLUME / 20.0
            } else {
                VOLUME / 5.0
            };
            near(reaction, DENSITY * integral * velocity[dof % 3] / dt, 2e-10);
        }
        for actual in step.motion.velocity_mps() {
            for i in 0..3 {
                near(actual[i], velocity[i], 2e-12);
            }
        }
        let kinetic = 0.5 * DENSITY * VOLUME * velocity.iter().map(|v| v * v).sum::<f64>();
        near(step.inertia.kinetic_after_j, kinetic, 2e-12);
        near(step.inertia.numerical_loss_j, kinetic, 2e-12);
        near(
            step.material_energy.support_endpoint_work_j,
            2.0 * kinetic,
            2e-11,
        );
        near(step.material_energy.stored_after_j, 0.0, 2e-12);
        near(step.material_energy.dissipated_j, 0.0, 2e-14);
        work_identity(&step);
    });
}
