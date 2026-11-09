//! Independent P2 mass, elastic-mode, inelastic-heat and weak-volume checks.
//! Synthetic coefficients test bounded backward-Euler motion, not calibrated
//! materials, general mesh stability or an energy-conserving time integrator.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::highorder::simplex::duffy_quadrature;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::plastic::finite::combined::FiniteNortonJ2;
use fs_material::plastic::finite::mixed::MixedFiniteInelasticMaterial;
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::mixed::plastic::dynamics::{
    FinitePlasticMotionControls, MixedPlasticMotionStep, MixedPlasticTetMotion,
};
use fs_solid::hyper3::mixed::plastic::{
    MixedPlasticLawField, MixedPlasticTetError, MixedPlasticTetSolid,
};
use fs_solid::hyper3::mixed::{MixedUpError, MixedUpSettings};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
use fs_solver::nonlinear::NewtonKrylovConfig;

const MU: f64 = 80.0;
const YIELD: f64 = 2.0;
const HARDENING: f64 = 8.0;
const VOLUME: f64 = 1.0 / 6.0;
const DENSITY: f64 = 6.0;
const DT: f64 = 0.02;
const ETA: f64 = 8.0;
const MASS_BYTES: usize = 128 * 1024;

fn with_gate(gate: &CancelGate, f: impl FnOnce(&Cx<'_>)) {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 65,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    with_gate(&CancelGate::new(), f);
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

fn law(bulk: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: bulk,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: HARDENING,
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

fn kinetics() -> FiniteOverstressParameters {
    FiniteOverstressParameters {
        viscosity_pa_s: ETA,
        stress_scale_pa: YIELD,
        exponent: 1.0,
    }
}

fn settings() -> MixedUpSettings {
    MixedUpSettings {
        load_steps: 1,
        maximum_newton_iterations: 40,
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

fn controls() -> FinitePlasticMotionControls {
    FinitePlasticMotionControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 2e-8,
            maximum_endpoint_work_defect_j: 10.0,
            yield_tolerance_pa: 2e-8,
        },
        maximum_numerical_kinetic_loss_j: 10.0,
    }
}

fn solid<'a, L: MixedFiniteInelasticMaterial>(
    nodes: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: &'a L,
    treatment: VolumetricTreatment,
    cx: &Cx<'_>,
) -> MixedPlasticTetSolid<'a, L> {
    MixedPlasticTetSolid::try_new(
        nodes,
        elements,
        MixedPlasticLawField::Uniform(material),
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

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn work_identity<L: MixedFiniteInelasticMaterial>(step: &MixedPlasticMotionStep<'_, L>) {
    let e = &step.material_energy;
    near(
        e.inertial_endpoint_work_j,
        step.inertia.inertial_work_j,
        0.0,
    );
    near(
        e.dead_load_endpoint_work_j + e.support_endpoint_work_j,
        e.stored_after_j - e.stored_before_j
            + e.inelastic_dissipation_j()
            + e.endpoint_work_defect_j
            + step.inertia.kinetic_change_j
            + step.inertia.numerical_loss_j
            + step.balance_error_j,
        2e-8,
    );
    near(step.inertia.balance_error_j, 0.0, 2e-11);
    near(step.balance_error_j, 0.0, 2e-8);
}

#[test]
fn free_body_preserves_rigid_motion_and_original_mass_momentum_in_both_volume_formulations() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let initial_velocity = [0.2, -0.1, 0.05];
        let origin = [-0.4, 0.15, -0.2];
        let mass = DENSITY * VOLUME;
        for (treatment, bulk) in [
            (VolumetricTreatment::OriginalMaterial, 160.0),
            (VolumetricTreatment::ExactIncompressible, 1e308),
        ] {
            let material = law(bulk);
            for acceleration in [[0.0; 3], [0.3, -0.2, 0.1]] {
                let mut motion = MixedPlasticTetMotion::try_new(
                    solid(&nodes, &elements, &material, treatment, cx),
                    &[DENSITY],
                    &vec![initial_velocity; nodes.len()],
                    MASS_BYTES,
                    cx,
                )
                .unwrap();
                let loads = motion
                    .reference_acceleration_load_n(&vec![acceleration; nodes.len()], cx)
                    .unwrap();
                for (dof, force) in loads {
                    let integral = if elements[0][..4].contains(&(dof / 3)) {
                        -VOLUME / 20.0
                    } else {
                        VOLUME / 5.0
                    };
                    near(force, DENSITY * integral * acceleration[dof % 3], 2e-13);
                }
                for interval in 1..=2 {
                    let n = f64::from(interval);
                    let velocity: [f64; 3] =
                        core::array::from_fn(|i| initial_velocity[i] + n * DT * acceleration[i]);
                    let displacement: [f64; 3] = core::array::from_fn(|i| {
                        n * DT * initial_velocity[i]
                            + 0.5 * n * (n + 1.0) * DT * DT * acceleration[i]
                    });
                    let step = motion
                        .advance(
                            DT,
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
                    near(step.motion.reference_mass_kg(), mass, 2e-14);
                    assert_eq!(step.motion.reference_density_kg_m3(), &[DENSITY]);
                    assert_eq!(step.motion.volumetric_treatment(), treatment);
                    for (u, v) in step
                        .motion
                        .solid()
                        .displacement_m()
                        .iter()
                        .zip(step.motion.velocity_mps())
                    {
                        for axis in 0..3 {
                            near(u[axis], displacement[axis], 3e-11);
                            near(v[axis], velocity[axis], 3e-9);
                        }
                    }
                    let momentum = step.motion.momentum(origin, cx).unwrap();
                    let angular = cross(
                        core::array::from_fn(|i| 0.25 + displacement[i] - origin[i]),
                        velocity.map(|v| mass * v),
                    );
                    for axis in 0..3 {
                        near(momentum.linear_kg_m_s[axis], mass * velocity[axis], 3e-9);
                        near(momentum.angular_kg_m2_s[axis], angular[axis], 3e-9);
                    }
                    let kinetic = 0.5 * mass * velocity.iter().map(|v| v * v).sum::<f64>();
                    near(step.inertia.kinetic_after_j, kinetic, 2e-10);
                    near(step.motion.kinetic_energy_j(cx).unwrap(), kinetic, 2e-10);
                    near(
                        step.inertia.numerical_loss_j,
                        0.5 * mass * DT * DT * acceleration.iter().map(|a| a * a).sum::<f64>(),
                        2e-11,
                    );
                    near(step.material_energy.stored_after_j, 0.0, 2e-12);
                    assert_eq!(step.material_energy.inelastic_dissipation_j(), 0.0);
                    for state in step.motion.solid().quadrature_states() {
                        assert_eq!(state, &material.initial_state());
                    }
                    for pressure in step.motion.solid().pressure_pa() {
                        near(*pressure, 0.0, 2e-8);
                    }
                    work_identity(&step);
                    motion = step.motion;
                }
            }
        }
    });
}

#[test]
fn elastic_quadratic_shear_mode_matches_independent_scalar_backward_euler() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(160.0);
        let edge = TET10_EDGES
            .iter()
            .position(|e| e.contains(&2) && e.contains(&3))
            .unwrap();
        let node = elements[0][4 + edge];
        let dof = 3 * node;
        let amplitude = 1e-5;
        let loaded: Vec<_> = (0..3 * nodes.len())
            .map(|d| (d, if d == dof { amplitude } else { 0.0 }))
            .collect();
        let fixed: Vec<_> = (0..3 * nodes.len())
            .filter(|&d| d != dof)
            .map(|d| (d, 0.0))
            .collect();
        let mut mode_settings = settings();
        mode_settings.newton.absolute_tolerance = 1e-13;
        let initial = solid(
            &nodes,
            &elements,
            &material,
            VolumetricTreatment::OriginalMaterial,
            cx,
        )
        .advance(
            &loaded,
            &[],
            [0.0; 3],
            None,
            mode_settings,
            controls().material,
            cx,
        )
        .unwrap()
        .solid;
        let mut motion = MixedPlasticTetMotion::try_new(
            initial,
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        // u_x=4*q*y*z. Independent simplex integrals give these P2 mass
        // and infinitesimal shear stiffness. Hencky corrections are O(q^3)
        // in force and negligible at the declared 1e-5 amplitude.
        let modal_mass = 8.0 * DENSITY * VOLUME / 105.0;
        let stiffness = 16.0 * MU * VOLUME / 5.0;
        let (mut q, mut velocity) = (amplitude, 0.0);
        let dt = 0.01;
        for _ in 0..12 {
            let old_q = q;
            q = (modal_mass * q + modal_mass * dt * velocity) / (modal_mass + dt * dt * stiffness);
            velocity = (q - old_q) / dt;
            let step = motion
                .advance(
                    dt,
                    &fixed,
                    &[],
                    [0.0; 3],
                    None,
                    mode_settings,
                    controls(),
                    cx,
                )
                .unwrap();
            near(step.motion.solid().displacement_m()[node][0], q, 5e-11);
            near(step.motion.velocity_mps()[node][0], velocity, 5e-9);
            near(
                step.material_energy.stored_after_j,
                0.5 * stiffness * q * q,
                2e-14,
            );
            near(
                step.inertia.kinetic_after_j,
                0.5 * modal_mass * velocity * velocity,
                2e-14,
            );
            assert_eq!(step.material_energy.inelastic_dissipation_j(), 0.0);
            assert!(step.inertia.numerical_loss_j >= 0.0);
            work_identity(&step);
            motion = step.motion;
        }
        assert!(motion.solid().displacement_m()[node][0] < 0.0);
    });
}

fn uniaxial_supports(nodes: &[[f64; 3]], axial: f64) -> Vec<(usize, f64)> {
    let mut fixed = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        fixed.push((3 * node, axial.exp_m1() * x[0]));
        for axis in 1..3 {
            if x[axis] == 0.0 {
                fixed.push((3 * node + axis, 0.0));
            }
        }
    }
    fixed
}

#[test]
fn combined_creep_and_plastic_motion_preserve_heat_and_accepted_weak_volume() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let quadrature = duffy_quadrature(4);
        let mut trajectories = Vec::new();
        for (treatment, bulk) in [
            (VolumetricTreatment::ExactIncompressible, 1e308),
            (VolumetricTreatment::OriginalMaterial, 1e8),
        ] {
            let material = FiniteNortonJ2::new(law(bulk), kinetics(), DT, 3.0).unwrap();
            let mut motion = MixedPlasticTetMotion::try_new(
                solid(&nodes, &elements, &material, treatment, cx),
                &[0.01 * DENSITY],
                &vec![[0.0; 3]; nodes.len()],
                MASS_BYTES,
                cx,
            )
            .unwrap();
            let mut trajectory = Vec::new();
            let (mut plastic_heat, mut creep_heat, mut numerical_loss) = (0.0, 0.0, 0.0);
            for axial in [0.04, 0.08, 0.08] {
                let before = motion.solid().quadrature_states().to_vec();
                let step = motion
                    .advance(
                        DT,
                        &uniaxial_supports(&nodes, axial),
                        &[],
                        [0.0; 3],
                        None,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                let after = step.motion.solid();
                let (mut expected_plastic, mut expected_creep) = (0.0, 0.0);
                let mut constraint = [0.0; 4];
                for (q, (lambda, weight)) in quadrature.iter().enumerate() {
                    let f = &after.deformation_gradients()[q];
                    let p = &after.first_piola_pa()[q];
                    let mut tau = [0.0; 9];
                    for i in 0..3 {
                        for j in 0..3 {
                            tau[3 * i + j] = (0..3).map(|a| p[3 * i + a] * f[3 * j + a]).sum();
                        }
                    }
                    let mean = (tau[0] + tau[4] + tau[8]) / 3.0;
                    for i in [0, 4, 8] {
                        tau[i] -= mean;
                    }
                    let equivalent = (1.5 * tau.iter().map(|x| x * x).sum::<f64>()).sqrt();
                    let state = &after.quadrature_states()[q];
                    let plastic = state.alpha - before[q].alpha;
                    let creep = state.creep_alpha - before[q].creep_alpha;
                    near(creep, DT * equivalent / ETA, 3e-10);
                    if plastic > 1e-10 {
                        near(equivalent, YIELD + HARDENING * state.alpha, 3e-8);
                    }
                    expected_plastic += weight * YIELD * plastic;
                    expected_creep += weight * equivalent * creep;
                    for a in 0..4 {
                        constraint[a] += weight * lambda[a] * (determinant(f) - 1.0);
                    }
                }
                near(
                    step.material_energy.plastic_dissipation_j,
                    expected_plastic,
                    2e-10,
                );
                near(
                    step.material_energy.creep_dissipation_j,
                    expected_creep,
                    2e-10,
                );
                near(
                    step.material_energy.flow_potential_j,
                    expected_plastic + 0.5 * expected_creep,
                    2e-10,
                );
                if treatment == VolumetricTreatment::ExactIncompressible {
                    for residual in constraint {
                        near(residual, 0.0, 2e-10);
                    }
                    assert_eq!(after.bulk_energy_j(), 0.0);
                }
                near(step.material_energy.pressure_constraint_work_j, 0.0, 2e-8);
                plastic_heat += step.material_energy.plastic_dissipation_j;
                creep_heat += step.material_energy.creep_dissipation_j;
                numerical_loss += step.inertia.numerical_loss_j;
                work_identity(&step);
                trajectory.push((
                    after.displacement_m().to_vec(),
                    step.motion.velocity_mps().to_vec(),
                ));
                motion = step.motion;
            }
            assert!(plastic_heat > 1e-3 && creep_heat > 1e-4 && numerical_loss > 1e-4);
            trajectories.push(trajectory);
        }
        for (exact, finite) in trajectories[0].iter().zip(&trajectories[1]) {
            for (a, b) in exact.0.iter().flatten().zip(finite.0.iter().flatten()) {
                near(*a, *b, 2e-5);
            }
            for (a, b) in exact.1.iter().flatten().zip(finite.1.iter().flatten()) {
                near(*a, *b, 2e-3);
            }
        }
    });
}

#[test]
fn kinetic_clock_mismatch_loss_refusal_and_cancellation_preserve_complete_motion() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = FiniteOverstressJ2::new(law(160.0), kinetics(), DT).unwrap();
        let diagonal = [0.08_f64.exp(), (-0.04_f64).exp(), (-0.04_f64).exp()];
        let fixed: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(node, x)| (0..3).map(move |i| (3 * node + i, (diagonal[i] - 1.0) * x[i])))
            .collect();
        let loaded = solid(
            &nodes,
            &elements,
            &material,
            VolumetricTreatment::OriginalMaterial,
            cx,
        )
        .advance(
            &fixed,
            &[],
            [0.0; 3],
            None,
            settings(),
            controls().material,
            cx,
        )
        .unwrap()
        .solid;
        assert!(
            loaded
                .quadrature_states()
                .iter()
                .all(|state| state.alpha > 0.0)
        );
        let motion = MixedPlasticTetMotion::try_new(
            loaded,
            &[DENSITY],
            &vec![[0.0; 3]; nodes.len()],
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let snapshot = || {
            (
                motion.solid().displacement_m().to_vec(),
                motion.solid().pressure_pa().to_vec(),
                motion.solid().quadrature_states().to_vec(),
                motion.solid().deformation_gradients().to_vec(),
                motion.solid().first_piola_pa().to_vec(),
                motion.solid().stored_energy_j(),
                motion.velocity_mps().to_vec(),
                motion.reference_mass_kg(),
            )
        };
        let before = snapshot();
        let moved: Vec<_> = fixed
            .iter()
            .map(|&(d, u)| (d, u + DT * [0.1, -0.05, 0.025][d % 3]))
            .collect();
        assert!(matches!(
            motion.advance(
                2.0 * DT,
                &moved,
                &[],
                [0.0; 3],
                None,
                settings(),
                controls(),
                cx
            ),
            Err(MixedPlasticTetError::Native(
                MixedUpError::InvalidInput { .. }
            ))
        ));
        let strict = FinitePlasticMotionControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            motion.advance(DT, &moved, &[], [0.0; 3], None, settings(), strict, cx),
            Err(MixedPlasticTetError::EnergyGate {
                quantity: "backward-Euler numerical kinetic loss",
                ..
            })
        ));
        let cancelled = CancelGate::new();
        cancelled.request();
        with_gate(&cancelled, |cx| {
            assert!(matches!(
                motion.advance(DT, &moved, &[], [0.0; 3], None, settings(), controls(), cx),
                Err(MixedPlasticTetError::Native(MixedUpError::Cancelled))
            ));
        });
        assert_eq!(snapshot(), before);
        assert_eq!(motion.solid().duration_s(), Some(DT));
        let retry = motion
            .advance(DT, &moved, &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        let repeat = motion
            .advance(DT, &moved, &[], [0.0; 3], None, settings(), controls(), cx)
            .unwrap();
        assert_eq!(
            retry.motion.solid().quadrature_states(),
            repeat.motion.solid().quadrature_states()
        );
        assert_eq!(
            retry.motion.solid().displacement_m(),
            repeat.motion.solid().displacement_m()
        );
        assert_eq!(
            retry.motion.solid().pressure_pa(),
            repeat.motion.solid().pressure_pa()
        );
        assert_eq!(retry.motion.velocity_mps(), repeat.motion.velocity_mps());
        assert_eq!(retry.material_energy, repeat.material_energy);
        assert_eq!(retry.inertia, repeat.inertia);
        assert!(retry.material_energy.plastic_dissipation_j > 0.0);
        assert!(retry.inertia.numerical_loss_j > 0.0);
        work_identity(&retry);
        assert_eq!(snapshot(), before);
    });
}
