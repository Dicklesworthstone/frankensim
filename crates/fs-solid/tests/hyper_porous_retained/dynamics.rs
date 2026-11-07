//! Independent finite-compaction and consistent-momentum oracles.

use super::*;
use fs_solid::hyper3::plastic::porous::dynamics::{
    FinitePlasticMotionControls, FinitePorousTetMotion,
};

const RHO: f64 = 6.0;
const VOLUME: f64 = 1.0 / 6.0;

fn motion_controls() -> FinitePlasticMotionControls {
    FinitePlasticMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn settings() -> HyperTetSettings {
    HyperTetSettings {
        relative_tolerance: 1e-11,
        ..HyperTetSettings::default()
    }
}

fn solid<'a>(
    nodes: &'a [[f64; 3]],
    material: &'a FinitePorous,
    cx: &Cx<'_>,
) -> FinitePorousTetSolid<'a> {
    FinitePorousTetSolid::try_new(
        nodes,
        &TETS,
        FinitePlasticLawField::Uniform(material),
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}

// Exact integration of rho*N_a*N_b on any linear tetrahedron:
// M_ab = mass/20 * (1 + delta_ab), independently of the production mass owner.
fn kinetic(velocity: &[[f64; 3]], mass: f64) -> f64 {
    0.5 * mass / 20.0
        * (0..3)
            .map(|i| {
                velocity.iter().map(|v| v[i] * v[i]).sum::<f64>()
                    + velocity.iter().map(|v| v[i]).sum::<f64>().powi(2)
            })
            .sum::<f64>()
}

#[test]
fn finite_hydrostatic_motion_matches_native_compaction_and_consistent_mass_work() {
    with_cx(|cx, _| {
        let material = law(A);
        let stretch = 0.9;
        let dt = 0.25;
        let mass = RHO * VOLUME;
        let (alpha, plastic_j, tau) = hydrostatic(stretch);
        let nominal = tau / stretch;
        let velocity: Vec<_> = NODES
            .iter()
            .map(|x| x.map(|component| (stretch - 1.0) * component / dt))
            .collect();
        let gradients = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let mut forces = Vec::new();
        for (node, gradient) in gradients.iter().enumerate() {
            for component in 0..3 {
                let acceleration_sum = velocity.iter().map(|v| v[component] / dt).sum::<f64>();
                let inertial = mass / 20.0 * (acceleration_sum + velocity[node][component] / dt);
                forces.push((
                    3 * node + component,
                    VOLUME * nominal * gradient[component] + inertial,
                ));
            }
        }
        let prescribed = affine(
            &NODES,
            [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch],
            [0.0; 3],
        );
        let expected_storage = VOLUME * (0.5 * tau * tau / K + 0.5 * H * alpha * alpha);
        let expected_heat = VOLUME * Y * alpha;
        let expected_kinetic = kinetic(&velocity, mass);
        let expected_internal_work = 3.0 * VOLUME * nominal * (stretch - 1.0);

        // The same endpoint is reached by a free dynamic solve and by prescribed
        // motion. The latter checks all off-diagonal inertial support reactions.
        for driven_by_supports in [false, true] {
            let motion = FinitePorousTetMotion::try_new(
                solid(&NODES, &material, cx),
                &[RHO],
                &[[0.0; 3]; 4],
                cx,
            )
            .unwrap();
            let step = motion
                .advance(
                    if driven_by_supports { &prescribed } else { &[] },
                    if driven_by_supports { &[] } else { &forces },
                    [0.0; 3],
                    dt,
                    settings(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            if driven_by_supports {
                for ((dof, reaction), (expected_dof, force)) in
                    step.equilibrium.reactions_n.iter().zip(&forces)
                {
                    assert_eq!(dof, expected_dof);
                    near(*reaction, *force, 3e-7);
                }
                assert_eq!(step.equilibrium.reactions_n.len(), forces.len());
            } else {
                assert!(step.equilibrium.reactions_n.is_empty());
                assert!(step.equilibrium.steps[0].linear_iterations > 0);
            }
            for (actual, expected) in step.motion.velocity_mps().iter().zip(&velocity) {
                for component in 0..3 {
                    near(actual[component], expected[component], 3e-8);
                }
            }
            for (component, value) in step.motion.solid().element_deformation()[0]
                .iter()
                .enumerate()
            {
                near(*value, if component % 4 == 0 { stretch } else { 0.0 }, 3e-8);
            }
            near(step.motion.solid().element_states()[0].alpha, alpha, 3e-8);
            near(
                step.motion.solid().element_plastic_jacobians()[0],
                plastic_j,
                3e-8,
            );
            near(
                step.motion.solid().element_jacobians()[0],
                stretch.powi(3),
                3e-8,
            );
            near(
                step.motion.current_density_kg_m3(cx).unwrap()[0],
                RHO / stretch.powi(3),
                3e-7,
            );
            assert!((RHO / plastic_j - RHO / stretch.powi(3)).abs() > 0.1);
            near(step.motion.reference_mass_kg(), mass, 1e-14);
            assert_eq!(step.motion.reference_density_kg_m3(), [RHO]);
            near(step.material_energy.stored_after_j, expected_storage, 3e-7);
            near(
                step.material_energy.plastic_dissipation_j,
                expected_heat,
                3e-7,
            );
            near(
                step.material_energy.internal_endpoint_work_j,
                expected_internal_work,
                3e-7,
            );
            near(step.inertia.kinetic_before_j, 0.0, 0.0);
            near(step.inertia.kinetic_after_j, expected_kinetic, 3e-8);
            near(step.inertia.numerical_loss_j, expected_kinetic, 3e-8);
            near(step.inertia.inertial_work_j, 2.0 * expected_kinetic, 3e-8);
            near(
                step.material_energy.dead_load_endpoint_work_j
                    + step.material_energy.support_endpoint_work_j,
                expected_storage
                    + expected_heat
                    + step.material_energy.endpoint_work_defect_j
                    + 2.0 * expected_kinetic,
                5e-7,
            );
            near(step.balance_error_j, 0.0, 2e-7);
            assert_eq!(motion.solid().element_states(), [material.initial_state()]);
            assert_eq!(motion.velocity_mps(), [[0.0; 3]; 4]);
        }
    });
}

#[test]
fn unrestrained_release_recovers_elastic_volume_and_keeps_permanent_compaction() {
    with_cx(|cx, _| {
        // A centered regular tetrahedron makes free dilation an exact mode.
        // Its gradients are X_a/4, V0=8/3, and sum X_a=0.
        let nodes = [
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ];
        let material = law(A);
        let initial_stretch = 0.9;
        let dt = 0.02;
        let (_, plastic_j, _) = hydrostatic(initial_stretch);
        let compacted = solid(&nodes, &material, cx)
            .advance(
                &affine(
                    &nodes,
                    [
                        initial_stretch,
                        0.0,
                        0.0,
                        0.0,
                        initial_stretch,
                        0.0,
                        0.0,
                        0.0,
                        initial_stretch,
                    ],
                    [0.0; 3],
                ),
                &[],
                [0.0; 3],
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let accepted_state = compacted.element_states().to_vec();
        let motion = FinitePorousTetMotion::try_new(compacted, &[RHO], &[[0.0; 3]; 4], cx).unwrap();

        // Original momentum reduces to rho*(s-s0)/(5*dt²)+P(s)=0.
        // Elastic unloading retains Jp, so P=K*(3 ln(s)-ln(Jp))/s.
        let mut lower = initial_stretch;
        let mut upper = plastic_j.cbrt();
        for _ in 0..64 {
            let middle = 0.5 * (lower + upper);
            let residual = RHO * (middle - initial_stretch) / (5.0 * dt * dt)
                + K * (3.0 * middle.ln() - plastic_j.ln()) / middle;
            if residual < 0.0 {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        let recovered = 0.5 * (lower + upper);
        assert!(recovered > initial_stretch && recovered < plastic_j.cbrt());
        let step = motion
            .advance(&[], &[], [0.0; 3], dt, settings(), motion_controls(), cx)
            .unwrap();
        assert!(step.equilibrium.reactions_n.is_empty());
        assert_eq!(step.motion.solid().element_states(), accepted_state);
        near(
            step.motion.solid().element_plastic_jacobians()[0],
            plastic_j,
            3e-8,
        );
        near(
            step.motion.solid().element_jacobians()[0],
            recovered.powi(3),
            3e-8,
        );
        for (x, (u, v)) in nodes.iter().zip(
            step.motion
                .solid()
                .displacement_m()
                .iter()
                .zip(step.motion.velocity_mps()),
        ) {
            for component in 0..3 {
                near(u[component], (recovered - 1.0) * x[component], 3e-8);
                near(
                    v[component],
                    (recovered - initial_stretch) * x[component] / dt,
                    3e-7,
                );
            }
        }
        let mass = RHO * 8.0 / 3.0;
        let current_density = step.motion.current_density_kg_m3(cx).unwrap()[0];
        near(current_density * recovered.powi(3) * 8.0 / 3.0, mass, 3e-7);
        near(step.motion.reference_mass_kg(), mass, 1e-13);
        near(step.inertia.reference_mass_kg, mass, 1e-13);
        assert!(step.inertia.kinetic_after_j > 0.0);
        near(
            step.inertia.kinetic_after_j,
            kinetic(step.motion.velocity_mps(), mass),
            3e-9,
        );
        near(
            step.inertia.numerical_loss_j,
            step.inertia.kinetic_after_j,
            3e-9,
        );
        near(step.material_energy.plastic_dissipation_j, 0.0, 0.0);
        near(step.material_energy.dead_load_endpoint_work_j, 0.0, 0.0);
        near(step.material_energy.support_endpoint_work_j, 0.0, 0.0);
        near(
            step.material_energy.stored_after_j - step.material_energy.stored_before_j
                + step.material_energy.endpoint_work_defect_j
                + step.inertia.kinetic_change_j
                + step.inertia.numerical_loss_j,
            0.0,
            2e-7,
        );
        near(step.balance_error_j, 0.0, 2e-7);
        assert_eq!(motion.solid().element_states(), accepted_state);
        assert_eq!(motion.velocity_mps(), [[0.0; 3]; 4]);
    });
}

#[test]
fn free_rigid_travel_and_acceleration_do_not_create_compaction_or_heat() {
    with_cx(|cx, _| {
        let material = law(A);
        let initial_velocity = [0.2, -0.1, 0.05];
        let dt = 0.02;
        let mass = RHO * VOLUME;
        for acceleration in [[0.0; 3], [1.0, -0.5, 0.25]] {
            let mut motion = FinitePorousTetMotion::try_new(
                solid(&NODES, &material, cx),
                &[RHO],
                &[initial_velocity; 4],
                cx,
            )
            .unwrap();
            for interval in 1..=3 {
                let n = f64::from(interval);
                let expected_v: [f64; 3] =
                    core::array::from_fn(|i| initial_velocity[i] + n * dt * acceleration[i]);
                let expected_u: [f64; 3] = core::array::from_fn(|i| {
                    n * dt * initial_velocity[i] + 0.5 * n * (n + 1.0) * dt * dt * acceleration[i]
                });
                let step = motion
                    .advance(
                        &[],
                        &[],
                        acceleration.map(|a| RHO * a),
                        dt,
                        settings(),
                        motion_controls(),
                        cx,
                    )
                    .unwrap();
                for (u, v) in step
                    .motion
                    .solid()
                    .displacement_m()
                    .iter()
                    .zip(step.motion.velocity_mps())
                {
                    for i in 0..3 {
                        near(u[i], expected_u[i], 1e-9);
                        near(v[i], expected_v[i], 1e-8);
                    }
                }
                assert_eq!(
                    step.motion.solid().element_states(),
                    [material.initial_state()]
                );
                near(step.motion.solid().element_jacobians()[0], 1.0, 1e-10);
                near(step.motion.solid().element_plastic_jacobians()[0], 1.0, 0.0);
                near(step.motion.reference_mass_kg(), mass, 1e-14);
                near(step.material_energy.plastic_dissipation_j, 0.0, 0.0);
                near(step.material_energy.stored_after_j, 0.0, 1e-12);
                near(
                    step.inertia.kinetic_after_j,
                    0.5 * mass * expected_v.iter().map(|v| v * v).sum::<f64>(),
                    1e-9,
                );
                near(
                    step.inertia.numerical_loss_j,
                    0.5 * mass * dt * dt * acceleration.iter().map(|a| a * a).sum::<f64>(),
                    1e-9,
                );
                near(step.balance_error_j, 0.0, 1e-9);
                motion = step.motion;
            }
        }
    });
}

#[test]
fn kinetic_loss_refusal_and_cancellation_keep_full_compacted_motion_retryable() {
    with_cx(|cx, gate| {
        let material = law(A);
        let held = affine(
            &NODES,
            [0.9, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.9],
            [0.0; 3],
        );
        let compacted = solid(&NODES, &material, cx)
            .advance(&held, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap()
            .solid;
        let original_history = compacted.element_states().to_vec();
        let original_jp = compacted.element_plastic_jacobians().to_vec();
        let original_u = compacted.displacement_m().to_vec();
        let original_storage = compacted.stored_energy_j();
        let velocity = [[0.2, -0.1, 0.05]; 4];
        let motion = FinitePorousTetMotion::try_new(compacted, &[RHO], &velocity, cx).unwrap();
        let mut strict = motion_controls();
        strict.maximum_numerical_kinetic_loss_j = 0.0;
        let refused = motion
            .advance(&held, &[], [0.0; 3], 0.1, settings(), strict, cx)
            .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate {
                quantity: "porous backward-Euler numerical kinetic loss",
                ..
            })
        ));
        let retry = motion
            .advance(&held, &[], [0.0; 3], 0.1, settings(), motion_controls(), cx)
            .unwrap();
        let repeat = motion
            .advance(&held, &[], [0.0; 3], 0.1, settings(), motion_controls(), cx)
            .unwrap();
        assert_eq!(retry.equilibrium, repeat.equilibrium);
        assert_eq!(retry.material_energy, repeat.material_energy);
        assert_eq!(retry.inertia, repeat.inertia);
        assert_eq!(retry.motion.solid().element_states(), original_history);
        assert_eq!(retry.motion.velocity_mps(), [[0.0; 3]; 4]);
        near(retry.material_energy.plastic_dissipation_j, 0.0, 0.0);
        near(
            retry.inertia.numerical_loss_j,
            kinetic(&velocity, RHO * VOLUME),
            1e-12,
        );
        assert!(
            motion
                .advance(&held, &[], [0.0; 3], 0.0, settings(), motion_controls(), cx)
                .is_err()
        );
        gate.request();
        let cancelled = motion
            .advance(&held, &[], [0.0; 3], 0.1, settings(), motion_controls(), cx)
            .unwrap_err();
        assert!(matches!(
            cancelled,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(HyperTetError::Cancelled))
        ));
        assert_eq!(motion.solid().element_states(), original_history);
        assert_eq!(motion.solid().element_plastic_jacobians(), original_jp);
        assert_eq!(motion.solid().displacement_m(), original_u);
        assert_eq!(motion.velocity_mps(), velocity);
        near(motion.solid().stored_energy_j(), original_storage, 0.0);
        near(motion.reference_mass_kg(), RHO * VOLUME, 1e-14);
    });
}

// A regular reference tet centered at zero has grad(N_a) = X_a/4 and V0 = 8/3.
// For uniform contraction, M*X = rho*V0*X/20, so the pressure needed to reach
// stretch s from rest is -tau(s)/s^3 + rho*(1-s)/(5*dt^2*s^2).
const PRESSURE_NODES: [[f64; 3]; 4] = [
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [-1.0, 1.0, -1.0],
];
const PRESSURE_VOLUME: f64 = 8.0 / 3.0;

#[test]
fn current_area_pressure_drives_free_compaction_with_consistent_inertial_work() {
    with_cx(|cx, _| {
        let material = law(A);
        let stretch = 0.9;
        let dt = 0.25;
        let (alpha, plastic_j, tau) = hydrostatic(stretch);
        let pressure =
            -tau / stretch.powi(3) + RHO * (1.0 - stretch) / (5.0 * dt * dt * stretch * stretch);
        let surface = ClosedTetPressureSurface::new(
            &PRESSURE_NODES,
            &TETS,
            TetAssemblyBudget::standard(),
            4,
            cx,
        )
        .unwrap();
        let mut face_pressures = vec![0.0; surface.faces().len()];
        let patch = surface
            .faces()
            .iter()
            .position(|face| !face.contains(&0))
            .unwrap();
        face_pressures[patch] = pressure;
        let expected_kinetic =
            3.0 * RHO * PRESSURE_VOLUME * (stretch - 1.0).powi(2) / (10.0 * dt * dt);
        let expected_storage = PRESSURE_VOLUME * (0.5 * tau * tau / K + 0.5 * H * alpha * alpha);
        let expected_heat = PRESSURE_VOLUME * Y * alpha;
        let expected_internal = 3.0 * PRESSURE_VOLUME * tau / stretch * (stretch - 1.0);

        // A one-face patch uses the nonsymmetric native pressure Jacobian. Its
        // independent manufactured dead loads restore the same affine endpoint,
        // without supports or using the production pressure force as an oracle.
        for patch_only in [false, true] {
            let mut dead_loads = Vec::new();
            if patch_only {
                let required = PRESSURE_VOLUME * tau / (4.0 * stretch)
                    + RHO * PRESSURE_VOLUME * (stretch - 1.0) / (20.0 * dt * dt);
                for (node, x) in PRESSURE_NODES.iter().enumerate() {
                    // Face opposite X_0 has outward reference area vector -2 X_0.
                    let patch_force = if node == 0 {
                        0.0
                    } else {
                        2.0 * pressure * stretch * stretch / 3.0
                    };
                    for (component, coordinate) in x.iter().enumerate() {
                        dead_loads
                            .push((3 * node + component, required * coordinate - patch_force));
                    }
                }
            }
            let fraction = if patch_only { 0.25 } else { 1.0 };
            let motion = FinitePorousTetMotion::try_new(
                solid(&PRESSURE_NODES, &material, cx),
                &[RHO],
                &[[0.0; 3]; 4],
                cx,
            )
            .unwrap();
            let step = motion
                .advance_with_pressure(
                    &[],
                    &dead_loads,
                    [0.0; 3],
                    dt,
                    &surface,
                    if patch_only {
                        TetPressure::PerFace(&face_pressures)
                    } else {
                        TetPressure::Uniform(pressure)
                    },
                    settings(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            assert!(step.equilibrium.reactions_n.is_empty());
            assert!(step.equilibrium.steps[0].linear_iterations > 0);
            for (node, x) in PRESSURE_NODES.iter().enumerate() {
                for (component, coordinate) in x.iter().enumerate() {
                    near(
                        step.motion.solid().displacement_m()[node][component],
                        (stretch - 1.0) * coordinate,
                        3e-8,
                    );
                    near(
                        step.motion.velocity_mps()[node][component],
                        (stretch - 1.0) * coordinate / dt,
                        1e-7,
                    );
                }
            }
            for (component, value) in step.motion.solid().element_deformation()[0]
                .iter()
                .enumerate()
            {
                near(*value, if component % 4 == 0 { stretch } else { 0.0 }, 3e-8);
            }
            near(step.motion.solid().element_states()[0].alpha, alpha, 3e-8);
            near(
                step.motion.solid().element_plastic_jacobians()[0],
                plastic_j,
                3e-8,
            );
            near(
                step.motion.current_density_kg_m3(cx).unwrap()[0],
                RHO / stretch.powi(3),
                3e-7,
            );
            near(
                step.motion.reference_mass_kg(),
                RHO * PRESSURE_VOLUME,
                1e-12,
            );
            near(step.material_energy.stored_after_j, expected_storage, 3e-7);
            near(
                step.material_energy.plastic_dissipation_j,
                expected_heat,
                3e-7,
            );
            near(
                step.material_energy.internal_endpoint_work_j,
                expected_internal,
                3e-7,
            );
            near(step.inertia.kinetic_change_j, expected_kinetic, 3e-7);
            near(step.inertia.numerical_loss_j, expected_kinetic, 3e-7);
            near(step.inertia.inertial_work_j, 2.0 * expected_kinetic, 5e-7);
            let account = step.pressure.unwrap();
            let endpoint =
                -fraction * 3.0 * pressure * PRESSURE_VOLUME * stretch * stretch * (stretch - 1.0);
            let exact = -fraction * pressure * PRESSURE_VOLUME * (stretch.powi(3) - 1.0);
            near(account.endpoint_work_j, endpoint, 3e-7);
            near(account.work_j, exact, 3e-7);
            near(account.endpoint_work_correction_j, exact - endpoint, 3e-7);
            assert!(account.endpoint_work_correction_j.abs() > 0.1);
            near(
                account.external_work_defect_j,
                step.material_energy.endpoint_work_defect_j + exact - endpoint,
                3e-7,
            );
            // The combined pressure defect already includes the material defect:
            // counting both again, or turning the kinetic loss into heat, fails.
            near(
                step.material_energy.dead_load_endpoint_work_j
                    + step.material_energy.support_endpoint_work_j
                    + account.work_j,
                step.material_energy.stored_after_j - step.material_energy.stored_before_j
                    + step.material_energy.plastic_dissipation_j
                    + step.inertia.kinetic_change_j
                    + step.inertia.numerical_loss_j
                    + account.external_work_defect_j
                    + step.balance_error_j,
                5e-7,
            );
            near(step.balance_error_j, 0.0, 2e-7);
            assert_eq!(motion.solid().element_states(), [material.initial_state()]);
            assert_eq!(motion.velocity_mps(), [[0.0; 3]; 4]);
        }
    });
}

#[test]
fn pressure_motion_supports_and_late_refusal_preserve_full_compacted_history() {
    with_cx(|cx, gate| {
        let material = law(A);
        let old_stretch = 0.95;
        let stretch = 0.9;
        let dt = 0.25;
        let pressure = 5.0;
        let old_supports = affine(
            &NODES,
            [
                old_stretch,
                0.0,
                0.0,
                0.0,
                old_stretch,
                0.0,
                0.0,
                0.0,
                old_stretch,
            ],
            [0.0; 3],
        );
        let supports = affine(
            &NODES,
            [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch],
            [0.0; 3],
        );
        let compacted = solid(&NODES, &material, cx)
            .advance(&old_supports, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap()
            .solid;
        let original_history = compacted.element_states().to_vec();
        let original_u = compacted.displacement_m().to_vec();
        let original_storage = compacted.stored_energy_j();
        let old_velocity = [0.2, -0.1, 0.05];
        let motion =
            FinitePorousTetMotion::try_new(compacted, &[RHO], &[old_velocity; 4], cx).unwrap();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let mut face_pressures = vec![0.0; surface.faces().len()];
        let patch = surface
            .faces()
            .iter()
            .position(|face| !face.contains(&0))
            .unwrap();
        face_pressures[patch] = pressure;
        let field = TetPressure::PerFace(&face_pressures);
        let (alpha, plastic_j, tau) = hydrostatic(stretch);
        let (old_alpha, _, _) = hydrostatic(old_stretch);
        let velocity: Vec<_> = NODES
            .iter()
            .map(|x| x.map(|value| (stretch - old_stretch) * value / dt))
            .collect();
        let acceleration: Vec<[f64; 3]> = velocity
            .iter()
            .map(|v| {
                core::array::from_fn(|component| (v[component] - old_velocity[component]) / dt)
            })
            .collect();
        let gradients = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let step = motion
            .advance_with_pressure(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &surface,
                field,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        for &(dof, actual) in &step.equilibrium.reactions_n {
            let node = dof / 3;
            let component = dof % 3;
            let inertia = RHO * VOLUME / 20.0
                * (acceleration.iter().map(|a| a[component]).sum::<f64>()
                    + acceleration[node][component]);
            // The sole loaded face has outward area vector [1,1,1]/2.
            let pressure_force = if node == 0 {
                0.0
            } else {
                -pressure * stretch * stretch / 6.0
            };
            let expected =
                VOLUME * tau / stretch * gradients[node][component] + inertia - pressure_force;
            near(actual, expected, 3e-7);
        }
        assert_eq!(step.equilibrium.reactions_n.len(), 12);
        near(step.motion.solid().element_states()[0].alpha, alpha, 3e-8);
        near(
            step.motion.solid().element_plastic_jacobians()[0],
            plastic_j,
            3e-8,
        );
        near(
            step.material_energy.plastic_dissipation_j,
            VOLUME * Y * (alpha - old_alpha),
            3e-8,
        );
        near(
            step.inertia.kinetic_change_j,
            kinetic(&velocity, RHO * VOLUME) - kinetic(&[old_velocity; 4], RHO * VOLUME),
            3e-8,
        );
        let account = step.pressure.unwrap();
        near(
            account.work_j,
            -pressure * VOLUME * (stretch.powi(3) - old_stretch.powi(3)),
            3e-8,
        );
        near(
            account.endpoint_work_j,
            -3.0 * pressure * VOLUME * stretch * stretch * (stretch - old_stretch),
            3e-8,
        );
        near(
            step.material_energy.support_endpoint_work_j + account.work_j,
            step.material_energy.stored_after_j - step.material_energy.stored_before_j
                + step.material_energy.plastic_dissipation_j
                + step.inertia.kinetic_change_j
                + step.inertia.numerical_loss_j
                + account.external_work_defect_j
                + step.balance_error_j,
            3e-7,
        );

        let mut strict = motion_controls();
        strict.maximum_numerical_kinetic_loss_j = 0.0;
        let refused = motion
            .advance_with_pressure(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &surface,
                field,
                settings(),
                strict,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate {
                quantity: "porous backward-Euler numerical kinetic loss",
                ..
            })
        ));
        let retry = motion
            .advance_with_pressure(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &surface,
                field,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.equilibrium, step.equilibrium);
        assert_eq!(retry.material_energy, step.material_energy);
        assert_eq!(retry.pressure, step.pressure);
        assert_eq!(retry.inertia, step.inertia);

        let dead = motion
            .advance(
                &supports,
                &[],
                [0.0; 3],
                dt,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        let zero = motion
            .advance_with_pressure(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &surface,
                0.0,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(zero.equilibrium, dead.equilibrium);
        assert_eq!(zero.material_energy, dead.material_energy);
        assert_eq!(zero.inertia, dead.inertia);
        assert!(dead.pressure.is_none());
        near(zero.pressure.unwrap().work_j, 0.0, 0.0);
        near(zero.pressure.unwrap().endpoint_work_correction_j, 0.0, 0.0);

        gate.request();
        let cancelled = motion
            .advance_with_pressure(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &surface,
                field,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            cancelled,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(HyperTetError::Cancelled))
        ));
        assert_eq!(motion.solid().element_states(), original_history);
        assert_eq!(motion.solid().displacement_m(), original_u);
        assert_eq!(motion.velocity_mps(), [old_velocity; 4]);
        near(motion.solid().stored_energy_j(), original_storage, 0.0);
        near(motion.reference_mass_kg(), RHO * VOLUME, 1e-14);
    });
}
