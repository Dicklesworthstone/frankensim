//! G0/G1/G3: actual thin-solid momentum, memory and physical energy accounts.
use super::*;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_solid::hyper2d::visco::dynamics::{
    FiniteViscoMotionControls, FiniteViscoPlaneMotion, FiniteViscoPlaneMotionStep,
};

const MASS_BYTES: usize = 1 << 20;

fn motion_controls() -> FiniteViscoMotionControls {
    FiniteViscoMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 10.0,
    }
}

#[derive(Debug, PartialEq)]
struct MotionSnapshot {
    solid: Snapshot,
    velocity: Vec<[f64; 2]>,
    density: Vec<f64>,
    mass: f64,
}

fn motion_snapshot(motion: &FiniteViscoPlaneMotion<'_>) -> MotionSnapshot {
    MotionSnapshot {
        solid: snapshot(motion.solid()),
        velocity: motion.velocity_mps().to_vec(),
        density: motion.reference_density_kg_m3().to_vec(),
        mass: motion.reference_mass_kg(),
    }
}

fn check_motion_balance(step: &FiniteViscoPlaneMotionStep<'_>) {
    let material = &step.material_energy;
    let inertia = &step.inertia;
    assert!(accounts(material).iter().all(|x| x.is_finite()));
    assert!(material.dissipated_j >= 0.0);
    assert!(inertia.numerical_loss_j >= 0.0);
    close(material.material_balance_error_j, 0.0, 1e-9);
    close(material.equilibrium_work_residual_j, 0.0, 1e-9);
    close(inertia.balance_error_j, 0.0, 1e-9);
    close(step.balance_error_j, 0.0, 1e-9);
    close(
        material.inertial_endpoint_work_j,
        inertia.inertial_work_j,
        1e-10,
    );
    close(
        inertia.inertial_work_j,
        inertia.kinetic_change_j + inertia.numerical_loss_j,
        1e-10,
    );
    close(
        material.stored_after_j - material.stored_before_j
            + inertia.kinetic_change_j
            + material.dissipated_j
            + inertia.numerical_loss_j,
        material.dead_load_endpoint_work_j + material.support_endpoint_work_j
            - material.endpoint_work_defect_j
            - step.balance_error_j,
        1e-9,
    );
}

#[test]
fn g1_free_translation_needs_no_support_and_preserves_velocity_and_reference_kinetic_energy() {
    with_cx(|cx, _| {
        let law = law();
        let rho = 2.3;
        let thickness = 0.03;
        let speed = [0.7, -0.3];
        for mesh in [
            Mesh2::triangles(1.0, 0.4, 1, 1),
            Mesh2::quads(1.0, 0.4, 1, 1),
        ] {
            let mass = rho * 1.0 * 0.4 * thickness;
            let kinetic = 0.5 * mass * (speed[0] * speed[0] + speed[1] * speed[1]);
            for exact in [false, true] {
                let solid = FiniteViscoPlaneSolid::try_new(
                    &mesh,
                    ViscoelasticLawField::Uniform(&law),
                    mode(exact),
                    thickness,
                    local(),
                    cx,
                )
                .unwrap();
                let density = vec![rho; mesh.elems.len()];
                let mut motion = FiniteViscoPlaneMotion::try_new(
                    solid,
                    &density,
                    &vec![speed; mesh.nodes.len()],
                    MASS_BYTES,
                    cx,
                )
                .unwrap();
                close(motion.reference_mass_kg(), mass, 2e-13);
                let mut time = 0.0;
                for dt in [0.013, 0.017, 0.011] {
                    let before = motion_snapshot(&motion);
                    let step = motion
                        .advance(dt, &[], &[], &[], newton(), motion_controls(), cx)
                        .unwrap();
                    assert_eq!(motion_snapshot(&motion), before);
                    time += dt;
                    for (u, velocity) in step
                        .motion
                        .solid()
                        .displacement_m()
                        .iter()
                        .zip(step.motion.velocity_mps())
                    {
                        for component in 0..2 {
                            close(u[component], time * speed[component], 2e-11);
                            close(velocity[component], speed[component], 2e-10);
                        }
                    }
                    assert!(step.reactions_n.is_empty());
                    close(step.motion.solid().time_s(), time, 1e-15);
                    close(step.inertia.kinetic_before_j, kinetic, 2e-12);
                    close(step.inertia.kinetic_after_j, kinetic, 2e-12);
                    close(step.inertia.numerical_loss_j, 0.0, 1e-18);
                    close(step.material_energy.dissipated_j, 0.0, 1e-12);
                    check_motion_balance(&step);
                    motion = step.motion;
                }
            }
        }
    });
}

#[test]
fn g1_consistent_reference_acceleration_matches_backward_euler_motion_and_kinetic_work() {
    with_cx(|cx, _| {
        let law = law();
        let rho = 2.6;
        let thickness = 0.04;
        let velocity = [0.25, -0.1];
        let acceleration = [0.4, 0.2];
        let dt = 0.05;
        let mass = rho * 1.2 * 0.5 * thickness;
        let next_v = [
            velocity[0] + dt * acceleration[0],
            velocity[1] + dt * acceleration[1],
        ];
        let expected_work = dt * mass * (acceleration[0] * next_v[0] + acceleration[1] * next_v[1]);
        let expected_loss = 0.5
            * mass
            * dt
            * dt
            * (acceleration[0] * acceleration[0] + acceleration[1] * acceleration[1]);
        for mesh in [
            Mesh2::triangles(1.2, 0.5, 1, 1),
            Mesh2::quads(1.2, 0.5, 1, 1),
        ] {
            for exact in [false, true] {
                let solid = FiniteViscoPlaneSolid::try_new(
                    &mesh,
                    ViscoelasticLawField::Uniform(&law),
                    mode(exact),
                    thickness,
                    local(),
                    cx,
                )
                .unwrap();
                let motion = FiniteViscoPlaneMotion::try_new(
                    solid,
                    &vec![rho; mesh.elems.len()],
                    &vec![velocity; mesh.nodes.len()],
                    MASS_BYTES,
                    cx,
                )
                .unwrap();
                let loads = motion
                    .reference_acceleration_load_n(&vec![acceleration; mesh.nodes.len()], cx)
                    .unwrap();
                for component in 0..2 {
                    let force: f64 = loads
                        .iter()
                        .filter(|(dof, _)| dof % 2 == component)
                        .map(|(_, f)| f)
                        .sum();
                    close(force, mass * acceleration[component], 2e-13);
                }
                let before = motion_snapshot(&motion);
                let step = motion
                    .advance(dt, &[], &[], &loads, newton(), motion_controls(), cx)
                    .unwrap();
                assert_eq!(motion_snapshot(&motion), before);
                for (u, v) in step
                    .motion
                    .solid()
                    .displacement_m()
                    .iter()
                    .zip(step.motion.velocity_mps())
                {
                    for component in 0..2 {
                        close(v[component], next_v[component], 2e-10);
                        // This is the admitted first-order BE position, not the
                        // continuous constant-acceleration position.
                        close(u[component], dt * next_v[component], 2e-11);
                    }
                }
                close(step.inertia.reference_mass_kg, mass, 2e-13);
                close(step.inertia.inertial_work_j, expected_work, 2e-11);
                close(step.inertia.numerical_loss_j, expected_loss, 2e-11);
                close(
                    step.material_energy.dead_load_endpoint_work_j,
                    expected_work,
                    2e-11,
                );
                close(step.material_energy.dissipated_j, 0.0, 1e-12);
                assert!(step.inertia.numerical_loss_j > 0.0);
                assert!(step.reactions_n.is_empty());
                check_motion_balance(&step);
            }
        }
    });
}

#[test]
fn g3_reference_event_retains_original_mass_action_velocity_and_nonuniform_material_memory() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.4, 0.6, 2, 1);
        let law = law();
        let mut target = law.clone().with_reference_deformation(
            ReferenceDeformation::try_new([1.04, 0.03, 0.01, 0.0, 0.98, -0.005, 0.0, 0.0, 1.02])
                .unwrap(),
        );
        target.branches[0].tau_s *= 0.5;
        let strain = |x: f64, y: f64| [0.1 * x + 0.05 * x * y, -0.025 * y + 0.03 * x * y];
        let initial: Vec<_> = mesh.nodes.iter().map(|[x, y]| strain(*x, *y)).collect();
        let velocity: Vec<_> = mesh
            .nodes
            .iter()
            .map(|[x, y]| [0.02 + 0.01 * x, -0.03 + 0.02 * y])
            .collect();
        let probe: Vec<_> = mesh
            .nodes
            .iter()
            .map(|[x, y]| [0.3 + 0.1 * x, -0.2 + 0.05 * y])
            .collect();
        let density = [1.5, 3.0];
        let thickness = 0.03;
        let expected_mass = 0.5 * 1.4 * 0.6 * thickness * (density[0] + density[1]);
        for exact in [false, true] {
            let solid = FiniteViscoPlaneSolid::try_new_from(
                &mesh,
                ViscoelasticLawField::Uniform(&law),
                mode(exact),
                thickness,
                &initial,
                local(),
                cx,
            )
            .unwrap()
            .advance(0.04, &boundary(&strain), &[], newton(), controls(), cx)
            .unwrap()
            .solid;
            let motion =
                FiniteViscoPlaneMotion::try_new(solid, &density, &velocity, MASS_BYTES, cx)
                    .unwrap();
            let before = motion_snapshot(&motion);
            let original_action = motion.reference_acceleration_load_n(&probe, cx).unwrap();
            let event = motion
                .change_reference(ViscoelasticLawField::Uniform(&target), controls(), cx)
                .unwrap();
            assert_eq!(motion_snapshot(&motion), before);
            assert_eq!(event.motion.velocity_mps(), velocity);
            assert_eq!(event.motion.reference_density_kg_m3(), density);
            assert_eq!(event.motion.reference_mass_kg(), motion.reference_mass_kg());
            close(event.motion.reference_mass_kg(), expected_mass, 2e-13);
            assert_eq!(
                event
                    .motion
                    .reference_acceleration_load_n(&probe, cx)
                    .unwrap(),
                original_action
            );
            assert_eq!(
                event.motion.solid().displacement_m(),
                motion.solid().displacement_m()
            );
            assert_eq!(
                event.motion.solid().material_points(),
                motion.solid().material_points()
            );
            assert_eq!(event.motion.solid().time_s(), 0.04);
            if exact {
                let previous = motion.solid().point_remainder_states().unwrap();
                let next = event.motion.solid().point_remainder_states().unwrap();
                for (a, b) in previous.iter().zip(next) {
                    assert_eq!(
                        a.inverse_viscous_cauchy_green,
                        b.inverse_viscous_cauchy_green
                    );
                    assert_eq!(a.dissipated, b.dissipated);
                }
            } else {
                let previous = motion.solid().point_states().unwrap();
                let next = event.motion.solid().point_states().unwrap();
                for (a, b) in previous.iter().zip(next) {
                    assert_eq!(
                        a.inverse_viscous_cauchy_green,
                        b.inverse_viscous_cauchy_green
                    );
                    assert_eq!(a.dissipated, b.dissipated);
                }
            }
        }
    });
}

#[test]
fn g0_loaded_strain_motion_has_correct_thickness_scaling_and_late_kinetic_gate_is_atomic() {
    with_cx(|cx, _| {
        let law = law();
        let mesh = Mesh2::quads(1.0, 0.4, 2, 1);
        let clamp = |_: f64, _: f64| [0.0, 0.0];
        let pull = |_: f64, _: f64| [0.8, 0.0];
        let supports = [(Patch::Left, &clamp as Field<'_>)];
        let tractions = [(Patch::Right, &pull as Field<'_>)];
        for exact in [false, true] {
            let make = |thickness| {
                let solid = FiniteViscoPlaneSolid::try_new(
                    &mesh,
                    ViscoelasticLawField::Uniform(&law),
                    mode(exact),
                    thickness,
                    local(),
                    cx,
                )
                .unwrap();
                FiniteViscoPlaneMotion::try_new(
                    solid,
                    &vec![2.5; mesh.elems.len()],
                    &vec![[0.0; 2]; mesh.nodes.len()],
                    MASS_BYTES,
                    cx,
                )
                .unwrap()
            };
            let thin = make(0.02);
            let thick = make(0.04);
            let before = motion_snapshot(&thin);
            let first = thin
                .advance(
                    0.04,
                    &supports,
                    &tractions,
                    &[],
                    newton(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            let scaled = thick
                .advance(
                    0.04,
                    &supports,
                    &tractions,
                    &[],
                    newton(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            assert_eq!(motion_snapshot(&thin), before);
            for (a, b) in first
                .motion
                .solid()
                .displacement_m()
                .iter()
                .flatten()
                .zip(scaled.motion.solid().displacement_m().iter().flatten())
            {
                close(*a, *b, 2e-10);
            }
            for (a, b) in first
                .motion
                .velocity_mps()
                .iter()
                .flatten()
                .zip(scaled.motion.velocity_mps().iter().flatten())
            {
                close(*a, *b, 2e-10);
            }
            for (a, b) in accounts(&first.material_energy)
                .into_iter()
                .zip(accounts(&scaled.material_energy))
            {
                close(b, 2.0 * a, 2e-10);
            }
            close(
                scaled.inertia.reference_mass_kg,
                2.0 * first.inertia.reference_mass_kg,
                2e-13,
            );
            close(
                scaled.inertia.kinetic_after_j,
                2.0 * first.inertia.kinetic_after_j,
                2e-10,
            );
            close(
                scaled.inertia.numerical_loss_j,
                2.0 * first.inertia.numerical_loss_j,
                2e-10,
            );
            assert!(first.material_energy.dissipated_j > 0.0);
            assert!(first.inertia.numerical_loss_j > 0.0);
            check_motion_balance(&first);
            check_motion_balance(&scaled);
            let strict = FiniteViscoMotionControls {
                maximum_numerical_kinetic_loss_j: 0.5 * first.inertia.numerical_loss_j,
                ..motion_controls()
            };
            match thin.advance(0.04, &supports, &tractions, &[], newton(), strict, cx) {
                Err(FiniteViscoPlaneError::EnergyGate {
                    quantity,
                    value_j,
                    maximum_abs_j,
                }) => {
                    assert_eq!(quantity, "backward-Euler numerical kinetic loss");
                    close(value_j, first.inertia.numerical_loss_j, 2e-12);
                    assert_eq!(maximum_abs_j, strict.maximum_numerical_kinetic_loss_j);
                }
                other => panic!("expected late kinetic-loss gate, got {other:?}"),
            }
            assert_eq!(motion_snapshot(&thin), before);
            let next = first
                .motion
                .advance(
                    0.04,
                    &supports,
                    &tractions,
                    &[],
                    newton(),
                    motion_controls(),
                    cx,
                )
                .unwrap();
            close(next.motion.solid().time_s(), 0.08, 1e-15);
            assert!(next.material_energy.dissipated_j > 0.0);
            assert_ne!(
                next.motion.solid().displacement_m(),
                first.motion.solid().displacement_m()
            );
            check_motion_balance(&next);
        }
    });
}
