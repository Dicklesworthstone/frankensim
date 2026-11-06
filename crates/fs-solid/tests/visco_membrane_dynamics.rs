//! Native finite-strain membrane histories coupled to original reference inertia.

#[path = "visco_membrane/fixture.rs"]
mod fixture;
use fixture::*;

use fs_solid::membrane::MembraneError;
use fs_solid::membrane::inertia::MembraneMassImage;
use fs_solid::membrane::visco::dynamics::{FiniteViscoMembraneMotion, FiniteViscoMotionControls};
use fs_solid::membrane::visco::{
    FiniteViscoMembrane, FiniteViscoMembraneError, ViscoelasticLawField,
};

fn motion_controls() -> FiniteViscoMotionControls {
    FiniteViscoMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

#[test]
fn free_translation_and_acceleration_retain_original_graded_inventory() {
    for image in [MembraneMassImage::Lumped, MembraneMassImage::Consistent] {
        let sheet = source(true);
        let law = material(&sheet);
        let patch = patch(triangle(), sheet)
            .with_reference_areal_density(vec![[0.7, 1.3, 2.0]], || false)
            .unwrap();
        let membrane =
            FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::Uniform(&law), || false)
                .unwrap();
        let v0 = [0.2, -0.1, 0.3];
        let motion =
            FiniteViscoMembraneMotion::try_new(membrane, &[v0; 3], image, || false).unwrap();
        let acceleration = [1.0, -2.0, 0.5];
        let load = patch
            .reference_mass_apply(image, &[acceleration; 3], || false)
            .unwrap();
        let dt = 0.004;
        let step = motion
            .advance(dt, &load, &[], settings(), motion_controls(), || false)
            .unwrap();
        for node in 0..3 {
            for axis in 0..3 {
                let v1 = v0[axis] + dt * acceleration[axis];
                close(step.motion.velocity_mps()[node][axis], v1, 2e-8);
                close(
                    step.motion.membrane().displacement_m()[node][axis],
                    dt * v1,
                    2e-10,
                );
            }
        }
        let mass = patch.reference_mass_kg();
        let expected_loss = 0.5 * mass * dt * dt * acceleration.iter().map(|a| a * a).sum::<f64>();
        close(step.inertia.reference_mass_kg, mass, 1e-14);
        close(step.inertia.numerical_loss_j, expected_loss, 1e-12);
        close(step.material_energy.dissipated_j, 0.0, 1e-16);
        close(step.material_energy.stored_after_j, 0.0, 1e-16);
        close(step.inertia.balance_error_j, 0.0, 1e-12);
        close(step.balance_error_j, 0.0, 1e-10);
        let momentum = step.motion.nodal_momentum_kg_m_s(|| false).unwrap();
        for axis in 0..3 {
            close(
                momentum.iter().map(|p| p[axis]).sum(),
                mass * (v0[axis] + dt * acceleration[axis]),
                2e-10,
            );
        }
        assert_eq!(motion.membrane().time_s(), 0.0);
        assert_eq!(motion.velocity_mps(), &[v0; 3]);
    }
}

#[test]
fn held_finite_strain_relaxes_material_without_numerical_kinetic_heat() {
    for exact in [false, true] {
        let sheet = source(exact);
        let law = material(&sheet);
        let patch = patch(triangle(), sheet);
        let u = affine(&patch, [[0.2, 0.08], [0.02, 0.1], [0.12, -0.04]]);
        let membrane = FiniteViscoMembrane::try_new_from(
            &patch,
            ViscoelasticLawField::Uniform(&law),
            &u,
            || false,
        )
        .unwrap();
        let stored = membrane.stored_energy_j();
        let motion = FiniteViscoMembraneMotion::try_new(
            membrane,
            &[[0.0; 3]; 3],
            MembraneMassImage::Consistent,
            || false,
        )
        .unwrap();
        let mut limits = motion_controls();
        limits.maximum_numerical_kinetic_loss_j = 0.0;
        let step = motion
            .advance(
                0.12,
                &[[0.0; 3]; 3],
                &prescribed(&u),
                settings(),
                limits,
                || false,
            )
            .unwrap();
        assert!(step.material_energy.dissipated_j > 0.0);
        assert!(step.material_energy.stored_after_j < stored);
        assert_eq!(step.inertia.numerical_loss_j, 0.0);
        assert_eq!(step.inertia.kinetic_after_j, 0.0);
        assert_eq!(step.motion.velocity_mps(), &[[0.0; 3]; 3]);
        close(step.material_energy.material_balance_error_j, 0.0, 1e-9);
        assert_eq!(step.motion.membrane().time_s(), 0.12);
    }
}

#[test]
fn moving_support_reactions_include_consistent_inertial_cross_terms() {
    let sheet = source(true);
    let law = material(&sheet);
    let patch = patch(triangle(), sheet)
        .with_reference_areal_density(vec![[0.7, 1.3, 2.0]], || false)
        .unwrap();
    let membrane =
        FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::Uniform(&law), || false)
            .unwrap();
    let endpoint = [[0.001, 0.0005, 0.0003], [0.0; 3], [0.0; 3]];
    let fixed = prescribed(&endpoint);
    let dt = 0.02;
    let quasistatic = membrane
        .advance(dt, &[[0.0; 3]; 3], &fixed, settings(), controls(), || false)
        .unwrap();
    let motion = FiniteViscoMembraneMotion::try_new(
        membrane,
        &[[0.0; 3]; 3],
        MembraneMassImage::Consistent,
        || false,
    )
    .unwrap();
    let step = motion
        .advance(
            dt,
            &[[0.0; 3]; 3],
            &fixed,
            settings(),
            motion_controls(),
            || false,
        )
        .unwrap();
    let acceleration = endpoint.map(|u| u.map(|x| x / (dt * dt)));
    let expected_inertia = patch
        .reference_mass_apply(MembraneMassImage::Consistent, &acceleration, || false)
        .unwrap();
    assert!(expected_inertia[1][0] > 0.0 && expected_inertia[2][0] > 0.0);
    for node in 0..3 {
        for axis in 0..3 {
            close(
                step.equilibrium.residual_n[node][axis]
                    - quasistatic.equilibrium.residual_n[node][axis],
                expected_inertia[node][axis],
                2e-10,
            );
            close(
                step.motion.velocity_mps()[node][axis],
                endpoint[node][axis] / dt,
                1e-13,
            );
        }
    }
    close(
        step.material_energy.dissipated_j,
        quasistatic.energy.dissipated_j,
        1e-12,
    );
    close(
        step.inertia.inertial_work_j,
        step.inertia.kinetic_change_j + step.inertia.numerical_loss_j,
        1e-12,
    );
    close(
        step.material_energy.inertial_endpoint_work_j,
        step.inertia.inertial_work_j,
        1e-12,
    );
}

#[test]
fn rejected_kinetic_budget_and_cancellation_publish_no_motion_or_memory() {
    let sheet = source(false);
    let law = material(&sheet);
    let patch = patch(triangle(), sheet);
    let membrane =
        FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::Uniform(&law), || false)
            .unwrap();
    let accepted = membrane.triangle_states().to_vec();
    let motion = FiniteViscoMembraneMotion::try_new(
        membrane,
        &[[0.0; 3]; 3],
        MembraneMassImage::Consistent,
        || false,
    )
    .unwrap();
    let load = patch
        .reference_mass_apply(MembraneMassImage::Consistent, &[[1.0, 0.0, 0.0]; 3], || {
            false
        })
        .unwrap();
    let mut limits = motion_controls();
    limits.maximum_numerical_kinetic_loss_j = 0.0;
    assert!(matches!(
        motion.advance(0.02, &load, &[], settings(), limits, || false),
        Err(FiniteViscoMembraneError::EnergyGate { .. })
    ));
    let mut calls = 0;
    assert!(matches!(
        motion.advance(0.02, &load, &[], settings(), motion_controls(), || {
            calls += 1;
            calls > 4
        }),
        Err(FiniteViscoMembraneError::Native(MembraneError::Cancelled))
    ));
    assert_eq!(motion.membrane().time_s(), 0.0);
    assert_eq!(motion.membrane().triangle_states(), accepted);
    assert_eq!(motion.velocity_mps(), &[[0.0; 3]; 3]);
    assert_eq!(motion.membrane().displacement_m(), &[[0.0; 3]; 3]);
}

#[test]
fn acceleration_trajectory_has_the_declared_first_order_time_limit() {
    let sheet = source(true);
    let law = material(&sheet);
    let patch = patch(triangle(), sheet);
    let acceleration = [0.0, 0.0, 1.0];
    let load = patch
        .reference_mass_apply(MembraneMassImage::Consistent, &[acceleration; 3], || false)
        .unwrap();
    let duration = 0.02;
    let mut errors = Vec::new();
    let mut numerical_losses = Vec::new();
    for intervals in [2, 4, 8] {
        let membrane =
            FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::Uniform(&law), || false)
                .unwrap();
        let mut motion = FiniteViscoMembraneMotion::try_new(
            membrane,
            &[[0.0; 3]; 3],
            MembraneMassImage::Consistent,
            || false,
        )
        .unwrap();
        let dt = duration / f64::from(intervals);
        let mut numerical_loss = 0.0;
        for _ in 0..intervals {
            let step = motion
                .advance(dt, &load, &[], settings(), motion_controls(), || false)
                .unwrap();
            numerical_loss += step.inertia.numerical_loss_j;
            close(step.material_energy.dissipated_j, 0.0, 1e-16);
            motion = step.motion;
        }
        close(motion.velocity_mps()[0][2], duration, 1e-9);
        errors.push((motion.membrane().displacement_m()[0][2] - 0.5 * duration * duration).abs());
        numerical_losses.push(numerical_loss);
    }
    for i in 0..2 {
        close(errors[i] / errors[i + 1], 2.0, 2e-5);
        close(numerical_losses[i] / numerical_losses[i + 1], 2.0, 2e-5);
    }
}
