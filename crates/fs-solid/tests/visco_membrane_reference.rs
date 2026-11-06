//! Natural-strain events retain the original spatial inventory and native memory.

#[path = "visco_membrane/fixture.rs"]
mod fixture;
use fixture::*;

use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::visco_hyper::membrane::ViscoMembraneState;
use fs_material::visco_hyper::{BranchEnergy, ViscoHyperelastic};
use fs_solid::membrane::inertia::MembraneMassImage;
use fs_solid::membrane::visco::dynamics::{FiniteViscoMembraneMotion, FiniteViscoMotionControls};
use fs_solid::membrane::visco::{
    FiniteViscoMembrane, FiniteViscoMembraneError, FiniteViscoMembraneReferenceEnergy,
    ViscoelasticLawField,
};
use fs_solid::membrane::{MembraneBudget, MembraneDomain, MembranePatch};

fn reference(law: ViscoHyperelastic, g: [f64; 3]) -> ViscoHyperelastic {
    law.with_reference_deformation(
        ReferenceDeformation::try_new([g[0], 0.0, 0.0, 0.0, g[1], 0.0, 0.0, 0.0, g[2]]).unwrap(),
    )
}

fn metrics(state: &ViscoMembraneState) -> &[[f64; 6]] {
    match state {
        ViscoMembraneState::Full(state) => &state.inverse_viscous_cauchy_green,
        ViscoMembraneState::Incompressible(state) => &state.inverse_viscous_cauchy_green,
    }
}

fn check_work(energy: &FiniteViscoMembraneReferenceEnergy) {
    assert_eq!(energy.dissipated_j, 0.0);
    close(energy.material_balance_error_j, 0.0, 1e-12);
    close(
        energy.reference_work_j + energy.native_transverse_loading_work_j,
        energy.stored_after_j - energy.stored_before_j,
        1e-12,
    );
    close(
        energy.reference_work_j
            + energy.constraint_reference_work_j
            + energy.transverse_endpoint_work_j
            - energy.endpoint_work_defect_j,
        energy.stored_after_j - energy.stored_before_j,
        1e-12,
    );
    assert!(energy.absolute_endpoint_work_defect_j >= energy.endpoint_work_defect_j.abs());
}

#[test]
fn reference_event_matches_native_transactions_after_finite_strain_relaxation() {
    for exact in [false, true] {
        let sheet = source(exact);
        let law = material(&sheet);
        let density = vec![[0.7, 1.0, 1.4]];
        let original = patch(triangle(), sheet.clone())
            .with_reference_areal_density(density.clone(), || false)
            .unwrap();
        let u = affine(&original, [[0.22, 0.08], [0.03, 0.11], [0.14, -0.02]]);
        let zero = vec![[0.0; 3]; u.len()];
        let initial = FiniteViscoMembrane::try_new_from(
            &original,
            ViscoelasticLawField::Uniform(&law),
            &u,
            || false,
        )
        .unwrap();
        let held = initial
            .advance(0.12, &zero, &prescribed(&u), settings(), controls(), || {
                false
            })
            .unwrap();
        assert!(held.energy.dissipated_j > 0.0);
        let before = held.membrane.triangle_states()[0].clone();
        let g = [1.05, 0.97, 1.03];
        let mut target = reference(law.clone(), g);
        target.branches[0].tau_s *= 0.5;
        let target_patch = patch(
            triangle(),
            sheet.with_hyperelastic_reference_stretches(g).unwrap(),
        )
        .with_reference_areal_density(density, || false)
        .unwrap();
        let event = held
            .membrane
            .change_reference(
                &target_patch,
                ViscoelasticLawField::Uniform(&target),
                controls(),
                || false,
            )
            .unwrap();

        // Independent native three-dimensional reference transactions own the
        // physical point, pressure actuation and work; original area/thickness
        // supply the spatial integration, without a reconstructed virgin state.
        let old_f = before.deformation_gradient();
        let (expected, piola, reference_work, constraint_work, loading_work, defect) = match &before
        {
            ViscoMembraneState::Full(state) => {
                let native = target
                    .change_plane_stress_reference_from(
                        &law,
                        state,
                        Some([old_f[2], old_f[5], old_f[8]]),
                        local(),
                    )
                    .unwrap();
                assert_eq!(native.material.dissipated_pa, 0.0);
                (
                    ViscoMembraneState::Full(native.material.state),
                    native.material.first_piola_pa,
                    native.reference_work_pa,
                    native.constraint_reference_work_pa,
                    native.material.loading_work_pa,
                    native.endpoint_work_defect_pa,
                )
            }
            ViscoMembraneState::Incompressible(state) => {
                let native = target
                    .change_incompressible_plane_stress_reference_from(
                        &law,
                        state,
                        Some([old_f[2], old_f[5]]),
                        local(),
                    )
                    .unwrap();
                assert_eq!(native.material.dissipated_pa, 0.0);
                (
                    ViscoMembraneState::Incompressible(native.material.state),
                    native.plane_stress.full_piola,
                    native.reference_work_pa,
                    native.constraint_reference_work_pa,
                    native.material.loading_work_pa,
                    native.endpoint_work_defect_pa,
                )
            }
        };
        assert_eq!(event.membrane.triangle_states()[0], expected);
        assert_eq!(event.membrane.triangle_first_piola_pa()[0], piola);
        assert_eq!(metrics(&expected), metrics(&before));
        assert_eq!(expected.dissipated_pa(), before.dissipated_pa());
        assert_eq!(
            event.membrane.displacement_m(),
            held.membrane.displacement_m()
        );
        assert_eq!(event.membrane.time_s(), held.membrane.time_s());
        assert_eq!(target_patch.nodal_mass_kg(), original.nodal_mass_kg());
        assert_eq!(
            target_patch.reference_corner_areal_density_kg_m2(),
            original.reference_corner_areal_density_kg_m2()
        );
        let volume = original.reference_triangle_areas_m2().next().unwrap() * THICKNESS_M;
        check_work(&event.energy);
        close(
            event.energy.reference_work_j,
            volume * reference_work,
            1e-14,
        );
        close(
            event.energy.constraint_reference_work_j,
            volume * constraint_work,
            1e-14,
        );
        close(
            event.energy.native_transverse_loading_work_j,
            volume * loading_work,
            1e-14,
        );
        close(event.energy.endpoint_work_defect_j, volume * defect, 1e-14);
        close(
            event.energy.stored_after_j,
            volume * expected.stored_pa(),
            1e-14,
        );
        assert!(event.energy.reference_work_j.abs() > 1e-8);
        if exact {
            let f = expected.deformation_gradient();
            close(f[0] * f[4] * f[8], g.iter().product(), 1e-12);
            assert!(event.energy.constraint_reference_work_j.abs() > 1e-8);
        }
        let continued = event
            .membrane
            .advance(0.08, &zero, &prescribed(&u), settings(), controls(), || {
                false
            })
            .unwrap();
        close(continued.membrane.time_s(), 0.2, 1e-15);
        assert!(continued.energy.dissipated_j > 0.0);
        assert!(continued.membrane.triangle_states()[0].dissipated_pa() > before.dissipated_pa());
        assert_eq!(held.membrane.triangle_states()[0], before);
    }
}

#[test]
fn released_membrane_moves_to_prescribed_natural_shape_without_inventing_heat() {
    for exact in [false, true] {
        let sheet = source(exact);
        let law = material(&sheet);
        let original = patch(triangle(), sheet.clone());
        let initial =
            FiniteViscoMembrane::try_new(&original, ViscoelasticLawField::Uniform(&law), || false)
                .unwrap();
        let g = [0.96, 0.98, 1.02];
        let target = reference(law.clone(), g);
        let target_patch = patch(
            triangle(),
            sheet.with_hyperelastic_reference_stretches(g).unwrap(),
        );
        let event = initial
            .change_reference(
                &target_patch,
                ViscoelasticLawField::Uniform(&target),
                controls(),
                || false,
            )
            .unwrap();
        assert!(event.membrane.stored_energy_j() > 0.0);
        assert_eq!(event.membrane.displacement_m(), initial.displacement_m());
        // Retain planar motion and remove only the three in-plane rigid modes.
        let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let released = event
            .membrane
            .advance(
                0.0,
                &[[0.0; 3]; 3],
                &supports,
                settings(),
                controls(),
                || false,
            )
            .unwrap();
        let expected = affine(
            &target_patch,
            [[g[0] - 1.0, 0.0], [0.0, g[1] - 1.0], [0.0, 0.0]],
        );
        for (actual, expected) in released
            .membrane
            .displacement_m()
            .iter()
            .flatten()
            .zip(expected.iter().flatten())
        {
            close(*actual, *expected, 2e-8);
        }
        close(released.membrane.stored_energy_j(), 0.0, 1e-11);
        assert_eq!(released.energy.dissipated_j, 0.0);
        assert_eq!(released.membrane.time_s(), 0.0);
        assert_eq!(
            released.membrane.patch().nodal_mass_kg(),
            original.nodal_mass_kg()
        );
    }
}

#[test]
fn regional_reference_and_kinetic_changes_preserve_original_point_order() {
    let sheet = source(true);
    let law = material(&sheet);
    let original = patch(diamond(), sheet.clone());
    let u = affine(&original, [[0.18, 0.06], [0.02, 0.09], [0.08, -0.03]]);
    let zero = vec![[0.0; 3]; u.len()];
    let initial = FiniteViscoMembrane::try_new_from(
        &original,
        ViscoelasticLawField::Uniform(&law),
        &u,
        || false,
    )
    .unwrap();
    let held = initial
        .advance(0.1, &zero, &prescribed(&u), settings(), controls(), || {
            false
        })
        .unwrap();
    let stretches = [
        [1.02, 0.98, 1.01],
        [0.97, 1.03, 1.04],
        [1.01, 1.01, 0.97],
        [1.0; 3],
    ];
    let laws: Vec<_> = stretches
        .iter()
        .zip([0.02, 0.2, 2.0, 20.0])
        .map(|(&g, tau)| {
            let mut target = reference(law.clone(), g);
            target.branches[0].tau_s = tau;
            target
        })
        .collect();
    let target_patch = MembranePatch::try_new(
        diamond(),
        stretches
            .iter()
            .map(|&g| {
                sheet
                    .clone()
                    .with_hyperelastic_reference_stretches(g)
                    .unwrap()
            })
            .collect(),
        vec![0.0; 4],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 2.0,
        },
        MembraneBudget {
            max_nodes: 32,
            max_triangles: 32,
        },
        || false,
    )
    .unwrap();
    let event = held
        .membrane
        .change_reference(
            &target_patch,
            ViscoelasticLawField::PerElement(&laws),
            controls(),
            || false,
        )
        .unwrap();
    check_work(&event.energy);
    for (i, g) in stretches.iter().enumerate() {
        let old = &held.membrane.triangle_states()[i];
        let new = &event.membrane.triangle_states()[i];
        assert_eq!(metrics(new), metrics(old));
        assert_eq!(new.dissipated_pa(), old.dissipated_pa());
        let f = new.deformation_gradient();
        close(f[0] * f[4] * f[8], g.iter().product(), 1e-12);
    }
    // The final region changed kinetics alone; the event has no physical work.
    let before = &held.membrane.triangle_states()[3];
    let after = &event.membrane.triangle_states()[3];
    assert_eq!(after.reference_work_pa(), before.reference_work_pa());
    close(after.stored_pa(), before.stored_pa(), 1e-12);
    let continued = event
        .membrane
        .advance(0.05, &zero, &prescribed(&u), settings(), controls(), || {
            false
        })
        .unwrap();
    assert!(continued.energy.dissipated_j > 0.0);
    assert_eq!(continued.membrane.triangle_states().len(), 4);
}

#[test]
fn motion_reference_event_retains_velocity_and_consistent_or_lumped_momentum() {
    let sheet = source(true);
    let law = material(&sheet);
    let density = vec![[0.7, 1.0, 1.4]];
    let original = patch(triangle(), sheet.clone())
        .with_reference_areal_density(density.clone(), || false)
        .unwrap();
    let g = [1.03, 0.98, 1.02];
    let target = reference(law.clone(), g);
    let target_patch = patch(
        triangle(),
        sheet.with_hyperelastic_reference_stretches(g).unwrap(),
    )
    .with_reference_areal_density(density, || false)
    .unwrap();
    let u = affine(&original, [[0.18, 0.02], [0.03, 0.11], [0.04, 0.02]]);
    let velocity = [[0.2, -0.1, 0.05], [-0.1, 0.15, 0.02], [0.04, -0.03, -0.08]];
    for image in [MembraneMassImage::Consistent, MembraneMassImage::Lumped] {
        let initial = FiniteViscoMembrane::try_new_from(
            &original,
            ViscoelasticLawField::Uniform(&law),
            &u,
            || false,
        )
        .unwrap();
        let motion =
            FiniteViscoMembraneMotion::try_new(initial, &velocity, image, || false).unwrap();
        let momentum = motion.nodal_momentum_kg_m_s(|| false).unwrap();
        let event = motion
            .change_reference(
                &target_patch,
                ViscoelasticLawField::Uniform(&target),
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(event.motion.mass_image(), image);
        assert_eq!(event.impulse_kg_m_s, [0.0; 3]);
        assert_eq!(event.kinetic_energy_change_j, 0.0);
        assert_eq!(event.motion.velocity_mps(), velocity);
        assert_eq!(
            event.motion.nodal_momentum_kg_m_s(|| false).unwrap(),
            momentum
        );
        assert_eq!(
            event.motion.membrane().displacement_m(),
            motion.membrane().displacement_m()
        );
        assert_eq!(event.motion.membrane().time_s(), motion.membrane().time_s());
        assert_eq!(
            metrics(&event.motion.membrane().triangle_states()[0]),
            metrics(&motion.membrane().triangle_states()[0])
        );
        check_work(&event.energy);
        let continued = event
            .motion
            .advance(
                0.04,
                &[[0.0; 3]; 3],
                &prescribed(&u),
                settings(),
                FiniteViscoMotionControls {
                    material: controls(),
                    maximum_numerical_kinetic_loss_j: 1.0,
                },
                || false,
            )
            .unwrap();
        assert_eq!(continued.motion.membrane().time_s(), 0.04);
        assert!(continued.material_energy.dissipated_j > 0.0);
        assert!(continued.inertia.numerical_loss_j > 0.0);
    }
}

#[test]
fn reference_refusals_and_late_cancellation_leave_the_accepted_body_unchanged() {
    let sheet = source(true);
    let law = material(&sheet);
    let original = patch(triangle(), sheet.clone());
    let u = affine(&original, [[0.2, 0.06], [0.02, 0.1], [0.05, -0.01]]);
    let initial = FiniteViscoMembrane::try_new_from(
        &original,
        ViscoelasticLawField::Uniform(&law),
        &u,
        || false,
    )
    .unwrap();
    let before = initial.triangle_states().to_vec();
    let g = [1.06, 0.97, 1.04];
    let target = reference(law.clone(), g);
    let target_sheet = sheet.with_hyperelastic_reference_stretches(g).unwrap();
    let target_patch = patch(triangle(), target_sheet.clone());
    let mut polls = 0;
    let proposed = initial
        .change_reference(
            &target_patch,
            ViscoelasticLawField::Uniform(&target),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(proposed.energy.absolute_endpoint_work_defect_j > 0.0);
    let mut strict = controls();
    strict.maximum_endpoint_work_defect_j = 0.0;
    assert!(matches!(
        initial.change_reference(
            &target_patch,
            ViscoelasticLawField::Uniform(&target),
            strict,
            || false,
        ),
        Err(FiniteViscoMembraneError::EnergyGate { .. })
    ));
    let mut count = 0;
    assert!(
        initial
            .change_reference(
                &target_patch,
                ViscoelasticLawField::Uniform(&target),
                controls(),
                || {
                    count += 1;
                    count == polls
                },
            )
            .is_err()
    );
    assert_eq!(
        count, polls,
        "cancellation must reach the final publication poll"
    );

    let mut changed_mesh = triangle();
    changed_mesh.nodes[1][0] *= 1.01;
    let changed_geometry = patch(changed_mesh, target_sheet.clone());
    let changed_inventory = patch(triangle(), target_sheet.clone())
        .with_reference_areal_density(vec![[0.7, 0.8, 0.9]], || false)
        .unwrap();
    let changed_volume = patch(
        triangle(),
        source(false)
            .with_hyperelastic_reference_stretches(g)
            .unwrap(),
    );
    for changed_patch in [&changed_geometry, &changed_inventory, &changed_volume] {
        assert!(
            initial
                .change_reference(
                    changed_patch,
                    ViscoelasticLawField::Uniform(&target),
                    controls(),
                    || false,
                )
                .is_err()
        );
    }
    let mut changed_branch = target.clone();
    changed_branch.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 13.0 };
    for refused in [&law, &changed_branch] {
        assert!(
            initial
                .change_reference(
                    &target_patch,
                    ViscoelasticLawField::Uniform(refused),
                    controls(),
                    || false,
                )
                .is_err()
        );
    }
    assert!(
        initial
            .change_reference(
                &target_patch,
                ViscoelasticLawField::PerElement(&[]),
                controls(),
                || false,
            )
            .is_err()
    );
    assert_eq!(initial.triangle_states(), before);
    assert_eq!(initial.displacement_m(), u);
    assert_eq!(initial.time_s(), 0.0);
}
