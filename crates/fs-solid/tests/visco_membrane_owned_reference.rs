//! G0/G3/G4: owned reference events preserve native work, memory and inventory.

#[path = "visco_membrane/fixture.rs"]
mod fixture;
use fixture::*;

use std::sync::Arc;

use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::visco_hyper::membrane::ViscoMembraneState;
use fs_material::visco_hyper::{BranchEnergy, ViscoHyperelastic};
use fs_solid::membrane::inertia::MembraneMassImage;
use fs_solid::membrane::visco::dynamics::{FiniteViscoMembraneMotion, FiniteViscoMotionControls};
use fs_solid::membrane::visco::{
    FiniteViscoMembrane, FiniteViscoMembraneError, OwnedFiniteViscoMembraneMotion,
    ViscoelasticLawField,
};
use fs_solid::membrane::{MembraneBudget, MembraneDomain, MembraneError, MembranePatch};

fn density() -> Vec<[f64; 3]> {
    vec![
        [0.7, 1.0, 1.4],
        [1.1, 0.8, 0.9],
        [0.6, 1.3, 1.2],
        [0.9, 0.95, 1.05],
    ]
}

fn motion_controls() -> FiniteViscoMotionControls {
    FiniteViscoMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn reference(law: ViscoHyperelastic, g: [f64; 3]) -> ViscoHyperelastic {
    law.with_reference_deformation(
        ReferenceDeformation::try_new([g[0], 0.0, 0.0, 0.0, g[1], 0.0, 0.0, 0.0, g[2]]).unwrap(),
    )
}

fn target_patch(exact: bool, stretches: &[[f64; 3]; 4]) -> Arc<MembranePatch> {
    Arc::new(
        MembranePatch::try_new(
            diamond(),
            stretches
                .iter()
                .map(|&g| {
                    source(exact)
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
        .unwrap()
        .with_reference_areal_density(density(), || false)
        .unwrap(),
    )
}

fn accepted(exact: bool, image: MembraneMassImage) -> OwnedFiniteViscoMembraneMotion {
    let sheet = source(exact);
    let law = material(&sheet);
    let patch = patch(diamond(), sheet)
        .with_reference_areal_density(density(), || false)
        .unwrap();
    let u1 = affine(&patch, [[0.2, 0.1], [0.03, 0.08], [0.12, -0.02]]);
    let u2 = affine(&patch, [[0.08, -0.07], [0.05, 0.23], [-0.03, 0.1]]);
    let u3 = affine(&patch, [[0.18, 0.04], [-0.02, 0.13], [0.07, 0.03]]);
    let zero = vec![[0.0; 3]; u1.len()];
    let initial =
        FiniteViscoMembrane::try_new_from(&patch, ViscoelasticLawField::Uniform(&law), &u1, || {
            false
        })
        .unwrap();
    let first = initial
        .advance(
            0.07,
            &zero,
            &prescribed(&u2),
            settings(),
            controls(),
            || false,
        )
        .unwrap();
    let second = first
        .membrane
        .advance(
            0.04,
            &zero,
            &prescribed(&u3),
            settings(),
            controls(),
            || false,
        )
        .unwrap();
    assert!(first.energy.dissipated_j > 0.0 && second.energy.dissipated_j > 0.0);
    let velocity = [
        [0.12, -0.03, 0.01],
        [0.02, 0.1, -0.04],
        [-0.08, 0.04, 0.06],
        [0.03, -0.05, 0.02],
        [0.07, 0.01, -0.03],
    ];
    let motion =
        FiniteViscoMembraneMotion::try_new(second.membrane, &velocity, image, || false).unwrap();
    OwnedFiniteViscoMembraneMotion::from_native(&motion)
}

fn metrics(state: &ViscoMembraneState) -> &[[f64; 6]] {
    match state {
        ViscoMembraneState::Full(state) => &state.inverse_viscous_cauchy_green,
        ViscoMembraneState::Incompressible(state) => &state.inverse_viscous_cauchy_green,
    }
}

fn same_motion(a: &FiniteViscoMembraneMotion<'_>, b: &FiniteViscoMembraneMotion<'_>) {
    assert_eq!(a.mass_image(), b.mass_image());
    assert_eq!(a.velocity_mps(), b.velocity_mps());
    let (a, b) = (a.membrane(), b.membrane());
    assert_eq!(a.displacement_m(), b.displacement_m());
    assert_eq!(a.triangle_states(), b.triangle_states());
    assert_eq!(a.triangle_first_piola_pa(), b.triangle_first_piola_pa());
    assert_eq!(a.triangle_pressure_pa(), b.triangle_pressure_pa());
    assert_eq!(a.triangle_resultants_n_m(), b.triangle_resultants_n_m());
    assert_eq!(a.stored_energy_j().to_bits(), b.stored_energy_j().to_bits());
    assert_eq!(a.time_s().to_bits(), b.time_s().to_bits());
}

#[test]
fn owned_regional_reference_events_match_native_work_and_continue_after_target_drop() {
    for exact in [false, true] {
        for image in [MembraneMassImage::Lumped, MembraneMassImage::Consistent] {
            let mut owned = accepted(exact, image);
            for stage in 0..2 {
                let before = owned.as_native();
                let u = before.membrane().displacement_m().to_vec();
                let zero = vec![[0.0; 3]; u.len()];
                let before_states = before.membrane().triangle_states().to_vec();
                let momentum = before.nodal_momentum_kg_m_s(|| false).unwrap();
                let (event, expected, energy, heat, inertia) = {
                    let stretches = if stage == 0 {
                        [
                            [1.04, 0.98, 1.02],
                            [0.98, 1.03, 1.01],
                            [1.02, 1.01, 0.97],
                            [1.01, 0.99, 1.03],
                        ]
                    } else {
                        [[0.99, 1.01, 1.02]; 4]
                    };
                    let target = target_patch(exact, &stretches);
                    let original = material(&source(exact));
                    let laws: Vec<_> = stretches
                        .into_iter()
                        .map(|g| reference(original.clone(), g))
                        .collect();
                    let field = if stage == 0 {
                        ViscoelasticLawField::PerElement(&laws)
                    } else {
                        ViscoelasticLawField::Uniform(&laws[0])
                    };
                    let native = before
                        .change_reference(&target, field, controls(), || false)
                        .unwrap();
                    let event = owned
                        .change_reference(Arc::clone(&target), field, controls(), || false)
                        .unwrap();
                    same_motion(&event.motion.as_native(), &native.motion);
                    assert_eq!(event.energy, native.energy);
                    assert_eq!(event.impulse_kg_m_s, native.impulse_kg_m_s);
                    assert_eq!(
                        event.kinetic_energy_change_j,
                        native.kinetic_energy_change_j
                    );
                    assert_eq!(event.energy.dissipated_j, 0.0);
                    assert!(event.energy.reference_work_j.abs() > 1e-8);
                    assert!(std::ptr::eq(event.motion.patch(), target.as_ref()));
                    let view = event.motion.as_native();
                    assert_eq!(view.velocity_mps().as_ptr(), before.velocity_mps().as_ptr());
                    assert_eq!(
                        view.membrane().displacement_m().as_ptr(),
                        before.membrane().displacement_m().as_ptr()
                    );
                    assert_eq!(view.nodal_momentum_kg_m_s(|| false).unwrap(), momentum);
                    assert_eq!(
                        event.motion.patch().reference_corner_areal_density_kg_m2(),
                        density()
                    );
                    for (old, new) in before_states.iter().zip(view.membrane().triangle_states()) {
                        assert_eq!(metrics(old), metrics(new));
                        assert_eq!(old.dissipated_pa(), new.dissipated_pa());
                    }
                    let step = native
                        .motion
                        .advance_with_pressure(
                            0.03,
                            &zero,
                            0.0,
                            &prescribed(&u),
                            settings(),
                            motion_controls(),
                            || false,
                        )
                        .unwrap();
                    (
                        event,
                        OwnedFiniteViscoMembraneMotion::from_native(&step.motion),
                        step.material_energy,
                        step.element_dissipated_j,
                        step.inertia,
                    )
                };
                // The temporary target Arc and law field are gone. The proposal
                // still owns the exact prepared state for the following solve.
                let continued = event
                    .motion
                    .advance_with_pressure(
                        0.03,
                        &zero,
                        0.0,
                        &prescribed(&u),
                        settings(),
                        motion_controls(),
                        || false,
                    )
                    .unwrap();
                same_motion(&continued.motion.as_native(), &expected.as_native());
                assert_eq!(continued.material_energy, energy);
                assert_eq!(continued.element_dissipated_j, heat);
                assert_eq!(continued.inertia, inertia);
                assert!(continued.material_energy.dissipated_j > 0.0);
                assert_eq!(before.membrane().triangle_states(), before_states);
                drop(before);
                owned = continued.motion;
            }
        }
    }
}

#[test]
fn native_reference_guards_and_final_cancellation_preserve_the_owned_input() {
    let owned = accepted(true, MembraneMassImage::Consistent);
    let before = owned.as_native();
    let g = [1.04, 0.98, 1.03];
    let target = target_patch(true, &[g; 4]);
    let law = reference(material(&source(true)), g);
    let valid = owned
        .change_reference(
            Arc::clone(&target),
            ViscoelasticLawField::Uniform(&law),
            controls(),
            || false,
        )
        .unwrap();
    assert!(valid.energy.absolute_endpoint_work_defect_j > 0.0);
    let mut strict = controls();
    strict.maximum_endpoint_work_defect_j = 0.0;
    assert!(matches!(
        owned.change_reference(
            Arc::clone(&target),
            ViscoelasticLawField::Uniform(&law),
            strict,
            || false,
        ),
        Err(FiniteViscoMembraneError::EnergyGate { .. })
    ));
    let mut wrong_density = density();
    wrong_density[0] = [0.8, 0.9, 1.4]; // Same triangle mass, different P1 inventory.
    let inventory = Arc::new(
        (*target)
            .clone()
            .with_reference_areal_density(wrong_density, || false)
            .unwrap(),
    );
    let wrong_volume = target_patch(false, &[g; 4]);
    let mut mesh = diamond();
    mesh.nodes[1][0] *= 1.01;
    let geometry = Arc::new(
        patch(
            mesh,
            source(true)
                .with_hyperelastic_reference_stretches(g)
                .unwrap(),
        )
        .with_reference_areal_density(density(), || false)
        .unwrap(),
    );
    for refused in [inventory, wrong_volume, geometry] {
        assert!(
            owned
                .change_reference(
                    refused,
                    ViscoelasticLawField::Uniform(&law),
                    controls(),
                    || false
                )
                .is_err()
        );
    }
    let wrong_reference = material(&source(true));
    let mut wrong_branch = law.clone();
    wrong_branch.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 13.0 };
    for refused in [&wrong_reference, &wrong_branch] {
        assert!(
            owned
                .change_reference(
                    Arc::clone(&target),
                    ViscoelasticLawField::Uniform(refused),
                    controls(),
                    || false
                )
                .is_err()
        );
    }
    assert!(
        owned
            .change_reference(
                Arc::clone(&target),
                ViscoelasticLawField::PerElement(&[]),
                controls(),
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    owned
        .change_reference(
            Arc::clone(&target),
            ViscoelasticLawField::Uniform(&law),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        owned.change_reference(
            target,
            ViscoelasticLawField::Uniform(&law),
            controls(),
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(FiniteViscoMembraneError::Native(MembraneError::Cancelled))
    ));
    assert_eq!(seen, polls);
    same_motion(&owned.as_native(), &before);
    assert_eq!(
        owned.as_native().membrane().triangle_states().as_ptr(),
        before.membrane().triangle_states().as_ptr()
    );
}
