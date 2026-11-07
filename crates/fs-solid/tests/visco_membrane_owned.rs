//! G0/G3/G4: owned native history and kinetic-only events retain physical state.

#[path = "visco_membrane/fixture.rs"]
mod fixture;
use fixture::*;

use fs_material::hyper::HyperelasticModel;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::sheet::ResolvedSheet;
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::membrane::MembraneError;
use fs_solid::membrane::inertia::MembraneMassImage;
use fs_solid::membrane::visco::dynamics::{FiniteViscoMembraneMotion, FiniteViscoMotionControls};
use fs_solid::membrane::visco::{
    FiniteViscoMembrane, FiniteViscoMembraneError, OwnedFiniteViscoMembraneMotion,
    ViscoelasticLawField,
};

fn motion_controls() -> FiniteViscoMotionControls {
    FiniteViscoMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn two_branch_law(sheet: &ResolvedSheet) -> ViscoHyperelastic {
    let original = material(sheet);
    ViscoHyperelastic::new(
        original.equilibrium,
        vec![
            original.branches[0].clone(),
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3.0 },
                tau_s: 1.5,
            },
        ],
    )
    .unwrap()
}

fn assert_same_motion(a: &FiniteViscoMembraneMotion<'_>, b: &FiniteViscoMembraneMotion<'_>) {
    assert_eq!(a.mass_image(), b.mass_image());
    assert_eq!(a.velocity_mps(), b.velocity_mps());
    assert_eq!(
        a.nodal_momentum_kg_m_s(|| false).unwrap(),
        b.nodal_momentum_kg_m_s(|| false).unwrap()
    );
    let (a, b) = (a.membrane(), b.membrane());
    assert_eq!(a.displacement_m(), b.displacement_m());
    assert_eq!(a.triangle_states(), b.triangle_states());
    assert_eq!(a.triangle_first_piola_pa(), b.triangle_first_piola_pa());
    assert_eq!(a.triangle_pressure_pa(), b.triangle_pressure_pa());
    assert_eq!(a.triangle_resultants_n_m(), b.triangle_resultants_n_m());
    assert_eq!(a.stored_energy_j().to_bits(), b.stored_energy_j().to_bits());
    assert_eq!(a.time_s().to_bits(), b.time_s().to_bits());
}

fn assert_shared_buffers(a: &FiniteViscoMembraneMotion<'_>, b: &FiniteViscoMembraneMotion<'_>) {
    assert_same_motion(a, b);
    assert_eq!(a.velocity_mps().as_ptr(), b.velocity_mps().as_ptr());
    let (a, b) = (a.membrane(), b.membrane());
    assert_eq!(a.displacement_m().as_ptr(), b.displacement_m().as_ptr());
    assert_eq!(a.triangle_states().as_ptr(), b.triangle_states().as_ptr());
    assert_eq!(
        a.triangle_first_piola_pa().as_ptr(),
        b.triangle_first_piola_pa().as_ptr()
    );
    assert_eq!(
        a.triangle_pressure_pa().as_ptr(),
        b.triangle_pressure_pa().as_ptr()
    );
    assert_eq!(
        a.triangle_resultants_n_m().as_ptr(),
        b.triangle_resultants_n_m().as_ptr()
    );
}

#[test]
fn owned_motion_outlives_sources_and_reproduces_the_native_pressure_interval() {
    for exact in [false, true] {
        for image in [MembraneMassImage::Lumped, MembraneMassImage::Consistent] {
            let (owned, expected, equilibrium, material_energy, heat, inertia, balance, u) = {
                let sheet = source(exact);
                let law = two_branch_law(&sheet);
                let patch = patch(triangle(), sheet)
                    .with_reference_areal_density(vec![[0.7, 1.0, 1.4]], || false)
                    .unwrap();
                let u = affine(&patch, [[0.2, 0.06], [0.02, 0.1], [0.12, -0.03]]);
                let initial = FiniteViscoMembrane::try_new_from(
                    &patch,
                    ViscoelasticLawField::Uniform(&law),
                    &u,
                    || false,
                )
                .unwrap();
                let held = initial
                    .advance(
                        0.03,
                        &[[0.0; 3]; 3],
                        &prescribed(&u),
                        settings(),
                        controls(),
                        || false,
                    )
                    .unwrap();
                assert!(held.energy.dissipated_j > 0.0);
                let velocity = [[0.1, -0.2, 0.05], [-0.05, 0.02, 0.03], [0.04, 0.1, -0.02]];
                let motion =
                    FiniteViscoMembraneMotion::try_new(held.membrane, &velocity, image, || false)
                        .unwrap();
                let owned = OwnedFiniteViscoMembraneMotion::from_native(&motion);
                assert_shared_buffers(&motion, &owned.as_native());
                assert_shared_buffers(&owned.as_native(), &owned.clone().as_native());
                assert!(matches!(
                    owned.native_laws(),
                    ViscoelasticLawField::Uniform(_)
                ));
                assert_eq!(
                    owned.patch().reference_corner_areal_density_kg_m2(),
                    patch.reference_corner_areal_density_kg_m2()
                );
                let step = motion
                    .advance_with_pressure(
                        0.08,
                        &[[0.0; 3]; 3],
                        0.4,
                        &prescribed(&u),
                        settings(),
                        motion_controls(),
                        || false,
                    )
                    .unwrap();
                (
                    owned,
                    OwnedFiniteViscoMembraneMotion::from_native(&step.motion),
                    step.equilibrium,
                    step.material_energy,
                    step.element_dissipated_j,
                    step.inertia,
                    step.balance_error_j,
                    u,
                )
            };
            // Neither the original patch nor either original law field remains
            // alive here. Ownership must preserve the real finite-strain history.
            let step = owned
                .advance_with_pressure(
                    0.08,
                    &[[0.0; 3]; 3],
                    0.4,
                    &prescribed(&u),
                    settings(),
                    motion_controls(),
                    || false,
                )
                .unwrap();
            assert_same_motion(&step.motion.as_native(), &expected.as_native());
            assert_eq!(step.equilibrium.displacement_m, equilibrium.displacement_m);
            assert_eq!(step.equilibrium.residual_n, equilibrium.residual_n);
            assert_eq!(step.material_energy, material_energy);
            assert_eq!(step.element_dissipated_j, heat);
            assert_eq!(step.inertia, inertia);
            assert_eq!(step.balance_error_j.to_bits(), balance.to_bits());
            assert!(step.material_energy.dissipated_j > 0.0);
            assert!(step.inertia.numerical_loss_j > 0.0);
            assert!(std::ptr::eq(owned.patch(), step.motion.patch()));
            assert_eq!(owned.as_native().membrane().time_s(), 0.03);
        }
    }
}

#[test]
fn regional_kinetics_retain_every_accepted_buffer_and_change_subsequent_relaxation() {
    for exact in [false, true] {
        let sheet = source(exact);
        let law = two_branch_law(&sheet);
        let patch = patch(diamond(), sheet);
        let u = affine(&patch, [[0.22, 0.06], [0.03, 0.12], [0.1, -0.02]]);
        let zero = vec![[0.0; 3]; u.len()];
        let initial = FiniteViscoMembrane::try_new_from(
            &patch,
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
        assert!(held.energy.dissipated_j > 0.0);
        let motion = FiniteViscoMembraneMotion::try_new(
            held.membrane,
            &zero,
            MembraneMassImage::Consistent,
            || false,
        )
        .unwrap();
        let owned = OwnedFiniteViscoMembraneMotion::from_native(&motion);
        let shifted = {
            let mut targets = vec![law.clone(); 4];
            for (target, factor) in targets.iter_mut().zip([0.3, 0.7, 1.5, 3.0]) {
                for branch in &mut target.branches {
                    branch.tau_s *= factor;
                }
            }
            owned.with_kinetics(&targets, || false).unwrap()
        };
        assert!(std::ptr::eq(owned.patch(), shifted.patch()));
        assert_shared_buffers(&owned.as_native(), &shifted.as_native());
        let ViscoelasticLawField::PerElement(targets) = shifted.native_laws() else {
            panic!("regional target spectra must be retained")
        };
        for (target, factor) in targets.iter().zip([0.3, 0.7, 1.5, 3.0]) {
            for (branch, original) in target.branches.iter().zip(&law.branches) {
                assert_eq!(branch.tau_s, original.tau_s * factor);
            }
        }
        let baseline = owned
            .advance_with_pressure(
                0.08,
                &zero,
                0.0,
                &prescribed(&u),
                settings(),
                motion_controls(),
                || false,
            )
            .unwrap();
        let continued = shifted
            .advance_with_pressure(
                0.08,
                &zero,
                0.0,
                &prescribed(&u),
                settings(),
                motion_controls(),
                || false,
            )
            .unwrap();
        let baseline_native = baseline.motion.as_native();
        let continued_native = continued.motion.as_native();
        let old = baseline_native.membrane().triangle_states();
        let new = continued_native.membrane().triangle_states();
        assert!(new[0].stored_pa() < old[0].stored_pa());
        assert!(new[3].stored_pa() > old[3].stored_pa());
        assert!(
            continued
                .element_dissipated_j
                .iter()
                .all(|heat| *heat > 0.0)
        );
        assert_eq!(continued.inertia.numerical_loss_j, 0.0);
        close(continued_native.membrane().time_s(), 0.18, 1e-15);
        assert_shared_buffers(&motion, &owned.as_native());
    }
}

#[test]
fn nonkinetic_changes_and_late_refusals_preserve_the_complete_owned_state() {
    let sheet = source(true);
    let law = two_branch_law(&sheet);
    let patch = patch(diamond(), sheet);
    let u = affine(&patch, [[0.2, 0.06], [0.02, 0.1], [0.1, -0.03]]);
    let membrane =
        FiniteViscoMembrane::try_new_from(&patch, ViscoelasticLawField::Uniform(&law), &u, || {
            false
        })
        .unwrap();
    let motion = FiniteViscoMembraneMotion::try_new(
        membrane,
        &[[0.1, -0.05, 0.03]; 5],
        MembraneMassImage::Consistent,
        || false,
    )
    .unwrap();
    let owned = OwnedFiniteViscoMembraneMotion::from_native(&motion);
    let mut refused = Vec::new();
    for tau in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        let mut target = law.clone();
        target.branches[0].tau_s = tau;
        refused.push(target);
    }
    let mut target = law.clone();
    target.equilibrium.model = HyperelasticModel::NeoHookean {
        mu: 12.0,
        lambda: 200.0,
    };
    refused.push(target);
    let mut target = law.clone();
    target.equilibrium.stretch_limit *= 1.01;
    refused.push(target);
    let mut target = law.clone();
    target.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 13.0 };
    refused.push(target);
    let mut target = law.clone();
    target.branches.swap(0, 1);
    refused.push(target);
    let mut target = law.clone();
    target.branches.pop();
    refused.push(target);
    refused.push(law.clone().with_reference_deformation(
        ReferenceDeformation::try_new([1.01, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap(),
    ));
    refused.push(
        law.clone()
            .with_local_newton(
                law.local_newton_tolerance() * 2.0,
                law.local_newton_max_iterations(),
            )
            .unwrap(),
    );
    refused.push(
        law.clone()
            .with_local_newton(
                law.local_newton_tolerance(),
                law.local_newton_max_iterations() + 1,
            )
            .unwrap(),
    );
    for target in refused {
        let mut targets = vec![law.clone(); 4];
        targets[3] = target;
        assert!(matches!(
            owned.with_kinetics(&targets, || false),
            Err(FiniteViscoMembraneError::Native(MembraneError::Material {
                element: 3,
                ..
            }))
        ));
        assert_shared_buffers(&motion, &owned.as_native());
    }
    assert!(
        owned
            .with_kinetics(&vec![law.clone(); 3], || false)
            .is_err()
    );
    let mut targets = vec![law; 4];
    for target in &mut targets {
        target.branches[0].tau_s *= 0.5;
    }
    let mut polls = 0;
    owned
        .with_kinetics(&targets, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        owned.with_kinetics(&targets, || {
            seen += 1;
            seen == polls
        }),
        Err(FiniteViscoMembraneError::Native(MembraneError::Cancelled))
    ));
    assert_eq!(
        seen, polls,
        "cancellation must reach the final publication check"
    );
    let mut strict = motion_controls();
    strict.maximum_numerical_kinetic_loss_j = 0.0;
    assert!(matches!(
        owned.advance_with_pressure(
            0.04,
            &[[0.0; 3]; 5],
            0.0,
            &prescribed(&u),
            settings(),
            strict,
            || false,
        ),
        Err(FiniteViscoMembraneError::EnergyGate { .. })
    ));
    assert_shared_buffers(&motion, &owned.as_native());
}
