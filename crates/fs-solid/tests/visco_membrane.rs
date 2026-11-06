//! G0/G1/G3/G4: native membrane spatial solves retain finite-strain memory.

#[path = "visco_membrane/fixture.rs"]
mod fixture;
use fixture::*;

use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::visco_hyper::membrane::ViscoMembraneState;
use fs_solid::membrane::visco::{
    FiniteViscoMembrane, FiniteViscoMembraneEnergy, FiniteViscoMembraneError, ViscoelasticLawField,
};

fn check_energy(energy: &FiniteViscoMembraneEnergy) {
    assert!(energy.dissipated_j >= 0.0);
    close(energy.material_balance_error_j, 0.0, 1e-9);
    close(energy.equilibrium_work_residual_j, 0.0, 1e-9);
    close(
        energy.native_loading_work_j,
        energy.stored_after_j - energy.stored_before_j + energy.dissipated_j,
        1e-9,
    );
    close(
        energy.endpoint_work_defect_j,
        energy.internal_endpoint_work_j - energy.native_loading_work_j,
        1e-9,
    );
    assert_eq!(energy.inertial_endpoint_work_j, 0.0);
}

#[test]
fn native_material_point_and_spatial_membrane_agree_through_loading_and_hold() {
    for exact in [false, true] {
        let sheet = source(exact);
        let law = material(&sheet);
        let patch = patch(triangle(), sheet);
        let original_mass = patch.reference_mass_kg();
        let membrane =
            FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::Uniform(&law), || false)
                .unwrap();
        let gradient = [[0.2, 0.07], [0.0, 0.1], [0.15, -0.03]];
        let u = affine(&patch, gradient);
        let fixed = prescribed(&u);
        let zero = vec![[0.0; 3]; u.len()];
        let load = membrane
            .advance(0.0, &zero, &fixed, settings(), controls(), || false)
            .unwrap();
        check_energy(&load.energy);
        assert!(load.membrane.stored_energy_j() > 0.0);
        assert_eq!(load.energy.dissipated_j, 0.0);
        assert_eq!(membrane.time_s(), 0.0);
        assert_eq!(membrane.stored_energy_j(), 0.0);

        // Independently form the full surface metric and its canonical upper
        // triangular F. The native 3-D plane-stress owner supplies the reference.
        let f1 = [1.0 + gradient[0][0], gradient[1][0], gradient[2][0]];
        let f2 = [gradient[0][1], 1.0 + gradient[1][1], gradient[2][1]];
        let a = f1.iter().map(|x| x * x).sum::<f64>().sqrt();
        let b = f1.iter().zip(f2).map(|(x, y)| x * y).sum::<f64>() / a;
        let d = (f2.iter().map(|x| x * x).sum::<f64>() - b * b).sqrt();
        let in_plane = [a, b, 0.0, d];
        let dt = 0.12;
        let (expected, storage, loss) = match &load.membrane.triangle_states()[0] {
            ViscoMembraneState::Full(state) => {
                let response = law
                    .step_plane_stress_response(state, &in_plane, dt, None, local())
                    .unwrap();
                (
                    ViscoMembraneState::Full(response.material.state.clone()),
                    response.material.state.stored,
                    response.material.dissipated_pa,
                )
            }
            ViscoMembraneState::Incompressible(state) => {
                let response = law
                    .step_incompressible_plane_stress_response(state, &in_plane, dt, None, local())
                    .unwrap();
                (
                    ViscoMembraneState::Incompressible(response.material.state.clone()),
                    response.material.state.stored,
                    response.material.dissipated_pa,
                )
            }
        };
        let hold = load
            .membrane
            .advance(dt, &zero, &fixed, settings(), controls(), || false)
            .unwrap();
        check_energy(&hold.energy);
        let volume = patch.reference_triangle_areas_m2().next().unwrap() * THICKNESS_M;
        close(hold.membrane.stored_energy_j(), storage * volume, 1e-10);
        close(hold.energy.dissipated_j, loss * volume, 1e-10);
        for (actual, expected) in hold.membrane.triangle_states()[0]
            .deformation_gradient()
            .iter()
            .zip(expected.deformation_gradient())
        {
            close(*actual, *expected, 2e-10);
        }
        assert!(hold.energy.dissipated_j > 0.0);
        assert!(hold.membrane.stored_energy_j() < load.membrane.stored_energy_j());
        assert_eq!(hold.membrane.time_s(), dt);
        assert_eq!(hold.membrane.patch().reference_mass_kg(), original_mass);
        if exact {
            let f = hold.membrane.triangle_states()[0].deformation_gradient();
            close(f[0] * f[4] * f[8], 1.0, 1e-12);
            assert!(hold.membrane.triangle_pressure_pa()[0].unwrap().is_finite());
        }
    }
}

#[test]
fn rigid_spatial_rotation_preserves_material_memory_and_rotates_reactions() {
    let sheet = source(true);
    let law = material(&sheet);
    let patch = patch(triangle(), sheet);
    let u = affine(&patch, [[0.2, 0.12], [0.04, 0.1], [0.15, -0.03]]);
    let initial =
        FiniteViscoMembrane::try_new_from(&patch, ViscoelasticLawField::Uniform(&law), &u, || {
            false
        })
        .unwrap();
    let zero = vec![[0.0; 3]; u.len()];
    let relaxed = initial
        .advance(0.15, &zero, &prescribed(&u), settings(), controls(), || {
            false
        })
        .unwrap();
    // Q(x,y,z)=(z,x,y) is a proper rigid rotation. Reference coordinates remain
    // fixed; the same material metrics must survive the changed spatial frame.
    let rotated: Vec<_> = patch
        .mesh()
        .nodes
        .iter()
        .zip(&u)
        .map(|(x, u)| {
            let y = [x[0] + u[0], x[1] + u[1], u[2]];
            [y[2] - x[0], y[0] - x[1], y[1]]
        })
        .collect();
    let turn = relaxed
        .membrane
        .advance(
            0.0,
            &zero,
            &prescribed(&rotated),
            settings(),
            controls(),
            || false,
        )
        .unwrap();
    check_energy(&turn.energy);
    close(
        turn.membrane.stored_energy_j(),
        relaxed.membrane.stored_energy_j(),
        1e-11,
    );
    assert_eq!(turn.energy.dissipated_j, 0.0);
    assert_eq!(turn.membrane.time_s(), relaxed.membrane.time_s());
    for (before, after) in relaxed
        .equilibrium
        .residual_n
        .iter()
        .zip(&turn.equilibrium.residual_n)
    {
        for (actual, expected) in after.iter().zip([before[2], before[0], before[1]]) {
            close(*actual, expected, 2e-9);
        }
    }
    let (ViscoMembraneState::Incompressible(old), ViscoMembraneState::Incompressible(new)) = (
        &relaxed.membrane.triangle_states()[0],
        &turn.membrane.triangle_states()[0],
    ) else {
        panic!("exact histories required")
    };
    for (a, b) in old
        .inverse_viscous_cauchy_green
        .iter()
        .flatten()
        .zip(new.inverse_viscous_cauchy_green.iter().flatten())
    {
        close(*a, *b, 1e-12);
    }
}

#[test]
fn pressure_loaded_membrane_creeps_with_free_center_and_original_mass() {
    let sheet = source(true);
    let law = material(&sheet);
    let patch = patch(diamond(), sheet);
    let u = affine(&patch, [[0.12, 0.0], [0.0, 0.12], [0.0, 0.0]]);
    let initial =
        FiniteViscoMembrane::try_new_from(&patch, ViscoelasticLawField::Uniform(&law), &u, || {
            false
        })
        .unwrap();
    let fixed: Vec<_> = prescribed(&u)
        .into_iter()
        .filter(|(dof, _)| *dof != 2)
        .collect();
    let zero = vec![[0.0; 3]; u.len()];
    let pressure = 0.5;
    let loaded = initial
        .advance_with_pressure(0.0, &zero, pressure, &fixed, settings(), controls(), || {
            false
        })
        .unwrap();
    let held = loaded
        .membrane
        .advance_with_pressure(0.2, &zero, pressure, &fixed, settings(), controls(), || {
            false
        })
        .unwrap();
    check_energy(&loaded.energy);
    check_energy(&held.energy);
    assert!(loaded.membrane.displacement_m()[0][2] > 0.0);
    assert!(held.membrane.displacement_m()[0][2] > loaded.membrane.displacement_m()[0][2]);
    assert!(held.energy.dissipated_j > 0.0);
    assert!(held.energy.pressure_path_work_j > 0.0);
    close(held.equilibrium.residual_n[0][2], 0.0, 1e-10);
    assert_eq!(
        held.membrane.patch().reference_mass_kg(),
        initial.patch().reference_mass_kg()
    );
    assert_eq!(
        held.membrane.patch().nodal_mass_kg(),
        initial.patch().nodal_mass_kg()
    );
}

#[test]
fn regional_relaxation_times_remain_distinct() {
    let sheet = source(true);
    let law = material(&sheet);
    let mut laws = vec![law; 4];
    for (law, tau) in laws.iter_mut().zip([0.02, 0.2, 2.0, 20.0]) {
        law.branches[0].tau_s = tau;
    }
    let patch = patch(diamond(), sheet);
    let u = affine(&patch, [[0.2, 0.0], [0.0, 0.1], [0.0, 0.0]]);
    let initial = FiniteViscoMembrane::try_new_from(
        &patch,
        ViscoelasticLawField::PerElement(&laws),
        &u,
        || false,
    )
    .unwrap();
    let held = initial
        .advance(
            0.1,
            &vec![[0.0; 3]; u.len()],
            &prescribed(&u),
            settings(),
            controls(),
            || false,
        )
        .unwrap();
    check_energy(&held.energy);
    let stresses: Vec<_> = held
        .membrane
        .triangle_resultants_n_m()
        .iter()
        .map(|s| s[0])
        .collect();
    assert!(
        stresses.windows(2).all(|pair| pair[0] < pair[1]),
        "{stresses:?}"
    );
    assert!(
        FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::PerElement(&laws[..3]), || {
            false
        })
        .is_err()
    );
}

#[test]
fn prescribed_natural_reference_sets_elastic_volume_without_changing_inventory() {
    let g = [1.04, 0.98, 1.02];
    let reference =
        ReferenceDeformation::try_new([g[0], 0.0, 0.0, 0.0, g[1], 0.0, 0.0, 0.0, g[2]]).unwrap();
    for exact in [false, true] {
        let sheet = source(exact)
            .with_hyperelastic_reference_stretches(g)
            .unwrap();
        let law = material(&sheet).with_reference_deformation(reference);
        let patch = patch(triangle(), sheet);
        let u = affine(&patch, [[g[0] - 1.0, 0.0], [0.0, g[1] - 1.0], [0.0, 0.0]]);
        let installed = FiniteViscoMembrane::try_new_from(
            &patch,
            ViscoelasticLawField::Uniform(&law),
            &u,
            || false,
        )
        .unwrap();
        close(installed.stored_energy_j(), 0.0, 1e-12);
        for resultant in installed.triangle_resultants_n_m()[0] {
            close(resultant, 0.0, 1e-11);
        }
        let f = installed.triangle_states()[0].deformation_gradient();
        close(f[0] * f[4] * f[8], g[0] * g[1] * g[2], 1e-12);
        close(patch.reference_mass_kg(), 0.005 * MASS_KG_M2, 1e-14);
        assert_eq!(patch.sheets()[0].thickness_m(), THICKNESS_M);
    }
}

#[test]
fn source_mismatch_work_refusal_and_cancellation_never_replace_accepted_history() {
    let sheet = source(false);
    let law = material(&sheet);
    let patch = patch(triangle(), sheet);
    let u = affine(&patch, [[0.2, 0.0], [0.0, 0.1], [0.0, 0.0]]);
    let initial =
        FiniteViscoMembrane::try_new_from(&patch, ViscoelasticLawField::Uniform(&law), &u, || {
            false
        })
        .unwrap();
    let snapshot = initial.triangle_states().to_vec();
    let zero = vec![[0.0; 3]; u.len()];
    let mut wrong = law.clone();
    wrong.equilibrium.stretch_limit *= 0.5;
    assert!(
        FiniteViscoMembrane::try_new(&patch, ViscoelasticLawField::Uniform(&wrong), || false)
            .is_err()
    );
    let target = affine(&patch, [[0.4, 0.1], [0.0, 0.2], [0.1, 0.0]]);
    let mut no_defect = controls();
    no_defect.maximum_endpoint_work_defect_j = 0.0;
    assert!(matches!(
        initial.advance(
            0.1,
            &zero,
            &prescribed(&target),
            settings(),
            no_defect,
            || false
        ),
        Err(FiniteViscoMembraneError::EnergyGate { .. })
    ));
    let mut polls = 0;
    assert!(
        initial
            .advance(
                0.1,
                &zero,
                &prescribed(&target),
                settings(),
                controls(),
                || {
                    polls += 1;
                    polls >= 8
                }
            )
            .is_err()
    );
    assert!(
        initial
            .advance(
                -0.1,
                &zero,
                &prescribed(&target),
                settings(),
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(initial.triangle_states(), snapshot);
    assert_eq!(initial.displacement_m(), u);
    assert_eq!(initial.time_s(), 0.0);
}
