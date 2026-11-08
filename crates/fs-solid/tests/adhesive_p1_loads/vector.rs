//! Full-stick adhesive vector fields on the original structural traces.
use super::*;
use fs_contact::normal_patch::adhesive::torsion::{JkrTorsionControls, SourceJkrTorsionContact};

fn mixed_source() -> SourceJkrMixedMode {
    SourceJkrMixedMode::new(
        sources::jkr(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap()
}

// At zero normal force the original Griffith balance gives this dimensional
// radius independently of the pressure/sector projection. On each edge node,
// N=positive(x)/L or positive(y)/L. The full-stick flat-punch field integrates
// to Q*a/(4L); the centre receives Q*(1-a/L).
fn expected_nodes(source: &SourceJkrMixedMode, shear: [f64; 2], l: f64) -> [[f64; 3]; 6] {
    let normal = source.normal_source();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let b = 3.0 * PI * w * r;
    let q2 = shear[0] * shear[0] + shear[1] * shear[1];
    let hertz = b + (b * b - e * q2 / (4.0 * source.reduced_shear_modulus_pa())).sqrt();
    let radius = (3.0 * r * hertz / (4.0 * e)).cbrt();
    let edge = radius / (4.0 * l);
    let normal_edge = radius * hertz / (16.0 * l);
    let mut result = [[0.0; 3]; 6];
    for node in 0..5 {
        let weight = if node == 0 { 1.0 - 4.0 * edge } else { edge };
        result[node] = [
            weight * shear[0],
            -weight * shear[1], // Original local y axis is world -y.
            if node == 0 {
                -4.0 * normal_edge
            } else {
                normal_edge
            },
        ];
    }
    result
}

fn sparse(forces: &[[f64; 3]; 6]) -> Vec<(usize, f64)> {
    forces
        .iter()
        .enumerate()
        .flat_map(|(node, row)| {
            row.iter().enumerate().filter_map(move |(axis, &force)| {
                (force != 0.0).then_some((3 * node + axis, force))
            })
        })
        .collect()
}

fn rotate(v: [f64; 3], rotated: bool) -> [f64; 3] {
    if rotated { [v[2], v[0], v[1]] } else { v }
}

#[test]
fn vector_projection_recovers_independent_nodal_shear_force_and_rotated_wrench() {
    let source = mixed_source();
    let critical = source.critical_tangential_force_n(0.0).unwrap();
    let shear = [0.3 * critical, -0.4 * critical];
    let profile = AdhesiveTractionProfile::jkr_mixed_mode_force(&source, 0.0, shear).unwrap();
    let a = profile.contact_radius_m();
    let l = 4.0 * a;
    for rotated in [false, true] {
        let nodes = nodes(l).map(|v| rotate(v, rotated));
        support::with_cx(|cx| {
            let geometry =
                TetLoadGeometry3::new(&nodes, &TETS, Default::default(), 16, cx).unwrap();
            let other = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let base_frame = frame();
            let selected_frame = AdhesiveP1Frame3 {
                centre_m: [0.0; 3],
                x_axis: rotate(base_frame.x_axis, rotated),
                y_axis: rotate(base_frame.y_axis, rotated),
            };
            let make = || {
                AdhesiveP1Patch3::new(
                    &geometry,
                    selected_frame,
                    grid(&geometry, vec![0.0, 0.37 * a, a]),
                    controls(&profile),
                    &gate(),
                    cx,
                )
                .unwrap()
            };
            let (_, normal_only) = make().assemble(&profile, &other, &gate(), cx).unwrap();
            assert!(normal_only.vector_source.is_none());
            assert!(
                normal_only
                    .force_world_n
                    .into_iter()
                    .all(|v| v.abs() < 1e-14)
            );
            let mut patch = make().with_full_stick_tractions();
            let (load, transfer) = patch.assemble(&profile, &other, &gate(), cx).unwrap();
            let expected = sparse(&expected_nodes(&source, shear, l).map(|v| rotate(v, rotated)));
            assert_eq!(load.nodal_forces_n().len(), expected.len());
            for (&(dof, actual), &(expected_dof, force)) in
                load.nodal_forces_n().iter().zip(&expected)
            {
                assert_eq!(dof, expected_dof);
                close(actual, force, 1e-17);
            }
            let resultant = rotate([shear[0], -shear[1], 0.0], rotated);
            for (actual, expected) in transfer.force_world_n.into_iter().zip(resultant) {
                close(actual, expected, 1e-17);
            }
            assert!(
                transfer
                    .moment_about_centre_n_m
                    .into_iter()
                    .all(|v| v.abs() < 1e-18)
            );
            let vector = transfer.vector_source.unwrap();
            assert_eq!(vector.normal.model_identity, source.identity());
            assert_eq!(vector.normal, transfer.source);
            assert_eq!(vector.normal.radial_evaluations, 0);
            assert!(
                vector
                    .tangential_force_residual_n
                    .into_iter()
                    .all(|v| v.abs() < 1e-14)
            );
        });
    }
}

#[test]
fn changing_vector_adhesion_drives_native_motion_work_and_atomic_retry() {
    let source = mixed_source();
    let critical = source.critical_tangential_force_n(0.0).unwrap();
    let shear = [
        [0.2 * critical, -0.1 * critical],
        [-0.35 * critical, 0.25 * critical],
        [0.1 * critical, 0.4 * critical],
    ];
    let profiles =
        shear.map(|q| AdhesiveTractionProfile::jkr_mixed_mode_force(&source, 0.0, q).unwrap());
    let outer = AdhesiveTractionProfile::jkr_mixed_mode_force(&source, 0.0, [0.0; 2])
        .unwrap()
        .contact_radius_m();
    let l = 4.0 * outer;
    let nodes = nodes(l);
    let material = material();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &[],
        budget: Default::default(),
    };
    support::with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
        let other = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let reference_load = |index| {
            geometry
                .assemble(
                    &[],
                    VolumeLoad3::None,
                    &sparse(&expected_nodes(&source, shear[index], l)),
                    cx,
                )
                .unwrap()
        };
        let patch = AdhesiveP1Patch3::new(
            &geometry,
            frame(),
            grid(&geometry, vec![0.0, 0.4 * outer, outer]),
            controls(&profiles[0]),
            &gate(),
            cx,
        )
        .unwrap()
        .with_full_stick_tractions();
        let mut actual = patch
            .begin_transient(
                &profiles[0],
                &other,
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                native_config(),
                &gate(),
                cx,
            )
            .unwrap();
        let mut reference = reference_load(0)
            .begin_transient(
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                native_config(),
                cx,
            )
            .unwrap();
        for index in [1, 2, 0, 1] {
            let checkpoint = actual.native().state().clone();
            let accepted = *actual.current_transfer();
            let rejected: Result<(), AdhesiveP1Error3> = actual.advance_with_admission(
                &profiles[index],
                &other,
                2e-7,
                &gate(),
                cx,
                |_, _| {
                    Err(AdhesiveP1Error3::Invalid(
                        "receiver rejects complete vector candidate",
                    ))
                },
            );
            assert!(rejected.is_err());
            assert_eq!(actual.native().state(), &checkpoint);
            assert_eq!(*actual.current_transfer(), accepted);
            let step = actual
                .advance(&profiles[index], &other, 2e-7, &gate(), cx)
                .unwrap()
                .clone();
            let expected = reference.advance(&reference_load(index), 2e-7, cx).unwrap();
            for (&a, &b) in actual
                .native()
                .state()
                .displacement
                .iter()
                .zip(&reference.state().displacement)
            {
                close(a, b, 1e-20);
            }
            for (&a, &b) in actual
                .native()
                .state()
                .velocity
                .iter()
                .zip(&reference.state().velocity)
            {
                close(a, b, 1e-14);
            }
            close(
                step.traction_work_on_solid_j,
                expected.mechanical_work_j,
                1e-24,
            );
            close(
                step.spatial_traction_work_j,
                expected.mechanical_work_j,
                1e-24,
            );
            assert!(step.work_residual_j.abs() <= step.work_tolerance_j);
            assert_eq!(step.quadrature_work_estimate_j, 0.0);
            assert!(step.transfer.vector_source.is_some());
        }
        assert!(
            actual
                .native()
                .state()
                .displacement
                .iter()
                .any(|u| *u != 0.0)
        );
        assert!(actual.native().state().strain_energy_j > 0.0);
    });
}

#[test]
fn prescribed_tangential_motion_accounts_for_adhesion_work_on_supported_nodes() {
    let source = mixed_source();
    let critical = source.critical_tangential_force_n(0.0).unwrap();
    let shear = [0.3 * critical, 0.4 * critical];
    let profile = AdhesiveTractionProfile::jkr_mixed_mode_force(&source, 0.0, shear).unwrap();
    let a = profile.contact_radius_m();
    let nodes = nodes(4.0 * a);
    let material = material();
    let fixed: Vec<_> = (0..18).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    support::with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
        let other = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let patch = AdhesiveP1Patch3::new(
            &geometry,
            frame(),
            grid(&geometry, vec![0.0, a]),
            controls(&profile),
            &gate(),
            cx,
        )
        .unwrap()
        .with_full_stick_tractions();
        let mut native = patch
            .begin_transient_with_supports(
                &profile,
                &other,
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                native_config(),
                &gate(),
                cx,
            )
            .unwrap();
        let acceleration: Vec<_> = (0..6)
            .flat_map(|node| [(3 * node, 0.8), (3 * node + 1, -0.6)])
            .collect();
        let h = 2e-7;
        let step = native
            .advance_with_supports(
                &profile,
                &other,
                AdhesiveP1SupportEndpoint3::Acceleration(&acceleration),
                h,
                &gate(),
                cx,
            )
            .unwrap();
        // With zero initial acceleration and velocity, Newmark travel is h²*a/4.
        // All loaded nodes translate alike, so only the original resultant enters.
        let expected = 0.25 * h * h * (0.8 * shear[0] + 0.6 * shear[1]);
        assert!(expected > 0.0);
        close(step.spatial_traction_work_j, expected, 1e-24);
        close(step.traction_work_on_solid_j, expected, 1e-24);
        close(step.native.mechanical_work_j, expected, 1e-24);
        assert_eq!(step.other_external_work_j, 0.0);
        assert!(step.work_residual_j.abs() <= step.work_tolerance_j);
        assert_eq!(step.quadrature_work_estimate_j, 0.0);
    });
}

fn torsion_profile() -> (SourceJkrTorsionContact, AdhesiveTractionProfile) {
    let normal = sources::jkr();
    let indentation = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let source = SourceJkrTorsionContact::new(
        normal,
        indentation,
        JkrTorsionControls {
            critical_twist_angle_rad: 1e-5,
            maximum_critical_twist_angle_rad: 0.01,
            maximum_mode_iii_energy_release_to_adhesion: 0.01,
        },
    )
    .unwrap();
    let profile = AdhesiveTractionProfile::jkr_torsion_elastic(&source, 4e-6).unwrap();
    (source, profile)
}

fn torsion_nodes(profile: &AdhesiveTractionProfile, moment: f64, l: f64) -> [[f64; 3]; 6] {
    let mut result = [[0.0; 3]; 6];
    for (dof, force) in jkr_nodal(profile, l) {
        result[dof / 3][dof % 3] = force;
    }
    // On each original edge N=max(+/-x,0)/L or max(+/-y,0)/L.
    // Integrating q_phi times N gives four forces of magnitude |M|/(4L),
    // with zero centre load and world couple -M*z because local normal=-z.
    result[1][1] = -moment / (4.0 * l);
    result[2][0] = moment / (4.0 * l);
    result[3][1] = moment / (4.0 * l);
    result[4][0] = -moment / (4.0 * l);
    result
}

#[test]
fn elastic_torsion_projects_original_nodal_couple_and_drives_native_deformation() {
    let (source, profile) = torsion_profile();
    let a = profile.contact_radius_m();
    let l = 4.0 * a;
    let moment = -(16.0 / 3.0) * source.law().reduced_shear_modulus_pa() * a.powi(3) * 4e-6;
    for rotated in [false, true] {
        let nodes = nodes(l).map(|v| rotate(v, rotated));
        let material = material();
        let problem = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &TETS,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &[],
            budget: Default::default(),
        };
        support::with_cx(|cx| {
            let geometry = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
            let other = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let base = frame();
            let selected_frame = AdhesiveP1Frame3 {
                centre_m: [0.0; 3],
                x_axis: rotate(base.x_axis, rotated),
                y_axis: rotate(base.y_axis, rotated),
            };
            let make = || {
                AdhesiveP1Patch3::new(
                    &geometry,
                    selected_frame,
                    grid(&geometry, vec![0.0, 0.37 * a, a]),
                    controls(&profile),
                    &gate(),
                    cx,
                )
                .unwrap()
            };
            let (_, normal_only) = make().assemble(&profile, &other, &gate(), cx).unwrap();
            assert!(normal_only.vector_source.is_none());
            for m in normal_only.moment_about_centre_n_m {
                close(m, 0.0, 1e-18);
            }
            let mut patch = make().with_full_stick_tractions();
            let (load, transfer) = patch.assemble(&profile, &other, &gate(), cx).unwrap();
            let expected = torsion_nodes(&profile, moment, l).map(|v| rotate(v, rotated));
            let mut actual = [[0.0; 3]; 6];
            for &(dof, force) in load.nodal_forces_n() {
                actual[dof / 3][dof % 3] = force;
            }
            for (actual, expected) in actual
                .into_iter()
                .flatten()
                .zip(expected.into_iter().flatten())
            {
                close(actual, expected, 1e-17);
            }
            for (actual, expected) in transfer
                .moment_about_centre_n_m
                .into_iter()
                .zip(rotate([0.0, 0.0, -moment], rotated))
            {
                close(actual, expected, 1e-18);
            }
            for residual in transfer.moment_residual_n_m {
                close(residual, 0.0, 1e-18);
            }
            assert_eq!(transfer.source.model_identity, source.identity());
            assert!(transfer.vector_source.is_some());

            let reference_load = geometry
                .assemble(&[], VolumeLoad3::None, &sparse(&expected), cx)
                .unwrap();
            let mut native = patch
                .begin_transient(
                    &profile,
                    &other,
                    &problem,
                    None,
                    0.0,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 6],
                    native_config(),
                    &gate(),
                    cx,
                )
                .unwrap();
            let mut reference = reference_load
                .begin_transient(
                    &problem,
                    None,
                    0.0,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 6],
                    native_config(),
                    cx,
                )
                .unwrap();
            let step = native
                .advance(&profile, &other, 2e-7, &gate(), cx)
                .unwrap()
                .clone();
            let expected_step = reference.advance(&reference_load, 2e-7, cx).unwrap();
            for (&actual, &expected) in native
                .native()
                .state()
                .displacement
                .iter()
                .zip(&reference.state().displacement)
            {
                close(actual, expected, 1e-20);
            }
            close(
                step.traction_work_on_solid_j,
                expected_step.mechanical_work_j,
                1e-24,
            );
            close(
                step.spatial_traction_work_j,
                expected_step.mechanical_work_j,
                1e-24,
            );
            assert!(native.native().state().strain_energy_j > 0.0);
            // A source normal pressure alone cannot drive this tangential mode.
            let tangential_dof = if rotated { 5 } else { 4 };
            assert_ne!(native.native().state().displacement[tangential_dof], 0.0);
        });
    }
}

#[test]
fn prescribed_p1_rotation_receives_original_torsional_work() {
    let (source, profile) = torsion_profile();
    let a = profile.contact_radius_m();
    let nodes = nodes(4.0 * a);
    let material = material();
    let fixed: Vec<_> = (0..18).collect();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    support::with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
        let other = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let patch = AdhesiveP1Patch3::new(
            &geometry,
            frame(),
            grid(&geometry, vec![0.0, a]),
            controls(&profile),
            &gate(),
            cx,
        )
        .unwrap()
        .with_full_stick_tractions();
        let mut native = patch
            .begin_transient_with_supports(
                &profile,
                &other,
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                native_config(),
                &gate(),
                cx,
            )
            .unwrap();
        let alpha = 4e6;
        let acceleration: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(node, p)| [(3 * node, alpha * p[1]), (3 * node + 1, -alpha * p[0])])
            .collect();
        let h = 2e-7;
        let step = native
            .advance_with_supports(
                &profile,
                &other,
                AdhesiveP1SupportEndpoint3::Acceleration(&acceleration),
                h,
                &gate(),
                cx,
            )
            .unwrap();
        // The fixed traction performs M*delta_phi work on a rigid affine
        // rotation. Native Newmark displacement is h²*alpha/4 from rest.
        let moment = -(16.0 / 3.0) * source.law().reduced_shear_modulus_pa() * a.powi(3) * 4e-6;
        let expected = moment * (0.25 * h * h * alpha);
        assert!(expected < 0.0);
        close(
            step.traction_work_on_solid_j,
            expected,
            expected.abs() * 1e-12,
        );
        close(
            step.spatial_traction_work_j,
            expected,
            expected.abs() * 1e-12,
        );
        close(
            step.native.mechanical_work_j,
            expected,
            expected.abs() * 1e-12,
        );
        assert!(step.work_residual_j.abs() <= step.work_tolerance_j);
        assert_eq!(step.quadrature_work_estimate_j, 0.0);
        assert_eq!(step.other_external_work_j, 0.0);
    });
}

#[test]
fn torsional_peeling_changes_native_normal_loads_and_twisting_couple_together() {
    let source = mixed_source();
    let normal = source.normal_source();
    let delta = 1e-6;
    let limit = source.torsional_displacement_limit(delta).unwrap().unwrap();
    let angles = [0.0, 0.6 * limit.twist_angle_rad];
    let profiles = angles
        .map(|phi| AdhesiveTractionProfile::jkr_torsion_equilibrium(&source, delta, phi).unwrap());
    let responses = angles.map(|phi| {
        source
            .evaluate_torsion_displacement_controlled(delta, phi)
            .unwrap()
    });
    assert!(profiles[1].contact_radius_m() < profiles[0].contact_radius_m());
    let a0 = profiles[0].contact_radius_m();
    let l = 4.0 * a0;
    let nodes = nodes(l);
    let material = material();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &[],
        budget: Default::default(),
    };
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let g = source.torsional_reduced_shear_modulus_pa();
    let expected = core::array::from_fn::<_, 2, _>(|i| {
        let a = responses[i].contact_radius_m;
        // Independent integrals of the actual Hertz and flat-punch pressures,
        // plus the circumferential field, against the true face barycentrics.
        let hertz = (4.0 / 3.0) * e * a.powi(3) / normal.effective_radius_m();
        let tensile = -2.0 * e * a * (a * a / normal.effective_radius_m() - delta);
        let edge = a * (3.0 * hertz / 8.0 + tensile / 2.0) / (2.0 * l);
        let moment = -(16.0 / 3.0) * g * a.powi(3) * angles[i];
        let mut forces = [[0.0; 3]; 6];
        forces[0][2] = -(hertz + tensile) + 4.0 * edge;
        for row in &mut forces[1..5] {
            row[2] = -edge;
        }
        forces[1][1] = -moment / (4.0 * l);
        forces[2][0] = moment / (4.0 * l);
        forces[3][1] = moment / (4.0 * l);
        forces[4][0] = -moment / (4.0 * l);
        forces
    });
    support::with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
        let other = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let start_load = geometry
            .assemble(&[], VolumeLoad3::None, &sparse(&expected[0]), cx)
            .unwrap();
        let end_load = geometry
            .assemble(&[], VolumeLoad3::None, &sparse(&expected[1]), cx)
            .unwrap();
        let patch = AdhesiveP1Patch3::new(
            &geometry,
            frame(),
            grid(&geometry, vec![0.0, 0.5 * a0, a0]),
            controls(&profiles[0]),
            &gate(),
            cx,
        )
        .unwrap()
        .with_full_stick_tractions();
        let mut native = patch
            .begin_transient(
                &profiles[0],
                &other,
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                native_config(),
                &gate(),
                cx,
            )
            .unwrap();
        let mut reference = start_load
            .begin_transient(
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 6],
                native_config(),
                cx,
            )
            .unwrap();
        let step = native
            .advance(&profiles[1], &other, 2e-7, &gate(), cx)
            .unwrap()
            .clone();
        let reference_step = reference.advance(&end_load, 2e-7, cx).unwrap();
        for (&actual, &expected) in native
            .native()
            .state()
            .displacement
            .iter()
            .zip(&reference.state().displacement)
        {
            close(actual, expected, 1e-20);
        }
        close(
            step.traction_work_on_solid_j,
            reference_step.mechanical_work_j,
            1e-24,
        );
        close(
            step.spatial_traction_work_j,
            reference_step.mechanical_work_j,
            1e-24,
        );
        let transfer = native.current_transfer();
        close(
            transfer.force_world_n[2],
            -responses[1].normal_force_n,
            1e-17,
        );
        close(
            transfer.moment_about_centre_n_m[2],
            responses[1].torsional_moment_n_m,
            1e-18,
        );
        assert!(transfer.vector_source.is_some());
        assert_ne!(native.native().state().displacement[4], 0.0);
        assert!(step.work_residual_j.abs() <= step.work_tolerance_j);
        assert_eq!(step.quadrature_work_estimate_j, 0.0);
    });
}
