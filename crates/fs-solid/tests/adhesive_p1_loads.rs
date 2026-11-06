//! Original adhesive fields on real exterior P1 traces and native dynamics.
#[path = "../../fs-contact/tests/normal_work_derivative/fixtures.rs"]
mod sources;
#[path = "tet10/support.rs"]
mod support;

use core::f64::consts::{FRAC_PI_2, PI, TAU};
use fs_blake3::ContentHash;
use fs_contact::normal_patch::{
    maugis::{MaugisEnvelope, SourceMaugisSphereContact},
    traction::{
        AdhesiveTractionProfile,
        sectors::{PolarTransferBudget, RadialMomentBudget},
    },
};
use fs_exec::CancelGate;
use fs_solid::linear3::{
    TetElasticMaterial, TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig,
    TetThermalStrainField, TetThermalStrainState,
};
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3, adhesive::*, transient::TransientConfig3};

const TETS: [[usize; 4]; 4] = [[0, 1, 2, 5], [0, 2, 3, 5], [0, 3, 4, 5], [0, 4, 1, 5]];
fn nodes(l: f64) -> [[f64; 3]; 6] {
    [
        [0.0; 3],
        [l, 0.0, 0.0],
        [0.0, l, 0.0],
        [-l, 0.0, 0.0],
        [0.0, -l, 0.0],
        [0.0, 0.0, -l],
    ]
}
fn gate() -> CancelGate {
    CancelGate::new_clock_free()
}
fn frame() -> AdhesiveP1Frame3 {
    AdhesiveP1Frame3 {
        centre_m: [0.0; 3],
        x_axis: [1.0, 0.0, 0.0],
        y_axis: [0.0, -1.0, 0.0],
    }
}
fn grid(g: &TetLoadGeometry3<'_>, radii_m: Vec<f64>) -> AdhesiveP1Grid3 {
    let one_ring: Vec<_> = [[0, 1, 4], [0, 3, 4], [0, 2, 3], [0, 1, 2]]
        .into_iter()
        .map(|nodes| g.faces().iter().position(|f| f.nodes == nodes).unwrap())
        .collect();
    let faces = (1..radii_m.len())
        .flat_map(|_| one_ring.iter().copied())
        .collect();
    AdhesiveP1Grid3 {
        radii_m,
        angles_rad: vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU],
        faces,
    }
}
fn controls(profile: &AdhesiveTractionProfile) -> AdhesiveP1Controls3 {
    let r = profile.interaction_radius_m();
    let whole = profile.annulus(0.0, r).unwrap();
    let f = whole.hertz_force_n - whole.adhesive_force_n;
    AdhesiveP1Controls3 {
        maximum_nodes: 32,
        geometry_tolerance_m: 1e-12 * r,
        maximum_frame_error: 1e-12,
        transfer: PolarTransferBudget {
            radial_moment: RadialMomentBudget {
                absolute_tolerance_n_m: 1e-12 * f * r,
                relative_tolerance: 1e-10,
                maximum_evaluations: 65_536,
            },
            force_absolute_tolerance_n: 1e-11 * f,
            moment_absolute_tolerance_n_m: 1e-11 * f * r,
        },
        projection_force_tolerance_n: 1e-10 * f,
        projection_moment_tolerance_n_m: 1e-10 * f * r,
        work_absolute_tolerance_j: 1e-24,
        work_relative_tolerance: 1e-10,
        maximum_quadrature_work_estimate_j: 1e-14,
    }
}
fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(1500.0, 2e7, 0.25, [71; 32]).unwrap()
}
fn native_config() -> TransientConfig3 {
    TransientConfig3 {
        absolute_force_tolerance_n: 1e-14,
        absolute_energy_tolerance_j: 1e-18,
        relative_energy_tolerance: 1e-9,
        maximum_displacement_m: 1e-3,
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-11,
            maximum_iterations: 512,
        },
        ..Default::default()
    }
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + 2e-9 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn zero_load_profile() -> AdhesiveTractionProfile {
    let source = sources::jkr();
    let normal = source.normal_state();
    let a = (4.5 * PI * normal.adhesion_energy_j_per_m2() * source.effective_radius_m().powi(2)
        / normal.elastic().reduced_modulus_pa())
    .cbrt();
    AdhesiveTractionProfile::jkr_attached(&source, a * a / (3.0 * source.effective_radius_m()))
        .unwrap()
}
// Independent closed-form integral on the four true top-face barycentric
// functions: edge N=|x|/L or |y|/L, centre N=1-(|x|+|y|)/L.
fn jkr_nodal(profile: &AdhesiveTractionProfile, l: f64) -> Vec<(usize, f64)> {
    let a = profile.contact_radius_m();
    let whole = profile.annulus(0.0, a).unwrap();
    let bending = 3.0 * whole.hertz_force_n / 8.0 + whole.adhesive_force_n / 2.0;
    let edge = a * bending / (2.0 * l);
    vec![
        (2, -whole.normal_force_n + 4.0 * edge),
        (5, -edge),
        (8, -edge),
        (11, -edge),
        (14, -edge),
    ]
}

#[test]
fn zero_resultant_patch_keeps_its_p1_load_pattern_and_drives_actual_native_deformation() {
    let profile = zero_load_profile();
    let a = profile.contact_radius_m();
    let l = 4.0 * a;
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
        let g = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
        let other = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut patch = AdhesiveP1Patch3::new(
            &g,
            frame(),
            grid(&g, vec![0.0, 0.3 * a, a]),
            controls(&profile),
            &gate(),
            cx,
        )
        .unwrap();
        let (load, transfer) = patch.assemble(&profile, &other, &gate(), cx).unwrap();
        let reference_load = g
            .assemble(&[], VolumeLoad3::None, &jkr_nodal(&profile, l), cx)
            .unwrap();
        for ((a, x), (b, y)) in load
            .nodal_forces_n()
            .iter()
            .zip(reference_load.nodal_forces_n())
        {
            assert_eq!(a, b);
            close(*x, *y, 1e-17);
        }
        assert_eq!(load.nodal_forces_n().len(), 5);
        assert!(load.nodal_forces_n()[0].1 < 0.0);
        assert!(load.nodal_forces_n()[1..].iter().all(|row| row.1 > 0.0));
        assert!(transfer.force_world_n.into_iter().all(|x| x.abs() < 1e-14));
        assert!(
            transfer
                .moment_about_centre_n_m
                .into_iter()
                .all(|x| x.abs() < 1e-18)
        );
        assert_eq!(transfer.source.model_identity, profile.model_identity());
        assert_eq!(transfer.source.radial_evaluations, 0);
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
        for step in 1..5 {
            let actual = native
                .advance(&profile, &other, 2e-7, &gate(), cx)
                .unwrap()
                .clone();
            let expected = reference.advance(&reference_load, 2e-7, cx).unwrap();
            assert_eq!(native.native().state().step, step);
            for (a, b) in native
                .native()
                .state()
                .displacement
                .iter()
                .zip(&reference.state().displacement)
            {
                close(*a, *b, 1e-20);
            }
            for (a, b) in native
                .native()
                .state()
                .velocity
                .iter()
                .zip(&reference.state().velocity)
            {
                close(*a, *b, 1e-14);
            }
            close(
                actual.traction_work_on_solid_j,
                expected.mechanical_work_j,
                1e-24,
            );
            assert!(actual.work_residual_j.abs() <= actual.work_tolerance_j);
            assert_eq!(actual.quadrature_work_estimate_j, 0.0);
        }
        assert!(native.native().nodal_motion(cx).unwrap()[0][0][2] < 0.0);
        assert!(native.native().state().strain_energy_j > 0.0);
    });
}

#[test]
fn changing_source_loads_retain_native_thermal_damping_and_supported_reactions() {
    let source = sources::jkr();
    let start = AdhesiveTractionProfile::jkr_attached(&source, 0.0).unwrap();
    let end = AdhesiveTractionProfile::jkr_attached(&source, 2e-6).unwrap();
    let r = end.contact_radius_m();
    let l = 4.0 * r;
    let nodes = nodes(l);
    let material = material();
    let fixed = [3, 4, 7, 15, 16, 17];
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &fixed,
        budget: Default::default(),
    };
    let thermal = TetThermalStrainState::try_new(
        301.0,
        300.0,
        [1e-7, 1e-7, 1e-7, 0.0, 0.0, 0.0],
        material.material_state_identity,
        ContentHash([72; 32]),
        0.01,
    )
    .unwrap();
    let cfg = TransientConfig3 {
        mass_damping_per_s: 10.0,
        stiffness_damping_s: 1e-7,
        ..native_config()
    };
    support::with_cx(|cx| {
        let g = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
        let other = g
            .assemble(
                &[],
                VolumeLoad3::Uniform([0.0, 0.0, -3.0]),
                &[(2, 1e-5)],
                cx,
            )
            .unwrap();
        let reference_load = |p: &AdhesiveTractionProfile| {
            let mut summed = [0.0; 18];
            for &(d, f) in other.nodal_forces_n() {
                summed[d] += f;
            }
            for (d, f) in jkr_nodal(p, l) {
                summed[d] += f;
            }
            let forces: Vec<_> = summed
                .into_iter()
                .enumerate()
                .filter(|row| row.1 != 0.0)
                .collect();
            g.assemble(&[], VolumeLoad3::None, &forces, cx).unwrap()
        };
        let mut velocity = [[0.0; 3]; 6];
        velocity[0][2] = 0.001;
        let mut reference = reference_load(&start)
            .begin_transient(
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                3.0,
                &[[0.0; 3]; 6],
                &velocity,
                cfg,
                cx,
            )
            .unwrap();
        let patch = AdhesiveP1Patch3::new(
            &g,
            frame(),
            grid(&g, vec![0.0, r]),
            controls(&end),
            &gate(),
            cx,
        )
        .unwrap();
        let mut native = patch
            .begin_transient(
                &start,
                &other,
                &problem,
                Some(TetThermalStrainField::Uniform(&thermal)),
                3.0,
                &[[0.0; 3]; 6],
                &velocity,
                cfg,
                &gate(),
                cx,
            )
            .unwrap();
        for p in [&end, &start, &end] {
            let a = native
                .advance(p, &other, 2e-7, &gate(), cx)
                .unwrap()
                .clone();
            let b = reference.advance(&reference_load(p), 2e-7, cx).unwrap();
            close(a.native.damping_loss_j, b.damping_loss_j, 1e-24);
            assert!(a.native.damping_loss_j > 0.0);
            close(
                a.native.thermal_equivalent_work_j,
                b.thermal_equivalent_work_j,
                1e-24,
            );
            close(
                a.traction_work_on_solid_j + a.other_external_work_j,
                b.mechanical_work_j,
                1e-23,
            );
            for (a, b) in native
                .native()
                .state()
                .displacement
                .iter()
                .zip(&reference.state().displacement)
            {
                close(*a, *b, 1e-20);
            }
            let response = native
                .native()
                .recover_response(ContentHash([72; 32]), 0.01, cx)
                .unwrap();
            let expected = reference
                .recover_response(ContentHash([72; 32]), 0.01, cx)
                .unwrap();
            for (a, b) in response
                .support_reactions_n
                .iter()
                .flatten()
                .zip(expected.support_reactions_n.iter().flatten())
            {
                close(*a, *b, 1e-12);
            }
        }
        // A completely fixed receiver still retains the distributed pressure
        // in the original static support reactions, with exactly zero work.
        let fixed: Vec<_> = (0..18).collect();
        let clamped = TetLinearElasticProblem {
            fixed_dofs: &fixed,
            ..problem
        };
        let zero_thermal = TetThermalStrainState::try_new(
            300.0,
            300.0,
            [0.0; 6],
            material.material_state_identity,
            ContentHash([72; 32]),
            0.01,
        )
        .unwrap();
        let mut patch = AdhesiveP1Patch3::new(
            &g,
            frame(),
            grid(&g, vec![0.0, r]),
            controls(&end),
            &gate(),
            cx,
        )
        .unwrap();
        let (load, _) = patch.assemble(&end, &other, &gate(), cx).unwrap();
        let solved = load
            .solve(
                &clamped,
                TetThermalStrainField::Uniform(&zero_thermal),
                Default::default(),
                ContentHash([72; 32]),
                0.01,
                cx,
            )
            .unwrap();
        for &(d, f) in load.nodal_forces_n() {
            close(solved.support_reactions_n[d / 3][d % 3], -f, 1e-17);
        }
    });
}

fn maugis() -> SourceMaugisSphereContact {
    let original = sources::jkr();
    let normal = original.normal_state().clone();
    let r = original.effective_radius_m();
    let a0 = (PI * normal.adhesion_energy_j_per_m2() * r * r
        / (4.0 * normal.elastic().reduced_modulus_pa() / 3.0))
        .cbrt();
    let d0 = a0 * a0 / r;
    SourceMaugisSphereContact::new(
        normal,
        r,
        2.0 * d0 / (PI * 0.2),
        MaugisEnvelope {
            maximum_indentation_m: original.envelope().maximum_indentation_m,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.5,
            maximum_adhesive_stress_to_yield: 0.9,
            maximum_interaction_range_to_radius: 0.01,
        },
    )
    .unwrap()
}

#[test]
fn maugis_contact_gap_and_detachment_keep_the_original_polar_source_and_quadrature_gate() {
    let source = maugis();
    let contact =
        AdhesiveTractionProfile::maugis(&source, &source.evaluate_contact_energy(0.0).unwrap())
            .unwrap();
    let bounds = source.noncontact_limits().unwrap().unwrap();
    let gap = AdhesiveTractionProfile::maugis(
        &source,
        &source
            .evaluate_noncontact(
                bounds
                    .minimum_indentation_m
                    .midpoint(bounds.first_contact_indentation_m),
            )
            .unwrap(),
    )
    .unwrap();
    let detached = AdhesiveTractionProfile::maugis(
        &source,
        &source
            .evaluate_detached(-2.0 * source.interaction_range_m())
            .unwrap(),
    )
    .unwrap();
    let r = 1.2
        * contact
            .interaction_radius_m()
            .max(gap.interaction_radius_m());
    let nodes = nodes(4.0 * r);
    let material = material();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &[],
        budget: Default::default(),
    };
    support::with_cx(|cx| {
        let g = TetLoadGeometry3::new(&nodes, &TETS, problem.budget, 16, cx).unwrap();
        let other = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let map = grid(
            &g,
            vec![
                0.0,
                0.4 * contact.contact_radius_m(),
                contact.contact_radius_m(),
                r,
            ],
        );
        let mut patch =
            AdhesiveP1Patch3::new(&g, frame(), map.clone(), controls(&contact), &gate(), cx)
                .unwrap();
        for p in [&contact, &gap, &detached] {
            let (load, transfer) = patch.assemble(p, &other, &gate(), cx).unwrap();
            let (force, moment) = load.resultant([0.0; 3], cx).unwrap();
            close(force[2], -p.normal_force_n(), 1e-15);
            assert!(moment.into_iter().all(|v| v.abs() < 1e-18));
            assert_eq!(transfer.source.model_identity, source.identity());
        }
        assert_eq!(gap.contact_radius_m(), 0.0);
        let mut velocity = [[0.0; 3]; 6];
        velocity[0][2] = 0.001;
        let mut native = patch
            .begin_transient(
                &contact,
                &other,
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &velocity,
                native_config(),
                &gate(),
                cx,
            )
            .unwrap();
        let first = native
            .advance(&contact, &other, 2e-7, &gate(), cx)
            .unwrap()
            .clone();
        assert!(first.transfer.source.radial_evaluations > 0);
        assert!(first.quadrature_work_estimate_j > 0.0);
        native.advance(&gap, &other, 2e-7, &gate(), cx).unwrap();
        native
            .advance(&detached, &other, 2e-7, &gate(), cx)
            .unwrap();
        let empty = native
            .advance(&detached, &other, 2e-7, &gate(), cx)
            .unwrap();
        assert_eq!(empty.spatial_traction_work_j, 0.0);
        assert_eq!(empty.traction_work_on_solid_j, 0.0);
        assert_eq!(empty.quadrature_work_estimate_j, 0.0);
        let mut exact = controls(&contact);
        exact.maximum_quadrature_work_estimate_j = 0.0;
        let patch = AdhesiveP1Patch3::new(&g, frame(), map, exact, &gate(), cx).unwrap();
        let mut refused = patch
            .begin_transient(
                &contact,
                &other,
                &problem,
                None,
                0.0,
                &[[0.0; 3]; 6],
                &velocity,
                native_config(),
                &gate(),
                cx,
            )
            .unwrap();
        let before = refused.native().state().clone();
        assert!(matches!(
            refused.advance(&contact, &other, 2e-7, &gate(), cx),
            Err(AdhesiveP1Error3::QuadratureWork { .. })
        ));
        assert_eq!(refused.native().state(), &before);
    });
}

#[test]
fn whole_sector_containment_and_reference_frame_refuse_real_geometric_misassignment() {
    let profile = AdhesiveTractionProfile::jkr_attached(&sources::jkr(), 0.0).unwrap();
    let a = profile.contact_radius_m();
    // Both radial endpoints and the centre are inside each top triangle. The
    // arc at 45 degrees crosses x+y=L: an endpoint-only test would be wrong.
    let bad_nodes = nodes(1.25 * a);
    let good_nodes = nodes(4.0 * a);
    support::with_cx(|cx| {
        let bad = TetLoadGeometry3::new(&bad_nodes, &TETS, Default::default(), 16, cx).unwrap();
        assert!(matches!(
            AdhesiveP1Patch3::new(
                &bad,
                frame(),
                grid(&bad, vec![0.0, a]),
                controls(&profile),
                &gate(),
                cx
            ),
            Err(AdhesiveP1Error3::Coverage { .. })
        ));
        let g = TetLoadGeometry3::new(&good_nodes, &TETS, Default::default(), 16, cx).unwrap();
        let mut wrong = frame();
        wrong.y_axis[1] = 1.0;
        assert!(
            AdhesiveP1Patch3::new(
                &g,
                wrong,
                grid(&g, vec![0.0, a]),
                controls(&profile),
                &gate(),
                cx
            )
            .is_err()
        );
        let mut wrong = frame();
        wrong.centre_m[2] = 0.01 * a;
        assert!(
            AdhesiveP1Patch3::new(
                &g,
                wrong,
                grid(&g, vec![0.0, a]),
                controls(&profile),
                &gate(),
                cx
            )
            .is_err()
        );
        let mut map = grid(&g, vec![0.0, a]);
        map.faces[0] = map.faces[2];
        assert!(matches!(
            AdhesiveP1Patch3::new(&g, frame(), map, controls(&profile), &gate(), cx),
            Err(AdhesiveP1Error3::Coverage { .. })
        ));
        let mut limited = controls(&profile);
        limited.maximum_nodes = 5;
        assert!(
            AdhesiveP1Patch3::new(&g, frame(), grid(&g, vec![0.0, a]), limited, &gate(), cx)
                .is_err()
        );

        // Rotate the full native geometry/frame and translate its contact
        // centre. Source pressure, true nodal shape functions, torque and
        // material virtual work must transform together.
        let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
        let centre = [2.0 * a, -3.0 * a, 4.0 * a];
        let rotated_nodes = good_nodes.map(|p| {
            let q = rotate(p);
            std::array::from_fn(|i| q[i] + centre[i])
        });
        let rotated =
            TetLoadGeometry3::new(&rotated_nodes, &TETS, Default::default(), 16, cx).unwrap();
        let other = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let rotated_other = rotated.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut original = AdhesiveP1Patch3::new(
            &g,
            frame(),
            grid(&g, vec![0.0, a]),
            controls(&profile),
            &gate(),
            cx,
        )
        .unwrap();
        let mut transformed = AdhesiveP1Patch3::new(
            &rotated,
            AdhesiveP1Frame3 {
                centre_m: centre,
                x_axis: rotate(frame().x_axis),
                y_axis: rotate(frame().y_axis),
            },
            grid(&rotated, vec![0.0, a]),
            controls(&profile),
            &gate(),
            cx,
        )
        .unwrap();
        let (original, _) = original.assemble(&profile, &other, &gate(), cx).unwrap();
        let (transformed, _) = transformed
            .assemble(&profile, &rotated_other, &gate(), cx)
            .unwrap();
        let mut f0 = [[0.0; 3]; 6];
        let mut f1 = [[0.0; 3]; 6];
        for &(d, f) in original.nodal_forces_n() {
            f0[d / 3][d % 3] = f;
        }
        for &(d, f) in transformed.nodal_forces_n() {
            f1[d / 3][d % 3] = f;
        }
        let mut work = [0.0; 2];
        for i in 0..6 {
            for (f, expected) in f1[i].into_iter().zip(rotate(f0[i])) {
                close(f, expected, 1e-17);
            }
            let p = good_nodes[i];
            let u = [
                0.2 * p[0],
                -0.1 * p[1],
                0.001 * (p[0].abs() + 0.5 * p[1].abs()),
            ];
            work[0] += (0..3).map(|c| f0[i][c] * u[c]).sum::<f64>();
            let u = rotate(u);
            work[1] += (0..3).map(|c| f1[i][c] * u[c]).sum::<f64>();
        }
        close(work[0], work[1], 1e-22);
        let (force, moment) = transformed.resultant([0.0; 3], cx).unwrap();
        close(force[0], -profile.normal_force_n(), 1e-17);
        close(moment[1], centre[2] * force[0], 1e-20);
        close(moment[2], -centre[1] * force[0], 1e-20);
    });
}

#[test]
fn late_receiver_cancellation_and_erased_background_load_preserve_the_native_checkpoint() {
    let source = sources::jkr();
    let a = AdhesiveTractionProfile::jkr_attached(&source, 0.0).unwrap();
    let b = AdhesiveTractionProfile::jkr_attached(&source, 1e-6).unwrap();
    let r = b.contact_radius_m();
    let nodes = nodes(4.0 * r);
    let material = material();
    let problem = TetLinearElasticProblem {
        nodes_m: &nodes,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &[],
        budget: Default::default(),
    };
    support::with_cx(|cx| {
        let g = TetLoadGeometry3::new(&nodes, &TETS, Default::default(), 16, cx).unwrap();
        let other = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let make = || {
            AdhesiveP1Patch3::new(
                &g,
                frame(),
                grid(&g, vec![0.0, r]),
                controls(&b),
                &gate(),
                cx,
            )
            .unwrap()
            .begin_transient(
                &a,
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
            .unwrap()
        };
        let mut native = make();
        let mut reference = make();
        native.advance(&a, &other, 2e-7, &gate(), cx).unwrap();
        reference.advance(&a, &other, 2e-7, &gate(), cx).unwrap();
        let state = native.native().state().clone();
        let frame = native.last_step().cloned();
        let transfer = *native.current_transfer();
        let refused: Result<(), AdhesiveP1Error3> =
            native.advance_with_admission(&b, &other, 2e-7, &gate(), cx, |_, _| {
                Err(AdhesiveP1Error3::Invalid("receiver refused"))
            });
        assert!(refused.is_err());
        let cancelled = gate();
        let refused: Result<(), AdhesiveP1Error3> =
            native.advance_with_admission(&b, &other, 2e-7, &cancelled, cx, |_, _| {
                cancelled.request();
                Ok(())
            });
        assert!(matches!(refused, Err(AdhesiveP1Error3::Cancelled)));
        assert_eq!(native.native().state(), &state);
        assert_eq!(native.last_step(), frame.as_ref());
        assert_eq!(native.current_transfer(), &transfer);
        assert_eq!(
            native.advance(&b, &other, 2e-7, &gate(), cx).unwrap(),
            reference.advance(&b, &other, 2e-7, &gate(), cx).unwrap()
        );
        assert_eq!(native.native().state(), reference.native().state());

        let background = g
            .assemble(
                &[],
                VolumeLoad3::None,
                &[(2, 1e20), (5, 1e20), (8, 1e20), (11, 1e20), (14, 1e20)],
                cx,
            )
            .unwrap();
        let mut patch = AdhesiveP1Patch3::new(
            &g,
            frame(),
            grid(&g, vec![0.0, r]),
            controls(&b),
            &gate(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            patch.assemble(&a, &background, &gate(), cx),
            Err(AdhesiveP1Error3::Projection { .. })
        ));
        patch.assemble(&a, &other, &gate(), cx).unwrap();
    });
}
