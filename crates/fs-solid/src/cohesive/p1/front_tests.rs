//! G0/G1/G3: matching crack fronts retain native work and material-point order.
use super::*;

fn mesh(shared: usize) -> (Vec<[f64; 3]>, Vec<[usize; 4]>, [[usize; 3]; 2]) {
    let mut nodes = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., -1.]];
    let b = if shared == 2 {
        nodes.extend([[0., 1., 0.], [0., 0., 1.]]);
        [0, 1, 4]
    } else {
        nodes.extend([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);
        [0, 4, 5]
    };
    let opposite = nodes.len() - 1;
    (
        nodes,
        vec![[0, 1, 2, 3], [b[0], b[1], b[2], opposite]],
        [[0, 1, 2], b],
    )
}

fn budget() -> TraceBudget {
    TraceBudget {
        maximum_nodes: 16,
        maximum_tetrahedra: 8,
        maximum_interfaces: 2,
        minimum_face_quality: 1e-10,
    }
}
fn law() -> BilinearCohesive {
    BilinearCohesive::new(1e9, 2e5, 100., 2e9, 0.01).unwrap()
}
fn pair(law: &BilinearCohesive, sides: [[usize; 3]; 2]) -> FacePair<'_> {
    FacePair {
        id: 23,
        side_a: sides[0],
        side_b: sides[1],
        law,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn partial_front_matches_restricted_analytic_triangle_matrix() {
    let law = law();
    for shared in [1, 2] {
        let (nodes, tets, [a, b]) = mesh(shared);
        let n = 3 * nodes.len();
        let faces = P1Interfaces::new_controlled(
            &nodes,
            &tets,
            &(0..n).collect::<Vec<_>>(),
            &[pair(&law, [a, b])],
            budget(),
            &mut || false,
        )
        .unwrap();
        assert!(
            faces
                .full_patches
                .iter()
                .all(|p| p.terms.len() == 6 * (3 - shared))
        );
        for j in 0..3 {
            for c in 0..3 {
                let mut u = vec![0.0; n];
                u[3 * b[j] + c] = 1.0;
                let jumps = faces.jumps(&u, &mut || false).unwrap();
                let force = faces.full_forces(&jumps, &mut || false).unwrap();
                let mut expected = vec![0.0; n];
                if a[j] != b[j] {
                    for i in 0..3 {
                        if a[i] != b[i] {
                            // Independent integral Ni*Nj over the reference triangle.
                            let entry = 0.5 * if i == j { 2.0 } else { 1.0 } / 12.0;
                            expected[3 * a[i] + c] -= entry;
                            expected[3 * b[i] + c] += entry;
                        }
                    }
                }
                for (&actual, expected) in force.iter().zip(expected) {
                    close(actual, expected, 2e-16);
                }
            }
        }
        let mut reordered = pair(&law, [a, b]);
        reordered.side_a.rotate_left(1);
        reordered.side_b.rotate_right(1);
        let again = P1Interfaces::new_controlled(
            &nodes,
            &tets,
            &(0..n).collect::<Vec<_>>(),
            &[reordered],
            budget(),
            &mut || false,
        )
        .unwrap();
        assert_eq!(faces.points(), again.points());
    }
}

#[test]
fn front_motion_cancels_before_arithmetic_and_maps_preserve_virtual_work() {
    let law = law();
    for shared in [1, 2] {
        let (nodes, tets, sides) = mesh(shared);
        let n = 3 * nodes.len();
        let free: Vec<_> = (0..n).filter(|dof| dof % 3 != 1).collect();
        let faces = P1Interfaces::new_controlled(
            &nodes,
            &tets,
            &free,
            &[pair(&law, sides)],
            budget(),
            &mut || false,
        )
        .unwrap();
        let u: Vec<_> = (0..n).map(|i| 0.001 * (0.7 * i as f64).sin()).collect();
        let traction = [[17., -31., 12.], [-20., 11., 30.], [6., 4., -7.]];
        let jumps = faces.jumps(&u, &mut || false).unwrap();
        let force = faces.full_forces(&traction, &mut || false).unwrap();
        close(
            force.iter().zip(&u).map(|(f, u)| f * u).sum(),
            jumps
                .iter()
                .zip(traction)
                .zip(faces.points())
                .map(|((g, t), p)| p.area_weight_m2 * dot(*g, t))
                .sum(),
            2e-16,
        );
        let mut large_front = u.clone();
        for i in 0..shared {
            for c in 0..3 {
                assert_eq!(force[3 * i + c], 0.0);
                large_front[3 * i + c] = 1e200;
            }
        }
        assert_eq!(faces.jumps(&large_front, &mut || false).unwrap(), jumps);
        let mut resultant = [0.0; 3];
        let mut moment = [0.0; 3];
        let mut rigid = Vec::new();
        for (x, f) in nodes.iter().zip(force.chunks_exact(3)) {
            let torque = cross(*x, [f[0], f[1], f[2]]);
            for c in 0..3 {
                resultant[c] += f[c];
                moment[c] += torque[c];
            }
            rigid.extend(cross([0.003, 0.001, -0.002], *x).map(|v| v + 0.01));
        }
        for value in resultant.into_iter().chain(moment) {
            close(value, 0.0, 1e-14);
        }
        for jump in faces.jumps(&rigid, &mut || false).unwrap() {
            for v in jump {
                close(v, 0., 1e-17);
            }
        }
        let reduced: Vec<_> = free.iter().map(|&dof| u[dof]).collect();
        let mut prescribed = u.clone();
        for &dof in &free {
            prescribed[dof] = 0.0;
        }
        let full = faces.expand(&reduced, &prescribed, &mut || false).unwrap();
        assert_eq!(full, u);
        let offsets = faces.prescribed_jumps(&prescribed, &mut || false).unwrap();
        for ((patch, offset), expected) in faces.patches().iter().zip(offsets).zip(jumps) {
            for (actual, expected) in patch
                .jump(&reduced, offset)
                .unwrap()
                .into_iter()
                .zip(expected)
            {
                close(actual, expected, 1e-17);
            }
        }
    }
}

#[test]
fn partial_front_native_fracture_work_and_compressive_reclosure_are_retained() {
    let law = law();
    let (nodes, tets, sides) = mesh(2);
    let n = 3 * nodes.len();
    let faces = P1Interfaces::new_controlled(
        &nodes,
        &tets,
        &(0..n).collect::<Vec<_>>(),
        &[pair(&law, sides)],
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut points = vec![law.bonded_point(); 3];
    let mut old_u = vec![0.0; n];
    let mut dissipated = 0.0;
    for opening in [0.00003, 0.0004, 0.0012, 0.003, 0.009, 0.003, -0.00003] {
        let mut u = vec![0.0; n];
        u[3 * sides[1][2] + 2] = opening;
        let jumps = faces.jumps(&u, &mut || false).unwrap();
        let mut means = Vec::new();
        let mut material_work = 0.0;
        for ((point, jump), location) in points.iter_mut().zip(jumps).zip(faces.points()) {
            let previous = *point;
            let (step, mean) = point.trial_with_mean_traction(jump, 1e-7).unwrap();
            assert!(step.point.damage() >= previous.damage());
            material_work += location.area_weight_m2 * step.energy.work;
            dissipated += location.area_weight_m2 * step.energy.dissipated;
            means.push(mean);
            *point = step.point;
        }
        let force = faces.full_forces(&means, &mut || false).unwrap();
        let nodal_work: f64 = force
            .iter()
            .zip(u.iter().zip(&old_u))
            .map(|(f, (u, old))| f * (u - old))
            .sum();
        close(nodal_work, material_work, 1e-8);
        for node in [0, 1] {
            assert_eq!(&force[3 * node..3 * node + 3], &[0.0; 3]);
        }
        old_u = u;
    }
    close(dissipated, 0.5 * law.fracture_energy_j_m2(), 1e-9);
    assert!(
        points
            .iter()
            .all(|point| point.failed() && point.traction_pa()[0] < 0.0)
    );
}

#[test]
fn ordered_trace_points_frames_and_actions_survive_global_node_renumbering() {
    let law = law();
    let (nodes, tets, mut sides) = mesh(1);
    sides[0].rotate_right(1);
    sides[1].rotate_left(1);
    let n = 3 * nodes.len();
    let free: Vec<_> = (0..n).collect();
    let before = P1Interfaces::new_ordered_controlled(
        &nodes,
        &tets,
        &free,
        &[pair(&law, sides)],
        budget(),
        &mut || false,
    )
    .unwrap();
    let rename = [5, 1, 6, 2, 4, 0, 3];
    let mut next_nodes = vec![[0.; 3]; nodes.len()];
    let u: Vec<_> = (0..n).map(|i| 0.001 * (i as f64).sin()).collect();
    let mut next_u = vec![0.0; n];
    for node in 0..nodes.len() {
        next_nodes[rename[node]] = nodes[node];
        next_u[3 * rename[node]..3 * rename[node] + 3].copy_from_slice(&u[3 * node..3 * node + 3]);
    }
    let next_tets: Vec<_> = tets
        .iter()
        .map(|tet| tet.map(|node| rename[node]))
        .collect();
    let next_sides = sides.map(|side| side.map(|node| rename[node]));
    let after = P1Interfaces::new_ordered_controlled(
        &next_nodes,
        &next_tets,
        &free,
        &[pair(&law, next_sides)],
        budget(),
        &mut || false,
    )
    .unwrap();
    for (a, b) in before.points().iter().zip(after.points()) {
        assert_eq!(a.interface_id, b.interface_id);
        assert_eq!(a.point_in_face, b.point_in_face);
        assert_eq!(a.location_m, b.location_m);
        assert_eq!(a.axes, b.axes);
        assert_eq!(a.area_weight_m2, b.area_weight_m2);
        assert_eq!(a.shape_values, b.shape_values);
        assert_eq!(
            a.paired_nodes.map(|side| side.map(|node| rename[node])),
            b.paired_nodes
        );
    }
    let canonical = P1Interfaces::new_controlled(
        &next_nodes,
        &next_tets,
        &free,
        &[pair(&law, next_sides)],
        budget(),
        &mut || false,
    )
    .unwrap();
    assert_ne!(canonical.points()[0].axes, after.points()[0].axes);
    let ga = before.jumps(&u, &mut || false).unwrap();
    let gb = after.jumps(&next_u, &mut || false).unwrap();
    for (a, b) in ga.iter().zip(gb) {
        for (a, b) in a.iter().zip(b) {
            close(*a, b, 1e-17);
        }
    }
    let traction = [[11., -2., 3.], [-5., 1., 3.], [4., 5., -2.]];
    let fa = before.full_forces(&traction, &mut || false).unwrap();
    let fb = after.full_forces(&traction, &mut || false).unwrap();
    for node in 0..nodes.len() {
        for c in 0..3 {
            close(fa[3 * node + c], fb[3 * rename[node] + c], 1e-15);
        }
    }
}

fn source() -> ResolvedBilinearCohesive {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
    };
    use fs_material::cohesive::source::*;
    use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
    use fs_qty::Pressure;
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 2e5),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured crack-front test; no measured claim".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        texture_frame: "fixed reference frame".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit cohesive seam".into(),
            third_body: None,
            environment: "fixed".into(),
            history: "virgin".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        0.01,
    )
    .unwrap()
}

#[test]
fn ordered_and_canonical_partial_fronts_retain_the_original_interface_source() {
    let source = source();
    let (nodes, tets, mut sides) = mesh(2);
    sides[0].rotate_right(1);
    let receiving = MaterialFacePair {
        id: 23,
        side_a: sides[0],
        side_b: sides[1],
        source: &source,
        surface_materials: [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
    };
    for ordered in [false, true] {
        let free: Vec<_> = (0..3 * nodes.len()).collect();
        let faces = if ordered {
            P1Interfaces::from_interfaces_ordered_controlled(
                &nodes,
                &tets,
                &free,
                &[receiving],
                budget(),
                &mut || false,
            )
        } else {
            P1Interfaces::from_interfaces_controlled(
                &nodes,
                &tets,
                &free,
                &[receiving],
                budget(),
                &mut || false,
            )
        }
        .unwrap();
        for patch in faces.patches().iter().chain(&faces.full_patches) {
            assert!(core::ptr::eq(patch.source().unwrap(), &source));
        }
        let mut wrong = receiving;
        wrong.surface_materials.swap(0, 1);
        assert!(
            P1Interfaces::from_interfaces_ordered_controlled(
                &nodes,
                &tets,
                &free,
                &[wrong],
                budget(),
                &mut || false
            )
            .is_err()
        );
    }
}

#[test]
fn full_ties_malformed_pairings_and_nonmatching_shared_dofs_still_refuse() {
    use super::nonmatching::{NonmatchingBudget, NonmatchingSurface};
    let law = law();
    let (nodes, tets, sides) = mesh(2);
    let free: Vec<_> = (0..3 * nodes.len()).collect();
    let valid = pair(&law, sides);
    for bad in [
        pair(&law, [sides[0], sides[0]]),
        pair(&law, [sides[0], [0, 4, 4]]),
    ] {
        assert!(
            P1Interfaces::new_ordered_controlled(
                &nodes,
                &tets,
                &free,
                &[bad],
                budget(),
                &mut || false
            )
            .is_err()
        );
    }
    assert!(
        P1Interfaces::new_ordered_controlled(
            &nodes,
            &tets,
            &free,
            &[valid, valid],
            budget(),
            &mut || false
        )
        .is_err()
    );
    let mut ambiguous = nodes.clone();
    ambiguous[4] = ambiguous[0];
    assert!(
        P1Interfaces::new_ordered_controlled(
            &ambiguous,
            &tets,
            &free,
            &[valid],
            budget(),
            &mut || false
        )
        .is_err()
    );
    for stop in [1, 5, 8] {
        let mut polls = 0;
        assert!(matches!(
            P1Interfaces::new_ordered_controlled(
                &nodes,
                &tets,
                &free,
                &[valid],
                budget(),
                &mut || {
                    polls += 1;
                    polls == stop
                }
            ),
            Err(Error::Cancelled)
        ));
    }
    let nonmatching = NonmatchingSurface {
        id: 23,
        side_a: &[sides[0]],
        side_b: &[sides[1]],
        law: &law,
    };
    let cap = NonmatchingBudget {
        mesh: budget(),
        maximum_surface_faces: 2,
        maximum_pair_checks: 4,
        maximum_quadrature_points: 3,
        relative_geometry_tolerance: 1e-10,
    };
    assert!(matches!(
        P1Interfaces::from_nonmatching_controlled(
            &nodes,
            &tets,
            &free,
            &[nonmatching],
            cap,
            &mut || false
        ),
        Err(Error::Invalid(
            "cohesive traces must have distinct vertices"
        ))
    ));
    P1Interfaces::new_ordered_controlled(&nodes, &tets, &free, &[valid], budget(), &mut || false)
        .unwrap();
}
