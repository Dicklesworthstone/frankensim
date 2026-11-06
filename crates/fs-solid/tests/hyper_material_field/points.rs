//! G0/G1: native quadrature owners retain independent strains and histories.

use super::*;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::{incompressible_plane_stress, plane_stress};
use fs_solid::hyper2d::{PlaneMaterial, material_point_kinematics};

const POINT_FORMS: [Formulation; 3] = [
    Formulation::PlaneStrain,
    Formulation::PlaneStress,
    Formulation::Incompressible,
];

#[derive(Clone)]
struct PointField<L> {
    elements: Vec<Vec<L>>,
    seeds: Option<Vec<Vec<[f64; 3]>>>,
    omitted: Option<(usize, usize)>,
    declare_elements: bool,
}

impl<L> PointField<L> {
    fn new(elements: Vec<Vec<L>>) -> Self {
        Self {
            elements,
            seeds: None,
            omitted: None,
            declare_elements: true,
        }
    }
}

impl<L: PlaneMaterial> HyperelasticMaterials for PointField<L> {
    fn material_count(&self) -> Option<usize> {
        self.declare_elements.then_some(self.elements.len())
    }

    fn material_at(&self, _: usize) -> Option<&dyn PlaneMaterial> {
        // A caller using the element's first history is a real wiring defect.
        None
    }

    fn material_point_count(&self, element: usize) -> Option<usize> {
        self.elements.get(element).map(Vec::len)
    }

    fn material_at_point(&self, element: usize, point: usize) -> Option<&dyn PlaneMaterial> {
        if self.omitted == Some((element, point)) {
            return None;
        }
        self.elements
            .get(element)?
            .get(point)
            .map(|law| law as &dyn PlaneMaterial)
    }

    fn transverse_seed_at(&self, element: usize, point: usize) -> Option<[f64; 3]> {
        self.seeds.as_ref()?.get(element)?.get(point).copied()
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn cheap_force(
    form: Formulation,
    p: &HyperProblem<'_>,
    u: &[f64],
    load: f64,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, SolidError> {
    let local = PlaneStressSettings::default();
    match form {
        Formulation::PlaneStrain => p.unconstrained_residual(u, load),
        Formulation::PlaneStress => p.plane_stress_unconstrained_residual(u, load, local, cx),
        Formulation::Incompressible => {
            p.incompressible_plane_stress_unconstrained_residual(u, load, local, cx)
        }
        Formulation::MeanDilatation => unreachable!(),
    }
}

#[test]
fn g0_q1_point_fields_share_force_energy_tangent_and_recovery_order() {
    with_cx(|cx| {
        let original = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 30.0,
            },
            4.0,
        )
        .unwrap();
        let mesh = Mesh2::mapped_quads(2, 1, &|x, y| [x + 0.17 * x * y + 0.05 * y, y + 0.07 * x]);
        let nodal: Vec<_> = mesh
            .nodes
            .iter()
            .map(|[x, y]| [0.12 * x + 0.07 * x * y, -0.03 * y + 0.02 * x * y])
            .collect();
        let u: Vec<_> = nodal.iter().flatten().copied().collect();
        let points = material_point_kinematics(&mesh, &nodal, cx).unwrap();
        assert_eq!(points.len(), 8);
        near(
            points.iter().map(|p| p.reference_area_m2).sum(),
            1.07555,
            2e-13,
        );
        let elements: Vec<_> = (0..2)
            .map(|element| {
                (0..4)
                    .map(|point| {
                        let k = (4 * element + point) as f64;
                        let reference = ReferenceDeformation::try_new([
                            1.0 + 0.012 * k,
                            0.005 * k,
                            0.02 - 0.003 * k,
                            0.0,
                            0.97 + 0.007 * k,
                            -0.01 + 0.002 * k,
                            0.0,
                            0.0,
                            1.03 - 0.003 * k,
                        ])
                        .unwrap();
                        ReferencedHyperelastic::new(&original, reference)
                    })
                    .collect()
            })
            .collect();
        let field = PointField::new(elements);
        let pull = |x: f64, y: f64| [0.1 + 0.04 * x, -0.03 * y];
        let mut p = problem(&mesh, &field);
        p.traction = vec![(Patch::Right, &pull)];
        for (index, point) in points.iter().enumerate() {
            assert_eq!((point.element, point.local_point), (index / 4, index % 4));
            assert!(point.reference_area_m2 > 0.0);
        }
        let local = PlaneStressSettings::default();
        for form in POINT_FORMS {
            let mut expected_energy = 0.0;
            let mut expected_states = Vec::new();
            for point in &points {
                let law = &field.elements[point.element][point.local_point];
                let a = point.in_plane;
                match form {
                    Formulation::PlaneStrain => {
                        expected_energy += point.reference_area_m2
                            * law
                                .energy_density(&[a[0], a[1], 0.0, a[2], a[3], 0.0, 0.0, 0.0, 1.0]);
                    }
                    Formulation::PlaneStress => {
                        let state = plane_stress::plane_stress_state(law, &a, local).unwrap();
                        expected_energy += point.reference_area_m2 * state.energy_density_pa;
                        expected_states.push((state.deformation_gradient, state.piola));
                    }
                    Formulation::Incompressible => {
                        let state = incompressible_plane_stress::incompressible_plane_stress_state(
                            law, &a, local,
                        )
                        .unwrap();
                        expected_energy += point.reference_area_m2 * state.energy_density_pa;
                        expected_states.push((state.deformation_gradient, state.piola));
                    }
                    Formulation::MeanDilatation => unreachable!(),
                }
            }
            near(
                form.stored_energy(&p, &u, cx).unwrap(),
                expected_energy,
                3e-10,
            );
            let actual_states: Vec<_> = match form {
                Formulation::PlaneStrain => vec![],
                Formulation::PlaneStress => p
                    .plane_stress_states(&u, local, cx)
                    .unwrap()
                    .into_iter()
                    .map(|s| (s.deformation_gradient, s.piola))
                    .collect(),
                Formulation::Incompressible => p
                    .incompressible_plane_stress_states(&u, local, cx)
                    .unwrap()
                    .into_iter()
                    .map(|s| (s.deformation_gradient, s.piola))
                    .collect(),
                Formulation::MeanDilatation => unreachable!(),
            };
            for ((actual_f, actual_p), (expected_f, expected_p)) in
                actual_states.iter().zip(&expected_states)
            {
                for (a, b) in actual_f.iter().zip(expected_f) {
                    near(*a, *b, 3e-9);
                }
                for (a, b) in actual_p.iter().zip(expected_p) {
                    near(*a, *b, 3e-9);
                }
            }
            let (force, tangent) = form.probe(&p, &u, 0.0, cx).unwrap();
            let h = 2e-6;
            for column in [0, 3, 7, 10] {
                let mut plus = u.clone();
                let mut minus = u.clone();
                plus[column] += h;
                minus[column] -= h;
                near(
                    (form.stored_energy(&p, &plus, cx).unwrap()
                        - form.stored_energy(&p, &minus, cx).unwrap())
                        / (2.0 * h),
                    force[column],
                    8e-7,
                );
                let rp = cheap_force(form, &p, &plus, 0.0, cx).unwrap();
                let rm = cheap_force(form, &p, &minus, 0.0, cx).unwrap();
                for row in 0..u.len() {
                    near(
                        (rp[row] - rm[row]) / (2.0 * h),
                        tangent.get(row, column),
                        8e-6,
                    );
                }
            }
            let full = form.probe(&p, &u, 0.37, cx).unwrap().0;
            for (a, b) in cheap_force(form, &p, &u, 0.37, cx)
                .unwrap()
                .iter()
                .zip(full)
            {
                near(*a, b, 3e-11);
            }
        }
        assert!(matches!(
            p.solve_mean_dilatation(cx),
            Err(SolidError::InvalidInput { .. })
        ));
    });
}

#[test]
fn g1_shared_material_kinematics_integrates_original_area_and_affine_strain() {
    with_cx(|cx| {
        for (mesh, area) in [
            (Mesh2::triangles(1.3, 0.7, 2, 2), 0.91),
            (
                Mesh2::mapped_quads(2, 2, &|x, y| [1.3 * x + 0.2 * y, 0.7 * y + 0.1 * x]),
                0.89,
            ),
        ] {
            let u: Vec<_> = mesh
                .nodes
                .iter()
                .map(|[x, y]| [0.12 * x + 0.04 * y + 2.0, -0.03 * x - 0.03 * y - 1.0])
                .collect();
            let before = u.clone();
            let points = material_point_kinematics(&mesh, &u, cx).unwrap();
            let count = if mesh.elems[0].len() == 3 { 1 } else { 4 };
            assert_eq!(points.len(), mesh.elems.len() * count);
            near(
                points.iter().map(|p| p.reference_area_m2).sum(),
                area,
                3e-13,
            );
            for (index, point) in points.iter().enumerate() {
                assert_eq!(
                    (point.element, point.local_point),
                    (index / count, index % count)
                );
                assert!(point.reference_area_m2 > 0.0);
                for (a, b) in point.in_plane.iter().zip([1.12, 0.04, -0.03, 0.97]) {
                    near(*a, b, 4e-13);
                }
            }
            assert!(material_point_kinematics(&mesh, &u[..u.len() - 1], cx).is_err());
            let mut invalid = u.clone();
            invalid[0][1] = f64::INFINITY;
            assert!(material_point_kinematics(&mesh, &invalid, cx).is_err());
            let inverted: Vec<_> = mesh.nodes.iter().map(|[x, _]| [-2.0 * x, 0.0]).collect();
            assert!(matches!(
                material_point_kinematics(&mesh, &inverted, cx),
                Err(SolidError::MaterialRefused { .. })
            ));
            let mut bad_mesh = mesh.clone();
            bad_mesh.elems[0].swap(0, 1);
            assert!(matches!(
                material_point_kinematics(&bad_mesh, &u, cx),
                Err(SolidError::InvalidInput { .. })
            ));
            assert_eq!(u, before);
        }
    });
}

#[test]
fn g0_point_coverage_and_mean_dilatation_cannot_silently_drop_histories() {
    with_cx(|cx| {
        let mesh = Mesh2::quads(1.0, 0.7, 2, 1);
        let field = PointField::new(vec![vec![law(2.0); 4]; 2]);
        let u = vec![0.0; 2 * mesh.nodes.len()];
        let p = problem(&mesh, &field);
        for form in POINT_FORMS {
            assert!(form.probe(&p, &u, 0.0, cx).is_ok());
        }
        assert!(matches!(
            p.solve_mean_dilatation(cx),
            Err(SolidError::InvalidInput { .. })
        ));
        for defect in 0..5 {
            let mut invalid = field.clone();
            match defect {
                0 => {
                    invalid.elements[1].pop();
                }
                1 => invalid.elements[0].push(law(2.0)),
                2 => {
                    invalid.elements.pop();
                }
                3 => invalid.omitted = Some((1, 2)),
                4 => invalid.declare_elements = false,
                _ => unreachable!(),
            }
            let p = problem(&mesh, &invalid);
            for form in POINT_FORMS {
                assert!(
                    matches!(
                        form.probe(&p, &u, 0.0, cx),
                        Err(SolidError::InvalidInput { .. })
                    ),
                    "{form:?}, defect {defect}"
                );
            }
        }
    });
}

#[test]
fn g0_point_transverse_seeds_preserve_admitted_shear_and_exact_volume() {
    with_cx(|cx| {
        let original = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 40.0,
            },
            1.15,
        )
        .unwrap();
        let mesh = Mesh2::quads(1.0, 0.7, 1, 1);
        let seeds: Vec<_> = (0..4)
            .map(|q| [0.8 + 0.03 * q as f64, -0.7 + 0.02 * q as f64, 1.0])
            .collect();
        let laws: Vec<_> = seeds
            .iter()
            .map(|seed| {
                ReferencedHyperelastic::new(
                    &original,
                    ReferenceDeformation::try_new([
                        1.0, 0.0, seed[0], 0.0, 1.0, seed[1], 0.0, 0.0, 1.0,
                    ])
                    .unwrap(),
                )
            })
            .collect();
        let mut field = PointField::new(vec![laws]);
        let u = vec![0.0; 2 * mesh.nodes.len()];
        let local = PlaneStressSettings::default();
        assert!(
            problem(&mesh, &field)
                .plane_stress_states(&u, local, cx)
                .is_err()
        );
        assert!(
            problem(&mesh, &field)
                .incompressible_plane_stress_states(&u, local, cx)
                .is_err()
        );
        field.seeds = Some(vec![seeds.clone()]);
        for state in problem(&mesh, &field)
            .plane_stress_states(&u, local, cx)
            .unwrap()
        {
            near(state.energy_density_pa, 0.0, 2e-12);
            assert!(state.piola.iter().all(|p| p.abs() < 2e-10));
        }
        // The retained thickness is only a finite-solve seed. Exact elastic
        // volume must recompute thickness from the new in-plane deformation.
        for seed in &mut field.seeds.as_mut().unwrap()[0] {
            seed[2] = 9.0;
        }
        for (q, state) in problem(&mesh, &field)
            .incompressible_plane_stress_states(&u, local, cx)
            .unwrap()
            .iter()
            .enumerate()
        {
            near(state.deformation_gradient[2], seeds[q][0], 2e-12);
            near(state.deformation_gradient[5], seeds[q][1], 2e-12);
            near(state.deformation_gradient[8], 1.0, 2e-12);
            assert!(state.piola.iter().all(|p| p.abs() < 2e-10));
        }
        field.seeds.as_mut().unwrap()[0][2][2] = -1.0;
        assert!(matches!(
            problem(&mesh, &field).incompressible_plane_stress_states(&u, local, cx),
            Err(SolidError::InvalidInput { .. })
        ));
    });
}
