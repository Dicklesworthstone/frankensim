//! G1/G3: exact density moments and native graded-inventory acceleration.
use super::*;

const CORNERS: [[f64; 4]; 1] = [[3.0, 15.0, 3.0, 3.0]];
type Polynomial = Vec<([usize; 4], f64)>;

fn multiply(a: &Polynomial, b: &Polynomial) -> Polynomial {
    a.iter()
        .flat_map(|(p, c)| {
            b.iter()
                .map(move |(q, d)| (core::array::from_fn(|i| p[i] + q[i]), c * d))
        })
        .collect()
}

fn integral(poly: &Polynomial) -> f64 {
    let factorial = |n: usize| (1..=n).map(|i| i as f64).product::<f64>();
    poly.iter()
        .map(|(p, c)| {
            c * p.iter().map(|&i| factorial(i)).product::<f64>()
                / factorial(3 + p.iter().sum::<usize>())
        })
        .sum()
}

fn exact_mass(curvature: [f64; 3]) -> [[f64; 10]; 10] {
    let shapes: Vec<Polynomial> = (0..10)
        .map(|i| {
            let mut p = [0; 4];
            if i < 4 {
                p[i] = 1;
                let mut q = p;
                q[i] = 2;
                vec![(q, 2.0), (p, -1.0)]
            } else {
                let [a, b] = TET10_EDGES[i - 4];
                p[a] = 1;
                p[b] = 1;
                vec![(p, 4.0)]
            }
        })
        .collect();
    let density: Polynomial = CORNERS[0]
        .into_iter()
        .enumerate()
        .map(|(i, rho)| {
            let mut p = [0; 4];
            p[i] = 1;
            (p, rho)
        })
        .collect();
    // X_i = lambda_(i+1) + c_i*lambda_(i+1)^2. Its Jacobian is diagonal,
    // so exact barycentric monomial integration is independent of Duffy rules.
    let weight = curvature.iter().enumerate().fold(density, |poly, (i, c)| {
        let mut p = [0; 4];
        p[i + 1] = 1;
        multiply(&poly, &vec![([0; 4], 1.0), (p, 2.0 * c)])
    });
    core::array::from_fn(|a| {
        core::array::from_fn(|b| integral(&multiply(&multiply(&shapes[a], &shapes[b]), &weight)))
    })
}

#[test]
fn g1_graded_inertia_matches_exact_affine_and_curved_mass_moments() {
    with_cx(|cx| {
        let (affine, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        for curvature in [[0.0; 3], [0.4, 0.3, 0.2]] {
            let nodes: Vec<[f64; 3]> = affine
                .iter()
                .map(|x| core::array::from_fn(|i| x[i] + curvature[i] * x[i] * x[i]))
                .collect();
            let points = if curvature == [0.0; 3] { 4_usize } else { 6 };
            let original = MixedViscoTetSolid::try_new(
                &nodes,
                &elements,
                MixedViscoLawField::Uniform(&material),
                VolumetricTreatment::OriginalMaterial,
                TetAssemblyBudget::standard(),
                points,
                points.pow(3),
                cx,
            )
            .unwrap();
            let motion = MixedViscoTetMotion::try_new_with_reference_density(
                original,
                &CORNERS,
                &[[0.0; 3]; 10],
                MASS_BYTES,
                cx,
            )
            .unwrap();
            let mass = exact_mass(curvature);
            near(
                motion.reference_mass_kg(),
                mass.iter().flatten().sum(),
                3e-13,
            );
            assert_eq!(
                motion.reference_corner_density_kg_m3(),
                Some(CORNERS.as_slice())
            );
            for b in 0..10 {
                let mut acceleration = [[0.0; 3]; 10];
                acceleration[elements[0][b]][1] = 1.0;
                let loads = motion
                    .reference_acceleration_load_n(&acceleration, cx)
                    .unwrap();
                let mut full = [0.0; 30];
                for (dof, value) in loads {
                    full[dof] = value;
                }
                for a in 0..10 {
                    near(full[3 * elements[0][a] + 1], mass[a][b], 2e-13);
                    assert_eq!(full[3 * elements[0][a]], 0.0);
                    assert_eq!(full[3 * elements[0][a] + 2], 0.0);
                }
            }
            assert!(
                mass[0][4] < 0.0,
                "signed P2 coupling must not become lumped inertia"
            );
        }
    });
}

#[test]
fn g1_graded_acceleration_has_original_centroid_and_native_momentum_balance() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let acceleration = [0.3, -0.2, 0.1];
        let dt = 0.02;
        let initial = [0.2, -0.1, 0.05];
        // rho=3+12*x on the unit tetrahedron: mass=1, centroid=(13/40,9/40,9/40).
        let centroid = [13.0 / 40.0, 9.0 / 40.0, 9.0 / 40.0];
        for treatment in [
            VolumetricTreatment::OriginalMaterial,
            VolumetricTreatment::ExactIncompressible,
        ] {
            let motion = MixedViscoTetMotion::try_new_with_reference_density(
                solid(&nodes, &elements, &material, treatment, cx),
                &CORNERS,
                &vec![initial; nodes.len()],
                MASS_BYTES,
                cx,
            )
            .unwrap();
            near(motion.reference_mass_kg(), 1.0, 2e-14);
            near(motion.reference_density_kg_m3()[0], 6.0, 2e-13);
            let loads = motion
                .reference_acceleration_load_n(&vec![acceleration; nodes.len()], cx)
                .unwrap();
            assert!(
                loads
                    .iter()
                    .any(|&(dof, force)| dof % 3 == 0 && force < 0.0)
            );
            let (step, momentum) = motion
                .advance_with_momentum(
                    dt,
                    &[],
                    &loads,
                    [0.0; 3],
                    None,
                    settings(),
                    controls(),
                    [0.0; 3],
                    cx,
                )
                .unwrap();
            let after: [f64; 3] = core::array::from_fn(|c| initial[c] + dt * acceleration[c]);
            for velocity in step.motion.velocity_mps() {
                for c in 0..3 {
                    near(velocity[c], after[c], 2e-9);
                }
            }
            for c in 0..3 {
                near(
                    momentum.applied.linear_kg_m_s[c],
                    dt * acceleration[c],
                    2e-11,
                );
                near(momentum.linear_residual_kg_m_s[c], 0.0, 3e-10);
                near(momentum.angular_residual_kg_m2_s[c], 0.0, 3e-10);
            }
            let arm: [f64; 3] = core::array::from_fn(|c| centroid[c] + dt * after[c]);
            let torque = [
                arm[1] * acceleration[2] - arm[2] * acceleration[1],
                arm[2] * acceleration[0] - arm[0] * acceleration[2],
                arm[0] * acceleration[1] - arm[1] * acceleration[0],
            ];
            for c in 0..3 {
                near(momentum.applied.angular_kg_m2_s[c], dt * torque[c], 2e-11);
            }
            assert_eq!(
                step.motion.reference_corner_density_kg_m3(),
                Some(CORNERS.as_slice())
            );
            assert_eq!(step.motion.reference_mass_kg(), motion.reference_mass_kg());
            work_identity(&step);
        }
    });
}

#[test]
fn g3_uniform_profile_retains_scalar_bits_and_invalid_density_is_refused() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(10.0 * MU, true);
        let make = || {
            solid(
                &nodes,
                &elements,
                &material,
                VolumetricTreatment::OriginalMaterial,
                cx,
            )
        };
        let v = [[0.0; 3]; 10];
        let scalar = MixedViscoTetMotion::try_new(make(), &[DENSITY], &v, MASS_BYTES, cx).unwrap();
        let uniform = MixedViscoTetMotion::try_new_with_reference_density(
            make(),
            &[[DENSITY; 4]],
            &v,
            MASS_BYTES,
            cx,
        )
        .unwrap();
        let acceleration: Vec<[f64; 3]> = nodes
            .iter()
            .map(|x| [1.0 + x[0] * x[1], x[2], -0.3])
            .collect();
        assert_eq!(scalar.reference_mass_kg(), uniform.reference_mass_kg());
        assert_eq!(
            scalar.reference_density_kg_m3(),
            uniform.reference_density_kg_m3()
        );
        assert_eq!(
            scalar
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            uniform
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );
        for profile in [
            vec![],
            vec![[0.0, 3.0, 3.0, 3.0]],
            vec![[f64::NAN; 4]],
            vec![[-1.0; 4]],
        ] {
            assert!(
                MixedViscoTetMotion::try_new_with_reference_density(
                    make(),
                    &profile,
                    &v,
                    MASS_BYTES,
                    cx
                )
                .is_err()
            );
        }
        assert!(
            MixedViscoTetMotion::try_new_with_reference_density(make(), &CORNERS, &v, 1, cx)
                .is_err()
        );
        assert!(scalar.reference_acceleration_load_n(&[], cx).is_err());
        assert!(
            scalar
                .reference_acceleration_load_n(&[[f64::NAN; 3]; 10], cx)
                .is_err()
        );
        let mut curved = nodes.clone();
        curved[4][0] += 0.02;
        let coarse = solid(
            &curved,
            &elements,
            &material,
            VolumetricTreatment::OriginalMaterial,
            cx,
        );
        assert!(
            MixedViscoTetMotion::try_new_with_reference_density(
                coarse, &CORNERS, &v, MASS_BYTES, cx
            )
            .is_err()
        );
    });
}
