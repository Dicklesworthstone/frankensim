//! G1/G3/G4: exact density moments, supported mobility and retained hot history.
use super::*;
use fs_solid::hyper3::plastic::hex::dynamics::{HexImpulsePoint, SolidImpulseControls};
use fs_solver::CheckedCgConfig;

fn impulse_controls() -> SolidImpulseControls {
    SolidImpulseControls {
        solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        maximum_workspace_components: 4096,
        linear_tolerance_kg_m_s: 1e-11,
        angular_tolerance_kg_m2_s: 1e-11,
        work_tolerance_j: 1e-11,
    }
}

fn graded<'a>(
    nodes: &'a [[f64; 3]],
    rho: [f64; 8],
    velocity: [[f64; 3]; 8],
    cx: &Cx<'_>,
) -> FiniteThermoplasticHexMotion<'a> {
    let solid = FiniteThermoplasticHexSolid::try_new(
        nodes,
        &HEX,
        vec![family()],
        &[300.0],
        FiniteHexFormulation::Standard,
        HexAssemblyBudget {
            maximum_tangent_entries: 768,
            ..Default::default()
        },
        cx,
    )
    .unwrap();
    FiniteThermoplasticHexMotion::try_new_with_reference_density(solid, &[rho], &velocity, cx)
        .unwrap()
}

// Integral of an explicitly expanded one-dimensional polynomial on [0,1].
// No production shape, quadrature or Jacobian routine enters these oracles.
fn integral(factors: &[&[f64]]) -> f64 {
    let mut p = vec![1.0];
    for factor in factors {
        let mut next = vec![0.0; p.len() + factor.len() - 1];
        for (i, a) in p.iter().enumerate() {
            for (j, b) in factor.iter().enumerate() {
                next[i + j] += a * b;
            }
        }
        p = next;
    }
    p.iter().enumerate().map(|(i, c)| c / (i + 1) as f64).sum()
}

#[test]
fn graded_mass_preserves_exact_affine_and_warped_inventory_and_kinetic_moments() {
    with_cx(|cx| {
        // x=s(1+a*t), y=t, z=u(1+b*t), J=(1+a*t)(1+b*t).
        // The quadratic Jacobian makes a two-point mass rule insufficient.
        for (a, b) in [(0.0, 0.0), (0.5, 0.25)] {
            let nodes = NODES.map(|[s, t, u]| [s * (1.0 + a * t), t, u * (1.0 + b * t)]);
            let rho = NODES.map(|[_, t, _]| 2.0 + 4.0 * t);
            let r = [2.0, 4.0];
            let x = [1.0, a];
            let z = [1.0, b];
            let mass = integral(&[&r, &x, &z]);
            let first = [
                0.5 * integral(&[&r, &x, &x, &z]),
                integral(&[&r, &x, &z, &[0.0, 1.0]]),
                0.5 * integral(&[&r, &x, &z, &z]),
            ];
            let moving = graded(&nodes, rho, [[0.0, 1.0, 0.0]; 8], cx);
            near(moving.reference_mass_kg(), mass, 3e-14, 2e-14);
            near(
                moving.reference_density_kg_m3()[0],
                mass / integral(&[&x, &z]),
                3e-14,
                2e-14,
            );
            let report = moving
                .impulse_candidate(&[], &[], [0.0; 3], impulse_controls(), cx)
                .unwrap()
                .impulse;
            near(report.before.linear_kg_m_s[1], mass, 3e-14, 2e-14);
            near(report.before.angular_kg_m2_s[0], -first[2], 3e-14, 2e-14);
            near(report.before.angular_kg_m2_s[2], first[0], 3e-14, 2e-14);
            near(report.kinetic_before_j, 0.5 * mass, 3e-14, 2e-14);

            let moving = graded(&nodes, rho, NODES.map(|[_, t, _]| [0.0, t, 0.0]), cx);
            let report = moving
                .impulse_candidate(&[], &[], [0.0; 3], impulse_controls(), cx)
                .unwrap()
                .impulse;
            near(report.before.linear_kg_m_s[1], first[1], 3e-14, 2e-14);
            near(
                report.kinetic_before_j,
                0.5 * integral(&[&r, &x, &z, &[0.0, 0.0, 1.0]]),
                3e-14,
                2e-14,
            );
            near(
                report.before.angular_kg_m2_s[2],
                0.5 * integral(&[&r, &x, &x, &z, &[0.0, 1.0]]),
                3e-14,
                2e-14,
            );
        }
    });
}

#[test]
fn equal_mass_density_redistribution_changes_native_point_mobility() {
    with_cx(|cx| {
        let held: Vec<_> = (0..24).filter(|dof| *dof != 3).collect();
        let point = HexImpulsePoint {
            hexahedron: 0,
            reference_coordinates: [1.0, -1.0, -1.0],
        };
        // M_11 = integral rho*x²*(1-y)²*(1-z)². Both profiles have mass 5,
        // but rho=2+6x gives 13/54 and rho=8-6x gives 7/54 kg.
        for (left, slope, diagonal) in [(2.0, 6.0, 13.0 / 54.0), (8.0, -6.0, 7.0 / 54.0)] {
            let rho = NODES.map(|[x, _, _]| left + slope * x);
            let motion = graded(&NODES, rho, [[0.0; 3]; 8], cx);
            near(motion.reference_mass_kg(), 5.0, 3e-14, 0.0);
            let mobility = motion
                .point_directional_inverse_mass_per_kg(
                    &held,
                    point,
                    [1.0, 0.0, 0.0],
                    impulse_controls(),
                    cx,
                )
                .unwrap();
            near(mobility.inverse_mass_per_kg, 1.0 / diagonal, 2e-12, 2e-13);
            let next = motion
                .point_impulse_candidate(
                    &held,
                    point,
                    [0.01, 0.0, 0.0],
                    [0.0; 3],
                    impulse_controls(),
                    cx,
                )
                .unwrap();
            near(
                next.motion.velocity_mps()[1][0],
                0.01 / diagonal,
                2e-13,
                2e-13,
            );
            assert_eq!(
                next.motion.reference_corner_density_kg_m3(),
                Some([rho].as_slice())
            );
            assert_eq!(motion.velocity_mps(), [[0.0; 3]; 8]);
        }
    });
}

#[test]
fn constant_profile_reproduces_exact_scalar_mass_on_affine_and_warped_bricks() {
    with_cx(|cx| {
        for (a, b) in [(0.0, 0.0), (0.5, 0.25)] {
            let nodes = NODES.map(|[s, t, u]| [s * (1.0 + a * t), t, u * (1.0 + b * t)]);
            let velocity = nodes.map(|[x, _, _]| [x, 0.0, 0.0]);
            let make = || {
                FiniteThermoplasticHexSolid::try_new(
                    &nodes,
                    &HEX,
                    vec![family()],
                    &[300.0],
                    FiniteHexFormulation::Standard,
                    HexAssemblyBudget {
                        maximum_tangent_entries: 768,
                        ..Default::default()
                    },
                    cx,
                )
                .unwrap()
            };
            let scalar =
                FiniteThermoplasticHexMotion::try_new(make(), &[4.0], &velocity, cx).unwrap();
            let exact = FiniteThermoplasticHexMotion::try_new_with_reference_density(
                make(),
                &[[4.0; 8]],
                &velocity,
                cx,
            )
            .unwrap();
            let measure = |motion: &FiniteThermoplasticHexMotion<'_>| {
                motion
                    .impulse_candidate(&[], &[], [0.0; 3], impulse_controls(), cx)
                    .unwrap()
                    .impulse
                    .kinetic_before_j
            };
            let energy = (2.0 / 3.0) * integral(&[&[1.0, a], &[1.0, a], &[1.0, a], &[1.0, b]]);
            near(measure(&exact), energy, 3e-14, 2e-14);
            assert_eq!(measure(&scalar).to_bits(), measure(&exact).to_bits());
            assert_eq!(
                scalar.reference_mass_kg().to_bits(),
                exact.reference_mass_kg().to_bits()
            );
        }
    });
}

#[test]
fn graded_inventory_and_hot_plastic_history_survive_refusal_impulse_and_motion() {
    with_cx(|cx| {
        let fixed: Vec<_> = NODES
            .iter()
            .flat_map(|x| {
                (0..3).map(move |i| {
                    let logarithm: f64 = if i == 0 { 0.16 } else { -0.08 };
                    (logarithm.exp() - 1.0) * x[i]
                })
            })
            .enumerate()
            .collect();
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            let solid = initial(family(), formulation, cx)
                .advance(
                    &[350.0],
                    &fixed,
                    &[],
                    0.2,
                    settings(),
                    controls().material,
                    cx,
                )
                .unwrap()
                .solid;
            let states = solid.quadrature_states().to_vec();
            assert!(states.iter().all(|state| state.alpha > 0.0));
            let rho = NODES.map(|[x, _, _]| 2.0 + 6.0 * x);
            let motion = FiniteThermoplasticHexMotion::try_new_with_reference_density(
                solid,
                &[rho],
                &[[0.1, 0.0, 0.0]; 8],
                cx,
            )
            .unwrap();
            let report = motion
                .impulse_candidate(&[], &[], [0.0; 3], impulse_controls(), cx)
                .unwrap();
            assert_eq!(report.motion.solid().quadrature_states(), states);
            assert_eq!(
                report.motion.reference_corner_density_kg_m3(),
                Some([rho].as_slice())
            );
            let mut refusal = controls();
            refusal.maximum_numerical_kinetic_loss_j = 0.0;
            assert!(
                report
                    .motion
                    .advance(&[400.0], &fixed, &[], 0.1, settings(), refusal, cx)
                    .is_err()
            );
            assert_eq!(report.motion.solid().temperatures_k(), [350.0]);
            assert_eq!(report.motion.solid().quadrature_states(), states);
            let step = report
                .motion
                .advance(&[400.0], &fixed, &[], 0.1, settings(), controls(), cx)
                .unwrap();
            assert_eq!(step.motion.solid().quadrature_states(), states);
            assert_eq!(
                step.motion.reference_corner_density_kg_m3(),
                Some([rho].as_slice())
            );
            near(step.motion.reference_mass_kg(), 5.0, 3e-14, 0.0);
            near(step.inertia.numerical_loss_j, 0.025, 3e-14, 2e-13);
        }
        for rho in [vec![], vec![[0.0; 8]], vec![[f64::NAN; 8]]] {
            assert!(
                FiniteThermoplasticHexMotion::try_new_with_reference_density(
                    initial(family(), FiniteHexFormulation::Standard, cx),
                    &rho,
                    &[[0.0; 3]; 8],
                    cx
                )
                .is_err()
            );
        }
    });
}
