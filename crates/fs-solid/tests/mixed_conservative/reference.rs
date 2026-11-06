//! Held natural strain in the existing conservative TET10 owner.
//! Independent invariant energies and exact reference mass distinguish total
//! deformation from elastic deformation without another time integrator.

use super::*;
use fs_material::hyper::discrete_gradient::bulk::DiscreteBulkSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_solid::hyper3::mixed::conservative::ConservativeReferenceField;

fn bulk_formulation() -> ConservativeFormulation {
    ConservativeFormulation::OriginalMaterial {
        bulk: DiscreteBulkSettings {
            minimum_volume_ratio: 0.8,
            maximum_volume_ratio: 1.2,
            maximum_iterations: 128,
            absolute_pressure_tolerance_pa: 1e-10,
            relative_pressure_tolerance: 1e-12,
            absolute_work_tolerance_pa: 1e-10,
            relative_work_tolerance: 1e-12,
        },
        maximum_relative_volume_residual: 1e-9,
        maximum_point_volume_work_defect: 1e-12,
    }
}

fn affine_displacement(nodes: &[[f64; 3]], f: [f64; 9]) -> Vec<[f64; 3]> {
    nodes
        .iter()
        .map(|x| core::array::from_fn(|i| (0..3).map(|j| f[3 * i + j] * x[j]).sum::<f64>() - x[i]))
        .collect()
}

#[test]
fn g1_large_reference_free_motion_preserves_natural_volume_and_original_mass() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(40.0);
        let densities = [DENSITY];
        let reference =
            ReferenceDeformation::try_new([5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0]).unwrap();
        let u0: Vec<_> = nodes.iter().map(|x| x.map(|v| 4.0 * v)).collect();
        let translation = [0.2, -0.1, 0.15];
        let v0 = vec![translation; nodes.len()];
        let base = definition(&nodes, &elements, &material, &densities);
        let (kinetic0, linear0, angular0) = independent_motion(&nodes, &elements[0], &u0, &v0);
        for formulation in [
            ConservativeFormulation::FullMaterial,
            base.formulation,
            bulk_formulation(),
        ] {
            let definition = ConservativeTetDefinition {
                formulation,
                ..base
            };
            // Total stretches five are outside the bare law's admitted domain;
            // elastic stretches are exactly one with the held natural map.
            assert!(definition.prepare(&u0, &v0, cx).is_err());
            let mut solid = definition
                .prepare_with_reference_deformations(
                    ConservativeReferenceField::Uniform(&reference),
                    &u0,
                    &v0,
                    &[],
                    cx,
                )
                .unwrap();
            near(solid.reference_mass_kg(), DENSITY * VOLUME, 2e-14);
            near(solid.stored_energy_j(), 0.0, WORK_ALLOWANCE);
            near(solid.kinetic_energy_j(), kinetic0, 2e-14);
            let before = state_bits(&solid);
            for _ in 0..2 {
                let step = solid
                    .advance(0.02, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                    .unwrap();
                near(step.inertia.numerical_loss_j, 0.0, 0.0);
                near(step.energy.material_work_j, 0.0, WORK_ALLOWANCE);
                near(step.energy.pressure_work_j, 0.0, WORK_ALLOWANCE);
                near(step.energy.balance_error_j, 0.0, WORK_ALLOWANCE);
                near(step.energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
                for residual in &step.pressure_residual_m3 {
                    near(*residual, 0.0, 2e-11);
                }
                if solid.time_s() == 0.0 {
                    assert_eq!(state_bits(&solid), before);
                }
                solid = step.solid;
                assert!(matches!(
                    solid.reference_deformations(),
                    ConservativeReferenceField::Uniform(retained)
                        if core::ptr::eq(retained, &reference)
                ));
                for ((u, v), initial) in solid
                    .displacement_m()
                    .iter()
                    .zip(solid.velocity_mps())
                    .zip(&u0)
                {
                    for c in 0..3 {
                        near(u[c], initial[c] + solid.time_s() * translation[c], 2e-9);
                        near(v[c], translation[c], 2e-9);
                    }
                }
                for volume in solid.bulk_volume_ratios() {
                    near(*volume, 125.0, 2e-9);
                }
                let (kinetic, linear, angular) = independent_motion(
                    &nodes,
                    &elements[0],
                    solid.displacement_m(),
                    solid.velocity_mps(),
                );
                near(kinetic, kinetic0, 2e-11);
                near(solid.reference_mass_kg(), DENSITY * VOLUME, 2e-14);
                for c in 0..3 {
                    near(linear[c], linear0[c], 2e-10);
                    near(angular[c], angular0[c], 2e-10);
                }
            }
        }
    });
}

fn upper_reference(lambda: [f64; 4]) -> [f64; 9] {
    [
        1.0 + 0.08 * lambda[0],
        0.03 * lambda[3],
        0.01 * lambda[1],
        0.0,
        1.0 - 0.06 * lambda[1],
        -0.02 * lambda[0],
        0.0,
        0.0,
        1.0 + 0.04 * lambda[2],
    ]
}

// Independent upper-triangular inverse and Mooney-Rivlin invariants. Neither
// the native reference pullback nor either discrete-gradient owner is used.
fn reference_energy(f: [f64; 9], g: [f64; 9]) -> f64 {
    let (a, b, c, s, t, u) = (g[0], g[4], g[8], g[1], g[5], g[2]);
    let jg = a * b * c;
    let h = [
        1.0 / a,
        -s / (a * b),
        (s * t - b * u) / jg,
        0.0,
        1.0 / b,
        -t / (b * c),
        0.0,
        0.0,
        1.0 / c,
    ];
    let fe: [f64; 9] = core::array::from_fn(|ij| {
        (0..3)
            .map(|k| f[3 * (ij / 3) + k] * h[3 * k + ij % 3])
            .sum()
    });
    let cof = [
        fe[4] * fe[8] - fe[5] * fe[7],
        fe[5] * fe[6] - fe[3] * fe[8],
        fe[3] * fe[7] - fe[4] * fe[6],
        fe[2] * fe[7] - fe[1] * fe[8],
        fe[0] * fe[8] - fe[2] * fe[6],
        fe[1] * fe[6] - fe[0] * fe[7],
        fe[1] * fe[5] - fe[2] * fe[4],
        fe[2] * fe[3] - fe[0] * fe[5],
        fe[0] * fe[4] - fe[1] * fe[3],
    ];
    let je = fe[0] * cof[0] + fe[1] * cof[1] + fe[2] * cof[2];
    let i1: f64 = fe.iter().map(|x| x * x).sum();
    let i2: f64 = cof.iter().map(|x| x * x).sum();
    let iso = je.powf(-2.0 / 3.0);
    jg * (4.0 * (iso * i1 - 3.0) + 2.0 * (iso * iso * i2 - 3.0) + 20.0 * (je - 1.0).powi(2))
}

fn reference_field_energy(element: &[usize; 10], u: &[[f64; 3]]) -> f64 {
    duffy_quadrature(4)
        .into_iter()
        .map(|(lambda, weight)| {
            // Analytic barycentric derivatives on the unit simplex, including
            // its quadratic edge modes after the initially affine body moves.
            let dl = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            let mut gradients = [[0.0; 3]; 10];
            for a in 0..4 {
                gradients[a] = dl[a].map(|v| (4.0 * lambda[a] - 1.0) * v);
            }
            for (edge, &[a, b]) in TET10_EDGES.iter().enumerate() {
                gradients[4 + edge] =
                    core::array::from_fn(|j| 4.0 * (lambda[a] * dl[b][j] + lambda[b] * dl[a][j]));
            }
            let f = core::array::from_fn(|ij| {
                let (i, j) = (ij / 3, ij % 3);
                (0..10).fold(if i == j { 1.0 } else { 0.0 }, |sum, a| {
                    sum + u[element[a]][i] * gradients[a][j]
                })
            });
            weight * reference_energy(f, upper_reference(lambda))
        })
        .sum()
}

#[test]
fn g1_nonuniform_reference_quadrature_work_matches_original_volume_invariant_energy() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(40.0);
        let densities = [DENSITY];
        let quadrature = duffy_quadrature(4);
        let references: Vec<_> = quadrature
            .iter()
            .map(|(lambda, _)| ReferenceDeformation::try_new(upper_reference(*lambda)).unwrap())
            .collect();
        let f0 = [1.08, 0.03, -0.01, 0.01, 0.96, 0.02, 0.0, -0.015, 1.04];
        let f1 = [
            1.083, 0.028, -0.009, 0.011, 0.958, 0.022, -0.001, -0.014, 1.042,
        ];
        let integrate = |f| {
            quadrature
                .iter()
                .map(|(lambda, weight)| weight * reference_energy(f, upper_reference(*lambda)))
                .sum::<f64>()
        };
        let (stored0, stored1) = (integrate(f0), integrate(f1));
        let u0 = affine_displacement(&nodes, f0);
        let u1 = affine_displacement(&nodes, f1);
        let zero = vec![[0.0; 3]; nodes.len()];
        let definition = ConservativeTetDefinition {
            formulation: ConservativeFormulation::FullMaterial,
            ..definition(&nodes, &elements, &material, &densities)
        };
        let accepted = definition
            .prepare_with_reference_deformations(
                ConservativeReferenceField::PerQuadrature(&references),
                &u0,
                &zero,
                &[],
                cx,
            )
            .unwrap();
        near(accepted.stored_energy_j(), stored0, 2e-13);
        let before = state_bits(&accepted);
        let prescribed: Vec<_> = u1.iter().flatten().copied().enumerate().collect();
        let step = accepted
            .advance(
                0.05,
                &prescribed,
                &[],
                [0.0; 3],
                settings(),
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        assert!(step.newton.is_none());
        assert_eq!(state_bits(&accepted), before);
        assert!(matches!(
            step.solid.reference_deformations(),
            ConservativeReferenceField::PerQuadrature(retained)
                if core::ptr::eq(retained.as_ptr(), references.as_ptr())
        ));
        near(step.energy.stored_before_j, stored0, 2e-13);
        near(step.energy.stored_after_j, stored1, 2e-13);
        near(step.solid.stored_energy_j(), stored1, 2e-13);
        near(step.energy.material_work_j, stored1 - stored0, 2e-13);
        near(step.energy.constitutive_work_error_j, 0.0, 2e-13);
        near(step.inertia.numerical_loss_j, 0.0, 0.0);
        let velocity: Vec<_> = u1
            .iter()
            .zip(&u0)
            .map(|(after, before)| core::array::from_fn(|c| 2.0 * (after[c] - before[c]) / 0.05))
            .collect();
        let (kinetic, _, _) = independent_motion(&nodes, &elements[0], &u1, &velocity);
        near(step.solid.kinetic_energy_j(), kinetic, 2e-14);
        near(
            step.energy.support_work_j,
            stored1 - stored0 + kinetic,
            2e-13,
        );
        near(step.energy.balance_error_j, 0.0, 2e-13);
        near(step.energy.cumulative_balance_error_j, 0.0, 2e-13);
        near(step.solid.reference_mass_kg(), DENSITY * VOLUME, 2e-14);

        // A separate released branch must deform through actual Newton/FGMRES
        // under the heterogeneous residual stress. Its total physical energy
        // and momenta are checked independently of the native work ledger.
        let free = accepted
            .advance(0.005, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
            .unwrap();
        assert!(free.newton.is_some());
        let stored = reference_field_energy(&elements[0], free.solid.displacement_m());
        let (kinetic, linear, angular) = independent_motion(
            &nodes,
            &elements[0],
            free.solid.displacement_m(),
            free.solid.velocity_mps(),
        );
        assert!(kinetic > 1e-8);
        assert!(stored < stored0 - 1e-8);
        near(free.solid.stored_energy_j(), stored, 2e-13);
        near(free.solid.kinetic_energy_j(), kinetic, 2e-14);
        near(stored + kinetic, stored0, WORK_ALLOWANCE);
        near(free.inertia.numerical_loss_j, 0.0, 0.0);
        near(free.energy.balance_error_j, 0.0, WORK_ALLOWANCE);
        for c in 0..3 {
            near(linear[c], 0.0, 2e-11);
            near(angular[c], 0.0, 2e-11);
        }
        assert_eq!(state_bits(&accepted), before);
    });
}

#[test]
fn g1_reference_coverage_identity_parity_and_graded_mass_compose_transactionally() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(40.0);
        let densities = [DENSITY];
        let base = definition(&nodes, &elements, &material, &densities);
        let identity = ReferenceDeformation::identity();
        let (u0, v0) = initial_state(&nodes);
        for formulation in [ConservativeFormulation::FullMaterial, base.formulation] {
            let definition = ConservativeTetDefinition {
                formulation,
                ..base
            };
            let legacy = definition.prepare(&u0, &v0, cx).unwrap();
            let explicit = definition
                .prepare_with_reference_deformations(
                    ConservativeReferenceField::Uniform(&identity),
                    &u0,
                    &v0,
                    &[],
                    cx,
                )
                .unwrap();
            assert_eq!(state_bits(&legacy), state_bits(&explicit));
            let legacy_step = legacy
                .advance(0.005, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            let explicit_step = explicit
                .advance(0.005, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            assert_eq!(
                state_bits(&legacy_step.solid),
                state_bits(&explicit_step.solid)
            );
        }
        let zero = vec![[0.0; 3]; nodes.len()];
        for count in [0, 2] {
            let references = vec![identity; count];
            assert!(
                base.prepare_with_reference_deformations(
                    ConservativeReferenceField::PerElement(&references),
                    &zero,
                    &zero,
                    &[],
                    cx,
                )
                .is_err()
            );
        }
        for count in [0, 63, 65] {
            let references = vec![identity; count];
            assert!(
                base.prepare_with_reference_deformations(
                    ConservativeReferenceField::PerQuadrature(&references),
                    &zero,
                    &zero,
                    &[],
                    cx,
                )
                .is_err()
            );
        }

        let reference =
            ReferenceDeformation::try_new([5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0]).unwrap();
        let references = [reference];
        let full = ConservativeTetDefinition {
            formulation: ConservativeFormulation::FullMaterial,
            ..base
        };
        assert!(
            full.prepare_with_reference_deformations(
                ConservativeReferenceField::PerElement(&references),
                &zero,
                &zero,
                &[],
                cx,
            )
            .is_err()
        );
        let u: Vec<_> = nodes.iter().map(|x| x.map(|v| 4.0 * v)).collect();
        let velocity = vec![[0.2, -0.1, 0.3]; nodes.len()];
        let corners = [[2.0, 4.0, 6.0, 12.0]];
        let accepted = full
            .prepare_with_reference_deformations_and_density(
                ConservativeReferenceField::PerElement(&references),
                &corners,
                &u,
                &velocity,
                &[],
                cx,
            )
            .unwrap();
        let original_mass = VOLUME * corners[0].iter().sum::<f64>() / 4.0;
        near(accepted.reference_mass_kg(), original_mass, 2e-14);
        near(
            accepted.kinetic_energy_j(),
            0.5 * original_mass * 0.14,
            2e-14,
        );
        assert!(matches!(
            accepted.reference_deformations(),
            ConservativeReferenceField::PerElement(retained)
                if core::ptr::eq(retained.as_ptr(), references.as_ptr())
        ));
        let before = state_bits(&accepted);
        let inverted = affine_displacement(&nodes, [-5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0]);
        let prescribed: Vec<_> = inverted.iter().flatten().copied().enumerate().collect();
        assert!(
            accepted
                .advance(
                    0.05,
                    &prescribed,
                    &[],
                    [0.0; 3],
                    settings(),
                    WORK_ALLOWANCE,
                    cx
                )
                .is_err()
        );
        assert_eq!(state_bits(&accepted), before);
        assert_eq!(accepted.displacement_m(), &u);
        assert_eq!(accepted.velocity_mps(), &velocity);
    });
}
