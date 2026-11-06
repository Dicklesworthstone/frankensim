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
fn g1_fixed_deformation_reference_shear_retains_work_beneath_large_bulk_storage() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e24);
        let density = [DENSITY];
        let deformation = [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let displacement = affine_displacement(&nodes, deformation);
        let velocity = vec![[0.0; 3]; nodes.len()];
        let definition = ConservativeTetDefinition {
            formulation: ConservativeFormulation::FullMaterial,
            ..definition(&nodes, &elements, &material, &density)
        };
        let accepted = definition.prepare(&displacement, &velocity, cx).unwrap();
        let original_bits = state_bits(&accepted);
        let shear = 2.0_f64.powi(-12);
        let reference =
            ReferenceDeformation::try_new([1.0, shear, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap();
        let changed = accepted
            .change_reference(
                MixedLawField::Uniform(&material),
                ConservativeReferenceField::Uniform(&reference),
                &[],
                1e-15,
                cx,
            )
            .unwrap();
        // Fe=diag(2,1,1)-2*shear e_x (x) e_y has unchanged J=2.
        // Both original Mooney-Rivlin invariants gain 4*shear^2 before
        // their respective J^(-2/3), J^(-4/3) factors.
        let expected = VOLUME
            * 4.0
            * shear
            * shear
            * (4.0 * 2.0_f64.powf(-2.0 / 3.0) + 2.0 * 2.0_f64.powf(-4.0 / 3.0));
        assert_eq!(
            changed.energy.stored_before_j,
            changed.energy.stored_after_j
        );
        near(changed.energy.stored_change_j, expected, 1e-16);
        near(changed.energy.reference_work_j, expected, 1e-16);
        near(changed.energy.bulk_stored_change_j, 0.0, 0.0);
        near(changed.energy.balance_error_j, 0.0, 0.0);
        near(changed.energy.cumulative_balance_error_j, 0.0, 0.0);
        assert_eq!(changed.solid.displacement_m(), accepted.displacement_m());
        assert_eq!(changed.solid.velocity_mps(), accepted.velocity_mps());
        assert_eq!(changed.solid.time_s(), accepted.time_s());
        assert_eq!(state_bits(&accepted), original_bits);
        let reversed = changed
            .solid
            .change_reference(
                MixedLawField::Uniform(&material),
                ConservativeReferenceField::Identity,
                &[],
                1e-15,
                cx,
            )
            .unwrap();
        near(reversed.energy.reference_work_j, -expected, 1e-16);
        near(reversed.energy.cumulative_balance_error_j, 0.0, 0.0);
        assert_eq!(state_bits(&reversed.solid), original_bits);
    });
}

#[test]
fn g1_bulk_reference_event_resolves_cancelling_modulus_work_without_resetting_volume() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let kappa = 2.0_f64.powi(80);
        let h = 2.0_f64.powi(-40);
        let next_kappa = kappa * (1.0 + 3.0 * h);
        let (material, target_material) = (law(kappa), law(next_kappa));
        let density = [DENSITY];
        let displacement =
            affine_displacement(&nodes, [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        let velocity = vec![[0.0; 3]; nodes.len()];
        let definition = ConservativeTetDefinition {
            formulation: ConservativeFormulation::OriginalMaterial {
                bulk: DiscreteBulkSettings {
                    minimum_volume_ratio: 0.8,
                    maximum_volume_ratio: 2.2,
                    maximum_iterations: 128,
                    absolute_pressure_tolerance_pa: 1e-12 * kappa,
                    relative_pressure_tolerance: 1e-12,
                    absolute_work_tolerance_pa: 1e-12,
                    relative_work_tolerance: 1e-12,
                },
                maximum_relative_volume_residual: 1e-9,
                maximum_point_volume_work_defect: 1e-12,
            },
            ..definition(&nodes, &elements, &material, &density)
        };
        let accepted = definition
            .prepare_with_pressure(&displacement, &velocity, &[kappa; 4], cx)
            .unwrap();
        let retained = accepted.bulk_volume_ratios().to_vec();
        let reference =
            ReferenceDeformation::try_new([1.0 + h, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
                .unwrap();
        let changed = accepted
            .change_reference(
                MixedLawField::Uniform(&target_material),
                ConservativeReferenceField::Uniform(&reference),
                &[],
                1e-12,
                cx,
            )
            .unwrap();
        // Integrate the exact rational source change at each retained q. The
        // pressure-space projection, including its original quadrature bits,
        // must not be replaced by the geometric value two for this audit.
        let expected = duffy_quadrature(4)
            .iter()
            .zip(&retained)
            .map(|((_, weight), q)| {
                let a = q - 1.0;
                weight * 0.5 * kappa * h * (2.0 * a * (q - 2.0) + h * (1.0 - 6.0 * a) + 3.0 * h * h)
                    / (1.0 + h)
            })
            .sum::<f64>();
        assert_eq!(
            changed.energy.bulk_stored_before_j,
            changed.energy.bulk_stored_after_j
        );
        near(changed.energy.bulk_stored_change_j, expected, 2e-12);
        near(changed.energy.bulk_reference_work_j, expected, 2e-12);
        assert_eq!(changed.solid.bulk_volume_ratios(), retained);
        assert_eq!(changed.solid.displacement_m(), accepted.displacement_m());
        near(changed.energy.cumulative_balance_error_j, 0.0, 0.0);
        let reversed = changed
            .solid
            .change_reference(
                MixedLawField::Uniform(&material),
                ConservativeReferenceField::Identity,
                &[],
                1e-12,
                cx,
            )
            .unwrap();
        near(reversed.energy.bulk_reference_work_j, -expected, 2e-12);
        near(reversed.energy.cumulative_balance_error_j, 0.0, 0.0);
        assert_eq!(reversed.solid.bulk_volume_ratios(), retained);
    });
}

#[test]
fn g1_material_stiffness_event_keeps_signed_work_beneath_unchanged_bulk_storage() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(1e24);
        let changed_material = Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 4.5,
                c01: 2.0,
                kappa: 1e24,
            },
            4.0,
        )
        .unwrap();
        let density = [DENSITY];
        let displacement =
            affine_displacement(&nodes, [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        let velocity = vec![[0.0; 3]; nodes.len()];
        let accepted = ConservativeTetDefinition {
            formulation: ConservativeFormulation::FullMaterial,
            ..definition(&nodes, &elements, &material, &density)
        }
        .prepare(&displacement, &velocity, cx)
        .unwrap();
        let changed = accepted
            .change_reference(
                MixedLawField::Uniform(&changed_material),
                ConservativeReferenceField::Identity,
                &[],
                1e-12,
                cx,
            )
            .unwrap();
        // Only c10 changes. At F=diag(2,1,1), I1=6 and J=2;
        // the original isochoric source work is Delta(c10)*(6 J^(-2/3)-3).
        let expected = VOLUME * 0.5 * (6.0 * 2.0_f64.powf(-2.0 / 3.0) - 3.0);
        assert_eq!(
            changed.energy.stored_after_j,
            changed.energy.stored_before_j
        );
        near(changed.energy.reference_work_j, expected, 2e-13);
        near(changed.energy.stored_change_j, expected, 2e-13);
        assert_eq!(changed.solid.displacement_m(), accepted.displacement_m());
        let held: Vec<_> = displacement
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .collect();
        let step = changed
            .solid
            .advance(0.01, &held, &[], [0.0; 3], settings(), 1e-12, cx)
            .unwrap();
        near(step.energy.stored_change_j, 0.0, 0.0);
        near(step.energy.constraint_reference_work_j, 0.0, 0.0);
        near(step.energy.cumulative_balance_error_j, 0.0, 1e-12);
        let reversed = step
            .solid
            .change_reference(
                MixedLawField::Uniform(&material),
                ConservativeReferenceField::Identity,
                &[],
                1e-12,
                cx,
            )
            .unwrap();
        near(reversed.energy.reference_work_j, -expected, 2e-13);
        near(reversed.energy.cumulative_balance_error_j, 0.0, 1e-12);
    });
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
    reference_energy_with_bulk(f, g, 40.0)
}

fn reference_energy_with_bulk(f: [f64; 9], g: [f64; 9], kappa: f64) -> f64 {
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
    jg * (4.0 * (iso * i1 - 3.0) + 2.0 * (iso * iso * i2 - 3.0) + 0.5 * kappa * (je - 1.0).powi(2))
}

fn reference_field_energy(element: &[usize; 10], u: &[[f64; 3]]) -> f64 {
    reference_field_storage(element, u, upper_reference, 40.0)
}

fn reference_field_storage(
    element: &[usize; 10],
    u: &[[f64; 3]],
    reference: impl Fn([f64; 4]) -> [f64; 9],
    kappa: f64,
) -> f64 {
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
            weight * reference_energy_with_bulk(f, reference(lambda), kappa)
        })
        .sum()
}

#[test]
fn g1_nonuniform_reference_event_keeps_motion_and_the_original_energy_ledger() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let (material, target_material) = (law(40.0), law(80.0));
        let densities = [DENSITY];
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let f0 = [1.04, 0.02, 0.0, 0.0, 0.98, 0.01, 0.0, 0.0, 1.03];
        let f1 = [1.05, 0.018, 0.003, 0.0, 0.985, 0.012, 0.0, 0.0, 1.025];
        let u0 = affine_displacement(&nodes, f0);
        let u1 = affine_displacement(&nodes, f1);
        let v0 = vec![[0.0; 3]; nodes.len()];
        let definition = ConservativeTetDefinition {
            formulation: ConservativeFormulation::FullMaterial,
            ..definition(&nodes, &elements, &material, &densities)
        };
        let initial = definition.prepare(&u0, &v0, cx).unwrap();
        let initial_energy = VOLUME * reference_energy(f0, identity);
        near(initial.total_energy_j(), initial_energy, 2e-13);
        // Begin with nonzero support work so a reset baseline cannot pass the
        // later accumulated-energy check accidentally.
        let prescribed: Vec<_> = u1.iter().flatten().copied().enumerate().collect();
        let loaded = initial
            .advance(
                0.04,
                &prescribed,
                &[],
                [0.0; 3],
                settings(),
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        assert!(loaded.energy.support_work_j.abs() > 1e-5);
        let accepted = loaded.solid;
        let before = state_bits(&accepted);
        let references: Vec<_> = duffy_quadrature(4)
            .iter()
            .map(|(lambda, _)| ReferenceDeformation::try_new(upper_reference(*lambda)).unwrap())
            .collect();
        let expected_before = VOLUME * reference_energy(f1, identity);
        let expected_after = reference_field_storage(&elements[0], &u1, upper_reference, 80.0);
        let changed = accepted
            .change_reference(
                MixedLawField::Uniform(&target_material),
                ConservativeReferenceField::PerQuadrature(&references),
                &[],
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        assert_eq!(state_bits(&accepted), before);
        assert_eq!(changed.solid.displacement_m(), accepted.displacement_m());
        assert_eq!(changed.solid.velocity_mps(), accepted.velocity_mps());
        assert_eq!(changed.solid.time_s(), accepted.time_s());
        assert_eq!(
            changed.solid.kinetic_energy_j(),
            accepted.kinetic_energy_j()
        );
        near(changed.energy.stored_before_j, expected_before, 2e-13);
        near(changed.energy.stored_after_j, expected_after, 2e-13);
        near(
            changed.energy.reference_work_j,
            expected_after - expected_before,
            2e-13,
        );
        near(changed.energy.bulk_reference_work_j, 0.0, 0.0);
        near(changed.energy.balance_error_j, 0.0, WORK_ALLOWANCE);
        near(
            changed.energy.cumulative_balance_error_j,
            0.0,
            WORK_ALLOWANCE,
        );
        near(changed.solid.reference_mass_kg(), DENSITY * VOLUME, 2e-14);
        let acceleration = vec![[0.2, -0.1, 0.4]; nodes.len()];
        assert_eq!(
            changed
                .solid
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            accepted
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );

        let continued = changed
            .solid
            .advance(0.005, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
            .unwrap();
        let stored = reference_field_storage(
            &elements[0],
            continued.solid.displacement_m(),
            upper_reference,
            80.0,
        );
        let (kinetic, _, _) = independent_motion(
            &nodes,
            &elements[0],
            continued.solid.displacement_m(),
            continued.solid.velocity_mps(),
        );
        near(continued.solid.stored_energy_j(), stored, 2e-13);
        near(continued.solid.kinetic_energy_j(), kinetic, 2e-13);
        near(
            stored + kinetic,
            initial_energy + loaded.energy.support_work_j + changed.energy.reference_work_j,
            WORK_ALLOWANCE,
        );
        near(continued.energy.constraint_reference_work_j, 0.0, 0.0);
        near(
            continued.energy.cumulative_balance_error_j,
            0.0,
            WORK_ALLOWANCE,
        );
    });
}

#[test]
fn g1_bulk_reference_event_preserves_retained_total_volumes_and_their_physical_storage() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let (material, target_material) = (law(40.0), law(80.0));
        let densities = [DENSITY];
        let (a, b) = (0.03, -0.02);
        let u: Vec<_> = nodes
            .iter()
            .map(|x| [a * x[0] * x[0], b * x[1] * x[1], 0.0])
            .collect();
        let velocity = vec![[0.1, -0.05, 0.02]; nodes.len()];
        // det(F)=1+2ax+2by+4abxy. Its exact L2 projection onto P1
        // uses xy -> (-1/30,2/15,2/15,-1/30) at the four vertices.
        // Retained q therefore differs from geometric det(F) pointwise.
        let projection = [-1.0 / 30.0, 2.0 / 15.0, 2.0 / 15.0, -1.0 / 30.0];
        let q_nodes: [f64; 4] = core::array::from_fn(|i| {
            1.0 + if i == 1 { 2.0 * a } else { 0.0 }
                + if i == 2 { 2.0 * b } else { 0.0 }
                + 4.0 * a * b * projection[i]
        });
        let pressure0 = q_nodes.map(|q| 40.0 * (q - 1.0));
        let definition = ConservativeTetDefinition {
            formulation: bulk_formulation(),
            ..definition(&nodes, &elements, &material, &densities)
        };
        let accepted = definition
            .prepare_with_pressure(&u, &velocity, &pressure0, cx)
            .unwrap();
        let before = state_bits(&accepted);
        let volumes = accepted.bulk_volume_ratios().to_vec();
        let quadrature = duffy_quadrature(4);
        let mut separation = 0.0_f64;
        for ((lambda, _), q) in quadrature.iter().zip(&volumes) {
            let analytic_q: f64 = (0..4).map(|i| lambda[i] * q_nodes[i]).sum();
            near(*q, analytic_q, 2e-14);
            let geometric = (1.0 + 2.0 * a * lambda[1]) * (1.0 + 2.0 * b * lambda[2]);
            separation = separation.max((q - geometric).abs());
        }
        assert!(separation > 1e-5);
        let g = [1.02, 0.01, 0.0, 0.0, 0.99, -0.02, 0.0, 0.0, 1.01];
        let reference = ReferenceDeformation::try_new(g).unwrap();
        let jg = g[0] * g[4] * g[8];
        let target_pressure = q_nodes.map(|q| 80.0 * (q / jg - 1.0));
        let bulk_energy = |values: &[f64], natural_volume: f64, kappa: f64| {
            quadrature
                .iter()
                .zip(values)
                .map(|((_, weight), q)| {
                    weight * natural_volume * 0.5 * kappa * (q / natural_volume - 1.0).powi(2)
                })
                .sum::<f64>()
        };
        let bulk_before = bulk_energy(&volumes, 1.0, 40.0);
        let bulk_after = bulk_energy(&volumes, jg, 80.0);
        let remainder_after = reference_field_storage(&elements[0], &u, |_| g, 0.0);
        let changed = accepted
            .change_reference(
                MixedLawField::Uniform(&target_material),
                ConservativeReferenceField::Uniform(&reference),
                &target_pressure,
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        assert_eq!(state_bits(&accepted), before);
        assert_eq!(changed.solid.bulk_volume_ratios(), &volumes);
        assert_eq!(changed.solid.displacement_m(), accepted.displacement_m());
        assert_eq!(changed.solid.velocity_mps(), accepted.velocity_mps());
        assert_eq!(changed.solid.algorithmic_pressure_pa(), target_pressure);
        near(changed.solid.reference_mass_kg(), DENSITY * VOLUME, 2e-14);
        near(changed.energy.bulk_stored_before_j, bulk_before, 2e-13);
        near(changed.energy.bulk_stored_after_j, bulk_after, 2e-13);
        near(
            changed.energy.bulk_reference_work_j,
            bulk_after - bulk_before,
            2e-13,
        );
        near(
            changed.energy.stored_after_j,
            remainder_after + bulk_after,
            2e-13,
        );
        near(
            changed.energy.reference_work_j,
            remainder_after + bulk_after - accepted.stored_energy_j(),
            2e-13,
        );
        let continued = changed
            .solid
            .advance(0.005, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
            .unwrap();
        let remainder =
            reference_field_storage(&elements[0], continued.solid.displacement_m(), |_| g, 0.0);
        let bulk = bulk_energy(continued.solid.bulk_volume_ratios(), jg, 80.0);
        let (kinetic, _, _) = independent_motion(
            &nodes,
            &elements[0],
            continued.solid.displacement_m(),
            continued.solid.velocity_mps(),
        );
        near(continued.solid.stored_energy_j(), remainder + bulk, 2e-13);
        near(
            remainder + bulk + kinetic,
            accepted.total_energy_j() + changed.energy.reference_work_j,
            WORK_ALLOWANCE,
        );
        near(continued.energy.constraint_reference_work_j, 0.0, 0.0);
        near(
            continued.energy.cumulative_balance_error_j,
            0.0,
            WORK_ALLOWANCE,
        );
    });
}

#[test]
fn g1_exact_natural_volume_change_supplies_signed_constraint_work_and_preserves_momenta() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(40.0);
        let densities = [DENSITY];
        let zero = vec![[0.0; 3]; nodes.len()];
        let accepted = definition(&nodes, &elements, &material, &densities)
            .prepare(&zero, &zero, cx)
            .unwrap();
        let before = state_bits(&accepted);
        for stretch in [1.01, 0.99] {
            let g = [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch];
            let reference = ReferenceDeformation::try_new(g).unwrap();
            let jg = stretch * stretch * stretch;
            // A zero-duration event cannot change the constrained total
            // volume at fixed F. The coupled interval must solve that motion.
            assert!(
                accepted
                    .change_reference(
                        MixedLawField::Uniform(&material),
                        ConservativeReferenceField::Uniform(&reference),
                        &[],
                        WORK_ALLOWANCE,
                        cx,
                    )
                    .is_err()
            );
            let changed = accepted
                .change_reference_and_advance(
                    MixedLawField::Uniform(&material),
                    ConservativeReferenceField::Uniform(&reference),
                    &[],
                    0.05,
                    &[],
                    &[],
                    [0.0; 3],
                    settings(),
                    WORK_ALLOWANCE,
                    cx,
                )
                .unwrap();
            let step = changed.step;
            let target_moment = VOLUME * (jg - 1.0) / 4.0;
            let actuation =
                -target_moment * step.solid.algorithmic_pressure_pa().iter().sum::<f64>();
            let (iso_storage, volume) =
                independent_fields(&elements[0], step.solid.displacement_m());
            let (kinetic, linear, angular) = independent_motion(
                &nodes,
                &elements[0],
                step.solid.displacement_m(),
                step.solid.velocity_mps(),
            );
            assert!(step.newton.is_some());
            assert!(kinetic > 1e-5);
            assert!(actuation > 1e-5);
            near(changed.reference.reference_work_j, 0.0, 2e-13);
            near(step.energy.constraint_reference_work_j, actuation, 2e-13);
            near(step.energy.pressure_work_j, -actuation, WORK_ALLOWANCE);
            near(step.solid.stored_energy_j(), jg * iso_storage, 2e-13);
            near(jg * iso_storage + kinetic, actuation, WORK_ALLOWANCE);
            near(step.inertia.numerical_loss_j, 0.0, 0.0);
            near(step.energy.balance_error_j, 0.0, WORK_ALLOWANCE);
            near(step.energy.cumulative_balance_error_j, 0.0, WORK_ALLOWANCE);
            near(step.solid.reference_mass_kg(), DENSITY * VOLUME, 2e-14);
            for moment in volume {
                near(moment, target_moment, 2e-11);
            }
            for c in 0..3 {
                near(linear[c], 0.0, 2e-11);
                near(angular[c], 0.0, 2e-11);
            }
            let continued = step
                .solid
                .advance(0.01, &[], &[], [0.0; 3], settings(), WORK_ALLOWANCE, cx)
                .unwrap();
            near(continued.energy.constraint_reference_work_j, 0.0, 0.0);
            near(continued.solid.total_energy_j(), actuation, WORK_ALLOWANCE);
            near(
                continued.energy.cumulative_balance_error_j,
                0.0,
                WORK_ALLOWANCE,
            );
            assert_eq!(state_bits(&accepted), before);
        }
    });
}

#[test]
fn g0_reference_event_admission_is_atomic_and_compatible_exact_changes_keep_inventory() {
    with_cx(|cx| {
        let (nodes, elements) = mesh(cx);
        let material = law(40.0);
        let densities = [DENSITY];
        let zero = vec![[0.0; 3]; nodes.len()];
        let velocity = vec![[0.1, -0.2, 0.3]; nodes.len()];
        let corners = [[2.0, 4.0, 6.0, 12.0]];
        let accepted = definition(&nodes, &elements, &material, &densities)
            .prepare_with_reference_density(&corners, &zero, &velocity, &[], cx)
            .unwrap();
        let before = state_bits(&accepted);
        let g = [1.1, 0.03, 0.0, 0.0, 1.0 / 1.1, -0.02, 0.0, 0.0, 1.0];
        let reference = ReferenceDeformation::try_new(g).unwrap();
        let changed = accepted
            .change_reference(
                MixedLawField::Uniform(&material),
                ConservativeReferenceField::Uniform(&reference),
                &[],
                WORK_ALLOWANCE,
                cx,
            )
            .unwrap();
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        near(
            changed.energy.reference_work_j,
            VOLUME * reference_energy_with_bulk(identity, g, 0.0),
            2e-13,
        );
        assert_eq!(
            changed.solid.reference_corner_density_kg_m3(),
            Some(corners.as_slice())
        );
        assert_eq!(
            changed.solid.reference_density_kg_m3(),
            accepted.reference_density_kg_m3()
        );
        assert_eq!(
            changed.solid.reference_mass_kg(),
            accepted.reference_mass_kg()
        );
        assert_eq!(
            changed.solid.kinetic_energy_j(),
            accepted.kinetic_energy_j()
        );
        assert_eq!(changed.solid.displacement_m(), accepted.displacement_m());
        assert_eq!(changed.solid.velocity_mps(), accepted.velocity_mps());
        assert_eq!(changed.solid.time_s(), accepted.time_s());
        let acceleration: Vec<_> = nodes
            .iter()
            .map(|x| [x[1], x[2] - x[0], 1.0 + x[0]])
            .collect();
        assert_eq!(
            changed
                .solid
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap(),
            accepted
                .reference_acceleration_load_n(&acceleration, cx)
                .unwrap()
        );
        assert!(
            accepted
                .change_reference(
                    MixedLawField::PerElement(&[]),
                    ConservativeReferenceField::Uniform(&reference),
                    &[],
                    WORK_ALLOWANCE,
                    cx,
                )
                .is_err()
        );
        assert!(
            accepted
                .change_reference(
                    MixedLawField::Uniform(&material),
                    ConservativeReferenceField::PerQuadrature(&[]),
                    &[],
                    WORK_ALLOWANCE,
                    cx,
                )
                .is_err()
        );
        for predictor in [&[0.0][..], &[0.0, 0.0, f64::NAN, 0.0][..]] {
            assert!(
                accepted
                    .change_reference(
                        MixedLawField::Uniform(&material),
                        ConservativeReferenceField::Uniform(&reference),
                        predictor,
                        WORK_ALLOWANCE,
                        cx,
                    )
                    .is_err()
            );
        }
        assert!(
            accepted
                .change_reference(
                    MixedLawField::Uniform(&material),
                    ConservativeReferenceField::Uniform(&reference),
                    &[],
                    f64::NAN,
                    cx,
                )
                .is_err()
        );
        let extreme =
            ReferenceDeformation::try_new([10.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 1.0]).unwrap();
        assert!(
            accepted
                .change_reference(
                    MixedLawField::Uniform(&material),
                    ConservativeReferenceField::Uniform(&extreme),
                    &[],
                    WORK_ALLOWANCE,
                    cx,
                )
                .is_err()
        );
        assert!(
            accepted
                .change_reference_and_advance(
                    MixedLawField::Uniform(&material),
                    ConservativeReferenceField::Uniform(&reference),
                    &[],
                    0.0,
                    &[],
                    &[],
                    [0.0; 3],
                    settings(),
                    WORK_ALLOWANCE,
                    cx,
                )
                .is_err()
        );
        assert_eq!(state_bits(&accepted), before);
    });
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
