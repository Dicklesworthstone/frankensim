//! Held-control shape work from the original tetrahedral constitutive energy.
use super::*;
use fs_em::electrostatic::nonlinear::spatial::SpatialNonlinearElectricField;

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn volume(positions: &[[f64; 3]], cell: [u32; 4]) -> f64 {
    let p = cell.map(|v| positions[v as usize]);
    let a: [f64; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
    let b: [f64; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
    let c: [f64; 3] = std::array::from_fn(|i| p[3][i] - p[0][i]);
    dot3(
        a,
        [
            b[1] * c[2] - b[2] * c[1],
            b[2] * c[0] - b[0] * c[2],
            b[0] * c[1] - b[1] * c[0],
        ],
    )
    .abs()
        / 6.0
}

fn mixed_functional(field: &SpatialNonlinearElectricField<'_>) -> f64 {
    let problem = field.problem();
    let source_pairing: f64 = problem
        .mesh()
        .tets
        .iter()
        .zip(field.free_charge_density_c_m3())
        .map(|(&cell, &rho)| {
            let phi: f64 = cell
                .iter()
                .map(|&v| field.relative_voltage_v()[v as usize])
                .sum();
            0.25 * rho * volume(problem.positions_m(), cell) * phi
        })
        .sum();
    let charge_pairing: f64 = field
        .electrode_controls()
        .iter()
        .zip(field.terminal_voltage_v())
        .filter_map(|(&control, &voltage)| match control {
            ElectrodeControl::Charge(q) => Some(q * voltage),
            ElectrodeControl::Voltage(_) => None,
        })
        .sum();
    field.coenergy_j() - source_pairing - charge_pairing
}

/// Independent finite-geometry oracle. Each cell's actual original q stays
/// fixed: the new density is q/new_volume, not the original density.
fn perturbed_functional(
    field: &SpatialNonlinearElectricField<'_>,
    displacement: &[[f64; 3]],
    step: f64,
) -> f64 {
    let original = field.problem();
    let positions: Vec<_> = original
        .positions_m()
        .iter()
        .zip(displacement)
        .map(|(&x, &u)| std::array::from_fn(|i| step.mul_add(u[i], x[i])))
        .collect();
    let cells = original.mesh().tets.clone();
    let rho: Vec<_> = cells
        .iter()
        .zip(field.free_charge_density_c_m3())
        .map(|(&cell, &density)| {
            density * volume(original.positions_m(), cell) / volume(&positions, cell)
        })
        .collect();
    let problem = SpatialNonlinearElectrostatics::new(
        positions,
        cells,
        original.materials().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let terminals: Vec<_> = field.terminal_faces().iter().map(Vec::as_slice).collect();
    let electrodes = problem
        .prepare_electrodes(&terminals, &mut || false)
        .unwrap();
    let result = electrodes
        .solve(
            &rho,
            field.electrode_controls(),
            field.relative_voltage_v(),
            controls(),
            &mut || false,
        )
        .unwrap();
    near(result.volume_charge_c(), field.volume_charge_c(), Q0);
    mixed_functional(&result)
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    let axes = [
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
        [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
    ];
    std::array::from_fn(|i| v[0].mul_add(axes[0][i], v[1].mul_add(axes[1][i], v[2] * axes[2][i])))
}

#[test]
fn nonlinear_shape_work_uses_original_storage_and_coenergy_under_voltage_and_charge_controls() {
    let (original, cells) = box_mesh(&[0.0, 0.5 * GAP, GAP], WIDTH, DEPTH, 1, 1);
    let laws = vec![SpatialDielectricMaterial::Reversible(nonlinear()); cells.len()];
    let positions = original.iter().copied().map(rotate).collect();
    let problem =
        SpatialNonlinearElectrostatics::new(positions, cells, laws, limits(), &mut || false)
            .unwrap();
    let terminals = [
        faces_at(&problem, &original, 0.0),
        faces_at(&problem, &original, GAP),
    ];
    let electrodes = problem
        .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
        .unwrap();
    let initial: Vec<_> = original.iter().map(|p| 2.0 * E0 * p[0]).collect();
    let rho = vec![0.0; problem.mesh().tets.len()];
    for right in [
        ElectrodeControl::Voltage(2.0 * V0),
        ElectrodeControl::Charge(2.5 * Q0),
    ] {
        let field = electrodes
            .solve(
                &rho,
                &[ElectrodeControl::Voltage(0.0), right],
                &initial,
                controls(),
                &mut || false,
            )
            .unwrap();
        near(field.energy_j(), 1.75 * U0, U0);
        near(field.coenergy_j(), 3.25 * U0, U0);
        let forces = field.shape_forces(&mut || false).unwrap();
        assert!(std::ptr::eq(forces.field(), &field));
        for axis in 0..3 {
            let u: Vec<_> = original
                .iter()
                .map(|p| {
                    let mut deformation = [0.0; 3];
                    deformation[axis] = p[axis];
                    rotate(deformation)
                })
                .collect();
            let expected = if axis == 0 { -1.75 * U0 } else { 3.25 * U0 };
            let work = field.virtual_work(&u, &mut || false).unwrap();
            near(work, expected, U0);
            let nodal: f64 = forces
                .nodal_force_n()
                .iter()
                .zip(&u)
                .map(|(&f, &u)| dot3(f, u))
                .sum();
            near(nodal, expected, U0);
            if axis == 0 {
                let h = 1e-5;
                let derivative = (perturbed_functional(&field, &u, h)
                    - perturbed_functional(&field, &u, -h))
                    / (2.0 * h);
                near(derivative, work, U0);
            }
        }
        // The cell-relative contraction must cancel common translations before
        // multiplying them by the force, even at a huge displacement origin.
        assert_eq!(
            field
                .virtual_work(&vec![[1e150, -2e150, 3e150]; original.len()], &mut || false)
                .unwrap(),
            0.0
        );
        let rotation: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| [-p[1], p[0], 0.0])
            .collect();
        near(
            field.virtual_work(&rotation, &mut || false).unwrap(),
            0.0,
            U0,
        );
    }
}

#[test]
fn materially_transported_charge_in_one_tetrahedron_has_the_original_shape_derivative() {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [GAP, 0.0, 0.0],
        [0.0, WIDTH, 0.0],
        [0.0, 0.0, DEPTH],
    ];
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        vec![[0, 1, 2, 3]],
        vec![SpatialDielectricMaterial::Reversible(nonlinear())],
        limits(),
        &mut || false,
    )
    .unwrap();
    let ground = faces_at(&problem, &positions, 0.0);
    let electrodes = problem
        .prepare_electrodes(&[&ground], &mut || false)
        .unwrap();
    let field = electrodes
        .solve(
            &[10.0 * D0 / GAP],
            &[ElectrodeControl::Voltage(0.0)],
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(field.energy_j(), 1.75 * U0 / 6.0, U0);
    near(mixed_functional(&field), -field.energy_j(), U0);
    let u: Vec<_> = positions
        .iter()
        .map(|p| {
            [
                0.2 * p[0] + 0.03 * p[1],
                -0.12 * p[1] + 0.07 * p[2],
                0.1 * p[2],
            ]
        })
        .collect();
    let work = field.virtual_work(&u, &mut || false).unwrap();
    // The axial integrated stress is U; both transverse stresses are -W*.
    near(work, (-0.2 * 1.75 - 0.02 * 3.25) * U0 / 6.0, U0);
    let forces = field.shape_forces(&mut || false).unwrap();
    let nodal: f64 = forces
        .nodal_force_n()
        .iter()
        .zip(&u)
        .map(|(&f, &u)| dot3(f, u))
        .sum();
    near(nodal, work, U0);
    let h = 1e-5;
    let derivative =
        (perturbed_functional(&field, &u, h) - perturbed_functional(&field, &u, -h)) / (2.0 * h);
    near(derivative, work, U0);
}

#[test]
fn charged_mixed_control_shape_work_keeps_cell_charge_and_source_voltage_pairings() {
    let (positions, cells) = box_mesh(&[0.0, 0.5 * GAP, GAP], WIDTH, DEPTH, 1, 1);
    let laws = vec![SpatialDielectricMaterial::Reversible(nonlinear()); cells.len()];
    let problem =
        SpatialNonlinearElectrostatics::new(positions.clone(), cells, laws, limits(), &mut || {
            false
        })
        .unwrap();
    let terminals = [
        faces_at(&problem, &positions, 0.0),
        faces_at(&problem, &positions, GAP),
    ];
    let electrodes = problem
        .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
        .unwrap();
    let rho = vec![0.2 * D0 / GAP; problem.mesh().tets.len()];
    let initial: Vec<_> = positions.iter().map(|p| 2.0 * E0 * p[0]).collect();
    let field = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(7.0 * V0),
                ElectrodeControl::Charge(2.5 * Q0),
            ],
            &initial,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(field.volume_charge_c() > 0.0);
    assert!(field.terminal_voltage_v()[1] * field.terminal_charge_c()[1] > U0);
    near(mixed_functional(&field), -field.energy_j(), U0);
    // Interior bending changes individual cell volumes differently. The oracle
    // must transport each cell's q, not merely preserve the global total.
    let u: Vec<_> = positions
        .iter()
        .map(|p| {
            [
                0.13 * p[0] + 0.015 * p[1],
                -0.07 * p[1] + 0.02 * p[2],
                0.05 * p[2] + 0.03 * p[0] * (GAP - p[0]) / GAP,
            ]
        })
        .collect();
    let work = field.virtual_work(&u, &mut || false).unwrap();
    let h = 1e-5;
    let derivative =
        (perturbed_functional(&field, &u, h) - perturbed_functional(&field, &u, -h)) / (2.0 * h);
    near(derivative, work, U0);
    let forces = field.shape_forces(&mut || false).unwrap();
    let nodal: f64 = forces
        .nodal_force_n()
        .iter()
        .zip(&u)
        .map(|(&f, &u)| dot3(f, u))
        .sum();
    near(nodal, work, U0);
}

fn almost_flat_source() -> ReversibleDielectricMaterial {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
                ELECTRIC_DISPLACEMENT_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "E".into(),
                abscissa_dims: ELECTRIC_FIELD_DIMS,
                knots: vec![(0.0, 0.0), (1e-6, D0), (4.0 * E0, D0 + 4e-18)],
                dims: ELECTRIC_DISPLACEMENT_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "E",
                QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
                0.0,
                4.0 * E0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored small-storage stress-cancellation fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic steep then shallow dielectric".into(),
            phase: "solid".into(),
            process: "3D shape test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let context = QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    ReversibleDielectricMaterial::resolve(
        &card,
        &context,
        "E",
        4.0 * E0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn longitudinal_shape_load_retains_small_original_storage_under_large_coenergy_bias() {
    let law = almost_flat_source();
    let response = law.evaluate(E0).unwrap();
    assert!(response.coenergy_density_j_m3 > 1e10 * response.energy_density_j_m3);
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        vec![[0, 1, 2, 3]],
        vec![SpatialDielectricMaterial::Reversible(law)],
        limits(),
        &mut || false,
    )
    .unwrap();
    let ground = faces_at(&problem, &positions, 0.0);
    let electrodes = problem
        .prepare_electrodes(&[&ground], &mut || false)
        .unwrap();
    let field = electrodes
        .solve(
            &[4.0 * response.displacement_c_m2],
            &[ElectrodeControl::Voltage(0.0)],
            &[0.0, E0, 0.0, 0.0],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(field.report().newton_iterations, 0);
    let u = [[0.0; 3], [1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    let exact = -response.energy_density_j_m3 / 6.0;
    near(field.virtual_work(&u, &mut || false).unwrap(), exact, exact);
    near(
        field.shape_forces(&mut || false).unwrap().nodal_force_n()[1][0],
        exact,
        exact,
    );
}

#[test]
fn undeclared_tensor_transport_invalid_motion_and_cancellation_preserve_the_field() {
    let (positions, cells) = box_mesh(&[0.0, GAP], WIDTH, DEPTH, 1, 1);
    let names = PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES;
    let material = scalar(&[
        (names[0], EPSILON),
        (names[1], EPSILON),
        (names[2], EPSILON),
    ]);
    let law =
        SpatialDielectricLaw::orthotropic([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .unwrap();
    // Even a numerically isotropic tensor source lacks a declared frame/strain
    // transport law; an explicitly isotropic source remains admitted.
    for material in [
        SpatialDielectricMaterial::Linear { material, law },
        SpatialDielectricMaterial::Reversible(nonlinear()),
    ] {
        let tensor = matches!(&material, SpatialDielectricMaterial::Linear { .. });
        let problem = SpatialNonlinearElectrostatics::new(
            positions.clone(),
            cells.clone(),
            vec![material; cells.len()],
            limits(),
            &mut || false,
        )
        .unwrap();
        let terminals = [
            faces_at(&problem, &positions, 0.0),
            faces_at(&problem, &positions, GAP),
        ];
        let electrodes = problem
            .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
            .unwrap();
        let field = electrodes
            .solve(
                &vec![0.0; cells.len()],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(2.0 * V0),
                ],
                &vec![0.0; positions.len()],
                controls(),
                &mut || false,
            )
            .unwrap();
        let energy = field.energy_j();
        let zero = vec![[0.0; 3]; positions.len()];
        if tensor {
            assert!(matches!(
                field.shape_forces(&mut || false),
                Err(NonlinearElectrostaticError::Electric(
                    ElectrostaticError::Invalid(_)
                ))
            ));
            assert!(matches!(
                field.virtual_work(&zero, &mut || false),
                Err(NonlinearElectrostaticError::Electric(
                    ElectrostaticError::Invalid(_)
                ))
            ));
        } else {
            assert!(
                field
                    .virtual_work(&zero[..zero.len() - 1], &mut || false)
                    .is_err()
            );
            let mut invalid = zero.clone();
            invalid[0][2] = f64::NAN;
            assert!(field.virtual_work(&invalid, &mut || false).is_err());
            assert!(matches!(
                field.shape_forces(&mut || true),
                Err(NonlinearElectrostaticError::Electric(
                    ElectrostaticError::Cancelled
                ))
            ));
            let mut polls = 0;
            assert!(matches!(
                field.virtual_work(&zero, &mut || {
                    polls += 1;
                    polls > 2
                }),
                Err(NonlinearElectrostaticError::Electric(
                    ElectrostaticError::Cancelled
                ))
            ));
            assert_eq!(field.virtual_work(&zero, &mut || false).unwrap(), 0.0);
        }
        assert_eq!(field.energy_j(), energy);
    }
}
