//! G1/G3/G4: source updates preserve real memory and continue native mechanics.
use super::*;
use fs_couple::tensor_material::{MaterialChangeBudget, solid::source::SolidMaterialSource};
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch};

fn changing_card(density_changes: bool) -> (MaterialCard, TensorPronyProperties) {
    let (mut claims, properties) = fixture(1000.0, 10.0, 0.4);
    for claim in &mut claims {
        let factor = match claim.key.name() {
            "young_modulus" => Some(0.4),
            "branch_tau" => Some(0.15),
            "density" if density_changes => Some(1.2),
            name if name.starts_with("branch_c") => Some(2.5),
            _ => None,
        };
        if let Some(factor) = factor {
            let (value, dims) = match &claim.value {
                PropertyValue::Scalar { value, dims } => (*value, *dims),
                _ => unreachable!(),
            };
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(300.0, value), (320.0, value * factor)],
                dims,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
    }
    (card(claims), properties)
}
fn bind(card: &MaterialCard, properties: TensorPronyProperties) -> ResolvedTensorMaterial {
    ResolvedTensorMaterial::resolve(
        card,
        &query(),
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn solid_projection() -> SolidProjectionBudget {
    SolidProjectionBudget {
        max_coordinates: 1,
        max_basis_scalars: 1,
        max_material_points: 1,
        max_scalar_terms: 100_000,
    }
}
fn build_elements(material: &ResolvedTensorMaterial, cx: &Cx<'_>) -> ViscoElements {
    let field = TetViscoMaterial::from_source(material, IDENTITY).unwrap();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &[[0, 1, 2, 3]],
        },
        TetViscoField::Uniform(&field),
        &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 100,
        },
        cx,
    )
    .unwrap()
}
fn body<'s>(
    elements: &'s ViscoElements,
    source: ResolvedTensorMaterial,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> TensorSolidBody<'s> {
    TensorSolidBody::from_material_sources(
        elements
            .project_basis(&[vec![1.0]], solid_projection(), cx)
            .unwrap(),
        vec![SolidMaterialSource {
            material: source,
            material_to_mesh: IDENTITY,
        }],
        dynamics(),
        gate,
    )
    .unwrap()
}
fn allowance(work: f64) -> MaterialChangeBudget {
    MaterialChangeBudget {
        max_absolute_parameter_work_j: work,
        relative_memory_tolerance: 1e-12,
    }
}

// Independently integrate the physical strain equations, without factor
// coordinates, FE assembly, a pH operator, or the production exponential.
// y=(displacement,velocity,Maxwell spring strain,physical heat,force work).
fn reference(mut y: [f64; 5], e: f64, branch: f64, tau: f64, force: f64, dt: f64) -> [f64; 5] {
    let volume = 1.0 / 6.0;
    let mass = 10.0 / 60.0;
    let rhs = |y: [f64; 5]| {
        [
            y[1],
            (force - volume * (e * y[0] + branch * y[2])) / mass,
            y[1] - y[2] / tau,
            volume * branch * y[2] * y[2] / tau,
            force * y[1],
        ]
    };
    let h = dt / 1024.0;
    for _ in 0..1024 {
        let a = rhs(y);
        let b = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| y[i] + h * c[i]));
        y = core::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn source_temperature_updates_preserve_spring_strain_and_continue_actual_tet_motion() {
    with_cx(|cx, gate| {
        let (card, properties) = changing_card(false);
        let material = bind(&card, properties);
        let elements = build_elements(&material, cx);
        let cold = body(&elements, material.clone(), cx, gate);
        let hot_query = query().with("T", 320.0).unwrap();
        let hot = cold
            .at_query_points(&[&card], &[hot_query], solid_projection(), dynamics(), gate)
            .unwrap();
        assert_eq!(cold.body().mass_matrix_kg(), hot.body().mass_matrix_kg());
        let old_source = &cold.material_sources().unwrap()[0].material;
        let new_source = &hot.material_sources().unwrap()[0].material;
        assert_eq!(
            old_source.state().card_identity(),
            new_source.state().card_identity()
        );
        assert_ne!(old_source.state().identity(), new_source.state().identity());
        for (a, b) in old_source
            .state()
            .properties()
            .iter()
            .zip(new_source.state().properties())
        {
            assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        }
        let mut history = cold
            .history(&[0.002], &[0.01], InitialMemory::Instantaneous, gate)
            .unwrap();
        let force = [[0.0; 3], [0.0; 3], [0.0; 3], [0.0, 0.0, 0.02]];
        history.advance(&force, 0.004, 1e-10, gate).unwrap();
        let expected = reference(
            [0.002, 0.01, 0.002, 0.0, 0.0],
            1000.0,
            280.0,
            0.02,
            0.02,
            0.004,
        );
        near(
            history.material_history().displacement_m()[0],
            expected[0],
            2e-12,
        );
        let before = history.material_history().state().to_vec();
        let old_heat = history.material_history().dissipated_heat_j();
        let old_work = history.material_history().supplied_work_j();
        let change = history
            .rebind_materials(&hot, allowance(1.0), gate)
            .unwrap();
        let expected_work = (1.0 / 12.0)
            * ((400.0 - 1000.0) * expected[0].powi(2) + (700.0 - 280.0) * expected[2].powi(2));
        near(change.parameter_work_j, expected_work, 2e-12);
        let mut continued = change.history;
        let actual = continued.material_history();
        assert_eq!(
            &actual.state()[..2],
            &before[..2],
            "configuration/momentum must not jump"
        );
        assert_eq!(actual.time_s(), 0.004);
        assert_eq!(actual.steps(), 1);
        assert_eq!(actual.dissipated_heat_j(), old_heat);
        assert_eq!(actual.supplied_work_j(), old_work);
        // Only the z normal branch factor is excited by this nodal strain.
        near(
            actual.state()[4] / (700.0_f64 / 6.0).sqrt(),
            expected[2],
            2e-12,
        );
        assert!(
            (actual.state()[4] - before[4]).abs() > 1e-4,
            "copying factor z would spuriously freeze energy and alter physical strain"
        );
        near(actual.balance_residual_j(), 0.0, 2e-11);
        continued.advance(&force, 0.005, 1e-10, gate).unwrap();
        let expected = reference(expected, 400.0, 700.0, 0.003, 0.02, 0.005);
        let actual = continued.material_history();
        near(actual.displacement_m()[0], expected[0], 2e-12);
        near(actual.velocity_m_s()[0], expected[1], 2e-11);
        near(
            actual.state()[4] / (700.0_f64 / 6.0).sqrt(),
            expected[2],
            2e-12,
        );
        near(actual.dissipated_heat_j(), expected[3], 2e-12);
        near(actual.supplied_work_j(), expected[4], 2e-12);
        near(actual.balance_residual_j(), 0.0, 2e-11);
        let back = continued
            .rebind_materials(&cold, allowance(1.0), gate)
            .unwrap();
        let expected_back = (1.0 / 12.0)
            * ((1000.0 - 400.0) * expected[0].powi(2) + (280.0 - 700.0) * expected[2].powi(2));
        near(back.parameter_work_j, expected_back, 2e-12);
        near(
            back.history.material_history().state()[4] / (280.0_f64 / 6.0).sqrt(),
            expected[2],
            2e-12,
        );
        near(
            back.history.material_history().parameter_work_j(),
            expected_work + expected_back,
            2e-12,
        );
        near(
            back.history.material_history().balance_residual_j(),
            0.0,
            2e-11,
        );
        assert_eq!(history.material_history().state(), before);
    });
}

#[test]
fn full_tensor_factor_change_preserves_hidden_shear_strains_and_refuses_singular_recovery() {
    let old_rows = [
        [2.0, 0.3, -0.2, 0.1, 0.0, 0.0],
        [0.0, 3.0, 0.2, 0.0, 0.4, 0.0],
        [0.0, 0.0, 4.0, 0.1, 0.0, 0.5],
        [0.0, 0.0, 0.0, 5.0, 0.2, 0.0],
        [0.0, 0.0, 0.0, 0.0, 6.0, -0.4],
        [0.0, 0.0, 0.0, 0.0, 0.0, 7.0],
    ];
    let new_rows = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 3.0],
        [0.0, 0.0, 0.0, 0.0, 2.0, 0.0],
        [0.0, 0.0, 0.0, 4.0, 0.0, 0.0],
        [0.0, 0.0, 5.0, 0.0, 0.0, 0.0],
        [0.0, 6.0, 0.0, 0.0, 0.0, 0.0],
        [7.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let project = |rows, tau| {
        let eq =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 10.0 } else { 0.0 }));
        TensorProny::new_mandel(
            eq,
            vec![TensorPronyBranch::from_gram_rows(rows, tau).unwrap()],
            0.1,
        )
        .unwrap()
        .project_strain(
            &[[0.2, 0.0, -0.1, 0.3, 0.1, 0.0]],
            0.7,
            projection(),
            || false,
        )
        .unwrap()
    };
    let old = project(old_rows, 0.2);
    let new = project(new_rows, 0.03);
    let strain = [0.003, -0.002, 0.004, 0.001, -0.006, 0.005];
    let memory: Vec<f64> = old_rows
        .iter()
        .map(|r| 0.7_f64.sqrt() * r.iter().zip(strain).map(|(a, b)| a * b).sum::<f64>())
        .collect();
    let expected: Vec<f64> = new_rows
        .iter()
        .map(|r| 0.7_f64.sqrt() * r.iter().zip(strain).map(|(a, b)| a * b).sum::<f64>())
        .collect();
    let result = old
        .transfer_factor_memory(&new, &memory, 1e-12, || false)
        .unwrap();
    for (a, b) in result.iter().zip(&expected) {
        near(*a, *b, 1e-13);
    }
    let back = new
        .transfer_factor_memory(&old, &result, 1e-12, || false)
        .unwrap();
    for (a, b) in back.iter().zip(&memory) {
        near(*a, *b, 1e-13);
    }
    let singular_rows = core::array::from_fn(|i| if i == 5 { [0.0; 6] } else { old_rows[i] });
    let singular = project(singular_rows, 0.2);
    let same_factors = project(singular_rows, 0.7);
    let mut singular_memory = memory.clone();
    singular_memory[5] = 0.0;
    assert!(
        singular
            .transfer_factor_memory(&new, &singular_memory, 1e-12, || false)
            .is_err()
    );
    assert!(
        singular
            .transfer_factor_memory(&new, &[0.0; 6], 1e-12, || false)
            .is_err()
    );
    assert!(
        old.transfer_factor_memory(&singular, &memory, 1e-12, || false)
            .is_err(),
        "a new nullspace would erase strains needed by a later source return"
    );
    assert!(
        old.transfer_factor_memory(&singular, &[0.0; 6], 1e-12, || false)
            .is_err()
    );
    assert_eq!(
        singular
            .transfer_factor_memory(&same_factors, &singular_memory, 1e-12, || false)
            .unwrap(),
        singular_memory
    );
    assert!(
        old.transfer_factor_memory(&new, &memory, 1e-12, || true)
            .is_err()
    );
    assert!(
        old.transfer_factor_memory(&new, &memory, 1.0, || false)
            .is_err()
    );
}

#[test]
fn source_domain_mass_work_and_cancellation_refusals_preserve_the_accepted_history() {
    with_cx(|cx, gate| {
        let (card, properties) = changing_card(false);
        let material = bind(&card, properties.clone());
        let elements = build_elements(&material, cx);
        let cold = body(&elements, material.clone(), cx, gate);
        let history = cold
            .history(&[0.002], &[0.01], InitialMemory::Instantaneous, gate)
            .unwrap();
        let before = history.material_history().state().to_vec();
        let hot_query = query().with("T", 320.0).unwrap();
        let hot = cold
            .at_query_points(
                &[&card],
                &[hot_query.clone()],
                solid_projection(),
                dynamics(),
                gate,
            )
            .unwrap();
        assert!(
            history
                .rebind_materials(&hot, allowance(0.0), gate)
                .is_err()
        );
        assert!(
            cold.at_query_points(
                &[&card],
                &[query().with("T", 340.0).unwrap()],
                solid_projection(),
                dynamics(),
                gate
            )
            .is_err()
        );
        assert!(
            cold.at_query_points(&[], &[], solid_projection(), dynamics(), gate)
                .is_err()
        );
        let (foreign, _) = changing_card(true);
        assert!(
            cold.at_query_points(
                &[&foreign],
                &[hot_query.clone()],
                solid_projection(),
                dynamics(),
                gate
            )
            .is_err()
        );
        let foreign_material = bind(&foreign, properties);
        let foreign_elements = build_elements(&foreign_material, cx);
        let foreign_body = body(&foreign_elements, foreign_material, cx, gate);
        assert!(
            foreign_body
                .at_query_points(
                    &[&foreign],
                    &[hot_query],
                    solid_projection(),
                    dynamics(),
                    gate
                )
                .is_err(),
            "changing density on unchanged geometry must not change closed specimen mass"
        );
        assert!(
            history
                .rebind_materials(&foreign_body, allowance(1.0), gate)
                .is_err()
        );
        let wrong = elements
            .project_basis(&[vec![1.0]], solid_projection(), cx)
            .unwrap();
        assert!(TensorSolidBody::from_material_sources(wrong, vec![], dynamics(), gate).is_err());
        gate.request();
        assert!(
            history
                .rebind_materials(&hot, allowance(1.0), gate)
                .is_err()
        );
        assert_eq!(history.material_history().state(), before);
        assert_eq!(history.material_history().time_s(), 0.0);
        assert_eq!(history.material_history().steps(), 0);
        assert_eq!(history.material_history().parameter_work_j(), 0.0);
        assert_eq!(history.material_history().dissipated_heat_j(), 0.0);
    });
}
