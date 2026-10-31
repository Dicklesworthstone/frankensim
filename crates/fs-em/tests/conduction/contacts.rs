//! Synthetic material and pressure substitutions in the native contact field.
//! These manufactured laws test the integration; they are not empirical claims.
use super::*;
use fs_em::conduction::{ConductionContact, motion::ConductionMotion};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    ResolvedElectricalInterface, resolve_electrical_interface_state,
};
use fs_qty::{
    Pressure, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

pub(super) fn mesh(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "synthetic matching material contact",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}

fn interface_card(
    first: MaterialStateId,
    second: MaterialStateId,
    scale: f64,
) -> InterfaceSystemCard {
    interface_card_with_domain(
        first,
        second,
        scale,
        ValidityDomain::unconstrained().with_quantity(
            "normal_pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
            2e5,
        ),
    )
}

fn interface_card_with_domain(
    first: MaterialStateId,
    second: MaterialStateId,
    scale: f64,
    validity: ValidityDomain,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "normal_pressure".into(),
                abscissa_dims: Pressure::DIMS,
                knots: vec![(1e5, 0.5 * scale), (2e5, 0.25 * scale)],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity,
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic pressure-dependent contact, not measured material data".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: first,
            texture_frame: "first synthetic surface".into(),
        },
        SurfaceSpec {
            material: second,
            texture_frame: "second synthetic surface".into(),
        },
        SystemContext {
            medium: "explicit synthetic interface".into(),
            third_body: None,
            environment: "fixed synthetic environment".into(),
            history: "fixed synthetic preparation".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn pressure_point(pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "normal_pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            pressure,
        )
        .unwrap()
}

pub(super) fn interface(
    left: &ResolvedConductorState,
    right: &ResolvedConductorState,
    pressure: f64,
    scale: f64,
) -> ResolvedElectricalInterface {
    resolve_electrical_interface_state(
        &interface_card(
            left.material().material().clone(),
            right.material().material().clone(),
            scale,
        ),
        &pressure_point(pressure),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

pub(super) fn binding(source: ResolvedElectricalInterface) -> ConductionContact {
    ConductionContact {
        first_vertices: [1, 2],
        second_vertices: [4, 7],
        source,
    }
}

fn specimen(
    left_sigma: f64,
    right_sigma: f64,
    pressure: f64,
    scale: f64,
    depth: f64,
) -> PlanarConduction {
    let left = conductor(left_sigma, ElectricalPropertyBasis::Conductivity);
    let right = conductor(right_sigma, ElectricalPropertyBasis::Resistivity);
    let contact = binding(interface(&left, &right, pressure, scale));
    PlanarConduction::new_with_contacts(
        mesh(depth),
        vec![left.clone(), left, right.clone(), right],
        vec![contact],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn source_pressure_and_bulk_substitutions_reach_field_current_heat_and_circuit() {
    for (left, right, pressure, scale, depth) in [
        (2.0, 4.0, 1e5, 1.0, 0.2),
        (2.0, 4.0, 2e5, 1.0, 0.2),
        (8.0, 3.0, 1e5, 1.0, 0.2),
        (2.0, 4.0, 1.5e5, 3.0, 0.6),
    ] {
        let problem = specimen(left, right, pressure, scale, depth);
        let r_specific = scale * (0.5 - 0.25 * (pressure - 1e5) / 1e5);
        let resistance = (1.0 / left + r_specific + 1.0 / right) / depth;
        let reduced = problem
            .resistance(&[5, 6], &[0, 3], controls(), &mut || false)
            .unwrap();
        near(reduced.resistance_ohm(), resistance);
        let Branch::Resistor { ohms } = reduced.circuit_branch() else {
            panic!("native resistor required")
        };
        near(ohms, resistance);
        assert_eq!(
            problem.contacts()[0].source.resolved().query_point(),
            &[("normal_pressure".into(), pressure)]
        );
        for voltage in [3.0, -2.0, 0.0] {
            let field = reduced
                .field_at_voltage(voltage, 1e-10, 1e-10, &mut || false)
                .unwrap();
            let current = voltage / resistance;
            let density = current / depth;
            let contact = &field.contact_response()[0];
            near(contact.current_a, -current);
            for endpoint in 0..2 {
                near(contact.voltage_jump_v[endpoint], -density * r_specific);
                near(contact.current_density_a_m2[endpoint], -density);
            }
            for (cell, power) in field.cell_power_w().iter().enumerate() {
                let sigma = if cell < 2 { left } else { right };
                near(*power, 0.5 * depth * density * density / sigma);
                near(field.current_density_a_m2()[cell][0], -density);
                near(field.current_density_a_m2()[cell][1], 0.0);
            }
            near(contact.power_w, depth * r_specific * density * density);
            near(field.power_w(), current * voltage);
            near(
                field.cell_power_w().iter().sum::<f64>() + contact.power_w,
                field.power_w(),
            );
            near(field.boundary_current_a(), 0.0);
            near(field.maximum_free_current_residual_a(), 0.0);
            near(field.power_balance_residual_w(), 0.0);
            let direct = problem
                .solve(
                    &[(0, 0.0), (3, 0.0), (5, voltage), (6, voltage)],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            for (a, b) in field.voltage_v().iter().zip(direct.voltage_v()) {
                near(*a, *b);
            }
        }
    }
}

#[test]
fn contact_connectivity_anchors_current_driven_otherwise_floating_body() {
    let problem = specimen(2.0, 4.0, 1e5, 1.0, 0.2);
    let currents = [ConductionCurrentTrace {
        vertices: [5, 6],
        inward_current_density_a_m2: 3.0,
    }];
    let field = problem
        .solve_with_current_traces(&[(0, 7.0), (3, 7.0)], &currents, controls(), &mut || false)
        .unwrap();
    near(field.voltage_v()[5], 7.0 + 3.0 * (0.5 + 0.5 + 0.25));
    near(field.voltage_v()[6], field.voltage_v()[5]);
    near(field.contact_response()[0].current_a, -0.6);
    near(field.power_w(), 0.2 * 9.0 * 1.25);
    let unconnected = PlanarConduction::new(
        problem.mesh().clone(),
        problem.conductors().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        unconnected.solve_with_current_traces(
            &[(0, 7.0), (3, 7.0)],
            &currents,
            controls(),
            &mut || false,
        ),
        Err(ConductionError::Invalid(_))
    ));
    for current in [-0.6, 0.0, 0.6] {
        let field = problem
            .solve_current(&[5, 6], &[0, 3], current, controls(), &mut || false)
            .unwrap();
        near(field.voltage_v()[5], current * 1.25 / 0.2);
        near(field.contact_response()[0].current_a, -current);
        near(field.power_w(), current * current * 1.25 / 0.2);
    }
}

#[test]
fn spatially_varying_contact_heat_matches_independent_gauss_integral() {
    let problem = specimen(2.0, 4.0, 1e5, 1.0, 0.3);
    let field = problem
        .solve(
            &[(0, 0.0), (3, 0.0), (5, 1.0), (6, 2.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let response = &field.contact_response()[0];
    let [d0, d1] = response.voltage_jump_v;
    assert!(
        (d0 - d1).abs() > 1e-6,
        "fixture must exercise nonuniform trace"
    );
    let offset = 1.0 / 3.0_f64.sqrt();
    let mut integral = 0.0;
    for x in [(1.0 - offset) / 2.0, (1.0 + offset) / 2.0] {
        let jump = (1.0 - x) * d0 + x * d1;
        integral += 0.5 * 0.3 * jump * jump / 0.5;
    }
    near(response.power_w, integral);
    near(
        field.power_w(),
        field.cell_power_w().iter().sum::<f64>() + integral,
    );
    near(field.power_balance_residual_w(), 0.0);
}

#[test]
fn contacts_refuse_duplicate_nonmatching_perfect_and_unsupported_sources() {
    let material = conductor(2.0, ElectricalPropertyBasis::Conductivity);
    let contact = binding(interface(&material, &material, 1e5, 1.0));
    let build = |contacts| {
        PlanarConduction::new_with_contacts(
            mesh(0.2),
            vec![material.clone(); 4],
            contacts,
            limits(),
            &mut || false,
        )
    };
    assert!(matches!(
        build(vec![contact.clone(), contact.clone()]),
        Err(ConductionError::Invalid(_))
    ));
    let mut crossed = contact.clone();
    crossed.second_vertices.reverse();
    assert!(matches!(
        build(vec![crossed]),
        Err(ConductionError::Invalid(_))
    ));
    let mut nonmatching = contact.clone();
    nonmatching.second_vertices = [5, 6];
    assert!(matches!(
        build(vec![nonmatching]),
        Err(ConductionError::Invalid(_))
    ));
    assert!(matches!(
        build(vec![binding(interface(&material, &material, 1e5, 0.0))]),
        Err(ConductionError::Invalid(_))
    ));

    let mut wrong_material = material.material().material().clone();
    wrong_material.revision += 1;
    let wrong_card = interface_card(wrong_material, material.material().material().clone(), 1.0);
    let wrong = resolve_electrical_interface_state(
        &wrong_card,
        &pressure_point(1e5),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        build(vec![binding(wrong)]),
        Err(ConductionError::Invalid(_))
    ));
    // The generic resolver can preserve an untyped pressure axis, but a static
    // field consumer cannot infer its physical meaning from its spelling.
    let untyped = resolve_electrical_interface_state(
        &interface_card_with_domain(
            material.material().material().clone(),
            material.material().material().clone(),
            1.0,
            ValidityDomain::unconstrained().with("normal_pressure", 1e5, 2e5),
        ),
        &QueryPoint::new().with("normal_pressure", 1e5).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        build(vec![binding(untyped)]),
        Err(ConductionError::Invalid(_))
    ));

    let mut limited = limits();
    limited.max_triplets = 4 * 9 + 15;
    assert!(matches!(
        PlanarConduction::new_with_contacts(
            mesh(0.2),
            vec![material.clone(); 4],
            vec![contact.clone()],
            limited,
            &mut || false,
        ),
        Err(ConductionError::Invalid(_))
    ));
    assert!(matches!(
        PlanarConduction::new_with_contacts(
            mesh(0.2),
            vec![material; 4],
            vec![contact.clone()],
            limits(),
            &mut || true,
        ),
        Err(ConductionError::Cancelled)
    ));

    let problem = specimen(2.0, 4.0, 1e5, 1.0, 0.2);
    assert!(matches!(
        problem.solve(
            &[(0, 0.0), (3, 0.0), (1, 1.0), (2, 1.0)],
            controls(),
            &mut || false,
        ),
        Err(ConductionError::Invalid(_))
    ));
    assert!(matches!(
        problem.solve_with_current_traces(
            &[(0, 0.0), (3, 0.0)],
            &[ConductionCurrentTrace {
                vertices: [4, 7],
                inward_current_density_a_m2: 1.0
            }],
            controls(),
            &mut || false,
        ),
        Err(ConductionError::Invalid(_))
    ));
    assert!(matches!(
        problem.solve_with_motion(
            &[(0, 0.0), (3, 0.0), (5, 1.0), (6, 1.0)],
            &[],
            ConductionMotion {
                velocity_m_s: &[[0.0; 2]; 4],
                normal_flux_density_t: &[0.0; 4]
            },
            controls(),
            &mut || false,
        ),
        Err(ConductionError::Invalid(_))
    ));
}
