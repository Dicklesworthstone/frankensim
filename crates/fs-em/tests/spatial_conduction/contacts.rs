//! Original bulk and interface data change the same native tetrahedral field.
//! These authored laws exercise integration; they are not empirical claims.
use super::*;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    ResolvedElectricalInterface, resolve_electrical_interface_state,
};
use fs_qty::{
    Pressure, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn conductor(resistivity: f64) -> ResolvedConductorState {
    ResolvedConductorState::try_from_material_with_basis(
        material(ElectricalPropertyBasis::Resistivity, &[resistivity]),
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap()
}

fn interface(
    first: &ResolvedConductorState,
    second: &ResolvedConductorState,
    pressure: f64,
    scale: f64,
) -> ResolvedElectricalInterface {
    let quantity =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
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
                knots: vec![(1e5, 3.0 * scale), (2e5, 1.5 * scale)],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "normal_pressure",
                quantity,
                1e5,
                2e5,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured pressure-dependent 3D contact; not measured data".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: first.material().material().clone(),
            texture_frame: "first surface".into(),
        },
        SurfaceSpec {
            material: second.material().material().clone(),
            texture_frame: "second surface".into(),
        },
        SystemContext {
            medium: "authored fixed interface".into(),
            third_body: None,
            environment: "authored environment".into(),
            history: "fixed preparation".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_electrical_interface_state(
        &card,
        &QueryPoint::new()
            .with_quantity("normal_pressure", quantity, pressure)
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn mesh(area: f64) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let (mut vertices, mut tets) = slab(&[0.0, 1.0]);
    let (other, other_tets) = slab(&[1.0, 2.0]);
    vertices.extend(other);
    tets.extend(other_tets.into_iter().map(|tet| tet.map(|v| v + 8)));
    for vertex in &mut vertices {
        vertex[1] *= area;
    }
    (vertices, tets)
}

fn contacts(source: ResolvedElectricalInterface) -> Vec<SpatialConductionContact> {
    vec![
        SpatialConductionContact {
            first_vertices: [4, 5, 7],
            second_vertices: [8, 9, 11],
            source: source.clone(),
        },
        SpatialConductionContact {
            first_vertices: [4, 6, 7],
            second_vertices: [8, 10, 11],
            source,
        },
    ]
}

fn specimen(rho: [f64; 2], pressure: f64, scale: f64, area: f64) -> SpatialConduction {
    let materials = rho.map(conductor);
    let (vertices, tets) = mesh(area);
    let source = interface(&materials[0], &materials[1], pressure, scale);
    SpatialConduction::new_with_contacts(
        vertices,
        tets,
        (0..12).map(|i| materials[i / 6].clone()).collect(),
        contacts(source),
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn pressure_material_area_and_polarity_change_bulk_and_surface_work_together() {
    for (rho, pressure, scale, area) in [
        ([2.0, 5.0], 1e5, 1.0, 1.0),
        ([2.0, 5.0], 2e5, 1.0, 1.0),
        ([7.0, 1.0], 1.5e5, 2.0, 0.4),
        ([2.0, 5.0], 1e5, 1.0, 1.6),
    ] {
        let problem = specimen(rho, pressure, scale, area);
        let resistance = scale * (3.0 - 1.5 * (pressure - 1e5) / 1e5);
        assert_eq!(problem.contacts().len(), 2);
        assert_eq!(problem.contact_cells().len(), 2);
        for voltage in [10.0, -4.0, 0.0] {
            let density = voltage / (rho[0] + resistance + rho[1]);
            let field = problem
                .solve(
                    &ends(&problem, 2.0, Some(voltage), 0.0),
                    controls(),
                    &mut || false,
                )
                .unwrap();
            for (vertex, p) in problem.positions_m().iter().enumerate() {
                let expected = if vertex < 8 {
                    density * rho[0] * p[0]
                } else {
                    density * (rho[0] + resistance + rho[1] * (p[0] - 1.0))
                };
                near(field.voltage_v[vertex], expected);
            }
            for cell in 0..12 {
                near(field.current_density_a_m2[cell][0], -density);
                near(field.current_density_a_m2[cell][1], 0.0);
                near(field.current_density_a_m2[cell][2], 0.0);
                near(
                    field.cell_power_w[cell],
                    area * density * density * rho[cell / 6] / 6.0,
                );
            }
            for contact in &field.contacts {
                for node in 0..3 {
                    near(contact.voltage_jump_v[node], -density * resistance);
                    near(contact.current_density_a_m2[node], -density);
                }
                near(contact.current_a, -0.5 * area * density);
                near(contact.power_w, 0.5 * area * resistance * density * density);
            }
            near(field.power_w, area * density * voltage);
            near(field.supplied_power_w, field.power_w);
            near(
                field.cell_power_w.iter().sum::<f64>()
                    + field.contacts.iter().map(|c| c.power_w).sum::<f64>(),
                field.power_w,
            );
            near(field.current_residual_a, 0.0);
        }
    }
}

#[test]
fn contact_connectivity_owns_one_current_gauge_and_refresh_retains_original_geometry() {
    let problem = specimen([2.0, 5.0], 1e5, 1.0, 1.0);
    for mut boundary in [
        ends(&problem, 2.0, None, 0.75),
        ends(&problem, 2.0, Some(0.0), 0.75),
    ] {
        for trace in &mut boundary.voltage {
            trace.voltage_v = [17.0; 3];
        }
        let field = problem.solve(&boundary, controls(), &mut || false).unwrap();
        let offset = if boundary.voltage.is_empty() {
            0.0
        } else {
            17.0
        };
        near(field.voltage_v[0], offset);
        near(field.voltage_v[15], offset + 7.5);
        near(field.power_w, 5.625);
        assert_eq!(field.gauges.len(), usize::from(boundary.voltage.is_empty()));
    }
    let (vertices, tets) = mesh(1.0);
    let disconnected = SpatialConduction::new(
        vertices,
        tets,
        problem.conductors().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        disconnected
            .solve(
                &ends(&disconnected, 2.0, None, 0.75),
                controls(),
                &mut || false
            )
            .is_err()
    );
    let changed_contacts = contacts(interface(&conductor(2.0), &conductor(5.0), 2e5, 1.0));
    let changed = problem
        .with_conductors_and_contacts(
            problem.conductors().to_vec(),
            changed_contacts,
            limits(),
            &mut || false,
        )
        .unwrap();
    near(
        changed
            .solve(
                &ends(&changed, 2.0, Some(10.0), 0.0),
                controls(),
                &mut || false,
            )
            .unwrap()
            .power_w,
        100.0 / 8.5,
    );
    let new_bulk = vec![conductor(1.0); 12];
    let changed = problem
        .with_conductors(new_bulk, limits(), &mut || false)
        .unwrap();
    assert_eq!(changed.contacts(), problem.contacts());
    near(
        changed
            .solve(
                &ends(&changed, 2.0, Some(10.0), 0.0),
                controls(),
                &mut || false,
            )
            .unwrap()
            .power_w,
        20.0,
    );
    near(
        problem
            .solve(
                &ends(&problem, 2.0, Some(10.0), 0.0),
                controls(),
                &mut || false,
            )
            .unwrap()
            .power_w,
        10.0,
    );
}

#[test]
fn nonuniform_triangle_jump_integrates_the_full_surface_mass_and_terminal_work() {
    let problem = specimen([2.0, 5.0], 1e5, 1.0, 1.0);
    let boundary = SpatialConductionBoundary {
        voltage: problem
            .exterior_faces()
            .map(|face| ConductionVoltageFace {
                face,
                voltage_v: problem.mesh().faces[face].map(|v| {
                    let p = problem.positions_m()[v as usize];
                    if v < 8 { 2.0 + p[1] + 2.0 * p[2] } else { 0.0 }
                }),
            })
            .collect(),
        current: vec![],
    };
    let field = problem.solve(&boundary, controls(), &mut || false).unwrap();
    // Integrals over the unit square: jump=2+y+2z, integral(jump)=7/2,
    // integral(jump^2)=38/3. The first bulk has |grad V|^2=5 and rho=2.
    near(field.contacts.iter().map(|c| c.current_a).sum(), 3.5 / 3.0);
    near(field.contacts.iter().map(|c| c.power_w).sum(), 38.0 / 9.0);
    near(field.cell_power_w.iter().sum(), 2.5);
    near(field.power_w, 2.5 + 38.0 / 9.0);
    near(field.supplied_power_w, field.power_w);
}

#[test]
fn invalid_trace_ownership_geometry_budget_and_cancel_refuse_without_replacing_field() {
    let problem = specimen([2.0, 5.0], 1e5, 1.0, 1.0);
    let (positions, tets) = mesh(1.0);
    let source = interface(&conductor(2.0), &conductor(5.0), 1e5, 1.0);
    for case in 0..3 {
        let mut positions = positions.clone();
        let mut contacts = contacts(source.clone());
        match case {
            0 => contacts[0].second_vertices.swap(0, 1),
            1 => contacts.push(contacts[0].clone()),
            _ => {
                for p in &mut positions[8..] {
                    p[0] = 2.0 - p[0];
                }
            }
        }
        assert!(
            SpatialConduction::new_with_contacts(
                positions,
                tets.clone(),
                problem.conductors().to_vec(),
                contacts,
                limits(),
                &mut || false
            )
            .is_err()
        );
    }
    let face = problem
        .mesh()
        .faces
        .iter()
        .position(|face| *face == [4, 5, 7])
        .unwrap();
    assert!(problem.exterior_normal(face).is_none());
    assert!(
        problem
            .solve(
                &SpatialConductionBoundary {
                    voltage: vec![ConductionVoltageFace {
                        face,
                        voltage_v: [0.0; 3]
                    }],
                    current: vec![]
                },
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        problem
            .with_conductors(
                problem.conductors().to_vec(),
                WeightedAssemblyLimits {
                    max_triplets: 12 * 36 + 2 * 36 - 1,
                    ..limits()
                },
                &mut || false
            )
            .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        problem.with_conductors(problem.conductors().to_vec(), limits(), &mut || {
            calls += 1;
            calls > 15
        }),
        Err(SpatialConductionError::Cancelled)
    ));
    near(
        problem
            .solve(
                &ends(&problem, 2.0, Some(10.0), 0.0),
                controls(),
                &mut || false,
            )
            .unwrap()
            .power_w,
        10.0,
    );
}

#[test]
fn total_current_uses_the_same_contact_connection_and_local_surface_heat() {
    for (pressure, area) in [(1e5, 0.4), (2e5, 1.6)] {
        let problem = specimen([2.0, 5.0], pressure, 1.0, area);
        let contact_resistivity = 3.0 - 1.5 * (pressure - 1e5) / 1e5;
        let expected_resistance = (7.0 + contact_resistivity) / area;
        let faces = super::terminals::terminal_faces(&problem, 2.0);
        let response = problem
            .resistance(&faces[0], &faces[1], controls(), &mut || false)
            .unwrap();
        near(response.resistance_ohm(), expected_resistance);
        for current in [-0.3, 0.0, 0.75] {
            let field = response
                .field_at_current(current, 1e-10, 1e-9, &mut || false)
                .unwrap();
            near(field.power_w, current * current * expected_resistance);
            near(
                field.cell_power_w.iter().sum(),
                current * current * 7.0 / area,
            );
            near(
                field.contacts.iter().map(|c| c.power_w).sum(),
                current * current * contact_resistivity / area,
            );
            for contact in &field.contacts {
                near(contact.current_a, -0.5 * current);
            }
            near(field.supplied_power_w, field.power_w);
        }
        let contact_face = problem
            .mesh()
            .faces
            .iter()
            .position(|f| *f == [4, 5, 7])
            .unwrap();
        assert!(
            problem
                .resistance(&[contact_face], &faces[1], controls(), &mut || false)
                .is_err()
        );
    }
}
