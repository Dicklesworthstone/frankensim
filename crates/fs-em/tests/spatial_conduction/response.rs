//! Actual 3D source response against independent one-dimensional field algebra.
use super::contacts::{contacts as contact_faces, interface};
use super::*;
use fs_em::spatial::conduction::response::*;
use fs_material::conductor::response::{ConductorMagnitude, ScalarConductorResponse};
use std::collections::BTreeSet;

fn law(coordinate: Option<ConductorMagnitude>, beta: f64, upper: f64) -> ScalarConductorResponse {
    let basis = ElectricalPropertyBasis::Conductivity;
    let key = PropertyKey::new(basis.property_name(), basis.dimensions());
    let mut claims = ClaimSet::new();
    let mut query = QueryPoint::new();
    let (value, validity, interpolation, axes) = if let Some(magnitude) = coordinate {
        query = query
            .with_quantity("operating", magnitude.quantity(), 0.0)
            .unwrap();
        (
            PropertyValue::Curve {
                abscissa: "operating".into(),
                abscissa_dims: magnitude.quantity().dims(),
                knots: vec![(0.0, 1.0), (upper, 1.0 + beta * upper)],
                dims: basis.dimensions(),
            },
            ValidityDomain::unconstrained().with_quantity(
                "operating",
                magnitude.quantity(),
                0.0,
                upper,
            ),
            InterpolationPolicy::LinearInside,
            vec![("operating", magnitude)],
        )
    } else {
        (
            PropertyValue::Scalar {
                value: 1.0,
                dims: basis.dimensions(),
            },
            ValidityDomain::unconstrained(),
            InterpolationPolicy::ConstantWithinValidity,
            vec![],
        )
    };
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured nonlinear 3D transport; not measured material".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid".into(),
            process: "bounded scalar field response".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &query,
        &[basis.requirement_for(&key).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = ResolvedConductorState::try_from_material_with_basis(material, basis).unwrap();
    ScalarConductorResponse::new(&card, original, &axes).unwrap()
}

fn response_controls() -> SpatialConductionResponseControls {
    SpatialConductionResponseControls {
        maximum_iterations: 64,
        maximum_field_solves: 128,
        absolute_conductivity_tolerance_s_m: 1e-12,
        relative_conductivity_tolerance: 1e-11,
        current_density_tolerance_a_m2: 1e-10,
    }
}

fn face_sets(problem: &SpatialConductionResponse, coordinates: &[f64]) -> Vec<Vec<usize>> {
    coordinates
        .iter()
        .map(|&z| {
            problem
                .exterior_faces()
                .filter(|&face| {
                    problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][2])
                        == [z; 3]
                })
                .collect()
        })
        .collect()
}

fn rotate(positions: &mut [[f64; 3]]) {
    // A proper cyclic rotation puts the entire axial field in its THIRD component.
    for p in positions {
        *p = [p[1], p[2], p[0]];
    }
}

#[test]
fn mixed_nonlinear_electrodes_and_contacts_close_actual_fields_and_local_power() {
    let nonlinear = law(Some(ConductorMagnitude::ElectricField), 0.2, 20.0);
    let linear = law(None, 0.0, 20.0);
    let source = interface(nonlinear.original(), nonlinear.original(), 1e5, 1.0 / 3.0);
    let (mut positions, mut tets) = slab(&[0.0, 0.5]);
    for xs in [&[0.5, 1.0][..], &[1.0, 1.5, 2.0][..]] {
        let (part, elements) = slab(xs);
        let offset = positions.len() as u32;
        positions.extend(part);
        tets.extend(elements.into_iter().map(|t| t.map(|v| v + offset)));
    }
    rotate(&mut positions);
    let problem = SpatialConductionResponse::new_with_contacts(
        positions,
        tets,
        (0..24)
            .map(|i| {
                if i < 12 {
                    nonlinear.clone()
                } else {
                    linear.clone()
                }
            })
            .collect(),
        contact_faces(source),
        limits(),
        &mut || false,
    )
    .unwrap();
    let faces = face_sets(&problem, &[0.0, 1.0, 2.0]);
    let terminals = faces.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let matrix = SpatialConductanceControls {
        maximum_terminals: 3,
        maximum_field_solves: 2,
        maximum_response_values: 1000,
        maximum_dense_work: 1000,
        reference_voltage_v: 0.5,
        coefficient_tolerance_s: 1e-10,
    };
    let mixed = SpatialMixedConductanceControls {
        maximum_current_terminals: 1,
        maximum_response_values: 100,
        maximum_dense_work: 100,
    };
    for voltage in [2.0_f64, 5.0] {
        for injection in [0.0, 0.5] {
            // A: J_A=E_A(1+.2 E_A), unit length/area plus R_contact=1.
            // B: J_B=E_B, unit length/area. Middle injection J_A-J_B=s.
            // V=E_A+J_A+J_B gives .4 E_A²+3 E_A=V+s.
            let e =
                2.0 * (voltage + injection) / (3.0 + (9.0 + 1.6 * (voltage + injection)).sqrt());
            let ja = e * (1.0 + 0.2 * e);
            let jb = ja - injection;
            let middle = e + ja;
            let power = ja * e + ja * ja + jb * jb;
            for sign in [-1.0, 1.0] {
                let voltages = [1e10, 1e10 + sign * voltage];
                let currents = [sign * injection];
                let drive = SpatialConductionResponseDrive::MixedTerminals {
                    terminals: &terminals,
                    current_terminals: &[1],
                    voltages_v: &voltages,
                    currents_a: &currents,
                    response: matrix,
                    mixed,
                };
                let mut polls = 0;
                let candidate = problem
                    .solve_candidate(
                        drive,
                        limits(),
                        controls(),
                        response_controls(),
                        &mut || {
                            polls += 1;
                            false
                        },
                    )
                    .unwrap();
                let report = candidate.report();
                let field = &report.field;
                for ((terminal, expected_v), expected_i) in terminals
                    .iter()
                    .zip([0.0, sign * middle, sign * voltage])
                    .zip([-sign * ja, sign * injection, sign * jb])
                {
                    let vertices = terminal
                        .iter()
                        .flat_map(|&face| problem.mesh().faces[face])
                        .collect::<BTreeSet<_>>();
                    near(
                        vertices
                            .iter()
                            .map(|&v| field.nodal_boundary_current_a[v as usize].unwrap())
                            .sum(),
                        expected_i,
                    );
                    for v in vertices {
                        near(field.voltage_v[v as usize], expected_v);
                    }
                }
                for (i, (field_e, field_j)) in field
                    .electric_field_v_m
                    .iter()
                    .zip(&field.current_density_a_m2)
                    .enumerate()
                {
                    near(field_e[0], 0.0);
                    near(field_e[1], 0.0);
                    near(field_e[2], -sign * if i < 12 { e } else { jb });
                    near(field_j[2], -sign * if i < 12 { ja } else { jb });
                }
                near(field.cell_power_w[..12].iter().sum(), ja * e);
                near(field.cell_power_w[12..].iter().sum(), jb * jb);
                near(field.contacts.iter().map(|c| c.power_w).sum(), ja * ja);
                near(field.power_w, power);
                near(field.supplied_power_w, voltage * jb + middle * injection);
                near(field.power_w, field.supplied_power_w);
                near(
                    report.evaluated[0].conductor.material().query().axes()["operating"],
                    e,
                );
                near(
                    candidate.problem().conductors()[0].resistivity_ohm_m(),
                    1.0 / (1.0 + 0.2 * e),
                );
                assert!(report.iterations > 1);
                assert_eq!(report.field_solves, 2 * report.iterations);
                assert_eq!(report.source_evaluations, 24 * report.iterations);
                assert!(
                    report
                        .evaluated
                        .iter()
                        .map(|s| s.power_defect_w.abs())
                        .sum::<f64>()
                        <= controls().power_tolerance_w
                );

                if voltage == 2.0 && injection == 0.5 && sign == 1.0 {
                    let bounded = SpatialConductionResponseControls {
                        maximum_field_solves: 3,
                        ..response_controls()
                    };
                    assert!(matches!(
                        problem
                            .solve_candidate(drive, limits(), controls(), bounded, &mut || false),
                        Err(SpatialConductionResponseError::FieldSolveBudget {
                            completed: 2,
                            required: 2,
                            limit: 3
                        })
                    ));
                    let mut count = 0;
                    assert!(matches!(
                        problem.solve_candidate(
                            drive,
                            limits(),
                            controls(),
                            response_controls(),
                            &mut || {
                                count += 1;
                                count == polls
                            }
                        ),
                        Err(SpatialConductionResponseError::Cancelled)
                    ));
                    let retry = problem
                        .solve_candidate(
                            drive,
                            limits(),
                            controls(),
                            response_controls(),
                            &mut || false,
                        )
                        .unwrap();
                    assert_eq!(retry.report().field.voltage_v, field.voltage_v);
                }
            }
        }
    }
    assert_eq!(
        problem.conductors()[0].material().query().axes()["operating"],
        0.0
    );
}

#[test]
fn current_density_dependent_source_uses_the_third_component_and_refuses_unclosed_states() {
    let source = law(Some(ConductorMagnitude::CurrentDensity), 0.1, 10.0);
    let (mut positions, tets) = slab(&[0.0, 1.0]);
    rotate(&mut positions);
    // Ordinary frozen-source admission still refuses an unclosed field law.
    assert!(
        SpatialConduction::new(
            positions.clone(),
            tets.clone(),
            vec![source.original().clone(); 6],
            limits(),
            &mut || false
        )
        .is_err()
    );
    let problem = SpatialConductionResponse::new(
        positions.clone(),
        tets.clone(),
        vec![source; 6],
        limits(),
        &mut || false,
    )
    .unwrap();
    let faces = face_sets(&problem, &[0.0, 1.0]);
    for current in [-2.0_f64, 0.0, 2.0] {
        let boundary = SpatialConductionBoundary {
            voltage: faces[0]
                .iter()
                .map(|&face| ConductionVoltageFace {
                    face,
                    voltage_v: [0.0; 3],
                })
                .collect(),
            current: faces[1]
                .iter()
                .map(|&face| ConductionCurrentFace {
                    face,
                    inward_current_density_a_m2: current,
                })
                .collect(),
        };
        let candidate = problem
            .solve_candidate(
                &boundary,
                limits(),
                controls(),
                response_controls(),
                &mut || false,
            )
            .unwrap();
        let field = &candidate.report().field;
        let sigma = 1.0 + 0.1 * current.abs();
        for (e, j) in field
            .electric_field_v_m
            .iter()
            .zip(&field.current_density_a_m2)
        {
            near(e[2], -current / sigma);
            near(j[2], -current);
        }
        near(field.power_w, current * current / sigma);
        near(field.supplied_power_w, current * current / sigma);
        near(
            candidate.report().evaluated[0]
                .conductor
                .material()
                .query()
                .axes()["operating"],
            current.abs(),
        );
        assert_eq!(
            candidate.report().field_solves,
            candidate.report().iterations
        );
        if current == 2.0 {
            // A loose sigma tolerance alone must not admit wrong J; even with
            // both loose, the independent local E.(J-J_source) gate must refuse.
            for current_density_tolerance_a_m2 in [1e-10, 1e3] {
                let once = SpatialConductionResponseControls {
                    maximum_iterations: 1,
                    absolute_conductivity_tolerance_s_m: 1e3,
                    current_density_tolerance_a_m2,
                    ..response_controls()
                };
                assert!(
                    matches!(problem.solve_candidate(&boundary, limits(), controls(), once, &mut || false),
                    Err(SpatialConductionResponseError::NotConverged { absolute_power_defect_w, .. })
                        if absolute_power_defect_w > 0.5)
                );
            }
            let limited = SpatialConductionResponse::new(
                positions.clone(),
                tets.clone(),
                vec![law(Some(ConductorMagnitude::CurrentDensity), 0.1, 1.0); 6],
                limits(),
                &mut || false,
            )
            .unwrap();
            assert!(matches!(
                limited.solve_candidate(
                    &boundary,
                    limits(),
                    controls(),
                    response_controls(),
                    &mut || false
                ),
                Err(SpatialConductionResponseError::Material(_))
            ));
            assert_eq!(
                limited.conductors()[0].material().query().axes()["operating"],
                0.0
            );
        }
    }
    assert_eq!(
        problem.conductors()[0].material().query().axes()["operating"],
        0.0
    );
}
