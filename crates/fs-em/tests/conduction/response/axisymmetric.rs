//! G1/G3/G4: nonlinear meridional fields, exact swept heat and source closure.
//! Synthetic laws test the numerical model, not measured material validity.
use super::*;
use fs_em::conduction::{
    ConductionContact,
    axisymmetric::{AxisymmetricConduction, response::AxisymmetricConductionResponse},
    response::contact::ContactResponseError,
};

#[path = "../../common/electrical_interface_response.rs"]
mod interface_response;

fn axial_mesh(angle: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "two half-length cylinder layers with unequal swept triangles",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 0.5],
            [1.0, 0.5],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        vec![[0, 1, 3], [0, 3, 2], [2, 3, 5], [2, 5, 4]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

fn axial_specimen(angle: f64) -> AxisymmetricConductionResponse {
    AxisymmetricConductionResponse::new(
        axial_mesh(angle),
        [1.0, 1.0, 3.0, 3.0]
            .into_iter()
            .zip([0.2, 0.2, 0.4, 0.4])
            .map(|(a, b)| {
                response(
                    ElectricalPropertyBasis::Conductivity,
                    ConductorMagnitude::ElectricField,
                    a,
                    b,
                    10.0,
                )
            })
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn contact_mesh(angle: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "two unit-length cylinder layers with a swept disk joint",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}

#[test]
fn actual_nonlinear_series_fields_keep_swept_heat_and_radial_current_loads() {
    let (density, expected_e) = series_current([1.0, 3.0], [0.2, 0.4], 1.5);
    for angle in [0.7, 2.0, core::f64::consts::TAU] {
        let area = angle / 2.0;
        let model = axial_specimen(angle);
        let before = model.conductors().to_vec();
        // A nonlinear source cannot escape through fixed Ohmic admission.
        assert!(
            AxisymmetricConduction::new(
                model.mesh().clone(),
                before.clone(),
                limits(),
                &mut || false,
            )
            .is_err()
        );
        for sign in [-1.0, 1.0] {
            let voltage = [(0, 0.0), (1, 0.0), (4, sign * 1.5), (5, sign * 1.5)];
            let traces = [ConductionCurrentTrace {
                vertices: [4, 5],
                inward_current_density_a_m2: sign * density,
            }];
            let voltage_result = model
                .solve_candidate(&voltage, &[], limits(), controls(), policy(), &mut || false)
                .unwrap();
            let current_result = model
                .solve_current_candidate(
                    &[4, 5],
                    &[0, 1],
                    sign * area * density,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            let trace_result = model
                .solve_candidate(
                    &[(0, 0.0), (1, 0.0)],
                    &traces,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            near(
                trace_result.report().nodal_applied_current_a[4].unwrap(),
                sign * area * density / 3.0,
            );
            near(
                trace_result.report().nodal_applied_current_a[5].unwrap(),
                sign * 2.0 * area * density / 3.0,
            );
            for result in [&voltage_result, &current_result, &trace_result] {
                let r = result.report();
                assert!(r.field_solves > 1);
                near(r.voltage_v[4], sign * 1.5);
                near(r.voltage_v[2], sign * 0.5 * expected_e[0]);
                for triangle in 0..4 {
                    let e = expected_e[triangle / 2];
                    near(r.electric_field_v_m[triangle][0], 0.0);
                    near(r.electric_field_v_m[triangle][1], -sign * e);
                    near(r.current_density_a_m2[triangle][1], -sign * density);
                    near(
                        r.cell_power_w[triangle],
                        model.mesh().face_measure(triangle).unwrap() * density * e,
                    );
                    near(
                        r.evaluated[triangle]
                            .conductor
                            .material()
                            .state_coordinate("operating")
                            .unwrap(),
                        e,
                    );
                    assert_eq!(
                        r.evaluated[triangle]
                            .conductor
                            .source_property()
                            .answer()
                            .receipt
                            .selected,
                        before[triangle].source_property().answer().receipt.selected,
                    );
                }
                near(r.power_w, area * density * 1.5);
                near(r.power_balance_residual_w, 0.0);
                near(r.boundary_current_a, 0.0);
            }
            let (continued, _) = current_result.into_parts();
            let next = continued
                .solve_current_candidate(
                    &[4, 5],
                    &[0, 1],
                    sign * area * density,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            near(next.report().power_w, area * density * 1.5);
        }
        assert_eq!(model.conductors(), before);
    }
}

#[test]
fn radial_nonlinear_current_converges_to_the_independent_coaxial_law() {
    let (angle, length, inner_radius, outer_radius) = (1.3, 0.8, 0.5, 1.5);
    let (a, b, flux) = (1.3_f64, 0.4_f64, 1.0_f64);
    let current = angle * length * flux;
    // Cylindrical continuity gives J(r)=flux/r and J=aE+bE^2.
    // Integrate the independent positive root, without a field interpolation.
    let field = |r: f64| {
        let density = flux / r;
        2.0 * density / (a + (a * a + 4.0 * b * density).sqrt())
    };
    let panels = 4096;
    let dr = (outer_radius - inner_radius) / panels as f64;
    let exact_voltage = (0..=panels)
        .map(|i| {
            let weight = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            weight * field(inner_radius + i as f64 * dr)
        })
        .sum::<f64>()
        * dr
        / 3.0;
    let law = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        a,
        b,
        10.0,
    );
    let cfg = ConductionControls {
        linear: fs_solver::CheckedCgConfig {
            tolerance: 1e-13,
            ..config()
        },
        current_tolerance_a: 1e-9,
        power_tolerance_w: 1e-9,
    };
    let mut previous = f64::INFINITY;
    for n in [4, 8, 16] {
        let mesh =
            super::super::axisymmetric::rectangle(inner_radius, outer_radius, length, n, n, angle);
        let count = mesh.faces().len();
        let model = AxisymmetricConductionResponse::new(
            mesh,
            vec![law.clone(); count],
            limits(),
            &mut || false,
        )
        .unwrap();
        let inner: Vec<_> = (0..=n).map(|j| j * (n + 1)).collect();
        let outer: Vec<_> = inner.iter().map(|v| v + n).collect();
        let solved = model
            .solve_current_candidate(
                &inner,
                &outer,
                current,
                limits(),
                cfg,
                policy(),
                &mut || false,
            )
            .unwrap();
        let voltage = solved.report().voltage_v[inner[0]];
        let error = (voltage / exact_voltage - 1.0).abs();
        assert!(
            error < 0.6 * previous,
            "n={n}, error={error}, previous={previous}"
        );
        previous = error;
        near(solved.report().power_w, current * voltage);
        near(
            solved.report().cell_power_w.iter().sum::<f64>(),
            solved.report().power_w,
        );
        near(solved.report().power_balance_residual_w, 0.0);
        assert!(
            solved
                .report()
                .evaluated
                .iter()
                .all(|sample| sample.conductivity_defect_s_m.abs() < 1e-10
                    && sample.power_defect_w.abs() < 1e-9)
        );
    }
    assert!(previous < 3e-3);
}

#[test]
fn actual_power_and_relaxed_current_close_original_laws_without_scaling_accepted_fields() {
    let angle = 1.3;
    let area = angle / 2.0;
    let bulk = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.2,
        20.0,
    );
    let model = AxisymmetricConductionResponse::new(
        axial_mesh(angle),
        vec![bulk; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    // Independent uniform law: E=2, sigma=1.4, J=2.8, P=A*E*J.
    let supply = ConductionPowerSource {
        power_w: area * 5.6,
        maximum_voltage_v: 2.01,
        maximum_current_a: area * 2.81,
    };
    for reverse in [false, true] {
        let (positive, negative, sign) = if reverse {
            (&[0, 1][..], &[4, 5][..], -1.0)
        } else {
            (&[4, 5][..], &[0, 1][..], 1.0)
        };
        let next = model
            .solve_power_candidate(
                positive,
                negative,
                supply,
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .unwrap();
        near(next.report().power_w, supply.power_w);
        near(
            next.report().voltage_v[positive[0]] - next.report().voltage_v[negative[0]],
            2.0,
        );
        for triangle in 0..4 {
            near(next.report().electric_field_v_m[triangle][1], -sign * 2.0);
            near(next.report().current_density_a_m2[triangle][1], -sign * 2.8);
            near(
                next.report().cell_power_w[triangle],
                model.mesh().face_measure(triangle).unwrap() * 5.6,
            );
        }
    }
    assert!(
        model
            .solve_power_candidate(
                &[4, 5],
                &[0, 1],
                ConductionPowerSource {
                    maximum_voltage_v: 1.5,
                    ..supply
                },
                limits(),
                controls(),
                policy(),
                &mut || false,
            )
            .is_err()
    );
    let steep = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        0.01,
        100.0,
        200.0,
    );
    let model = AxisymmetricConductionResponse::new(
        axial_mesh(angle),
        vec![steep; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let relaxed = model.clone().with_relaxation(0.5).unwrap();
    let solve = |model: &AxisymmetricConductionResponse| {
        model.solve_current_candidate(
            &[4, 5],
            &[0, 1],
            area,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
    };
    assert!(matches!(
        solve(&model),
        Err(ConductionResponseError::NotConverged { .. })
    ));
    let expected_e = 2.0 / (0.01 + (0.0001_f64 + 400.0).sqrt());
    let solved = solve(&relaxed).unwrap();
    near(solved.report().voltage_v[4], expected_e);
    near(solved.report().power_w, area * expected_e);
    assert_eq!(solved.problem().relaxation(), 0.5);
    assert!(solved.report().field_solves < 32);
    let adaptive = AdaptiveConductionRelaxation {
        initial_weight: 0.5,
        minimum_weight: 0.02,
        maximum_weight: 1.0,
    };
    let adapted = relaxed.with_adaptive_relaxation(adaptive).unwrap();
    let solved = solve(&adapted).unwrap();
    near(solved.report().voltage_v[4], expected_e);
    assert_eq!(solved.problem().adaptive_relaxation(), Some(adaptive));
}

#[test]
fn nonlinear_swept_contact_uses_actual_local_sources_and_separate_interface_heat() {
    let beta = 0.75;
    for angle in [0.7, core::f64::consts::TAU] {
        let area = angle / 2.0;
        let bulk = response(
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            1.0,
            0.0,
            20.0,
        );
        let law = interface_response::law(bulk.original(), bulk.original(), beta);
        let contact = ConductionContact {
            first_vertices: [2, 3],
            second_vertices: [5, 4],
            source: law.original().clone(),
        };
        assert!(
            AxisymmetricConductionResponse::new_with_phase_contacts(
                contact_mesh(angle),
                vec![bulk.clone(); 4],
                vec![contact.clone()],
                vec![[None, None]],
                limits(),
                &mut || false,
            )
            .is_err()
        );
        let model = AxisymmetricConductionResponse::new_with_contact_response(
            contact_mesh(angle),
            vec![bulk; 4],
            vec![contact],
            vec![law],
            interface_response::controls(),
            limits(),
            &mut || false,
        )
        .unwrap()
        .with_relaxation(0.5)
        .unwrap();
        for temperature in [300.0, 350.0] {
            let law = &model.contact_response_sources().unwrap().laws()[0];
            let mut contacts = model.contacts().to_vec();
            contacts[0].source = law
                .at_query_point(
                    &contacts[0]
                        .source
                        .resolved()
                        .query()
                        .with("T", temperature)
                        .unwrap(),
                )
                .unwrap();
            for sign in [-1.0, 1.0] {
                let density = interface_response::density(2.0, temperature, beta);
                let rc = 0.5 + 0.002 * (temperature - 300.0) + beta * density;
                let voltage = [(0, 0.0), (1, 0.0), (6, 2.0 * sign), (7, 2.0 * sign)];
                let current = model
                    .solve_current_with_materials_and_contacts_candidate(
                        model.conductors().to_vec(),
                        contacts.clone(),
                        &[6, 7],
                        &[0, 1],
                        sign * area * density,
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                    .unwrap();
                let voltage_driven = model
                    .solve_with_materials_and_contacts_candidate(
                        model.conductors().to_vec(),
                        contacts.clone(),
                        &voltage,
                        &[],
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                    .unwrap();
                for solved in [&current, &voltage_driven] {
                    let r = solved.report();
                    near(r.voltage_v[6], 2.0 * sign);
                    near(r.power_w, 2.0 * area * density);
                    near(r.contact_response[0].current_a, -sign * area * density);
                    near(r.contact_response[0].power_w, area * density * density * rc);
                    for jump in r.contact_response[0].voltage_jump_v {
                        near(jump, -sign * density * rc);
                    }
                    for triangle in 0..4 {
                        near(
                            r.cell_power_w[triangle],
                            model.mesh().face_measure(triangle).unwrap() * density * density,
                        );
                    }
                    near(
                        r.cell_power_w.iter().sum::<f64>() + r.contact_response[0].power_w,
                        r.power_w,
                    );
                    for evaluated in &r.evaluated_contacts[0].evaluated {
                        near(
                            evaluated
                                .resolved()
                                .state_coordinate("normal_current")
                                .unwrap(),
                            density,
                        );
                        near(
                            evaluated.resolved().state_coordinate("T").unwrap(),
                            temperature,
                        );
                        assert_eq!(
                            evaluated.resolved().properties()[0]
                                .answer()
                                .receipt
                                .selected,
                            model.contacts()[0].source.resolved().properties()[0]
                                .answer()
                                .receipt
                                .selected,
                        );
                    }
                }
            }
        }
        assert!(matches!(
            model.solve_candidate(
                &[(0, 0.0), (1, 0.0), (6, 3.0), (7, 1.0)],
                &[],
                limits(),
                controls(),
                policy(),
                &mut || false,
            ),
            Err(ConductionResponseError::Contact(
                ContactResponseError::NonuniformTrace { .. }
            ))
        ));
    }
}

#[test]
fn native_bounds_refuse_metric_axis_domain_and_cancelled_candidates_without_publication() {
    let model = axial_specimen(2.0);
    assert!(
        PlanarConductionResponse::new(
            model.mesh().clone(),
            model.sources().to_vec(),
            limits(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        AxisymmetricConductionResponse::new(
            layered(0.5, 1.0, 1.0),
            model.sources().to_vec(),
            limits(),
            &mut || false,
        )
        .is_err()
    );
    let boundary = [(0, 0.0), (1, 0.0), (4, 1.5), (5, 1.5)];
    let mut polls = 0;
    let baseline = model
        .solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let before = model.conductors().to_vec();
    let mut interrupted = 0;
    assert!(matches!(
        model.solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            interrupted += 1;
            interrupted == polls / 2
        }),
        Err(ConductionResponseError::Cancelled)
    ));
    assert!(matches!(
        model.solve_candidate(
            &boundary,
            &[],
            limits(),
            controls(),
            ConductionResponseControls {
                maximum_iterations: 1,
                ..policy()
            },
            &mut || false,
        ),
        Err(ConductionResponseError::NotConverged { .. })
    ));
    let retry = model
        .solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(baseline.report().voltage_v, retry.report().voltage_v);
    assert_eq!(baseline.report().cell_power_w, retry.report().cell_power_w);
    assert_eq!(model.conductors(), before);
    let bounded = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.2,
        0.1,
    );
    let bounded = AxisymmetricConductionResponse::new(
        axial_mesh(2.0),
        vec![bounded; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        bounded.solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            false
        }),
        Err(ConductionResponseError::Material(_))
    ));
    // Two physical disk edges join these two terminal sets. The zero-area
    // axis must not borrow either disk from the distinct outer electrode.
    let mesh = TriComplex2::from_indexed_triangles(
        "axis-only terminal beside positive-area disk edges",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        vec![[0, 1, 2], [0, 2, 3]],
        Metric2::axisymmetric(2.0).unwrap(),
    )
    .unwrap();
    let axis = AxisymmetricConductionResponse::new(
        mesh,
        vec![model.sources()[0].clone(); 2],
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(
        axis.solve_current_with_materials_and_contacts_candidate(
            axis.conductors().to_vec(),
            vec![],
            &[0, 3],
            &[1, 2],
            1.0,
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .is_err()
    );
    assert!(
        axis.solve_power_with_materials_and_contacts_candidate(
            axis.conductors().to_vec(),
            vec![],
            &[0, 3],
            &[1, 2],
            ConductionPowerSource {
                power_w: 1.0,
                maximum_voltage_v: 10.0,
                maximum_current_a: 10.0,
            },
            limits(),
            controls(),
            policy(),
            &mut || false,
        )
        .is_err()
    );
}
