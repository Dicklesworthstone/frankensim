//! Prescribed baths act through the same live-contact/material/feedback dual.
use super::*;
use fs_conduction::adjoint::{MeanFluxFeedback, RobinResponse};
use fs_conduction::interface::temperature::steady::SteadyContactNewtonConfig;
use fs_conduction::{ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBc};

fn source_card() -> InterfaceSystemCard {
    card_with_response(PropertyValue::Surface2D {
        axes: ["TA".into(), "TB".into()],
        axis_dims: [fs_conduction::TEMPERATURE_DIMS; 2],
        samples: [[256., 256.], [384., 256.], [256., 384.], [384., 384.]]
            .into_iter().map(|p| (p, 0.5 + (p[0] - 320.) / 1024.
                + 3. * (p[1] - 320.) / 1024.)).collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    }, two_sided_validity())
}

fn boundary(mesh: &ConductionMesh, bath: f64, row: [f64; 3], bath_side_a: bool) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("bath", |f| if bath_side_a {
            f.centroid[2] < 0. && f.centroid[1] == 0.
        } else { f.centroid[2] > 0. && f.centroid[0] == 0. },
            ThermalBc::dirichlet(bath).unwrap()).unwrap()
        .region("film", |f| f.centroid[2] < 0. && f.centroid[0] == 0.,
            ThermalBc::robin(row[0], row[1]).unwrap()).unwrap()
        .region("loss", |f| f.centroid[2] > 0. && f.centroid[1] == 0.,
            ThermalBc::neumann(row[2]).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-12, max_iterations: 1024, restart: 8 }
}

fn config() -> SteadyContactNewtonConfig {
    let mut config = SteadyContactNewtonConfig::default();
    config.newton.absolute_tolerance = 1e-11;
    config.newton.forcing_minimum = 1e-12;
    config.newton.forcing_maximum = 1e-10;
    config.newton.linear_restart = 8;
    config.residual_tolerance_w = 1e-9;
    config.energy_tolerance_w = 1e-9;
    config.max_face_pairs = 1;
    config
}

// Independently solve physical boundary VALUES. Neither the new lift nor a
// derivative owner participates in these contact/material Newton solves.
fn physical(cx: &Cx<'_>, mesh: &ConductionMesh, material: &ConductivityModel,
    source: &ScalarField, card: &InterfaceSystemCard, bath: f64, feedback: bool,
    bath_side_a: bool,
) -> (ThermalBoundary, Vec<f64>) {
    let mut row = [2., 300., 0.];
    let mut initial = vec![320.; mesh.vertex_count()];
    for _ in 0..128 {
        let boundary = boundary(mesh, bath, row, bath_side_a);
        let set = two_sided_contacts(cx, mesh, &boundary, card, &two_sided_point());
        let problem = ConductionProblem { mesh, boundary: &boundary, material,
            element_materials: None, source };
        let solved = set.solve_steady_newton(cx, problem, &initial, config()).unwrap();
        if !feedback { return (boundary, solved.temperature); }
        // These exact faces are film={0,2,3}, loss={4,5,7}. A bath on
        // either side fixes contact nodes and one feedback patch's nodes.
        let film = (solved.temperature[0] + solved.temperature[2] + solved.temperature[3]) / 3.;
        let loss = (solved.temperature[4] + solved.temperature[5] + solved.temperature[7]) / 3.;
        let next = [2. + 0.01 * (film - 320.), 300. + 0.2 * (film - 320.),
            0.15 * (film - 320.) + 0.1 * (loss - 320.)];
        if row.iter().zip(next).all(|(a, b)| (a - b).abs() <= 1e-12) {
            return (boundary, solved.temperature);
        }
        for i in 0..3 { row[i] = 0.5 * (row[i] + next[i]); }
        initial = solved.temperature;
    }
    panic!("physical boundary fixed point did not converge");
}

#[test]
fn g3_prescribed_contact_lift_matches_fresh_bath_solves_with_all_boundary_feedback() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(256., 1.4), (384., 4.0)]).unwrap());
        let source = ScalarField::Nodal(vec![0., 0., 0., 0., 30., 30., 30., 30.]);
        let card = source_card();
        for bath_side_a in [false, true] { for with_feedback in [false, true] {
            let (boundary, temperature) = physical(cx, &mesh, &material, &source,
                &card, 310., with_feedback, bath_side_a);
            let before = temperature.clone();
            let set = two_sided_contacts(cx, &mesh, &boundary, &card, &two_sided_point());
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, element_materials: None, source: &source };
            let response = set.linearize_steady_endpoint(cx, problem, &temperature,
                linear(), 1e-9, 1e-9, &["film"], 1).unwrap();
            let flux = with_feedback.then(|| MeanFluxFeedback::new(cx, &mesh, &boundary,
                &["film", "loss"], &[0., 0., 0.15, 0.1], 32).unwrap());
            let h = if with_feedback { 0.01 } else { 0. };
            // Local and cross-reference arrays contribute once each.
            let reference = if with_feedback { 0.07 } else { 0. };
            let cross_reference = if with_feedback { 0.13 } else { 0. };
            let goal = if bath_side_a { 7 } else { 2 };
            let mut weights = [0.; 8]; weights[goal] = 1.;
            let full = response.pullback_with_robin_and_flux_feedback_and_resistance_controls(cx,
                &[h], &[reference], &[cross_reference], flux.as_ref(), &weights, 128, 1).unwrap();
            let lift = response.prescribed_temperature_pullback_with_flux_feedback(cx,
                &weights, &full.transport.nodal_load, &[h], &[reference], &[cross_reference],
                flux.as_ref(), 4096).unwrap();
            let derivative: f64 = lift.iter().sum();
            assert!(full.transport.relative_residual <= linear().tolerance);
            for epsilon in [0.001, 0.0005] {
                let plus = physical(cx, &mesh, &material, &source, &card,
                    310. + epsilon, with_feedback, bath_side_a).1;
                let minus = physical(cx, &mesh, &material, &source, &card,
                    310. - epsilon, with_feedback, bath_side_a).1;
                near(derivative, (plus[goal] - minus[goal]) / (2. * epsilon), 2e-6);
            }
            // Same complete dual, but deliberately omit ONLY the source-law
            // fixed columns. A frozen-contact lift is physically different.
            let frozen = RobinResponse::prescribed_temperature_pullback_with_flux_feedback(cx,
                problem, Some(response.contact().physical_contacts().frozen_interfaces()),
                &temperature, &weights, &full.transport.nodal_load, &["film"], &[h],
                &[reference], &[cross_reference], flux.as_ref(), 4096).unwrap();
            assert!((derivative - frozen.iter().sum::<f64>()).abs() > 1e-4,
                "the prescribed trace must exercise dR/dT, not only the full free-node dual");
            if with_feedback {
                let no_flux = response.prescribed_temperature_pullback_with_flux_feedback(cx,
                    &weights, &full.transport.nodal_load, &[h], &[reference], &[cross_reference],
                    None, 4096).unwrap();
                assert!((derivative - no_flux.iter().sum::<f64>()).abs() > 1e-3,
                    "flux means must retain their fixed-node weights");
                if bath_side_a {
                    let no_robin = response.prescribed_temperature_pullback_with_flux_feedback(cx,
                        &weights, &full.transport.nodal_load, &[0.], &[0.], &[0.],
                        flux.as_ref(), 4096).unwrap();
                    assert!((derivative - no_robin.iter().sum::<f64>()).abs() > 1e-4,
                        "Robin means must retain their fixed-node weights");
                }
            }
            for vertex in 0..mesh.vertex_count() {
                if !boundary.dirichlet().iter().any(|&(v, _)| v == vertex) {
                    assert_eq!(lift[vertex], 0.);
                }
            }
            assert_eq!(temperature, before, "a lift never replaces the accepted field");
        } }
    });
}

#[test]
fn g0_g1_prescribed_contact_goal_keeps_direct_term_and_original_refusals() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let card = source_card();
        let (boundary, temperature) = physical(cx, &mesh, &material, &source, &card, 310., false, true);
        let set = two_sided_contacts(cx, &mesh, &boundary, &card, &two_sided_point());
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source };
        let response = set.linearize_steady_endpoint(cx, problem, &temperature,
            linear(), 1e-9, 1e-9, &["film"], 1).unwrap();
        let fixed = boundary.dirichlet()[0].0;
        let mut weights = [0.; 8]; weights[fixed] = 1.;
        let full = response.pullback_with_robin_and_flux_feedback_and_resistance_controls(cx,
            &[0.], &[0.], &[0.], None, &weights, 128, 1).unwrap();
        assert!(full.transport.nodal_load.iter().all(|&v| v == 0.));
        let lift = response.prescribed_temperature_pullback_with_flux_feedback(cx,
            &weights, &full.transport.nodal_load, &[0.], &[0.], &[0.], None, 4096).unwrap();
        assert_eq!(lift, weights);
        assert_eq!(temperature[fixed], 310.);
        let epsilon = 0.001;
        let plus = physical(cx, &mesh, &material, &source, &card, 310. + epsilon, false, true).1;
        let minus = physical(cx, &mesh, &material, &source, &card, 310. - epsilon, false, true).1;
        near(lift[fixed], (plus[fixed] - minus[fixed]) / (2. * epsilon), 1e-9);
        let mut wrong_dual = full.transport.nodal_load.clone(); wrong_dual[fixed] = 1.;
        assert!(response.prescribed_temperature_pullback_with_flux_feedback(cx,
            &weights, &wrong_dual, &[0.], &[0.], &[0.], None, 4096).is_err());
        assert!(response.prescribed_temperature_pullback_with_flux_feedback(cx,
            &weights, &full.transport.nodal_load, &[0.], &[0.], &[0.], None, 1).is_err());
        assert!(response.prescribed_temperature_pullback_with_flux_feedback(cx,
            &weights[..7], &full.transport.nodal_load, &[0.], &[0.], &[0.], None, 4096).is_err());
        let gate = CancelGate::new_clock_free(); gate.request();
        context(&gate, |cancelled| {
            assert!(matches!(response.prescribed_temperature_pullback_with_flux_feedback(cancelled,
                &weights, &full.transport.nodal_load, &[0.], &[0.], &[0.], None, 4096),
                Err(ConductionError::Cancelled { .. })));
        });
    });
}
