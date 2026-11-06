//! Native bulk/contact current and heat against an independently solved joint.
use super::*;
#[path = "../../../common/electrical_interface_response.rs"]
mod interface_response;
use fs_em::conduction::response::contact::ContactResponseError;

fn specimen() -> PlanarConductionResponse {
    let bulk = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        1.,
        0.,
        20.,
    );
    let law = interface_response::law(bulk.original(), bulk.original(), 0.75);
    let contact = ConductionContact {
        first_vertices: [1, 2],
        second_vertices: [4, 7],
        source: law.original().clone(),
    };
    // A source with unclosed electrical axes is never silently admitted as a
    // static interface merely because one coefficient was already resolved.
    assert!(
        PlanarConductionResponse::new_with_phase_contacts(
            mesh(0.4),
            vec![bulk.clone(); 4],
            vec![contact.clone()],
            vec![[None, None]],
            limits(),
            &mut || false,
        )
        .is_err()
    );
    PlanarConductionResponse::new_with_contact_response(
        mesh(0.4),
        vec![bulk; 4],
        vec![contact],
        vec![law],
        interface_response::controls(),
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap()
}

#[test]
fn nonlinear_contact_refreshes_actual_current_and_temperature_with_one_surface_heat_account() {
    let model = specimen();
    for temperature in [300., 350.] {
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
        for sign in [-1., 1.] {
            let voltage = 2. * sign;
            let j = interface_response::density(voltage, temperature, 0.75);
            let rc = 0.5 + 0.002 * (temperature - 300.) + 0.75 * j;
            let boundary = [(0, 0.), (3, 0.), (5, voltage), (6, voltage)];
            let solved = model
                .solve_with_materials_and_contacts_candidate(
                    model.conductors().to_vec(),
                    contacts.clone(),
                    &boundary,
                    &[],
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            let r = solved.report();
            assert!(r.field_solves > 1);
            near(r.power_w, 0.4 * j * voltage.abs());
            near(r.contact_response[0].power_w, 0.4 * j * j * rc);
            near(r.contact_response[0].current_a, -sign * 0.4 * j);
            for v in r.contact_response[0].voltage_jump_v {
                near(v, -sign * j * rc);
            }
            near(r.cell_power_w.iter().sum::<f64>(), 0.8 * j * j);
            near(r.power_balance_residual_w, 0.);
            assert_eq!(r.evaluated_contacts.len(), 1);
            for source in &r.evaluated_contacts[0].evaluated {
                near(
                    source
                        .resolved()
                        .state_coordinate("normal_current")
                        .unwrap(),
                    j,
                );
                assert_eq!(source.resolved().state_coordinate("T"), Some(temperature));
                assert_eq!(
                    source.resolved().properties()[0].answer().receipt.selected,
                    model.contacts()[0].source.resolved().properties()[0]
                        .answer()
                        .receipt
                        .selected
                );
            }
            let current = model
                .solve_current_with_materials_and_contacts_candidate(
                    model.conductors().to_vec(),
                    contacts.clone(),
                    &[5, 6],
                    &[0, 3],
                    sign * 0.4 * j,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            near(current.report().voltage_v[5], voltage);
            near(current.report().power_w, r.power_w);
            let expected_power = r.power_w;
            let (native, sources) = solved.into_parts().0.into_frozen_parts();
            assert_eq!(sources.relaxation(), 0.5);
            assert_eq!(sources.contacts().unwrap().laws().len(), 1);
            let resumed =
                PlanarConductionResponse::from_frozen_parts(native, sources, &mut || false)
                    .unwrap();
            let again = resumed
                .solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
                    false
                })
                .unwrap();
            near(again.report().power_w, expected_power);
        }
    }
    assert_eq!(
        model.contacts()[0]
            .source
            .resolved()
            .state_coordinate("normal_current"),
        Some(0.)
    );
}

#[test]
fn nonuniform_contact_never_becomes_an_averaged_law_and_cancellation_preserves_source_state() {
    let model = specimen();
    let varying = [(0, 0.), (3, 0.), (5, 1.), (6, 3.)];
    assert!(matches!(
        model.solve_candidate(&varying, &[], limits(), controls(), policy(), &mut || false),
        Err(ConductionResponseError::Contact(
            ContactResponseError::NonuniformTrace { .. }
        ))
    ));
    let boundary = [(0, 0.), (3, 0.), (5, 2.), (6, 2.)];
    let mut calls = 0;
    let accepted = model
        .solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut polled = 0;
    assert!(matches!(
        model.solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            polled += 1;
            polled == calls / 2
        }),
        Err(ConductionResponseError::Cancelled)
    ));
    let retry = model
        .solve_candidate(&boundary, &[], limits(), controls(), policy(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(accepted.report().voltage_v, retry.report().voltage_v);
    assert_eq!(
        accepted.report().contact_response,
        retry.report().contact_response
    );
    let limited = ConductionResponseControls {
        maximum_iterations: 1,
        ..policy()
    };
    assert!(
        matches!(model.solve_candidate(&boundary,&[],limits(),controls(),limited,&mut ||false),
        Err(ConductionResponseError::NotConverged{maximum_contact_resistance_defect_ohm_m2,..})
            if maximum_contact_resistance_defect_ohm_m2>0.)
    );
    let original = &model.contacts()[0].source;
    let law = &model.contact_response_sources().unwrap().laws()[0];
    assert!(law.at_operating_point(original, 1., 11.).is_err());
    let warmer = law
        .at_query_point(&original.resolved().query().with("T", 350.).unwrap())
        .unwrap();
    assert!(law.relax_operating_point(original, &warmer, 0.5).is_err());
}

#[test]
fn fixed_contact_bindings_retain_nonuniform_native_p1_current_and_heat() {
    let bulk = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        1.,
        0.,
        20.,
    );
    let law = interface_response::fixed_law(bulk.original(), bulk.original());
    let contact = ConductionContact {
        first_vertices: [1, 2],
        second_vertices: [4, 7],
        source: law.original().clone(),
    };
    let fixed = PlanarConductionResponse::new_with_phase_contacts(
        mesh(0.4),
        vec![bulk.clone(); 4],
        vec![contact.clone()],
        vec![[None, None]],
        limits(),
        &mut || false,
    )
    .unwrap();
    let bound = PlanarConductionResponse::new_with_contact_response(
        mesh(0.4),
        vec![bulk; 4],
        vec![contact],
        vec![law],
        interface_response::controls(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let varying = [(0, 0.), (3, 0.), (5, 1.), (6, 3.)];
    let original = fixed
        .solve_candidate(&varying, &[], limits(), controls(), policy(), &mut || false)
        .unwrap();
    let recovered = bound
        .solve_candidate(&varying, &[], limits(), controls(), policy(), &mut || false)
        .unwrap();
    assert_eq!(original.report().voltage_v, recovered.report().voltage_v);
    assert_eq!(
        original.report().contact_response,
        recovered.report().contact_response
    );
    let jump = recovered.report().contact_response[0].voltage_jump_v;
    assert!((jump[1] - jump[0]).abs() > 0.01);
    near(recovered.report().power_balance_residual_w, 0.);
}
