//! Original pinned material claims coupled through the resolved thermal node.
use super::*;
use fs_couple::tribothermal::partial_slip::resolved::*;

fn resolved_controls() -> ResolvedCaloricContactSolve {
    let mut contact = controls();
    contact.newton.linear_restart = 3;
    ResolvedCaloricContactSolve {
        contact,
        temperature_scale_j_per_k: 1.0,
        interface_energy_tolerance_j: 2e-8,
    }
}
fn resolved_source(
    state: &DepositionSimulation<'_>,
    _: f64,
) -> Result<ResolvedCaloricContactSource<()>, Infallible> {
    Ok(ResolvedCaloricContactSource {
        conductivity_w_m_k: [3.0, 1.0],
        transport: other(state)?,
    })
}

#[test]
fn pinned_mu_uses_actual_body_temperatures_and_network_derived_heat() {
    let (a, ca) = common::material(1.0, 0.05, 0.0);
    let (b, cb) = common::material(1.0, 0.1, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = source(&cards, [(200.0, 0.1), (600.0, 0.5)], 600.0);
    let thermal = initial(&curves, [300.0; 2], 1.5);
    let original = loaded(
        &thermal,
        &curves,
        data(&cards, &curves, &card, &point),
        0.2,
        [[0.0; 2]; 2],
        1.0,
    );
    let old = original
        .clone()
        .with_resolved_interface(&thermal, 7, [3.0, 1.0])
        .unwrap();
    assert_eq!(original.checkpoint(), old.checkpoint());
    let step = old
        .step_resolved_candidate(
            &thermal,
            1.0,
            frame(),
            motion(1.0),
            &work("resolved material"),
            resolved_controls(),
            &[0.0; 2],
            300.0,
            resolved_source,
            || false,
        )
        .unwrap();
    // Independent three-node conduction reduction plus the original two-slider
    // work polynomial. Source mu=.2+.001*(Ta-300); its Ta axis remains Ta.
    let y = 6.0 * 0.05 / (0.05 + 6.0) + 2.0 * 0.1 / (0.1 + 2.0);
    let rise_a_per_j = (6.0 / (0.05 + 6.0)) / y;
    let heat = root(0.0, 4.0, |d| {
        let mu = 0.2 + 0.001 * rise_a_per_j * d;
        d - 10.0 * mu + 0.292 * mu * mu
    });
    let mu = 0.2 + 0.001 * rise_a_per_j * heat;
    let r = step.report();
    near(r.caloric.contact.dissipation.heat_j, heat);
    near(r.interface_temperature_k, 300.0 + heat / y);
    near(
        step.contact().temperature_k()[0],
        300.0 + rise_a_per_j * heat,
    );
    near(
        step.contact().temperature_k()[1],
        300.0 + (2.0 / 2.1) * heat / y,
    );
    near(r.caloric.contact.tangent_force_n[0], -10.0 * mu);
    near(
        r.caloric.contact.generalized_work.numerical_dissipation_j,
        0.5 * 0.292 * mu * mu,
    );
    near(r.caloric.thermal.stored_energy_change_j, heat);
    near(r.caloric.total_energy_residual_j, 0.0);
    let sample = r.caloric.material.as_ref().unwrap();
    near(sample.partition.surface_a(), 0.75);
    near(sample.partition.surface_b(), 0.25);
    near(
        sample.resolved.state_coordinate("Ta").unwrap(),
        step.contact().temperature_k()[0],
    );
    near(
        sample.resolved.state_coordinate("Tb").unwrap(),
        step.contact().temperature_k()[1],
    );
    assert!((r.interface_temperature_k - step.contact().temperature_k()[0]).abs() > 0.01);
    let MaterialPropertySelection::PinnedByProperty(pins) =
        &old.material_friction_binding().unwrap().selection
    else {
        panic!("original claims unpinned")
    };
    for property in sample.resolved.properties() {
        assert_eq!(
            pins.iter()
                .find(|(name, _)| name == property.requirement().name())
                .unwrap()
                .1,
            property.answer().receipt.selected
        );
    }
    near(sample.speed_mps, 1.0);
    near(sample.accumulated_distance_m, 1.0);
    assert_eq!(step.contact().checkpoint().state().accepted_steps(), 8);
}

#[test]
fn resolved_source_domain_refusal_preserves_thermal_and_loaded_contact() {
    let (a, ca) = common::material(1.0, 0.05, 0.0);
    let (b, cb) = common::material(1.0, 0.1, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = source(&cards, [(200.0, 0.1), (600.0, 0.5)], 305.0);
    let thermal = initial(&curves, [300.0; 2], 1.5);
    let old = loaded(
        &thermal,
        &curves,
        data(&cards, &curves, &card, &point),
        0.2,
        [[0.0; 2]; 2],
        1.0,
    )
    .with_resolved_interface(&thermal, 7, [3.0, 1.0])
    .unwrap();
    let checkpoint = old.checkpoint().clone();
    let cells = thermal.cells().clone();
    let calls = std::cell::Cell::new(0usize);
    let refused = old.step_resolved_candidate(
        &thermal,
        1.0,
        frame(),
        motion(1.0),
        &work("resolved material domain"),
        resolved_controls(),
        &[0.0; 2],
        300.0,
        |state, t| {
            calls.set(calls.get() + 1);
            resolved_source(state, t)
        },
        || false,
    );
    assert!(
        calls.get() > 0,
        "initial source query must pass before later refusal"
    );
    assert!(matches!(
        refused,
        Err(CaloricContactError::Material(SlidingError::Material(_)))
    ));
    assert_eq!(old.checkpoint(), &checkpoint);
    assert_eq!(thermal.cells(), &cells);
    near(old.time_s(), 0.0);
    near(thermal.time_s(), 0.0);
}
