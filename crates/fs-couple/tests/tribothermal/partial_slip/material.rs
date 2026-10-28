//! Original pinned interface queries drive the loaded native return map.
use super::super::{binding, context, interface_with_static};
use super::*;
use fs_couple::tribothermal::material::{KINETIC_FRICTION_PROPERTY, STATIC_FRICTION_PROPERTY};
use fs_couple::tribothermal::{MaterialSlidingContact, SlidingError, SlidingInterface};
use fs_matdb::{ClaimSet, InterfaceSystemCard, MaterialCard, QueryPoint, SystemContext};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Temperature;
use fs_tribo::InputAuthority;

fn source(cards: &[MaterialCard; 2], knots: [(f64, f64); 2], maximum: f64) -> InterfaceSystemCard {
    let base = interface_with_static(cards, 1.0, "Ta", Temperature::DIMS, knots, maximum);
    let kinetic = base
        .claims()
        .claims_ordered()
        .find(|(_, c)| c.key.name() == KINETIC_FRICTION_PROPERTY)
        .unwrap()
        .1
        .value
        .clone();
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == STATIC_FRICTION_PROPERTY {
            claim.value = kinetic.clone();
            claim.interpolation = fs_matdb::InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic thermal friction".into(),
            history: base.history().into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn data<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    source: &'a InterfaceSystemCard,
    point: &'a QueryPoint,
) -> SlidingInterface<'a> {
    let mut data = binding(
        cards,
        curves,
        source,
        point,
        HeatPartition::new(0.25, 0.75, 0.0).unwrap(),
    );
    data.accumulated_distance_axis = None;
    data.outside_destination = None;
    data
}
fn loaded<'a>(
    thermal: &DepositionSimulation<'_>,
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    data: SlidingInterface<'a>,
    mu: f64,
    strains: [[f64; 2]; 2],
    initial_speed: f64,
) -> CaloricPartialSlip<'a> {
    let sampler = MaterialSlidingContact::new(
        thermal,
        data.clone(),
        10.0,
        initial_speed,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let bound = sampler.interface().interface();
    let checkpoint = PartialSlipCheckpoint::new(
        NormalPatchView::new(
            "patch",
            "held patch",
            "synthetic patch",
            NormalPatchAuthority::SyntheticFixture,
            10.0,
            0.1,
            0.1,
            0.005,
        )
        .unwrap(),
        PartialSlipInterface::new(
            bound.ordered_system_id(),
            bound.history_id(),
            bound.provenance().source_id(),
            NormalPatchAuthority::SyntheticFixture,
        )
        .unwrap(),
        model(mu, zero_slopes()).reference_law().clone(),
        PartialSlipState::from_region_checkpoint(strains, [0.0; 2], [false; 2], Some(frame()), 7)
            .unwrap(),
    )
    .unwrap();
    CaloricPartialSlip::new(
        thermal,
        model(mu, zero_slopes()),
        checkpoint,
        [curves[0].curve(), curves[1].curve()],
        DESTINATIONS,
        HeatPartition::new(0.25, 0.75, 0.0).unwrap(),
    )
    .unwrap()
    .with_material_friction(
        thermal,
        data,
        initial_speed,
        InputAuthority::SyntheticFixture,
    )
    .unwrap()
}
fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (lo + hi);
        if residual(middle) > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn source_mu_temperature_feedback_matches_independent_two_slider_equation() {
    let (a, ca) = common::material(1.0, 0.05, 0.0);
    let (b, cb) = common::material(1.0, 0.1, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = source(&cards, [(200.0, 0.1), (600.0, 0.5)], 600.0);
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = loaded(
        &thermal,
        &curves,
        data(&cards, &curves, &card, &point),
        0.2,
        [[0.0; 2]; 2],
        1.0,
    );
    let mut drive = motion(1.0);
    drive.creepage = [0.3, 0.4];
    drive.rolling_speed_mps = 2.0;
    let step = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            drive,
            &work("source slide"),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    // Ci=N*mu*(0.15,0.85), ki=(750,250), unit travel. Physical heat is
    // D=10*mu-0.292*mu^2; thermal/source feedback gives mu=0.2+0.005*D.
    let heat = root(0.0, 4.0, |d| {
        let mu = 0.2 + 0.005 * d;
        d - 10.0 * mu + 0.292 * mu * mu
    });
    let mu = 0.2 + 0.005 * heat;
    near(step.report().contact.dissipation.heat_j, heat);
    near(step.contact().temperature_k()[0], 300.0 + 5.0 * heat);
    near(step.contact().temperature_k()[1], 300.0 + 7.5 * heat);
    near(step.report().contact.tangent_force_n[0], -0.6 * 10.0 * mu);
    near(step.report().contact.tangent_force_n[1], -0.8 * 10.0 * mu);
    near(
        step.report()
            .contact
            .generalized_work
            .numerical_dissipation_j,
        0.5 * 0.292 * mu * mu,
    );
    near(step.report().total_energy_residual_j, 0.0);
    let sample = step.report().material.as_ref().unwrap();
    near(sample.speed_mps, 1.0);
    near(sample.accumulated_distance_m, 1.0);
    near(sample.resolved.state_coordinate("speed").unwrap(), 1.0);
    near(
        sample.resolved.state_coordinate("Ta").unwrap(),
        step.contact().temperature_k()[0],
    );
    near(
        sample.resolved.state_coordinate("Tb").unwrap(),
        step.contact().temperature_k()[1],
    );
    let MaterialPropertySelection::PinnedByProperty(pins) =
        &old.material_friction_binding().unwrap().selection
    else {
        panic!("original claims not pinned")
    };
    for p in sample.resolved.properties() {
        assert_eq!(
            pins.iter()
                .find(|(name, _)| name == p.requirement().name())
                .unwrap()
                .1,
            p.answer().receipt.selected
        );
    }
    assert_eq!(step.contact().checkpoint().state().accepted_steps(), 8);
    assert_eq!(old.checkpoint().state().accepted_steps(), 7);
    assert_eq!(
        old.checkpoint().law().source_id(),
        step.contact().checkpoint().law().source_id()
    );
    assert_ne!(
        old.checkpoint().law().source_id(),
        old.checkpoint().interface().source_id()
    );
}

#[test]
fn source_thermal_weakening_releases_retained_strain_as_physical_heat_at_zero_motion() {
    let (a, ca) = common::material(1.0, 0.1, 0.0);
    let (b, cb) = common::material(1.0, 0.2, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = source(&cards, [(200.0, 0.9), (600.0, 0.5)], 600.0);
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let mu0 = 0.8;
    let strain = [[1.5 * mu0 / 750.0, 0.0], [8.5 * mu0 / 250.0, 0.0]];
    let old = loaded(
        &thermal,
        &curves,
        data(&cards, &curves, &card, &point),
        mu0,
        strain,
        0.0,
    );
    let step = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.0),
            &work("heated relaxation"),
            controls(),
            &[0.0, 10.0],
            |_| {
                Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::from([(DESTINATIONS[0], 10.0)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || false,
        )
        .unwrap();
    // Heating alone would give Ta=400 K and mu=.7. Released slip heat adds
    // Ta-400=.25D/.1; hence mu=.7-.0025D and D=.292*mu*(.8-mu).
    let heat = root(0.0, 0.1, |d| {
        let mu = 0.7 - 0.0025 * d;
        d - 0.292 * mu * (mu0 - mu)
    });
    let mu = 0.7 - 0.0025 * heat;
    near(step.report().contact.dissipation.heat_j, heat);
    near(
        step.report().contact.generalized_work.work_into_interface_j,
        0.0,
    );
    near(
        step.report().contact.generalized_work.coefficient_work_j,
        0.0,
    );
    near(
        step.report()
            .contact
            .generalized_work
            .numerical_dissipation_j,
        0.5 * 0.292 * (mu0 - mu).powi(2),
    );
    near(step.report().thermal.stored_energy_change_j, 10.0 + heat);
    near(
        step.report().contact_internal_energy_change_j,
        0.5 * 0.292 * (mu * mu - mu0 * mu0),
    );
    near(step.contact().temperature_k()[0], 400.0 + 2.5 * heat);
    near(step.contact().temperature_k()[1], 300.0 + 3.75 * heat);
    let next = step.contact().checkpoint().state();
    near(next.region_elastic_displacement_m()[0][0], 1.5 * mu / 750.0);
    near(next.region_elastic_displacement_m()[1][0], 8.5 * mu / 250.0);
    assert_eq!(next.sliding_regions(), [true; 2]);
    assert_eq!(
        old.checkpoint().state().region_elastic_displacement_m(),
        strain
    );
    near(step.report().total_energy_residual_j, 0.0);
}

#[test]
fn original_source_domain_failure_rolls_back_loaded_and_thermal_histories() {
    let (a, ca) = common::material(1.0, 0.05, 0.0);
    let (b, cb) = common::material(1.0, 0.1, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = source(&cards, [(200.0, 0.1), (600.0, 0.5)], 305.0);
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = loaded(
        &thermal,
        &curves,
        data(&cards, &curves, &card, &point),
        0.2,
        [[0.0; 2]; 2],
        1.0,
    );
    let history = old.checkpoint().clone();
    let cells = thermal.cells().clone();
    let calls = std::cell::Cell::new(0);
    let refused = old.step_candidate(
        &thermal,
        1.0,
        frame(),
        motion(1.0),
        &work("unsupported temperature"),
        controls(),
        &[0.0; 2],
        |state| {
            calls.set(calls.get() + 1);
            other(state)
        },
        || false,
    );
    assert!(
        calls.get() > 0,
        "source must pass at initial state before the late trial refusal"
    );
    assert!(matches!(
        refused,
        Err(CaloricContactError::Material(SlidingError::Material(_)))
    ));
    assert_eq!(old.checkpoint(), &history);
    assert_eq!(thermal.cells(), &cells);
    near(old.time_s(), 0.0);
    near(thermal.time_s(), 0.0);
    let mut twist = motion(0.0);
    twist.torsional_spin_rad_per_s = 1.0;
    assert!(
        old.step_candidate(
            &thermal,
            1.0,
            frame(),
            twist,
            &work("unsupported torsion"),
            controls(),
            &[0.0; 2],
            other,
            || false
        )
        .is_err()
    );
    let data = data(&cards, &curves, &card, &point);
    assert!(
        old.clone()
            .with_material_friction(&thermal, data, 1.0, InputAuthority::SyntheticFixture)
            .is_err(),
        "already-pinned source must not be rebound"
    );
}

#[path = "resolved_material.rs"]
mod resolved;
