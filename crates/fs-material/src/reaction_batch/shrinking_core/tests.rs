use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_thermochem::stoichiometry::minerals::MineralReaction;
use std::collections::BTreeMap;

fn source(resistances: [f64; 3]) -> ConstitutiveModelCard {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims) in AXES {
        let lo = if [TEMPERATURE_AXIS, VOLUME_AXIS, RADIUS_AXIS].contains(&axis) {
            1e-9
        } else {
            0.0
        };
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), lo, 1e9);
    }
    ConstitutiveModelCard {
        law: LawId(LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            (FILM_RESISTANCE, resistances[0], TransferResistance::DIMS),
            (SURFACE_RESISTANCE, resistances[1], TransferResistance::DIMS),
            (
                DIFFUSION_RESISTANCE,
                resistances[2],
                DiffusionResistance::DIMS,
            ),
            (SOLID_MOLAR_DENSITY, 10_000.0, AmountConcentration::DIMS),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect::<BTreeMap<_, _>>(),
        validity,
        sources: vec![hash_domain(
            "shrinking-core-test",
            b"authored analytic fixture, no ore calibration",
        )],
        provenance: Provenance {
            source: "synthetic spherical fixture; not measured calcite kinetics".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn model(resistances: [f64; 3]) -> ShrinkingCoreReaction {
    ShrinkingCoreReaction::from_card(
        MineralReaction::CalciteBicarbonateNeutralization
            .build()
            .unwrap(),
        SpeciesId::new("CaCO3(s)").unwrap(),
        SpeciesId::new("H+(aq)").unwrap(),
        source(resistances),
    )
    .unwrap()
}
fn inventory(model: &ShrinkingCoreReaction, solid: f64, reagent: f64) -> Vec<(SpeciesId, Amount)> {
    model
        .reaction()
        .terms()
        .iter()
        .map(|term| {
            (
                term.species.clone(),
                Amount::new(if term.species == model.solid {
                    solid
                } else if term.species == model.reagent {
                    reagent
                } else {
                    0.0
                }),
            )
        })
        .collect()
}
fn run(model: &ShrinkingCoreReaction, solid: f64, reagent: f64, t: f64) -> ShrinkingCoreAdvance {
    model
        .advance(
            &inventory(model, solid, reagent),
            Volume::new(1.0),
            Length::new(0.001),
            Time::new(t),
            Temperature::new(298.15),
            &QueryPoint::new(),
            LeachNumerics::default(),
            || false,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn film_limited_finite_reagent_matches_exponential_and_conserves_all_species() {
    let m = model([10.0, 0.0, 0.0]);
    // rho R = 10 mol/m2, inverse k_f=10 s/m. dX/dt=.3(.5-X).
    let t = core::f64::consts::LN_2 / 0.3;
    let out = run(&m, 10.0, 5.0, t);
    close(out.conversion.value(), 0.25, 2e-9);
    close(out.extent.value(), 2.5, 2e-8);
    close(out.reactants[0].remaining.value(), 7.5, 2e-8);
    close(out.reactants[1].remaining.value(), 2.5, 2e-8);
    for product in ["Ca[2+](aq)", "HCO3[-](aq)"] {
        close(
            stock(&out.final_inventory, &SpeciesId::new(product).unwrap()),
            2.5,
            2e-8,
        );
    }
    assert_eq!(out.model().source(), m.source());
    assert_eq!(out.model().reaction().terms(), m.reaction().terms());
    assert!(out.estimated_time_error.value() <= LeachNumerics::default().tolerance(t).unwrap());
    assert!(out.evaluations <= LeachNumerics::default().max_evaluations);
    // Four half-lives cannot cross the finite-reagent ceiling X=.5.
    let later = run(&m, 10.0, 5.0, 4.0 * t);
    close(later.conversion.value(), 0.5 * 15.0 / 16.0, 2e-9);
    assert!(later.reactants[1].remaining.value() > 0.0);
}

#[test]
fn surface_diffusion_and_combined_matched_capacity_known_answers() {
    let ln2 = core::f64::consts::LN_2;
    // n_A0=n_B0 => C_A=C_A0*y^3. At y=r/R=.5 the three integrals
    // are ln2, (1/y^2-1)/2=1.5, and 1/y+ln(y)-1=1-ln2.
    for (r, t) in [
        ([10.0, 0.0, 0.0], 10.0 * ln2),
        ([0.0, 20.0, 0.0], 30.0),
        ([0.0, 0.0, 30_000.0], 30.0 * (1.0 - ln2)),
        ([10.0, 20.0, 30_000.0], 60.0 - 20.0 * ln2),
    ] {
        let out = run(&model(r), 10.0, 10.0, t);
        close(out.conversion.value(), 0.875, 4e-9);
        close(out.core_radius.value(), 0.0005, 3e-12);
        close(out.reactants[1].remaining.value(), 1.25, 4e-8);
    }
}

#[test]
fn near_stoichiometric_excess_resolves_short_time_without_integrating_completion() {
    // rho_B R r_s = 1, C_A0 = 1 and demand = .9999 mol/m3. The
    // complete-dissolution integral has a sharp late-time tail, but the
    // requested early-time state must remain solvable on a modest budget.
    let m = model([0.0, 0.1, 0.0]);
    let target = 0.01;
    let settings = LeachNumerics {
        max_evaluations: 4096,
        ..LeachNumerics::default()
    };
    let out = m
        .advance(
            &inventory(&m, 0.9999, 1.0),
            Volume::new(1.0),
            Length::new(0.001),
            Time::new(target),
            Temperature::new(298.15),
            &QueryPoint::new(),
            settings,
            || false,
        )
        .unwrap();
    // Independent analytical bounds: exact stoichiometric loading gives
    // y=(1+2t)^(-1/2); replenishing the bath gives u=t. Slight excess
    // reagent must advance farther than the first, less than the second.
    let matched_radius = 1.0 / (1.0 + 2.0 * target).sqrt();
    let matched_conversion = 1.0 - matched_radius.powi(3);
    let maintained_conversion = 1.0 - (1.0 - target).powi(3);
    assert!(out.conversion.value() >= matched_conversion - 5e-9);
    assert!(out.conversion.value() < maintained_conversion);
    assert!(out.reactants[0].remaining.value() > 0.0);
    assert!(out.reactants[1].remaining.value() > 0.0);
    close(out.extent.value(), 0.9999 * out.conversion.value(), 1e-15);
    assert!(out.estimated_time_error.value() <= settings.tolerance(target).unwrap());
    assert!(out.evaluations <= settings.max_evaluations);
}

#[test]
fn maintained_concentration_has_correct_radius_scaling_and_series_addition() {
    let time = |m: &ShrinkingCoreReaction, r| {
        m.constant_concentration_time(
            Dimensionless::new(0.5),
            Length::new(r),
            AmountConcentration::new(10.0),
            Volume::new(1.0),
            AmountConcentration::new(10.0),
            Temperature::new(298.15),
            &QueryPoint::new(),
        )
        .unwrap()
        .value()
    };
    let f = model([10.0, 0.0, 0.0]);
    let s = model([0.0, 20.0, 0.0]);
    let d = model([0.0, 0.0, 30_000.0]);
    close(time(&f, 0.002), 2.0 * time(&f, 0.001), 1e-13);
    close(time(&s, 0.002), 2.0 * time(&s, 0.001), 1e-13);
    close(time(&d, 0.002), 4.0 * time(&d, 0.001), 1e-13);
    let combined = model([10.0, 20.0, 30_000.0]);
    close(
        time(&combined, 0.001),
        time(&f, 0.001) + time(&s, 0.001) + time(&d, 0.001),
        1e-13,
    );
    let closed = run(&combined, 10.0, 10.0, time(&combined, 0.001));
    assert!(
        closed.conversion.value() < 0.875,
        "depletion must slow leaching before any ceiling is reached"
    );
}

#[test]
fn excess_reagent_allows_finite_completion_and_zero_controls_preserve_inventory() {
    let m = model([10.0, 20.0, 30_000.0]);
    let complete = run(&m, 10.0, 100.0, 100.0);
    assert_eq!(complete.conversion.value(), 1.0);
    assert_eq!(complete.core_radius.value(), 0.0);
    assert_eq!(complete.reactants[0].remaining.value(), 0.0);
    assert_eq!(complete.reactants[1].remaining.value(), 90.0);
    for (solid, reagent, time) in [(10.0, 10.0, 0.0), (0.0, 10.0, 10.0), (10.0, 0.0, 10.0)] {
        let out = run(&m, solid, reagent, time);
        assert_eq!(out.extent.value(), 0.0);
        assert_eq!(out.initial_inventory, out.final_inventory);
        assert_eq!(out.evaluations, 0);
    }
}

#[test]
fn stoichiometric_scaling_and_batch_scale_preserve_physical_conversion() {
    let m = model([10.0, 20.0, 30_000.0]);
    let a = run(&m, 10.0, 10.0, 10.0);
    let mut stock = inventory(&m, 100.0, 100.0);
    stock.reverse();
    let b = m
        .advance(
            &stock,
            Volume::new(10.0),
            Length::new(0.001),
            Time::new(10.0),
            Temperature::new(298.15),
            &QueryPoint::new(),
            LeachNumerics::default(),
            || false,
        )
        .unwrap();
    assert_eq!(a.conversion, b.conversion);
    close(b.extent.value(), 10.0 * a.extent.value(), 1e-13);
    // Two acid moles per calcite change both particle flux and bath demand.
    let acid = ShrinkingCoreReaction::from_card(
        MineralReaction::CalciteAcidNeutralization.build().unwrap(),
        SpeciesId::new("CaCO3(s)").unwrap(),
        SpeciesId::new("H+(aq)").unwrap(),
        source([10.0, 20.0, 30_000.0]),
    )
    .unwrap();
    let out = run(&acid, 10.0, 20.0, 10.0);
    assert_eq!(out.conversion, a.conversion);
    close(
        out.reactants[1].consumed.value(),
        2.0 * out.extent.value(),
        1e-13,
    );
}

#[test]
fn source_support_and_model_parameters_are_not_inferred_or_overridden() {
    let m = model([10.0, 20.0, 30_000.0]);
    let advance = |model: &ShrinkingCoreReaction, conditions: &QueryPoint| {
        model.advance(
            &inventory(model, 10.0, 10.0),
            Volume::new(1.0),
            Length::new(0.001),
            Time::new(10.0),
            Temperature::new(298.15),
            conditions,
            LeachNumerics::default(),
            || false,
        )
    };
    assert!(matches!(
        advance(
            &m,
            &QueryPoint::new().with(CONCENTRATION_AXIS, 1e6).unwrap()
        ),
        Err(ShrinkingCoreError::Invalid(_))
    ));
    let mut limited = m.clone();
    limited.source.validity = limited.source.validity.clone().with_quantity(
        CONCENTRATION_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        9.99,
        10.0,
    );
    assert!(matches!(
        advance(&limited, &QueryPoint::new()),
        Err(ShrinkingCoreError::OutsideValidity)
    ));
    let mut conditional = m.clone();
    conditional.source.validity = conditional.source.validity.clone().with_quantity(
        "stirrer-speed",
        QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
        1.0,
        2.0,
    );
    assert!(matches!(
        advance(&conditional, &QueryPoint::new()),
        Err(ShrinkingCoreError::OutsideValidity)
    ));
    assert!(
        advance(
            &conditional,
            &QueryPoint::new()
                .with_quantity(
                    "stirrer-speed",
                    QuantitySpec::dimensional(fs_qty::Frequency::DIMS),
                    1.5
                )
                .unwrap()
        )
        .is_ok()
    );
    for mutate in [0, 1, 2, 3] {
        let mut card = source([10.0, 20.0, 30_000.0]);
        match mutate {
            0 => {
                card.parameters.get_mut(DIFFUSION_RESISTANCE).unwrap().dims =
                    TransferResistance::DIMS
            }
            1 => card.parameters.get_mut(SOLID_MOLAR_DENSITY).unwrap().value = 0.0,
            2 => card.sources.clear(),
            _ => card.parameters.get_mut(FILM_RESISTANCE).unwrap().value = -1.0,
        }
        assert!(
            ShrinkingCoreReaction::from_card(
                m.reaction.clone(),
                m.solid.clone(),
                m.reagent.clone(),
                card
            )
            .is_err()
        );
    }
}

#[test]
fn work_exhaustion_cancellation_and_inventory_failures_are_atomic() {
    let m = model([10.0, 20.0, 30_000.0]);
    let input = inventory(&m, 10.0, 10.0);
    let original = input.clone();
    let advance = |n, settings, cancel: &mut dyn FnMut() -> bool| {
        m.advance(
            n,
            Volume::new(1.0),
            Length::new(0.001),
            Time::new(10.0),
            Temperature::new(298.15),
            &QueryPoint::new(),
            settings,
            cancel,
        )
    };
    assert!(matches!(
        advance(
            &input,
            LeachNumerics {
                max_evaluations: 5,
                ..LeachNumerics::default()
            },
            &mut || false
        ),
        Err(ShrinkingCoreError::Budget)
    ));
    let mut calls = 0;
    assert_eq!(
        advance(&input, LeachNumerics::default(), &mut || {
            calls += 1;
            calls > 10
        })
        .unwrap_err(),
        ShrinkingCoreError::Cancelled
    );
    let mut malformed = input.clone();
    malformed.push(input[0].clone());
    assert!(matches!(
        advance(&malformed, LeachNumerics::default(), &mut || false),
        Err(ShrinkingCoreError::Reaction(_))
    ));
    assert_eq!(input, original);
    assert_eq!(
        advance(&input, LeachNumerics::default(), &mut || false)
            .unwrap()
            .conversion,
        run(&m, 10.0, 10.0, 10.0).conversion
    );
}
