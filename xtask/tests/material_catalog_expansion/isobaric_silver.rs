//! Source Cp reaches the original heat timeline without fabricating rho(T).
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, PropertyClaim, PropertyValue, Provenance,
};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink,
    ThermalStorageBasis, TimedEvent,
};

const SLUG: &str = "silver-pure-smith-fickett-nist1995";
const CP_NAME: &str = "specific-heat-capacity-cp";
const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const CONFIG: StepConfig = StepConfig {
    max_dt_s: 0.05,
    energy_tolerance_j: 1e-8,
};

fn point(t: f64, pressure: f64) -> QueryPoint {
    temperature(flags(&[("engineering-reference-model", 1.0)]), "T", t)
        .with_quantity("p", QuantitySpec::dimensional(PA), pressure)
        .unwrap()
}
fn state(name: &str) -> MaterialStateId {
    MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "isobaric single-phase caloric reference".into(),
        revision: 0,
    }
}
fn law(
    card: &MaterialCard,
    lower: &QueryPoint,
    upper: &QueryPoint,
    datum: f64,
) -> Result<IsobaricHeatCapacityLaw, fs_material::thermal::transport::ThermalTransportError> {
    IsobaricHeatCapacityLaw::try_new(
        card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key(CP_NAME, CP),
            temperature_axis: "T",
            pressure_axis: "p",
            lower,
            upper,
            selection: SELECTION,
        },
        datum,
    )
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:e} != {expected:e}"
    );
}
fn analytic_card(scale: f64, pressure_known: bool) -> MaterialCard {
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        point(300.0, 101000.0).axis_quantities()["T"],
        280.0,
        340.0,
    );
    if pressure_known {
        validity = validity.with_quantity("p", QuantitySpec::dimensional(PA), 101000.0, 101000.0);
    }
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key(CP_NAME, CP),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: TEMP,
                knots: [280.0, 300.0, 310.0, 340.0]
                    .map(|t| (t, scale * (2.0 * t - 100.0)))
                    .to_vec(),
                dims: CP,
            },
            validity,
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored analytic Cp fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(state("analytic caloric material"), claims, vec![]).unwrap()
}
fn cell(law: &IsobaricHeatCapacityLaw, id: u64, t: f64, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(id),
        mass_kg: 0.1,
        specific_h_j_kg: law.specific_enthalpy_j_kg(t).unwrap(),
        geometry: CellGeometry {
            start_m: [0.0, 0.0, id as f64],
            end_m: [0.01, 0.0, id as f64],
            volume_m3: 1e-5,
            surface_area_m2: 0.01,
        },
        boundary: Boundary {
            ambient_k: 293.15,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: power,
        },
    }
}
fn timeline<'a>(
    cells: Vec<DepositedCell>,
    laws: &[(MaterialId, &'a dyn ThermalLaw)],
    link: Option<ThermalLink>,
) -> DepositionSimulation<'a> {
    let mut schedule: Vec<_> = cells
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    if let Some(link) = link {
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(link),
        });
    }
    DepositionSimulation::new(
        schedule,
        laws,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 20000,
        },
    )
    .unwrap()
    .advance_to(0.0, CONFIG, || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_g3_isobaric_cp_uses_the_shared_integral_inverse_and_source_conditions() {
    let card = analytic_card(1.0, true);
    let lower = point(280.0, 101000.0);
    let upper = point(340.0, 101000.0);
    for datum in [280.0, 300.0, 340.0] {
        let law = law(&card, &lower, &upper, datum).unwrap();
        assert_eq!(law.material_card_identity(), card.content_hash());
        assert_eq!(law.storage_basis(), ThermalStorageBasis::Enthalpy);
        assert_eq!(law.source_line().knots.len(), 4);
        assert_eq!(law.pressure_pa(), 101000.0);
        assert_eq!(law.query_endpoints(), &[lower.clone(), upper.clone()]);
        assert_eq!(law.maximum_temperature_slope(), (1.0_f64 / 460.0).next_up());
        for (_, answer) in &law.source_line().knots {
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
        for t in [280.0, 297.25, 300.0, 310.0, 337.5, 340.0] {
            // Independent primitive for Cp=2T-100; no production sampler.
            let expected_h = (t - datum) * (t + datum - 100.0);
            let h = law.specific_enthalpy_j_kg(t).unwrap();
            near(h, expected_h, 2e-10);
            near(law.temperature_k(h).unwrap(), t, 1e-12);
            near(
                law.specific_entropy_change_at_temperatures_j_kg_k(datum, t)
                    .unwrap(),
                2.0 * (t - datum) - 100.0 * (t / datum).ln(),
                2e-13,
            );
        }
        let [lo, hi] = law.specific_enthalpy_span_j_kg();
        assert!(law.temperature_k(lo.next_down()).is_err());
        assert!(law.temperature_k(hi.next_up()).is_err());
        assert!(law.temperature_after_specific_heat(340.0, 1.0).is_err());
        assert!(law.temperature_after_specific_heat(280.0, -1.0).is_err());
        assert!(law.specific_heat_cp_j_kg_k(f64::NAN).is_err());
    }
    let a = law(&card, &lower, &upper, 300.0).unwrap();
    let b = law(&card, &lower, &upper, 301.0).unwrap();
    assert_ne!(a.identity(), b.identity());
    assert!(law(&card, &lower, &upper, 279.0).is_err());
    assert!(law(&card, &lower, &point(340.0, 102000.0), 300.0).is_err());
    assert!(
        law(&analytic_card(1.0, false), &lower, &upper, 300.0).is_err(),
        "an unsupported pressure query must not acquire source authority"
    );
    assert!(law(&analytic_card(0.0, true), &lower, &upper, 300.0).is_err());
    let wrong_units = lower
        .clone()
        .with_quantity("p", QuantitySpec::dimensional(Dims::NONE), 101000.0)
        .unwrap();
    assert!(law(&card, &wrong_units, &upper, 300.0).is_err());

    // Material substitution changes only source data, on the same native owner.
    let doubled = analytic_card(2.0, true);
    let dense_law = law(&doubled, &lower, &upper, 300.0).unwrap();
    for (law, scale) in [(&a, 1.0), (&dense_law, 2.0)] {
        let initial = timeline(
            vec![cell(law, 1, 300.0, 1.0)],
            &[(MaterialId(1), law)],
            None,
        );
        let (heated, audit) = initial.advance_to(10.0, CONFIG, || Ok(())).unwrap();
        let delta = 200.0 / scale / (500.0 + (250000.0 + 400.0 / scale).sqrt());
        near(
            heated.temperature_k(CellId(1)).unwrap(),
            300.0 + delta,
            2e-12,
        );
        near(audit.external_heat_j, 10.0, 1e-10);
        near(heated.energy_j() - initial.energy_j(), 10.0, 1e-10);
    }
}

#[test]
fn g1_g3_silver_caloric_columns_preserve_original_values_and_domains() {
    let pack = compile(SLUG, 16);
    let energy = Dims([2, 0, -2, 0, 0, 0]);
    for [t, cp, h, s] in [
        [1.0, 0.00758, 0.00340, 0.00654],
        [10.0, 1.696, 4.2829, 0.5896],
        [100.0, 186.34, 9659.78, 158.71],
        [273.2, 232.87, 47413.2, 374.06],
        [290.0, 234.36, 51338.1, 388.06],
        [295.0, 234.73, 52518.6, 392.05],
        [298.2, 235.01, 53261.8, 394.55],
        [300.0, 235.19, 53695.2, 396.03],
    ] {
        sample(&pack, CP_NAME, CP, &point(t, 101000.0), cp);
        sample(
            &pack,
            "specific-enthalpy-relative-to-0-k",
            energy,
            &point(t, 101000.0),
            h,
        );
        sample(&pack, "specific-entropy", CP, &point(t, 101000.0), s);
    }
    for invalid in [
        point(0.0, 101000.0),
        point(300.01, 101000.0),
        point(295.0, 101325.0),
        point(295.0, 101000.0)
            .with("engineering-reference-model", 0.0)
            .unwrap(),
    ] {
        refuses(&pack, CP_NAME, CP, &invalid);
    }
    // The added thermal observations do not extend the old density domain.
    refuses(&pack, "density", RHO, &point(295.0, 101000.0));
    let normalized =
        NormalizedMaterialCardPack::new(state("Ag; NIST caloric reference"), pack).unwrap();
    let full = law(
        normalized.card(),
        &point(1.0, 101000.0),
        &point(300.0, 101000.0),
        293.15,
    )
    .unwrap();
    assert_eq!(full.source_line().knots.len(), 78);
    for (_, answer) in &full.source_line().knots {
        normalized
            .card()
            .claims()
            .verify_receipt(&answer.receipt)
            .unwrap();
        assert!(answer.receipt.observation_backed);
    }
    // Independently integrate the 290/295 K Cp row. Rounded H is a separate source.
    near(
        full.specific_enthalpy_change_j_kg(290.0, 295.0).unwrap(),
        1172.725,
        1e-10,
    );
    assert!((1172.725 - (52518.6 - 51338.1)).abs() > 1.0);
}

#[test]
fn g1_g3_silver_heating_and_finite_body_cooling_keep_native_energy_and_entropy() {
    let normalized =
        NormalizedMaterialCardPack::new(state("Ag; NIST caloric reference"), compile(SLUG, 16))
            .unwrap();
    let law = law(
        normalized.card(),
        &point(280.0, 101000.0),
        &point(300.0, 101000.0),
        290.0,
    )
    .unwrap();
    // 0.1kg * integral_290^295 Cp dT = 117.2725J, with no rho(T) assumption.
    let initial = timeline(
        vec![cell(&law, 1, 290.0, 11.72725)],
        &[(MaterialId(1), &law)],
        None,
    );
    let (middle, first) = initial.advance_to(5.0, CONFIG, || Ok(())).unwrap();
    let (heated, second) = middle.advance_to(10.0, CONFIG, || Ok(())).unwrap();
    assert!(middle.temperature_k(CellId(1)).unwrap() > 290.0);
    near(heated.temperature_k(CellId(1)).unwrap(), 295.0, 2e-11);
    near(
        first.external_heat_j + second.external_heat_j,
        117.2725,
        1e-9,
    );
    near(heated.energy_j() - initial.energy_j(), 117.2725, 1e-9);
    near(heated.mass_kg(), 0.1, 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert!(heated.advance_to(30.0, CONFIG, || Ok(())).is_err());
    assert!(
        initial
            .advance_to(1.0, CONFIG, || Err(DepositionError::Cancelled))
            .is_err()
    );

    let cold = cell(&law, 1, 290.0, 0.0);
    let hot = cell(&law, 2, 299.0, 0.0);
    let initial = timeline(
        vec![cold, hot],
        &[(MaterialId(1), &law), (MaterialId(2), &law)],
        Some(ThermalLink {
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: 0.001,
            conductance_w_k: 1.0,
        }),
    );
    let (final_state, audit) = initial.advance_to(20.0, CONFIG, || Ok(())).unwrap();
    let cold_end = final_state.temperature_k(CellId(1)).unwrap();
    let hot_end = final_state.temperature_k(CellId(2)).unwrap();
    assert!(cold_end > 290.0 && hot_end < 299.0 && cold_end < hot_end);
    near(final_state.energy_j(), initial.energy_j(), 1e-8);
    near(audit.external_heat_j, 0.0, 0.0);
    let entropy = 0.1
        * law
            .specific_entropy_change_at_temperatures_j_kg_k(290.0, cold_end)
            .unwrap()
        + 0.1
            * law
                .specific_entropy_change_at_temperatures_j_kg_k(299.0, hot_end)
                .unwrap();
    assert!(entropy > 0.0);
    for id in [CellId(1), CellId(2)] {
        let from = initial.cells()[&id].specific_h_j_kg;
        let to = final_state.cells()[&id].specific_h_j_kg;
        let native_entropy = law
            .specific_entropy_change_j_kg_k(from, to)
            .unwrap()
            .unwrap();
        near(
            native_entropy,
            law.specific_entropy_change_at_temperatures_j_kg_k(
                initial.temperature_k(id).unwrap(),
                final_state.temperature_k(id).unwrap(),
            )
            .unwrap(),
            1e-12,
        );
    }
}
