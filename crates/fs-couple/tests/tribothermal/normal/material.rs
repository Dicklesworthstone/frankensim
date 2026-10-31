//! G1 finite-temperature references with original pinned elastic/caloric data.
//! These authored curves are implementation fixtures, not measured calibration.
use super::*;
use fs_couple::tribothermal::normal::material::MaterialNormalModulus;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::elastic::thermal::contact::{ThermalHertzError, ThermalHertzPair};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY,
    SPECIFIC_HEAT_CAPACITY_DIMS, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
    resolve_isotropic_solid_state_point,
};
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const CP: [f64; 2] = [0.1, 0.2];
const INITIAL: [f64; 2] = [300.0, 360.0];
fn point(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            pressure,
        )
        .unwrap()
}
fn card(index: usize, density_segments: &[(f64, f64, f64)]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims: Dims, low: f64, high: f64, values: [f64; 2]| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(low, values[0]), (high, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), low, high)
                    .with_quantity("p", point(300.0, 1e5).axis_quantities()["p"], 1e5, 2e5),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored coupled normal contact".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for &(lo, hi, density) in density_segments {
        add(DENSITY_PROPERTY, Density::DIMS, lo, hi, [density; 2]);
    }
    let young = if index == 0 { [2.2e6, 1.4e6] } else { [3e6; 2] };
    add(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 200.0, 600.0, young);
    add(
        POISSON_RATIO_PROPERTY,
        Dims::NONE,
        200.0,
        600.0,
        [if index == 0 { 0.25 } else { 0.3 }; 2],
    );
    add(
        YIELD_STRESS_PROPERTY,
        Pressure::DIMS,
        200.0,
        600.0,
        [1e12; 2],
    );
    add(
        "specific_heat_capacity",
        SPECIFIC_HEAT_CAPACITY_DIMS,
        200.0,
        600.0,
        [CP[index]; 2],
    );
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("normal-reference-{index}"),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard, pressure: f64, range: [f64; 2]) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(range[0], pressure),
            upper: &point(range[1], pressure),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap()
}
fn data() -> ([MaterialCard; 2], [ResolvedHeatCapacityCurve; 2]) {
    let cards = [
        card(0, &[(200.0, 600.0, 1.0)]),
        card(1, &[(200.0, 600.0, 1.0)]),
    ];
    let curves = [
        caloric(&cards[0], 1e5, [200.0, 600.0]),
        caloric(&cards[1], 1e5, [200.0, 600.0]),
    ];
    (cards, curves)
}
fn source(cards: &[MaterialCard; 2], ranges: [[f64; 2]; 2]) -> ThermalHertzPair {
    let state = |i| {
        resolve_isotropic_solid_state_point(
            &cards[i],
            &point(INITIAL[i], 1e5),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    };
    ThermalHertzPair::new([&cards[0], &cards[1]], [&state(0), &state(1)], ranges).unwrap()
}
fn normal_request(source: &ThermalHertzPair, indentation: f64) -> NormalPatchRequest {
    let mut request = request(indentation, 0.0);
    request.law = NormalPatchLaw::HertzSpherePlane {
        effective_radius_m: 1.0,
        reduced_modulus_pa: source
            .at_temperatures(source.reference_temperatures_k())
            .unwrap()
            .reduced_modulus
            .value_pa,
    };
    request.applicability.yield_strength_pa = source.limiting_yield_stress_pa();
    request
}
fn initial(curves: &[ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'_> {
    let events = CELLS
        .iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: *id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: CP[i] * (INITIAL[i] - 200.0),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: common::boundary(),
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 4,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn bind_source<'a>(
    thermal: &DepositionSimulation<'_>,
    cards: &[MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    indentation: f64,
    ranges: [[f64; 2]; 2],
) -> CaloricNormalContact<'a, MaterialNormalModulus<'a>> {
    let source = source(cards, ranges);
    let model = MaterialNormalModulus::new(
        normal_request(&source, indentation),
        source,
        [&curves[0], &curves[1]],
    )
    .unwrap();
    CaloricNormalContact::new(
        thermal,
        model,
        [curves[0].curve(), curves[1].curve()],
        CELLS,
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
    )
    .unwrap()
}
fn bisect(mut lo: f64, mut hi: f64, fun: impl Fn(f64) -> f64) -> f64 {
    assert!(fun(lo) < 0.0 && fun(hi) > 0.0);
    for _ in 0..70 {
        let middle = 0.5 * (lo + hi);
        if fun(middle) > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    0.5 * (lo + hi)
}
// Independent rational formula and its elementary derivative, not calls to
// either production material or contact owner.
fn reference(t: f64) -> (f64, f64) {
    let ea = 2e6 - 2000.0 * (t - 300.0);
    let eb = 3e6;
    let a = 1.0 - 0.25 * 0.25;
    let b = 1.0 - 0.3 * 0.3;
    let denominator = a * eb + b * ea;
    (
        ea * eb / denominator,
        -2000.0 * a * eb * eb / (denominator * denominator),
    )
}
fn residual(t: f64, d0: f64, d1: f64, heat: f64) -> f64 {
    let (e, g) = reference(t);
    let (e0, g0) = reference(INITIAL[0]);
    CP[0] * (t - INITIAL[0]) + basis(d0) * (e - e0) - t * basis(d1) * g
        + INITIAL[0] * basis(d0) * g0
        - heat
}

#[test]
fn source_hertz_loading_matches_independent_scalar_caloric_root() {
    let (cards, curves) = data();
    let thermal = initial(&curves);
    let contact = bind_source(&thermal, &cards, &curves, 0.002, [[200.0, 600.0]; 2]);
    let candidate = contact
        .step_candidate(
            &thermal,
            1.0,
            0.01,
            controls(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let t = bisect(200.0, 300.0, |t| residual(t, 0.002, 0.01, 0.0));
    let (e, g) = reference(t);
    near(candidate.contact().temperature_k()[0], t);
    near(candidate.contact().temperature_k()[1], 360.0);
    // The solved energy tolerance propagates through g=dE*/dT; compare Pa
    // relatively instead of applying the small-energy absolute tolerance.
    assert!((candidate.report().modulus.modulus_pa - e).abs() < 1e-11 * e);
    near(
        candidate.report().endpoint.normal_force_n,
        (4.0 / 3.0) * e * 0.01_f64.powf(1.5),
    );
    near(
        candidate.report().thermodynamics.internal_energy_j,
        basis(0.01) * (e - t * g),
    );
    near(candidate.report().reversible_heat_j[0], CP[0] * (t - 300.0));
    near(candidate.report().reversible_heat_j[1], 0.0);
    near(
        candidate.report().work.work_into_contact_j,
        e * (basis(0.01) - basis(0.002)),
    );
    assert_eq!(candidate.report().work.physical_dissipation_j, 0.0);
    near(candidate.report().total_energy_residual_j, 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(contact.temperature_k(), INITIAL);
}

#[test]
fn source_hertz_held_heating_includes_nonlinear_internal_energy() {
    let (cards, curves) = data();
    let thermal = initial(&curves);
    let contact = bind_source(&thermal, &cards, &curves, 0.01, [[200.0, 600.0]; 2]);
    let heat = 3.0;
    let candidate = contact
        .step_candidate(
            &thermal,
            1.0,
            0.01,
            controls(),
            &guess(&thermal),
            |_| {
                Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::from([(CELLS[0], heat)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || false,
        )
        .unwrap();
    let t = bisect(300.0, 400.0, |t| residual(t, 0.01, 0.01, heat));
    let (e, g) = reference(t);
    let u = basis(0.01) * (e - t * g);
    near(candidate.contact().temperature_k()[0], t);
    near(candidate.report().thermodynamics.internal_energy_j, u);
    near(
        candidate.report().reversible_heat_j[0],
        -(u - contact.thermodynamics().internal_energy_j),
    );
    near(candidate.report().work.work_into_contact_j, 0.0);
    near(candidate.report().external_heat_j, heat);
    assert_eq!(candidate.report().work.physical_dissipation_j, 0.0);
    assert!(candidate.report().work.internal_energy_change_j > 0.0);
    assert!(candidate.report().work.coefficient_free_energy_change_j < 0.0);
    assert!(
        (candidate.report().work.internal_energy_change_j
            - candidate.report().work.coefficient_free_energy_change_j)
            .abs()
            > 0.1
    );
    near(
        CP[0] * (t - 300.0) + u - contact.thermodynamics().internal_energy_j,
        heat,
    );
    near(candidate.report().total_energy_residual_j, 0.0);
}

#[test]
fn source_caloric_mismatch_and_late_domain_failure_preserve_normal_state() {
    let (cards, curves) = data();
    let thermal = initial(&curves);
    let pair = source(&cards, [[200.0, 600.0]; 2]);
    let request = normal_request(&pair, 0.002);
    assert!(
        MaterialNormalModulus::new(request.clone(), pair.clone(), [&curves[1], &curves[0]])
            .is_err()
    );
    let foreign_context = caloric(&cards[0], 2e5, [200.0, 600.0]);
    assert!(MaterialNormalModulus::new(request, pair, [&foreign_context, &curves[1]]).is_err());
    // Both independently valid sources use the same card, but density changes
    // claim between disjoint source domains. The wrapper must retain its pin.
    let split = [
        card(0, &[(200.0, 400.0, 1.0), (450.0, 600.0, 2.0)]),
        cards[1].clone(),
    ];
    let different_density = caloric(&split[0], 1e5, [450.0, 600.0]);
    let split_source = source(&split, [[200.0, 400.0], [200.0, 600.0]]);
    let mismatch = MaterialNormalModulus::new(
        normal_request(&split_source, 0.002),
        split_source,
        [&different_density, &curves[1]],
    );
    assert!(matches!(
        mismatch,
        Err(NormalThermalError::Invalid(
            "caloric and elastic sources must share complete external conditions and original density claim"
        ))
    ));
    let contact = bind_source(
        &thermal,
        &cards,
        &curves,
        0.002,
        [[299.5, 600.0], [200.0, 600.0]],
    );
    let failed = contact.step_candidate(
        &thermal,
        1.0,
        0.01,
        controls(),
        &guess(&thermal),
        other,
        || false,
    );
    assert!(
        matches!(failed,Err(NormalThermalError::MaterialModel(error))
        if matches!(error.downcast_ref::<ThermalHertzError>(),Some(ThermalHertzError::Invalid(
            "temperature leaves the pinned solid source interval"))))
    );
    let run = || {
        contact
            .step_candidate(
                &thermal,
                1.0,
                0.0021,
                controls(),
                &guess(&thermal),
                other,
                || false,
            )
            .unwrap()
    };
    let first = run();
    let replay = run();
    assert_eq!(first.report().work, replay.report().work);
    assert_eq!(
        first.contact().temperature_k(),
        replay.contact().temperature_k()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(contact.temperature_k(), INITIAL);
}

// Independent K/G sources with deliberately noncanonical density names.
// No Young/Poisson source exists on these cards.
fn independent_pair_card(index: usize) -> MaterialCard {
    let base = card(index, &[(200.0, 600.0, 1.0)]);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        if ![
            DENSITY_PROPERTY,
            YOUNG_MODULUS_PROPERTY,
            POISSON_RATIO_PROPERTY,
        ]
        .contains(&original.key.name())
        {
            claims.insert_claim(original.clone()).unwrap();
        }
    }
    let (bulk, shear) = if index == 0 {
        ([1.6e6, 1.2e6], [0.85e6, 0.65e6])
    } else {
        ([2e6; 2], [1.2e6; 2])
    };
    for (name, dims, values) in [
        ("rho-source", Density::DIMS, [1.0; 2]),
        ("other-rho-source", Density::DIMS, [1.0; 2]),
        ("contact-bulk", Pressure::DIMS, bulk),
        ("contact-shear", Pressure::DIMS, shear),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, values[0]), (600.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        200.0,
                        600.0,
                    )
                    .with_quantity("p", point(300.0, 1e5).axis_quantities()["p"], 1e5, 2e5),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored independent contact parameters".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn independent_pair_caloric(card: &MaterialCard, density: &str) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new(density, Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0, 1e5),
            upper: &point(600.0, 1e5),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap()
}

fn independent_pair_source(cards: &[MaterialCard; 2]) -> ThermalHertzPair {
    use fs_material::state_point::{
        IsotropicElasticParameters, IsotropicSolidStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    };
    let requirement = |name, dims| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    };
    let density = requirement("rho-source", Density::DIMS);
    let yield_requirement = requirement(YIELD_STRESS_PROPERTY, Pressure::DIMS);
    let parameters = IsotropicElasticParameters::BulkShear {
        bulk_modulus: requirement("contact-bulk", Pressure::DIMS),
        shear_modulus: requirement("contact-shear", Pressure::DIMS),
    };
    let mut requirements = vec![
        density.clone(),
        yield_requirement.clone(),
        requirement("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
    ];
    requirements.extend(parameters.requirements().into_iter().cloned());
    let states: Vec<_> = cards
        .iter()
        .enumerate()
        .map(|(i, card)| {
            let resolved = resolve_material_state_point(
                card,
                &point(INITIAL[i], 1e5),
                &requirements,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            IsotropicSolidStatePoint::try_from_resolved_with_parameters(
                resolved,
                &density,
                &parameters,
                &yield_requirement,
            )
            .unwrap()
        })
        .collect();
    ThermalHertzPair::new(
        [&cards[0], &cards[1]],
        [&states[0], &states[1]],
        [[200.0, 600.0]; 2],
    )
    .unwrap()
}

// Independent reference through E(K,G), nu(K,G) and their quotient derivatives,
// rather than the production compliance decomposition or derivative owner.
fn independent_pair_reference(t: f64) -> (f64, f64) {
    let k = 1.5e6 - 1000.0 * (t - 300.0);
    let g = 0.8e6 - 500.0 * (t - 300.0);
    let q = 3.0 * k + g;
    let q_prime = -3500.0;
    let e = 9.0 * k * g / q;
    let e_prime = 9.0 * ((-1000.0 * g - 500.0 * k) * q - k * g * q_prime) / (q * q);
    let n = 3.0 * k - 2.0 * g;
    let nu = n / (2.0 * q);
    let nu_prime = (-2000.0 * q - n * q_prime) / (2.0 * q * q);
    let compliance = (1.0 - nu * nu) / e + (1.0 - 0.25 * 0.25) / 3e6;
    let compliance_prime = -2.0 * nu * nu_prime / e - (1.0 - nu * nu) * e_prime / (e * e);
    let modulus = 1.0 / compliance;
    (modulus, -modulus * modulus * compliance_prime)
}

#[test]
fn independent_bulk_shear_heating_uses_original_caloric_owner_and_conserves_energy() {
    let cards = [independent_pair_card(0), independent_pair_card(1)];
    let curves = [
        independent_pair_caloric(&cards[0], "rho-source"),
        independent_pair_caloric(&cards[1], "rho-source"),
    ];
    let thermal = initial(&curves);
    let source = independent_pair_source(&cards);
    let model = MaterialNormalModulus::new(
        normal_request(&source, 0.01),
        source,
        [&curves[0], &curves[1]],
    )
    .unwrap();
    let contact = CaloricNormalContact::new(
        &thermal,
        model,
        [curves[0].curve(), curves[1].curve()],
        CELLS,
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
    )
    .unwrap();
    let heat = 3.0;
    let candidate = contact
        .step_candidate(
            &thermal,
            1.0,
            0.01,
            controls(),
            &guess(&thermal),
            |_| {
                Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::from([(CELLS[0], heat)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
            || false,
        )
        .unwrap();
    let (m0, d0) = independent_pair_reference(300.0);
    let initial_u = basis(0.01) * (m0 - 300.0 * d0);
    let expected_t = bisect(300.0, 400.0, |t| {
        let (m, d) = independent_pair_reference(t);
        CP[0] * (t - 300.0) + basis(0.01) * (m - t * d) - initial_u - heat
    });
    let (modulus, derivative) = independent_pair_reference(expected_t);
    let expected_u = basis(0.01) * (modulus - expected_t * derivative);
    near(candidate.contact().temperature_k()[0], expected_t);
    near(candidate.contact().temperature_k()[1], INITIAL[1]);
    near(
        candidate.report().endpoint.normal_force_n,
        (4.0 / 3.0) * modulus * 0.01_f64.powf(1.5),
    );
    near(
        candidate.report().thermodynamics.internal_energy_j,
        expected_u,
    );
    near(CP[0] * (expected_t - 300.0) + expected_u - initial_u, heat);
    near(candidate.report().work.work_into_contact_j, 0.0);
    near(candidate.report().total_energy_residual_j, 0.0);
    assert_eq!(candidate.report().work.physical_dissipation_j, 0.0);
    let refreshed = candidate
        .contact()
        .model()
        .source()
        .at_temperatures(candidate.contact().temperature_k())
        .unwrap();
    for solid in &refreshed.surfaces {
        assert!(solid.resolved().property(YOUNG_MODULUS_PROPERTY).is_none());
        assert!(solid.resolved().property(POISSON_RATIO_PROPERTY).is_none());
        assert!(
            solid
                .resolved()
                .property("specific_heat_capacity")
                .is_some()
        );
        assert_eq!(
            solid
                .elastic_state()
                .derivation()
                .unwrap()
                .density_requirement()
                .name(),
            "rho-source"
        );
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(contact.temperature_k(), INITIAL);
}

#[test]
fn independent_pair_caloric_binding_rejects_equal_density_from_a_different_source_role() {
    let cards = [independent_pair_card(0), independent_pair_card(1)];
    let source = independent_pair_source(&cards);
    let right = independent_pair_caloric(&cards[1], "rho-source");
    let wrong = independent_pair_caloric(&cards[0], "other-rho-source");
    assert!(matches!(
        MaterialNormalModulus::new(normal_request(&source, 0.01), source, [&wrong, &right]),
        Err(NormalThermalError::Invalid(
            "caloric and elastic sources must share complete external conditions and original density claim"
        ))
    ));
}
