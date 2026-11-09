use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    SelectionPolicy, UncertaintyModel,
};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositionBudget, DepositionSimulation, Event, MaterialId,
    StepConfig, TimedEvent,
    advection::{Advected, MassBoundary},
};

fn key() -> PropertyKey { PropertyKey::new("authored_cv", SPECIFIC_HEAT) }

fn card(density_supported: bool, capacities: [f64; 3]) -> MaterialCard {
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        "temperature", QuantitySpec::dimensional(Temperature::DIMS), 250.0, 500.0,
    );
    if density_supported {
        validity = validity.with_quantity("rho", QuantitySpec::dimensional(DENSITY), 2.0, 4.0);
    }
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: key(),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots: [250.0, 350.0, 500.0].into_iter().zip(capacities).collect(),
            dims: SPECIFIC_HEAT,
        },
        validity,
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic Cv calorimetry oracle".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId {
        chemistry: "numerical single phase".into(),
        phase: "fluid".into(),
        process: "test".into(),
        revision: 0,
    }, claims, vec![]).unwrap()
}

fn point(t: f64, rho: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
        .with_quantity("rho", QuantitySpec::dimensional(DENSITY), rho).unwrap()
}

fn prepare(
    card: &MaterialCard, lower: &QueryPoint, upper: &QueryPoint, datum: f64,
) -> Result<IsochoricHeatCapacityLaw, ThermalTransportError> {
    IsochoricHeatCapacityLaw::try_new(card, IsochoricHeatCapacityQuery {
        specific_heat_cv: &key(),
        temperature_axis: "temperature",
        density_axis: "rho",
        lower,
        upper,
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }, datum)
}

fn law(capacities: [f64; 3]) -> IsochoricHeatCapacityLaw {
    prepare(&card(true, capacities), &point(250.0, 2.0), &point(500.0, 2.0), 300.0).unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-12 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn g1_affine_cv_uses_exact_internal_energy_entropy_and_inverse_across_source_knots() {
    // Independent closed forms for Cv(T)=2T: delta u=T1²-T0²; delta s=2(T1-T0).
    let law = law([500.0, 700.0, 1000.0]);
    assert_eq!(law.storage_basis(), ThermalStorageBasis::InternalEnergy);
    for (a, b) in [(250.0, 500.0), (300.0, 425.0), (490.0, 275.0), (350.0, 350.0)] {
        let energy = (b - a) * (a + b);
        close(law.specific_internal_energy_change_j_kg(a, b).unwrap(), energy);
        close(law.specific_entropy_change_at_temperatures_j_kg_k(a, b).unwrap(), 2.0 * (b - a));
        close(law.temperature_after_specific_heat(a, energy).unwrap(), b);
        let ua = law.specific_internal_energy_j_kg(a).unwrap();
        let ub = law.specific_internal_energy_j_kg(b).unwrap();
        close(law.temperature_k(ub).unwrap(), b);
        close(law.specific_entropy_change_j_kg_k(ua, ub).unwrap().unwrap(), 2.0 * (b - a));
    }
    let t = 400.0_f64;
    let delta = t.next_up() - t;
    close(law.specific_internal_energy_change_j_kg(t, t.next_up()).unwrap() / delta, 2.0 * t);
    close(law.specific_heat_cv_j_kg_k(t).unwrap(), 2.0 * t);
    let bounds = law.specific_internal_energy_span_j_kg();
    assert_eq!(law.temperature_k(bounds[0]).unwrap(), 250.0);
    assert_eq!(law.temperature_k(bounds[1]).unwrap(), 500.0);
    assert!(law.temperature_k(bounds[0].next_down()).is_err());
    assert!(law.temperature_k(bounds[1].next_up()).is_err());
}

#[test]
fn g1_indexed_cv_heat_preserves_relative_changes_and_original_support() {
    let card = card(true, [500.0, 2000.0, 300.0]);
    for datum in [250.0_f64, 301.25, 350.0, 425.25, 500.0] {
        let indexed = prepare(&card, &point(250.0, 2.0), &point(500.0, 2.0), datum).unwrap();
        let mut direct = indexed.clone();
        direct.inverse = None;
        for target in [250.0, 290.25, 350.0, 450.75, 500.0,
            datum.next_down(), datum, datum.next_up()]
        {
            if !(250.0..=500.0).contains(&target) { continue; }
            let expected = direct.specific_internal_energy_j_kg(target).unwrap();
            let u = indexed.specific_internal_energy_j_kg(target).unwrap();
            assert!((u - expected).abs() <= 8.0 * f64::EPSILON * expected.abs());
            assert_eq!(indexed.specific_internal_energy_change_j_kg(datum, target).unwrap(), u);
            assert_eq!(indexed.specific_internal_energy_change_j_kg(target, datum).unwrap(), -u);
            assert!((indexed.temperature_k(u).unwrap() - target).abs()
                <= 8.0 * f64::EPSILON * target);
            if target == 250.0 || target == 500.0 || (target - datum).abs() < 1e-10 {
                assert_eq!(u, expected);
            }
        }
        let a = 400.0_f64;
        let du = indexed.specific_internal_energy_change_j_kg(a, a.next_up()).unwrap();
        assert!(du > 0.0);
        assert_eq!(du, direct.specific_internal_energy_change_j_kg(a, a.next_up()).unwrap());
        for outside in [249.0, 501.0, f64::NAN, f64::INFINITY] {
            assert!(indexed.specific_internal_energy_j_kg(outside).is_err());
            assert!(indexed.specific_internal_energy_change_j_kg(outside, datum).is_err());
        }
    }
}

#[test]
fn g0_source_density_context_is_required_and_remains_in_the_law_identity() {
    let card = card(true, [700.0; 3]);
    let lower = point(250.0, 2.0);
    let upper = point(500.0, 2.0);
    let a = prepare(&card, &lower, &upper, 300.0).unwrap();
    assert_eq!(a.material_card_identity(), card.content_hash());
    assert_eq!(a.query_endpoints(), &[lower.clone(), upper.clone()]);
    assert_eq!(a.caloric_identity(), Some(*a.identity().as_bytes()));
    assert_eq!(a.source_line().knots.len(), 3);
    let b = prepare(&card, &point(250.0, 3.0), &point(500.0, 3.0), 300.0).unwrap();
    let datum = prepare(&card, &lower, &upper, 350.0).unwrap();
    assert_ne!(a.identity(), b.identity());
    assert_ne!(a.identity(), datum.identity());
    close(a.specific_internal_energy_change_j_kg(280.0, 450.0).unwrap(),
        datum.specific_internal_energy_change_j_kg(280.0, 450.0).unwrap());
    assert!(prepare(&card, &lower, &point(500.0, 3.0), 300.0).is_err());
    assert!(prepare(&card, &point(250.0, 1.0), &point(500.0, 1.0), 300.0).is_err());
    assert!(prepare(&card, &lower, &upper, 200.0).is_err());
    let with_pressure = |p: QueryPoint| p.with_quantity(
        "p", QuantitySpec::dimensional(Pressure::DIMS), 1e5,
    ).unwrap();
    assert!(prepare(&card, &with_pressure(lower), &with_pressure(upper), 300.0).is_err());
}

#[test]
fn g0_unqualified_density_and_nonpositive_cv_cannot_supply_a_thermal_law() {
    // Adding density to a query cannot add density support to the source.
    for card in [
        card(false, [700.0; 3]),
        card(true, [700.0, 0.0, 700.0]),
        card(true, [700.0, -1.0, 700.0]),
    ] {
        assert!(prepare(&card, &point(250.0, 2.0), &point(500.0, 2.0), 300.0).is_err());
    }
}

fn cell(id: u64, mass: f64, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: mass,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, 0.0],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: mass / 2.0,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: power,
        },
    }
}

fn simulation(law: &IsochoricHeatCapacityLaw, cells: Vec<DepositedCell>) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        cells.into_iter().map(|cell| TimedEvent {
            time_s: 0.0, event: Event::Birth(cell),
        }).collect(),
        &[(MaterialId(1), law)],
        DepositionBudget { max_cells: 4, max_links: 4, max_events: 8, max_substeps: 1024 },
    ).unwrap()
}

fn config() -> StepConfig {
    StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 }
}

#[test]
fn g1_native_fixed_volume_heat_uses_cv_and_conserves_internal_energy() {
    let law = law([500.0, 700.0, 1000.0]);
    let initial = simulation(&law, vec![cell(0, 2.0, 20_000.0)]);
    let (heated, energy, entropy) = initial
        .advance_to_with_entropy(7.0, config(), || Ok(())).unwrap();
    let value = heated.cells()[&CellId(0)].specific_h_j_kg;
    // Q=140000 J, m=2 kg, Cv=2T gives T=sqrt(300²+Q/m)=400 K.
    close(law.temperature_k(value).unwrap(), 400.0);
    close(energy.external_heat_j, 140_000.0);
    close(energy.residual_j, 0.0);
    close(entropy.thermal_storage_change_j_k, 400.0);
    assert_eq!(heated.material_storage_basis(MaterialId(1)).unwrap(),
        ThermalStorageBasis::InternalEnergy);
    assert_eq!(initial.time_s(), 0.0);
    assert!(initial.cells().is_empty());
    // A frozen initial Cv would give 416.67 K, violating this exact calorimetry.
    assert!((law.temperature_k(value).unwrap() - (300.0 + 70_000.0 / 600.0)).abs() > 10.0);
}

#[test]
fn g0_native_cell_admission_rejects_incompatible_volume_without_publishing_birth() {
    let law = law([700.0; 3]);
    let mut bad = cell(0, 2.0, 0.0);
    bad.geometry.volume_m3 *= 2.0;
    let initial = simulation(&law, vec![bad]);
    assert!(matches!(initial.advance_to(0.0, config(), || Ok(())),
        Err(DepositionError::InvalidInput("cell mass/volume differs from its Cv source density"))));
    assert!(initial.cells().is_empty());
}

#[test]
fn g0_fixed_volume_internal_energy_is_refused_as_open_stream_enthalpy() {
    let law = law([700.0; 3]);
    let (initial, _) = simulation(&law, vec![cell(0, 2.0, 0.0)])
        .advance_to(0.0, config(), || Ok(())).unwrap();
    let flow = Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
        MassBoundary::Inlet {
            id: 1, cell: CellId(0), kg_s: 0.1,
            material_id: MaterialId(1), specific_h_j_kg: 700.0,
        },
        MassBoundary::Outlet { id: 2, cell: CellId(0), kg_s: 0.1 },
    ], 2).unwrap();
    assert!(matches!(initial.advance_observed(&flow, 1.0, config(), || Ok(())),
        Err(DepositionError::InvalidInput(
            "mass flow requires an explicitly enthalpy-based storage law"
        ))));
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 0.0);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g0_relative_retained_internal_energy_keeps_sub_temperature_ulp_heat_and_source_limits() {
    let law = law([500.0, 700.0, 1000.0]);
    let reference = 400.0_f64;
    let tiny =
        law.specific_heat_cv_j_kg_k(reference).unwrap() * (reference.next_up() - reference) * 0.125;
    assert!(tiny > 0.0);
    assert_eq!(
        law.temperature_at_relative_internal_energy(reference, tiny)
            .unwrap(),
        reference
    );
    assert!(
        law.temperature_after_specific_heat(reference, tiny)
            .is_err()
    );
    for target in [250.0, 300.0, 350.0, 425.0, 500.0] {
        let energy = law
            .specific_internal_energy_change_j_kg(reference, target)
            .unwrap();
        close(
            law.temperature_at_relative_internal_energy(reference, energy)
                .unwrap(),
            target,
        );
    }
    let upper = law
        .specific_internal_energy_change_j_kg(reference, 500.0)
        .unwrap();
    let lower = law
        .specific_internal_energy_change_j_kg(reference, 250.0)
        .unwrap();
    assert!(
        law.temperature_at_relative_internal_energy(reference, upper.next_up())
            .is_err()
    );
    assert!(
        law.temperature_at_relative_internal_energy(reference, lower.next_down())
            .is_err()
    );
    assert!(
        law.temperature_at_relative_internal_energy(249.0, 0.0)
            .is_err()
    );
    assert!(
        law.temperature_at_relative_internal_energy(reference, f64::NAN)
            .is_err()
    );
}
