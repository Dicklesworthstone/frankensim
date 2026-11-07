//! G1/G4: actual source contacts on original nonlinear sensible storage.
//! Synthetic sources; the heat oracle uses independent quadratic primitives.

use fs_conduction::transient::deposition::contact::implicit::{
    PreparedCaloricContact, PreparedContactError,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard,
    MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::{
    conduction::{LayerConductionConfig, LayerGeometry},
    interface::ThermalInterfacePoint,
    transport::{
        CONDUCTIVITY, DENSITY, SPECIFIC_HEAT,
        conduction::{CaloricConductionLayer, IntegratedCaloricLink},
        isobaric::{IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery, conduction::IsobaricConductionLaw},
        isochoric::{IsochoricHeatCapacityLaw, IsochoricHeatCapacityQuery, conduction::IsochoricConductionLaw},
        storage::{
            HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw, PreparedHeatStorage,
            conduction::HeatStorageConductionLaw,
        },
        volumetric::VOLUMETRIC_HEAT_CAPACITY_DIMS,
    },
};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink,
    ThermalStorageBasis, TimedEvent,
    conduction::{
        ConductionObserver, ConductivityLaw,
        implicit::{
            ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
            ImplicitConductionWorkspace, link::ImplicitLinkLaw,
        },
    },
};

const AREA: f64 = 0.02;
const MASSES: [f64; 2] = [0.4, 0.2];
const INITIAL_T: [f64; 2] = [650.0, 330.0];
const DISTANCES: [f64; 2] = [0.01, 0.005];
const K: [f64; 2] = [12.0, 8.0];
const SLOPES: [f64; 2] = [1e-5, 2e-5];
const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn energy(side: usize, t: f64) -> f64 {
    let x = t - 300.0;
    [1.0, 2.0][side] * (350.0 * x + 0.25 * x * x)
}
fn temperature(side: usize, e: f64) -> f64 {
    let e = e / [1.0, 2.0][side];
    300.0 + 2.0 * e / (350.0 + (350.0 * 350.0 + e).sqrt())
}
fn capacity(side: usize, t: f64) -> f64 {
    [1.0, 2.0][side] * (200.0 + 0.5 * t)
}
fn resistance(t: [f64; 2]) -> f64 {
    0.004 + SLOPES[0] * (t[0] - 300.0) + SLOPES[1] * (t[1] - 300.0)
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}");
}
fn bulk_point(side: usize, t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
        .with_quantity(
            if side == 0 { "rho" } else { "p" },
            QuantitySpec::dimensional(if side == 0 { DENSITY } else { Pressure::DIMS }),
            if side == 0 { 1000.0 } else { 100000.0 },
        ).unwrap()
}
fn caloric_sources() -> ([MaterialCard; 2], IsochoricHeatCapacityLaw, IsobaricHeatCapacityLaw) {
    let cards = [0, 1].map(|side| {
        let mut claims = ClaimSet::new();
        let validity = ValidityDomain::unconstrained()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0, 700.0)
            .with_quantity(
                if side == 0 { "rho" } else { "p" },
                QuantitySpec::dimensional(if side == 0 { DENSITY } else { Pressure::DIMS }),
                if side == 0 { 1000.0 } else { 100000.0 },
                if side == 0 { 1000.0 } else { 100000.0 },
            );
        for (name, dims, values) in [
            ("capacity", SPECIFIC_HEAT, [capacity(side, 300.0), capacity(side, 475.0), capacity(side, 700.0)]),
            ("k", CONDUCTIVITY, [K[side]; 3]),
            ("volume-capacity", VOLUMETRIC_HEAT_CAPACITY_DIMS,
                [1000.0 * capacity(side, 300.0), 1000.0 * capacity(side, 475.0), 1000.0 * capacity(side, 700.0)]),
        ] {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: [300.0, 475.0, 700.0].into_iter().zip(values).collect(), dims,
                },
                validity: validity.clone(), interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated, observations: vec![],
                provenance: Provenance { source: "authored nonlinear contact fixture".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
        }
        MaterialCard::assemble(MaterialStateId {
            chemistry: if side == 0 { "Cv body".into() } else { "Cp body".into() },
            phase: "single phase".into(), process: "fixed sensible body".into(), revision: 0,
        }, claims, vec![]).unwrap()
    });
    let key = PropertyKey::new("capacity", SPECIFIC_HEAT);
    let cv = IsochoricHeatCapacityLaw::try_new(&cards[0], IsochoricHeatCapacityQuery {
        specific_heat_cv: &key, temperature_axis: "T", density_axis: "rho",
        lower: &bulk_point(0, 300.0), upper: &bulk_point(0, 700.0), selection: SELECTION,
    }, 300.0).unwrap();
    let cp = IsobaricHeatCapacityLaw::try_new(&cards[1], IsobaricHeatCapacityQuery {
        specific_heat_cp: &key, temperature_axis: "T", pressure_axis: "p",
        lower: &bulk_point(1, 300.0), upper: &bulk_point(1, 700.0), selection: SELECTION,
    }, 300.0).unwrap();
    (cards, cv, cp)
}
fn contact_source(cards: &[MaterialCard; 2], complete: bool) -> (InterfaceSystemCard, ClaimId) {
    let mut claims = ClaimSet::new();
    let pin = claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(RP, RD),
        value: PropertyValue::Surface2D {
            axes: ["aT".into(), "bT".into()], axis_dims: [Temperature::DIMS; 2],
            samples: [[300.0, 300.0], [700.0, 300.0], [300.0, 700.0], [700.0, 700.0]]
                .into_iter().map(|t| (t, resistance(t))).collect(),
            triangles: if complete { vec![[0, 1, 2], [1, 3, 2]] } else { vec![[0, 1, 2]] },
            dims: RD,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("aT", QuantitySpec::dimensional(Temperature::DIMS), 300.0, 700.0)
            .with_quantity("bT", QuantitySpec::dimensional(Temperature::DIMS), 300.0, 700.0)
            .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 200000.0, 200000.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated, observations: vec![],
        provenance: Provenance { source: "authored ordered two-temperature contact".into(),
            license: "test-only".into(), artifact: None },
    }).unwrap();
    (super::common::interface(cards, claims), pin)
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("aT", QuantitySpec::dimensional(Temperature::DIMS), 400.0).unwrap()
        .with_quantity("bT", QuantitySpec::dimensional(Temperature::DIMS), 400.0).unwrap()
        .with_quantity("mean", QuantitySpec::dimensional(Temperature::DIMS), 400.0).unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 200000.0).unwrap()
}
fn point<'a>(cards: &'a [MaterialCard; 2], source: &'a InterfaceSystemCard, context: &'a QueryPoint)
    -> ThermalInterfacePoint<'a> {
    ThermalInterfacePoint { interface: source, materials: [&cards[0], &cards[1]], context,
        temperature_axes: ["aT", "bT"], liquid_fraction_axes: [None; 2],
        mean_temperature_axis: Some("mean") }
}
fn timeline<'a>(laws: [&'a dyn ThermalLaw; 2], energies: [f64; 2], masses: [f64; 2])
    -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    let mut start = 0.0;
    for side in 0..2 {
        let length = masses[side] / (1000.0 * AREA);
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(side as u64), segment_id: side as u64, body_id: side as u64,
            material_id: MaterialId(side as u64), mass_kg: masses[side],
            specific_h_j_kg: energies[side], geometry: CellGeometry {
                start_m: [start, 0.0, 0.0], end_m: [start + length, 0.0, 0.0],
                volume_m3: masses[side] / 1000.0,
                surface_area_m2: 2.0 * AREA + 4.0 * AREA.sqrt() * length,
            },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }) });
        start += length;
    }
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 7, a: CellId(0), b: CellId(1), area_m2: AREA, conductance_w_k: 0.0,
    }) });
    DepositionSimulation::new(events, &[(MaterialId(0), laws[0]), (MaterialId(1), laws[1])],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 32 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 }, || Ok(())).unwrap().0
}

// A single transferred heat q closes both original energy inventories. The
// reference never invokes a production source, caloric inverse or layer kernel.
fn exact_transfer(dt: f64, frozen_contact: bool) -> f64 {
    let e0 = [energy(0, INITIAL_T[0]), energy(1, INITIAL_T[1])];
    let mut lo = 0.0;
    let mut hi = MASSES[0] * (e0[0] - energy(0, INITIAL_T[1]));
    for _ in 0..80 {
        let q = (lo + hi) / 2.0;
        let t = [temperature(0, e0[0] - q / MASSES[0]), temperature(1, e0[1] + q / MASSES[1])];
        let r = resistance(if frozen_contact { INITIAL_T } else { t });
        let power = AREA * (t[0] - t[1]) / (DISTANCES[0] / K[0] + DISTANCES[1] / K[1] + r);
        if q - dt * power < 0.0 { lo = q; } else { hi = q; }
    }
    (lo + hi) / 2.0
}

#[test]
fn original_cv_cp_contact_closes_native_heat_and_entropy_with_exact_source_tangents() {
    let (cards, cv, cp) = caloric_sources();
    let key = PropertyKey::new("k", CONDUCTIVITY);
    let a = IsochoricConductionLaw::try_new(&cards[0], &cv, &key, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&cards[1], &cp, &key, SELECTION).unwrap();
    let layers = [CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap()];
    let e0 = [energy(0, INITIAL_T[0]), energy(1, INITIAL_T[1])];
    let old = timeline([&cv, &cp], e0, MASSES);
    let (source, pin) = contact_source(&cards, true);
    let context = context();
    let contact = PreparedCaloricContact::new(&old, 7, [&layers[0], &layers[1]],
        point(&cards, &source, &context), ClaimSelection::Pinned(pin),
        layers.each_ref().map(|layer| layer.energy_bounds_j_kg()), || Ok(())).unwrap();
    assert_eq!(contact.reference().temperatures_k, INITIAL_T);
    assert_eq!(contact.reference().specific_energies_j_kg, e0);
    assert_eq!(contact.reference().interface.properties()[0].answer().receipt.selected, pin);
    assert_eq!(contact.reference().interface.state_coordinate("mean"), Some(490.0));
    assert!(contact.source_surface().is_some() && contact.source_line().is_none());
    assert!(std::ptr::addr_eq(contact.caloric_laws()[0], &cv as &dyn ThermalLaw));
    assert_eq!(context.axes()["aT"], 400.0);
    let state = contact.state_at_specific_energies(e0).unwrap();
    for side in 0..2 {
        near(state.resistance_enthalpy_derivatives_m2_kg_k_w_j[side],
            SLOPES[side] / capacity(side, INITIAL_T[side]), 1e-18);
    }
    let link = IntegratedCaloricLink::new_with_contact(&old, 7, [&layers[0], &layers[1]],
        LayerGeometry::planar(DISTANCES, AREA).unwrap(), &contact,
        LayerConductionConfig { max_iterations: 64, absolute_heat_tolerance_w: 1e-10,
            relative_heat_tolerance: 1e-12 }, || Ok(())).unwrap();
    assert_eq!(link.native_geometry().contact_resistance_m2_k_w, 0.0);
    let flux = link.state_at_specific_enthalpies(e0, link.native_geometry(),
        old.contact_patch(7).unwrap(), &mut || Ok(())).unwrap();
    let d = DISTANCES[0] / K[0] + DISTANCES[1] / K[1] + resistance(INITIAL_T);
    let jump = INITIAL_T[0] - INITIAL_T[1];
    for side in 0..2 {
        let expected = AREA * ([d, -d][side] - jump * SLOPES[side])
            / (d * d * capacity(side, INITIAL_T[side]));
        near(flux.power_enthalpy_derivatives_w_kg_j[side], expected, 1e-12);
    }
    let observer = ConductionObserver::new(&old,
        &[(MaterialId(0), &a as &dyn ConductivityLaw), (MaterialId(1), &b as &dyn ConductivityLaw)],
        &[link.native_geometry()], &[], || Ok(())).unwrap();
    let chain = ImplicitConductionChain::new(&old, &observer,
        &[(MaterialId(0), &a as &dyn ImplicitConductionLaw), (MaterialId(1), &b as &dyn ImplicitConductionLaw)],
        || Ok(())).unwrap().with_link_laws(&old, &[(7, &link)], || Ok(())).unwrap();
    let dt = 120.0;
    let controls = ImplicitConductionControls { maximum_dt_s: dt, maximum_halvings: 0,
        cell_tolerance_j: 1e-6, total_tolerance_j: 1e-6, ..ImplicitConductionControls::default() };
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, report) = chain.advance_to(&old, dt, controls, &mut workspace, || Ok(())).unwrap();
    let q = exact_transfer(dt, false);
    let e = [next.cells()[&CellId(0)].specific_h_j_kg, next.cells()[&CellId(1)].specific_h_j_kg];
    near(e[0], e0[0] - q / MASSES[0], 5e-6);
    near(e[1], e0[1] + q / MASSES[1], 5e-6);
    near(next.energy_j(), old.energy_j(), 1e-6);
    assert!((q - exact_transfer(dt, true)).abs() > 500.0);
    let t = [temperature(0, e[0]), temperature(1, e[1])];
    let entropy = MASSES[0] * (200.0 * (t[0] / INITIAL_T[0]).ln() + 0.5 * (t[0] - INITIAL_T[0]))
        + MASSES[1] * (400.0 * (t[1] / INITIAL_T[1]).ln() + t[1] - INITIAL_T[1]);
    near(report.entropy.thermal_storage_change_j_k, entropy, 2e-8);
    near(report.entropy.conduction_generation_j_k, q * (1.0 / t[1] - 1.0 / t[0]), 2e-8);
    assert!(entropy > 0.0 && report.energy.residual_j.abs() <= 1e-6);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(next.material_storage_basis(MaterialId(0)).unwrap(), ThermalStorageBasis::InternalEnergy);
    assert_eq!(next.material_storage_basis(MaterialId(1)).unwrap(), ThermalStorageBasis::Enthalpy);
    assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, e0[0]);
}

#[test]
fn sensible_contact_keeps_original_identities_source_holes_and_transactional_preparation() {
    let (cards, cv, cp) = caloric_sources();
    let key = PropertyKey::new("k", CONDUCTIVITY);
    let a = IsochoricConductionLaw::try_new(&cards[0], &cv, &key, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&cards[1], &cp, &key, SELECTION).unwrap();
    let layers = [CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap()];
    let old = timeline([&cv, &cp], [energy(0, INITIAL_T[0]), energy(1, INITIAL_T[1])], MASSES);
    let (source, pin) = contact_source(&cards, false);
    let context = context();
    let point = point(&cards, &source, &context);
    let bounds = layers.each_ref().map(|layer| layer.energy_bounds_j_kg());
    let contact = PreparedCaloricContact::new(&old, 7, [&layers[0], &layers[1]], point,
        ClaimSelection::Pinned(pin), bounds, || Ok(())).unwrap();
    assert!(contact.state_at_specific_energies([energy(0, 650.0), energy(1, 650.0)]).is_err());
    assert!(matches!(contact.state_at_specific_energies([-1.0, 0.0]),
        Err(PreparedContactError::OutsideEnergyDomain { side: 0, .. })));
    let foreign = MaterialCard::assemble(cards[0].id().clone(), ClaimSet::new(), vec![]).unwrap();
    for bad in [
        ThermalInterfacePoint { materials: [&foreign, &cards[1]], ..point },
        ThermalInterfacePoint { liquid_fraction_axes: [Some("fa"), None], ..point },
        ThermalInterfacePoint { temperature_axes: ["aT", "aT"], ..point },
    ] {
        assert!(PreparedCaloricContact::new(&old, 7, [&layers[0], &layers[1]], bad,
            ClaimSelection::Pinned(pin), bounds, || Ok(())).is_err());
    }
    assert!(PreparedCaloricContact::new(&old, 7, [&layers[1], &layers[0]], point,
        ClaimSelection::Pinned(pin), bounds, || Ok(())).is_err());
    let copied = cv.clone();
    let copied_conduction = IsochoricConductionLaw::try_new(&cards[0], &copied, &key, SELECTION).unwrap();
    let copied_layer = CaloricConductionLayer::from_isochoric(&copied_conduction).unwrap();
    assert!(PreparedCaloricContact::new(&old, 7, [&copied_layer, &layers[1]], point,
        ClaimSelection::Pinned(pin), bounds, || Ok(())).is_err());
    let before = old.cells().clone();
    let mut polls = 0;
    let cancelled = PreparedCaloricContact::new(&old, 7, [&layers[0], &layers[1]], point,
        ClaimSelection::Pinned(pin), bounds, || {
            polls += 1;
            if polls == 5 { Err(DepositionError::Cancelled) } else { Ok(()) }
        });
    assert!(matches!(cancelled, Err(PreparedContactError::Thermal(DepositionError::Cancelled))));
    let retry = PreparedCaloricContact::new(&old, 7, [&layers[0], &layers[1]], point,
        ClaimSelection::Pinned(pin), bounds, || Ok(())).unwrap();
    assert_eq!(retry.reference(), contact.reference());
    assert_eq!(old.cells(), &before);
}

#[test]
fn assembled_contact_uses_total_mass_and_volume_capacity_in_its_native_energy_coordinate() {
    let (cards, _, cp) = caloric_sources();
    let capacity_key = PropertyKey::new("capacity", SPECIFIC_HEAT);
    let volume_key = PropertyKey::new("volume-capacity", VOLUMETRIC_HEAT_CAPACITY_DIMS);
    let (lower, upper) = (bulk_point(0, 300.0), bulk_point(0, 700.0));
    let component = |id, capacity, amount| HeatStorageComponentQuery {
        id, card: &cards[0], capacity, amount, temperature_axis: "T",
        lower: &lower, upper: &upper, selection: SELECTION,
    };
    // C_total(T) = (0.4 kg + 0.0002 m³ * 1000 kg/m³) * Cv(T).
    // The declared accounting mass is 2 kg, not the sum of constituent masses.
    let storage = PreparedHeatStorage::try_new(&[
        component(1, &capacity_key, HeatStorageAmount::MassKg(0.4)),
        component(2, &volume_key, HeatStorageAmount::ReferenceVolumeM3(0.0002)),
    ]).unwrap();
    let assembly = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0).unwrap();
    let key = PropertyKey::new("k", CONDUCTIVITY);
    let a = HeatStorageConductionLaw::try_new(&cards[0], &assembly, &key, "T", &lower, &upper,
        SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&cards[1], &cp, &key, SELECTION).unwrap();
    let layers = [CaloricConductionLayer::from_heat_storage(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap()];
    let energies = [0.3 * energy(0, INITIAL_T[0]), energy(1, INITIAL_T[1])];
    let old = timeline([&assembly, &cp], energies, [2.0, MASSES[1]]);
    let (source, pin) = contact_source(&cards, true);
    let context = context();
    let contact = PreparedCaloricContact::new(&old, 7, [&layers[0], &layers[1]],
        point(&cards, &source, &context), ClaimSelection::Pinned(pin),
        layers.each_ref().map(|layer| layer.energy_bounds_j_kg()), || Ok(())).unwrap();
    let response = contact.state_at_specific_energies(energies).unwrap();
    near(response.resistance_m2_k_w, resistance(INITIAL_T), 1e-15);
    near(response.resistance_enthalpy_derivatives_m2_kg_k_w_j[0],
        SLOPES[0] * 2.0 / (0.6 * capacity(0, INITIAL_T[0])), 1e-18);
    assert_eq!(contact.caloric_laws()[0].storage_basis(), ThermalStorageBasis::DeclaredSensibleStorage);
    assert_eq!(contact.energy_bounds_j_kg()[0], assembly.specific_storage_span_j_kg());
    assert!(layers[0].supports_material_card(cards[0].content_hash()));
    assert!(!layers[0].supports_material_card(cards[1].content_hash()));
}
