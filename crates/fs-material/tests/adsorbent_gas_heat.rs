//! Independent manufactured caloric/relaxation references; no empirical fit.
use std::{collections::BTreeMap, sync::Arc};

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    boundary::GasSpeciesBoundarySolve,
    diffusion::dusty_gas::{DustyGasFaceSolve, DustyGasMedium, DustyGasPores},
    network::{
        SpeciesNetworkSolve,
        transport::{
            GasTransportPaths, GasTransportSolve,
            pore_domain::{
                PoreCellGeometry, PoreDomain, PoreDomainError, PoreDomainSolve, PoreFace,
                adsorption::AdsorbingPoreDomain,
            },
        },
    },
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_material::{
    adsorption::{
        AdsorptionCaloricBinding, AdsorptionSource, LANGMUIR_AFFINITY_DIMS, LANGMUIR_CAPACITY_DIMS,
        LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS,
        LANGMUIR_REFERENCE_AFFINITY_PROPERTY, LangmuirAdsorbent,
        competitive::{
            CompetitiveAdsorbent, CompetitiveSolid, CompetitiveSolve,
            network::{
                self, AdsorptionNetworkBinding, AdsorptionNetworkError,
                gas_heat::{
                    AdsorbentGasHeatLink, AdsorbentGasHeatSource,
                    GAS_SOLID_THERMAL_CONDUCTANCE_DIMS,
                },
            },
        },
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    porous::permeability::{
        DarcyGasBarrier, IntrinsicPermeability, PorousGasLayer, PorousGasSupport,
    },
    state_point::MaterialPropertySelection,
    thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel},
};
use fs_qty::{
    Dims, DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, mixture::FrozenIdealGasMixtureModelV1,
    required_nasa9_parameter_dims,
};

const HEAT: &str = "gas_solid_thermal_conductance";
const MASS: f64 = 0.028;
const CP_R: f64 = 3.5;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p)
        .unwrap()
}
fn heat_point(tg: f64, ts: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "gas_temperature",
            typed(QuantityKind::AbsoluteTemperature),
            tg,
        )
        .unwrap()
        .with_quantity(
            "solid_temperature",
            typed(QuantityKind::AbsoluteTemperature),
            ts,
        )
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p)
        .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured closed sensible heat reference".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                200.0,
                800.0,
            )
            .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: provenance(),
    }
}
fn card(name: &str, claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "manufactured".into(),
            process: "analytic thermal/adsorption fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}
fn gas_model() -> ThermochemicalGasModel {
    let mut coefficients = [0.0; 9];
    coefficients[2] = CP_R;
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let law = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 800.0),
        sources: vec![],
        provenance: provenance(),
    };
    let species = SpeciesId::new("A").unwrap();
    let nasa = Nasa9StandardStateModelV1::new(
        species.clone(),
        MolarMass::new(MASS),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("analytic-reference").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(law).unwrap()],
    )
    .unwrap();
    let transport = SpeciesTransportMaterial {
        species,
        material: card(
            "gas A",
            vec![
                claim("dynamic_viscosity", DynViscosity::DIMS, 1.7e-5),
                claim("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), 0.025),
            ],
        ),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![nasa],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport],
        None,
    )
    .unwrap()
}
fn chart(name: &str, cp: f64) -> ResolvedHeatCapacityCurve {
    let c = Dims([2, 0, -2, -1, 0, 0]);
    let r = Dims([-3, 1, 0, 0, 0, 0]);
    let source = card(
        name,
        vec![
            claim("specific_heat_capacity", c, cp),
            claim("density", r, 1000.0),
        ],
    );
    ResolvedHeatCapacityCurve::try_new(
        &source,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", c),
            density: &PropertyKey::new("density", r),
            temperature_axis: "temperature",
            lower: &point(200.0, 1e5),
            upper: &point(800.0, 1e5),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap()
}
fn fixture(
    ts: f64,
    rate: f64,
    loading: f64,
) -> (CompetitiveSolid, ClosedGasParcel, ThermochemicalGasModel) {
    let catalog = gas_model();
    let species = SpeciesId::new("A").unwrap();
    let dry = chart("dry skeleton", 500.0);
    let bound = chart("bound A", CP_R * R / MASS);
    let offset = catalog
        .species_specific_enthalpy_j_kg(&species, 300.0)
        .unwrap()
        - bound
            .state_at_temperature(300.0)
            .unwrap()
            .specific_enthalpy_j_kg()
        - 6000.0 / MASS;
    let law = LangmuirAdsorbent::try_new(
        &catalog,
        species.clone(),
        AdsorptionSource {
            material: card(
                "sites",
                vec![
                    claim(LANGMUIR_CAPACITY_PROPERTY, LANGMUIR_CAPACITY_DIMS, 0.5),
                    claim(
                        LANGMUIR_REFERENCE_AFFINITY_PROPERTY,
                        LANGMUIR_AFFINITY_DIMS,
                        2e-5,
                    ),
                    claim(LANGMUIR_DESORPTION_RATE_PROPERTY, LANGMUIR_RATE_DIMS, rate),
                ],
            ),
            point: point(300.0, 1e5),
            selection: MaterialPropertySelection::SingleClaimOnly,
            reference_temperature_k: 300.0,
        },
        AdsorptionCaloricBinding {
            dry_skeleton: dry,
            bound_adsorbate: bound,
            bound_enthalpy_offset_j_kg: offset,
        },
    )
    .unwrap();
    let solid = CompetitiveSolid::try_new(
        Arc::new(CompetitiveAdsorbent::try_new(vec![law]).unwrap()),
        0.1,
        &[(species, loading)],
        ts,
    )
    .unwrap();
    let gas = ClosedGasParcel::new(catalog.clone(), point(300.0, 1e5), 0.001).unwrap();
    (solid, gas, catalog)
}
fn local() -> CompetitiveSolve {
    CompetitiveSolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 800.0,
                energy_tolerance_j_kg: 1e-8,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        maximum_gas_gross_fraction: 0.1,
        maximum_site_turnover: 0.25,
        maximum_gas_surface_temperature_difference_k: 5.0,
        maximum_temperature_change_k: 2.0,
        entropy_tolerance_j_k: 1e-10,
        maximum_iterations: 128,
    }
}
fn transport() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 2e-9,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}
fn heat_claim(g: f64, solid_range: [f64; 2]) -> PropertyClaim {
    let mut source = claim(HEAT, GAS_SOLID_THERMAL_CONDUCTANCE_DIMS, g);
    source.validity = ValidityDomain::unconstrained()
        .with_quantity(
            "gas_temperature",
            typed(QuantityKind::AbsoluteTemperature),
            290.0,
            310.0,
        )
        .with_quantity(
            "solid_temperature",
            typed(QuantityKind::AbsoluteTemperature),
            solid_range[0],
            solid_range[1],
        )
        .with_quantity("pressure", typed(QuantityKind::Pressure), 90000.0, 110000.0);
    source
}
fn source(g: f64) -> AdsorbentGasHeatSource {
    AdsorbentGasHeatSource::try_new(
        card("exposed path", vec![heat_claim(g, [290.0, 310.0])]),
        &PropertyKey::new(HEAT, GAS_SOLID_THERMAL_CONDUCTANCE_DIMS),
        &heat_point(300.0, 300.0, 1e5),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn sensible_heat_uses_original_gas_cv_and_dry_plus_bound_storage_in_both_directions() {
    let source = source(0.5);
    for ts in [299.0, 301.0] {
        let (solid, gas, catalog) = fixture(ts, 0.0, 0.005);
        let saved = (solid.clone(), gas.clone());
        let parcels = BTreeMap::from([(4, gas)]);
        let endpoints = BTreeMap::from([(4, local().gas.gas)]);
        let bindings = [AdsorptionNetworkBinding {
            id: 7,
            parcel: 4,
            solid: &solid,
            solve: local(),
        }];
        let links = [AdsorbentGasHeatLink {
            id: 2,
            solid: 7,
            source: &source,
        }];
        let candidate = network::step_candidate_with_heat(
            &parcels,
            GasTransportPaths::default(),
            &bindings,
            &catalog,
            0.125,
            &endpoints,
            transport(),
            &[],
            &links,
            || false,
        )
        .unwrap();
        let q = 0.5 * (300.0 - ts) * 0.125;
        let cg = parcels[&4].mass_kg() / MASS * (CP_R - 1.0) * R;
        let cs = 0.1 * 500.0 + 0.005 * CP_R * R;
        near(
            candidate.parcels()[&4].state().state().temperature,
            300.0 - q / cg,
            2e-10,
        );
        near(candidate.solids()[&7].temperature_k(), ts + q / cs, 2e-10);
        assert_eq!(candidate.parcels()[&4].mass_kg(), parcels[&4].mass_kg());
        assert_eq!(
            candidate.solids()[&7].adsorbed_moles(&SpeciesId::new("A").unwrap()),
            Some(0.005)
        );
        near(
            candidate.step().exchange.parcels[0].net_additional_heat_j,
            -q,
            0.0,
        );
        near(
            candidate.step().exchange.parcels[0].net_enthalpy_j,
            0.0,
            0.0,
        );
        near(candidate.step().solids[0].gas_sensible_heat_j, q, 0.0);
        near(candidate.step().solids[0].conductive_heat_j, 0.0, 0.0);
        near(
            candidate.step().exchange.parcels[0].energy_change_j
                + candidate.step().solids[0].solid_energy_change_j,
            0.0,
            2e-9,
        );
        near(candidate.step().total_energy_residual_j, 0.0, 2e-9);
        assert!(candidate.step().entropy_change_j_k > 0.0);
        assert_eq!(candidate.step().exchange.external_heat.len(), 1);
        assert_eq!((solid, parcels[&4].clone()), saved);
    }
}

#[test]
fn material_coefficients_refresh_actual_temperatures_and_preserve_pinned_claims() {
    let mut c = heat_claim(0.0, [290.0, 310.0]);
    c.value = PropertyValue::Curve {
        abscissa: "gas_temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(290.0, 0.1), (310.0, 0.3)],
        dims: GAS_SOLID_THERMAL_CONDUCTANCE_DIMS,
    };
    c.interpolation = InterpolationPolicy::LinearInside;
    let source = AdsorbentGasHeatSource::try_new(
        card("temperature-dependent path", vec![c]),
        &PropertyKey::new(HEAT, GAS_SOLID_THERMAL_CONDUCTANCE_DIMS),
        &heat_point(290.0, 300.0, 1e5),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let (solid, gas, catalog) = fixture(301.0, 0.0, 0.0);
    let parcels = BTreeMap::from([(0, gas)]);
    let endpoints = BTreeMap::from([(0, local().gas.gas)]);
    let bindings = [AdsorptionNetworkBinding {
        id: 7,
        parcel: 0,
        solid: &solid,
        solve: local(),
    }];
    let links = [
        AdsorbentGasHeatLink {
            id: 8,
            solid: 7,
            source: &source,
        },
        AdsorbentGasHeatLink {
            id: 3,
            solid: 7,
            source: &source,
        },
    ];
    let run = |links: &[AdsorbentGasHeatLink<'_>]| {
        network::step_candidate_with_heat(
            &parcels,
            GasTransportPaths::default(),
            &bindings,
            &catalog,
            0.125,
            &endpoints,
            transport(),
            &[],
            links,
            || false,
        )
        .unwrap()
    };
    let a = run(&links);
    assert_eq!(a, run(&[links[1], links[0]]));
    assert_eq!(
        a.step().gas_heat.iter().map(|s| s.id).collect::<Vec<_>>(),
        vec![3, 8]
    );
    for sample in &a.step().gas_heat {
        near(sample.conductance_w_k, 0.2, 1e-16);
        assert_eq!(
            sample.source[0].state_coordinate("gas_temperature"),
            Some(300.0)
        );
        assert_eq!(
            sample.source[0].state_coordinate("solid_temperature"),
            Some(301.0)
        );
        assert_eq!(
            sample.source[1].state_coordinate("gas_temperature"),
            Some(a.parcels()[&0].state().state().temperature)
        );
        assert_eq!(
            sample.source[1].state_coordinate("pressure"),
            Some(a.parcels()[&0].state().state().pressure)
        );
        assert!(sample.source[1].property(HEAT).unwrap().value_si() > sample.conductance_w_k);
        assert_eq!(
            sample.source[0]
                .property(HEAT)
                .unwrap()
                .answer()
                .receipt
                .selected,
            sample.source[1]
                .property(HEAT)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    // No invented default heat transfer: zero G recovers the preexisting owner.
    let zero = super_source_zero();
    let zero_links = [AdsorbentGasHeatLink {
        id: 1,
        solid: 7,
        source: &zero,
    }];
    let old = network::step_candidate(
        &parcels,
        GasTransportPaths::default(),
        &bindings,
        &catalog,
        0.125,
        &endpoints,
        transport(),
        || false,
    )
    .unwrap();
    let new = run(&zero_links);
    assert_eq!(old.parcels(), new.parcels());
    assert_eq!(old.solids(), new.solids());
    assert!(new.step().exchange.external_heat.is_empty());
}
fn super_source_zero() -> AdsorbentGasHeatSource {
    source(0.0)
}

fn medium() -> DustyGasMedium {
    let layer = PorousGasLayer::try_new(
        IntrinsicPermeability::isotropic(2e-14).unwrap(),
        0.01,
        0.0,
        "manufactured rigid pore geometry",
    )
    .unwrap();
    let barrier = DarcyGasBarrier::try_new(
        vec![layer],
        [1.0, 0.0, 0.0],
        PorousGasSupport {
            temperature_k: (200.0, 800.0),
            pressure_pa: (1000.0, 1e7),
            maximum_superficial_speed_m_s: 100.0,
        },
    )
    .unwrap();
    DustyGasMedium::try_new(
        barrier,
        DustyGasPores::try_new(1e-7, 0.4, 2.0, "manufactured voids").unwrap(),
    )
    .unwrap()
}
fn domain(medium: &DustyGasMedium, ts: f64, rate: f64, loading: f64) -> AdsorbingPoreDomain {
    let (solid, gas, catalog) = fixture(ts, rate, loading);
    AdsorbingPoreDomain::try_new(
        PoreDomain::try_new(
            0.25,
            vec![PoreCellGeometry::try_new(medium.clone()).unwrap()],
            gas.clone(),
            vec![gas.clone()],
            gas,
        )
        .unwrap(),
        BTreeMap::from([(0, solid)]),
        catalog,
    )
    .unwrap()
}
fn faces<'a>(medium: &'a DustyGasMedium, context: &'a QueryPoint, area: f64) -> [PoreFace<'a>; 2] {
    [PoreFace {
        medium,
        area_m2: area,
        pairs: &[],
        gas_context: context,
        solve: DustyGasFaceSolve {
            maximum_temperature_difference_k: 5.0,
            maximum_relative_pressure_difference: 0.25,
            gradient_tolerance_per_m: 1e-8,
        },
    }; 2]
}
fn window(dt: f64) -> PoreDomainSolve {
    PoreDomainSolve {
        maximum_step_s: dt,
        minimum_step_s: 1e-8,
        maximum_substeps: 1024,
        maximum_attempts: 4096,
        maximum_throat_volume_fraction: 0.02,
        transport: transport(),
    }
}

#[test]
fn closed_pore_heat_converges_to_independent_two_body_exponential_and_replays() {
    let medium = medium();
    let original = domain(&medium, 301.0, 0.0, 0.005);
    let saved = original.clone();
    let context = point(300.0, 1e5);
    let faces = faces(&medium, &context, 0.0);
    let source = source(0.5);
    let links = [AdsorbentGasHeatLink {
        id: 5,
        solid: 0,
        source: &source,
    }];
    let local = BTreeMap::from([(0, local())]);
    let endpoints = [local[&0].gas.gas; 3];
    let gas = original.pore_domain().pore_gas(0).unwrap();
    let cg = gas.mass_kg() / MASS * (CP_R - 1.0) * R;
    let cs = 50.0 + 0.005 * CP_R * R;
    let teq = (cg * 300.0 + cs * 301.0) / (cg + cs);
    let exact = teq + (300.0 - teq) * (-0.5 * (1.0 / cg + 1.0 / cs)).exp();
    let mut errors = Vec::new();
    for dt in [0.25, 0.125, 0.0625] {
        let candidate = original
            .step_candidate_with_heat(
                1.0,
                &faces,
                &endpoints,
                window(dt),
                &local,
                &[],
                &links,
                || false,
            )
            .unwrap();
        let t = candidate
            .domain()
            .pore_domain()
            .pore_gas(0)
            .unwrap()
            .state()
            .state()
            .temperature;
        errors.push((t - exact).abs());
        assert!(candidate.step().pore.entropy_change_j_k > 0.0);
        near(candidate.step().pore.total_energy_residual_j, 0.0, 2e-9);
        near(
            candidate.step().gas_heat[0].gross_solid_to_gas_j,
            cg * (t - 300.0),
            2e-9,
        );
        assert_eq!(candidate.step().gas_heat[0].gross_gas_to_solid_j, 0.0);
        assert_eq!(candidate.step().pore.accepted_substeps, (1.0 / dt) as usize);
    }
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!((1.8..2.3).contains(&ratio));
    }
    let whole = original
        .step_candidate_with_heat(
            1.0,
            &faces,
            &endpoints,
            window(0.125),
            &local,
            &[],
            &links,
            || false,
        )
        .unwrap();
    let mut split = original.clone();
    split
        .advance_with_heat(
            0.5,
            &faces,
            &endpoints,
            window(0.125),
            &local,
            &[],
            &links,
            || false,
        )
        .unwrap();
    split
        .advance_with_heat(
            0.5,
            &faces,
            &endpoints,
            window(0.125),
            &local,
            &[],
            &links,
            || false,
        )
        .unwrap();
    assert_eq!(&split, whole.domain());
    assert_eq!(original, saved);
}

#[test]
fn sensible_heat_changes_later_adsorption_and_actual_porous_transport() {
    let medium = medium();
    let original = domain(&medium, 301.0, 0.4, 0.0);
    let context = point(300.0, 1e5);
    let faces = faces(&medium, &context, 0.001);
    let source = source(200.0);
    let links = [AdsorbentGasHeatLink {
        id: 5,
        solid: 0,
        source: &source,
    }];
    let local = BTreeMap::from([(0, local())]);
    let endpoints = [local[&0].gas.gas; 3];
    let dt = 0.00025;
    let species = SpeciesId::new("A").unwrap();
    let first = original
        .step_candidate_with_heat(
            dt,
            &faces,
            &endpoints,
            window(dt),
            &local,
            &[],
            &links,
            || false,
        )
        .unwrap();
    let first_old = original
        .step_candidate(dt, &faces, &endpoints, window(dt), &local, || false)
        .unwrap();
    // Both mechanisms see the same original gas; this interval's kinetic and
    // throat rates are unchanged, with all sensible heat in the same inversion.
    assert_eq!(
        first.step().pore.last_substep.transport,
        first_old.step().pore.last_substep.transport
    );
    assert_eq!(
        first.domain().solids()[&0].adsorbed_moles(&species),
        first_old.domain().solids()[&0].adsorbed_moles(&species)
    );
    assert!(
        first
            .domain()
            .pore_domain()
            .pore_gas(0)
            .unwrap()
            .state()
            .state()
            .temperature
            > first_old
                .domain()
                .pore_domain()
                .pore_gas(0)
                .unwrap()
                .state()
                .state()
                .temperature
    );
    let next = original
        .step_candidate_with_heat(
            8.0 * dt,
            &faces,
            &endpoints,
            window(dt),
            &local,
            &[],
            &links,
            || false,
        )
        .unwrap();
    let old = original
        .step_candidate(8.0 * dt, &faces, &endpoints, window(dt), &local, || false)
        .unwrap();
    assert_ne!(
        next.domain().solids()[&0].adsorbed_moles(&species),
        old.domain().solids()[&0].adsorbed_moles(&species)
    );
    assert_ne!(next.step().pore.faces, old.step().pore.faces);
    assert!(next.step().pore.entropy_change_j_k > 0.0);
    near(next.step().pore.total_energy_residual_j, 0.0, 2e-9);
    near(next.step().pore.total_mass_residual_kg, 0.0, 1e-13);
    assert!(next.step().gas_heat[0].gross_solid_to_gas_j > 0.0);
}

#[test]
fn excessive_heat_retries_typed_entropy_failure_with_bounded_atomic_work() {
    let medium = medium();
    let original = domain(&medium, 301.0, 0.0, 0.0);
    let mut mutable = original.clone();
    let context = point(300.0, 1e5);
    let faces = faces(&medium, &context, 0.0);
    let source = source(1.0);
    let links = [AdsorbentGasHeatLink {
        id: 5,
        solid: 0,
        source: &source,
    }];
    let mut c = local();
    c.maximum_temperature_change_k = 10.0;
    c.maximum_gas_surface_temperature_difference_k = 10.0;
    let local = BTreeMap::from([(0, c)]);
    let endpoints = [c.gas.gas; 3];
    let mut controls = window(3.0);
    controls.maximum_attempts = 1;
    assert!(matches!(
        mutable.advance_with_heat(
            3.0,
            &faces,
            &endpoints,
            controls,
            &local,
            &[],
            &links,
            || false
        ),
        Err(PoreDomainError::WorkLimit { .. })
    ));
    assert_eq!(mutable, original);
    let mut polls = 0;
    let successful = original
        .step_candidate_with_heat(
            3.0,
            &faces,
            &endpoints,
            window(3.0),
            &local,
            &[],
            &links,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(successful.step().pore.attempts > successful.step().pore.accepted_substeps);
    assert!(successful.step().pore.entropy_change_j_k >= -1e-10);
    near(successful.step().pore.total_energy_residual_j, 0.0, 2e-9);
    let last = polls;
    let mut calls = 0;
    assert!(matches!(
        mutable.advance_with_heat(
            3.0,
            &faces,
            &endpoints,
            window(3.0),
            &local,
            &[],
            &links,
            || {
                calls += 1;
                calls >= last
            }
        ),
        Err(PoreDomainError::Cancelled)
    ));
    assert_eq!(mutable, original);
}

#[test]
fn unowned_source_axes_and_late_original_domain_refusal_publish_nothing() {
    let key = PropertyKey::new(HEAT, GAS_SOLID_THERMAL_CONDUCTANCE_DIMS);
    for (alias, quantity, value) in [
        ("T", typed(QuantityKind::AbsoluteTemperature), 300.0),
        ("p", typed(QuantityKind::Pressure), 1e5),
        ("coverage", QuantitySpec::dimensional(Dims::NONE), 0.25),
    ] {
        let mut c = heat_claim(1.0, [290.0, 310.0]);
        c.validity = c.validity.with_quantity(alias, quantity, value, value);
        let query = heat_point(300.0, 300.0, 1e5)
            .with_quantity(alias, quantity, value)
            .unwrap();
        assert!(
            AdsorbentGasHeatSource::try_new(
                card(alias, vec![c]),
                &key,
                &query,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    assert!(
        AdsorbentGasHeatSource::try_new(
            card("negative", vec![heat_claim(-1.0, [290.0, 310.0])]),
            &key,
            &heat_point(300.0, 300.0, 1e5),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    assert!(
        AdsorbentGasHeatSource::try_new(
            card("wrong units", vec![heat_claim(1.0, [290.0, 310.0])]),
            &PropertyKey::new(HEAT, Dims([1, 1, -3, -1, 0, 0])),
            &heat_point(300.0, 300.0, 1e5),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let medium = medium();
    let mut original = domain(&medium, 299.0, 0.0, 0.0);
    let saved = original.clone();
    let source = AdsorbentGasHeatSource::try_new(
        card(
            "narrow actual solid support",
            vec![heat_claim(1.0, [298.0, 299.0005])],
        ),
        &key,
        &heat_point(300.0, 299.0, 1e5),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let context = point(300.0, 1e5);
    let faces = faces(&medium, &context, 0.0);
    let links = [AdsorbentGasHeatLink {
        id: 5,
        solid: 0,
        source: &source,
    }];
    let local = BTreeMap::from([(0, local())]);
    let result = original.advance_with_heat(
        0.125,
        &faces,
        &[local[&0].gas.gas; 3],
        window(0.125),
        &local,
        &[],
        &links,
        || false,
    );
    assert!(
        matches!(
            result,
            Err(PoreDomainError::Adsorption {
                error: AdsorptionNetworkError::Conduction { id: 5, .. },
                ..
            })
        ),
        "{result:?}"
    );
    assert_eq!(original, saved);
}
