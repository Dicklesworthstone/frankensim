//! Manufactured original-source gas fixtures; no empirical material claim.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::sublimation::gas::{
    GasSurfaceBinding, GasSurfaceError, GasSurfaceSolve,
    chamber::{GasChamberSolve, coating::GasSupportedChamberSolve},
};
use fs_material::sublimation::gas::chamber::coating::resistance::{
    advance::GasResistiveAdvanceError,
    growing::advance_growing_integrated_chamber,
};
use fs_material::thermal::conduction::{ConductivityCurve, LayerConductionConfig};
use fs_material::thermochemical_gas::{
    SpeciesTransportMaterial,
    energy::{
        GasEnergySolve,
        parcel::transfer::species::boundary::GasSpeciesBoundarySolve,
    },
};
use fs_qty::{DynViscosity, MolarMass};
use fs_thermal::enthalpy::DepositionError;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::sublimation::{
    SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw, SublimationProperties,
    coating::resistance::{
        ResistiveSurfaceState, geometry::ConductiveCoatingGeometry,
    },
};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1,
    NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1, required_nasa9_parameter_dims,
};
use crate::graph::thermal::medium::mass_transfer::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, evaluate_laminar_plate_gas_film,
};

const LATENT: f64 = 200000.0;
const MOLAR_MASS: f64 = 0.018015;
const AREA: f64 = 1e-4;

fn provenance() -> Provenance {
    Provenance {
        source: "manufactured live binary-film fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p).unwrap()
}
fn material(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims { set.insert_claim(claim).unwrap(); }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored binary pair".into(),
            phase: "gas".into(),
            process: "manufactured original-source fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    ).unwrap()
}
fn scalar(
    property: &str, dims: fs_qty::Dims, axis: &str, knots: Vec<(f64, f64)>,
    temperatures: [f64; 2],
) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(property, dims),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: if axis == "temperature" { Temperature::DIMS } else { Pressure::DIMS },
            knots, dims,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature),
                temperatures[0], temperatures[1])
            .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e6),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: provenance(),
    }
}
fn transport(name: &str, mu: f64, k: f64) -> SpeciesTransportMaterial {
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: material(vec![
            scalar("dynamic_viscosity", DynViscosity::DIMS, "temperature",
                vec![(200.0, mu), (900.0, 2.0 * mu)], [200.0, 900.0]),
            scalar(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, "temperature",
                vec![(200.0, k), (900.0, 2.0 * k)], [200.0, 900.0]),
        ]),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn species(name: &str, mass: f64, cp_r: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    let mut parameters = BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter().take(9).zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert("reference_pressure".into(),
        LawParameter { value: 100000.0, dims: Pressure::DIMS });
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 900.0),
        sources: vec![],
        provenance: provenance(),
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(), MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas, ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("manufactured-reference").unwrap(),
        ).unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    ).unwrap()
}
fn catalog() -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("A", 0.02897, 3.5), species("B", MOLAR_MASS, 4.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![0.97, 0.03]).unwrap(),
        ).unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.03)],
        Some(&SpeciesId::new("B").unwrap()),
    ).unwrap()
}
fn gas(catalog: &ThermochemicalGasModel, t: f64, p: f64, vapor: f64) -> ClosedGasParcel {
    let model = catalog.with_composition(
        Composition::new(CompositionBasis::MoleFraction, vec![1.0-vapor, vapor]).unwrap(),
    ).unwrap();
    ClosedGasParcel::new(
        model, point(t, p).with("local-lot", 9.0).unwrap(), 0.01,
    ).unwrap()
}
#[derive(Clone, Copy)]
struct Phase;
impl SublimationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] { [200.0, 900.0] }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let r = R / MOLAR_MASS;
        let pressure = 1000.0 * fs_math::det::exp(LATENT/r * (1.0/400.0 - 1.0/t));
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * LATENT/(r*t*t),
            latent_heat_j_kg: LATENT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 917.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}
struct HighSaturation;
impl SublimationLaw for HighSaturation {
    fn temperature_range_k(&self) -> [f64; 2] { [200.0, 900.0] }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let mut phase = Phase.properties(t)?;
        phase.saturation_pressure_pa *= 10.0;
        phase.saturation_pressure_derivative_pa_k *= 10.0;
        Ok(phase)
    }
}
fn surface_binding() -> GasSurfaceBinding {
    GasSurfaceBinding {
        species: SpeciesId::new("B").unwrap(),
        accommodation_coefficient: 0.1,
        solid_enthalpy_offset_j_kg: 0.0,
    }
}
fn chart(catalog: &ThermochemicalGasModel) -> (f64, f64) {
    let vapor = &catalog.thermodynamics().components()[1];
    let cp = vapor.evaluate(Temperature::new(400.0)).unwrap().properties().cp().value()
        / MOLAR_MASS;
    let h = catalog.species_specific_enthalpy_j_kg(
        &SpeciesId::new("B").unwrap(), 400.0,
    ).unwrap() - LATENT;
    (h, cp)
}
fn support(h: f64) -> Result<f64, SublimationError> { Ok(400.0 + h/500.0) }
fn state(chart: (f64, f64), t: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: 500.0*(350.0-400.0),
        deposit: SolidSurfaceState {
            mass_kg: 2e-6,
            specific_enthalpy_j_kg: chart.0+chart.1*(t-400.0),
            area_m2: AREA,
        },
    }
}
fn source<'m>(card: &'m MaterialCard, context: &'m QueryPoint) -> LaminarPlateFilmBinding<'m> {
    LaminarPlateFilmBinding {
        source: BinaryDiffusionMaterial {
            species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            material: card, context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        geometry: LaminarPlateGeometry { start_m: 0.02, end_m: 0.1, speed_m_s: 1.0 },
        property_temperature: ChamberFilmTemperature::MeanSurfaceAndBulk,
        maximum_vapor_mole_fraction: 0.05,
    }
}
fn diffusion() -> MaterialCard {
    material(vec![scalar(
        BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS, "temperature",
        vec![(200.0, 1e-5), (900.0, 3e-5)], [200.0, 900.0],
    )])
}
struct ConstantConductivity;
impl ConductivityLaw for ConstantConductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> { Ok(1e-3) }
}
fn controls() -> GasGrowingCoatingSolve {
    GasGrowingCoatingSolve {
        advance: GasResistiveAdvanceSolve {
            step: GasSupportedChamberSolve {
                chamber: GasChamberSolve {
                    surface: GasSurfaceSolve {
                        gas: GasSpeciesBoundarySolve {
                            gas: GasEnergySolve {
                                minimum_temperature_k: 200.0, maximum_temperature_k: 900.0,
                                energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
                            },
                            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-7,
                        },
                        maximum_gross_inventory_fraction: 0.25,
                        latent_compatibility_tolerance_j_kg: 1e-5,
                    },
                    maximum_surfaces: 8,
                },
                storage: SolverConfig {
                    relative_tolerance: 1e-12,
                    absolute_energy_tolerance_j: 1e-11, max_iterations: 128,
                },
            },
            absolute_mass_tolerance_kg: 1e-10,
            absolute_temperature_tolerance_k: 0.02,
            absolute_specific_enthalpy_tolerance_j_kg: 40.0,
            absolute_heat_tolerance_j: 1e-5,
            relative_tolerance: 1e-5,
            minimum_dt_s: 1e-9, maximum_dt_s: 0.25, maximum_trials: 2048,
        },
        maximum_deposit_volume_fraction: 0.01,
        absolute_extent_tolerance_m: 1e-10,
        absolute_conductance_tolerance_w_k: 1e-8,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn live_film_uses_actual_pressure_composition_and_declared_temperature_convention() {
    let catalog = catalog();
    let initial_gas = gas(&catalog, 400.0, 100000.0, 0.03);
    let changed_gas = gas(&catalog, 450.0, 130000.0, 0.04);
    let chart = chart(&catalog);
    let deposit = |h| Ok(400.0+(h-chart.0)/chart.1);
    let surface_binding = surface_binding();
    let row = GasResistiveSurface {
        id: 4, law: &Phase, binding: &surface_binding, state: state(chart, 350.0),
        conductance_w_k: 0.004, substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let card = material(vec![scalar(
        BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS, "pressure",
        vec![(50000.0, 2e-5), (150000.0, 2e-5/3.0)], [200.0, 900.0],
    )]);
    let context = point(250.0, 50000.0);
    for convention in [ChamberFilmTemperature::Bulk, ChamberFilmTemperature::MeanSurfaceAndBulk] {
        let binding = LaminarPlateFilmBinding {
            property_temperature: convention, ..source(&card, &context)
        };
        let prepared = SourceLaminarPlateChamberFilms::new(
            &BTreeMap::from([(4, binding.clone())]), &[row], &initial_gas, 8, || false,
        ).unwrap();
        let mut coefficients = BTreeMap::from([(4, f64::NAN)]);
        prepared.resolve(&[row], &initial_gas, &mut coefficients, &mut || false).unwrap();
        let first = coefficients[&4];
        prepared.resolve(&[row], &changed_gas, &mut coefficients, &mut || false).unwrap();
        let temperature = match convention {
            ChamberFilmTemperature::Bulk => 450.0,
            ChamberFilmTemperature::MeanSurfaceAndBulk => 400.0,
        };
        let fluid = changed_gas.model().evaluate(&point(temperature, 130000.0)).unwrap();
        let direct = evaluate_laminar_plate_gas_film(
            &fluid, &binding.source, &surface_binding.species,
            LaminarPlateFilmConditions {
                start: Length::new(0.02), end: Length::new(0.1), speed_m_s: 1.0,
                temperature: Temperature::new(temperature),
                total_pressure: Pressure::new(130000.0),
            },
            || false,
        ).unwrap();
        near(coefficients[&4], direct.mass_transfer_coefficient_m_s(), 1e-15);
        assert!((coefficients[&4]-first).abs() > 1e-5);
        let initial = prepared.selected_diffusivity(4).unwrap();
        assert_eq!(initial.state_coordinate("pressure"), Some(100000.0));
        assert_eq!(initial.state_coordinate("temperature"),
            Some(if convention == ChamberFilmTemperature::Bulk { 400.0 } else { 375.0 }));
        assert!(core::ptr::eq(prepared.binding(4).unwrap().source.material, &card));
    }
    assert_eq!(initial_gas.point().axes().get("local-lot"), Some(&9.0));
    assert_eq!(changed_gas.point().axes().get("local-lot"), Some(&9.0));
}

#[test]
fn growing_history_refreshes_original_film_and_closes_mass_and_source_energy() {
    let catalog = catalog();
    let gas = gas(&catalog, 400.0, 100000.0, 0.03);
    let original_gas = gas.clone();
    let chart = chart(&catalog);
    let deposit = |h| Ok(400.0+(h-chart.0)/chart.1);
    let binding = surface_binding();
    let row = GasResistiveSurface {
        id: 4, law: &Phase, binding: &binding, state: state(chart, 350.0),
        conductance_w_k: 0.004, substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let card = diffusion();
    let context = point(250.0, 50000.0);
    let prepared = SourceLaminarPlateChamberFilms::new(
        &BTreeMap::from([(4, source(&card, &context))]), &[row], &gas, 8, || false,
    ).unwrap();
    let mut frozen_film = BTreeMap::new();
    prepared.resolve(&[row], &gas, &mut frozen_film, &mut || false).unwrap();
    let curve = ConductivityCurve::new(vec![[200.0, 1e-3], [900.0, 1e-3]]).unwrap();
    let layers = BTreeMap::from([(4, GasIntegratedGrowingCoating {
        geometry: ConductiveCoatingGeometry::planar(917.0, AREA).unwrap(),
        conductivity: &curve, contact_resistance_m2_k_w: 2e-4,
        interface_config: LayerConductionConfig {
            max_iterations: 96, absolute_heat_tolerance_w: 1e-11,
            relative_heat_tolerance: 1e-12,
        },
    })]);
    let powers = BTreeMap::from([(4, GasResistivePower {
        substrate_power_w: -1.0, deposit_power_w: 0.05, interface: None,
    })]);
    let result = prepared.advance_growing_integrated(
        &[row], &layers, None, Some(&powers), &gas, &catalog, 1.0, controls(), || false,
    ).unwrap();
    let frozen = advance_growing_integrated_chamber(
        &[row], &layers, None, Some(&frozen_film), Some(&powers),
        &gas, &catalog, 1.0, controls(), || false,
    ).unwrap();
    let constant = ConstantConductivity;
    let authored_layers = BTreeMap::from([(4, GasGrowingCoating {
        geometry: layers[&4].geometry,
        conductivity: &constant,
        contact_resistance_m2_k_w: 2e-4,
    })]);
    let authored = prepared.advance_growing_resistive(
        &[row], &authored_layers, None, Some(&powers), &gas, &catalog, 1.0,
        controls(), || false,
    ).unwrap();
    near(
        authored.history().endpoint().surfaces()[0].state.deposit.mass_kg,
        result.history().endpoint().surfaces()[0].state.deposit.mass_kg,
        1e-15,
    );
    for geometry in [
        ConductiveCoatingGeometry::spherical(917.0, 0.01).unwrap(),
        ConductiveCoatingGeometry::cylindrical(917.0, 0.01, 0.1).unwrap(),
    ] {
        let curved = BTreeMap::from([(4, GasIntegratedGrowingCoating {
            geometry, ..layers[&4]
        })]);
        assert!(matches!(
            prepared.advance_growing_integrated(
                &[row], &curved, None, None, &gas, &catalog, 0.0,
                controls(), || false,
            ),
            Err(GasResistiveFilmError::Source(GasFilmError::Invalid(_)))
        ));
        let curved = BTreeMap::from([(4, GasGrowingCoating {
            geometry, ..authored_layers[&4]
        })]);
        assert!(matches!(
            prepared.advance_growing_resistive(
                &[row], &curved, None, None, &gas, &catalog, 0.0,
                controls(), || false,
            ),
            Err(GasResistiveFilmError::Source(GasFilmError::Invalid(_)))
        ));
    }
    let history = result.history();
    let end = history.endpoint().surfaces()[0].state;
    assert!(end.deposit.mass_kg > row.state.deposit.mass_kg);
    assert!((end.deposit.mass_kg-frozen.history().endpoint().surfaces()[0].state.deposit.mass_kg)
        .abs() > 1e-16);
    near(history.audit().thermal[0].supplied_substrate_heat_j, -1.0, 1e-12);
    near(history.audit().thermal[0].supplied_deposit_heat_j, 0.05, 1e-12);
    near(history.audit().chamber.energy_residual_j, 0.0, 1e-7);
    near(history.audit().chamber.mass_residual_kg, 0.0, 1e-14);
    near(
        end.deposit.mass_kg-row.state.deposit.mass_kg
            + result.gas().species_mass_kg(&binding.species).unwrap()
            - gas.species_mass_kg(&binding.species).unwrap(),
        0.0, 1e-14,
    );
    near(result.gas().species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(),
        gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(), 1e-14);
    let final_row = GasResistiveSurface { state: end, ..row };
    let mut final_film = BTreeMap::new();
    prepared.resolve(&[final_row], result.gas(), &mut final_film, &mut || false).unwrap();
    assert!((final_film[&4]-frozen_film[&4]).abs() > 1e-8);
    assert_eq!(gas, original_gas);
}

#[test]
fn binary_identity_carrier_and_dynamic_context_are_admitted_without_trace_seeds() {
    let catalog = catalog();
    let gas = gas(&catalog, 400.0, 100000.0, 0.03);
    let chart = chart(&catalog);
    let deposit = |h| Ok(400.0+(h-chart.0)/chart.1);
    let binding = surface_binding();
    let row = GasResistiveSurface {
        id: 4, law: &Phase, binding: &binding, state: state(chart, 350.0),
        conductance_w_k: 0.004, substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let card = diffusion();
    let context = point(250.0, 50000.0);
    let bindings = BTreeMap::from([(4, source(&card, &context))]);
    let carrier = super::tests::gas(&catalog, 400.0, 100000.0, 0.0);
    let prepared = SourceLaminarPlateChamberFilms::new(
        &bindings, &[row], &carrier, 8, || false,
    ).unwrap();
    assert_eq!(carrier.species_mass_kg(&binding.species).unwrap(), 0.0);
    let rich = super::tests::gas(&catalog, 400.0, 100000.0, 0.06);
    assert!(matches!(
        prepared.advance_resistive(
            &[row], None, &BTreeMap::new(), &rich, &catalog, 0.0,
            controls().advance, || false,
        ),
        Err(GasResistiveFilmError::Source(GasFilmError::OutsideDilute { .. }))
    ));
    let high_saturation = GasResistiveSurface { law: &HighSaturation, ..row };
    assert!(matches!(
        prepared.advance_resistive(
            &[high_saturation], None, &BTreeMap::new(), &gas, &catalog, 0.0,
            controls().advance, || false,
        ),
        Err(GasResistiveFilmError::Source(GasFilmError::OutsideDilute { .. }))
    ));
    let mut narrow = bindings.clone();
    narrow.get_mut(&4).unwrap().maximum_vapor_mole_fraction = 0.01;
    let narrow = SourceLaminarPlateChamberFilms::new(
        &narrow, &[row], &carrier, 8, || false,
    ).unwrap();
    let hot_carrier = super::tests::gas(&catalog, 900.0, 100000.0, 0.0);
    assert!(matches!(
        narrow.resolve(&[row], &hot_carrier, &mut BTreeMap::new(), &mut || false),
        Err(GasFilmError::OutsideDilute {
            saturation_pressure_fraction, kinetic_equilibrium_pressure_fraction,
            maximum_vapor_mole_fraction, ..
        }) if saturation_pressure_fraction <= maximum_vapor_mole_fraction
            && kinetic_equilibrium_pressure_fraction > maximum_vapor_mole_fraction
    ));
    let mut unrestricted = bindings.clone();
    unrestricted.get_mut(&4).unwrap().maximum_vapor_mole_fraction = 0.5;
    assert!(matches!(
        SourceLaminarPlateChamberFilms::new(&unrestricted, &[row], &gas, 8, || false),
        Err(GasFilmError::Invalid(_))
    ));
    let pure_vapor = super::tests::gas(&catalog, 400.0, 100000.0, 1.0);
    assert!(matches!(
        prepared.resolve(&[row], &pure_vapor, &mut BTreeMap::new(), &mut || false),
        Err(GasFilmError::Invalid(_))
    ));
    let third_model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species("A", 0.02897, 3.5), species("B", MOLAR_MASS, 4.0),
                species("C", 0.032, 3.5)],
            Composition::new(CompositionBasis::MoleFraction, vec![0.96, 0.03, 0.01]).unwrap(),
        ).unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.03),
            transport("C", 1.8e-5, 0.026)],
        Some(&binding.species),
    ).unwrap();
    let third = ClosedGasParcel::new(third_model.clone(), point(400.0, 100000.0), 0.01).unwrap();
    assert!(matches!(
        prepared.resolve(&[row], &third, &mut BTreeMap::new(), &mut || false),
        Err(GasFilmError::Invalid(_))
    ));
    // The original gas owner retains the source bank while omitting explicit
    // zero fractions from its active model. No trace of C may be invented.
    let zero_third = ClosedGasParcel::new(
        third_model.with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![0.97, 0.03, 0.0]).unwrap(),
        ).unwrap(),
        point(400.0, 100000.0), 0.01,
    ).unwrap();
    assert_eq!(zero_third.species_mass_kg(&SpeciesId::new("C").unwrap()).unwrap(), 0.0);
    let mut coefficients = BTreeMap::new();
    prepared.resolve(&[row], &zero_third, &mut coefficients, &mut || false).unwrap();
    assert!(coefficients[&4] > 0.0);
    let stale = context.clone().with("carrier_fraction", 0.97).unwrap();
    assert!(matches!(
        SourceLaminarPlateChamberFilms::new(
            &BTreeMap::from([(4, source(&card, &stale))]), &[row], &gas, 8, || false,
        ),
        Err(GasFilmError::Invalid(_))
    ));
    assert!(matches!(
        SourceLaminarPlateChamberFilms::new(&bindings, &[row], &gas, 0, || false),
        Err(GasFilmError::Invalid(_))
    ));
    assert!(matches!(
        prepared.resolve(&[GasResistiveSurface { id: 9, ..row }], &gas,
            &mut BTreeMap::new(), &mut || false),
        Err(GasFilmError::Invalid(_))
    ));
}

#[test]
fn selected_source_support_and_late_cancellation_refuse_without_physical_mutation() {
    let catalog = catalog();
    let gas = gas(&catalog, 350.0, 100000.0, 0.03);
    let original_gas = gas.clone();
    let chart = chart(&catalog);
    let deposit = |h| Ok(400.0+(h-chart.0)/chart.1);
    let binding = surface_binding();
    let row = GasResistiveSurface {
        id: 4, law: &Phase, binding: &binding, state: state(chart, 350.0),
        conductance_w_k: 0.004, substrate_temperature: &support,
        deposit_temperature: &deposit,
    };
    let card = material(vec![
        scalar(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS, "temperature",
            vec![(300.0, 1e-5), (400.0, 1.5e-5)], [300.0, 400.0]),
        scalar(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS, "temperature",
            vec![(450.0, 2e-5), (550.0, 2.5e-5)], [450.0, 550.0]),
    ]);
    let context = point(250.0, 50000.0);
    let bindings = BTreeMap::from([(4, source(&card, &context))]);
    let prepared = SourceLaminarPlateChamberFilms::new(
        &bindings, &[row], &gas, 8, || false,
    ).unwrap();
    let selected = prepared.selected_diffusivity(4).unwrap().clone();
    let hot_gas = super::tests::gas(&catalog, 500.0, 100000.0, 0.03);
    let hot = GasResistiveSurface { state: state(chart, 500.0), ..row };
    assert!(SourceLaminarPlateChamberFilms::new(
        &bindings, &[hot], &hot_gas, 8, || false,
    ).is_ok());
    assert!(matches!(
        prepared.advance_resistive(
            &[hot], None, &BTreeMap::new(), &hot_gas, &catalog, 0.0,
            controls().advance, || false,
        ),
        Err(GasResistiveFilmError::Source(GasFilmError::Source(_)))
    ));
    let mut polls = 0;
    let mut solve = controls().advance;
    solve.maximum_dt_s = 0.025;
    prepared.advance_resistive(
        &[row], None, &BTreeMap::new(), &gas, &catalog, 0.1, solve,
        || { polls += 1; false },
    ).unwrap();
    assert!(polls > 20);
    let at = polls-1;
    let mut replay = 0;
    let cancelled = prepared.advance_resistive(
        &[row], None, &BTreeMap::new(), &gas, &catalog, 0.1, solve,
        || { replay += 1; replay == at },
    );
    assert!(matches!(
        cancelled,
        Err(GasResistiveFilmError::Source(GasFilmError::Cancelled))
            | Err(GasResistiveFilmError::Advance(GasResistiveAdvanceError::Surface(
                GasSurfaceError::Surface(SublimationError::Cancelled)
            )))
            | Err(GasResistiveFilmError::Advance(GasResistiveAdvanceError::Thermal(
                DepositionError::Cancelled
            )))
    ));
    assert_eq!(replay, at);
    assert_eq!(prepared.selected_diffusivity(4), Some(&selected));
    assert_eq!(gas, original_gas);
}
