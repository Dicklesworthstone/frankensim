//! Original material/NASA-9 fixtures and independent liquid-column identities.
use super::super::gas::TrappedGasPocket;
use super::*;
use fs_duct::dg::DgMesh;
use fs_duct::{Duct, Segment};
use fs_evidence::ValidityDomain;
use fs_flux::cavitation::BubbleGas;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};

const REFERENCE: f64 = 100_000.0;
const VAPOR: f64 = 2000.0;
const PG: f64 = REFERENCE - VAPOR;
const INITIAL_T: f64 = 400.0;
const VOLUME: f64 = 1e-6;
const RADIUS: f64 = 0.01;
const MOLAR_MASS: f64 = 0.028;
fn area() -> f64 {
    core::f64::consts::PI * RADIUS * RADIUS
}
fn gas_constant() -> f64 {
    fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance(label: &str) -> Provenance {
    Provenance {
        source: label.into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn parcel(slope: f64, formation: f64, transport_maximum: f64) -> ClosedGasParcel {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[3] = slope;
    coefficients[7] = formation;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let species = SpeciesId::new("caloric-pocket-fixture").unwrap();
    let caloric = Nasa9StandardStateModelV1::new(
        species.clone(),
        MolarMass::new(MOLAR_MASS),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic").unwrap(),
        )
        .unwrap(),
        vec![
            Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: NASA9_LAW_VERSION_V1,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
                sources: vec![],
                provenance: provenance("synthetic variable-Cv identity"),
            })
            .unwrap(),
        ],
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (2000.0, 2.0 * value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        transport_maximum,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance("synthetic supported transport"),
            })
            .unwrap();
    }
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![caloric],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![SpeciesTransportMaterial {
            species,
            material: MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "synthetic gas".into(),
                    phase: "gas".into(),
                    process: "test".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            INITIAL_T,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), PG)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(model, point, VOLUME).unwrap()
}
fn material() -> CavitationState {
    CavitationState::try_new(293.15, 1000.0, 0.0, 0.0, VAPOR, 1500.0).unwrap()
}
fn native(outlet: DgBoundary) -> DgDuct<DgFluidMaterial<CavitationState>> {
    DgDuct::from_material(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: 1.0,
            }],
        },
        &material(),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(8)),
        outlet,
    )
    .unwrap()
}
fn limits() -> DuctCavityLimits {
    DuctCavityLimits {
        reference_pressure_pa: REFERENCE,
        maximum_cavity_volume_m3: 3e-6,
        maximum_liquid_strain: 0.001,
        maximum_interface_mach: 0.01,
        maximum_vapor_density_ratio: 0.01,
    }
}
fn numerics() -> DuctCavityNumerics {
    DuctCavityNumerics {
        relative_tolerance: 1e-10,
        absolute_energy_tolerance_j: 5e-11,
        absolute_pressure_tolerance_pa: 2e-3,
        absolute_volume_tolerance_m3: 1e-14,
        maximum_step_s: 1e-5,
        event_time_tolerance_s: 1e-13,
        event_scan_substeps: 8,
        maximum_attempts: 30_000,
        maximum_events: 32,
    }
}
fn specification(port: u32, volume: f64) -> DuctPocket {
    DuctPocket {
        port: PortId(port),
        initial_volume_m3: volume,
        maximum_volume_m3: 2e-6,
    }
}
fn caloric(gas: ClosedGasParcel) -> CaloricGasPocket {
    CaloricGasPocket {
        port: PortId(8),
        gas,
        minimum_volume_m3: 2e-7,
        solve: GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 1600.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
    }
}
fn model(gas: ClosedGasParcel, flow: f64) -> CavitatingDuctPockets<CavitationState> {
    let duct = native(DgBoundary::IdealOpen);
    let initial = duct.project_initial_state(|_| 0.0, |_| flow);
    CavitatingDuctPockets::new_with_caloric_gas(
        duct,
        initial,
        0.0,
        &[specification(8, VOLUME)],
        &[],
        &[caloric(gas)],
        limits(),
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e} ({tolerance:e})");
}

#[test]
fn original_caloric_zero_heat_recovers_adiabatic_polytrope_and_source_inventory() {
    let gas = parcel(0.0, 0.0, 2000.0);
    let original = gas.clone();
    let gamma = gas.state().state().gamma;
    let q0 = area() * 0.1;
    let mut thermal = model(gas, q0);
    let duct = native(DgBoundary::IdealOpen);
    let initial = duct.project_initial_state(|_| 0.0, |_| q0);
    let mut polytropic = CavitatingDuctPockets::new_with_trapped_gas(
        duct,
        initial,
        0.0,
        &[specification(8, VOLUME)],
        &[TrappedGasPocket {
            port: PortId(8),
            gas: BubbleGas::from_reference_volume_m3(VOLUME, PG, gamma).unwrap(),
            minimum_volume_m3: 2e-7,
        }],
        limits(),
    )
    .unwrap();
    let report = thermal.advance(5e-4, &numerics(), &mut || false).unwrap();
    let reference = polytropic
        .advance(5e-4, &numerics(), &mut || false)
        .unwrap();
    let endpoint = thermal.caloric_gas(PortId(8)).unwrap();
    let ratio = VOLUME / endpoint.volume_m3();
    near(
        endpoint.state().state().temperature,
        INITIAL_T * ratio.powf(gamma - 1.0),
        2e-6,
    );
    near(
        endpoint.state().state().pressure,
        PG * ratio.powf(gamma),
        2e-3,
    );
    near(
        report.pockets[0].volume_m3,
        reference.pockets[0].volume_m3,
        2e-13,
    );
    assert_eq!(endpoint.mass_kg(), original.mass_kg());
    assert_eq!(endpoint.model(), original.model());
    assert_eq!(endpoint.point().axes().get("specimen-lot"), Some(&42.0));
    assert_eq!(report.supplied_heat_j, 0.0);
    assert!(report.caloric_gas[0].internal_energy_change_j < 0.0);
    assert!(report.energy_balance_residual_j.abs() < 1e-10);
}

#[test]
fn variable_cv_heat_drives_exact_isobaric_column_and_preserves_large_formation_reference() {
    let slope = 0.002;
    let duration = 1e-3;
    for (sign, formation) in [(1.0, 0.0), (-1.0, 1e12)] {
        let gas = parcel(slope, formation, 2000.0);
        let mass = gas.mass_kg();
        let q0 = sign * area() * 0.1;
        let mut coupled = model(gas, q0);
        let original = coupled.state().duct().clone();
        // Independent exact solution: p_g=constant, V'=Q0, T=T0 V/V0.
        // Ideal enthalpy delivery supplies both actual Cv(T) storage and p dV.
        let report = coupled
            .advance_with_gas_heat(
                duration,
                &numerics(),
                &|_, sample| Ok(PG * q0 * (1.0 + sample.cv_j_kg_k / gas_constant())),
                &mut || false,
            )
            .unwrap();
        let gas = coupled.caloric_gas(PortId(8)).unwrap();
        let v1 = VOLUME + q0 * duration;
        let t1 = INITIAL_T * v1 / VOLUME;
        let expected_u = mass
            * gas_constant()
            * (2.5 * (t1 - INITIAL_T) + 0.5 * slope * (t1 * t1 - INITIAL_T * INITIAL_T));
        let expected_heat = expected_u + PG * (v1 - VOLUME);
        near(gas.volume_m3(), v1, 2e-13);
        near(gas.state().state().temperature, t1, 2e-6);
        near(gas.state().state().pressure, PG, 2e-3);
        near(
            report.caloric_gas[0].internal_energy_change_j,
            expected_u,
            1e-10,
        );
        near(
            report.caloric_gas[0].work_by_gas_j,
            PG * (v1 - VOLUME),
            1e-10,
        );
        near(report.supplied_heat_j, expected_heat, 1e-10);
        let view = coupled.duct().continuous_operator();
        let before = view.energy_coordinates(&original).unwrap();
        let after = view.energy_coordinates(coupled.state().duct()).unwrap();
        for (&a, &b) in after.iter().zip(before) {
            near(a, b, 5e-9);
        }
        assert!(report.energy_balance_residual_j.abs() < 1e-10);
        assert_eq!(gas.mass_kg(), mass);
    }
}

#[test]
fn caloric_gas_and_neighboring_pure_vapor_share_the_same_collapse_event() {
    let duct = native(DgBoundary::Port(PortId(17)));
    let initial = duct.project_initial_state(|x| -PG * x, |_| area() * 0.01);
    let mut coupled = CavitatingDuctPockets::new_with_caloric_gas(
        duct,
        initial,
        0.0,
        &[specification(17, 1e-10), specification(8, VOLUME)],
        &[],
        &[caloric(parcel(0.002, 0.0, 2000.0))],
        limits(),
    )
    .unwrap();
    let duration = 5e-5;
    let report = coupled
        .advance_with_gas_heat(duration, &numerics(), &|_, _| Ok(20.0), &mut || false)
        .unwrap();
    assert_eq!(report.pockets[0].port, PortId(8));
    assert_eq!(report.pockets[0].collapse_events, 0);
    assert_eq!(report.pockets[1].collapse_events, 1);
    assert_eq!(report.pockets[1].mode, DuctCavityMode::LiquidFilled);
    near(report.supplied_heat_j, 20.0 * duration, 1e-11);
    assert!(report.energy_balance_residual_j.abs() < 1e-10);
}

#[test]
fn late_cancel_and_original_transport_refusal_preserve_liquid_gas_and_clock() {
    let q0 = area() * 0.1;
    let heat =
        |_: f64, s: &CaloricGasObservation| Ok(PG * q0 * (1.0 + s.cv_j_kg_k / gas_constant()));
    let mut probe = model(parcel(0.002, 0.0, 2000.0), q0);
    let mut polls = 0;
    probe
        .advance_with_gas_heat(1e-3, &numerics(), &heat, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut actual = model(parcel(0.002, 0.0, 2000.0), q0);
    let before = actual.state().clone();
    let before_gas = actual.caloric_gas(PortId(8)).unwrap().clone();
    let mut seen = 0;
    assert_eq!(
        actual
            .advance_with_gas_heat(1e-3, &numerics(), &heat, &mut || {
                seen += 1;
                seen == polls
            })
            .unwrap_err(),
        DuctCavityError::Cancelled
    );
    assert_eq!(actual.state(), &before);
    assert_eq!(actual.caloric_gas(PortId(8)).unwrap(), &before_gas);

    // Caloric source remains valid, but the final transport source is not.
    let mut unsupported = model(parcel(0.002, 0.0, 401.0), q0);
    let before = unsupported.state().clone();
    let before_gas = unsupported.caloric_gas(PortId(8)).unwrap().clone();
    assert!(matches!(
        unsupported.advance_with_gas_heat(1e-3, &numerics(), &heat, &mut || false),
        Err(DuctCavityError::Gas(_))
    ));
    assert_eq!(unsupported.state(), &before);
    assert_eq!(unsupported.caloric_gas(PortId(8)).unwrap(), &before_gas);
}
