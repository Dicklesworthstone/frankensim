//! Actual NASA-9/material integration; manufactured fixtures, no empirical claims.
#[path = "sealed/advance.rs"]
mod advance;
#[path = "sealed/contact_line.rs"]
mod contact_line;

use super::*;
use fs_flux::capillary::material::sealed::{
    SealedCapillaryError, SealedCapillaryGeometry, SealedCapillaryInput, SealedCapillarySolve,
};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured sealed capillary gas".into(),
        license: "test-only".into(),
        artifact: None,
    }
}
fn gas(volume: f64, cp_over_r: f64) -> ClosedGasParcel {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_over_r;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("A").unwrap(),
        MolarMass::new(0.028),
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
                provenance: provenance(),
            })
            .unwrap(),
        ],
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 1.8e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.026),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        typed(QuantityKind::AbsoluteTemperature),
                        200.0,
                        2000.0,
                    )
                    .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e8),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    let transport = SpeciesTransportMaterial {
        species: SpeciesId::new("A").unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: "A".into(),
                phase: "gas".into(),
                process: "manufactured".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport],
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            typed(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), 100000.0)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(model, point, volume).unwrap()
}
fn geometry() -> SealedCapillaryGeometry {
    SealedCapillaryGeometry {
        radius_m: 1e-4,
        tube_length_m: 0.1,
        wet_length_m: 0.02,
        dead_volume_m3: 0.0,
    }
}
fn volume(g: SealedCapillaryGeometry) -> f64 {
    g.dead_volume_m3 + PI * g.radius_m.powi(2) * (g.tube_length_m - g.wet_length_m)
}
fn controls() -> SealedCapillarySolve {
    SealedCapillarySolve {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1000.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 80,
        },
        maximum_iterations: 80,
        pressure_absolute_tolerance_pa: 1e-6,
        length_absolute_tolerance_m: 1e-12,
        energy_absolute_tolerance_j: 1e-13,
        energy_relative_tolerance: 1e-9,
    }
}
fn boundary() -> SealedCapillaryInput {
    SealedCapillaryInput {
        inlet_pressure_pa: 100000.0,
        duration_s: 1.0,
        gas_heat_power_w: 0.0,
        opposing_gravity_m_s2: 0.0,
    }
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + 2e-9 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn actual_gas_midpoint_drives_compression_drainage_heat_and_gravity_work() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = geometry();
    let original_gas = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(original_gas.clone(), g)
        .unwrap();
    for displacement in [0.001, -0.001] {
        for heat in [0.0, 1e-5] {
            for gravity in [0.0, 9.81, -9.81] {
                let x1 = g.wet_length_m + displacement;
                let xm = g.wet_length_m + 0.5 * displacement;
                let transferred = PI * g.radius_m.powi(2) * displacement;
                let v0 = original_gas.volume_m3();
                let v1 = v0 - transferred;
                // Independent closed-form perfect-gas MIDPOINT law, not a finite isentrope.
                let pm = (2.0 * 100000.0 * v0 + 0.4 * heat) / (1.4 * v1 + 0.6 * v0);
                let t1 = 2.0 * pm * (0.5 * (v0 + v1)) / (100000.0 * v0) * 300.0 - 300.0;
                let pin = pm - 300.0
                    + 8.0 * xm * displacement / g.radius_m.powi(2)
                    + 1000.0 * gravity * xm;
                let candidate = owner
                    .step_candidate(
                        SealedCapillaryInput {
                            inlet_pressure_pa: pin,
                            gas_heat_power_w: heat,
                            opposing_gravity_m_s2: gravity,
                            ..boundary()
                        },
                        controls(),
                        || false,
                    )
                    .unwrap();
                let step = candidate.step();
                near(step.wet_length_m, x1, 2e-12);
                near(step.midpoint_gas_pressure_pa, pm, 1e-6);
                near(
                    candidate.endpoint().gas().state().state().temperature,
                    t1,
                    1e-8,
                );
                near(
                    candidate.endpoint().gas().state().state().pressure,
                    100000.0 * (t1 / 300.0) * (v0 / v1),
                    1e-5,
                );
                near(step.transferred_liquid_mass_kg, 1000.0 * transferred, 1e-16);
                near(
                    step.liquid_viscous_dissipation_j,
                    8.0 * PI * xm * displacement.powi(2),
                    1e-15,
                );
                near(
                    step.gas.stored_energy_change_j,
                    heat + pm * transferred,
                    1e-13,
                );
                near(
                    step.gravitational_energy_change_j,
                    1000.0 * gravity * xm * transferred,
                    1e-15,
                );
                assert!(step.energy_balance_residual_j.abs() < 1e-13);
                assert_eq!(step.gas.supplied_heat_j, heat);
                assert_eq!(candidate.endpoint().gas().mass_kg(), original_gas.mass_kg());
                assert_eq!(candidate.endpoint().gas().model(), original_gas.model());
                assert_eq!(
                    candidate
                        .endpoint()
                        .gas()
                        .point()
                        .axes()
                        .get("specimen-lot"),
                    Some(&42.0)
                );
                assert!(std::ptr::eq(candidate.endpoint().tube().liquid(), &liquid));
                assert!(std::ptr::eq(
                    candidate.endpoint().tube().wetting(),
                    &wetting
                ));
                assert_eq!(owner.gas(), &original_gas);
                assert_eq!(owner.geometry(), g);
            }
        }
    }
}

#[test]
fn gas_caloric_source_changes_coupled_motion_without_replacing_pressure() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    let mut lengths = Vec::new();
    for cp in [3.5, 2.5] {
        let owner = tube.seal(gas(volume(geometry()), cp), geometry()).unwrap();
        let candidate = owner
            .step_candidate(
                SealedCapillaryInput {
                    inlet_pressure_pa: 125000.0,
                    gas_heat_power_w: 1e-5,
                    ..boundary()
                },
                controls(),
                || false,
            )
            .unwrap();
        lengths.push(candidate.step().wet_length_m);
        assert!(candidate.endpoint().gas().state().state().temperature > 300.0);
    }
    assert!(lengths[0] > lengths[1] + 1e-8, "{lengths:?}");
}

#[test]
fn cancellation_domains_equilibrium_overshoot_and_geometry_refuse_atomically() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    let original = gas(volume(geometry()), 3.5);
    let owner = tube.seal(original.clone(), geometry()).unwrap();
    assert!(matches!(
        owner.step_candidate(boundary(), controls(), || true),
        Err(SealedCapillaryError::Cancelled)
    ));
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                duration_s: 1e5,
                ..boundary()
            },
            controls(),
            || false
        ),
        Err(SealedCapillaryError::StepTooLarge(_))
    ));
    let narrow = SealedCapillarySolve {
        gas: GasEnergySolve {
            maximum_temperature_k: 300.0001,
            ..controls().gas
        },
        ..controls()
    };
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                gas_heat_power_w: 1e-4,
                ..boundary()
            },
            narrow,
            || false
        ),
        Err(SealedCapillaryError::Gas(_))
    ));
    assert!(matches!(
        tube.seal(
            original.clone(),
            SealedCapillaryGeometry {
                dead_volume_m3: original.volume_m3(),
                ..geometry()
            }
        ),
        Err(SealedCapillaryError::Invalid(_))
    ));
    assert_eq!(owner.gas(), &original);
    assert_eq!(owner.geometry(), geometry());
}

#[test]
fn finite_dead_volume_allows_full_tube_but_crossing_its_end_never_publishes() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = SealedCapillaryGeometry {
        wet_length_m: 0.1,
        dead_volume_m3: 1e-8,
        ..geometry()
    };
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(gas(volume(g), 3.5), g)
        .unwrap();
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                inlet_pressure_pa: 110000.0,
                ..boundary()
            },
            controls(),
            || false
        ),
        Err(SealedCapillaryError::StepTooLarge(_))
    ));
    assert_eq!(owner.geometry(), g);
}

#[path = "sealed/inertia.rs"]
mod inertia;
