//! Manufactured public-API checks for source-backed finite molecular transport.
#![cfg(all(feature = "material", feature = "experimental-euler"))]

use fs_evidence::ValidityDomain;
use fs_gas::material::transport::parcel::{
    ParcelMolecularCandidate, ParcelMolecularError, ParcelMolecularLink, ParcelMolecularSolve,
    ParcelMolecularTransport,
};
use fs_gas::material::transport::{
    MaterialSoretDufourTransport, MaterialTransport, SpeciesHeatOfTransportMaterial,
    ThermalConductivityConvention,
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, DissipationError, MaxwellStefanSolveV1, PlanarTransport,
    TransportCoefficients, TransportDomain,
};
use fs_gas::multicomponent::{CaloricDomain, MixtureGas};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, network::snapshot::SnapshotNozzleLink,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError, SpeciesNetworkSolve,
    transport::{self, GasTransportError, GasTransportPaths, GasTransportSolve},
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{
    CompositionBasis as CoordinateBasis, QuantityKind, SemanticType, ValueForm,
};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{
    ChargeVector, Composition, CompositionBasis, ElementId, ElementalMatrix,
    ElementalReferenceIdV1, MassSpecificEnergyQuantityV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

const MASS: f64 = 0.01;
const DIFFUSIVITY: f64 = 1e-4;
const CONDUCTIVITY: f64 = 0.05;
type Law = MaterialSoretDufourTransport<MaterialTransport<ConstantTransport>>;
type Operator = ParcelMolecularTransport<Law>;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured finite gas transport".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn species(name: &str, mass: f64, formation: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 4.0;
    coefficients[7] = formation / R;
    let mut parameters: BTreeMap<_, _> = required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
        .map(|((key, dims), value)| (key.into(), LawParameter { value, dims }))
        .collect();
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("manufactured-finite-gas").unwrap(),
        )
        .unwrap(),
        vec![
            Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: NASA9_LAW_VERSION_V1,
                parameters,
                state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
                sources: vec![],
                provenance: provenance(),
            })
            .unwrap(),
        ],
    )
    .unwrap()
}
fn card(name: &str, claims: impl IntoIterator<Item = PropertyClaim>) -> MaterialCard {
    let mut source = ClaimSet::new();
    for claim in claims {
        source.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "gas".into(),
            process: "manufactured-finite-gas".into(),
            revision: 0,
        },
        source,
        vec![],
    )
    .unwrap()
}
fn domain(maximum: f64) -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            200.0,
            maximum,
        )
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100.0, 1e8)
}
fn scalar(name: &str, dims: fs_qty::Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: domain(1000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: provenance(),
    }
}
fn transport_card(name: &str) -> SpeciesTransportMaterial {
    SpeciesTransportMaterial {
        species: id(name),
        material: card(
            name,
            [
                scalar("dynamic_viscosity", DynViscosity::DIMS, 1e-5),
                scalar(
                    "thermal_conductivity",
                    ThermalConductivityQuantityV1::DIMS,
                    CONDUCTIVITY,
                ),
            ],
        ),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), p)
        .unwrap()
}
fn heat_source(name: &str, shift: f64) -> SpeciesHeatOfTransportMaterial {
    let dims = MassSpecificEnergyQuantityV1::DIMS;
    let value = if name == "A" {
        PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 40000.0 + shift), (650.0, 85000.0 + shift)],
            dims,
        }
    } else {
        PropertyValue::Scalar {
            value: -40000.0 + shift,
            dims,
        }
    };
    let interpolation = if name == "A" {
        InterpolationPolicy::LinearInside
    } else {
        InterpolationPolicy::ConstantWithinValidity
    };
    let key = PropertyKey::new("kinetic_heat_of_transport", dims);
    SpeciesHeatOfTransportMaterial {
        species: id(name),
        material: card(
            name,
            [PropertyClaim {
                key: key.clone(),
                value,
                validity: domain(650.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation,
                observations: vec![],
                provenance: provenance(),
            }],
        ),
        property: key,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn fixture(mass_b: f64, formation_b: f64, heat_shift: f64) -> (ThermochemicalGasModel, Operator) {
    // Reversed authored source order must be mapped to the original canonical axis.
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![species("B", mass_b, formation_b), species("A", MASS, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
    )
    .unwrap();
    let models = mixture.components().to_vec();
    let ids: Vec<_> = models.iter().map(|s| s.species().clone()).collect();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, if mass_b == MASS { 1 } else { 2 }]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 0]).unwrap(),
        models,
    )
    .unwrap();
    let bank = ThermochemicalGasModel::new(
        mixture,
        vec![transport_card("B"), transport_card("A")],
        None,
    )
    .unwrap();
    let binary = ConstantTransport::new(
        source.clone(),
        "manufactured-pair",
        1,
        TransportDomain {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1000.0,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.0),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.0),
        },
        vec![0.0, DIFFUSIVITY, DIFFUSIVITY, 0.0],
    )
    .unwrap();
    let material = MaterialTransport::new(
        bank.clone(),
        point(400.0, 1e5),
        binary,
        DynViscosity::new(0.0),
        || false,
    )
    .unwrap();
    let mole = quantity(QuantityKind::Composition(CoordinateBasis::MoleFraction));
    let reference = point(400.0, 1e5)
        .with_quantity("x_A", mole, 0.5)
        .unwrap()
        .with_quantity("x_B", mole, 0.5)
        .unwrap();
    let law = MaterialSoretDufourTransport::new(
        material,
        ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        vec![heat_source("B", heat_shift), heat_source("A", heat_shift)],
        vec![(id("B"), "x_B".into()), (id("A"), "x_A".into())],
        reference,
        || false,
    )
    .unwrap();
    let gas = MixtureGas::new(
        source,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap();
    (
        bank,
        ParcelMolecularTransport::new(PlanarTransport::new(gas, law, || false).unwrap()),
    )
}
fn parcel(bank: &ThermochemicalGasModel, x_a: f64, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        bank.with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![x_a, 1.0 - x_a]).unwrap(),
        )
        .unwrap(),
        point(t, p),
        1e-3,
    )
    .unwrap()
}
fn endpoint() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 250.0,
        maximum_temperature_k: 850.0,
        energy_tolerance_j_kg: 1e-6,
        maximum_iterations: 100,
    }
}
fn controls() -> ParcelMolecularSolve {
    ParcelMolecularSolve {
        face: MaxwellStefanSolveV1::default(),
        physical: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.4,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-8,
            },
            entropy_tolerance_j_k: 1e-10,
        },
    }
}
fn link() -> ParcelMolecularLink {
    ParcelMolecularLink {
        id: 7,
        parcels: [1, 2],
        area_m2: 0.01,
        distance_m: 0.1,
    }
}
fn solve(
    op: &mut Operator,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[ParcelMolecularLink],
    dt: f64,
) -> Result<ParcelMolecularCandidate, ParcelMolecularError> {
    let endpoints = parcels.keys().map(|&id| (id, endpoint())).collect();
    op.step_candidate(
        parcels,
        links,
        GasTransportPaths::default(),
        dt,
        &endpoints,
        controls(),
        || false,
    )
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance={tolerance:e}"
    );
}
fn mass(parcel: &ClosedGasParcel, species: &str) -> f64 {
    let thermodynamics = parcel.model().thermodynamics();
    thermodynamics
        .components()
        .iter()
        .zip(thermodynamics.mole_composition().fractions())
        .find(|(source, _)| source.species() == &id(species))
        .map_or(0.0, |(source, &x)| {
            parcel.mass_kg() * x * source.molar_mass().value()
                / parcel
                    .state()
                    .thermodynamics()
                    .properties()
                    .molar_mass()
                    .value()
        })
}
fn analytic(t: f64, p: f64, x: f64, grad_x: f64, grad_t: f64) -> (f64, f64) {
    let heat_difference = 60000.0 + 100.0 * t;
    let drive = grad_x + x * (1.0 - x) * MASS * heat_difference / (R * t) * grad_t / t;
    let j = -p / (R * t) * DIFFUSIVITY * MASS * drive;
    (j, -CONDUCTIVITY * grad_t + heat_difference * j)
}

#[test]
fn material_soret_and_fourier_dufour_drive_actual_finite_temperature_and_species() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.5, 390.0, 1e5)),
        (2, parcel(&bank, 0.5, 410.0, 1e5)),
    ]);
    let next = solve(&mut op, &original, &[link()], 0.1).unwrap();
    let (j, q) = analytic(400.0, 1e5, 0.5, 0.0, 200.0);
    assert!(j < 0.0 && q < 0.0);
    let report = &next.step().links[0];
    close(report.mass_left_to_right_kg[0], j * 0.001, 1e-18);
    close(report.face.reduced_heat_flux_w_m2, q, 1e-10);
    assert_eq!(next.step().species, vec![id("A"), id("B")]);
    for (id, sign) in [(1, -1.0), (2, 1.0)] {
        let old = &original[&id];
        let new = &next.parcels()[&id];
        close(mass(new, "A") - mass(old, "A"), sign * j * 0.001, 1e-18);
        close(new.mass_kg(), old.mass_kg(), 1e-18);
        let temperature =
            old.state().state().temperature + sign * q * 0.001 / (old.mass_kg() * 3.0 * R / MASS);
        close(new.state().state().temperature, temperature, 2e-10);
    }
    assert!(next.step().transport.entropy_change_j_k > 0.0);
    let law = op.face().law();
    let last_temperature = next.parcels()[&2].state().state().temperature;
    close(
        law.last_resolved().temperature().value(),
        last_temperature,
        1e-12,
    );
    close(
        law.last_resolved().heats()[0].value(),
        20000.0 + 100.0 * last_temperature,
        1e-8,
    );
    assert_eq!(
        law.last_resolved().components()[0].card_identity(),
        law.sources()[0].material.content_hash()
    );
    let next = next.into_parts().0;
    solve(&mut op, &next, &[link()], 0.1).unwrap();
    assert_ne!(
        op.face().law().last_resolved().temperature().value(),
        last_temperature
    );
}

#[test]
fn isothermal_composition_diffusion_produces_dufour_heat_without_double_enthalpy() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.8, 400.0, 1e5)),
        (2, parcel(&bank, 0.2, 400.0, 1e5)),
    ]);
    let next = solve(&mut op, &original, &[link()], 0.1).unwrap();
    let (j, q) = analytic(400.0, 1e5, 0.5, -6.0, 0.0);
    let report = &next.step().links[0];
    close(report.mass_left_to_right_kg[0], j * 0.001, 1e-18);
    close(report.heat_correction_left_to_right_j, q * 0.001, 1e-12);
    for (id, sign) in [(1, -1.0), (2, 1.0)] {
        close(
            next.parcels()[&id].state().state().temperature,
            400.0 + sign * q * 0.001 / (original[&id].mass_kg() * 3.0 * R / MASS),
            2e-10,
        );
    }
    let total_before: f64 = original
        .values()
        .map(|p| p.mass_kg() * p.state().specific_internal_energy_j_kg())
        .sum();
    let total_after: f64 = next
        .parcels()
        .values()
        .map(|p| p.mass_kg() * p.state().specific_internal_energy_j_kg())
        .sum();
    close(total_after, total_before, 1e-9);
}

#[test]
fn formation_and_common_kinetic_heat_datums_do_not_change_the_endpoint() {
    let run = |formation, shift| {
        let (bank, mut op) = fixture(MASS, formation, shift);
        let old = BTreeMap::from([
            (1, parcel(&bank, 0.65, 390.0, 1e5)),
            (2, parcel(&bank, 0.35, 410.0, 1e5)),
        ]);
        solve(&mut op, &old, &[link()], 0.1).unwrap()
    };
    let baseline = run(0.0, 0.0);
    let shifted = run(1e14, 2.0_f64.powi(40));
    for id in [1, 2] {
        close(
            shifted.parcels()[&id].state().state().temperature,
            baseline.parcels()[&id].state().state().temperature,
            1e-9,
        );
        close(
            mass(&shifted.parcels()[&id], "A"),
            mass(&baseline.parcels()[&id], "A"),
            1e-18,
        );
    }
    close(
        shifted.step().links[0].heat_correction_left_to_right_j,
        baseline.step().links[0].heat_correction_left_to_right_j,
        1e-12,
    );
    assert!(
        shifted
            .step()
            .transport
            .exchange
            .total_energy_residual_j
            .abs()
            <= controls().physical.exchange.energy_tolerance_j
    );
}

#[test]
fn unequal_masses_use_the_mass_frame_and_reversing_a_face_preserves_the_update() {
    let (bank, mut op) = fixture(2.0 * MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.7, 390.0, 1e5)),
        (2, parcel(&bank, 0.3, 410.0, 1e5)),
    ]);
    let next = solve(&mut op, &original, &[link()], 0.1).unwrap();
    let mass = &next.step().links[0].mass_left_to_right_kg;
    close(mass.iter().sum(), 0.0, 1e-18);
    assert!((mass[0] / MASS + mass[1] / (2.0 * MASS)).abs() > 1e-10);
    let reversed = solve(
        &mut op,
        &original,
        &[ParcelMolecularLink {
            parcels: [2, 1],
            ..link()
        }],
        0.1,
    )
    .unwrap();
    for id in [1, 2] {
        close(
            next.parcels()[&id].mass_kg(),
            original[&id].mass_kg(),
            1e-18,
        );
        close(
            next.parcels()[&id].state().state().temperature,
            reversed.parcels()[&id].state().state().temperature,
            1e-10,
        );
    }
}

#[test]
fn native_nozzle_and_thermal_faces_share_one_original_gross_withdrawal_budget() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 0.8, 400.0, 1e5)),
        (2, parcel(&bank, 0.2, 400.0, 1e5)),
        (3, parcel(&bank, 0.5, 400.0, 8e4)),
    ]);
    let endpoints: BTreeMap<_, _> = original.keys().map(|&id| (id, endpoint())).collect();
    let nozzles = [SnapshotNozzleLink {
        id: 91,
        parcels: [1, 3],
        nozzle: GasNozzle {
            area_m2: 1e-8,
            discharge_coefficient: 1.0,
        },
        nozzle_solve: endpoint(),
    }];
    let paths = GasTransportPaths {
        nozzles: &nozzles,
        ..GasTransportPaths::default()
    };
    let native = transport::step_candidate_with_paths(
        &original,
        paths,
        0.1,
        &endpoints,
        controls().physical,
        || false,
    )
    .unwrap();
    let molecular = solve(&mut op, &original, &[link()], 0.1).unwrap();
    let joint = op
        .step_candidate(
            &original,
            &[link()],
            paths,
            0.1,
            &endpoints,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(joint.step().transport.nozzles, native.step().nozzles);
    assert_eq!(joint.step().links, molecular.step().links);
    let native_mass: f64 = native
        .step()
        .exchange
        .fluxes
        .iter()
        .filter(|row| row.transfer.donor == 1 && row.transfer.species == id("A"))
        .map(|row| row.transfer.mass_kg)
        .sum();
    let molecular_mass = molecular.step().links[0].mass_left_to_right_kg[0];
    assert!(native_mass > 0.0 && molecular_mass > 0.0);
    let mut tight = controls();
    tight.physical.exchange.maximum_species_fraction = (native_mass.max(molecular_mass)
        + 0.5 * native_mass.min(molecular_mass))
        / mass(&original[&1], "A");
    let result = op.step_candidate(&original, &[link()], paths, 0.1, &endpoints, tight, || {
        false
    });
    assert!(matches!(
        result,
        Err(ParcelMolecularError::Transport(
            GasTransportError::Exchange(SpeciesNetworkError::Withdrawal { id: 1, .. })
        ))
    ));
}

#[test]
fn source_refusal_late_cancellation_and_missing_species_are_atomic() {
    let (bank, mut op) = fixture(MASS, 0.0, 0.0);
    let original = BTreeMap::from([
        (1, parcel(&bank, 1.0, 390.0, 1e5)),
        (2, parcel(&bank, 0.0, 410.0, 1e5)),
    ]);
    let before = original.clone();
    let next = solve(&mut op, &original, &[link()], 0.001).unwrap();
    assert!(mass(&next.parcels()[&1], "B") > 0.0 && mass(&next.parcels()[&2], "A") > 0.0);
    assert_eq!(original, before);
    // Its midpoint is supported (550 K), but the actual right gas is not.
    let hot = BTreeMap::from([
        (1, parcel(&bank, 0.5, 400.0, 1e5)),
        (2, parcel(&bank, 0.5, 700.0, 1e5)),
    ]);
    let saved_hot = hot.clone();
    assert!(matches!(
        solve(&mut op, &hot, &[link()], 0.1),
        Err(ParcelMolecularError::Face(DissipationError::Material(_)))
    ));
    assert_eq!(hot, saved_hot);
    close(
        op.face().law().last_resolved().temperature().value(),
        400.0,
        1e-12,
    );
    // Every initial state and the original face are inside support. Dufour
    // heats the right endpoint beyond 650 K, which must refuse publication.
    let edge = BTreeMap::from([
        (1, parcel(&bank, 0.8, 649.99, 1e5)),
        (2, parcel(&bank, 0.2, 649.99, 1e5)),
    ]);
    let saved_edge = edge.clone();
    assert!(matches!(
        solve(&mut op, &edge, &[link()], 0.1),
        Err(ParcelMolecularError::Face(DissipationError::Material(_)))
    ));
    assert_eq!(edge, saved_edge);
    assert!(op.face().law().last_resolved().temperature().value() < 649.99);
    let endpoints = original.keys().map(|&id| (id, endpoint())).collect();
    let mut calls = 0;
    let cancelled = op.step_candidate(
        &original,
        &[link()],
        GasTransportPaths::default(),
        0.001,
        &endpoints,
        controls(),
        || {
            calls += 1;
            calls > 20
        },
    );
    assert!(matches!(cancelled, Err(ParcelMolecularError::Cancelled)));
    assert!(calls > 20);
    assert_eq!(original, before);
    let (other, _) = fixture(MASS, 1234.0, 0.0);
    let mismatched = BTreeMap::from([
        (1, parcel(&other, 0.5, 400.0, 1e5)),
        (2, original[&2].clone()),
    ]);
    assert!(matches!(
        solve(&mut op, &mismatched, &[link()], 0.001),
        Err(ParcelMolecularError::SourceMismatch(_))
    ));
}
