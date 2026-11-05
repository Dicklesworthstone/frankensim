//! G1/G3/G4: physical cross diffusion in the original complete live assembly.
//! Synthetic source cards; no measured-material or engine validation claim.
#[path = "molecular_coupled.rs"]
mod coupled_tests;
#[path = "molecular_endpoint.rs"]
mod endpoint_tests;
#[path = "molecular_gas_law.rs"]
mod gas_law_tests;

use super::*;
use fs_gas::material::transport::parcel::{
    ParcelMolecularError, ParcelMolecularLink, ParcelMolecularSolve, ParcelMolecularTransport,
};
use fs_gas::material::transport::{
    MaterialSoretDufourTransport, SpeciesHeatOfTransportMaterial, ThermalConductivityConvention,
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, MaxwellStefanSolveV1, PlanarTransport, TransportCoefficients,
    TransportDomain,
};
use fs_gas::multicomponent::{CaloricDomain, MixtureGas};
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::network::snapshot::SnapshotNozzleLink,
    species::network::{
        SpeciesNetworkSolve,
        transport::{GasTransportPaths, GasTransportSolve},
    },
};
use fs_qty::DynViscosity;
use fs_qty::semantic::CompositionBasis as CoordinateBasis;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix, MassSpecificEnergyQuantityV1};

const DT: f64 = 0.00005;
const MASS: f64 = 0.028;
const DIFFUSIVITY: f64 = 2e-5;
const CONDUCTIVITY: f64 = 0.03;
type MolecularLaw = MaterialSoretDufourTransport<ConstantTransport>;
type Molecular = ParcelMolecularTransport<MolecularLaw>;

fn species_id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn gas_point(t: f64, p: f64) -> QueryPoint {
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
fn heat_source(name: &str, scale: f64) -> SpeciesHeatOfTransportMaterial {
    let dims = MassSpecificEnergyQuantityV1::DIMS;
    let property = PropertyKey::new("kinetic_heat_of_transport", dims);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: property.clone(),
            value: if name == "A" {
                PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, 40000.0 * scale), (1800.0, 200000.0 * scale)],
                    dims,
                }
            } else {
                PropertyValue::Scalar {
                    value: -40000.0 * scale,
                    dims,
                }
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    200.0,
                    1800.0,
                )
                .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if name == "A" {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            observations: vec![],
            provenance: Provenance {
                source: "synthetic molecular chamber cross coefficient".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    SpeciesHeatOfTransportMaterial {
        species: species_id(name),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "molecular assembly fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        property,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn molecular_fixture(
    mass_b: f64,
    cross_scale: f64,
    formation_b: f64,
) -> (ThermochemicalGasModel, Molecular) {
    let (a, ta) = gas_fixture::component("A", MASS, 0.0);
    let (b, tb) = gas_fixture::component("B", mass_b, formation_b);
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![b, a],
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
    let bank = ThermochemicalGasModel::new(mixture, vec![tb, ta], None).unwrap();
    let binary = ConstantTransport::new(
        source.clone(),
        "synthetic-assembly-pair",
        1,
        TransportDomain {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            minimum_pressure_pa: 1000.0,
            maximum_pressure_pa: 1e7,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.0),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(CONDUCTIVITY),
        },
        vec![0.0, DIFFUSIVITY, DIFFUSIVITY, 0.0],
    )
    .unwrap();
    let mole = quantity(QuantityKind::Composition(CoordinateBasis::MoleFraction));
    let reference = gas_point(400.0, 1e5)
        .with_quantity("x_A", mole, 0.5)
        .unwrap()
        .with_quantity("x_B", mole, 0.5)
        .unwrap();
    let law = MaterialSoretDufourTransport::new(
        binary,
        ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        vec![heat_source("B", cross_scale), heat_source("A", cross_scale)],
        vec![
            (species_id("B"), "x_B".into()),
            (species_id("A"), "x_A".into()),
        ],
        reference,
        || false,
    )
    .unwrap();
    let gas = MixtureGas::new(
        source,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1800.0)),
        || false,
    )
    .unwrap();
    (
        bank,
        ParcelMolecularTransport::new(PlanarTransport::new(gas, law, || false).unwrap()),
    )
}
fn molecular_gases(
    bank: &ThermochemicalGasModel,
    states: [(f64, f64, f64); 2],
    fixed: bool,
) -> BTreeMap<u64, GasPiston> {
    states
        .into_iter()
        .enumerate()
        .map(|(i, (x, t, p))| {
            let model = bank
                .with_composition(
                    Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
                )
                .unwrap();
            let gas = GasPiston::new(
                ClosedGasParcel::new(model, gas_point(t, p), 0.001).unwrap(),
                1.0,
                0.01,
                0.0,
            )
            .unwrap();
            (
                10 + 10 * i as u64,
                if fixed {
                    gas.with_fixed_support().unwrap()
                } else {
                    gas
                },
            )
        })
        .collect()
}
fn molecular_controls() -> ParcelMolecularSolve {
    ParcelMolecularSolve {
        face: MaxwellStefanSolveV1::default(),
        physical: GasTransportSolve {
            exchange: SpeciesNetworkSolve {
                maximum_species_fraction: 0.05,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-7,
            },
            entropy_tolerance_j_k: 1e-10,
        },
    }
}
fn molecular_link() -> ParcelMolecularLink {
    ParcelMolecularLink {
        id: 27,
        parcels: [10, 20],
        area_m2: 0.1,
        distance_m: 0.0001,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn species_mass(gas: &ClosedGasParcel, species: &str) -> f64 {
    let thermo = gas.model().thermodynamics();
    let mean = gas
        .state()
        .thermodynamics()
        .properties()
        .molar_mass()
        .value();
    thermo
        .components()
        .iter()
        .zip(thermo.mole_composition().fractions())
        .find(|(s, _)| s.species() == &species_id(species))
        .map_or(0.0, |(s, x)| {
            gas.mass_kg() * x * s.molar_mass().value() / mean
        })
}

#[test]
fn g1_soret_and_reciprocal_heat_reach_the_complete_assembly_and_pure_heat_has_stream_ports() {
    for scale in [1.0, 0.0] {
        fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
            with_cx(|cx, _| {
                let (bank, mut molecular) = molecular_fixture(MASS, scale, 0.0);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let state = MaterialExchangeState::from_parts(
                    thermal,
                    circuit,
                    molecular_gases(&bank, [(0.5, 390.0, 1e5), (0.5, 410.0, 1e5)], true),
                );
                let drives = drives(state.chambers());
                let (next, report) = state
                    .molecular_transport_step_candidate(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &[],
                        &mut molecular,
                        &[molecular_link()],
                        GasTransportPaths::default(),
                        molecular_controls(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                // Equal-mass binary law, independent of the Maxwell--Stefan solve.
                let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
                let heat_difference = 100000.0 * scale;
                let gradient_t = 20.0 / molecular_link().distance_m;
                let drive = 0.25 * MASS * heat_difference / (r * 400.0) * gradient_t / 400.0;
                let j = -1e5 / (r * 400.0) * DIFFUSIVITY * MASS * drive;
                let q = -CONDUCTIVITY * gradient_t + heat_difference * j;
                let area_time = molecular_link().area_m2 * DT;
                let link = &report.flow.links[0];
                close(link.mass_left_to_right_kg[0], j * area_time, 1e-16);
                close(link.mass_left_to_right_kg[1], -j * area_time, 1e-16);
                close(link.face.reduced_heat_flux_w_m2, q, 1e-7);
                for (id, sign) in [(10, -1.0), (20, 1.0)] {
                    let old = state.chambers()[&id].gas();
                    let new = next.chambers()[&id].gas();
                    close(new.mass_kg(), old.mass_kg(), 1e-16);
                    close(
                        species_mass(new, "A") - species_mass(old, "A"),
                        sign * j * area_time,
                        1e-16,
                    );
                    let temperature = old.state().state().temperature
                        + sign * q * area_time / (old.mass_kg() * 2.5 * r / MASS);
                    close(new.state().state().temperature, temperature, 2e-8);
                    assert!(matches!(
                        next.chambers()[&id].audit_last_window().unwrap(),
                        ModuleWindowAudit::Unavailable { .. }
                    ));
                }
                if scale == 0.0 {
                    assert!(report.flow.transport.exchange.fluxes.is_empty());
                    assert!(!report.flow.transport.exchange.heat.is_empty());
                }
                assert!(
                    report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j
                );
                assert!(report.flow.transport.entropy_change_j_k > 0.0);
                let last = molecular.face().law().last_resolved();
                close(
                    last.temperature().value(),
                    next.chambers()[&20].gas().state().state().temperature,
                    2e-8,
                );
                close(
                    last.heats()[0].value(),
                    scale * (20000.0 + 100.0 * last.temperature().value()),
                    1e-8,
                );
                assert_ne!(state.circuit(), next.circuit());
                assert_eq!(state.time_s(), 0.0);
            })
        });
    }
}

#[test]
fn g1_isothermal_dufour_flux_is_counted_once_and_zero_cross_recovers_the_original_face() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.8, 400.0, 1e5), (0.2, 400.0, 1e5)], true),
            );
            let drives = drives(state.chambers());
            let run = |op: &mut Molecular| {
                state
                    .molecular_transport_step_candidate(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &[],
                        op,
                        &[molecular_link()],
                        GasTransportPaths::default(),
                        molecular_controls(),
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap()
            };
            let (cross, report) = run(&mut molecular);
            let j = 1e5 / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0) * DIFFUSIVITY * MASS * 0.6
                / molecular_link().distance_m;
            let q = 100000.0 * j;
            close(
                report.flow.links[0].mass_left_to_right_kg[0],
                j * molecular_link().area_m2 * DT,
                1e-16,
            );
            close(
                report.flow.links[0].heat_correction_left_to_right_j,
                q * molecular_link().area_m2 * DT,
                1e-12,
            );
            assert!(cross.chambers()[&10].gas().state().state().temperature < 400.0);
            assert!(cross.chambers()[&20].gas().state().state().temperature > 400.0);
            let (_, mut zero) = molecular_fixture(MASS, 0.0, 0.0);
            let mut bare = ParcelMolecularTransport::new(
                PlanarTransport::new(
                    zero.face().gas().clone(),
                    zero.face().law().law().clone(),
                    || false,
                )
                .unwrap(),
            );
            let (zero_state, zero_report) = run(&mut zero);
            let (bare_state, bare_report) = state
                .molecular_transport_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &[],
                    &mut bare,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(zero_report.flow, bare_report.flow);
            assert_eq!(zero_state.chambers(), bare_state.chambers());
            assert_eq!(zero_state.thermal().cells(), bare_state.thermal().cells());
            assert_eq!(zero_state.circuit(), bare_state.circuit());
            assert_eq!(
                zero_report.flow.links[0].heat_correction_left_to_right_j,
                0.0
            );
            close(
                zero_report.flow.links[0].mass_left_to_right_kg[0],
                report.flow.links[0].mass_left_to_right_kg[0],
                1e-16,
            );
        })
    });
}

#[test]
fn g3_unequal_masses_and_distinct_nozzle_paths_share_one_inventory_and_feedback_solve() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(2.0 * MASS, 1.0, -5000.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.7, 390.0, 1e5), (0.3, 410.0, 1.2e5)], false),
            );
            let drives = drives(state.chambers());
            let nozzles = [SnapshotNozzleLink {
                id: 901,
                parcels: [10, 20],
                nozzle: GasNozzle {
                    area_m2: 2e-8,
                    discharge_coefficient: 0.8,
                },
                nozzle_solve: chemical().gas,
            }];
            let paths = || GasTransportPaths {
                nozzles: &nozzles,
                ..Default::default()
            };
            let parcels = state
                .chambers()
                .iter()
                .map(|(&id, g)| (id, g.gas().clone()))
                .collect();
            let endpoints = drives
                .iter()
                .map(|d| (d.chamber, d.piston_solve.gas))
                .collect();
            let (transferred, original) = molecular
                .step_candidate(
                    &parcels,
                    &[molecular_link()],
                    paths(),
                    DT,
                    &endpoints,
                    molecular_controls(),
                    || false,
                )
                .unwrap()
                .into_parts();
            let staged = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                transferred
                    .into_iter()
                    .map(|(id, g)| (id, GasPiston::new(g, 1.0, 0.01, 0.0).unwrap()))
                    .collect(),
            );
            let (expected, _) = staged
                .step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let inventory = |state: &MaterialExchangeState<'_>, name| {
                state
                    .chambers()
                    .values()
                    .map(|g| species_mass(g.gas(), name))
                    .sum::<f64>()
            };
            let before = [inventory(&state, "A"), inventory(&state, "B")];
            let report = state
                .advance_with_molecular_transport(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &mut molecular,
                    &[molecular_link()],
                    paths(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(report.flow, original);
            assert_eq!(state.chambers(), expected.chambers());
            assert_eq!(state.circuit(), expected.circuit());
            assert_eq!(state.thermal().cells(), expected.thermal().cells());
            assert!(report.flow.transport.nozzles[0].transferred_mass_kg > 0.0);
            close(
                report.flow.links[0].mass_left_to_right_kg.iter().sum(),
                0.0,
                1e-15,
            );
            for (i, name) in ["A", "B"].into_iter().enumerate() {
                close(inventory(&state, name), before[i], 1e-13);
            }
            assert!(report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
            assert!(
                report.total_mass_residual_kg.abs()
                    <= molecular_controls().physical.exchange.mass_tolerance_kg
            );
            let first_temperature = state.chambers()[&10].gas().state().state().temperature;
            let next = state
                .advance_with_molecular_transport(
                    &prepared,
                    cx,
                    2.0 * DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &mut molecular,
                    &[molecular_link()],
                    paths(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_ne!(
                next.flow.links[0].face.enthalpy_temperature_k,
                report.flow.links[0].face.enthalpy_temperature_k
            );
            assert_ne!(
                state.chambers()[&10].gas().state().state().temperature,
                first_temperature
            );
            assert_eq!(state.time_s(), 2.0 * DT);
        })
    });
}

#[test]
fn g4_source_and_late_feedback_refusals_and_cancellation_preserve_every_physical_owner() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, gate| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.5, 390.0, 1e5), (0.5, 410.0, 1e5)], false),
            );
            let drives = drives(state.chambers());
            let gases = state.chambers().clone();
            let cells = state.thermal().cells().clone();
            let electrical = state.circuit().clone();
            let (_, mut wrong) = molecular_fixture(2.0 * MASS, 1.0, 0.0);
            assert!(matches!(
                state.advance_with_molecular_transport(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &mut wrong,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2]
                ),
                Err(MaterialFlowError::Molecular(
                    ParcelMolecularError::SourceMismatch(_)
                ))
            ));
            let mut limited = controls();
            limited.network.maximum_circuit_trials = 1;
            assert!(
                state
                    .advance_with_molecular_transport(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &gas_fixture::patches(),
                        &mut molecular,
                        &[molecular_link()],
                        GasTransportPaths::default(),
                        molecular_controls(),
                        limited,
                        &[0.0; 2]
                    )
                    .is_err()
            );
            assert_eq!(state.chambers(), &gases);
            assert_eq!(state.thermal().cells(), &cells);
            assert_eq!(state.circuit(), &electrical);
            assert_eq!(state.time_s(), 0.0);
            for (&id, gas) in state.chambers() {
                assert_eq!(
                    gas.last_window_transfers(),
                    gases[&id].last_window_transfers()
                );
                assert_eq!(
                    gas.last_window_duration_s(),
                    gases[&id].last_window_duration_s()
                );
            }
            let (preview, _) = state
                .molecular_transport_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            state
                .advance_with_molecular_transport(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_eq!(state.chambers(), preview.chambers());
            assert_eq!(state.circuit(), preview.circuit());
            assert_eq!(state.thermal().cells(), preview.thermal().cells());
            gate.request();
            assert!(matches!(
                state.advance_with_molecular_transport(
                    &prepared,
                    cx,
                    2.0 * DT,
                    &[],
                    &drives,
                    &[],
                    &gas_fixture::patches(),
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2]
                ),
                Err(MaterialFlowError::Cancelled)
            ));
            assert_eq!(state.chambers(), preview.chambers());
            assert_eq!(state.circuit(), preview.circuit());
            assert_eq!(state.thermal().cells(), preview.thermal().cells());
            assert_eq!(state.time_s(), DT);
            for (&id, gas) in state.chambers() {
                assert_eq!(
                    gas.last_window_transfers(),
                    preview.chambers()[&id].last_window_transfers()
                );
                assert_eq!(
                    gas.last_window_duration_s(),
                    preview.chambers()[&id].last_window_duration_s()
                );
            }
        })
    });
}
