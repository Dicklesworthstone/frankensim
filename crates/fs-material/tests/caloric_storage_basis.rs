//! G1 source-capacity heating and G0 refusal of unsupported open-mass semantics.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{
    storage::{
        HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw, PreparedHeatStorage,
    },
    volumetric::{
        PreparedVolumetricThermalTransport, VOLUMETRIC_HEAT_CAPACITY_DIMS,
        VolumetricThermalPathQuery, timeline::FixedVolumeThermalLaw,
    },
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalStorageBasis,
    TimedEvent,
    advection::{Advected, MassBoundary, MassFlow},
    evaporation::LiquidCellExchange,
};
use fs_thermal::evaporation::LiquidSurfaceState;

fn sources() -> (PreparedVolumetricThermalTransport, PreparedHeatStorage) {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let capacity = PropertyKey::new("selected-sensible-capacity", Dims([2, 0, -2, -1, 0, 0]));
    let volumetric = PropertyKey::new("volumetric-capacity", VOLUMETRIC_HEAT_CAPACITY_DIMS);
    let conductivity = PropertyKey::new("conductivity", Dims([1, 1, -3, -1, 0, 0]));
    let mut claims = ClaimSet::new();
    for (key, dims, values) in [
        (capacity.clone(), Dims([2, 0, -2, -1, 0, 0]), [2.0, 6.0]),
        (
            volumetric.clone(),
            VOLUMETRIC_HEAT_CAPACITY_DIMS,
            [20.0, 60.0],
        ),
        (conductivity.clone(), Dims([1, 1, -3, -1, 0, 0]), [1.0, 1.0]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (320.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    300.0,
                    320.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored sensible-capacity integration fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal body".into(),
            phase: "solid".into(),
            process: "fixed amount; no pressure or flow-work claim".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", temperature, t)
            .unwrap()
    };
    let (lower, upper) = (point(300.0), point(320.0));
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let path = PreparedVolumetricThermalTransport::try_new(
        &card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &volumetric,
            conductivity: &conductivity,
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selections: [selection; 2],
        },
    )
    .unwrap();
    let storage = PreparedHeatStorage::try_new(&[
        HeatStorageComponentQuery {
            id: 1,
            card: &card,
            capacity: &capacity,
            amount: HeatStorageAmount::MassKg(1.0),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection,
        },
        HeatStorageComponentQuery {
            id: 2,
            card: &card,
            capacity: &volumetric,
            amount: HeatStorageAmount::ReferenceVolumeM3(0.1),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection,
        },
    ])
    .unwrap();
    (path, storage)
}

fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0,
        energy_tolerance_j: 1e-10,
    }
}

fn initial(law: &dyn ThermalLaw, cells: usize, source_w: f64) -> DepositionSimulation<'_> {
    let events = (0..cells)
        .map(|id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id as u64),
                segment_id: id as u64,
                body_id: id as u64,
                material_id: MaterialId(1),
                mass_kg: 2.0,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64 + 0.5, 0.0, 0.0],
                    volume_m3: 0.2,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: cells,
            max_links: 0,
            max_events: cells,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(0.0, step(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_declared_sensible_capacity_keeps_nonlinear_heat_and_entropy_on_the_native_timeline() {
    let (path, storage) = sources();
    let volume_law = FixedVolumeThermalLaw::try_new(&path, 10.0, 300.0).unwrap();
    let assembly_law = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0).unwrap();
    for law in [&volume_law as &dyn ThermalLaw, &assembly_law] {
        assert_eq!(
            law.storage_basis(),
            ThermalStorageBasis::DeclaredSensibleStorage
        );
        assert_eq!(law.specific_entropy_change_j_kg_k(0.0, 30.0).unwrap(), None);
    }
    let volume_law = volume_law.with_sensible_entropy();
    let assembly_law = assembly_law.with_sensible_entropy();
    for law in [&volume_law as &dyn ThermalLaw, &assembly_law] {
        assert_eq!(
            law.storage_basis(),
            ThermalStorageBasis::DeclaredSensibleStorage
        );
        let state = initial(law, 1, 60.0);
        let (next, energy, entropy) = state
            .advance_to_with_entropy(1.0, step(), || Ok(()))
            .unwrap();
        // Both inventories have C(T)=4+0.4*(T-300) J/K. Integral from 300 to
        // 310 K is 60 J; freezing the initial capacity would incorrectly give 315 K.
        assert!((next.temperature_k(CellId(0)).unwrap() - 310.0).abs() < 1e-12);
        assert_eq!(next.cells()[&CellId(0)].specific_h_j_kg, 30.0);
        assert_eq!(energy.stored_energy_change_j, 60.0);
        assert_eq!(energy.external_heat_j, 60.0);
        assert_eq!(energy.residual_j, 0.0);
        let exact_entropy = 4.0 - 116.0 * (310.0_f64 / 300.0).ln();
        assert!((entropy.thermal_storage_change_j_k - exact_entropy).abs() < 2e-14);
        assert_eq!(state.cells()[&CellId(0)].specific_h_j_kg, 0.0);
        let replay = state
            .advance_to_with_entropy(1.0, step(), || Ok(()))
            .unwrap();
        assert_eq!(replay.0.cells(), next.cells());
        assert_eq!(replay.1, energy);
        assert_eq!(replay.2, entropy);
    }
}

#[test]
fn g0_declared_capacity_cannot_supply_stream_or_evaporation_enthalpy() {
    let (path, storage) = sources();
    let volume_law = FixedVolumeThermalLaw::try_new(&path, 10.0, 300.0)
        .unwrap()
        .with_sensible_entropy();
    let assembly_law = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0)
        .unwrap()
        .with_sensible_entropy();
    for law in [&volume_law as &dyn ThermalLaw, &assembly_law] {
        let mut state = initial(law, 2, 0.0);
        let before = state.cells().clone();
        for rate in [0.1, -0.1] {
            let flows = vec![
                MassFlow {
                    id: 1,
                    a: CellId(0),
                    b: CellId(1),
                    kg_s: rate,
                },
                MassFlow {
                    id: 2,
                    a: CellId(1),
                    b: CellId(0),
                    kg_s: rate,
                },
            ];
            let observer = Advected::new((), flows, 2).unwrap();
            assert!(matches!(
                state.advance_observed(&observer, 1.0, step(), || Ok(())),
                Err(DepositionError::InvalidInput(
                    "mass flow requires an explicitly enthalpy-based storage law"
                ))
            ));
            assert_eq!(observer.audit().substeps, 0);
        }
        let boundary = Advected::new((), vec![], 0)
            .unwrap()
            .with_boundaries(
                vec![
                    MassBoundary::Inlet {
                        id: 1,
                        cell: CellId(0),
                        kg_s: 0.1,
                        material_id: MaterialId(1),
                        specific_h_j_kg: 30.0,
                    },
                    MassBoundary::Outlet {
                        id: 2,
                        cell: CellId(0),
                        kg_s: 0.1,
                    },
                ],
                2,
            )
            .unwrap();
        assert!(matches!(
            state.advance_observed(&boundary, 1.0, step(), || Ok(())),
            Err(DepositionError::InvalidInput(
                "mass flow requires an explicitly enthalpy-based storage law"
            ))
        ));
        let exchange = LiquidCellExchange {
            expected: before[&CellId(0)].clone(),
            interval_exposed_area_m2: 1.0,
            replacement: LiquidSurfaceState {
                mass_kg: 1.9,
                specific_enthalpy_j_kg: 0.0,
                area_m2: 1.0,
            },
            pure_liquid_enthalpy_range_j_kg: [0.0, 80.0],
        };
        assert!(matches!(
            state.apply_liquid_exchanges(&[exchange], || Ok(())),
            Err(DepositionError::InvalidInput(
                "evaporation requires enthalpy storage and an explicit open-mass model"
            ))
        ));
        assert_eq!(state.cells(), &before);
        assert_eq!(state.time_s(), 0.0);
    }
}
