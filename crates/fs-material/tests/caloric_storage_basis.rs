//! G1 source-capacity heating and G0 refusal of unsupported open-mass semantics.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::conduction::LayerConductionConfig;
use fs_material::thermal::transport::{
    conduction::{CaloricConductionLayer, IntegratedCaloricBoundary},
    storage::{
        HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw, PreparedHeatStorage,
        conduction::HeatStorageConductionLaw,
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
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, ThermalLink, ThermalStorageBasis,
    TimedEvent,
    advection::{Advected, MassBoundary, MassFlow},
    conduction::{ConductionObserver, ConductiveLink, ConductiveReservoir, ConductivityLaw},
    conduction::implicit::{
        ImplicitConductionControls, ImplicitConductionError, ImplicitConductionGraph,
        ImplicitConductionLaw, ImplicitConductionLinearControls, ImplicitConductionWorkspace,
        boundary::{ImplicitBoundaryGeometry, ImplicitBoundaryLaw},
    },
    evaporation::LiquidCellExchange,
};
use fs_thermal::evaporation::LiquidSurfaceState;

fn sources() -> (MaterialCard, PreparedVolumetricThermalTransport, PreparedHeatStorage) {
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
        (conductivity.clone(), Dims([1, 1, -3, -1, 0, 0]), [1.0, 3.0]),
    ] {
        let knots = if key == conductivity {
            // An independent source breakpoint must survive the assembled C(T).
            vec![(300.0, 1.0), (307.0, 2.0), (320.0, 3.0)]
        } else {
            [300.0, 304.25, 310.0, 317.5, 320.0].into_iter()
                .map(|t| (t, values[0] + (t - 300.0) / 20.0 * (values[1] - values[0])))
                .collect()
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
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
    (card, path, storage)
}

#[test]
fn g1_indexed_sensible_storage_keeps_mass_scaling_and_small_heat() {
    let (_, path, storage) = sources();
    for datum in [300.0_f64, 303.125, 310.0, 319.0, 320.0] {
        let volume = FixedVolumeThermalLaw::try_new(&path, 10.0, datum).unwrap();
        let assembly = HeatStorageThermalLaw::try_new(&storage, 2.0, datum).unwrap();
        for target in [300.0, 304.25, 309.75, 310.0, 317.5, 319.999, 320.0,
            datum.next_down(), datum, datum.next_up()]
        {
            if !(300.0..=320.0).contains(&target) { continue; }
            // Both amount bindings have c(T)=2+0.2*(T-300) J/(kg K).
            let expected = (2.0 + 0.1 * ((datum - 300.0) + (target - 300.0))) * (target - datum);
            for (law, h, direct) in [
                (&volume as &dyn ThermalLaw, volume.specific_storage_j_kg(target).unwrap(),
                    path.storage_change_j_m3(datum, target).unwrap() / 10.0),
                (&assembly as &dyn ThermalLaw, assembly.specific_storage_j_kg(target).unwrap(),
                    storage.energy_change_j(datum, target).unwrap() / 2.0),
            ] {
                assert!((h - expected).abs() <= 16.0 * f64::EPSILON * expected.abs());
                assert!((law.temperature_k(h).unwrap() - target).abs()
                    <= 8.0 * f64::EPSILON * target);
                if target != datum { assert_ne!(h, 0.0); }
                if target == 300.0 || target == 320.0 || (target - datum).abs() < 1e-10 {
                    assert_eq!(h, direct);
                }
            }
        }
        for outside in [299.0, 321.0, f64::NAN, f64::INFINITY] {
            assert!(volume.specific_storage_j_kg(outside).is_err());
            assert!(assembly.specific_storage_j_kg(outside).is_err());
        }
    }
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
    let (_, path, storage) = sources();
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
    let (_, path, storage) = sources();
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

fn temperature_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            t,
        )
        .unwrap()
}

#[test]
fn g1_assembly_conduction_uses_total_capacity_and_original_source_tangents() {
    let (card, _, storage) = sources();
    let caloric = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0)
        .unwrap()
        .with_sensible_entropy();
    let conductivity = PropertyKey::new("conductivity", Dims([1, 1, -3, -1, 0, 0]));
    let lower = temperature_point(300.0);
    let upper = temperature_point(320.0);
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let law = HeatStorageConductionLaw::try_new(
        &card, &caloric, &conductivity, "T", &lower, &upper, selection,
    )
    .unwrap();
    assert!(core::ptr::eq(law.caloric(), &caloric));
    assert_eq!(law.conductivity_card_identity(), card.content_hash());
    assert_eq!(
        law.caloric_law().storage_basis(),
        ThermalStorageBasis::DeclaredSensibleStorage
    );
    assert_eq!(law.enthalpy_bounds_j_kg(), caloric.specific_storage_span_j_kg());
    assert_eq!(
        law.source_line().knots.iter().map(|knot| knot.0).collect::<Vec<_>>(),
        vec![300.0, 307.0, 320.0]
    );
    for t in [300.0, 304.25, 306.999, 307.001, 310.0, 317.5, 320.0] {
        // Independent amount-weighted capacity: 1 kg * c(T) + 0.1 m³ * c_vol(T).
        // C_total=4+0.4*(T-300) J/K, while this body has fixed mass 2 kg.
        let theta = t - 300.0;
        let e = (4.0 * theta + 0.2 * theta * theta) / 2.0;
        let d_t_de = 2.0 / (4.0 + 0.4 * theta);
        let (k, d_k_d_t) = if t < 307.0 {
            (1.0 + theta / 7.0, 1.0 / 7.0)
        } else {
            (2.0 + (t - 307.0) / 13.0, 1.0 / 13.0)
        };
        let state = law.state_at_specific_enthalpy(e).unwrap();
        assert!((state.temperature_k - t).abs() < 2e-12);
        assert!((state.temperature_slope_kg_k_j - d_t_de).abs() < 2e-14);
        assert!((state.conductivity_w_m_k - k).abs() < 2e-13);
        assert!((state.conductivity_slope_kg_w_m_k_j - d_k_d_t * d_t_de).abs() < 2e-14);
        assert_eq!(law.conductivity_w_m_k(e).unwrap(), state.conductivity_w_m_k);
        let observable = law.time_observables(e).unwrap();
        assert_eq!(observable.temperature_k, state.temperature_k);
        assert_eq!(observable.phase_mass_fraction, None);
    }
    for outside in [-1.0, 81.0, f64::NAN, f64::INFINITY] {
        assert!(law.state_at_specific_storage(outside).is_err());
        assert!(law.conductivity_w_m_k(outside).is_err());
    }
    assert!(HeatStorageConductionLaw::try_new(
        &card, &caloric, &conductivity, "T", &temperature_point(301.0), &upper, selection,
    ).is_err());
    let wrong_role = PropertyKey::new("selected-sensible-capacity", Dims([2, 0, -2, -1, 0, 0]));
    assert!(HeatStorageConductionLaw::try_new(
        &card, &caloric, &wrong_role, "T", &lower, &upper, selection,
    ).is_err());

    // Register the very same caloric owner supplied to the implicit adapter;
    // native event admission must still reject rescaling its fixed inventory.
    let mut changed = initial(law.caloric_law(), 1, 0.0).cells()[&CellId(0)].clone();
    changed.mass_kg = 3.0;
    let pending = DepositionSimulation::new(
        vec![TimedEvent { time_s: 0.0, event: Event::Birth(changed) }],
        &[(MaterialId(1), law.caloric_law())],
        DepositionBudget { max_cells: 1, max_links: 0, max_events: 1, max_substeps: 16 },
    ).unwrap();
    assert!(pending.advance_to(0.0, step(), || Ok(())).is_err());
    assert!(pending.cells().is_empty());
}

#[test]
fn g1_assembly_native_graph_conserves_nonlinear_storage_across_a_conductivity_knot() {
    let (card, _, storage) = sources();
    let caloric = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0)
        .unwrap()
        .with_sensible_entropy();
    let law = HeatStorageConductionLaw::try_new(
        &card,
        &caloric,
        &PropertyKey::new("conductivity", Dims([1, 1, -3, -1, 0, 0])),
        "T",
        &temperature_point(300.0),
        &temperature_point(320.0),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    ).unwrap();
    let mut events: Vec<_> = initial(&caloric, 4, 0.0).cells().values().cloned()
        .map(|mut cell| {
            cell.specific_h_j_kg = if cell.id == CellId(0) { 136.8 / 2.0 } else { 8.8 / 2.0 };
            cell.geometry.surface_area_m2 = 4.0;
            TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
        })
        .collect();
    for id in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id, a: CellId(0), b: CellId(id + 1), area_m2: 1.0, conductance_w_k: 0.0,
            }),
        });
    }
    let old = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &caloric)],
        DepositionBudget { max_cells: 4, max_links: 3, max_events: 7, max_substeps: 16 },
    ).unwrap().advance_to(0.0, step(), || Ok(())).unwrap().0;
    let geometry: Vec<_> = (0..3).map(|id| ConductiveLink {
        link_id: id, half_distances_m: [0.5, 0.5], contact_resistance_m2_k_w: 0.2,
    }).collect();
    let observer = ConductionObserver::new(
        &old, &[(MaterialId(1), &law)], &geometry, &[], || Ok(()),
    ).unwrap();
    let graph = ImplicitConductionGraph::new(
        &old, &observer, &[(MaterialId(1), &law)],
        ImplicitConductionLinearControls::default(), || Ok(()),
    ).unwrap();
    let dt = 20.0;
    let controls = ImplicitConductionControls {
        maximum_dt_s: dt, cell_tolerance_j: 1e-10, total_tolerance_j: 1e-9,
        ..ImplicitConductionControls::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(4).unwrap();
    let (next, report) = graph.advance_to(&old, dt, controls, &mut workspace, || Ok(())).unwrap();

    // Independent scalar energy balance for this symmetric three-leaf graph.
    // Solve q=dt*G(T_hot,T_cold)*(T_hot-T_cold) for heat to EACH cold body.
    // Invert E=4*theta+0.2*theta² without using the material law or native graph.
    let temperature = |energy: f64| 300.0 + 2.0 * energy / (4.0 + (16.0 + 0.8 * energy).sqrt());
    let conductivity = |t: f64| if t < 307.0 {
        1.0 + (t - 300.0) / 7.0
    } else {
        2.0 + (t - 307.0) / 13.0
    };
    let mut bracket = [0.0, 32.0];
    for _ in 0..80 {
        let heat = 0.5 * (bracket[0] + bracket[1]);
        let hot = temperature(136.8 - 3.0 * heat);
        let cold = temperature(8.8 + heat);
        let conductance = 1.0 / (0.5 / conductivity(hot) + 0.2 + 0.5 / conductivity(cold));
        if heat - dt * conductance * (hot - cold) > 0.0 {
            bracket[1] = heat;
        } else {
            bracket[0] = heat;
        }
    }
    let heat = 0.5 * (bracket[0] + bracket[1]);
    assert!(temperature(8.8 + heat) > 307.0);
    for cell in next.cells().values() {
        let energy = if cell.id == CellId(0) { 136.8 - 3.0 * heat } else { 8.8 + heat };
        assert_eq!(cell.mass_kg, 2.0);
        assert!((cell.mass_kg * cell.specific_h_j_kg - energy).abs() < 1e-8);
        assert!((next.temperature_k(cell.id).unwrap() - temperature(energy)).abs() < 1e-9);
    }
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.energy.external_heat_j, 0.0);
    assert!(report.energy.stored_energy_change_j.abs() < 1e-9);
    assert!(report.energy.residual_j.abs() < 1e-9);
    assert!(report.entropy.thermal_storage_change_j_k > 0.0);
    assert!(report.entropy.conduction_generation_j_k > 0.0);
    let replay = graph.advance_to(&old, dt, controls, &mut workspace, || Ok(())).unwrap();
    assert_eq!(replay.0.cells(), next.cells());
    assert_eq!(replay.1, report);
    assert!(matches!(
        graph.advance_to(&old, dt, controls, &mut workspace, || Err(DepositionError::Cancelled)),
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, 136.8 / 2.0);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 8.8 / 2.0);
}

#[test]
fn g1_integrated_storage_cools_through_source_conductivity_and_a_supported_film() {
    let (card, path, storage) = sources();
    let volume = FixedVolumeThermalLaw::try_new(&path, 10.0, 300.0)
        .unwrap()
        .with_sensible_entropy();
    let assembly = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0)
        .unwrap()
        .with_sensible_entropy();
    let assembly_conduction = HeatStorageConductionLaw::try_new(
        &card,
        &assembly,
        &PropertyKey::new("conductivity", Dims([1, 1, -3, -1, 0, 0])),
        "T",
        &temperature_point(300.0),
        &temperature_point(320.0),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    ).unwrap();
    let layers = [
        CaloricConductionLayer::from_fixed_volume(&volume).unwrap(),
        CaloricConductionLayer::from_heat_storage(&assembly_conduction).unwrap(),
    ];
    let caloric: [&dyn ThermalLaw; 2] = [&volume, &assembly];
    let conduction: [&dyn ConductivityLaw; 2] = [&volume, &assembly_conduction];
    let implicit: [&dyn ImplicitConductionLaw; 2] = [&volume, &assembly_conduction];
    let geometry = ImplicitBoundaryGeometry {
        area_m2: 0.2, half_distance_m: 0.5, surface_resistance_m2_k_w: 0.16,
    };
    let flux_controls = LayerConductionConfig {
        max_iterations: 128, absolute_heat_tolerance_w: 1e-11, relative_heat_tolerance: 1e-12,
    };
    let (initial_energy, dt, bath) = (136.8, 5.0, 295.0);

    // Independent coupled storage/film solution. K is the exact piecewise
    // integral of the fixture k(T), and T(E) inverts total C=4+0.4*(T-300).
    // The 0.4 m * 0.5 m * 1 m body has V=0.2 m³, A=0.2 m² and L_half=0.5 m.
    // In the surface coordinate, P=(Ts-Tbath)/R and E_new=E_old-dt*P,
    // with total film resistance R=0.16/A=0.8 K/W.
    let integrated_k = |t: f64| {
        if t <= 307.0 {
            let x = t - 300.0;
            x + x * x / 14.0
        } else {
            let x = t - 307.0;
            10.5 + 2.0 * x + x * x / 26.0
        }
    };
    let temperature = |energy: f64| 300.0 + 2.0 * energy / (4.0 + (16.0 + 0.8 * energy).sqrt());
    let k = |t: f64| if t < 307.0 { 1.0 + (t - 300.0) / 7.0 } else { 2.0 + (t - 307.0) / 13.0 };
    let mut bracket = [300.0, 307.0];
    for _ in 0..80 {
        let surface = 0.5 * (bracket[0] + bracket[1]);
        let power = (surface - bath) / 0.8;
        let bulk = temperature(initial_energy - dt * power);
        if integrated_k(bulk) - integrated_k(surface) - 2.5 * power > 0.0 {
            bracket[0] = surface;
        } else {
            bracket[1] = surface;
        }
    }
    let expected_surface = 0.5 * (bracket[0] + bracket[1]);
    let expected_power = (expected_surface - bath) / 0.8;
    let expected_energy = initial_energy - dt * expected_power;
    let expected_temperature = temperature(expected_energy);
    assert!(expected_surface > 300.0 && expected_surface < 307.0);
    assert!(expected_temperature > 307.0);
    let point_k_power = (expected_temperature - bath) / (2.5 / k(expected_temperature) + 0.8);
    assert!((point_k_power - expected_power).abs() > 0.5);

    for (side, layer) in layers.iter().enumerate() {
        assert!(core::ptr::addr_eq(
            layer.caloric_law() as *const dyn ThermalLaw,
            caloric[side] as *const dyn ThermalLaw,
        ));
        assert_eq!(layer.caloric_law().storage_basis(), ThermalStorageBasis::DeclaredSensibleStorage);
        assert_eq!(layer.energy_bounds_j_kg(), implicit[side].enthalpy_bounds_j_kg());
        assert_eq!(
            layer.source_line().knots.iter().map(|knot| knot.0).collect::<Vec<_>>(),
            vec![300.0, 307.0, 320.0],
        );
        let mut cell = initial(caloric[side], 1, 0.0).cells()[&CellId(0)].clone();
        cell.specific_h_j_kg = initial_energy / cell.mass_kg;
        cell.geometry.end_m = [1.0, 0.0, 0.0];
        cell.geometry.surface_area_m2 = 2.2;
        cell.boundary.reservoir = Some(Reservoir {
            area_m2: 0.2, conductance_w_k: 999.0, temperature_k: bath,
        });
        let old = DepositionSimulation::new(
            vec![TimedEvent { time_s: 0.0, event: Event::Birth(cell.clone()) }],
            &[(MaterialId(1), caloric[side])],
            DepositionBudget { max_cells: 1, max_links: 0, max_events: 1, max_substeps: 16 },
        ).unwrap().advance_to(0.0, step(), || Ok(())).unwrap().0;
        // No integrated layer may replace the original amount/volume admission.
        cell.mass_kg = 3.0;
        assert!(layer.caloric_law().validate_cell(&cell).is_err());
        let boundary = IntegratedCaloricBoundary::new(layer, geometry, flux_controls).unwrap();
        let observer = ConductionObserver::new(
            &old, &[(MaterialId(1), conduction[side])], &[],
            &[ConductiveReservoir {
                cell_id: CellId(0), half_distance_m: 0.5, surface_resistance_m2_k_w: 0.16,
            }], || Ok(()),
        ).unwrap();
        let graph = ImplicitConductionGraph::new(
            &old, &observer, &[(MaterialId(1), implicit[side])],
            ImplicitConductionLinearControls::default(), || Ok(()),
        ).unwrap().with_boundary_laws(&old, &[(CellId(0), &boundary)], &[], || Ok(())).unwrap();
        let controls = ImplicitConductionControls {
            maximum_dt_s: dt, maximum_halvings: 0,
            cell_tolerance_j: 1e-8, total_tolerance_j: 1e-7,
            ..ImplicitConductionControls::default()
        };
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = graph.advance_to(&old, dt, controls, &mut workspace, || Ok(())).unwrap();
        let energy = next.cells()[&CellId(0)].specific_h_j_kg;
        assert!((next.temperature_k(CellId(0)).unwrap() - expected_temperature).abs() < 2e-8);
        assert!((2.0 * energy - expected_energy).abs() < 2e-7);
        let response = boundary.response_at_specific_energy(energy, bath, &mut || Ok(())).unwrap();
        assert!((response.material_face_temperature_k - expected_surface).abs() < 2e-8);
        assert!((response.power_to_reservoir_w - expected_power).abs() < 2e-8);
        let tangent = boundary.state_at_specific_enthalpy(energy, bath, geometry, &mut || Ok(())).unwrap();
        let expected_tangent = 0.2 * k(expected_temperature) / (0.5 + 0.16 * k(expected_surface))
            * 2.0 / (4.0 + 0.4 * (expected_temperature - 300.0));
        assert!((tangent.power_enthalpy_derivative_w_kg_j - expected_tangent).abs() < 2e-9);
        assert_eq!(report.energy.substeps, 1);
        assert!((report.energy.external_heat_j + dt * expected_power).abs() < 2e-7);
        assert!((report.energy.stored_energy_change_j + dt * expected_power).abs() < 2e-7);
        assert!(report.energy.residual_j.abs() < 1e-7);
        let entropy = 0.4 * (expected_temperature - 318.0)
            - 116.0 * (expected_temperature / 318.0).ln();
        assert!((report.entropy.thermal_storage_change_j_k - entropy).abs() < 2e-9);
        assert!((report.entropy.reservoir_entropy_input_j_k + dt * expected_power / bath).abs() < 2e-9);
        assert!(entropy - report.entropy.reservoir_entropy_input_j_k > 0.0);

        // A bath outside source support is admissible only while the actual
        // material surface stays supported. Removing the film must refuse.
        let perfect = IntegratedCaloricBoundary::new(
            layer, ImplicitBoundaryGeometry { surface_resistance_m2_k_w: 0.0, ..geometry }, flux_controls,
        ).unwrap();
        assert!(perfect.response_at_specific_energy(energy, bath, &mut || Ok(())).is_err());
        let replay = graph.advance_to(&old, dt, controls, &mut workspace, || Ok(())).unwrap();
        assert_eq!(replay.0.cells(), next.cells());
        assert_eq!(replay.1, report);
        assert!(matches!(
            graph.advance_to(&old, dt, controls, &mut workspace, || Err(DepositionError::Cancelled)),
            Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, initial_energy / 2.0);
    }
}
