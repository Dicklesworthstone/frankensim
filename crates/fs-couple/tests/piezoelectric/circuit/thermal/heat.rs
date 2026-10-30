//! Synthetic ceramic/support sources with independent two-temperature balances.
//! These numerical fixtures do not establish measured ceramic calibration.
use super::*;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::transient::deposition::{
    contact::{ContactCell, MaterialContactLink},
    radiation::{MaterialRadiation, RadiatingCell},
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, GrayDiffuseEnclosure, RadiationSurface, SURFACE_EMISSIVITY_PROPERTY,
    SurfaceEmissivity, ViewFactorMatrix,
};
use fs_couple::circuit_enthalpy::contacts::{
    heat::MaterialHeatNetwork, material::ElectricalContactNetwork,
};
use fs_couple::circuit_enthalpy::material::ConductionLink;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SelectionPolicy, SurfaceSpec, SystemContext};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS;
use fs_material::thermal::interface::ThermalInterfacePoint;

const SUPPORT: CellId = CellId(8);
const IDS: [CellId; 2] = [CELL, SUPPORT];
const MIDS: [MaterialId; 2] = [MATERIAL, MaterialId(4)];
const CONTACT_AREA: f64 = 0.1;
const RADIATING_AREA: f64 = 0.01;
const INITIAL: [f64; 2] = [300.0, 340.0];
const OLD: [f64; 3] = [0.02, 0.01, 0.1];
const FORCE: f64 = 0.03;
const DT: f64 = 0.1;

fn heat_card(emissivity: f64, support: bool) -> MaterialCard {
    let original = source_card(0.2, SLOPES, 0.0);
    let mut claims = original.claims().clone();
    for (name, dims, value, slope) in [
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            0.2,
            0.001,
        ),
        (
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            emissivity,
            0.002 * emissivity,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, value - 20.0 * slope), (360.0, value + 60.0 * slope)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored ceramic/support heat fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let mut id = original.id().clone();
    if support {
        id.chemistry = "synthetic thermal support".into();
    }
    MaterialCard::assemble(id, claims, vec![]).unwrap()
}

fn interface(cards: &[MaterialCard; 2], resistance: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Curve {
                abscissa: "aT".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 0.9 * resistance), (360.0, 1.3 * resistance)],
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("aT", temperature_quantity(), 280.0, 360.0)
                .with_quantity("bT", temperature_quantity(), 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored two-sided thermal contact".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "ceramic".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "support".into(),
        },
        SystemContext {
            medium: "fixed interface".into(),
            third_body: None,
            environment: "authored".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn with_case<R>(
    emissivity: f64,
    resistance: f64,
    extra_paths: bool,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        PiezoelectricCircuit,
        &PiezoelectricThermalBinding,
    ) -> R,
) -> R {
    let cards = [heat_card(emissivity, false), heat_card(0.6, true)];
    let curves = cards.each_ref().map(|card| caloric_curve(card, false));
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: IDS[i],
        material_id: MIDS[i],
        curve: &curves[i],
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [resistor()];
    let links = [ConductionLink {
        link: 19,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    }];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let bases = [(1, ElectricalPropertyBasis::Resistivity)];
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &bases,
    };
    let contact_card = interface(&cards, resistance);
    let contact_context = QueryPoint::new()
        .with_quantity("aT", temperature_quantity(), INITIAL[0])
        .unwrap()
        .with_quantity("bT", temperature_quantity(), INITIAL[1])
        .unwrap();
    let contacts = [MaterialContactLink {
        link: 19,
        sides: [0, 1].map(|i| ContactCell {
            cell: IDS[i],
            material_id: MIDS[i],
            curve: &curves[i],
        }),
        point: ThermalInterfacePoint {
            interface: &contact_card,
            materials: [&cards[0], &cards[1]],
            context: &contact_context,
            temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let (complex, positions) = box_grid([1, 1, 1], [0.1; 3]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let surfaces = [0, 1]
        .into_iter()
        .map(|i| {
            let phase = curves[i]
                .state_at_specific_enthalpy((INITIAL[i] - 280.0) * CAPACITY / DENSITY)
                .unwrap();
            let emissivity = SurfaceEmissivity::from_thermal_state(
                "fixture",
                cells[i].point,
                phase,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
            .0;
            RadiationSurface::new(
                &mesh,
                format!("body-{i}"),
                |face| on_box_face(face.centroid[0], i as f64 * 0.1),
                emissivity,
            )
            .unwrap()
        })
        .collect();
    // Analytic parallel-plate exchange; mesh faces supply area only, not a claim
    // that these are finite-cube geometric view factors.
    let enclosure = GrayDiffuseEnclosure::new(
        surfaces,
        ViewFactorMatrix::infinite_parallel_plates(RADIATING_AREA).unwrap(),
    )
    .unwrap();
    let radiating = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: IDS[i],
        material_id: MIDS[i],
        curve: &curves[i],
        material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &radiating,
        power_tolerance_w: 1e-10,
        radiosity_tolerance_w_m2: 1e-10,
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: if extra_paths { &contacts } else { &[] },
        radiation: extra_paths.then_some(&radiation),
    };
    let mut events: Vec<_> = [0, 1]
        .into_iter()
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MIDS[i],
                mass_kg: DENSITY,
                specific_h_j_kg: (INITIAL[i] - 280.0) * CAPACITY / DENSITY,
                geometry: CellGeometry {
                    start_m: [i as f64 * LENGTH, 0.0, 0.0],
                    end_m: [(i + 1) as f64 * LENGTH, 0.0, 0.0],
                    volume_m3: LENGTH * AREA,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 19,
            a: CELL,
            b: SUPPORT,
            area_m2: CONTACT_AREA,
            conductance_w_k: 999.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MIDS[0], &curves[0] as &dyn ThermalLaw),
            (MIDS[1], &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let body = model(&cards[0], &bulk, &thermal, OLD[0], OLD[1], OLD[2]);
    let binding = thermal_binding(&cards[0], [280.0, 360.0]);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 81,
                kernel_id: 53,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &assembly, thermal, body, &binding)
    })
}

// Direct series resistance and two diffuse parallel-surface Stefan-Boltzmann
// equations. No native source report or resolved coefficient defines the oracle.
fn transfer(t: [f64; 2], e: f64, resistance: f64) -> (f64, f64, f64) {
    let k = t.map(|t| 0.2 + 0.001 * (t - 300.0));
    let contact = resistance * (1.0 + 0.005 * (t[0] - 300.0));
    let g = CONTACT_AREA / (0.5 / k[0] + 0.5 / k[1] + contact);
    let eps = [
        e * (1.0 + 0.002 * (t[0] - 300.0)),
        0.6 * (1.0 + 0.002 * (t[1] - 300.0)),
    ];
    let radiation = RADIATING_AREA * 5.670_374_419e-8 * (t[1].powi(4) - t[0].powi(4))
        / (1.0 / eps[0] + 1.0 / eps[1] - 1.0);
    (g * (t[1] - t[0]) + radiation, g, radiation)
}
fn expected(e: f64, resistance: f64) -> ([f64; 2], [f64; 5]) {
    let support = |t| {
        root(280.0, 360.0, |s| {
            CAPACITY * (s - INITIAL[1]) + DT * transfer([t, s], e, resistance).0
        })
    };
    let before = internal(INITIAL[0], OLD[0], OLD[2], 0.2, SLOPES) + 0.5 * OLD[1].powi(2);
    let t = root(280.0, 360.0, |t| {
        let [x, v, q, _, _] = motion(t, OLD, 0.2, SLOPES, FORCE, DT);
        CAPACITY * (t - INITIAL[0]) + internal(t, x, q, 0.2, SLOPES) + 0.5 * v * v
            - before
            - FORCE * (x - OLD[0])
            - DT * transfer([t, support(t)], e, resistance).0
    });
    ([t, support(t)], motion(t, OLD, 0.2, SLOPES, FORCE, DT))
}

#[test]
fn g1_g3_contact_and_emissivity_sources_reach_ceramic_current_motion_and_energy() {
    let mut endpoints = Vec::new();
    for (e, resistance) in [(0.2, 0.5), (0.8, 0.5), (0.2, 8.0)] {
        let (temperatures, motion) = expected(e, resistance);
        let endpoint = with_case(
            e,
            resistance,
            true,
            |cx, _, assembly, thermal, body, binding| {
                let before_cells = thermal.cells().clone();
                let before_model = body.clone();
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let next = prepared
                    .step_piezoelectric_candidate(
                        cx,
                        &thermal,
                        DT,
                        &body,
                        &[],
                        FORCE,
                        binding,
                        thermal_controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let r = next.report();
                let sample = &r.coupled.source.network;
                for (i, id) in IDS.into_iter().enumerate() {
                    near(
                        next.thermal().temperature_k(id).unwrap(),
                        temperatures[i],
                        2e-8,
                    );
                }
                near(next.model().extension_m(), motion[0], 2e-10);
                near(next.model().velocity_m_s(), motion[1], 2e-9);
                near(next.model().charge_c(), motion[2], 2e-10);
                near(r.coupled.electrical.resistors[0].voltage_v, motion[3], 2e-9);
                near(r.coupled.deposited_joule_j, motion[4], 2e-10);
                let (power, g, radiation) = transfer(temperatures, e, resistance);
                near(sample.thermal_contacts[0].conductance_w_k, g, 2e-10);
                let photons = sample.radiation.as_ref().unwrap();
                near(photons.sources_w()[&CELL], radiation, 2e-9);
                near(photons.sources_w()[&SUPPORT], -radiation, 2e-9);
                near(
                    CAPACITY * (temperatures[1] - INITIAL[1]),
                    -DT * power,
                    2e-10,
                );
                near(r.external_heat_j, 0.0, 2e-10);
                near(r.closed_radiation_net_j, DT * photons.net_heat_w(), 0.0);
                near(r.total_energy_residual_j, 0.0, 2e-10);
                for surface in photons.surfaces() {
                    let bulk = sample
                        .material
                        .cells
                        .iter()
                        .find(|c| c.cell == surface.cell)
                        .unwrap();
                    assert_eq!(bulk.phase, surface.phase);
                    assert_eq!(
                        surface.material.card_identity(),
                        assembly
                            .electrical
                            .bulk
                            .cells
                            .iter()
                            .find(|c| c.cell == surface.cell)
                            .unwrap()
                            .point
                            .material
                            .content_hash()
                    );
                }
                let before = thermal.energy_j()
                    + internal(300.0, OLD[0], OLD[2], 0.2, SLOPES)
                    + 0.5 * OLD[1].powi(2);
                let after = next.thermal().energy_j()
                    + internal(temperatures[0], motion[0], motion[2], 0.2, SLOPES)
                    + 0.5 * motion[1].powi(2);
                near(after - before, FORCE * (motion[0] - OLD[0]), 2e-9);
                assert_eq!(thermal.cells(), &before_cells);
                assert_eq!(body, before_model);
                [
                    temperatures[0],
                    temperatures[1],
                    next.model().extension_m(),
                    next.model().charge_c(),
                ]
            },
        );
        endpoints.push(endpoint);
    }
    for changed in &endpoints[1..] {
        for (original, changed) in endpoints[0].iter().zip(changed) {
            assert!((original - changed).abs() > 1e-8);
        }
    }
}

#[test]
fn g0_no_extra_paths_recover_original_owner_and_cancellation_preserves_histories() {
    with_case(
        0.2,
        0.5,
        false,
        |cx, gate, assembly, thermal, body, binding| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let next = prepared
                .step_piezoelectric_candidate(
                    cx,
                    &thermal,
                    DT,
                    &body,
                    &[],
                    FORCE,
                    binding,
                    thermal_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let original = heated::step_candidate(
                assembly.electrical.bulk,
                &thermal,
                DT,
                &body,
                &[],
                FORCE,
                binding,
                assembly.electrical.electrical_bases,
                thermal_controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
            assert_eq!(next.model(), original.model());
            assert_eq!(next.thermal().cells(), original.thermal().cells());
            assert_eq!(next.report().closed_radiation_net_j, 0.0);
            let before = thermal.cells().clone();
            let before_model = body.clone();
            gate.request();
            assert!(matches!(
                prepared.step_piezoelectric_candidate(
                    cx,
                    &thermal,
                    DT,
                    &body,
                    &[],
                    FORCE,
                    binding,
                    thermal_controls(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.cells(), &before);
            assert_eq!(body, before_model);
            assert_eq!(thermal.time_s(), 0.0);
        },
    );
}
