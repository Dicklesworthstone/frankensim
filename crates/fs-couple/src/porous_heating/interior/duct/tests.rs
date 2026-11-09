//! Native molecular losses, original caloric cells and accepted thermal work.
//! Synthetic reservoirs/local wall isolate mechanisms; these are not empirical
//! material claims. No replacement acoustic or thermal solver is used.

use super::*;
use fs_duct::molecular::{MolecularCompliantWall, MolecularWallResolution};
use fs_duct::network::field::{NetworkFieldBudget, NetworkFieldProbe, NetworkFieldSection};
use fs_duct::network::molecular::{
    ConeAbsorptionControls, MolecularEdge, MolecularNetwork, MolecularNetworkFieldRequest,
};
use fs_duct::network::{
    DuctNetwork, JunctionCorrections, NetworkEdge, NetworkMedia, NetworkNode, NodeBoundary,
    PortExcitation, PortSource,
};
use fs_duct::{ConeSlices, Duct, Segment, WallLaw};
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::input::GasRecord;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, MaterialId, Reservoir, ThermalLaw,
    ThermalLink, TimedEvent,
};
use fs_thermal::observers::CoupledObservers;
use fs_thermochem::SpeciesId;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

const IDS: [CellId; 2] = [CellId(91), CellId(7)];
const LENGTHS: [f64; 2] = [0.03, 0.05];
const OMEGA: f64 = 2800.0;

fn cv(gas: GasState) -> f64 {
    gas.specific_heat_cp - gas.specific_gas_constant
}

struct Surface(AtomicUsize);
impl WallLaw for Surface {
    fn specific_impedance(&self, radius: f64, omega: f64) -> Result<C64, DuctError> {
        self.0.fetch_add(1, Ordering::Relaxed);
        Ok(C64::new(
            500.0,
            3e6 * (0.002 / radius).powi(2) / omega - 0.02 * omega,
        ))
    }
}

fn network() -> (MolecularNetwork, Arc<Surface>, GasState) {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let id = SpeciesId::new("synthetic").unwrap();
    let state = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes: (0..2)
            .map(|j| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * j as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (2e-5 + j as f64 * 3e-4),
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap();
    let graph = DuctNetwork {
        nodes: vec![
            NetworkNode {
                id: "left".into(),
                boundary: NodeBoundary::Port {
                    label: "left".into(),
                },
            },
            NetworkNode {
                id: "joint".into(),
                boundary: NodeBoundary::Junction,
            },
            NetworkNode {
                id: "right".into(),
                boundary: NodeBoundary::Port {
                    label: "right".into(),
                },
            },
        ],
        edges: [
            (
                "tube",
                0,
                1,
                Segment::Cylinder {
                    radius: 0.0015,
                    length: LENGTHS[0],
                },
            ),
            (
                "cone",
                1,
                2,
                Segment::Cone {
                    inlet_radius: 0.0015,
                    outlet_radius: 0.0025,
                    length: LENGTHS[1],
                },
            ),
        ]
        .into_iter()
        .map(|(id, from, to, segment)| NetworkEdge {
            id: id.into(),
            from,
            to,
            duct: Duct {
                segments: vec![segment],
            },
            corrections: JunctionCorrections::Disabled,
        })
        .collect(),
    };
    let media = NetworkMedia::new(
        &graph,
        graph
            .edges
            .iter()
            .map(|edge| (edge.id.clone(), GasRecord::Legacy(gas)))
            .collect(),
        GasRecord::Legacy(gas),
    )
    .unwrap();
    let edges = graph
        .edges
        .iter()
        .map(|edge| MolecularEdge {
            edge_id: edge.id.clone(),
            state: state.clone(),
            cones: ConeSlices::Fixed(1),
            wall: Some(MolecularWallResolution {
                poles_per_field: 8,
                band_rad_s: [0.1, 7000.0],
                maximum_relative_error: 0.02,
            }),
        })
        .collect();
    let wall = Arc::new(Surface(AtomicUsize::new(0)));
    let walls = graph
        .edges
        .iter()
        .map(|edge| (edge.id.clone(), wall.clone() as MolecularCompliantWall))
        .collect();
    (
        MolecularNetwork::new_with_compliant_walls(graph, media, edges, walls, 64).unwrap(),
        wall,
        gas,
    )
}

fn sources() -> [PortSource; 2] {
    [
        ("left", C64::new(80.0, 20.0)),
        ("right", C64::new(-40.0, 60.0)),
    ]
    .map(|(label, value)| PortSource {
        label: label.into(),
        excitation: PortExcitation::Pressure {
            open_circuit_pressure_pa: value,
            impedance_pa_s_m3: C64::from_re(2e7),
        },
    })
}
fn probes() -> [NetworkFieldProbe; 4] {
    ["tube", "cone"]
        .into_iter()
        .flat_map(|edge| {
            [
                NetworkFieldProbe {
                    edge_id: edge.into(),
                    section: NetworkFieldSection::Inlet { segment: 0 },
                },
                NetworkFieldProbe {
                    edge_id: edge.into(),
                    section: NetworkFieldSection::Outlet { segment: 0 },
                },
            ]
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}
fn request<'a>(
    sources: &'a [PortSource],
    probes: &'a [NetworkFieldProbe],
) -> MolecularNetworkFieldRequest<'a> {
    MolecularNetworkFieldRequest {
        omega: OMEGA,
        source_clock: "native-molecular-thermal-clock",
        sources,
        source_bank: None,
        loads: &[],
        coupled_load: None,
        probes,
        budget: NetworkFieldBudget::default(),
    }
}
fn volumes() -> [f64; 2] {
    [
        core::f64::consts::PI * 0.0015_f64.powi(2) * LENGTHS[0],
        core::f64::consts::PI
            * LENGTHS[1]
            * (0.0015_f64.powi(2) + 0.0015 * 0.0025 + 0.0025_f64.powi(2))
            / 3.0,
    ]
}

#[derive(Clone, Copy)]
struct Caloric {
    gas: GasState,
    basis: ThermalStorageBasis,
    identity: Option<[u8; 32]>,
}
impl Caloric {
    fn new(gas: GasState) -> Self {
        Self {
            gas,
            basis: ThermalStorageBasis::InternalEnergy,
            identity: Some([41; 32]),
        }
    }
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.identity
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        let temperature = self.gas.temperature + u / cv(self.gas);
        if !(250.0..=400.0).contains(&temperature) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(temperature)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / cv(self.gas)
    }
}
fn boundary(gas: GasState, source_w: f64) -> Boundary {
    Boundary {
        ambient_k: gas.temperature,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn schedule(gas: GasState, sources: [f64; 2]) -> Vec<TimedEvent> {
    let volume = volumes();
    (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: 10 + i as u64,
                body_id: 2,
                material_id: MaterialId(5),
                mass_kg: gas.density * volume[i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [if i == 0 { 0.0 } else { LENGTHS[0] }, 0.0, 0.0],
                    end_m: [
                        if i == 0 {
                            LENGTHS[0]
                        } else {
                            LENGTHS[0] + LENGTHS[1]
                        },
                        0.0,
                        0.0,
                    ],
                    volume_m3: volume[i],
                    surface_area_m2: 0.002,
                },
                boundary: boundary(gas, sources[i]),
            }),
        })
        .collect()
}
fn initial<'t>(events: Vec<TimedEvent>, law: &'t Caloric) -> DepositionSimulation<'t> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(5), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 2048,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn routes(law: &Caloric) -> [InteriorHeatCell<'_>; 2] {
    IDS.map(|cell| InteriorHeatCell { cell, caloric: law })
}
fn config(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-12,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * expected.abs().max(1e-12),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn check_energy(energy: DuctHeatingEnergy, powers: &[MolecularNetworkCylinderAbsorption], dt: f64) {
    let expected: [f64; 4] =
        std::array::from_fn(|i| dt * powers.iter().map(|p| channels(*p)[i]).sum::<f64>());
    for (actual, expected) in [
        energy.viscous_heat_j,
        energy.thermal_availability_loss_j,
        energy.molecular_availability_loss_j,
        energy.wall_work_j,
    ]
    .into_iter()
    .zip(expected)
    {
        close(actual, expected);
    }
    close(energy.deposited_j, expected[0]);
}

#[test]
fn mixed_coherent_cells_heat_original_internal_energy_once_across_events() {
    let (network, wall, gas) = network();
    let (sources, points) = (sources(), probes());
    let fields = network
        .prepare_fields(&request(&sources, &points), || false)
        .unwrap();
    let cylinders = fields
        .prepare_cylinder_absorption(&[[0, 1]], 1, || false)
        .unwrap();
    let cones = fields
        .prepare_cone_absorption(&[[2, 3]], ConeAbsorptionControls::default(), || false)
        .unwrap();
    let mean = fields.fields().port_response().default_excitations();
    let factor = [C64::new(30.0, -10.0), C64::new(20.0, 5.0)];
    let factors = [&factor[..]];
    let expected: Vec<_> = cylinders
        .observe_coherence(Some(mean), &factors, || false)
        .unwrap()
        .into_iter()
        .chain(
            cones
                .observe_coherence(Some(mean), &factors, || false)
                .unwrap(),
        )
        .collect();
    let capacity = volumes().map(|v| gas.density * v * cv(gas));
    let law = Caloric::new(gas);
    let mut events = schedule(gas, [0.1 * capacity[0], -0.03 * capacity[1]]);
    events.push(TimedEvent {
        time_s: 0.375,
        event: Event::Boundary(IDS[0], boundary(gas, 0.2 * capacity[0])),
    });
    let thermal = initial(events, &law);
    let heating = MolecularDuctViscousHeating::try_new(
        Some(&cylinders),
        Some(&cones),
        Some(mean),
        &factors,
        &thermal,
        &routes(&law),
        2,
        2.0,
        || Ok(()),
    )
    .unwrap();
    assert!(core::ptr::eq(heating.fields(), &fields));
    assert!(core::ptr::eq(heating.cylinders().unwrap(), &cylinders));
    assert!(core::ptr::eq(heating.cones().unwrap(), &cones));
    assert_eq!(heating.mean(), Some(mean));
    assert_eq!(heating.factors(), factors);
    assert_eq!(wall.0.load(Ordering::Relaxed), 2);
    let accepted = heating
        .advance_to(&thermal, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let authored = [
        (0.1 * 0.375 + 0.2 * 0.625) * capacity[0],
        -0.03 * capacity[1],
    ];
    for (i, cell) in heating.cells().enumerate() {
        assert_eq!(cell.cell, IDS[i]);
        assert_eq!(cell.power, expected[i]);
        close(cell.mass_kg, thermal.cells()[&IDS[i]].mass_kg);
        let original = &thermal.cells()[&IDS[i]];
        let next = &accepted.thermal.cells()[&IDS[i]];
        close(
            next.mass_kg * (next.specific_h_j_kg - original.specific_h_j_kg),
            authored[i] + expected[i].viscous_power_w,
        );
        assert_eq!(next.geometry, original.geometry);
        assert_eq!(accepted.cell_energy_j[i].0, IDS[i]);
        close(accepted.cell_energy_j[i].1, expected[i].viscous_power_w);
    }
    check_energy(accepted.acoustic, &expected, 1.0);
    assert!(accepted.acoustic.thermal_availability_loss_j > 0.0);
    assert!(accepted.acoustic.molecular_availability_loss_j > 0.0);
    assert!(accepted.acoustic.wall_work_j > 0.0);
    close(
        accepted.heat_report.external_heat_j,
        authored.iter().sum::<f64>() + accepted.acoustic.deposited_j,
    );
    assert!(accepted.heat_report.residual_j.abs() < 1e-14);
    assert_eq!(accepted.heat_report.events, 1);
    assert_eq!(wall.0.load(Ordering::Relaxed), 2);
}

#[test]
fn native_transport_composition_split_and_late_cancellation_preserve_both_accounts() {
    let (network, wall, gas) = network();
    let (sources, points) = (sources(), probes());
    let fields = network
        .prepare_fields(&request(&sources, &points), || false)
        .unwrap();
    let cylinders = fields
        .prepare_cylinder_absorption(&[[0, 1]], 1, || false)
        .unwrap();
    let cones = fields
        .prepare_cone_absorption(&[[2, 3]], ConeAbsorptionControls::default(), || false)
        .unwrap();
    let mean = fields.fields().port_response().default_excitations();
    let expected: Vec<_> = cylinders
        .sample(mean, || false)
        .unwrap()
        .into_iter()
        .chain(cones.sample(mean, || false).unwrap())
        .collect();
    let law = Caloric::new(gas);
    let capacity = volumes().map(|v| gas.density * v * cv(gas));
    let mut events = schedule(gas, [0.1 * capacity[0], 0.0]);
    if let Event::Birth(cell) = &mut events[1].event {
        cell.boundary.reservoir = Some(Reservoir {
            area_m2: 1e-5,
            conductance_w_k: capacity[1],
            temperature_k: gas.temperature - 0.2,
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 3,
            a: IDS[0],
            b: IDS[1],
            area_m2: 1e-5,
            conductance_w_k: 4.0 * capacity[0],
        }),
    });
    events.push(TimedEvent {
        time_s: 0.375,
        event: Event::Boundary(IDS[0], boundary(gas, 0.2 * capacity[0])),
    });
    let thermal = initial(events.clone(), &law);
    let heating = MolecularDuctViscousHeating::try_new(
        Some(&cylinders),
        Some(&cones),
        Some(mean),
        &[],
        &thermal,
        &routes(&law),
        2,
        2.0,
        || Ok(()),
    )
    .unwrap();
    let observers = CoupledObservers::new((), heating.observer());
    let (middle, middle_account, _) = thermal
        .advance_observed(&observers, 0.5, config(0.125), || Ok(()))
        .unwrap();
    let (actual, accepted, _) = middle
        .advance_observed(&middle_account, 1.0, config(0.125), || Ok(()))
        .unwrap();
    let mut polls = 0;
    let (whole, whole_account, whole_report) = thermal
        .advance_observed(&observers, 1.0, config(0.125), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    // Independent injection of native viscous watts on the SAME thermal owner,
    // including the authored source replacement event. No availability credit.
    for event in &mut events {
        match &mut event.event {
            Event::Birth(cell) => {
                let i = IDS.iter().position(|&id| id == cell.id).unwrap();
                cell.boundary.source_w += expected[i].viscous_power_w;
            }
            Event::Boundary(id, boundary) => {
                let i = IDS.iter().position(|candidate| candidate == id).unwrap();
                boundary.source_w += expected[i].viscous_power_w;
            }
            _ => {}
        }
    }
    let reference = initial(events, &law)
        .advance_to(1.0, config(0.125), || Ok(()))
        .unwrap();
    for id in IDS {
        close(
            actual.cells()[&id].specific_h_j_kg,
            reference.0.cells()[&id].specific_h_j_kg,
        );
        close(
            actual.cells()[&id].specific_h_j_kg,
            whole.cells()[&id].specific_h_j_kg,
        );
        assert_eq!(actual.cells()[&id].geometry, thermal.cells()[&id].geometry);
    }
    close(whole_report.external_heat_j, reference.1.external_heat_j);
    assert_eq!(whole_report.substeps, reference.1.substeps);
    assert_eq!(actual.links(), thermal.links());
    check_energy(accepted.second.energy().unwrap(), &expected, 1.0);
    check_energy(middle_account.second.energy().unwrap(), &expected, 0.5);
    assert_eq!(
        accepted.second.energy().unwrap(),
        whole_account.second.energy().unwrap()
    );
    let before = thermal.cells().clone();
    let mut attempted = 0;
    assert!(
        thermal
            .advance_observed(&observers, 1.0, config(0.125), || {
                attempted += 1;
                if attempted >= polls * 3 / 4 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &before);
    assert_eq!(
        observers.second.energy().unwrap(),
        DuctHeatingEnergy::default()
    );
    assert_eq!(middle.time_s(), 0.5);
    check_energy(middle_account.second.energy().unwrap(), &expected, 0.5);
    let (retry, retry_account, _) = thermal
        .advance_observed(&observers, 1.0, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), whole.cells());
    assert_eq!(
        retry_account.second.energy().unwrap(),
        whole_account.second.energy().unwrap()
    );
    assert_eq!(wall.0.load(Ordering::Relaxed), 2);
}

#[test]
fn physical_mapping_coherence_identity_and_temperature_excursion_are_required() {
    let (network, _, gas) = network();
    let (sources, points) = (sources(), probes());
    let fields = network
        .prepare_fields(&request(&sources, &points), || false)
        .unwrap();
    let cylinders = fields
        .prepare_cylinder_absorption(&[[0, 1]], 1, || false)
        .unwrap();
    let cones = fields
        .prepare_cone_absorption(&[[2, 3]], ConeAbsorptionControls::default(), || false)
        .unwrap();
    let mean = fields.fields().port_response().default_excitations();
    let law = Caloric::new(gas);
    let thermal = initial(schedule(gas, [0.0; 2]), &law);
    let make = |simulation: &DepositionSimulation<'_>, destinations: &[InteriorHeatCell<'_>]| {
        MolecularDuctViscousHeating::try_new(
            Some(&cylinders),
            Some(&cones),
            Some(mean),
            &[],
            simulation,
            destinations,
            2,
            2.0,
            || Ok(()),
        )
        .is_ok()
    };
    assert!(make(&thermal, &routes(&law)));
    let equal_but_foreign = law;
    assert!(!make(&thermal, &routes(&equal_but_foreign)));
    let duplicate = [routes(&law)[0]; 2];
    assert!(!make(&thermal, &duplicate));
    for missing in [true, false] {
        let mut invalid = law;
        if missing {
            invalid.identity = None;
        } else {
            invalid.basis = ThermalStorageBasis::Enthalpy;
        }
        let timeline = initial(schedule(gas, [0.0; 2]), &invalid);
        assert!(!make(&timeline, &routes(&invalid)));
    }
    for changed in 0..4 {
        let mut events = schedule(gas, [0.0; 2]);
        if let Event::Birth(cell) = &mut events[0].event {
            match changed {
                0 => cell.mass_kg *= 1.01,
                1 => cell.geometry.volume_m3 *= 1.01,
                2 => cell.geometry.end_m[0] *= 1.01,
                _ => cell.specific_h_j_kg = cv(gas),
            }
        }
        assert!(!make(&initial(events, &law), &routes(&law)));
    }
    let other_fields = network
        .prepare_fields(&request(&sources, &points), || false)
        .unwrap();
    let other_cones = other_fields
        .prepare_cone_absorption(&[[2, 3]], ConeAbsorptionControls::default(), || false)
        .unwrap();
    assert!(
        MolecularDuctViscousHeating::try_new(
            Some(&cylinders),
            Some(&other_cones),
            Some(mean),
            &[],
            &thermal,
            &routes(&law),
            2,
            2.0,
            || Ok(())
        )
        .is_err()
    );
    let overlap = fields
        .prepare_cylinder_absorption(&[[0, 1], [0, 1]], 2, || false)
        .unwrap();
    let mut same_volume = schedule(gas, [0.0; 2]);
    let first = match &same_volume[0].event {
        Event::Birth(cell) => cell.clone(),
        _ => unreachable!(),
    };
    if let Event::Birth(cell) = &mut same_volume[1].event {
        cell.mass_kg = first.mass_kg;
        cell.geometry = first.geometry;
    }
    let overlapping = initial(same_volume, &law);
    assert!(matches!(
        MolecularDuctViscousHeating::try_new(
            Some(&overlap),
            None,
            Some(mean),
            &[],
            &overlapping,
            &routes(&law),
            2,
            2.0,
            || Ok(())
        ),
        Err(DuctHeatingError::Binding(InteriorHeatingError::Invalid(
            "duct heat intervals overlap on the same original edge and segment"
        )))
    ));
    let zero = fields
        .prepare_cylinder_absorption(&[[0, 0]], 1, || false)
        .unwrap();
    assert!(
        MolecularDuctViscousHeating::try_new(
            Some(&zero),
            None,
            Some(mean),
            &[],
            &thermal,
            &routes(&law)[..1],
            1,
            2.0,
            || Ok(())
        )
        .is_err()
    );
    let power = cylinders.sample(mean, || false).unwrap()[0].viscous_power_w;
    let change = power / (thermal.cells()[&IDS[0]].mass_kg * cv(gas));
    assert!(change > 0.0);
    let strict = MolecularDuctViscousHeating::try_new(
        Some(&cylinders),
        None,
        Some(mean),
        &[],
        &thermal,
        &routes(&law)[..1],
        1,
        0.5 * change,
        || Ok(()),
    )
    .unwrap();
    let observer = strict.observer();
    let before = thermal.cells().clone();
    assert!(matches!(
        thermal.advance_observed(&observer, 1.0, config(0.25), || Ok(())),
        Err(DuctHeatingError::TemperatureExcursion { .. })
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(observer.energy().unwrap(), DuctHeatingEnergy::default());
    let relaxed = MolecularDuctViscousHeating::try_new(
        Some(&cylinders),
        None,
        Some(mean),
        &[],
        &thermal,
        &routes(&law)[..1],
        1,
        2.0 * change,
        || Ok(()),
    )
    .unwrap();
    close(
        relaxed
            .advance_to(&thermal, 1.0, config(0.25), || Ok(()))
            .unwrap()
            .acoustic
            .deposited_j,
        power,
    );
}
