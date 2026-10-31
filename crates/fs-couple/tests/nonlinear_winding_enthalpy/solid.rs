//! G1/G3/G4: resolved solid currents, spatial temperatures and the original
//! circuit/field/caloric histories. All material values are synthetic fixtures.

use super::*;
#[path = "../../../fs-em/tests/winding_ohmic/common.rs"]
mod linear;

use fs_couple::circuit_enthalpy::{
    contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    eddy::thermal::MagneticThermalLaw,
    material::{ConductionLink, MaterialNetwork},
};
use fs_em::{
    planar::{
        eddy::{EddyConductor, nonlinear::PlanarNonlinearEddy, solid::SolidWindingRegion},
        nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
    },
    winding::nonlinear::{
        NonlinearWindingCircuit, NonlinearWindingInitial, NonlinearWindingStep, NonlinearWindings,
    },
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{InterpolationPolicy, PropertyClaim, PropertyValue, Provenance, UncertaintyModel};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::thermal::magnetic::ThermalMagneticMaterial,
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_rep_mesh::{Metric2, TriComplex2};

const STEP: f64 = 0.1;
const VOLTAGE: f64 = 1.0;
const END_H: f64 = 0.05;
const CELL_VOLUME: f64 = 0.5 * DEPTH;
const HEAT_CAPACITY: f64 = 4.0 * CELL_VOLUME;
const HEATER_W: f64 = 6.0;
const SOLID_MATERIAL: MaterialId = MaterialId(30);
const CELL_IDS: [CellId; 4] = [CellId(20), CellId(21), CellId(22), CellId(23)];
const TURN_COUNTS: [f64; 2] = [6.0, -6.0];

fn config() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 4,
        ..controls()
    }
}

fn source(slope: f64) -> MaterialCard {
    let base = if slope == 0.0 {
        linear::card(1.0 / 18.0, 0.005, 1.0, ElectricalPropertyBasis::Resistivity)
    } else {
        linear::absolute_affine_magnetic_card(
            1.0 / 18.0,
            0.005,
            1.0,
            slope,
            ElectricalPropertyBasis::Resistivity,
        )
    };
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 1.0,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: base
                .claims()
                .claims_ordered()
                .next()
                .unwrap()
                .1
                .validity
                .clone(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic spatial solid winding heat fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn mesh() -> TriComplex2 {
    // Two unit-area conductor sections separated by a unit-width insulating
    // gap. Multiple free nodes make redistribution change the field load.
    let vertices = (0..=6)
        .flat_map(|x| (0..=2).map(move |y| [0.5 * x as f64, 0.5 * y as f64]))
        .collect();
    let mut faces = Vec::new();
    for x in 0..6 {
        for y in 0..2 {
            let a = (3 * x + y) as u32;
            faces.extend([[a, a + 3, a + 4], [a, a + 4, a + 1]]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "spatial-solid-winding-thermal-cells",
        vertices,
        faces,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap()
}

fn receiver(mesh: &TriComplex2, triangle: usize) -> Option<usize> {
    let center: [f64; 2] = [0, 1].map(|axis| {
        mesh.faces()[triangle]
            .iter()
            .map(|&v| mesh.vertices()[v as usize][axis] / 3.0)
            .sum()
    });
    let section = if center[0] < 1.0 {
        0
    } else if center[0] > 2.0 {
        1
    } else {
        return None;
    };
    Some(2 * section + usize::from(center[1] > 0.5))
}

fn circuit(problem: &PlanarNonlinearMagnetostatics) -> NonlinearWindingCircuit<'_> {
    let mesh = problem.mesh();
    let boundary: Vec<_> = mesh
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0 || p[0] == 3.0 || p[1] == 0.0 || p[1] == 1.0)
        .map(|(vertex, _)| (vertex, 0.0))
        .collect();
    let mut conductors = Vec::new();
    let mut insulators = Vec::new();
    let mut profile = vec![0.0; mesh.faces().len()];
    for (triangle, law) in problem.materials().iter().enumerate() {
        if let Some(cell) = receiver(mesh, triangle) {
            profile[triangle] = TURN_COUNTS[cell / 2];
            conductors.push(EddyConductor {
                triangle,
                region: if cell < 2 { 10 } else { 20 },
                conductor: ResolvedConductorState::try_from_material_with_basis(
                    law.material().clone(),
                    ElectricalPropertyBasis::Resistivity,
                )
                .unwrap(),
            });
        } else {
            insulators.push(triangle);
        }
    }
    let zero = vec![0.0; mesh.faces().len()];
    let initial = problem
        .solve(
            &zero,
            &boundary,
            &vec![0.0; mesh.vertices().len()],
            field::material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let field = PlanarNonlinearEddy::new_with_impressed_current(
        &initial,
        conductors,
        &insulators,
        &boundary,
        &zero,
        0.0,
        field::material::controls(),
        &mut || false,
    )
    .unwrap();
    NonlinearWindingCircuit::new_with_solid_regions(
        field,
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 123.0 }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profile],
            end_inductance_h: vec![END_H],
        },
        vec![],
        [10, 20]
            .into_iter()
            .enumerate()
            .map(|(section, region)| SolidWindingRegion {
                region,
                winding: 0,
                signed_turns: TURN_COUNTS[section],
            })
            .collect(),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![VOLTAGE],
        },
        field::config(),
        &mut || false,
    )
    .unwrap()
}

fn with_system<R>(
    slope: f64,
    linked: bool,
    run: impl FnOnce(
        &MaterialCard,
        DepositionSimulation<'_>,
        &PlanarNonlinearMagnetostatics,
        &[MaterialCell<'_>],
        &[MagneticRegion],
        &[ConductionLink],
    ) -> R,
) -> R {
    let card = source(slope);
    let air = source(0.0);
    let context = linear::absolute_point(300.0);
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &linear::absolute_point(200.0),
            upper: &linear::absolute_point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap();
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = CELL_IDS.map(|cell| MaterialCell {
        cell,
        material_id: SOLID_MATERIAL,
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
    });
    let mut events: Vec<_> = CELL_IDS
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: SOLID_MATERIAL,
                mass_kg: CELL_VOLUME,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [2.0 * (i / 2) as f64, 0.5 * (i % 2) as f64, 0.0],
                    end_m: [2.0 * (i / 2) as f64 + 1.0, 0.5 * (i % 2) as f64, 0.0],
                    volume_m3: CELL_VOLUME,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: if i == 0 { HEATER_W } else { 0.0 },
                },
            }),
        })
        .collect();
    let links = if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 17,
                a: CELL_IDS[0],
                b: CELL_IDS[1],
                area_m2: DEPTH,
                conductance_w_k: 999.0,
            }),
        });
        vec![ConductionLink {
            link: 17,
            distance_a_m: 0.25,
            distance_b_m: 0.25,
        }]
    } else {
        vec![]
    };
    let thermal = DepositionSimulation::new(
        events,
        &[(SOLID_MATERIAL, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 1,
            max_events: 5,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: STEP,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let mesh = mesh();
    let regions: Vec<_> = (0..mesh.faces().len())
        .map(|triangle| match receiver(&mesh, triangle) {
            Some(cell) => MagneticRegion::Cell {
                triangle,
                cell: CELL_IDS[cell],
            },
            None => MagneticRegion::Fixed { triangle },
        })
        .collect();
    let materials = (0..mesh.faces().len())
        .map(|triangle| {
            IsotropicMagneticMaterial::Linear(linear::bundle_at(
                if receiver(&mesh, triangle).is_some() {
                    &card
                } else {
                    &air
                },
                &context,
                ElectricalPropertyBasis::Resistivity,
            ))
        })
        .collect();
    let problem =
        PlanarNonlinearMagnetostatics::new(mesh, materials, linear::limits(), &mut || false)
            .unwrap();
    run(&card, thermal, &problem, &cells, &regions, &links)
}

fn temperatures(thermal: &DepositionSimulation<'_>) -> [f64; 4] {
    CELL_IDS.map(|cell| thermal.temperature_k(cell).unwrap())
}

// Independent exact P1 quadrature and constitutive relations. This uses the
// returned A only as a field to check, not the driver's recovered E/J or heat.
fn check_receipt(
    card: &MaterialCard,
    old: &NonlinearWindingCircuit<'_>,
    next: &NonlinearWindingCircuit<'_>,
    temperature: [f64; 4],
    step: &NonlinearWindingStep,
) -> [f64; 4] {
    let mesh = next.field().problem().mesh();
    let delta: Vec<_> = next
        .field()
        .potential_wb_m()
        .iter()
        .zip(old.field().potential_wb_m())
        .map(|(a, b)| a - b)
        .collect();
    let sigma = temperature.map(|t| 18.0 / (1.0 + 0.005 * (t - 300.0)));
    let mut weight = [0.0; 2];
    let mut mean_delta = [0.0; 2];
    for triangle in 0..mesh.faces().len() {
        if let Some(cell) = receiver(mesh, triangle) {
            let volume = mesh.face_measure(triangle).unwrap();
            weight[cell / 2] += sigma[cell] * volume;
            mean_delta[cell / 2] += sigma[cell]
                * volume
                * mesh.faces()[triangle]
                    .iter()
                    .map(|&v| delta[v as usize] / 3.0)
                    .sum::<f64>();
        }
    }
    for region in 0..2 {
        mean_delta[region] /= weight[region];
    }
    let current = next.currents_a()[0];
    let mut cell_heat = [0.0; 4];
    let mut net_current = [0.0; 2];
    let mut pairing = 0.0;
    let mut profile_motion = 0.0;
    assert_eq!(step.field.cells.len(), 16);
    assert!(step.cells.is_empty());
    assert!(step.graph_resistors.is_empty());
    for loss in &step.field.cells {
        let triangle = loss.source.triangle;
        let cell = receiver(mesh, triangle).unwrap();
        let region = cell / 2;
        let volume = mesh.face_measure(triangle).unwrap();
        let face = mesh.faces()[triangle];
        let unit_drive = TURN_COUNTS[region] * DEPTH / weight[region];
        let profile = sigma[cell] * unit_drive;
        let electric =
            face.map(|v| unit_drive * current - (delta[v as usize] - mean_delta[region]) / STEP);
        let heat = STEP * sigma[cell] * volume / 12.0
            * (electric.iter().sum::<f64>().powi(2) + electric.iter().map(|e| e * e).sum::<f64>());
        near(next.profiles_per_m2()[0][triangle], profile);
        near(loss.source.conductor.resistivity_ohm_m(), 1.0 / sigma[cell]);
        assert_eq!(
            loss.source.conductor.material().card_identity(),
            card.content_hash()
        );
        let source_t = loss
            .source
            .conductor
            .material()
            .query_point()
            .iter()
            .find(|(name, _)| name == "T")
            .unwrap()
            .1;
        assert!((source_t - temperature[cell]).abs() < 1e-10);
        for vertex in 0..3 {
            near(loss.electric_field_v_m[vertex], electric[vertex]);
            near(
                loss.current_density_a_m2[vertex],
                sigma[cell] * electric[vertex],
            );
        }
        near(loss.volume_m3, volume);
        near(loss.heat_j, heat);
        near(loss.power_w * STEP, heat);
        cell_heat[cell] += heat;
        net_current[region] += sigma[cell] * electric.iter().sum::<f64>() / 3.0 * volume / DEPTH;
        pairing += profile * volume * face.iter().map(|&v| delta[v as usize] / 3.0).sum::<f64>();
        profile_motion += (profile - old.profiles_per_m2()[0][triangle])
            * volume
            * face
                .iter()
                .map(|&v| old.field().potential_wb_m()[v as usize] / 3.0)
                .sum::<f64>();
    }
    let resistance: f64 = (0..2)
        .map(|region| (TURN_COUNTS[region] * DEPTH).powi(2) / weight[region])
        .sum();
    near(next.resistance_ohm()[0], resistance);
    for (region, &turns) in TURN_COUNTS.iter().enumerate() {
        near(net_current[region], turns * current);
        near(step.field.region_currents_a[region].1, turns * current);
    }
    let end_increment = END_H * (current - old.currents_a()[0]);
    assert!((step.inductive_voltage_v[0] - (pairing + end_increment) / STEP).abs() < 1e-9);
    near(
        step.terminal_voltage_v[0],
        resistance * current + (pairing + end_increment) / STEP,
    );
    near(step.terminal_voltage_v[0], VOLTAGE);
    near(
        next.flux_linkages_wb()[0] - old.flux_linkages_wb()[0],
        pairing + end_increment + profile_motion,
    );
    if old.time_s() > 0.0 {
        // A moving conductivity profile must not turn the old retained A into
        // an artificial voltage. This quantity is nonzero in this geometry.
        assert!(profile_motion.abs() > 1e-9);
    }
    near(
        step.field.uniform_drive_work_j,
        STEP * resistance * current * current,
    );
    near(
        step.field.supplied_work_j,
        current * pairing + step.field.uniform_drive_work_j,
    );
    near(step.field.heat_j, cell_heat.iter().sum());
    near(step.heat_j, cell_heat.iter().sum());
    assert!(step.field.heat_j > step.field.uniform_drive_work_j + 1e-6);
    near(step.supplied_work_j, STEP * VOLTAGE * current);
    cell_heat
}

#[test]
fn spatial_solid_heat_and_profiles_follow_the_same_finite_temperature_endpoint() {
    with_system(
        0.0,
        false,
        |card, mut thermal, problem, cells, regions, links| {
            let layout = WindingThermalNetwork {
                cells,
                links,
                magnetic_regions: regions,
                graph_resistors: &[],
                volume_tolerance_m3: 1e-12,
            };
            let mut circuit = circuit(problem);
            for interval in 1..=2 {
                let before = temperatures(&thermal);
                let old_a = circuit.field().potential_wb_m().to_vec();
                let old_state = circuit.state().to_vec();
                let target = interval as f64 * STEP;
                let next = layout
                    .step_nonlinear_candidate(
                        &thermal,
                        target,
                        &circuit,
                        &[VOLTAGE],
                        config(),
                        field::controls(),
                        &[0.0; 4],
                        || false,
                    )
                    .unwrap();
                let t = temperatures(next.thermal());
                let report = next.report();
                let heat = check_receipt(card, &circuit, next.circuit(), t, &report.electrical);
                for cell in 0..4 {
                    near(
                        HEAT_CAPACITY * (t[cell] - before[cell]),
                        heat[cell] + if cell == 0 { STEP * HEATER_W } else { 0.0 },
                    );
                }
                assert!(t[0] > t[1] + 0.5);
                near(report.deposited_joule_j, heat.iter().sum());
                near(report.non_electrical_heat_j, STEP * HEATER_W);
                near(
                    report.thermal.stored_energy_change_j,
                    heat.iter().sum::<f64>() + STEP * HEATER_W,
                );
                assert!(report.electrical.numerical_dissipation_j > 1e-5);
                assert!(report.total_energy_residual_j.abs() < config().total_energy_tolerance_j);
                assert!(report.maximum_cell_residual_j <= config().cell_energy_tolerance_j);
                assert!(report.circuit_trials > report.electrical.field_trials);
                assert_eq!(temperatures(&thermal), before);
                assert_eq!(circuit.field().potential_wb_m(), old_a);
                assert_eq!(circuit.state(), old_state);
                near(circuit.time_s(), target - STEP);
                let (accepted_thermal, accepted_circuit, _) = next.into_parts();
                thermal = accepted_thermal;
                circuit = accepted_circuit;
            }
        },
    );
}

#[test]
fn full_heat_assembly_preserves_solid_loss_and_magnetic_caloric_owners() {
    let mut outcomes = Vec::new();
    for (linked, slope) in [(false, 0.0), (true, 0.0), (true, 0.004)] {
        outcomes.push(with_system(
            slope,
            linked,
            |card, mut thermal, problem, cells, regions, links| {
                let bulk = MaterialNetwork {
                    cells,
                    resistors: &[],
                    links,
                };
                let electrical = ElectricalContactNetwork {
                    bulk: &bulk,
                    contacts: &[],
                    electrical_bases: &[],
                };
                let assembly = MaterialHeatNetwork {
                    electrical: &electrical,
                    thermal_contacts: &[],
                    radiation: None,
                };
                let laws: Vec<_> = problem
                    .materials()
                    .iter()
                    .enumerate()
                    .filter(|(triangle, _)| receiver(problem.mesh(), *triangle).is_some())
                    .map(|(triangle, original)| MagneticThermalLaw {
                        triangle,
                        law: ThermalMagneticMaterial::new(
                            card,
                            original.material(),
                            "T",
                            [290.0, 360.0],
                        )
                        .unwrap(),
                        maximum_flux_density_t: 4.0,
                    })
                    .collect();
                let gate = CancelGate::new();
                let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
                pool.scope(|arena| {
                    let cx = Cx::new(
                        &gate,
                        arena,
                        fs_exec::StreamKey {
                            seed: 43,
                            kernel_id: 15,
                            tile: 0,
                            iteration: 0,
                        },
                        fs_exec::Budget::INFINITE,
                        fs_exec::ExecMode::Deterministic,
                    );
                    let prepared = assembly.prepare(&cx, &thermal).unwrap();
                    let mut circuit = circuit(problem);
                    for interval in 1..=2 {
                        let before = temperatures(&thermal);
                        let next = prepared
                            .step_nonlinear_winding_candidate(
                                &cx,
                                &thermal,
                                interval as f64 * STEP,
                                &circuit,
                                &[VOLTAGE],
                                regions,
                                &[],
                                1e-12,
                                Some((&laws, linear::limits())),
                                config(),
                                field::controls(),
                                &[0.0; 4],
                            )
                            .unwrap();
                        let t = temperatures(next.thermal());
                        let report = next.report();
                        let heat =
                            check_receipt(card, &circuit, next.circuit(), t, &report.electrical);
                        let caloric = report.source.magnetic_caloric.as_ref().unwrap();
                        let balance = next.caloric_balance().unwrap();
                        let mut reversible = [0.0; 4];
                        let mut internal = [0.0; 2];
                        for magnetic in &report.electrical.field.magnetic_cells {
                            let owner = receiver(problem.mesh(), magnetic.triangle);
                            let ts = owner.map_or([300.0; 2], |cell| [before[cell], t[cell]]);
                            let local_slope = if owner.is_some() { slope } else { 0.0 };
                            let mu =
                                ts.map(|temperature| 1.0 + local_slope * (temperature - 300.0));
                            let b2 = magnetic.flux_density_t.map(|b| b[0] * b[0] + b[1] * b[1]);
                            let volume = problem.mesh().face_measure(magnetic.triangle).unwrap();
                            let entropy: [f64; 2] =
                                [0, 1].map(|i| 0.5 * volume * b2[i] * local_slope / mu[i].powi(2));
                            for i in 0..2 {
                                internal[i] += 0.5 * volume * b2[i] / mu[i] + ts[i] * entropy[i];
                            }
                            if let Some(cell) = owner {
                                // The nonlinear backward-Euler owner changes
                                // the material at the retained old B first.
                                let parameter = 0.5 * volume * b2[0] * (1.0 / mu[1] - 1.0 / mu[0]);
                                let q = -parameter - (ts[1] * entropy[1] - ts[0] * entropy[0]);
                                let sample = caloric
                                    .fields
                                    .iter()
                                    .find(|s| s.triangle == magnetic.triangle)
                                    .unwrap();
                                near(sample.parameter_work_j, parameter);
                                near(sample.reversible_heat_j, q);
                                near(sample.temperature_k[1], t[cell]);
                                reversible[cell] += q;
                            }
                        }
                        let transfer = if linked {
                            STEP * 0.5 * (t[0] - t[1])
                        } else {
                            0.0
                        };
                        if linked {
                            let link = &report.source.network.material.links[0];
                            near(link.conductance_w_k, 0.5);
                            assert_eq!(link.cells, [CELL_IDS[0], CELL_IDS[1]]);
                        }
                        for cell in 0..4 {
                            let other = match cell {
                                0 => STEP * HEATER_W - transfer,
                                1 => transfer,
                                _ => 0.0,
                            };
                            near(
                                HEAT_CAPACITY * (t[cell] - before[cell]),
                                heat[cell] + reversible[cell] + other,
                            );
                        }
                        near(report.deposited_joule_j, heat.iter().sum());
                        near(caloric.reversible_heat_j, reversible.iter().sum());
                        near(caloric.field_internal_energy_j[0], internal[0]);
                        near(caloric.field_internal_energy_j[1], internal[1]);
                        near(balance.external_heat_j, STEP * HEATER_W);
                        assert!(
                            balance.total_energy_residual_j.abs()
                                < config().total_energy_tolerance_j
                        );
                        let delta_end = 0.5
                            * END_H
                            * (next.circuit().currents_a()[0].powi(2)
                                - circuit.currents_a()[0].powi(2));
                        let delta_caloric: f64 = t
                            .iter()
                            .zip(before)
                            .map(|(t, before)| HEAT_CAPACITY * (t - before))
                            .sum();
                        near(
                            internal[1] - internal[0]
                                + delta_end
                                + delta_caloric
                                + report.electrical.numerical_dissipation_j,
                            STEP * VOLTAGE * next.circuit().currents_a()[0] + STEP * HEATER_W,
                        );
                        if slope != 0.0 {
                            assert!(caloric.reversible_heat_j.abs() > 1e-7);
                            for source in &report.source.magnetic {
                                if let Some(cell) = source.cell {
                                    let source_t = source
                                        .material
                                        .query_point()
                                        .iter()
                                        .find(|(name, _)| name == "T")
                                        .unwrap()
                                        .1;
                                    near(source_t, next.thermal().temperature_k(cell).unwrap());
                                }
                            }
                        }
                        let (accepted_thermal, accepted_circuit, _) = next.into_parts();
                        thermal = accepted_thermal;
                        circuit = accepted_circuit;
                    }
                    (temperatures(&thermal), circuit.currents_a()[0])
                })
            },
        ));
    }
    assert!(outcomes[1].0[0] < outcomes[0].0[0] - 1e-3);
    assert!(outcomes[1].0[1] > outcomes[0].0[1] + 1e-3);
    assert!((outcomes[1].1 - outcomes[0].1).abs() > 1e-9);
    assert!((outcomes[2].1 - outcomes[1].1).abs() > 1e-7);
}

#[test]
fn actual_nested_field_budget_and_late_cancellation_preserve_solid_histories() {
    with_system(0.0, false, |_, thermal, problem, cells, regions, links| {
        let layout = WindingThermalNetwork {
            cells,
            links,
            magnetic_regions: regions,
            graph_resistors: &[],
            volume_tolerance_m3: 1e-12,
        };
        let circuit = circuit(problem);
        let original_a = circuit.field().potential_wb_m().to_vec();
        let original_state = circuit.state().to_vec();
        let original_profile = circuit.profiles_per_m2().to_vec();
        let original_t = temperatures(&thermal);
        let mut checkpoints = 0usize;
        let complete = layout
            .step_nonlinear_candidate(
                &thermal,
                STEP,
                &circuit,
                &[VOLTAGE],
                config(),
                field::controls(),
                &[0.0; 4],
                || {
                    checkpoints += 1;
                    false
                },
            )
            .unwrap();
        assert!(complete.report().circuit_trials > complete.report().electrical.field_trials);
        for maximum_circuit_trials in [1, complete.report().circuit_trials - 1] {
            let limited = CircuitEnthalpySolve {
                maximum_circuit_trials,
                ..config()
            };
            assert!(matches!(
                layout.step_nonlinear_candidate(
                    &thermal,
                    STEP,
                    &circuit,
                    &[VOLTAGE],
                    limited,
                    field::controls(),
                    &[0.0; 4],
                    || false,
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
        }
        let mut visited = 0usize;
        assert!(checkpoints > 2);
        assert!(matches!(
            layout.step_nonlinear_candidate(
                &thermal,
                STEP,
                &circuit,
                &[VOLTAGE],
                config(),
                field::controls(),
                &[0.0; 4],
                || {
                    visited += 1;
                    visited >= checkpoints - 1
                },
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(circuit.field().potential_wb_m(), original_a);
        assert_eq!(circuit.state(), original_state);
        assert_eq!(circuit.profiles_per_m2(), original_profile);
        assert_eq!(temperatures(&thermal), original_t);
        assert_eq!(circuit.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        let exact_budget = CircuitEnthalpySolve {
            maximum_circuit_trials: complete.report().circuit_trials,
            ..config()
        };
        let retry = layout
            .step_nonlinear_candidate(
                &thermal,
                STEP,
                &circuit,
                &[VOLTAGE],
                exact_budget,
                field::controls(),
                &[0.0; 4],
                || false,
            )
            .unwrap();
        assert_eq!(
            retry.circuit().field().potential_wb_m(),
            complete.circuit().field().potential_wb_m()
        );
        assert_eq!(retry.circuit().state(), complete.circuit().state());
        assert_eq!(
            temperatures(retry.thermal()),
            temperatures(complete.thermal())
        );
        assert_eq!(
            retry.report().circuit_trials,
            complete.report().circuit_trials
        );
    });
}
