//! G1/G3/G4: original three-dimensional terminal currents respond to finite heat.
//! Manufactured cards check equations and conservation, not measured devices.
#[path = "../../fs-em/tests/winding_ohmic/common.rs"]
mod source;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, ResistorCell,
    contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    magnetic::caloric::MagneticConstitutiveThermalLaw,
    material::{CellResistor, MaterialCell, MaterialNetwork},
    spatial_induction::{
        SpatialMagneticRegion, SpatialMagneticThermalLaw, SpatialMagneticThermalNetwork,
        ThermalMagneticMaterial,
        nonlinear::winding::{SpatialWindingThermalNetwork, heat_assembly_candidate},
    },
};
use fs_couple::material_circuit::MaterialCircuitSolve;
use fs_em::{
    spatial::{
        nonlinear::SpatialNonlinearMaterial,
        vector_potential::{
            VectorMagneticBoundary, VectorMagneticControls, VectorPotentialFace,
            induction::{
                SpatialInductionControls, SpatialInductionLimits,
                nonlinear::{
                    SpatialNonlinearInduction, SpatialNonlinearInductionControls,
                    winding::{
                        SpatialNonlinearWindingCircuit, SpatialNonlinearWindingConfig,
                        SpatialNonlinearWindings, SpatialWindingConductor,
                    },
                },
            },
            nonlinear::{NonlinearVectorControls, SpatialNonlinearVectorMagnetostatics},
        },
    },
    winding::nonlinear::{
        NonlinearWindingControls, NonlinearWindingInitial, resistors::GraphResistorMaterial,
    },
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{MaterialCard, PropertyKey, SelectionPolicy};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS},
    thermal::source::ThermalMaterialPoint,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::Density;
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};
use fs_thermal::enthalpy::*;

const BASIS: ElectricalPropertyBasis = ElectricalPropertyBasis::Resistivity;
const CELL: CellId = CellId(7);
const MAT: MaterialId = MaterialId(9);
const END_L: f64 = 0.03;
const GRAPH_R: f64 = 0.7;
const ROUTES: [ResistorCell; 1] = [ResistorCell {
    branch: 1,
    cell: CELL,
}];

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}
fn problem(card: &MaterialCard) -> SpatialNonlinearVectorMagnetostatics {
    let mut p = Vec::new();
    for k in 0..=2 {
        for j in 0..=2 {
            for i in 0..=2 {
                p.push([f64::from(i) / 2.0, f64::from(j) / 2.0, f64::from(k) / 2.0]);
            }
        }
    }
    let index = |i, j, k| ((k * 3 + j) * 3 + i) as u32;
    let mut tets = Vec::new();
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let (a, b, c, d) = (
                    index(i, j, k),
                    index(i + 1, j, k),
                    index(i, j + 1, k),
                    index(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    index(i, j, k + 1),
                    index(i + 1, j, k + 1),
                    index(i, j + 1, k + 1),
                    index(i + 1, j + 1, k + 1),
                );
                tets.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    let material = source::bundle_at(card, &source::absolute_point(300.0), BASIS);
    SpatialNonlinearVectorMagnetostatics::new(
        p,
        tets,
        vec![SpatialNonlinearMaterial::Linear(material); 48],
        WeightedAssemblyLimits {
            max_cells: 100,
            max_dofs: 2000,
            max_triplets: 100_000,
        },
        &mut || false,
    )
    .unwrap()
}
// Junit=2 curl(N_center e_z), a closed one-ampere loop with a complete return.
fn profile(p: &SpatialNonlinearVectorMagnetostatics) -> Vec<[f64; 3]> {
    let center = 13;
    p.mesh()
        .tets
        .iter()
        .map(|tet| {
            if !tet.contains(&center) {
                return [0.0; 3];
            }
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != center).collect();
            let origin = p.positions_m()[face[0] as usize];
            let normal = cross(
                sub(p.positions_m()[face[1] as usize], origin),
                sub(p.positions_m()[face[2] as usize], origin),
            );
            let distance = dot(normal, sub(p.positions_m()[center as usize], origin));
            cross(normal.map(|v| v / distance), [0.0, 0.0, 2.0])
        })
        .collect()
}
fn controls() -> NonlinearWindingControls<SpatialNonlinearInductionControls> {
    NonlinearWindingControls {
        field: SpatialNonlinearInductionControls {
            induction: SpatialInductionControls {
                relative_tolerance: 2e-11,
                absolute_residual_a: 1e-11,
                restart: 64,
                maximum_cycles: 100,
                maximum_krylov_scalars: 2_000_000,
                current_tolerance_a: 2e-8,
                flux_tolerance_wb: 2e-9,
                energy_tolerance_j: 2e-8,
                maximum_cell_magnetic_peclet: 100.0,
            },
            maximum_newton_iterations: 50,
            maximum_backtracks: 35,
        },
        maximum_field_trials: 150,
        maximum_newton_iterations: 25,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 2e-7,
        power_tolerance_w: 2e-7,
        energy_tolerance_j: 3e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn thermal_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 3000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 2e-7,
        thermal_energy_tolerance_j: 2e-7,
        total_energy_tolerance_j: 5e-7,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 3e-7,
            power_tolerance_w: 3e-7,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn circuit(
    p: &SpatialNonlinearVectorMagnetostatics,
    fill: f64,
) -> SpatialNonlinearWindingCircuit<'_> {
    let j = profile(p);
    let cg = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 3000,
        max_corrections: 2,
    };
    let initial = p
        .solve(
            &vec![[0.0; 3]; 48],
            &VectorMagneticBoundary {
                potential: p
                    .exterior_faces()
                    .into_iter()
                    .map(|face| VectorPotentialFace {
                        face,
                        edge_integrals_wb: [0.0; 3],
                    })
                    .collect(),
            },
            &vec![0.0; p.mesh().edges.len()],
            NonlinearVectorControls {
                field: VectorMagneticControls {
                    field: cg,
                    gauge: cg,
                    maximum_gauge_product_terms: 2_000_000,
                    current_tolerance_a: 2e-8,
                    flux_tolerance_wb: 2e-9,
                    energy_tolerance_j: 2e-8,
                },
                relative_residual_tolerance: 2e-11,
                maximum_newton_iterations: 50,
                maximum_backtracks: 35,
            },
            &mut || false,
        )
        .unwrap();
    let electrical = p
        .materials()
        .iter()
        .map(|m| {
            ResolvedConductorState::try_from_material_with_basis(m.material().clone(), BASIS)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let eddy = j
        .iter()
        .enumerate()
        .map(|(cell, &j)| (j == [0.0; 3]).then(|| electrical[cell].clone()))
        .collect();
    let sources = j
        .iter()
        .enumerate()
        .filter(|(_, j)| **j != [0.0; 3])
        .map(|(cell, _)| SpatialWindingConductor {
            winding: 0,
            cell,
            fill_fraction: fill,
            conductor: electrical[cell].clone(),
        })
        .collect();
    let field = SpatialNonlinearInduction::new(
        &initial,
        eddy,
        SpatialInductionLimits {
            maximum_unknowns: 20000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 1_000_000,
        },
        2e-9,
        &mut || false,
    )
    .unwrap();
    SpatialNonlinearWindingCircuit::new(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: GRAPH_R }),
                (2, 0, Branch::Inductor { henries: 7.0 }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![2],
            profiles_per_m2: vec![j],
            end_inductance_h: vec![END_L],
        },
        sources,
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![10.0],
        },
        SpatialNonlinearWindingConfig {
            max_windings: 2,
            max_profile_values: 100_000,
            current_balance_tolerance_turns: 1e-10,
        },
        &mut || false,
    )
    .unwrap()
}
fn curve(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &source::absolute_point(200.0),
            upper: &source::absolute_point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap()
}
fn thermal(curve: &ResolvedHeatCapacityCurve, cooling: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 0,
                body_id: 0,
                material_id: MAT,
                mass_kg: 1.0,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 290.0,
                    htc_w_m2_k: cooling,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MAT, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn regions() -> Vec<SpatialMagneticRegion> {
    (0..48)
        .map(|tetrahedron| SpatialMagneticRegion::Cell {
            tetrahedron,
            cell: CELL,
        })
        .collect()
}

#[test]
fn native_voltage_current_eddy_and_strand_heat_follow_original_temperature_over_two_steps() {
    let run = |alpha| {
        let card = source::absolute_affine_magnetic_card(0.2, alpha, 2.0, 0.0, BASIS);
        let caloric = curve(&card);
        let mut heat = thermal(&caloric, 0.0);
        let p = problem(&card);
        let mut electrical = circuit(&p, 1.0);
        let query = source::absolute_point(300.0);
        let cells = [MaterialCell {
            cell: CELL,
            material_id: MAT,
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &query,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        }];
        let regions = regions();
        let network = SpatialWindingThermalNetwork {
            magnetic: SpatialMagneticThermalNetwork {
                cells: &cells,
                links: &[],
                magnetic_regions: &regions,
                volume_tolerance_m3: 1e-12,
            },
            graph_resistors: &ROUTES,
            fill_tolerance: 1e-12,
        };
        let geometric_r = profile(&p).iter().map(|&j| dot(j, j) / 48.0).sum::<f64>();
        let mut temperatures = Vec::new();
        let mut currents = Vec::new();
        for interval in 1..=2 {
            let old_temperature = heat.temperature_k(CELL).unwrap();
            let next = network
                .step_candidate(
                    &heat,
                    0.2 * f64::from(interval),
                    &electrical,
                    &[10.0],
                    None,
                    thermal_controls(),
                    controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let t = next.thermal().temperature_k(CELL).unwrap();
            let report = next.report();
            let step = &report.coupled.electrical;
            let i = step.currents_a[0];
            let rho = 0.2 * (1.0 + alpha * (t - 300.0));
            let strand = step.cells.iter().map(|c| c.heat_j).sum::<f64>();
            near(strand, 0.2 * i * i * rho * geometric_r, 3e-7);
            near(step.graph_resistors[0].heat_j, 0.2 * i * i * GRAPH_R, 3e-7);
            near(
                step.heat_j,
                strand + step.field.heat_j + step.graph_resistors[0].heat_j,
                3e-7,
            );
            near(4.0 * (t - old_temperature), step.heat_j, 5e-7);
            near(report.coupled.deposited_joule_j, step.heat_j, 5e-7);
            near(report.total_energy_residual_j, 0.0, 5e-7);
            assert!(step.field.numerical_dissipation_j > 0.0);
            assert!(report.coupled.circuit_trials >= step.field_trials);
            for source in next.circuit().conductors() {
                near(source.conductor.resistivity_ohm_m(), rho, 1e-12);
            }
            for source in next.circuit().field().conductors().iter().flatten() {
                near(source.resistivity_ohm_m(), rho, 1e-12);
            }
            temperatures.push(t);
            currents.push(i);
            let (next_heat, next_electrical, _) = next.into_parts();
            heat = next_heat;
            electrical = next_electrical;
        }
        (temperatures, currents)
    };
    let fixed = run(0.0);
    let responsive = run(0.005);
    assert!(responsive.0[1] > responsive.0[0]);
    assert!(responsive.1[0] < fixed.1[0]);
    assert!(responsive.1[1] < fixed.1[1]);
}

#[test]
fn complete_heat_assembly_closes_magnetic_entropy_native_circuit_and_cooling() {
    let slope = 0.0004;
    let card = source::absolute_affine_magnetic_card(0.2, 0.005, 2.0, slope, BASIS);
    let caloric = curve(&card);
    let heat = thermal(&caloric, 3.0);
    let p = problem(&card);
    let circuit = circuit(&p, 1.0);
    let query = source::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CELL,
        material_id: MAT,
        curve: caloric.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &query,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions = regions();
    let laws = p
        .materials()
        .iter()
        .enumerate()
        .map(|(tetrahedron, m)| SpatialMagneticThermalLaw {
            tetrahedron,
            law: MagneticConstitutiveThermalLaw::Linear(
                ThermalMagneticMaterial::new(&card, m.material(), "T", [280.0, 350.0]).unwrap(),
            ),
            maximum_flux_density_t: 10.0,
        })
        .collect::<Vec<_>>();
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &[],
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
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 47,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let prepared = assembly.prepare(&cx, &heat).unwrap();
        let next = heat_assembly_candidate(
            &prepared,
            &cx,
            &heat,
            0.2,
            &circuit,
            &[10.0],
            &regions,
            &ROUTES,
            1e-12,
            1e-12,
            Some(&laws),
            thermal_controls(),
            controls(),
            &[0.0],
        )
        .unwrap();
        let t = next.thermal().temperature_k(CELL).unwrap();
        let r = next.report();
        let step = &r.coupled.electrical;
        let mu = 2.0 + slope * (t - 300.0);
        let field_u = step
            .field
            .endpoint_flux_density_t
            .iter()
            .map(|&b| dot(b, b) / (96.0 * mu) * (1.0 + t * slope / mu))
            .sum::<f64>();
        let end = 0.5 * END_L * step.currents_a[0].powi(2);
        let cooling = -0.2 * 3.0 * (t - 290.0);
        near(r.external_heat_j, cooling, 5e-7);
        near(r.coupled.source.field_internal_energy_j[1], field_u, 5e-7);
        near(
            field_u + end + 4.0 * (t - 300.0) + step.numerical_dissipation_j
                - step.supplied_work_j
                - cooling,
            0.0,
            8e-7,
        );
        near(
            4.0 * (t - 300.0),
            step.heat_j + r.coupled.source.reversible_heat_j + cooling,
            5e-7,
        );
        assert!(!r.coupled.source.fields.is_empty());
        assert!(r.coupled.source.reversible_heat_j.abs() > 1e-10);
        near(r.total_energy_residual_j, 0.0, 5e-7);
        assert_eq!(circuit.time_s(), 0.0);
        near(heat.temperature_k(CELL).unwrap(), 300.0, 1e-12);
    });
}

#[test]
fn finite_receivers_budget_source_and_late_cancellation_preserve_all_histories() {
    let card = source::absolute_affine_magnetic_card(0.2, 0.005, 2.0, 0.0, BASIS);
    let caloric = curve(&card);
    let heat = thermal(&caloric, 0.0);
    let p = problem(&card);
    let original = circuit(&p, 1.0);
    let query = source::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CELL,
        material_id: MAT,
        curve: caloric.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &query,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions = regions();
    let network = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &ROUTES,
        fill_tolerance: 1e-12,
    };
    let partial = circuit(&p, 0.5);
    assert!(
        network
            .step_candidate(
                &heat,
                0.2,
                &partial,
                &[10.0],
                None,
                thermal_controls(),
                controls(),
                &[0.0],
                || false
            )
            .is_err()
    );
    let foreign = source::absolute_affine_magnetic_card(0.3, 0.005, 2.0, 0.0, BASIS);
    let changed_cells = [MaterialCell {
        cell: CELL,
        material_id: MAT,
        curve: caloric.curve(),
        point: ThermalMaterialPoint {
            material: &foreign,
            context: &query,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let changed = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &changed_cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &ROUTES,
        fill_tolerance: 1e-12,
    };
    assert!(
        changed
            .step_candidate(
                &heat,
                0.2,
                &original,
                &[10.0],
                None,
                thermal_controls(),
                controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let mut tiny = thermal_controls();
    tiny.maximum_circuit_trials = 1;
    assert!(matches!(
        network.step_candidate(
            &heat,
            0.2,
            &original,
            &[10.0],
            None,
            tiny,
            controls(),
            &[0.0],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let mut polls = 0;
    let good = network
        .step_candidate(
            &heat,
            0.2,
            &original,
            &[10.0],
            None,
            thermal_controls(),
            controls(),
            &[0.0],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    let aborted = network.step_candidate(
        &heat,
        0.2,
        &original,
        &[10.0],
        None,
        thermal_controls(),
        controls(),
        &[0.0],
        || {
            seen += 1;
            seen >= polls - 1
        },
    );
    assert!(matches!(aborted, Err(CircuitEnthalpyError::Cancelled)));
    let retry = network
        .step_candidate(
            &heat,
            0.2,
            &original,
            &[10.0],
            None,
            thermal_controls(),
            controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(good.circuit().state(), retry.circuit().state());
    assert_eq!(
        good.circuit().field().potential_wb(),
        retry.circuit().field().potential_wb()
    );
    assert_eq!(
        good.thermal().temperature_k(CELL).unwrap().to_bits(),
        retry.thermal().temperature_k(CELL).unwrap().to_bits()
    );
    assert_eq!(original.time_s(), 0.0);
    near(heat.temperature_k(CELL).unwrap(), 300.0, 1e-12);
}

#[test]
fn separate_material_graph_body_changes_current_and_retains_its_source_in_complete_heat() {
    const RES_CELL: CellId = CellId(8);
    const RES_MAT: MaterialId = MaterialId(10);
    let run = |alpha| {
        let coil_card = source::absolute_affine_magnetic_card(0.2, 0.0, 2.0, 0.0, BASIS);
        let resistor_card = source::absolute_affine_magnetic_card(GRAPH_R, alpha, 2.0, 0.0, BASIS);
        let coil_curve = curve(&coil_card);
        let resistor_curve = curve(&resistor_card);
        let heat = DepositionSimulation::new(
            [(CELL, MAT), (RES_CELL, RES_MAT)]
                .into_iter()
                .map(|(id, material_id)| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(DepositedCell {
                        id,
                        segment_id: 0,
                        body_id: id.0,
                        material_id,
                        mass_kg: 1.0,
                        specific_h_j_kg: 400.0,
                        geometry: CellGeometry {
                            start_m: [0.0; 3],
                            end_m: [1.0, 0.0, 0.0],
                            volume_m3: 1.0,
                            surface_area_m2: 1.0,
                        },
                        boundary: Boundary {
                            ambient_k: 300.0,
                            htc_w_m2_k: 0.0,
                            reservoir: None,
                            source_w: 0.0,
                        },
                    }),
                })
                .collect(),
            &[
                (MAT, coil_curve.curve() as &dyn ThermalLaw),
                (RES_MAT, resistor_curve.curve() as &dyn ThermalLaw),
            ],
            DepositionBudget {
                max_cells: 2,
                max_links: 2,
                max_events: 2,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let p = problem(&coil_card);
        let original = circuit(&p, 1.0);
        let query = source::absolute_point(300.0);
        let cells = [
            MaterialCell {
                cell: CELL,
                material_id: MAT,
                curve: coil_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &coil_card,
                    context: &query,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
            MaterialCell {
                cell: RES_CELL,
                material_id: RES_MAT,
                curve: resistor_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &resistor_card,
                    context: &query,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
        ];
        let routes = [ResistorCell {
            branch: 1,
            cell: RES_CELL,
        }];
        let resistors = [CellResistor {
            destination: routes[0],
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 1.0,
            area_m2: 1.0,
        }];
        let bases = [(1, BASIS)];
        let regions = regions();
        let layout = SpatialWindingThermalNetwork {
            magnetic: SpatialMagneticThermalNetwork {
                cells: &cells,
                links: &[],
                magnetic_regions: &regions,
                volume_tolerance_m3: 1e-12,
            },
            graph_resistors: &routes,
            fill_tolerance: 1e-12,
        };
        // Explicitly reconfiguring the electrical source does not establish
        // agreement with the accepted caloric state; the coupled owner checks it.
        let wrong = ResolvedConductorState::try_from_material_with_basis(
            source::bundle_at(&resistor_card, &source::absolute_point(301.0), BASIS),
            BASIS,
        )
        .unwrap()
        .into_uniform_conductor(1.0, 1.0)
        .unwrap();
        let mismatched = original
            .with_graph_resistors(
                vec![GraphResistorMaterial {
                    branch: 1,
                    conductor: wrong,
                }],
                &mut || false,
            )
            .unwrap();
        assert!(
            layout
                .step_material_candidate(
                    &heat,
                    0.2,
                    &mismatched,
                    &[10.0],
                    &resistors,
                    &bases,
                    None,
                    thermal_controls(),
                    controls(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        let first = layout
            .step_material_candidate(
                &heat,
                0.2,
                &original,
                &[10.0],
                &resistors,
                &bases,
                None,
                thermal_controls(),
                controls(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        let verify=|before:&DepositionSimulation<'_>,after:&DepositionSimulation<'_>,
            step:&fs_em::spatial::vector_potential::induction::nonlinear::winding::SpatialNonlinearWindingStep| {
            let tr=after.temperature_k(RES_CELL).unwrap();
            let r=GRAPH_R*(1.0+alpha*(tr-300.0));
            let i=step.currents_a[0];
            assert_eq!(step.graph_materials.len(),1);
            near(step.graph_materials[0].conductor.resistance_ohm(),r,1e-12);
            near(step.graph_resistors[0].heat_j,0.2*r*i*i,3e-7);
            near(4.0*(tr-before.temperature_k(RES_CELL).unwrap()),step.graph_resistors[0].heat_j,5e-7);
            near(4.0*(after.temperature_k(CELL).unwrap()-before.temperature_k(CELL).unwrap()),
                step.field.heat_j+step.cells.iter().map(|cell|cell.heat_j).sum::<f64>(),5e-7);
        };
        verify(&heat, first.thermal(), &first.report().coupled.electrical);
        assert!(
            layout
                .step_candidate(
                    first.thermal(),
                    0.4,
                    first.circuit(),
                    &[10.0],
                    None,
                    thermal_controls(),
                    controls(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        let material = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &material,
            contacts: &[],
            electrical_bases: &bases,
        };
        let assembly = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 91,
                    kernel_id: 3,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let prepared = assembly.prepare(&cx, first.thermal()).unwrap();
            let second = heat_assembly_candidate(
                &prepared,
                &cx,
                first.thermal(),
                0.4,
                first.circuit(),
                &[10.0],
                &regions,
                &routes,
                1e-12,
                1e-12,
                None,
                thermal_controls(),
                controls(),
                &[0.0; 2],
            )
            .unwrap();
            verify(
                first.thermal(),
                second.thermal(),
                &second.report().coupled.electrical,
            );
            near(second.report().total_energy_residual_j, 0.0, 5e-7);
            assert_eq!(original.time_s(), 0.0);
            near(heat.temperature_k(RES_CELL).unwrap(), 300.0, 1e-12);
            (
                first.circuit().currents_a()[0],
                second.circuit().currents_a()[0],
            )
        })
    };
    let fixed = run(0.0);
    let responsive = run(0.005);
    assert!(responsive.0 < fixed.0);
    assert!(responsive.1 < fixed.1);
}
