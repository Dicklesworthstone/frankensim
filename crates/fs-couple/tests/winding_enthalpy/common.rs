#![allow(dead_code)]
use super::field;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell};
use fs_couple::circuit_enthalpy::winding::{MagneticRegion, WindingThermalNetwork};
use fs_couple::circuit_enthalpy::{CircuitEnthalpySolve, ResistorCell};
use fs_couple::material_circuit::MaterialCircuitSolve;
use fs_em::winding::ohmic::OhmicWindingCircuit;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

pub fn cards(scale: f64, alpha: f64) -> [MaterialCard; 2] {
    [2e-5 * scale, 3e-5].map(|rho| {
        let source = field::card(rho, alpha, 2e-4, ElectricalPropertyBasis::Resistivity);
        let mut claims = source.claims().clone();
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
                value: PropertyValue::Scalar {
                    value: 1.0,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    2000.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic winding thermal transfer fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic winding".into(),
                phase: "solid".into(),
                process: "test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap()
    })
}
pub fn curves(cards: &[MaterialCard; 2]) -> [ResolvedHeatCapacityCurve; 2] {
    cards.each_ref().map(|card| {
        ResolvedHeatCapacityCurve::try_new(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(
                    "specific_heat_capacity",
                    SPECIFIC_HEAT_CAPACITY_DIMS,
                ),
                density: &PropertyKey::new("density", Density::DIMS),
                temperature_axis: "T",
                lower: &field::point(200.0),
                upper: &field::point(2000.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            0.0,
            1e-6,
        )
        .unwrap()
    })
}
pub fn config() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 8,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 1e-6,
        thermal_energy_tolerance_j: 1e-6,
        total_energy_tolerance_j: 2e-6,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-8,
            power_tolerance_w: 1e-8,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
pub fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
pub fn initial<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    first_volume: f64,
    event: Option<f64>,
) -> DepositionSimulation<'a> {
    initial_with_core(curves, first_volume, event, None)
}
pub fn initial_with_core<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    first_volume: f64,
    event: Option<f64>,
    core: Option<&'a ResolvedHeatCapacityCurve>,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (index, (id, material_id, volume, source_w)) in [
        (CellId(30), MaterialId(0), first_volume, 0.0),
        (CellId(10), MaterialId(1), 0.25, 0.0),
        (CellId(50), MaterialId(2), 0.1, 0.7),
    ]
    .into_iter()
    .enumerate()
    {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: index as u64,
                body_id: index as u64,
                material_id,
                mass_kg: volume,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [index as f64, 0.0, 0.0],
                    end_m: [index as f64 + 1.0, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 1.0,
                },
                boundary: boundary(source_w),
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(30),
            b: CellId(10),
            area_m2: 0.5,
            conductance_w_k: 999.0,
        }),
    });
    if let Some(time_s) = event {
        events.push(TimedEvent {
            time_s,
            event: Event::Boundary(CellId(30), boundary(7.0)),
        });
    }
    let mut laws = vec![
        (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
        (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        (MaterialId(2), curves[1].curve() as &dyn ThermalLaw),
    ];
    if let Some(curve) = core {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(40),
                segment_id: 3,
                body_id: 3,
                material_id: MaterialId(3),
                mass_kg: 0.015625,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [3.0, 0.0, 0.0],
                    end_m: [4.0, 0.0, 0.0],
                    volume_m3: 0.015625,
                    surface_area_m2: 1.0,
                },
                boundary: boundary(1.0),
            }),
        });
        laws.push((MaterialId(3), curve.curve() as &dyn ThermalLaw));
    }
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
pub fn cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> Vec<MaterialCell<'a>> {
    [
        (
            CellId(30),
            MaterialId(0),
            0,
            Some((key, MaterialPropertySelection::SingleClaimOnly)),
        ),
        (
            CellId(10),
            MaterialId(1),
            1,
            Some((key, MaterialPropertySelection::SingleClaimOnly)),
        ),
        (CellId(50), MaterialId(2), 1, None),
    ]
    .into_iter()
    .map(|(cell, material_id, index, conductivity)| MaterialCell {
        cell,
        material_id,
        curve: curves[index].curve(),
        point: ThermalMaterialPoint {
            material: &cards[index],
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity,
    })
    .collect()
}
pub fn routes(circuit: &OhmicWindingCircuit<'_, '_>) -> Vec<MagneticRegion> {
    circuit
        .conductors()
        .iter()
        .zip(circuit.conductor_windings())
        .map(|(s, &w)| MagneticRegion::Cell {
            triangle: s.triangle,
            cell: if w == 0 { CellId(30) } else { CellId(10) },
        })
        .collect()
}
pub const GRAPH: [ResistorCell; 3] = [
    ResistorCell {
        branch: 1,
        cell: CellId(30),
    },
    ResistorCell {
        branch: 3,
        cell: CellId(50),
    },
    ResistorCell {
        branch: 5,
        cell: CellId(10),
    },
];
pub const LINKS: [ConductionLink; 1] = [ConductionLink {
    link: 7,
    distance_a_m: 0.125,
    distance_b_m: 0.125,
}];
pub fn network<'a>(
    cells: &'a [MaterialCell<'a>],
    routes: &'a [MagneticRegion],
) -> WindingThermalNetwork<'a> {
    WindingThermalNetwork {
        cells,
        links: &LINKS,
        magnetic_regions: routes,
        graph_resistors: &GRAPH,
        volume_tolerance_m3: 1e-12,
    }
}
pub fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-6, "{a} != {b}");
}
