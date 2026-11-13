//! Shared original-source fixtures for layered thermal shell integration tests.
#![allow(dead_code)]

use fs_couple::laminate_thermal::{LaminateShellDiscretization, LaminateShellLimits};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY as CP_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY as K_PROPERTY,
};
use fs_material::thermal::conduction::LayerConductionConfig;
use fs_material::thermal::transport::stack::LayerStackSourceQuery;
use fs_material::thermal::transport::stack::transient::{
    BoundLayeredThermalBody, LayeredThermalIds, PreparedLayeredThermalBody,
    TransientLayerSourceQuery,
    implicit::{LayeredThermalSolveControls, LayeredThermalWorkspace},
};
use fs_material::thermal::transport::storage::{HeatStorageAmount, HeatStorageComponentQuery};
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, NodalDirectors, ShellConstraint, ShellMass,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{DepositionBudget, DepositionSimulation, Event, StepConfig};
use std::sync::Arc;

pub const AREA: f64 = 0.6;
pub const H: f64 = 0.02;
pub const RHO: f64 = 1000.0;
pub const CP: f64 = 100.0;
pub const K: f64 = 10.0;
pub const E: f64 = 2e8;

pub fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
pub fn at(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
pub fn young(t: f64) -> f64 {
    E * (1.0 - 0.002 * (t - 300.0))
}
pub fn free(t: f64) -> f64 {
    let dt = t - 300.0;
    1e-5 * dt + 0.5e-7 * dt * dt
}
pub fn card(name: &str) -> Arc<MaterialCard> {
    card_with_coefficients(name, -0.002 * E, 1e-5, 1e-7, CP)
}
pub fn card_with_coefficients(
    name: &str,
    young_slope_pa_k: f64,
    alpha_reference_per_k: f64,
    alpha_slope_per_k2: f64,
    specific_heat_j_kg_k: f64,
) -> Arc<MaterialCard> {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (
            DENSITY_PROPERTY,
            Density::DIMS,
            vec![(280.0, 1040.0), (300.0, RHO), (360.0, 880.0)],
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            vec![
                (280.0, E - 20.0 * young_slope_pa_k),
                (300.0, E),
                (360.0, E + 60.0 * young_slope_pa_k),
            ],
        ),
        ("poisson_ratio", Dims::NONE, vec![(300.0, 0.0)]),
        (
            CTE,
            INVERSE_TEMPERATURE_DIMS,
            vec![
                (280.0, alpha_reference_per_k - 20.0 * alpha_slope_per_k2),
                (320.0, alpha_reference_per_k + 20.0 * alpha_slope_per_k2),
                (360.0, alpha_reference_per_k + 60.0 * alpha_slope_per_k2),
            ],
        ),
        (K_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, vec![(300.0, K)]),
        (
            CP_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            vec![(300.0, specific_heat_j_kg_k)],
        ),
    ] {
        let constant = values.len() == 1;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if constant {
                    PropertyValue::Scalar {
                        value: values[0].1,
                        dims,
                    }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: values,
                        dims,
                    }
                },
                interpolation: if constant {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature(),
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored coupled laminate oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "solid".into(),
                process: "numerical thermal shell".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}
pub fn mechanical(cards: &[Arc<MaterialCard>]) -> Vec<LaminaSourceInput> {
    cards
        .iter()
        .map(|card| LaminaSourceInput {
            card: Arc::clone(card),
            reference: at(300.0),
            thickness_m: H / cards.len() as f64,
            orientation_rad: 0.0,
            elasticity: LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungPoisson {
                young_modulus: ScalarPropertyRequirement::try_new(
                    "young_modulus",
                    Pressure::DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap(),
                poisson_ratio: ScalarPropertyRequirement::try_new(
                    "poisson_ratio",
                    Dims::NONE,
                    ScalarAdmissibility::Finite,
                )
                .unwrap(),
            }),
            expansion: LaminaExpansionSource::Isotropic,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect()
}
pub fn prepare_body(
    cards: &[Arc<MaterialCard>],
    counts: &[usize],
    mass_scale: f64,
) -> PreparedLayeredThermalBody {
    let lower = at(280.0);
    let upper = at(360.0);
    let conductivity = PropertyKey::new(K_PROPERTY, THERMAL_CONDUCTIVITY_DIMS);
    let capacity = PropertyKey::new(CP_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let h = H / cards.len() as f64;
    let mass = mass_scale * AREA * h * RHO;
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let queries: Vec<_> = cards
        .iter()
        .zip(counts)
        .enumerate()
        .map(|(i, (card, &count))| TransientLayerSourceQuery {
            conduction: LayerStackSourceQuery {
                id: i as u64,
                card,
                conductivity: &conductivity,
                thickness_m: h,
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            capacity: HeatStorageComponentQuery {
                id: i as u64,
                card,
                capacity: &capacity,
                amount: HeatStorageAmount::MassKg(mass),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection,
            },
            mass_kg: mass,
            cell_count: count,
        })
        .collect();
    PreparedLayeredThermalBody::try_new(&queries, &vec![0.0; cards.len() - 1], AREA, 300.0).unwrap()
}
pub fn ids() -> LayeredThermalIds {
    LayeredThermalIds {
        first_cell: 10,
        first_material: 30,
        first_link: 50,
        body_id: 70,
    }
}
pub fn conduction() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12,
    }
}
pub fn initial<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    temperatures: &[f64],
) -> DepositionSimulation<'a> {
    initial_with_sources(bound, temperatures, &vec![0.0; temperatures.len()])
}
pub fn initial_with_sources<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    temperatures: &[f64],
    sources_w: &[f64],
) -> DepositionSimulation<'a> {
    assert_eq!(temperatures.len(), sources_w.len());
    let mut events = bound
        .birth_events(0.0, temperatures, conduction(), || Ok(()))
        .unwrap();
    let mut cell = 0;
    for event in &mut events {
        if let Event::Birth(born) = &mut event.event {
            born.boundary.source_w = sources_w[cell];
            cell += 1;
        }
    }
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: temperatures.len(),
            max_links: temperatures.len() - 1,
            max_events: 2 * temperatures.len() - 1,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
pub fn discretization() -> LaminateShellDiscretization {
    let mesh = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.6, 0.0],
            [1.0, 0.6, 0.0],
        ],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 4], 0.9).unwrap();
    LaminateShellDiscretization {
        mesh,
        directors,
        constraints: vec![ShellConstraint {
            node: 0,
            mask: [true; 6],
            frame: ConstraintFrame::Global,
        }],
        section_angles_rad: vec![0.0; 2],
        mass: ShellMass::Consistent,
        budget: MitcBudget {
            max_nodes: 4,
            max_triangles: 2,
        },
        drilling_alpha: 1e-3,
    }
}
pub fn limits() -> LaminateShellLimits {
    LaminateShellLimits {
        maximum_translation_m: 0.2,
        maximum_rotation_rad: 0.5,
        maximum_reduced_elastic_strain: 0.01,
    }
}
pub fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}
