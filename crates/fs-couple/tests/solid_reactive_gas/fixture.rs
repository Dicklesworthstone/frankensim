//! Manufactured isomers and two spatial tetrahedra, using the original owners.
use super::*;
use fs_couple::cavity_phs::piston::reactive::material_rate::{
    ACTIVITY_RATE_SCALE_DIMS, ACTIVITY_RATE_SCALE_PROPERTY, MaterialRateScale,
};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::*;

// The existing fixture supplies actual NASA-9 models and source transport cards.
#[allow(dead_code)]
#[path = "../material_exchange/gas.rs"]
mod gas_data;

pub const NODES: [[f64; 3]; 5] = [
    [0., 0., 0.],
    [1., 0., 0.],
    [0., 1., 0.],
    [0., 0., 1.],
    [0., 0., -1.],
];
pub const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 2, 1, 4]];
pub const FREE: [usize; 6] = [0, 1, 2, 4, 5, 8];
pub const VOLUME: f64 = 0.1;

pub fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 3,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
pub fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
pub fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with("lot", 7.0)
        .unwrap()
}

pub fn reaction(rate_at_300: f64, formation_release: f64, product_k_factor: f64) -> GasReaction {
    let (a, ta) = gas_data::component("A", 0.028, 0.0);
    let (b, mut tb) = gas_data::component("B", 0.028, -formation_release);
    let mut claims = ClaimSet::new();
    for (_, original) in tb.material.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Scalar {
                value: 0.03 * product_k_factor,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    tb.material = MaterialCard::assemble(tb.material.id().clone(), claims, vec![]).unwrap();
    let species = vec![a.species().clone(), b.species().clone()];
    let id = ReactionId::new("manufactured-isomerization").unwrap();
    let matrix =
        StoichiometricMatrix::new(species.clone(), vec![id.clone()], vec![vec![-1], vec![1]])
            .unwrap();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        species.clone(),
        vec![vec![1, 1]],
    )
    .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let certificate = verify_conservation(&elements, &matrix, &charge).unwrap();
    let equilibrium =
        IdealGasReactionEquilibriumV1::new(matrix, certificate, id, vec![a.clone(), b.clone()])
            .unwrap();
    let kinetics = DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(
            equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(1.0)).unwrap(),
        ),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    );
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![ta, tb],
        None,
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (200.0, rate_at_300 * 200.0 / 300.0),
                    (2000.0, rate_at_300 * 2000.0 / 300.0),
                ],
                dims: ACTIVITY_RATE_SCALE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    200.0,
                    2000.0,
                )
                .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7)
                .with("catalyst-lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured spatial reaction rate".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic catalyst".into(),
            phase: "gas".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    GasReaction::new(catalog, kinetics)
        .unwrap()
        .with_material_rate(
            MaterialRateScale::new(
                card,
                MaterialPropertySelection::SingleClaimOnly,
                Pressure::new(100000.0),
            )
            .unwrap(),
        )
        .unwrap()
}
pub fn gas(reaction: &GasReaction) -> IsochoricGasParcel {
    let model = reaction
        .catalog()
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap()
        .with("catalyst-lot", 7.0)
        .unwrap();
    IsochoricGasParcel::new(model, point, VOLUME).unwrap()
}
pub fn wall(k_at_300: f64, varying_k: bool) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, value) in [
        (heat.clone(), 10.0),
        (density.clone(), 6.0),
        (PropertyKey::new("young_modulus", Pressure::DIMS), 1000.0),
        (PropertyKey::new("poisson_ratio", Dims([0; 6])), 0.0),
        (
            PropertyKey::new(
                "linear_thermal_expansion_coefficient",
                INVERSE_TEMPERATURE_DIMS,
            ),
            1e-5,
        ),
        (
            PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            k_at_300,
        ),
    ] {
        let dims = key.quantity().dims();
        let varying = varying_k && key.name() == "thermal_conductivity";
        let data = if varying {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (280.0, value * 280.0 / 300.0),
                    (600.0, value * 600.0 / 300.0),
                ],
                dims,
            }
        } else {
            PropertyValue::Scalar { value, dims }
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value: data,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        280.0,
                        600.0,
                    )
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if varying {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured spatial solid".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored specimen".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(600.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-10,
    )
    .unwrap();
    (card, caloric)
}
pub fn initial(curve: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events = [30, 10]
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 77,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: 200.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64],
                    end_m: [1.0, 0.0, i as f64],
                    volume_m3: 1.0 / 6.0,
                    surface_area_m2: 3.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(30),
            b: CellId(10),
            area_m2: 0.5,
            conductance_w_k: 0.5,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve as &dyn ThermalLaw)],
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
            max_dt_s: 0.125,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
pub fn solid<'a>(
    thermal: &DepositionSimulation<'_>,
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    cx: &Cx<'_>,
) -> SolidEnthalpy<'a> {
    let cells = [30, 10].map(|id| TetCellMaterial {
        cell: CellId(id),
        curve,
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
        stiffness_model: ThermalStiffnessModel::Isotropic,
        expansion_model: ThermalExpansionModel::Isotropic,
        expansion_selection: MaterialPropertySelection::SingleClaimOnly,
    });
    SolidEnthalpy::try_new(
        thermal,
        &NODES,
        &TETS,
        cells.to_vec(),
        TetAssemblyBudget::standard(),
        0.1,
        config().mesh_frame,
        cx,
    )
    .unwrap()
}
pub fn contacts() -> [SolidGasContact; 2] {
    [(99, 30, 0.1), (8, 10, 0.2)].map(|(id, cell, area_m2)| SolidGasContact {
        id,
        cell: CellId(cell),
        geometry: GasConductionGeometry {
            area_m2,
            gas_distance_m: 0.001,
            cell_distance_m: 0.02,
        },
        conductivity: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        selection: MaterialPropertySelection::SingleClaimOnly,
    })
}
pub fn config() -> SolidEquilibriumConfig {
    SolidEquilibriumConfig {
        static_solve: TetStaticSolveConfig::default(),
        mesh_frame: ContentHash([1; 32]),
        maximum_nodal_displacement_m: 0.1,
    }
}
pub fn solve() -> SolidGasSolve {
    SolidGasSolve {
        maximum_cells: 2,
        maximum_contacts: 4,
        maximum_trials: 500,
        maximum_newton_iterations: 20,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 1e-7,
        thermal_energy_tolerance_j: 1e-7,
        gas_energy_tolerance_j: 1e-7,
        total_energy_tolerance_j: 2e-7,
        maximum_gas_temperature_increment_k: 200.0,
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-9,
            maximum_iterations: 128,
        },
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..Default::default()
        },
    }
}
pub fn chemical() -> GasReactionSolve {
    GasReactionSolve {
        gas: solve().gas,
        maximum_consumed_fraction: 0.2,
        species_tolerance_mol: 1e-12,
        entropy_tolerance_j_k: 1e-10,
    }
}

// Independent three-reservoir balance for constant capacities/conductivities.
// This is an analytical test oracle; production only uses its existing owners.
pub fn reference(chemical_t: f64, gas_capacity: f64, gas_k: f64) -> [f64; 3] {
    let [g1, g2] = [0.1, 0.2].map(|area| area / (0.001 / gas_k + 0.02));
    let (c, h) = (10.0, 0.5);
    let (a, b) = (c + h + g1, c + h + g2);
    let determinant = a * b - h * h;
    let a1 = (b + h) * c * 300.0 / determinant;
    let a2 = (h + a) * c * 300.0 / determinant;
    let b1 = (b * g1 + h * g2) / determinant;
    let b2 = (h * g1 + a * g2) / determinant;
    let gas_t = (gas_capacity * chemical_t + g1 * a1 + g2 * a2)
        / (gas_capacity + g1 + g2 - g1 * b1 - g2 * b2);
    [a1 + b1 * gas_t, a2 + b2 * gas_t, gas_t]
}
