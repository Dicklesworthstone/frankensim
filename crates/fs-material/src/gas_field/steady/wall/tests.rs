//! Manufactured source laws exercise the real gas/thermal owners, not fitted audio.
use super::*;
use crate::fluid_transport::FluidTransport;
use crate::gas_field::PrimitiveGasState;
use crate::gas_field::steady::network::{GasNetworkCell, GasOutlet};
use crate::gas_field::steady::{AxialGasCell, GasInjection};
use crate::gas_species::input::ThermochemicalPhase;
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use crate::state_point::THERMAL_CONDUCTIVITY_DIMS;
use crate::thermal::source::ThermalMaterialPoint;
use crate::thermochemical_gas::{energy::GasEnergySolve, tests::model};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::Temperature;
use fs_thermal::enthalpy::CellGeometry;
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const FLOW: f64 = 5e-5;
const CP: f64 = 3.5 * R / 0.028;

fn controls() -> SteadyGasSolve {
    SteadyGasSolve {
        caloric: GasEnergySolve {
            minimum_temperature_k: 275.0,
            maximum_temperature_k: 350.0,
            energy_tolerance_j_kg: 1e-8,
            maximum_iterations: 64,
        },
        maximum_dt_s: 0.05,
        maximum_relaxations: 2000,
        maximum_substeps: 100000,
        steady_power_tolerance_w: 1e-8,
        energy_tolerance_j: 1e-8,
        maximum_mean_mach: 0.05,
    }
}

fn wall_controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-12,
        relative_heat_tolerance: 1e-13,
    }
}

fn source(knots: Vec<(f64, f64)>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured wall conductivity; numerical verification only".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn caloric(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [(0.0, 200.0), (20_000.0, 400.0)]
            .into_iter()
            .map(|(h, t)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            })
            .collect(),
    )
    .unwrap()
}

fn layer<'a>(
    card: &'a MaterialCard,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
    range: [f64; 2],
) -> ResolvedConductionLayer<'a> {
    ResolvedConductionLayer::new(
        caloric,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        SolidLiquidPhase::Solid,
        range,
        || Ok(()),
    )
    .unwrap()
}

fn fixture() -> SteadyGasNetwork {
    SteadyGasNetwork {
        source: model(vec![0.75, 0.25], 0.0),
        source_conditions: QueryPoint::new(),
        static_pressure_pa: 101325.0,
        phase: ThermochemicalPhase::FrozenGas,
        cells: vec![GasNetworkCell {
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [0.1, 0.0, 0.0],
                volume_m3: 2e-5,
                surface_area_m2: 0.005,
            },
            flow_area_m2: 2e-4,
            wall_conductance_w_k: 0.777,
            wall_temperature_k: 290.0,
            initial_temperature_k: 310.0,
        }],
        flows: vec![],
        inlets: vec![GasInjection {
            cell: 0,
            mass_flow_kg_s: FLOW,
            state: PrimitiveGasState {
                temperature_k: 330.0,
                composition: Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0])
                    .unwrap(),
            },
        }],
        outlets: vec![GasOutlet {
            cell: 0,
            mass_flow_kg_s: FLOW,
        }],
    }
}

#[test]
fn changing_only_wall_source_propagates_through_heat_transport_and_acoustics() {
    let problem = fixture();
    let before = problem.clone();
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let geometry = ReservoirLayerGeometry::planar(0.005, 0.1).unwrap();
    let mut states = Vec::new();
    for k in [0.2, 2.0] {
        let card = source(vec![(200.0, k), (400.0, k)]);
        let caloric = caloric(&card);
        let layer = layer(&card, &caloric, &context, [200.0, 400.0]);
        let wall = MaterialGasWall::new(0, &layer, geometry, 10.0).unwrap();
        assert!(core::ptr::eq(wall.layer(), &layer));
        let solved = problem
            .solve_with_material_walls(&[wall], wall_controls(), controls(), || false)
            .unwrap();
        let g = 1.0 / (10.0 + 0.1 / (0.005 * k));
        let expected_t = (FLOW * CP * 330.0 + g * 290.0) / (FLOW * CP + g);
        let state = &solved.states()[0];
        assert!((state.state().temperature - expected_t).abs() < 2e-5);
        state.state().validate_coherence().unwrap();
        assert_eq!(
            state.transport_state().unwrap().density_kg_m3(),
            state.state().density
        );
        assert!(solved.maximum_power_residual_bound_w() <= controls().steady_power_tolerance_w);
        assert!(solved.maximum_energy_residual_j() <= controls().energy_tolerance_j);
        states.push(state.clone());
    }
    assert!(states[0].state().temperature > states[1].state().temperature + 5.0);
    assert!(states[0].state().sound_speed > states[1].state().sound_speed);
    assert!(states[0].state().density < states[1].state().density);
    assert_ne!(
        states[0].state().dynamic_viscosity,
        states[1].state().dynamic_viscosity
    );
    assert_ne!(
        states[0].state().thermal_conductivity,
        states[1].state().thermal_conductivity
    );
    assert_eq!(problem, before);
}

#[test]
fn nonlinear_source_and_film_close_exact_heat_balance_and_publish_existing_axial_field() {
    let card = source(vec![(200.0, 0.1), (300.0, 0.3), (400.0, 0.5)]);
    let caloric = caloric(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&card, &caloric, &context, [200.0, 400.0]);
    let wall = MaterialGasWall::new(
        0,
        &layer,
        ReservoirLayerGeometry::planar(0.005, 0.1).unwrap(),
        10.0,
    )
    .unwrap();
    // k(T)=-0.3+0.002*T. A 310 K material face and 290 K reservoir transfer
    // exactly 0.3 W; the 10 K/W film gives 313 K bulk gas. Choose its inlet
    // enthalpy analytically to supply this exact heat.
    let mut problem = fixture();
    problem.inlets[0].state.temperature_k = 313.0 + 0.3 / (FLOW * CP);
    let solved = problem
        .solve_with_material_walls(&[wall], wall_controls(), controls(), || false)
        .unwrap();
    let t = solved.states()[0].state().temperature;
    assert!((t - 313.0).abs() < 2e-5);
    let response = wall
        .evaluate([t, 290.0], wall_controls(), || Ok(()))
        .unwrap();
    assert!((response.material_face_temperature_k - 310.0).abs() < 2e-5);
    assert!(
        (response.power_to_reservoir_w - FLOW * CP * (problem.inlets[0].state.temperature_k - t))
            .abs()
            < 2e-8
    );

    // A loose constitutive solve cannot hide its flux defect in a tight gas
    // steady-residual claim, even if the approximate gas iteration converges.
    assert!(matches!(
        problem.solve_with_material_walls(
            &[wall],
            LayerConductionConfig {
                absolute_heat_tolerance_w: 1.0,
                ..wall_controls()
            },
            controls(),
            || false
        ),
        Err(SteadyGasError::PublishedBalance { .. })
    ));

    let axial = SteadyAxialGas {
        source: problem.source.clone(),
        source_conditions: problem.source_conditions.clone(),
        static_pressure_pa: problem.static_pressure_pa,
        phase: problem.phase,
        inlet: problem.inlets[0].state.clone(),
        inlet_mass_flow_kg_s: FLOW,
        cells: vec![AxialGasCell {
            end_m: 0.1,
            volume_m3: 2e-5,
            wall_area_m2: 0.005,
            wall_conductance_w_k: 0.777,
            wall_temperature_k: 290.0,
            initial_temperature_k: 310.0,
        }],
        injections: vec![],
    };
    let result = axial
        .solve_with_material_walls(&[wall], wall_controls(), controls(), || false)
        .unwrap();
    let emitted = result.field().state_at(0.1).unwrap();
    assert_eq!(emitted, solved.states()[0]);
    emitted.state().validate_coherence().unwrap();
    assert!(result.maximum_power_residual_bound_w() <= controls().steady_power_tolerance_w);
}

#[test]
fn source_refusal_and_cancellation_publish_nothing_and_empty_bindings_preserve_fixed_g() {
    let problem = fixture();
    let before = problem.clone();
    let card = source(vec![(200.0, 0.1), (400.0, 0.1)]);
    let caloric = caloric(&card);
    let context = QueryPoint::new().with("T", 300.0).unwrap();
    // Initial data fit this source, but its true warm steady state does not.
    let layer = layer(&card, &caloric, &context, [290.0, 315.0]);
    let geometry = ReservoirLayerGeometry::planar(0.005, 0.1).unwrap();
    let wall = MaterialGasWall::new(0, &layer, geometry, 0.0).unwrap();
    assert!(matches!(
        problem.solve_with_material_walls(&[wall], wall_controls(), controls(), || false),
        Err(SteadyGasError::Thermal(DepositionError::MaterialDomain))
    ));
    assert!(matches!(
        problem.solve_with_material_walls(&[wall], wall_controls(), controls(), || true),
        Err(SteadyGasError::Thermal(DepositionError::Cancelled))
    ));
    assert!(
        problem
            .solve_with_material_walls(&[wall, wall], wall_controls(), controls(), || false)
            .is_err()
    );
    let missing = MaterialGasWall::new(1, &layer, geometry, 0.0).unwrap();
    assert!(
        problem
            .solve_with_material_walls(&[missing], wall_controls(), controls(), || false)
            .is_err()
    );
    let excess = MaterialGasWall::new(
        0,
        &layer,
        ReservoirLayerGeometry::planar(0.006, 0.1).unwrap(),
        0.0,
    )
    .unwrap();
    assert!(
        problem
            .solve_with_material_walls(&[excess], wall_controls(), controls(), || false)
            .is_err()
    );
    assert!(MaterialGasWall::new(0, &layer, geometry, -1.0).is_err());
    assert_eq!(problem, before);
    assert_eq!(
        problem
            .solve_with_material_walls(&[], wall_controls(), controls(), || false)
            .unwrap(),
        problem.solve(controls(), || false).unwrap()
    );
}
