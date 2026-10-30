//! Native pressure/heat/material continuation with an independent one-DOF oracle.
use super::*;
use fs_couple::solid_enthalpy::gas::{GasConductionGeometry, SolidGasContact};
use fs_couple::solid_enthalpy::thermoelastic::dynamic::cavity::{CavitySolve, ThermoelasticCavity};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_solid::loads3::SurfaceLoad3;
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::{Composition, CompositionBasis};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

#[allow(dead_code)]
#[path = "../../../material_exchange/gas.rs"]
mod gas_fixture;

const GAS_VOLUME: f64 = 0.05;
const PRESSURE: f64 = 2000.0;
const PATCH_AREA: f64 = 0.01;
const CP: f64 = 0.2;

fn gas(t: f64, formation: f64) -> ClosedGasParcel {
    let (species, transport) = gas_fixture::component("cavity-fixture", 0.028, formation);
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport], None,
    ).unwrap();
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), t).unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), PRESSURE).unwrap()
        .with("lot", 42.0).unwrap();
    ClosedGasParcel::new(model, point, GAS_VOLUME).unwrap()
}

fn wall(young: f64, slope: f64, alpha: f64, rho: f64, conductivity: f64)
    -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, value, derivative) in [
        ("young_modulus", Pressure::DIMS, young, slope),
        ("poisson_ratio", Dims([0; 6]), 0.0, 0.0),
        ("linear_thermal_expansion_coefficient", INVERSE_TEMPERATURE_DIMS, alpha, 0.0),
        ("density", Density::DIMS, rho, 0.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, CP, 0.0),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, conductivity, 0.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * derivative), (500.0, value + 200.0 * derivative)], dims },
            validity: ValidityDomain::unconstrained().with_quantity("T", quantity(), 280.0, 500.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance {
                source: "synthetic deformable cavity material".into(), license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic cavity wall".into(), phase: "solid".into(),
        process: "numerical fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(), SolidLiquidPhase::Solid, [280.0, 500.0].map(|t| EnthalpyPhaseKnot {
            temperature_k: t, specific_enthalpy_j_kg: CP * (t - 280.0),
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: rho,
        }).to_vec(),
    ).unwrap();
    (card, curve)
}

fn heat(curve: &EquilibriumEnthalpyPhaseCurve, rho: f64, reservoir: bool) -> DepositionSimulation<'_> {
    DepositionSimulation::new(vec![TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
        id: IDS[0], segment_id: 0, body_id: 3, material_id: MaterialId(0),
        mass_kg: rho * VOLUME, specific_h_j_kg: 20.0 * CP,
        geometry: CellGeometry { start_m: [0.0; 3], end_m: [1.0, 0.0, 0.0],
            volume_m3: VOLUME, surface_area_m2: 3.0 },
        boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, source_w: 0.0,
            reservoir: reservoir.then_some(Reservoir {
                area_m2: 0.1, conductance_w_k: 0.0, temperature_k: 400.0,
            }) },
    }) }], &[(MaterialId(0), curve as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 1, max_links: 0, max_events: 1, max_substeps: 1000,
    }).unwrap().advance_to(0.0, thermal_controls(), || Ok(())).unwrap().0
}
fn contacts(k_area: f64) -> [SolidGasContact; 1] {
    [SolidGasContact { id: 17, cell: IDS[0],
        geometry: GasConductionGeometry { area_m2: k_area, gas_distance_m: 0.001, cell_distance_m: 0.02 },
        conductivity: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        selection: MaterialPropertySelection::SingleClaimOnly }]
}
fn controls() -> CavitySolve {
    CavitySolve {
        gas: GasEnergySolve { minimum_temperature_k: 280.0, maximum_temperature_k: 400.0,
            energy_tolerance_j_kg: 1e-8, maximum_iterations: 128 },
        solid: solve(), thermal: StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-10 },
        maximum_trials: 400, maximum_newton_iterations: 20, difference_step_j: 1e-5,
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-10, relative_tolerance: 1e-12,
            linear_restart: 4, forcing_minimum: 1e-4, ..NewtonKrylovConfig::default() },
        volume_interval_m3: [0.04, 0.06], volume_tolerance_m3: 1e-12, energy_tolerance_j: 2e-9,
    }
}
fn launch<'s, 'm>(
    source: &'s SolidEnthalpy<'m>, thermal: &DepositionSimulation<'_>,
    gas: ClosedGasParcel, speed: f64, reference_force: f64, cx: &fs_exec::Cx<'_>,
) -> ThermoelasticCavity<'s, 'm> {
    // The sloping x+y+z=1 face has unit-pressure coefficient -1/6
    // on the sole free x displacement at node 1.
    let face = source.geometry().faces().iter().position(|face|
        face.nodes.contains(&1) && face.nodes.contains(&2) && face.nodes.contains(&3)).unwrap();
    let faces = [face];
    let reference = [(3, reference_force)];
    let loads = source.geometry().assemble(
        &[SurfaceLoad3::Pressure { faces: &faces, pascals: gas.state().state().pressure - PRESSURE }],
        VolumeLoad3::None, &reference, cx,
    ).unwrap();
    let body = ThermoelasticSolid::try_new(source, thermal, &loads, &fixed(),
        &[[280.0, 360.0]], mechanical(), solve(), cx).unwrap();
    body.begin_cavity(thermal, gas, &faces, PRESSURE, &reference,
        &contacts(PATCH_AREA), &velocity(speed), dynamic_controls(), cx).unwrap()
}

// Independent analytic scalar Newmark plus constant-capacity ideal-gas equation.
// NASA source has Cp=3.5 R, constant k_g=0.03; consistent tetra mass is rho V/10.
// Neither the native material solver nor its pressure/volume helpers enter here.
fn reference(speed: f64, reference_force: f64, young: f64, rho: f64, kb: f64,
    initial_gas_t: f64, dt: f64, external_heat: f64) -> [f64; 5] {
    let b = -1.0 / 6.0;
    let mass = rho * VOLUME / 10.0;
    let k = young * VOLUME;
    let u0 = reference_force / k;
    let capacity_g = 2.5 * PRESSURE * GAS_VOLUME / initial_gas_t;
    let g = PATCH_AREA / (0.001 / 0.03 + 0.02 / kb);
    let a = 2.0 * mass / dt.powi(2) + 0.5 * k;
    let rhs = (2.0 * mass / dt.powi(2) - 0.5 * k) * u0
        + 2.0 * mass / dt * speed + reference_force;
    let (mut t, mut u) = (initial_gas_t, u0);
    for _ in 0..30 {
        let dv = b * (u - u0);
        let volume = GAS_VOLUME + dv;
        let p = PRESSURE * t / initial_gas_t * GAS_VOLUME / volume;
        let r1 = a * u - rhs - 0.5 * b * (p - PRESSURE);
        let r2 = capacity_g * (t - initial_gas_t)
            + 0.5 * (PRESSURE + p) * dv + dt * g * (t - 300.0) - external_heat;
        if r1.abs() < 1e-11 && r2.abs() < 1e-12 {
            let tw = 300.0 + dt * g * (t - 300.0) / (rho * VOLUME * CP);
            return [t, u, 2.0 * (u - u0) / dt - speed, tw, p];
        }
        let j11 = a + 0.5 * b * b * p / volume;
        let j12 = -0.5 * b * p / t;
        let j21 = 0.5 * b * (PRESSURE + p) - 0.5 * p * b * dv / volume;
        let j22 = capacity_g + 0.5 * p * dv / t + dt * g;
        let determinant = j11 * j22 - j12 * j21;
        u -= (j22 * r1 - j12 * r2) / determinant;
        t -= (j11 * r2 - j21 * r1) / determinant;
    }
    panic!("independent gas/elastic reference failed to converge");
}

#[test]
fn native_pressure_work_and_film_heat_match_independent_deformable_cavity() {
    with_cx(|cx| {
        for (speed, formation) in [(0.02, 0.0), (-0.02, 0.0), (0.02, 1e15)] {
            let (card, curve) = wall(10000.0, 0.0, 0.0, 6.0, 1.0);
            let thermal = heat(&curve, 6.0, false);
            let context = context(300.0);
            let source = SolidEnthalpy::try_new(&thermal, &NODES, &TETS,
                vec![source(30, &card, &curve, &context)], TetAssemblyBudget::standard(),
                0.05, mechanical().mesh_frame, cx).unwrap();
            let cavity = launch(&source, &thermal, gas(340.0, formation), speed, 0.3, cx);
            let before = cavity.solid().nodal_motion(cx).unwrap();
            let gas_before = cavity.gas().clone();
            let expected = reference(speed, 0.3, 10000.0, 6.0, 1.0, 340.0, 0.005, 0.01);
            let next = cavity.advance(&thermal, &(), 0.005, 0.01, controls(), [0.0; 2], cx).unwrap();
            let actual = next.cavity.solid().nodal_motion(cx).unwrap();
            near(next.cavity.gas().state().state().temperature, expected[0], 1e-7);
            near(actual[0][1][0], expected[1], 2e-11);
            near(actual[1][1][0], expected[2], 1e-8);
            near(next.thermal.temperature_k(IDS[0]).unwrap(), expected[3], 1e-8);
            near(next.cavity.gas().state().state().pressure, expected[4], 1e-6);
            let dv = -(expected[1] - before[0][1][0]) / 6.0;
            near(next.balance.swept_volume_m3, dv, 1e-12);
            near(next.balance.gas_pressure_work_j, 0.5 * (PRESSURE + expected[4]) * dv, 2e-9);
            near(next.balance.exterior_pressure_work_j, PRESSURE * dv, 2e-9);
            near(next.balance.reference_load_work_j, 0.3 * (expected[1] - before[0][1][0]), 1e-11);
            assert!(next.balance.gas_to_solid_heat_j > 0.0);
            assert!(next.balance.pressure_work_residual_j.abs() < 2e-9);
            assert!(next.balance.combined_residual_j.abs() < 2e-9);
            assert_eq!(next.cavity.gas().mass_kg().to_bits(), gas_before.mass_kg().to_bits());
            assert_eq!(next.cavity.gas().point().axes().get("lot"), Some(&42.0));
            assert_eq!(next.cavity.solid().step(), 1);
            assert_eq!(cavity.solid().step(), 0);
            assert_eq!(cavity.solid().nodal_motion(cx).unwrap(), before);
            assert_eq!(*cavity.gas(), gas_before);
            assert_eq!(thermal.time_s(), 0.0);
        }
    });
}

#[test]
fn wall_density_conductivity_and_stiffness_data_change_continued_native_motion() {
    with_cx(|cx| {
        let mut observed = Vec::new();
        for (young, rho, conductivity) in [
            (10000.0, 6.0, 1.0), (20000.0, 6.0, 1.0),
            (10000.0, 12.0, 1.0), (10000.0, 6.0, 2.0),
        ] {
            let (card, curve) = wall(young, 20.0, 0.0001, rho, conductivity);
            let mut thermal = heat(&curve, rho, false);
            let context = context(300.0);
            let source = SolidEnthalpy::try_new(&thermal, &NODES, &TETS,
                vec![source(30, &card, &curve, &context)], TetAssemblyBudget::standard(),
                0.05, mechanical().mesh_frame, cx).unwrap();
            let mut cavity = launch(&source, &thermal, gas(340.0, 0.0), 0.02, 0.3, cx);
            for step in 1..=8 {
                let next = cavity.advance(&thermal, &(), step as f64 * 0.002,
                    0.001, controls(), [0.0; 2], cx).unwrap();
                assert_eq!(next.cavity.solid().step(), step);
                assert_eq!(next.cavity.time_s(), next.thermal.time_s());
                assert!(next.balance.combined_residual_j.abs() < 2e-9);
                near(next.cavity.gas().volume_m3(), next.balance.geometric_volume_m3, 1e-12);
                cavity = next.cavity;
                thermal = next.thermal;
            }
            let accepted = cavity.solid().nodal_motion(cx).unwrap();
            observed.push([accepted[0][1][0], cavity.gas().state().state().pressure,
                thermal.temperature_k(IDS[0]).unwrap()]);
            let bad = CavitySolve { volume_interval_m3: [0.1, 0.2], ..controls() };
            assert!(cavity.advance(&thermal, &(), 0.02, 0.001, bad, [0.0; 2], cx).is_err());
            assert_eq!(cavity.solid().nodal_motion(cx).unwrap(), accepted);
            assert_eq!(cavity.solid().step(), 8);
        }
        for other in &observed[1..] {
            assert!((other[0] - observed[0][0]).abs() > 1e-10);
            assert!((other[1] - observed[0][1]).abs() > 1e-7);
        }
        assert!((observed[3][2] - observed[0][2]).abs() > 1e-5);
    });
}

#[derive(Clone)]
struct ReservoirObserver { steps: usize, reservoir_calls: usize, fail_at: Option<usize> }
impl DepositionObserver for ReservoirObserver {
    type Error = DepositionError;
    fn reservoir_conductances(&mut self, _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        self.reservoir_calls += 1;
        Ok(BTreeMap::from([(IDS[0], 0.2)]))
    }
    fn requires_step_interval(&self) -> bool { true }
    fn step_sources(&mut self, _: &DepositionSimulation<'_>, _: f64, _: f64, _: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        Ok(BTreeMap::from([(IDS[0], 0.4)]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>, _: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        check()?;
        self.steps += 1;
        if self.fail_at == Some(self.steps) { return Err(DepositionError::InvalidInput("late observer refusal")); }
        Ok(())
    }
}

#[test]
fn existing_reservoir_and_interval_heat_hooks_survive_cavity_and_rollback() {
    with_cx(|cx| {
        let (card, curve) = wall(10000.0, 0.0, 0.0, 6.0, 1.0);
        let thermal = heat(&curve, 6.0, true);
        let context = context(300.0);
        let source = SolidEnthalpy::try_new(&thermal, &NODES, &TETS,
            vec![source(30, &card, &curve, &context)], TetAssemblyBudget::standard(),
            0.05, mechanical().mesh_frame, cx).unwrap();
        let cavity = launch(&source, &thermal, gas(300.0, 0.0), 0.0, 0.0, cx);
        let observer = ReservoirObserver { steps: 0, reservoir_calls: 0, fail_at: None };
        let next = cavity.advance(&thermal, &observer, 0.005, 0.0, controls(), [0.0; 2], cx).unwrap();
        near(next.thermal.temperature_k(IDS[0]).unwrap(), 300.51, 1e-8);
        near(next.balance.external_solid_heat_j, 0.102, 1e-10);
        near(next.balance.gas_to_solid_heat_j, 0.0, 1e-10);
        assert_eq!(next.observer.steps, 1);
        assert_eq!(next.observer.reservoir_calls, 1);
        assert_eq!(observer.steps, 0);
        assert_eq!(observer.reservoir_calls, 0);
        let before = cavity.solid().nodal_motion(cx).unwrap();
        let failing = ReservoirObserver { fail_at: Some(1), ..observer.clone() };
        assert!(cavity.advance(&thermal, &failing, 0.005, 0.0, controls(), [0.0; 2], cx).is_err());
        let budget = CavitySolve { maximum_trials: 1, ..controls() };
        assert!(cavity.advance(&thermal, &observer, 0.005, 0.0, budget, [0.0; 2], cx).is_err());
        assert_eq!(cavity.solid().nodal_motion(cx).unwrap(), before);
        assert_eq!(cavity.solid().step(), 0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(failing.steps, 0);
    });
}
