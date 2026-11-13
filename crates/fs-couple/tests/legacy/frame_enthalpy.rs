//! G0/G1: native heated frame histories, independent scalar work and real heat cells.
//! Manufactured coefficients exercise the coupled equations, not alloy calibration.

use fs_couple::frame_enthalpy::{
    FrameEnthalpyConfig, FrameEnthalpyError, FrameFiberCell, FrameThermalCell,
    FrameThermalMechanics, ThermoplasticFrame,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::creep::CreepTolerances;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::source::{HardeningPropertyBasis, KINEMATIC_HARDENING_MODULUS_PROPERTY};
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS, CAUCHY_YIELD_CURVE_PROPERTY,
};
use fs_material::plastic::thermal::tabulated::ThermalTabulatedSource;
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, Section};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::frame2d::creep::FrameCreepTolerances;
use fs_solid::{Frame2d, Geometry, SectionModel};
use fs_thermal::enthalpy::*;
use std::collections::BTreeMap;

const E: f64 = 1.2e9;
const DE: f64 = -5e6;
const HK: f64 = 1.6e8;
const DHK: f64 = -1e6;
const EXPANSION: f64 = 2e-5;
const CP: f64 = 500.0;
const RHO: f64 = 1000.0;
const AREA: f64 = 5e-5;
const ARM: f64 = 0.3;
const V: f64 = 2.0 * AREA;
const STRESS: f64 = 4.25e6;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const KNOTS: [(f64, f64); 5] = [
    (0.0, 3e6),
    (0.002, 4e6),
    (0.006, 4e6),
    (0.012, 7e6),
    (0.04, 8.4e6),
];
const SLOPES: [f64; 5] = [-10e3, -20e3, -20e3, -50e3, -80e3];

fn tq() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", tq(), t).unwrap()
}
fn claim(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if slope == 0.0 {
            PropertyValue::Scalar { value, dims }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * slope), (340.0, value + 40.0 * slope)],
                dims,
            }
        },
        validity: ValidityDomain::unconstrained().with_quantity("T", tq(), 280.0, 340.0),
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic coupled thermal frame fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new(revision: u32) -> Self {
        let mut claims = ClaimSet::new();
        for (name, dims, value, slope) in [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E, DE),
            (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2, 0.0),
            (YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6, -10e3),
            (BASIS.property_name(), Pressure::DIMS, 0.0, 0.0),
            (
                KINEMATIC_HARDENING_MODULUS_PROPERTY,
                Pressure::DIMS,
                HK,
                DHK,
            ),
            ("specific_heat_capacity", CP_DIMS, CP, 0.0),
            (DENSITY_PROPERTY, Density::DIMS, RHO, 0.0),
            (
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                EXPANSION,
                0.0,
            ),
            (
                OVERSTRESS_VISCOSITY_PROPERTY,
                OVERSTRESS_VISCOSITY_DIMS,
                1e7,
                -50e3,
            ),
            (OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 2e6, 0.0),
            (OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, 2.0, 0.0),
        ] {
            claims
                .insert_claim(claim(name, dims, value, slope))
                .unwrap();
        }
        let mut curve = claim(CAUCHY_YIELD_CURVE_PROPERTY, Pressure::DIMS, 3e6, 0.0);
        curve.value = PropertyValue::Surface2D {
            axes: [ACCUMULATED_PLASTIC_STRAIN_AXIS.into(), "T".into()],
            axis_dims: [Dims::NONE, Temperature::DIMS],
            samples: [280.0, 300.0, 340.0]
                .into_iter()
                .flat_map(|t| {
                    KNOTS
                        .into_iter()
                        .zip(SLOPES)
                        .map(move |((p, y), slope)| ([p, t], y + slope * (t - 300.0)))
                })
                .collect(),
            triangles: vec![],
            dims: Pressure::DIMS,
        };
        curve.validity = curve.validity.with_quantity(
            ACCUMULATED_PLASTIC_STRAIN_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            0.04,
        );
        curve.interpolation = InterpolationPolicy::TabulatedOnly;
        claims.insert_claim(curve).unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "thermal frame test".into(),
                phase: "solid".into(),
                process: "manufactured source".into(),
                revision,
            },
            claims,
            vec![],
        )
        .unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &point(280.0),
                upper: &point(340.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * CP,
            1e-6,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn source(&self, rate: bool) -> ThermalTabulatedSource<'_> {
        let base = ThermoplasticSource::try_new(
            SourceQuery {
                card: &self.card,
                caloric: &self.caloric,
                context: point(300.0),
                temperature_axis: "T".into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
                hardening_basis: BASIS,
                thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
                strain_limit: 0.04,
            },
            0.0,
        )
        .unwrap();
        let base = if rate {
            base.with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
                .unwrap()
        } else {
            base
        };
        base.with_tabulated_hardening(ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly), 15)
            .unwrap()
    }
    fn raw_frame(&self, axial: bool, rate: bool, mass: FrameMassFormulation) -> Frame2d {
        let law = self
            .source(rate)
            .at_uniaxial_temperatures(300.0, 300.0, 1.0)
            .unwrap();
        let material = FiberMaterial::thermal_tabulated(law).unwrap();
        let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, mass);
        frame.add_node(0.0, 0.0);
        frame.add_node(1.0, 0.0);
        frame.fix(0, [true; 3]);
        if axial {
            frame.fix(1, [false, true, true]);
        }
        frame
            .add_element(0, 1, 3, RHO * V, &|| {
                SectionModel::Fiber(Section::new(
                    [-ARM, ARM]
                        .into_iter()
                        .map(|y| Fiber {
                            y,
                            area: AREA,
                            law: material.law().clone(),
                        })
                        .collect(),
                ))
            })
            .unwrap();
        frame
    }
    fn thermal(&self, cells: usize, mass_factor: f64) -> DepositionSimulation<'_> {
        let events = (0..cells)
            .map(|i| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64),
                    segment_id: i as u64,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: RHO * V / cells as f64 * mass_factor,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [0.0; 3],
                        end_m: [1.0, 0.0, 0.0],
                        volume_m3: V / cells as f64,
                        surface_area_m2: 0.02,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        source_w: 0.0,
                        reservoir: None,
                    },
                }),
            })
            .collect();
        DepositionSimulation::new(
            events,
            &[(MaterialId(1), self.caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 2,
                max_links: 4,
                max_events: 4,
                max_substeps: 1000,
            },
        )
        .unwrap()
        .advance_to(0.0, heat_config(), || Ok(()))
        .unwrap()
        .0
    }
    fn mappings(&self, cells: usize) -> Vec<FrameFiberCell> {
        (0..3)
            .flat_map(|section| {
                (0..2).map(move |fiber| FrameFiberCell {
                    element: 0,
                    section,
                    fiber,
                    cell: CellId(if cells == 1 { 0 } else { fiber as u64 }),
                })
            })
            .collect()
    }
    fn bind<'a>(
        &'a self,
        thermal: &DepositionSimulation<'_>,
        axial: bool,
        rate: bool,
        mass: FrameMassFormulation,
        cx: &Cx<'_>,
    ) -> ThermoplasticFrame<'a> {
        let cells = thermal.cells().len();
        ThermoplasticFrame::try_new(
            self.raw_frame(axial, rate, mass),
            thermal,
            (0..cells)
                .map(|i| FrameThermalCell {
                    cell: CellId(i as u64),
                    source: self.source(rate).into(),
                })
                .collect(),
            self.mappings(cells),
            6,
            cx,
        )
        .unwrap()
    }
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 0x4654484d,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn config() -> FrameEnthalpyConfig {
    FrameEnthalpyConfig {
        maximum_iterations: 64,
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-7,
        enthalpy_relaxation: 1.0,
    }
}
fn mechanical() -> FrameCreepTolerances {
    FrameCreepTolerances {
        material: CreepTolerances {
            stress_pa: 1e-3,
            energy_j_m3: 1e-5,
        },
        energy_j: 1e-8,
        force_n: 1e-5,
        moment_n_m: 1e-7,
        maximum_fiber_points: 6,
    }
}
fn heat_config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.01,
        energy_tolerance_j: 1e-9,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "actual={a:.17e}, expected={b:.17e}, tolerance={tol:e}"
    );
}
fn axial_load(stress: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, V * stress, 0.0, 0.0]
}

#[derive(Clone, Debug, Default)]
struct Heating {
    rates: BTreeMap<CellId, f64>,
    delivered_j: f64,
    steps: usize,
    fail: bool,
}
impl Heating {
    fn with_heat(mut self, heat: &[(CellId, f64)], dt: f64) -> Self {
        self.rates = heat.iter().map(|(id, q)| (*id, q / dt)).collect();
        self
    }
}
impl DepositionObserver for Heating {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        _: f64,
        _: f64,
        _: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        Ok(self.rates.clone())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        self.delivered_j += dt * self.rates.values().sum::<f64>();
        if self.fail {
            Err(DepositionError::InvalidInput("injected candidate failure"))
        } else {
            Ok(())
        }
    }
}

// Independent stress-controlled free-sided first law, including the full
// original elastic/Prager/table internal energy and native endpoint loss.
#[derive(Clone, Copy, Debug)]
struct Point {
    t: f64,
    stress: f64,
    p: f64,
}
impl Point {
    fn initial() -> Self {
        Self {
            t: 300.0,
            stress: 0.0,
            p: 0.0,
        }
    }
    fn strain(self) -> f64 {
        self.stress / youngs(self.t) + self.p + EXPANSION * (self.t - 300.0)
    }
}
fn youngs(t: f64) -> f64 {
    E + DE * (t - 300.0)
}
fn hk(t: f64) -> f64 {
    HK + DHK * (t - 300.0)
}
fn ordinates(t: f64) -> [f64; 5] {
    core::array::from_fn(|i| KNOTS[i].1 + SLOPES[i] * (t - 300.0))
}
fn value(p: f64, values: [f64; 5]) -> f64 {
    let i = KNOTS.windows(2).position(|w| p <= w[1].0).unwrap();
    values[i] + (values[i + 1] - values[i]) * (p - KNOTS[i].0) / (KNOTS[i + 1].0 - KNOTS[i].0)
}
fn primitive(p: f64, values: [f64; 5]) -> f64 {
    let mut total = 0.0;
    for i in 0..KNOTS.len() - 1 {
        if p <= KNOTS[i].0 {
            break;
        }
        let width = p.min(KNOTS[i + 1].0) - KNOTS[i].0;
        total += width
            * (values[i] - values[0]
                + 0.5 * width * (values[i + 1] - values[i]) / (KNOTS[i + 1].0 - KNOTS[i].0));
    }
    total
}
fn energy(point: Point) -> (f64, f64) {
    let elastic = point.stress * point.stress / (2.0 * youngs(point.t));
    let psi = elastic
        + 0.5 * hk(point.t) * point.p * point.p
        + primitive(point.p.abs(), ordinates(point.t));
    let entropy = EXPANSION * point.stress
        - DE / youngs(point.t) * elastic
        - 0.5 * DHK * point.p * point.p
        - primitive(point.p.abs(), SLOPES);
    (psi, entropy)
}
fn root(mut lower: f64, mut upper: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lower) <= 0.0 && residual(upper) >= 0.0);
    for _ in 0..100 {
        let middle = lower + 0.5 * (upper - lower);
        if middle == lower || middle == upper {
            break;
        }
        if residual(middle) > 0.0 {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    lower + 0.5 * (upper - lower)
}
fn stress_endpoint(old: Point, stress: f64, temperature: f64) -> Point {
    let p0 = old.p.abs();
    let drive = stress.abs() - hk(temperature) * p0 - value(p0, ordinates(temperature));
    let p = if stress != 0.0 && drive > 0.0 {
        stress.signum()
            * root(p0, 0.02, |a| {
                hk(temperature) * a + value(a, ordinates(temperature)) - stress.abs()
            })
    } else {
        old.p
    };
    Point {
        t: temperature,
        stress,
        p,
    }
}
fn numerical(old: Point, next: Point) -> f64 {
    let dp = next.p.abs() - old.p.abs();
    let elastic_change = next.stress / youngs(next.t) - old.stress / youngs(old.t);
    let y = ordinates(next.t);
    0.5 * youngs(next.t) * elastic_change.powi(2)
        + 0.5 * hk(next.t) * (next.p - old.p).powi(2)
        + (value(next.p.abs(), y) - y[0]) * dp
        - (primitive(next.p.abs(), y) - primitive(old.p.abs(), y))
}
fn reference(old: Point, stress: f64, external_j: f64, volume: f64) -> Point {
    let old_energy = energy(old);
    let residual = |t| {
        let next = stress_endpoint(old, stress, t);
        let end_energy = energy(next);
        RHO * CP * (t - old.t) + end_energy.0 + t * end_energy.1
            - old_energy.0
            - old.t * old_energy.1
            + numerical(old, next)
            - stress * (next.strain() - old.strain())
            - external_j / volume
    };
    stress_endpoint(old, stress, root(280.0, 339.9, residual))
}
fn check_point(
    frame: &ThermoplasticFrame<'_>,
    thermal: &DepositionSimulation<'_>,
    fiber: usize,
    cell: CellId,
    expected: Point,
) {
    let t = 300.0 + thermal.cells()[&cell].specific_h_j_kg / CP;
    near(t, expected.t, 5e-8);
    for section in 0..3 {
        let (_, state) = frame
            .frame()
            .thermal_tabulated_fiber(0, section, fiber)
            .unwrap();
        near(state.mechanical.plastic_strain(), expected.p, 1e-10);
        near(
            state.mechanical.equivalent_plastic_strain(),
            expected.p.abs(),
            1e-10,
        );
        assert_eq!(
            state.temperature_k(),
            frame
                .frame()
                .thermal_tabulated_fiber(0, section, fiber)
                .unwrap()
                .0
                .temperature_k()
        );
        near(state.temperature_k(), t, 1e-12);
    }
}

#[test]
fn loaded_member_follows_actual_heating_cooling_and_release_with_retained_shape() {
    let fixture = Fixture::new(0);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut thermal = fixture.thermal(1, 1.0);
        let mut frame = fixture.bind(
            &thermal,
            true,
            false,
            FrameMassFormulation::ConsistentReference,
            cx,
        );
        let mut observer = Heating::default();
        let mut expected = Point::initial();
        let mut target = 0.0;
        for (dt, stress, q) in [
            (0.13, STRESS, 0.0),
            (0.27, STRESS, RHO * V * CP * 35.0),
            (0.19, STRESS, -RHO * V * CP * 20.0),
            (0.11, 0.0, 0.0),
        ] {
            target += dt;
            let old = expected;
            expected = reference(old, stress, q, V);
            let before_u = frame.frame().displacements().to_vec();
            let before_t = thermal.time_s();
            let before_energy = observer.delivered_j;
            observer = observer.with_heat(&[(CellId(0), q)], dt);
            let step = frame
                .advance(
                    &thermal,
                    &observer,
                    target,
                    &axial_load(stress),
                    FrameThermalMechanics::Quasistatic,
                    mechanical(),
                    heat_config(),
                    config(),
                    cx,
                )
                .unwrap();
            assert_eq!(frame.frame().displacements(), before_u);
            assert_eq!(thermal.time_s(), before_t);
            check_point(step.frame(), step.thermal(), 0, CellId(0), expected);
            check_point(step.frame(), step.thermal(), 1, CellId(0), expected);
            near(
                step.frame().frame().displacements()[3],
                expected.strain(),
                1e-10,
            );
            let report = step.report();
            near(report.external_heat_j, q, 1e-7);
            near(
                report.mechanical_energy.applied_work_j,
                V * stress * (expected.strain() - old.strain()),
                1e-7,
            );
            near(
                report.mechanical_energy.numerical_dissipation_j,
                V * numerical(old, expected),
                1e-7,
            );
            near(
                report.cell_heat_j[&CellId(0)],
                RHO * V * CP * (expected.t - old.t) - q,
                1e-7,
            );
            near(report.thermal_residual_j, 0.0, 1e-7);
            near(report.combined_residual_j, 0.0, 1e-7);
            near(step.observer().delivered_j - before_energy, q, 1e-7);
            assert!(report.iterations > 1);
            let (next_t, next_f, next_o, _) = step.into_parts();
            thermal = next_t;
            frame = next_f;
            observer = next_o;
            if q > 0.0 {
                assert!(expected.p > 0.006);
            }
        }
        assert!(expected.p > 0.006);
        near(
            frame.frame().displacements()[3],
            expected.p + EXPANSION * (expected.t - 300.0),
            1e-10,
        );
        assert!(
            frame
                .frame()
                .element_forces(0)
                .iter()
                .all(|value| value.abs() < 1e-5)
        );
    });
}

#[test]
fn nonuniform_fiber_temperatures_change_loaded_curvature_through_original_cells() {
    let fixture = Fixture::new(0);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let mut thermal = fixture.thermal(2, 1.0);
        let mut frame = fixture.bind(
            &thermal,
            false,
            false,
            FrameMassFormulation::ConsistentReference,
            cx,
        );
        let mut observer = Heating::default();
        let mut expected = [Point::initial(); 2];
        let mut target = 0.0;
        let moment = V * ARM * STRESS;
        let load = [0.0, 0.0, 0.0, 0.0, 0.0, moment];
        for (dt, heat) in [(0.1, 0.0), (0.23, RHO * AREA * CP * 35.0)] {
            target += dt;
            expected[0] = reference(expected[0], STRESS, heat, AREA);
            expected[1] = reference(expected[1], -STRESS, 0.0, AREA);
            observer = observer.with_heat(&[(CellId(0), heat)], dt);
            let step = frame
                .advance(
                    &thermal,
                    &observer,
                    target,
                    &load,
                    FrameThermalMechanics::Quasistatic,
                    mechanical(),
                    heat_config(),
                    config(),
                    cx,
                )
                .unwrap();
            for f in 0..2 {
                check_point(
                    step.frame(),
                    step.thermal(),
                    f,
                    CellId(f as u64),
                    expected[f],
                );
            }
            let extension = 0.5 * (expected[0].strain() + expected[1].strain());
            let curvature = (expected[0].strain() - expected[1].strain()) / (2.0 * ARM);
            assert!(curvature.abs() < 0.04);
            for &(eps, kappa) in step.frame().frame().section_deformations(0) {
                near(eps, extension, 1e-10);
                near(kappa, curvature, 1e-8);
            }
            near(step.frame().frame().displacements()[3], extension, 1e-10);
            near(
                step.frame().frame().displacements()[4],
                0.5 * curvature,
                1e-8,
            );
            near(step.frame().frame().displacements()[5], curvature, 1e-8);
            near(
                step.report().cell_heat_j.values().sum(),
                step.report().mechanical_energy.caloric_heat_j(),
                1e-7,
            );
            near(step.report().combined_residual_j, 0.0, 1e-7);
            let (t, f, o, _) = step.into_parts();
            thermal = t;
            frame = f;
            observer = o;
        }
        assert!(expected[0].p > 0.006 && expected[1].p.abs() < 0.002);
        assert!(frame.frame().displacements()[3] > 0.002);
    });
}

#[test]
fn coupled_rate_motion_keeps_numerical_inertia_out_of_heat_for_both_mass_forms() {
    let fixture = Fixture::new(0);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        for mass in [
            FrameMassFormulation::Lumped,
            FrameMassFormulation::ConsistentReference,
        ] {
            let mut thermal = fixture.thermal(1, 1.0);
            let mut frame = fixture.bind(&thermal, true, true, mass, cx);
            let mut observer = Heating::default();
            let mut target = 0.0;
            let free_mass = RHO
                * V
                * if mass == FrameMassFormulation::ConsistentReference {
                    1.0 / 3.0
                } else {
                    0.5
                };
            for (dt, q) in [(0.0011, RHO * V * CP * 5.0), (0.0017, RHO * V * CP * 2.0)] {
                target += dt;
                let actual_dt = target - frame.time_s();
                let old_v = frame.frame().velocities()[3];
                let old_h = thermal.cells()[&CellId(0)].specific_h_j_kg;
                observer = observer.with_heat(&[(CellId(0), q)], actual_dt);
                let step = frame
                    .advance(
                        &thermal,
                        &observer,
                        target,
                        &axial_load(STRESS),
                        FrameThermalMechanics::BackwardEuler,
                        mechanical(),
                        heat_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let dynamic = step.report().dynamic_energy.unwrap();
                let v = step.frame().frame().velocities()[3];
                near(
                    dynamic.kinetic_change_j,
                    0.5 * free_mass * (v * v - old_v * old_v),
                    1e-8,
                );
                near(
                    dynamic.inertial_numerical_dissipation_j,
                    0.5 * free_mass * (v - old_v).powi(2),
                    1e-8,
                );
                assert!(dynamic.inertial_numerical_dissipation_j > 0.0);
                assert!(dynamic.material.overstress_dissipation_j > 0.0);
                let caloric =
                    RHO * V * (step.thermal().cells()[&CellId(0)].specific_h_j_kg - old_h);
                near(caloric - q, dynamic.material.caloric_heat_j(), 1e-7);
                near(step.report().combined_residual_j, 0.0, 1e-7);
                let t = fixture
                    .caloric
                    .curve()
                    .state_at_specific_enthalpy(step.thermal().cells()[&CellId(0)].specific_h_j_kg)
                    .unwrap()
                    .temperature_k();
                for s in 0..3 {
                    for f in 0..2 {
                        let (law, state) = step
                            .frame()
                            .frame()
                            .thermal_tabulated_fiber(0, s, f)
                            .unwrap();
                        assert_eq!(state.temperature_k(), t);
                        assert_eq!(law.temperature_k(), t);
                        assert_eq!(law.duration_s(), actual_dt);
                    }
                }
                let (t, f, o, _) = step.into_parts();
                thermal = t;
                frame = f;
                observer = o;
            }
        }
    });
}

#[test]
fn mismatched_inventory_and_failed_feedback_cannot_change_accepted_histories() {
    let fixture = Fixture::new(0);
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let wrong_mass = fixture.thermal(1, 1.01);
        assert!(
            ThermoplasticFrame::try_new(
                fixture.raw_frame(true, false, FrameMassFormulation::ConsistentReference),
                &wrong_mass,
                vec![FrameThermalCell {
                    cell: CellId(0),
                    source: fixture.source(false).into()
                }],
                fixture.mappings(1),
                6,
                cx
            )
            .is_err()
        );
        let thermal = fixture.thermal(1, 1.0);
        let frame = fixture.bind(
            &thermal,
            true,
            false,
            FrameMassFormulation::ConsistentReference,
            cx,
        );
        let loaded = frame
            .advance(
                &thermal,
                &Heating::default(),
                0.1,
                &axial_load(STRESS),
                FrameThermalMechanics::Quasistatic,
                mechanical(),
                heat_config(),
                config(),
                cx,
            )
            .unwrap();
        let (thermal, frame, observer, _) = loaded.into_parts();
        let saved_u = frame.frame().displacements().to_vec();
        let saved_h = thermal.cells()[&CellId(0)].specific_h_j_kg;
        let saved_state = frame
            .frame()
            .thermal_tabulated_fiber(0, 0, 0)
            .unwrap()
            .1
            .clone();
        let heating = observer
            .clone()
            .with_heat(&[(CellId(0), RHO * V * CP * 10.0)], 0.2);
        let mut limited = config();
        limited.maximum_iterations = 1;
        assert!(matches!(
            frame.advance(
                &thermal,
                &heating,
                0.3,
                &axial_load(STRESS),
                FrameThermalMechanics::Quasistatic,
                mechanical(),
                heat_config(),
                limited,
                cx
            ),
            Err(FrameEnthalpyError::NotConverged { .. })
        ));
        let mut failed = heating.clone();
        failed.fail = true;
        assert!(matches!(
            frame.advance(
                &thermal,
                &failed,
                0.3,
                &axial_load(STRESS),
                FrameThermalMechanics::Quasistatic,
                mechanical(),
                heat_config(),
                config(),
                cx
            ),
            Err(FrameEnthalpyError::Observer(_))
        ));
        assert_eq!(frame.frame().displacements(), saved_u);
        assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, saved_h);
        assert_eq!(
            frame.frame().thermal_tabulated_fiber(0, 0, 0).unwrap().1,
            &saved_state
        );
        assert_eq!(observer.steps, heating.steps);
        let retry = frame
            .advance(
                &thermal,
                &heating,
                0.3,
                &axial_load(STRESS),
                FrameThermalMechanics::Quasistatic,
                mechanical(),
                heat_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.frame().time_s(), 0.3);
        assert_eq!(retry.thermal().time_s(), 0.3);
        assert!(retry.frame().frame().displacements()[3] > saved_u[3]);
    });
}

#[path = "frame_enthalpy/chaboche.rs"]
mod chaboche;
