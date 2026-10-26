//! G1/G3: independent energy reduction and continuing two-temperature motion.
//! These deliberately simple material sources are synthetic, not instrument data.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section_temperature::ConcentricSectionTemperatureLaw;
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, SPECIFIC_HEAT_CAPACITY_DIMS,
    THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS};
use fs_qty::{Density, Dims, DynViscosity, Pressure, Temperature};
use fs_thermal::enthalpy::{Boundary, CellGeometry, DepositedCell, DepositionBudget,
    Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};

const LENGTH: f64 = 0.502;
const REFERENCE: f64 = 0.5;
const RATE: u32 = 8000;
const RADII: [f64; 2] = [0.001, 0.0015];
const VISCOSITY: [f64; 2] = [1e6, 3e6];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, fs_exec::StreamKey {
            seed: 9, kernel_id: 11, tile: 0, iteration: 0,
        }, fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic);
        f(&cx, &gate)
    })
}

struct Sources {
    cards: Vec<MaterialCard>,
    caloric: Vec<ResolvedHeatCapacityCurve>,
    reference: ConcentricElasticSection,
}
impl Sources {
    fn new(outer_alpha: f64) -> Self {
        Self::with_resistivity(outer_alpha, 1.0)
    }
    fn with_resistivity(outer_alpha: f64, outer_resistivity_scale: f64) -> Self {
        let mut cards = Vec::new();
        let mut requirements = Vec::new();
        for i in 0..2 {
            let mut claims = ClaimSet::new();
            for (name, dims, value) in [
                (DENSITY_PROPERTY, Density::DIMS, scalar((i+1) as f64*1000.0, Density::DIMS)),
                ("young_modulus", Pressure::DIMS, PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, 2.02e9), (300.0, 2e9), (400.0, 1.9e9)], dims: Pressure::DIMS }),
                (ALPHA, INVERSE_TEMPERATURE_DIMS,
                    scalar(if i == 0 { 2e-5 } else { outer_alpha }, INVERSE_TEMPERATURE_DIMS)),
                (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS,
                    scalar(VISCOSITY[i], DynViscosity::DIMS)),
                (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS,
                    scalar(10.0, SPECIFIC_HEAT_CAPACITY_DIMS)),
                (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS,
                    scalar(20.0, THERMAL_CONDUCTIVITY_DIMS)),
                (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS,
                    PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                        knots: [280.0, 300.0, 400.0].into_iter().map(|t| (t,
                            (i+1) as f64*1e-6*(if i == 1 { outer_resistivity_scale } else { 1.0 })
                                *(1.0+0.002*(t-300.0)))).collect(), dims: ELECTRICAL_RESISTIVITY_DIMS }),
            ] {
                let interpolation = if matches!(value, PropertyValue::Curve { .. }) {
                    InterpolationPolicy::LinearInside
                } else { InterpolationPolicy::ConstantWithinValidity };
                claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims), value,
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 400.0)
                        .with("omega", 1.0, 20000.0),
                    interpolation, uncertainty: UncertaintyModel::Unstated, observations: vec![],
                    provenance: Provenance { source: "synthetic continuing layered string".into(),
                        license: "test-only".into(), artifact: None },
                }).unwrap();
                if i == 0 { requirements.push(ScalarPropertyRequirement::try_new(
                    name, dims, ScalarAdmissibility::NonNegative).unwrap()); }
            }
            cards.push(MaterialCard::assemble(MaterialStateId { chemistry: format!("layer-{i}"),
                phase: "solid".into(), process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap());
        }
        let caloric = cards.iter().map(caloric).collect();
        let reference = ConcentricElasticSection::try_new(0.0,
            cards.iter().zip(RADII).map(|(card, radius)| ConcentricMaterialLayer {
                outer_radius_m: radius,
                material: resolve_material_state_point(card, &point(300.0), &requirements,
                    MaterialPropertySelection::SingleClaimOnly).unwrap(),
            }).collect()).unwrap();
        Self { cards, caloric, reference }
    }
    fn runtime<'a>(&'a self, cx: &Cx<'_>)
        -> (ThermalConcentricStringRuntime<'a>, DepositionSimulation<'a>)
    {
        let (runtime, thermal) = self.runtime_with(cx, template());
        (runtime.unwrap(), thermal)
    }
    fn runtime_with<'a>(&'a self, cx: &Cx<'_>, template: PrestressedString)
        -> (Result<ThermalConcentricStringRuntime<'a>, ConcentricStringError>, DepositionSimulation<'a>)
    {
        let mut events = Vec::new();
        let mut materials: Vec<(MaterialId, &'a dyn ThermalLaw)> = Vec::new();
        for (i, (geometry, chart)) in self.reference.layer_integrals().iter().zip(&self.caloric).enumerate() {
            materials.push((MaterialId(i as u64), chart.curve()));
            events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
                id: CellId(i as u64), segment_id: 0, body_id: 0, material_id: MaterialId(i as u64),
                mass_kg: geometry.line_density_kg_m()*LENGTH,
                specific_h_j_kg: 200.0,
                geometry: CellGeometry { start_m: [0.0; 3], end_m: [0.0, 0.0, LENGTH],
                    volume_m3: geometry.area_m2()*LENGTH, surface_area_m2: 0.01 },
                boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
            }) });
        }
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 0, a: CellId(0), b: CellId(1), area_m2: 0.001, conductance_w_k: 0.5,
        }) });
        let thermal = DepositionSimulation::new(events, &materials, DepositionBudget {
            max_cells: 2, max_links: 1, max_events: 3, max_substeps: 128,
        }).unwrap().advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 0.0 }, || Ok(())).unwrap().0;
        let cards: Vec<_> = self.cards.iter().collect();
        let law = ConcentricThermoelasticLaw::try_new(ConcentricSectionTemperatureLaw::try_new(
            self.reference.clone(), &cards, &["T", "T"]).unwrap()).unwrap();
        let bindings = self.caloric.iter().enumerate().map(|(i, caloric)| ConcentricStringThermalCell {
            cell: CellId(i as u64), caloric,
        }).collect();
        let runtime = ThermalConcentricStringRuntime::try_new(cx, template, law, bindings,
            &thermal, REFERENCE, 0.02, Some(Pluck { station_frac: 0.4, height_m: 0.0005 }),
            super::super::gas_state(fs_scenario::AmbientGas::sea_level()).unwrap(), 1.0, RATE);
        (runtime, thermal)
    }
}

#[test]
fn a_tension_loss_factor_refuses_because_no_layer_owns_its_heat() {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        for eta in [1e-4, f64::NAN] {
            let lossy = PrestressedString { tension_loss_factor: Some(eta), ..template() };
            assert!(matches!(sources.runtime_with(cx, lossy).0,
                Err(ConcentricStringError::Invalid("a layered thermoelastic string has no owner for tension-loss heat"))),
                "eta_F {eta} was admitted");
        }
        // Zero is the absent factor: the same modes and states.
        let zero = PrestressedString { tension_loss_factor: Some(0.0), ..template() };
        let (runtime, _) = sources.runtime_with(cx, zero);
        let runtime = runtime.unwrap();
        let (base, _) = sources.runtime(cx);
        assert_eq!(runtime.modes(), base.modes());
        assert_eq!(runtime.states(), base.states());
    });
}
fn scalar(value: f64, dims: Dims) -> PropertyValue { PropertyValue::Scalar { value, dims } }
fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T", t).unwrap().with("omega", 1.0).unwrap() }
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS), temperature_axis: "T",
        lower: &point(280.0), upper: &point(400.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0.0, 1e-9).unwrap()
}
fn template() -> PrestressedString {
    PrestressedString { length_m: LENGTH, tension_n: 1.0, lin_density_kg_m: 1.0,
        axial_stiffness_n: 0.0, width_m: 0.001, n_modes: 2, damping_ratio: 0.0,
        rayleigh: None, bending_stiffness_n_m2: 0.0, kelvin_voigt_bending: None,
        relaxation_bending: None, polarization_detune: 0.0, moving_end: false,
        tension_loss_factor: None }
}
fn controls() -> ConcentricStringThermalSolve {
    ConcentricStringThermalSolve { maximum_cells: 2, maximum_thermal_trials: 256,
        maximum_newton_iterations: 16, difference_step_j: 1e-7,
        maximum_temperature_increment_k: 2.0, cell_energy_tolerance_j: 2e-11,
        thermal_energy_tolerance_j: 4e-11, total_energy_tolerance_j: 1e-10,
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-11, relative_tolerance: 0.0,
            linear_restart: 2, ..NewtonKrylovConfig::default() } }
}
fn heat(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients { sources_w: BTreeMap::from([(CellId(0), 16.0), (CellId(1), 48.0)]),
        conductances_w_k: BTreeMap::new(), payload: () })
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!((actual-expected).abs() <= relative*actual.abs().max(expected.abs()).max(1e-15),
        "{actual} != {expected}");
}

// Independent energy using fixture geometry/coefficients, not runtime helpers.
fn independent_internal(runtime: &ThermalConcentricStringRuntime<'_>, outer_alpha: f64) -> f64 {
    let pi = core::f64::consts::PI;
    let mut result: f64 = runtime.states().iter().zip(runtime.modes())
        .map(|(q, m)| 0.5*(q.velocity_m_sqrt_kg_per_s.powi(2)
            +(m.angular_frequency_rad_s*q.displacement_m_sqrt_kg).powi(2))).sum();
    for i in 0..2 {
        let inside = if i == 0 { 0.0_f64 } else { RADII[i-1] };
        let area = pi*(RADII[i].powi(2)-inside.powi(2));
        let moment = pi/4.0*(RADII[i].powi(4)-inside.powi(4));
        let t = runtime.source().layers[i].temperature_k;
        let alpha = if i == 0 { 2e-5 } else { outer_alpha };
        let young = 2e9-1e6*(t-300.0);
        let e = (LENGTH-REFERENCE)/REFERENCE-alpha*(t-300.0);
        let mut entropy = REFERENCE*area*(young*alpha*e+0.5e6*e*e);
        for (n, q) in runtime.states().iter().enumerate() {
            let k = (n+1) as f64*pi/LENGTH;
            let slope = (-1e6*moment*k.powi(4)+area*(-1e6*e-young*alpha)*k*k)
                /runtime.string().lin_density_kg_m;
            entropy -= 0.5*slope*q.displacement_m_sqrt_kg.powi(2);
        }
        result += 0.5*REFERENCE*area*young*e*e+t*entropy;
    }
    result
}

#[test]
fn g1_layer_heat_changes_live_tension_sound_and_each_thermodynamic_reservoir() {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        let (mut runtime, mut thermal) = sources.runtime(cx);
        let initial_frequency = runtime.modes()[0].angular_frequency_rad_s;
        let mut emitted = false;
        for sample in 0..32 {
            let before = independent_internal(&runtime, 4e-5);
            let before_h = thermal.energy_j();
            let mut frozen = runtime.member.clone();
            frozen.model.step(&[0.0; 2]).unwrap();
            let frame = runtime.step_with_thermal_network(cx, &mut thermal, controls(),
                &[0.002, 0.006], &[0.0; 2], heat).unwrap();
            assert_eq!(runtime.states(), frozen.model.states());
            assert_eq!(runtime.accepted_samples(), sample+1);
            assert_eq!(frame.vibration.epoch, sample+1);
            let after = independent_internal(&runtime, 4e-5);
            near(runtime.mechanical_internal_energy_j(cx).unwrap(), after, 5e-13);
            let loss = &frame.vibration.dissipation;
            let balance = frame.vibration.acoustic.input_work_j+frame.external_heat_j
                -(after-before+thermal.energy_j()-before_h+loss.air_loss_j
                    +loss.authored_loss_j+loss.roundoff_residual_j);
            assert!(balance.abs() <= controls().total_energy_tolerance_j+frame.roundoff_tolerance_j);
            let weights = runtime.section().layer_viscous_stiffness_n_m2_s();
            near(frame.layers[0].material_heat_j*weights[1], frame.layers[1].material_heat_j*weights[0], 5e-13);
            near(frame.layers.iter().map(|l| l.material_heat_j).sum(), loss.material_heat_j, 5e-13);
            assert!(frame.layers.iter().any(|l| l.mechanical_entropy_energy_change_j != 0.0));
            emitted |= frame.vibration.acoustic.observer_pressure_pa != 0.0;
            assert_eq!(runtime.source().section.section().layers()[0].material.card_identity(), sources.cards[0].content_hash());
        }
        assert!(emitted);
        assert!(thermal.temperature_k(CellId(1)).unwrap() > thermal.temperature_k(CellId(0)).unwrap());
        assert!(runtime.modes()[0].angular_frequency_rad_s < initial_frequency);
    });
}

#[test]
fn g3_one_layer_data_swap_changes_force_and_pressure_with_identical_solver_and_geometry() {
    with_cx(|cx, _| {
        let run = |alpha| {
            let source = Sources::new(alpha);
            let (mut runtime, mut thermal) = source.runtime(cx);
            let mut pressure = Vec::new();
            for _ in 0..16 {
                pressure.push(runtime.step_with_thermal_network(cx, &mut thermal, controls(),
                    &[0.002, 0.006], &[0.0; 2], heat).unwrap().vibration.acoustic.observer_pressure_pa);
            }
            (runtime.string().tension_n, runtime.modes()[0].angular_frequency_rad_s, pressure)
        };
        let low = run(0.0);
        let high = run(4e-5);
        assert!(high.0 < low.0 && high.1 < low.1);
        assert_ne!(high.2, low.2);
    });
}

#[test]
fn failed_source_or_step_budget_preserves_all_histories_and_retry() {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        let (mut runtime, mut thermal) = sources.runtime(cx);
        let mut replay = runtime.clone();
        let mut replay_thermal = thermal.clone();
        let states = runtime.states().to_vec();
        let cells = thermal.cells().clone();
        let source = runtime.source().clone();
        let mut small = controls();
        small.maximum_temperature_increment_k = 1e-6;
        assert!(runtime.step_with_thermal_network(cx, &mut thermal, small,
            &[0.002, 0.006], &[0.0; 2], heat).is_err());
        struct SourceRefusal;
        let refused = runtime.step_with_thermal_network(cx, &mut thermal, controls(),
            &[0.0; 2], &[0.0; 2], |_| Err::<StepCoefficients<()>, _>(SourceRefusal));
        assert!(matches!(refused, Err(ConcentricStringError::Source(_))));
        assert_eq!(runtime.states(), states);
        assert_eq!(runtime.source(), &source);
        assert_eq!(runtime.accepted_samples(), 0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), &cells);
        let actual = runtime.step_with_thermal_network(cx, &mut thermal, controls(),
            &[0.002, 0.006], &[0.0; 2], heat).unwrap();
        let expected = replay.step_with_thermal_network(cx, &mut replay_thermal, controls(),
            &[0.002, 0.006], &[0.0; 2], heat).unwrap();
        assert_eq!(actual.vibration, expected.vibration);
        assert_eq!(actual.layers, expected.layers);
        assert_eq!(thermal.cells(), replay_thermal.cells());
    });
}

#[test]
fn source_resistance_drives_layer_current_heat_and_sound_in_the_same_endpoint_solve() {
    with_cx(|cx, _| {
        let run = |outer_scale| {
            let sources = Sources::with_resistivity(4e-5, outer_scale);
            let (mut runtime, mut thermal) = sources.runtime(cx);
            let mut supplied = 0.0;
            let mut pressure = Vec::new();
            for _ in 0..16 {
                let old_u = independent_internal(&runtime, 4e-5);
                let old_h = thermal.energy_j();
                let frame = runtime.step_with_ohmic_heating(cx, &mut thermal, controls(),
                    &[0.003, 0.002], &[0.0; 2], OhmicDrive::Voltage(2.0), |_| {
                        Ok::<_, Infallible>(StepCoefficients { sources_w: BTreeMap::new(),
                            conductances_w_k: BTreeMap::new(), payload: () })
                    }).unwrap();
                let electric = &frame.source.electrical;
                near(electric.voltage_drop_v, 2.0, 4e-13);
                near(electric.joule_power_w, 2.0*electric.current_a, 4e-13);
                near(frame.source.joule_heat_j,
                    electric.joule_power_w/f64::from(RATE), 4e-13);
                assert_eq!(frame.source.other_external_heat_j, 0.0);
                near(frame.external_heat_j, frame.source.joule_heat_j, 4e-13);
                for i in 0..2 {
                    let inside = if i == 0 { 0.0_f64 } else { RADII[i-1] };
                    let area = core::f64::consts::PI*(RADII[i].powi(2)-inside.powi(2));
                    let t = thermal.temperature_k(CellId(i as u64)).unwrap();
                    let rho_e = (i+1) as f64*1e-6*(if i == 1 { outer_scale } else { 1.0 })
                        *(1.0+0.002*(t-300.0));
                    near(electric.layers[i].current_a, 2.0*area/(rho_e*LENGTH), 4e-13);
                    near(frame.source.layer_joule_heat_j[i],
                        2.0*electric.layers[i].current_a/f64::from(RATE), 4e-13);
                }
                let loss = &frame.vibration.dissipation;
                let balance = frame.source.joule_heat_j
                    -(independent_internal(&runtime, 4e-5)-old_u+thermal.energy_j()-old_h
                        +loss.air_loss_j+loss.authored_loss_j+loss.roundoff_residual_j);
                assert!(balance.abs() <= controls().total_energy_tolerance_j+frame.roundoff_tolerance_j);
                supplied += frame.source.joule_heat_j;
                pressure.push(frame.vibration.acoustic.observer_pressure_pa);
            }
            (supplied, thermal.temperature_k(CellId(1)).unwrap(), runtime.string().tension_n, pressure)
        };
        let a = run(1.0);
        let b = run(2.0);
        assert!(b.0 < a.0 && b.1 < a.1 && b.2 > a.2);
        assert_ne!(a.3, b.3);
    });
}

mod finite_gas;
