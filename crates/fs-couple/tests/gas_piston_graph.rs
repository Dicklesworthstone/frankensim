//! The hand-written gas piston against a declared graph.
//!
//! `fs_couple::cavity_phs::piston` advances a gas and its wall with its own
//! impulse solve. Here the same device is declared from parts: the fs-gas
//! cavity, a mass, and constant sources, joined at the wall's velocity, each
//! through its own face area. Nothing is written for the pair.
//!
//! The module is the reference and stays. Both sides run the gas, controls
//! and clock of the module's own tests (`src/cavity_phs/piston/tests.rs`,
//! copied here: a two-species mixture with a temperature-dependent heat
//! capacity) and must agree at every step. The module's finite thermal body
//! (`piston::thermal::discrete_gradient`) is the phase chart's own body on
//! the contract, the kernel's conductance and the same cavity. The graph has
//! no equivalent of the module's time-resolution admission, of a conductance
//! that depends on state, of its midpoint thermal step, or of its real-gas,
//! reactive and nozzle variants.

use std::convert::Infallible;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
};
use fs_couple::cavity_phs::piston::thermal::discrete_gradient::DiscreteGradientThermalCandidate;
use fs_couple::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, ThermalPistonError,
};
use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Vec3;
use fs_gas::component::GasCavity;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::phase::component::UniformEnthalpyBody;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, SolidLiquidPhase,
    UniformEnthalpyStepInput,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_mbd::component::RigidBody;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion,
    Vec3 as BodyVec3,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_ports::PortKind;
use fs_ports::component::{Class, Component};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_system::elements::{
    Conductance, Force, HeatCapacity, Mass, Motion, SourceAccount, Waveform,
};
use fs_system::{NodeId, SolveConfig, System, SystemBuilder};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};

/// The module tests' clock [s].
const DT_S: f64 = 0.00005;

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn provenance(source: &str) -> Provenance {
    Provenance {
        source: source.into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}

/// One authored species with `Cp/R = cp_r + slope T`.
fn species(name: &str, molar_mass: f64, cp_r: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[3] = slope;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        sources: vec![],
        provenance: provenance("synthetic caloric gas"),
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(molar_mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100_000.0),
            ElementalReferenceIdV1::new("synthetic").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn transport(name: &str) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2.0e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (2000.0, 2.0 * value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        2000.0,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1.0e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance("synthetic transport"),
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "synthetic".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

/// A litre at 400 K and 1 bar: mole fraction `x` of a diatomic-like species
/// whose heat capacity rises with temperature by `slope`, the rest monatomic.
fn gas(x: f64, slope: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species("A", 0.028, 3.5, slope),
                species("B", 0.004, 2.5, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        vec![transport("A"), transport("B")],
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            400.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(model, point, 0.001).unwrap()
}

fn controls() -> PistonSolve {
    PistonSolve {
        volume_bracket_m3: [0.0007, 0.0015],
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1.0e-7,
            maximum_iterations: 128,
        },
        momentum_tolerance_n_s: 2.0e-7 * DT_S,
        energy_tolerance_j: 1.0e-7,
        maximum_iterations: 128,
    }
}

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 29,
            kernel_id: 5,
            tile: 0,
            iteration: 0,
        };
        f(&Cx::new(
            &gate,
            arena,
            key,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

/// One face of the wall: its gas, its area with the sign of the side it is
/// on, and the heat rate [W] it receives as a function of the step's start.
/// A face that nothing heats is the cavity behind an insulating wall, the
/// one-port cavity whose class is a scalar monotone relation.
struct Face {
    gas: ClosedGasParcel,
    signed_area_m2: f64,
    heat_w: Option<fn(f64) -> f64>,
}

struct Wall {
    system: System,
    mass: NodeId,
    gases: Vec<NodeId>,
    exterior: Option<NodeId>,
    force: Option<NodeId>,
}

/// A mass whose velocity is shared by every face's gas through that face's
/// area, by the exterior pressure through minus the first face's area, and
/// by a constant force.
fn wall(
    mass_kg: f64,
    velocity_m_s: f64,
    faces: Vec<Face>,
    exterior_pressure_pa: f64,
    external_force_n: f64,
) -> Wall {
    let mut builder = SystemBuilder::new();
    let first_area = faces[0].signed_area_m2;
    let mass = builder
        .add_with_tangent(
            "wall",
            Mass::new(mass_kg).unwrap(),
            Motion {
                position_m: faces[0].gas.volume_m3() / first_area,
                velocity_m_s,
            },
        )
        .unwrap();
    let mut members = Vec::new();
    let mut gases = Vec::new();
    for (name, face) in ["primary", "opposing"].into_iter().zip(faces) {
        let closed = GasCavity::new(face.gas.clone(), controls().gas).unwrap();
        let cavity = if let Some(rate) = face.heat_w {
            let cavity = builder.add_with_tangent(name, closed, face.gas).unwrap();
            let heat = builder.port(cavity, "heat").unwrap();
            builder
                .boundary(heat, Box::new(move |t0, _, values| values[0] = rate(t0)))
                .unwrap();
            cavity
        } else {
            let insulated = closed.insulated().unwrap();
            assert_eq!(insulated.class(), Class::Monotone1d);
            builder.add_with_tangent(name, insulated, face.gas).unwrap()
        };
        members.push((builder.port(cavity, "p").unwrap(), face.signed_area_m2));
        gases.push(cavity);
    }
    let exterior = (exterior_pressure_pa != 0.0).then(|| {
        let node = builder
            .add_with_tangent(
                "exterior",
                Force::in_domain(
                    PortKind::FluidPressureFlux,
                    Waveform::constant(exterior_pressure_pa).unwrap(),
                ),
                SourceAccount::default(),
            )
            .unwrap();
        members.push((builder.port(node, "p").unwrap(), -first_area));
        node
    });
    let force = (external_force_n != 0.0).then(|| {
        let node = builder
            .add_with_tangent(
                "force",
                Force::new(Waveform::constant(external_force_n).unwrap()),
                SourceAccount::default(),
            )
            .unwrap();
        members.push((builder.port(node, "p").unwrap(), 1.0));
        node
    });
    let velocity = builder.port(mass, "p").unwrap();
    builder.common_flow_ratios(velocity, &members).unwrap();
    let mut system = builder.seal().unwrap();
    system.set_config(SolveConfig {
        relative_tolerance: 1.0e-13,
        ..SolveConfig::default()
    });
    Wall {
        system,
        mass,
        gases,
        exterior,
        force,
    }
}

impl Wall {
    fn velocity_m_s(&self) -> f64 {
        self.system.state::<Motion>(self.mass).unwrap().velocity_m_s
    }
    fn position_m(&self) -> f64 {
        self.system.state::<Motion>(self.mass).unwrap().position_m
    }
    fn gas(&self, face: usize) -> &ClosedGasParcel {
        self.system
            .state::<ClosedGasParcel>(self.gases[face])
            .unwrap()
    }
    /// Work the source has received [J]; zero when the graph has no such
    /// source because its value is zero.
    fn received_j(&self, node: Option<NodeId>) -> f64 {
        node.map_or(0.0, |node| {
            self.system
                .state::<SourceAccount>(node)
                .unwrap()
                .energy_in_j
        })
    }
}

fn temperature_k(gas: &ClosedGasParcel) -> f64 {
    gas.state().state().temperature
}

/// Largest differences seen over a run, each in its own unit.
#[derive(Debug, Default)]
struct Worst {
    velocity_m_s: f64,
    volume_m3: f64,
    temperature_k: f64,
}

impl Worst {
    fn see(&mut self, wall: &Wall, face: usize, velocity_m_s: f64, gas: &ClosedGasParcel) {
        self.velocity_m_s = self
            .velocity_m_s
            .max((wall.velocity_m_s() - velocity_m_s).abs());
        self.volume_m3 = self
            .volume_m3
            .max((wall.gas(face).volume_m3() - gas.volume_m3()).abs());
        self.temperature_k = self
            .temperature_k
            .max((temperature_k(wall.gas(face)) - temperature_k(gas)).abs());
    }

    /// The module stops its impulse solve at `momentum_tolerance_n_s`, which
    /// on the unit mass is 1e-11 m/s a step; the graph solves to 1e-13 of its
    /// unknowns. The largest differences measured over these runs are
    /// 1e-11 m/s, 2e-16 m3 and 4e-11 K.
    fn check(&self) {
        assert!(self.velocity_m_s <= 1.0e-10, "{self:?}");
        assert!(self.volume_m3 <= 1.0e-15, "{self:?}");
        assert!(self.temperature_k <= 1.0e-9, "{self:?}");
    }
}

/// Heat for the first eight steps of the module's clock: 20 J a step.
fn early_heat_w(t0_s: f64) -> f64 {
    if t0_s < 7.5 * DT_S { 400_000.0 } else { 0.0 }
}

#[test]
fn a_wall_on_one_gas_with_exterior_pressure_force_and_heat_matches_the_module() {
    let dt = DT_S;
    let heated: Option<fn(f64) -> f64> = Some(early_heat_w);
    for (velocity, exterior_pa, force_n, heat_w) in [
        // Released into vacuum from rest, behind an insulating wall.
        (0.0, 0.0, 0.0, None),
        // Pushed inwards against the atmosphere while the gas is heated.
        (-0.2, 100_000.0, -0.5, heated),
        // Thrown outwards with nothing but the atmosphere on it, insulated.
        (0.2, 100_000.0, 0.0, None),
    ] {
        let mut module = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, velocity).unwrap();
        let initial = module.gas().clone();
        let mut graph = wall(
            1.0,
            velocity,
            vec![Face {
                gas: gas(0.75, 0.001),
                signed_area_m2: 0.01,
                heat_w,
            }],
            exterior_pa,
            force_n,
        );
        let mut worst = Worst::default();
        let (mut exterior_work_j, mut force_work_j) = (0.0, 0.0);
        with_cx(|cx| {
            for step in 0..48 {
                let input = PistonInput {
                    duration_s: dt,
                    heat_j: heat_w.map_or(0.0, |rate| dt * rate(f64::from(step) * dt)),
                    exterior_pressure_pa: exterior_pa,
                    external_force_n: force_n,
                };
                let report = module
                    .step_discrete_gradient(input, controls(), || false)
                    .unwrap();
                exterior_work_j += report.exterior_pressure_work_j;
                force_work_j += report.external_force_work_j;
                graph.system.step(dt, cx).unwrap();
                worst.see(&graph, 0, module.velocity_m_s(), module.gas());
            }
        });
        // The sources' accounts are the module's work terms: work done on
        // the exterior is received by it, work done by the force leaves it.
        let exterior_gap_j = graph.received_j(graph.exterior) - exterior_work_j;
        let force_gap_j = graph.received_j(graph.force) + force_work_j;
        let ledger_j = graph.system.ledger().energy_residual_j();
        eprintln!(
            "one gas, v0 {velocity}: {worst:?}; exterior work {exterior_work_j:.6} J \
             (gap {exterior_gap_j:.2e}), force work {force_work_j:.6} J (gap {force_gap_j:.2e}), \
             ledger residual {ledger_j:.2e} J"
        );
        // The run did something: the wall moved and the gas changed.
        assert!((module.gas().volume_m3() - initial.volume_m3()).abs() > 1.0e-8);
        assert!((temperature_k(module.gas()) - temperature_k(&initial)).abs() > 1.0e-3);
        worst.check();
        // The graph keeps a second position; it stays the gas's volume.
        assert!((0.01 * graph.position_m() - graph.gas(0).volume_m3()).abs() <= 1.0e-15);
        // Measured: 2e-11 J of about 1 J of exterior work, 8e-15 J of force
        // work, and a ledger residual of 2e-14 J.
        assert!(exterior_gap_j.abs() <= 2.0e-10, "{exterior_gap_j}");
        assert!(force_gap_j.abs() <= 1.0e-12, "{force_gap_j}");
        assert!(ledger_j.abs() <= 1.0e-12, "{ledger_j}");
    }
}

#[test]
fn a_wall_between_two_gases_matches_the_double_acting_module() {
    let dt = DT_S;
    // A rod on the opposing side: its face is smaller than the primary one.
    for (velocity, opposing_area, force_n) in [(0.3, 0.01, 0.0), (-0.2, 0.008, 1.5)] {
        let primary = GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, velocity).unwrap();
        let mut module =
            DoubleActingGasPiston::new(primary, gas(0.25, 0.001), opposing_area).unwrap();
        let initial = module.clone();
        let mut graph = wall(
            1.0,
            velocity,
            vec![
                Face {
                    gas: gas(0.75, 0.001),
                    signed_area_m2: 0.01,
                    heat_w: Some(early_heat_w),
                },
                // Nothing heats the opposing gas: the insulated cavity.
                Face {
                    gas: gas(0.25, 0.001),
                    signed_area_m2: -opposing_area,
                    heat_w: None,
                },
            ],
            0.0,
            force_n,
        );
        // The module tries both ends of its volume bracket. At the far end
        // of the one-gas bracket the opposing gas would be halved in one
        // step, which its owner refuses for this mixture, so the bracket
        // here is one both gases admit. The graph has no bracket.
        let solve = DoublePistonSolve {
            piston: PistonSolve {
                volume_bracket_m3: [0.0009, 0.0011],
                ..controls()
            },
            opposing_gas: controls().gas,
            volume_tolerance_m3: 1.0e-18,
        };
        let mut worst = Worst::default();
        let mut force_work_j = 0.0;
        with_cx(|cx| {
            for step in 0..48 {
                let input = DoublePistonInput {
                    duration_s: dt,
                    heat_j: [dt * early_heat_w(f64::from(step) * dt), 0.0],
                    external_force_n: force_n,
                };
                let report = module
                    .step_discrete_gradient(input, solve, || false)
                    .unwrap_or_else(|error| panic!("module refused step {step}: {error:?}"));
                force_work_j += report.external_force_work_j;
                graph.system.step(dt, cx).unwrap();
                let velocity = module.primary().velocity_m_s();
                worst.see(&graph, 0, velocity, module.primary().gas());
                worst.see(&graph, 1, velocity, module.opposing());
            }
        });
        let force_gap_j = graph.received_j(graph.force) + force_work_j;
        let ledger_j = graph.system.ledger().energy_residual_j();
        eprintln!(
            "two gases, v0 {velocity}: {worst:?}; force work {force_work_j:.6} J \
             (gap {force_gap_j:.2e}), ledger residual {ledger_j:.2e} J"
        );
        for (now, before) in [
            (module.primary().gas(), initial.primary().gas()),
            (module.opposing(), initial.opposing()),
        ] {
            assert!((now.volume_m3() - before.volume_m3()).abs() > 1.0e-8);
            assert!((temperature_k(now) - temperature_k(before)).abs() > 1.0e-3);
        }
        worst.check();
        assert!(force_gap_j.abs() <= 1.0e-12, "{force_gap_j}");
        assert!(ledger_j.abs() <= 1.0e-12, "{ledger_j}");
    }
}

/// A 0.1 m square plate, 10 mm thick, of an authored solid with the given
/// density: the fixture of the module's material-wall tests. Its mass is
/// what the rigid-body compiler makes of the density claim and the geometry.
fn plate(density_kg_m3: f64) -> ResolvedBody {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", Density::DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, density_kg_m3), (400.0, density_kg_m3)],
                dims: Density::DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance("synthetic specimen density"),
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic specimen".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    resolve_body(
        &[RegionBinding {
            key: "wall",
            material: &card,
            point: &point,
            geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.01)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}

/// A free fs-mbd rigid body at `z = V / area` moving along `z`, with the
/// gas on its face through a 1 x 3 power map of the area, no torque, and
/// the early heat on the gas.
fn plate_on_gas(
    properties: MassProperties,
    velocity_m_s: f64,
    area_m2: f64,
    gas: &ClosedGasParcel,
) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let start = RigidBodyState::new(
        Pose::new(
            BodyVec3::new(0.0, 0.0, gas.volume_m3() / area_m2),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        BodyVec3::new(0.0, 0.0, properties.mass() * velocity_m_s),
        BodyVec3::ZERO,
    )
    .unwrap();
    let wall = builder
        .add(
            "wall",
            RigidBody::new(RigidBodyIntegrator::new(Gravity::ZERO), properties).unwrap(),
            start,
        )
        .unwrap();
    let cavity = builder
        .add_with_tangent(
            "gas",
            GasCavity::new(gas.clone(), controls().gas).unwrap(),
            gas.clone(),
        )
        .unwrap();
    builder
        .power_map(
            builder.port(wall, "force").unwrap(),
            builder.port(cavity, "p").unwrap(),
            vec![0.0, 0.0, area_m2],
        )
        .unwrap();
    let torque = builder.port(wall, "torque").unwrap();
    builder
        .boundary(torque, Box::new(|_, _, values| values.fill(0.0)))
        .unwrap();
    let heat = builder.port(cavity, "heat").unwrap();
    builder
        .boundary(
            heat,
            Box::new(|t0, _, values| values[0] = early_heat_w(t0)),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    system.set_config(SolveConfig {
        relative_tolerance: 1.0e-13,
        ..SolveConfig::default()
    });
    (system, wall, cavity)
}

#[test]
fn a_material_plate_as_a_free_rigid_body_on_the_gas_matches_the_module() {
    // The module's material wall takes its mass from the resolved body and
    // otherwise is `GasPiston`. Here the same body is the fs-mbd rigid body
    // in three dimensions, and the gas pushes it along the plate's normal
    // through a 1 x 3 power map of the face area. Nothing constrains the
    // body: it stays on the axis because nothing pushes it off.
    let area = 0.01;
    for (density, velocity) in [(1000.0, 0.0), (8000.0, -0.2)] {
        let body = plate(density);
        let inertia = body.inertia();
        let tensor = inertia.inertia_com().m;
        // The plate's own frame is principal and centred, so its tensor
        // goes to the rigid-body owner as three moments.
        let offset = inertia.center_of_mass();
        assert!(offset.x.abs() + offset.y.abs() + offset.z.abs() <= 0.0);
        assert!(
            [1, 2, 3, 5, 6, 7]
                .iter()
                .all(|&entry| tensor[entry].abs() <= 0.0),
            "{tensor:?}"
        );
        let mass_kg = inertia.mass();
        let properties = MassProperties::new(
            mass_kg,
            BodyVec3::ZERO,
            BodyVec3::new(tensor[0], tensor[4], tensor[8]),
        )
        .unwrap();
        let mut module = GasPiston::new(gas(0.75, 0.001), mass_kg, area, velocity).unwrap();
        let initial = module.gas().clone();
        let (mut system, wall, cavity) = plate_on_gas(properties, velocity, area, &initial);

        let (mut worst_velocity, mut worst_volume, mut worst_temperature) = (0.0_f64, 0.0_f64, 0.0_f64);
        let mut off_axis = 0.0_f64;
        with_cx(|cx| {
            for step in 0..48 {
                let input = PistonInput {
                    duration_s: DT_S,
                    heat_j: DT_S * early_heat_w(f64::from(step) * DT_S),
                    exterior_pressure_pa: 0.0,
                    external_force_n: 0.0,
                };
                module
                    .step_discrete_gradient(input, controls(), || false)
                    .unwrap_or_else(|error| panic!("module refused step {step}: {error:?}"));
                system.step(DT_S, cx).unwrap();
                let state = system.state::<RigidBodyState>(wall).unwrap();
                let momentum = state.linear_momentum_world();
                let position = state.pose().position_world();
                let spin = state.angular_momentum_body();
                let parcel = system.state::<ClosedGasParcel>(cavity).unwrap();
                worst_velocity =
                    worst_velocity.max((momentum.z / mass_kg - module.velocity_m_s()).abs());
                worst_volume =
                    worst_volume.max((parcel.volume_m3() - module.gas().volume_m3()).abs());
                worst_temperature = worst_temperature
                    .max((temperature_k(parcel) - temperature_k(module.gas())).abs());
                off_axis = off_axis
                    .max(momentum.x.abs())
                    .max(momentum.y.abs())
                    .max(position.x.abs())
                    .max(position.y.abs())
                    .max(spin.dot(spin).sqrt())
                    .max((area * position.z - parcel.volume_m3()).abs());
            }
        });
        let ledger_j = system.ledger().energy_residual_j();
        eprintln!(
            "plate of {mass_kg} kg, v0 {velocity}: velocity {worst_velocity:.2e} m/s, volume \
             {worst_volume:.2e} m3, temperature {worst_temperature:.2e} K, off axis \
             {off_axis:.2e}, ledger residual {ledger_j:.2e} J"
        );
        assert!((module.gas().volume_m3() - initial.volume_m3()).abs() > 1.0e-8);
        // Measured: 9e-12 m/s, 2e-16 m3, 4e-11 K, 9e-19 off the axis and a
        // ledger residual of 4e-14 J.
        assert!(worst_velocity <= 1.0e-10, "{worst_velocity}");
        assert!(worst_volume <= 1.0e-15, "{worst_volume}");
        assert!(worst_temperature <= 1.0e-9, "{worst_temperature}");
        assert!(off_axis <= 1.0e-16, "{off_axis}");
        assert!(ledger_j.abs() <= 1.0e-12, "{ledger_j}");
    }
}

/// Mass [kg] of the finite body of the module's thermal tests
/// (`src/cavity_phs/piston/thermal/discrete_gradient/tests.rs`), whose
/// fixtures follow, copied.
const BODY_KG: f64 = 0.001;

/// A solid of 500 J/(kg K) from 200 K to 1000 K.
fn sensible() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        fs_blake3::hash_bytes(b"synthetic entropy contact cp=500"),
        SolidLiquidPhase::Solid,
        [200.0, 1000.0]
            .into_iter()
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 500.0 * (temperature_k - 200.0),
                temperature_k,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 4000.0,
            })
            .collect(),
    )
    .unwrap()
}

/// A solid from 200 K that melts at `transition` between 100 and 160 kJ/kg
/// and is liquid to 900 K at 400 kJ/kg.
fn latent(transition: f64) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, t, f| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: 4000.0,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        fs_blake3::hash_bytes(format!("synthetic latent contact T={transition}").as_bytes()),
        vec![
            knot(0.0, 200.0, 0.0),
            knot(100_000.0, transition, 0.0),
            knot(160_000.0, transition, 1.0),
            knot(400_000.0, 900.0, 1.0),
        ],
    )
    .unwrap()
}

/// What stores the heat on the body's side of the contact.
enum Store<'a> {
    /// The chart's own body: the phase owner on the composition contract.
    /// Its heat port returns the chart's entropy-conjugate mean, which is
    /// the temperature the module's contact law uses.
    Chart(&'a EquilibriumEnthalpyPhaseCurve),
    /// The kernel's ideal heat capacity [J/K] in its place. It returns the
    /// arithmetic mean of its two end temperatures: another discrete law,
    /// put here to show that the comparison tells the two apart.
    Arithmetic(f64),
}

struct Contact {
    system: System,
    gas: NodeId,
    body: NodeId,
    wall: Option<NodeId>,
}

/// A gas in contact with a finite body through a constant conductance
/// [W/K]: `body.heat`, the kernel's `Conductance`, `gas.heat`. With a
/// velocity the gas is behind a 1 kg wall of 0.01 m2 under the atmosphere,
/// as in [`wall`]; without one its walls are fixed.
fn contact(
    gas: &ClosedGasParcel,
    store: &Store<'_>,
    start: EquilibriumPhaseState,
    conductance_w_k: f64,
    velocity_m_s: Option<f64>,
) -> Contact {
    let mut builder = SystemBuilder::new();
    let cavity = builder
        .add_with_tangent(
            "gas",
            GasCavity::new(gas.clone(), controls().gas).unwrap(),
            gas.clone(),
        )
        .unwrap();
    let link = builder
        .add_with_tangent("contact", Conductance::new(conductance_w_k).unwrap(), ())
        .unwrap();
    let body = match store {
        Store::Chart(curve) => {
            let body = UniformEnthalpyBody::new((*curve).clone(), BODY_KG).unwrap();
            assert_eq!(body.class(), Class::Monotone1d);
            builder.add_with_tangent("body", body, start).unwrap()
        }
        Store::Arithmetic(capacity_j_k) => builder
            .add_with_tangent(
                "body",
                HeatCapacity::new(*capacity_j_k).unwrap(),
                start.temperature_k(),
            )
            .unwrap(),
    };
    let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
    let (a, b) = (port(&builder, link, "a"), port(&builder, link, "b"));
    builder
        .common_effort(port(&builder, body, "heat"), &[a])
        .unwrap();
    builder
        .common_effort(port(&builder, cavity, "heat"), &[b])
        .unwrap();
    let pressure = port(&builder, cavity, "p");
    let wall = if let Some(velocity_m_s) = velocity_m_s {
        let mass = builder
            .add_with_tangent(
                "wall",
                Mass::new(1.0).unwrap(),
                Motion {
                    position_m: gas.volume_m3() / 0.01,
                    velocity_m_s,
                },
            )
            .unwrap();
        let exterior = builder
            .add_with_tangent(
                "exterior",
                Force::in_domain(
                    PortKind::FluidPressureFlux,
                    Waveform::constant(100_000.0).unwrap(),
                ),
                SourceAccount::default(),
            )
            .unwrap();
        builder
            .common_flow_ratios(
                port(&builder, mass, "p"),
                &[(pressure, 0.01), (port(&builder, exterior, "p"), -0.01)],
            )
            .unwrap();
        Some(mass)
    } else {
        builder
            .boundary(pressure, Box::new(|_, _, values| values[0] = 0.0))
            .unwrap();
        None
    };
    let mut system = builder.seal().unwrap();
    system.set_config(SolveConfig {
        relative_tolerance: 1.0e-13,
        ..SolveConfig::default()
    });
    // The conductance's two temperature inputs are in the loop and need a
    // start.
    assert!(system.guess(a, &[start.temperature_k()]));
    assert!(system.guess(b, &[temperature_k(gas)]));
    Contact {
        system,
        gas: cavity,
        body,
        wall,
    }
}

/// `|a - b|`, which must be a number: a NaN would pass a `max`.
fn gap(a: f64, b: f64) -> f64 {
    let difference = (a - b).abs();
    assert!(difference.is_finite(), "{a} against {b}");
    difference
}

/// Largest differences between the graph and the module, each in its own
/// unit: after every step the two bodies' temperature [K] and specific
/// enthalpy [J/kg], the two gases' temperature [K] and volume [m3] and the
/// walls' velocity [m/s]; within every step the heat exchanged [J] and the
/// two temperatures the contact law was given [K]; over the run the entropy
/// the contact produced [J/K].
#[derive(Debug, Default)]
struct Apart {
    body_k: f64,
    body_j_kg: f64,
    gas_k: f64,
    volume_m3: f64,
    velocity_m_s: f64,
    heat_j: f64,
    mean_k: f64,
    entropy_j_k: f64,
}

impl Apart {
    /// What the two may differ by after `steps` steps of a contact of
    /// `conductance_w_k` between `gas` and a body whose smallest heat
    /// capacity is `capacity_j_kg_k`, with the module's heat bracket at
    /// `bracket_j`. All of it is from the two solves' own tolerances and
    /// none of it from a run.
    ///
    /// - The module accepts a heat `q` with `|q - dt G (theta_b - theta_g)|`
    ///   up to its source tolerance; the graph's unknowns converge to 1e-13
    ///   of their size, a tenth of that here, and are allowed as much again.
    /// - Each side solves the gas's energy equation to its tolerance per
    ///   kilogram.
    /// - The exact step does not amplify a difference between two runs:
    ///   contact damps the difference of the two stores' energies and
    ///   keeps their sum. A heat error moves energy between the stores and
    ///   is damped by later steps; a gas error changes their sum and stays,
    ///   and on a plateau all of it ends in the body. So after `steps`
    ///   steps either store's energy differs by at most that many errors of
    ///   both kinds, and each temperature by its energy over its heat
    ///   capacity (the gas's taken at the start, with a tenth more for its
    ///   change with temperature).
    /// - A step's heat follows the two means, which differ at most as the
    ///   states do; its entropy is `q (1 / theta_g - 1 / theta_b)`, no
    ///   temperature being below the controls' 200 K.
    /// - The module stops its impulse solve at its momentum tolerance, on a
    ///   wall of 1 kg and 0.01 m2. With the wall `moving`, the two first
    ///   order couplings are added, each from the other's bound without
    ///   it: a volume error changes the gas temperature by at most
    ///   `(R / cv) T dV / V`, and a gas temperature error changes the
    ///   pressure by `(m R / V) dT`, which pushes the wall by that times
    ///   `A dt / m_wall` a step; `T` below the controls' maximum and `V`
    ///   above the lower end of their volume bracket.
    fn allowed(
        gas: &ClosedGasParcel,
        capacity_j_kg_k: f64,
        conductance_w_k: f64,
        bracket_j: f64,
        moving: bool,
        steps: u32,
    ) -> Self {
        let n = f64::from(steps);
        let heat_error_j = 2.0 * exchange(bracket_j).source_tolerance_j;
        let gas_error_j = 2.0 * gas.mass_kg() * controls().gas.energy_tolerance_j_kg;
        let energy_j = n * (heat_error_j + gas_error_j);
        let cv = gas
            .state()
            .thermodynamics()
            .properties()
            .mass_specific_cv()
            .value();
        let gas_constant = gas.state().state().specific_gas_constant;
        let (wall_kg, area_m2) = (1.0, 0.01);
        let least_m3 = controls().volume_bracket_m3[0];
        let body_j_kg = energy_j / BODY_KG;
        let body_k = body_j_kg / capacity_j_kg_k;
        let mut gas_k = 1.1 * energy_j / (gas.mass_kg() * cv);
        let mut velocity_m_s = n * controls().momentum_tolerance_n_s / wall_kg;
        if moving {
            let volume_m3 = area_m2 * DT_S * n * velocity_m_s;
            let pressure_pa = gas.mass_kg() * gas_constant / least_m3 * gas_k;
            velocity_m_s += n * pressure_pa * area_m2 * DT_S / wall_kg;
            gas_k += gas_constant / cv * controls().gas.maximum_temperature_k * volume_m3 / least_m3;
        }
        let mean_k = body_k + gas_k;
        let heat_j = heat_error_j + DT_S * conductance_w_k * mean_k;
        let floor_k = controls().gas.minimum_temperature_k;
        Self {
            body_k,
            body_j_kg,
            gas_k,
            volume_m3: area_m2 * DT_S * n * velocity_m_s,
            velocity_m_s,
            heat_j,
            mean_k,
            entropy_j_k: n * (heat_j / floor_k + bracket_j * mean_k / (floor_k * floor_k)),
        }
    }

    /// After a step: the graph's states against the module's gas, body and
    /// wall velocity.
    fn see_states(
        &mut self,
        graph: &Contact,
        store: &Store<'_>,
        gas: &ClosedGasParcel,
        body: EquilibriumPhaseState,
        velocity_m_s: f64,
    ) {
        let held = graph.system.state::<ClosedGasParcel>(graph.gas).unwrap();
        self.gas_k = self.gas_k.max(gap(temperature_k(held), temperature_k(gas)));
        self.volume_m3 = self.volume_m3.max(gap(held.volume_m3(), gas.volume_m3()));
        match store {
            Store::Chart(_) => {
                let state = graph
                    .system
                    .state::<EquilibriumPhaseState>(graph.body)
                    .unwrap();
                self.body_k = self
                    .body_k
                    .max(gap(state.temperature_k(), body.temperature_k()));
                self.body_j_kg = self.body_j_kg.max(gap(
                    state.specific_enthalpy_j_kg(),
                    body.specific_enthalpy_j_kg(),
                ));
            }
            Store::Arithmetic(_) => {
                let kelvin = *graph.system.state::<f64>(graph.body).unwrap();
                self.body_k = self.body_k.max(gap(kelvin, body.temperature_k()));
            }
        }
        if let Some(wall) = graph.wall {
            let velocity = graph.system.state::<Motion>(wall).unwrap().velocity_m_s;
            self.velocity_m_s = self.velocity_m_s.max(gap(velocity, velocity_m_s));
        }
    }

    /// Within the step just taken: the heat that went into the graph's gas
    /// and the two temperatures its contact law was given, against the heat
    /// the module exchanged [J] and its two means [K].
    fn see_step(&mut self, graph: &Contact, heat_j: f64, body_mean_k: f64, gas_mean_k: f64) {
        let body_heat = graph.system.port(graph.body, "heat").unwrap();
        let gas_heat = graph.system.port(graph.gas, "heat").unwrap();
        let into_gas_j = DT_S * graph.system.last_input(gas_heat)[0];
        self.heat_j = self.heat_j.max(gap(into_gas_j, heat_j));
        self.mean_k = self
            .mean_k
            .max(gap(graph.system.last_output(body_heat)[0], body_mean_k))
            .max(gap(graph.system.last_output(gas_heat)[0], gas_mean_k));
    }

    fn within(&self, allowed: &Self) -> bool {
        self.body_k <= allowed.body_k
            && self.body_j_kg <= allowed.body_j_kg
            && self.gas_k <= allowed.gas_k
            && self.volume_m3 <= allowed.volume_m3
            && self.velocity_m_s <= allowed.velocity_m_s
            && self.heat_j <= allowed.heat_j
            && self.mean_k <= allowed.mean_k
            && self.entropy_j_k <= allowed.entropy_j_k
    }
}

/// The module's heat bracket in its own thermal fixture [J], either way.
const FIXTURE_BRACKET_J: f64 = 100.0;

/// The module's exchange controls with the caller's heat bracket; the
/// tolerances are the fixture's. The module's solve tries the far end of
/// its bracket on the side of the root, and on some of these runs the gas
/// or the body cannot give or take the fixture's 100 J inside its source's
/// range a few steps on. A run has a narrower bracket only where
/// [`contact_against_the_module`] has first seen the module refuse its
/// own, and it prints the step and the refusal.
fn exchange(bracket_j: f64) -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-bracket_j, bracket_j],
        source_tolerance_j: 1.0e-10,
        energy_tolerance_j: 1.0e-7,
        maximum_iterations: 128,
    }
}

/// The module's piston for these runs: `gas` behind a 1 kg wall of 0.01 m2
/// thrown at `velocity_m_s`, or, without a velocity, with fixed walls.
fn module_piston(gas: ClosedGasParcel, velocity_m_s: Option<f64>) -> GasPiston {
    let free = GasPiston::new(gas, 1.0, 0.01, velocity_m_s.unwrap_or(0.0)).unwrap();
    if velocity_m_s.is_some() {
        free
    } else {
        free.with_fixed_support().unwrap()
    }
}

/// The module's body fixture on `curve`, at `start`.
fn module_body(
    curve: &EquilibriumEnthalpyPhaseCurve,
    start: EquilibriumPhaseState,
) -> UniformEnthalpyStepInput<'_> {
    UniformEnthalpyStepInput {
        curve,
        initial: start,
        mass_kg: BODY_KG,
        volume_m3: BODY_KG / 4000.0,
        surface_area_m2: 0.02,
        internal_heat_j: 0.0,
        duration_s: DT_S,
    }
}

/// One step of the module's finite-body contact under the atmosphere, with
/// a constant conductance [W/K] and the heat bracket `bracket_j`.
fn module_step(
    module: &GasPiston,
    thermal: UniformEnthalpyStepInput<'_>,
    conductance_w_k: f64,
    bracket_j: f64,
) -> Result<DiscreteGradientThermalCandidate<()>, ThermalPistonError<Infallible>> {
    let input = PistonInput {
        duration_s: DT_S,
        heat_j: 0.0,
        exterior_pressure_pa: 100_000.0,
        external_force_n: 0.0,
    };
    module.thermal_discrete_gradient_candidate(
        input,
        thermal,
        controls(),
        exchange(bracket_j),
        |_, _| {
            Ok::<_, Infallible>(ConductanceSample {
                conductance_w_k,
                source: (),
            })
        },
        || false,
    )
}

/// What [`contact_against_the_module`] found.
struct ContactRun {
    apart: Apart,
    /// The module's body, gas temperature [K] and gas volume [m3] at the
    /// end.
    body: EquilibriumPhaseState,
    gas_k: f64,
    gas_m3: f64,
    /// Entropy the module's contact produced [J/K] and heat it moved, in
    /// either direction [J].
    produced_j_k: f64,
    moved_j: f64,
    /// The graph's energy residual [J] and Newton iterations.
    ledger_j: f64,
    iterations: u32,
}

/// The module's discrete-gradient thermal step and the graph of [`contact`]
/// from the same start, a step at a time: `gas` and a body of `curve` at
/// `enthalpy_j_kg`, a constant conductance, fixed walls or a wall thrown at
/// `velocity_m_s`.
#[allow(clippy::too_many_arguments)]
fn contact_against_the_module(
    gas: ClosedGasParcel,
    curve: &EquilibriumEnthalpyPhaseCurve,
    enthalpy_j_kg: f64,
    store: &Store<'_>,
    conductance_w_k: f64,
    bracket_j: f64,
    velocity_m_s: Option<f64>,
    steps: u32,
) -> ContactRun {
    let start = curve.state_at_specific_enthalpy(enthalpy_j_kg).unwrap();
    if bracket_j < FIXTURE_BRACKET_J {
        // A bracket narrower than the module's fixture is used only where
        // the module, alone, refuses a step of this run with its own.
        let mut module = module_piston(gas.clone(), velocity_m_s);
        let mut thermal = module_body(curve, start);
        let refusal = (0..steps).find_map(|step| {
            match module_step(&module, thermal, conductance_w_k, FIXTURE_BRACKET_J) {
                Ok(candidate) => {
                    let (next, accepted, _) = candidate.into_parts();
                    module = next;
                    thermal.initial = accepted.state;
                    None
                }
                Err(error) => Some((step, format!("{error:?}"))),
            }
        });
        let (step, why) =
            refusal.expect("the module takes its fixture's bracket on this run: use it");
        eprintln!(
            "from {:.0} K with its fixture's {FIXTURE_BRACKET_J} J either way the module refuses \
             step {step}: {why}",
            start.temperature_k()
        );
    }
    let mut graph = contact(&gas, store, start, conductance_w_k, velocity_m_s);
    let mut module = module_piston(gas, velocity_m_s);
    let mut thermal = module_body(curve, start);
    let mut apart = Apart::default();
    let (mut produced_j_k, mut moved_j, mut iterations) = (0.0, 0.0, 0_u32);
    with_cx(|cx| {
        for step in 0..steps {
            let candidate = module_step(&module, thermal, conductance_w_k, bracket_j)
                .unwrap_or_else(|error| panic!("module refused step {step}: {error:?}"));
            let (next, accepted, report) = candidate.into_parts();
            module = next;
            thermal.initial = accepted.state;
            produced_j_k += report.contact_entropy_production_j_k;
            moved_j += report.exchanged_heat_j.abs();

            iterations += graph
                .system
                .step(DT_S, cx)
                .unwrap_or_else(|error| panic!("graph refused step {step}: {error}"))
                .iterations;
            apart.see_states(
                &graph,
                store,
                module.gas(),
                thermal.initial,
                module.velocity_m_s(),
            );
            apart.see_step(
                &graph,
                report.exchanged_heat_j,
                report.thermal_mean_temperature_k,
                report.piston.thermodynamic_mean.state().temperature,
            );
        }
    });
    apart.entropy_j_k = gap(
        graph.system.ledger().entropy_generated_j_per_k(),
        produced_j_k,
    );
    ContactRun {
        apart,
        body: thermal.initial,
        gas_k: temperature_k(module.gas()),
        gas_m3: module.gas().volume_m3(),
        produced_j_k,
        moved_j,
        ledger_j: graph.system.ledger().energy_residual_j(),
        iterations,
    }
}

#[test]
fn a_gas_in_contact_with_a_finite_body_matches_the_modules_thermal_step_in_both_directions() {
    // The module's sensible fixture: fixed walls, the gas at 400 K, the
    // body at 300 K and at 650 K, 4000 W/K. These are no small steps: in
    // the twelve the hot body gives 66 J and cools by 132 K (measured).
    let curve = sensible();
    // From 650 K the module's own 100 J bracket serves; from 300 K it does
    // not (the run prints where the module refuses it).
    let (conductance, steps) = (4000.0, 12);
    for (body_k, bracket) in [(300.0, 60.0), (650.0, FIXTURE_BRACKET_J)] {
        let enthalpy = 500.0 * (body_k - 200.0);
        let allowed = Apart::allowed(&gas(0.75, 0.0), 500.0, conductance, bracket, false, steps);
        let run = contact_against_the_module(
            gas(0.75, 0.0),
            &curve,
            enthalpy,
            &Store::Chart(&curve),
            conductance,
            bracket,
            None,
            steps,
        );
        eprintln!(
            "contact, body from {body_k} K to {:.3} K, gas to {:.3} K, {:.3} J moved, \
             {:.6e} J/K produced, {} iterations, ledger residual {:.2e} J: {:?}; allowed {allowed:?}",
            run.body.temperature_k(),
            run.gas_k,
            run.moved_j,
            run.produced_j_k,
            run.iterations,
            run.ledger_j,
            run.apart,
        );
        // It did something: the body went more than 40 K towards the gas,
        // the gas more than 20 K towards the body, and entropy was made.
        assert!((run.body.temperature_k() - body_k) * (400.0 - body_k) > 0.0);
        assert!((run.body.temperature_k() - body_k).abs() > 40.0);
        assert!((run.gas_k - 400.0).abs() > 20.0);
        assert!(run.produced_j_k > 1.0e-3);
        assert!(run.apart.within(&allowed), "{:?} against {allowed:?}", run.apart);
        assert!(run.ledger_j.abs() <= 1.0e-10, "{}", run.ledger_j);
        // The states above do not depend on the tangents the graph is
        // given, only its iterations do: with the owners' own, two or
        // three a step (measured 27 and 28 in the 12).
        assert!(run.iterations <= 3 * steps, "{}", run.iterations);

        // The same graph with the kernel's ideal heat capacity for the
        // body, 0.5 J/K: the same energy, but the arithmetic mean of the
        // end temperatures where the module's law has the entropy-conjugate
        // one. The two means differ by about dT^2 / (12 T) for a step that
        // changes the body by dT, and the comparison must say so. Measured:
        // 0.19 K and 0.61 K in the temperature the contact is given, 0.065 K
        // and 0.21 K in the body.
        let other = contact_against_the_module(
            gas(0.75, 0.0),
            &curve,
            enthalpy,
            &Store::Arithmetic(BODY_KG * 500.0),
            conductance,
            bracket,
            None,
            steps,
        );
        eprintln!(
            "the arithmetic mean in its place, body from {body_k} K: {:?}",
            other.apart
        );
        assert!(!other.apart.within(&allowed));
        assert!(
            other.apart.body_k > 1.0e6 * allowed.body_k
                && other.apart.heat_j > 1.0e6 * allowed.heat_j
                && other.apart.mean_k > 1.0e6 * allowed.mean_k,
            "{:?} against {allowed:?}",
            other.apart
        );
    }
}

#[test]
fn a_moving_wall_a_gas_and_a_finite_body_in_contact_match_the_modules_thermal_step() {
    // The module's moving fixture: the wall thrown at 2 m/s either way, a
    // gas whose heat capacity varies, the body at 650 K, 2000 W/K.
    let curve = sensible();
    let (conductance, bracket, steps) = (2000.0, FIXTURE_BRACKET_J, 12);
    for velocity in [-2.0, 2.0] {
        let allowed = Apart::allowed(&gas(0.75, 0.001), 500.0, conductance, bracket, true, steps);
        let run = contact_against_the_module(
            gas(0.75, 0.001),
            &curve,
            225_000.0,
            &Store::Chart(&curve),
            conductance,
            bracket,
            Some(velocity),
            steps,
        );
        eprintln!(
            "contact behind a wall at {velocity} m/s, body to {:.3} K, gas to {:.3} K, {:.3} J \
             moved, {:.6e} J/K produced, {} iterations, ledger residual {:.2e} J: {:?}; allowed \
             {allowed:?}",
            run.body.temperature_k(),
            run.gas_k,
            run.moved_j,
            run.produced_j_k,
            run.iterations,
            run.ledger_j,
            run.apart,
        );
        // Heat crossed and the wall moved, the way it was thrown.
        assert!(run.body.temperature_k() < 600.0 && run.gas_k > 430.0);
        assert!((run.gas_m3 - 0.001) * velocity > 2.0e-6, "{}", run.gas_m3);
        assert!(run.apart.within(&allowed), "{:?} against {allowed:?}", run.apart);
        assert!(run.ledger_j.abs() <= 1.0e-10, "{}", run.ledger_j);
        // Measured 25 in the 12 steps.
        assert!(run.iterations <= 3 * steps, "{}", run.iterations);
    }
}

#[test]
fn a_body_that_melts_or_freezes_in_contact_with_the_gas_matches_the_modules_thermal_step() {
    // The module's latent fixture: fixed walls, the gas at 400 K, the body
    // half molten on its plateau at 350 K or at 500 K. Here the conductance
    // is a constant 2000 W/K (the module's own test makes it depend on the
    // liquid fraction, which the kernel's element does not do), and the run
    // is long enough for the body to leave its plateau: it melts through at
    // 350 K and freezes through at 500 K. Across the plateau its
    // temperature does not move, so its enthalpy is what is compared.
    let (conductance, bracket, steps) = (2000.0, 20.0, 40);
    for transition in [350.0_f64, 500.0] {
        let curve = latent(transition);
        let capacity = (100_000.0 / (transition - 200.0)).min(240_000.0 / (900.0 - transition));
        let allowed =
            Apart::allowed(&gas(0.75, 0.001), capacity, conductance, bracket, false, steps);
        let run = contact_against_the_module(
            gas(0.75, 0.001),
            &curve,
            130_000.0,
            &Store::Chart(&curve),
            conductance,
            bracket,
            None,
            steps,
        );
        eprintln!(
            "contact with a body on its plateau at {transition} K: body to {:.3} K at liquid \
             fraction {}, gas to {:.3} K, {:.3} J moved, {:.6e} J/K produced, {} iterations, \
             ledger residual {:.2e} J: {:?}; allowed {allowed:?}",
            run.body.temperature_k(),
            run.body.liquid_mass_fraction(),
            run.gas_k,
            run.moved_j,
            run.produced_j_k,
            run.iterations,
            run.ledger_j,
            run.apart,
        );
        // It left the plateau on the side the gas pushed it to: all liquid
        // above 350 K, or all solid below 500 K. The plateau held 30 J
        // either way.
        if transition < 400.0 {
            assert!(run.body.liquid_mass_fraction() >= 1.0 && run.body.temperature_k() > transition);
        } else {
            assert!(run.body.liquid_mass_fraction() <= 0.0 && run.body.temperature_k() < transition);
        }
        assert!(run.moved_j > 30.0 && run.produced_j_k > 1.0e-3);
        assert!(run.apart.within(&allowed), "{:?} against {allowed:?}", run.apart);
        assert!(run.ledger_j.abs() <= 1.0e-10, "{}", run.ledger_j);
        // Measured 77 and 69 in the 40 steps.
        assert!(run.iterations <= 3 * steps, "{}", run.iterations);
    }
}
