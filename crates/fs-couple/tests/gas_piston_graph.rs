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
//! capacity) and must agree at every step. The graph has no equivalent of
//! the module's time-resolution admission, of its fixed support, or of its
//! real-gas, reactive, nozzle and thermal-body variants.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
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
use fs_system::elements::{Force, Mass, Motion, SourceAccount, Waveform};
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
