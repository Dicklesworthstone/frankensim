//! G1/G3/G4: retained local microslip, shared native bodies and finite heat.
//! Explicit synthetic constant mechanical claims; caloric/conductive references
//! use independent temperature-coordinate equations, not production residuals.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use fs_contact::interface_binding::{
    NormalContactModelSelection,
    mindlin::{SourceMindlinContact, SourceMindlinState},
};
use fs_contact::normal_patch::source::{PointContactGeometry, state::PointContactMaterialCards};
use fs_couple::friction_heat::mindlin::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_la::prepared_lu::LuWorkspace;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::tangential_friction::mindlin::multiple::{
    PreparedRigidMindlinIsland, RigidMindlinIslandContact, RigidMindlinIslandControls,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_qty::Dims;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::{matrix::MatrixTangentialConfig, multiple::MultipleTangentialConfig};
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 1e-3;
const DEPTH: f64 = 2e-5;
const HALF: f64 = 0.04;
const LENGTHS: [f64; 3] = [2.0 * HALF, 0.03, 0.05];
const DENSITY: [f64; 3] = [1000.0, 1500.0, 2000.0];
const K: [f64; 3] = [100.0, 15.0, 40.0];
const IDS: [CellId; 3] = [CellId(11), CellId(29), CellId(47)];
const LINKS: [u64; 2] = [73, 31];
const NAMES: [&str; 2] = ["upper retained source", "lower retained source"];
const SIGN: [f64; 2] = [1.0, -1.0];

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let t = self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b - a) / t).ln_1p()))
    }
}
static SENSIBLE: Sensible = Sensible;

struct Fixture {
    cards: [MaterialCard; 3],
    interfaces: [InterfaceSystemCard; 2],
}
impl Fixture {
    fn new() -> Self {
        let cards = [0, 1, 2].map(|i| {
            common::solid(
                ["shared center", "top pad", "bottom pad"][i],
                DENSITY[i],
                [1e7, 2e7, 3e7][i],
            )
        });
        let interfaces = [0, 1].map(|i| {
            let normal = common::bind(&cards[0], &cards[i + 1], None, 0.0);
            let original = normal.normal_state().elastic().interface().card();
            let mut claims = original.claims().clone();
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new("local-friction", Dims::NONE),
                    value: PropertyValue::Scalar {
                        value: [0.6, 0.45][i],
                        dims: Dims::NONE,
                    },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic invariant-mechanics heat fixture".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
            let mut context = original.context().clone();
            context.history = NAMES[i].into();
            InterfaceSystemCard::assemble(
                original.surface_a().clone(),
                original.surface_b().clone(),
                context,
                claims,
                original.models().to_vec(),
            )
            .unwrap()
        });
        Self { cards, interfaces }
    }
    fn law(&self, i: usize, bulk: [f64; 2], face: f64) -> SourceMindlinContact {
        let queries = bulk.map(|t| common::query().with("T", t).unwrap());
        let interface = common::query().with("T", face).unwrap();
        let selection = MaterialPropertySelection::SingleClaimOnly;
        let normal = PointContactMaterialCards {
            surface_a: &self.cards[0],
            surface_b: &self.cards[i + 1],
            interface: &self.interfaces[i],
            bulk_selection: [selection.clone(), selection.clone()],
            interface_selection: selection.clone(),
            model_selection: NormalContactModelSelection::SingleSupported,
            authority: InputAuthority::SyntheticFixture,
        }
        .resolve_at_temperatures(
            [&queries[0], &queries[1]],
            &interface,
            PointContactGeometry::SpherePlane {
                effective_radius_m: 0.1,
            },
            common::envelope(),
        )
        .unwrap();
        SourceMindlinContact::resolve(
            normal,
            DEPTH,
            ScalarPropertyRequirement::try_new(
                "local-friction",
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            selection,
            64,
        )
        .unwrap()
    }
    fn laws(&self, thermal: &DepositionSimulation<'_>) -> [SourceMindlinContact; 2] {
        [0, 1].map(|i| {
            self.law(
                i,
                [
                    thermal.temperature_k(IDS[0]).unwrap(),
                    thermal.temperature_k(IDS[i + 1]).unwrap(),
                ],
                300.0,
            )
        })
    }
    fn source(
        &self,
        thermal: &DepositionSimulation<'_>,
        links: &[u64],
        faces: &[f64],
        outside: [f64; 3],
    ) -> Result<RigidMindlinIslandHeatSource<Vec<f64>>, &'static str> {
        Ok(RigidMindlinIslandHeatSource {
            interfaces: links
                .iter()
                .zip(faces)
                .map(|(&id, &face)| {
                    let i = LINKS.iter().position(|&x| x == id).unwrap();
                    RigidMindlinIslandHeatInterface {
                        link_id: id,
                        law: self.law(
                            i,
                            [
                                thermal.temperature_k(IDS[0]).unwrap(),
                                thermal.temperature_k(IDS[i + 1]).unwrap(),
                            ],
                            face,
                        ),
                        conductivity_w_m_k: [K[0], K[i + 1]],
                    }
                })
                .collect(),
            transport: StepCoefficients {
                sources_w: IDS.into_iter().zip(outside).collect(),
                conductances_w_k: BTreeMap::new(),
                payload: faces.to_vec(),
            },
        })
    }
}

fn area() -> f64 {
    core::f64::consts::PI * 0.1 * DEPTH
}
fn thermal(old: [f64; 3]) -> DepositionSimulation<'static> {
    let bounds = [
        [-HALF, HALF],
        [HALF, HALF + LENGTHS[1]],
        [-HALF - LENGTHS[2], -HALF],
    ];
    let mut events = Vec::new();
    for i in 0..3 {
        // The circular contact occupies half of the finite prismatic face.
        let face_area = 2.0 * area();
        let volume = face_area * (bounds[i][1] - bounds[i][0]);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: IDS[i].0,
                body_id: i as u64,
                material_id: MaterialId(i as u64 + 1),
                mass_kg: DENSITY[i] * volume,
                specific_h_j_kg: old[i] - 300.0,
                geometry: CellGeometry {
                    start_m: [0.0, bounds[i][0], 0.0],
                    end_m: [0.0, bounds[i][1], 0.0],
                    volume_m3: volume,
                    surface_area_m2: 2.0 * face_area + 4.0 * face_area.sqrt() * LENGTHS[i],
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINKS[i],
                a: IDS[0],
                b: IDS[i + 1],
                area_m2: area(),
                conductance_w_k: series_conductance_w_k(
                    area(),
                    [LENGTHS[0] / 2.0, LENGTHS[i + 1] / 2.0],
                    [K[0], K[i + 1]],
                )
                .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(1), &SENSIBLE as &dyn ThermalLaw),
            (MaterialId(2), &SENSIBLE as &dyn ThermalLaw),
            (MaterialId(3), &SENSIBLE as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 6,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn bodies() -> ([RigidBodyState; 3], [MassProperties; 3]) {
    let mass = [0.02, 0.03, 0.04];
    let properties = [0, 1, 2]
        .map(|i| MassProperties::new(mass[i], Vec3::ZERO, Vec3::new(3e-4, 4e-4, 5e-4)).unwrap());
    let positions = [0.0, HALF, -HALF];
    let bodies = [0, 1, 2].map(|i| {
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, positions[i], 0.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap()
    });
    (bodies, properties)
}
fn routes(t: &DepositionSimulation<'_>) -> [FrictionHeatRoute; 2] {
    [0, 1].map(|i| FrictionHeatRoute::from_link(t, LINKS[i], [K[0], K[i + 1]]).unwrap())
}
fn contacts<'a>(
    laws: &'a [SourceMindlinContact; 2],
    states: &'a [SourceMindlinState; 2],
    routes: &'a [FrictionHeatRoute; 2],
) -> [RigidMindlinIslandHeatContact<'a>; 2] {
    [0, 1].map(|i| RigidMindlinIslandHeatContact {
        contact: RigidMindlinIslandContact {
            contact_id: NAMES[i],
            law: &laws[i],
            history: &states[i],
            body_indices: [0, i + 1],
            arms_body: [Vec3::new(0.0, SIGN[i] * HALF, 0.0), Vec3::ZERO],
            frame: TangentFrame::new([0.0, SIGN[i], 0.0], [1.0, 0.0, 0.0]).unwrap(),
        },
        route: &routes[i],
        temperature_guess_k: 300.0,
    })
}
fn controls() -> RigidMindlinIslandHeatControls {
    RigidMindlinIslandHeatControls {
        mechanical: RigidMindlinIslandControls {
            solve: TangentialSolveConfig {
                max_iterations: 96,
                maximum_force_n: 100.0,
                force_absolute_tolerance_n: 1e-13,
                force_relative_tolerance: 1e-13,
            },
            multiple: MultipleTangentialConfig {
                maximum_contacts: 2,
                maximum_matrix_terms: 16,
                maximum_sweeps: 128,
                maximum_law_evaluations: 8192,
                local: MatrixTangentialConfig {
                    maximum_evaluations: 128,
                    displacement_tolerance_m: 1e-14,
                    energy_tolerance_j: 1e-12,
                    matrix_relative_tolerance: 1e-12,
                },
                total_energy_tolerance_j: 1e-12,
            },
            maximum_bodies: 3,
            maximum_step_s: DT,
            maximum_slip_speed_m_s: 1.0,
            geometry_tolerance_m: 1e-11,
            velocity_tolerance_m_s: 1e-10,
            angular_velocity_tolerance_rad_s: 1e-10,
            linear_momentum_tolerance_ns: 1e-12,
            angular_momentum_tolerance_nms: 1e-12,
            energy_tolerance_j: 1e-12,
        },
        maximum_cells: 3,
        maximum_unknowns: 5,
        maximum_trials: 256,
        maximum_newton_iterations: 20,
        difference_step_j: 1e-8,
        temperature_scale_j_per_k: 1e-3,
        heat: FrictionHeatTolerances {
            cell_j: 2e-11,
            thermal_j: 2e-11,
            interface_j: 2e-11,
            total_j: 3e-11,
        },
        combined_energy_tolerance_j: 4e-11,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 2e-13,
            relative_tolerance: 1e-12,
            linear_restart: 5,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
        audit_entropy: true,
    }
}
fn drive() -> [Wrench; 3] {
    [
        Wrench {
            force_world: Vec3::new(1.2, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.003),
        },
        Wrench {
            force_world: Vec3::new(-0.3, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ]
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:e} != {b:e}, tol={t:e}");
}
fn prescribed(law: &SourceMindlinContact, du: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: law.normal_response().normal_force_n,
        slip_increment_m: [du, 0.0],
        dt_s: DT,
        contact_temperature_k: Some(law.temperature_kelvin()),
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: controls().mechanical.solve,
    }
}

// Eliminate each massless interface analytically, then solve the independent
// three-capacity temperature matrix. This reference does not use heat residuals.
fn reference(
    old: &DepositionSimulation<'_>,
    heat: [f64; 2],
    outside: [f64; 3],
) -> ([f64; 3], [f64; 2]) {
    let mut matrix = [0.0; 9];
    let mut rhs = [0.0; 3];
    for i in 0..3 {
        let c = old.cells()[&IDS[i]].mass_kg / DT;
        matrix[3 * i + i] = c;
        rhs[i] = c * (old.temperature_k(IDS[i]).unwrap() - 300.0) + outside[i];
    }
    for i in 0..2 {
        let ga = area() * K[0] / (LENGTHS[0] / 2.0);
        let gb = area() * K[i + 1] / (LENGTHS[i + 1] / 2.0);
        let g = ga * gb / (ga + gb);
        let j = i + 1;
        matrix[0] += g;
        matrix[3 * j + j] += g;
        matrix[j] -= g;
        matrix[3 * j] -= g;
        rhs[0] += ga / (ga + gb) * heat[i] / DT;
        rhs[j] += gb / (ga + gb) * heat[i] / DT;
    }
    let mut lu = LuWorkspace::new(3).unwrap();
    lu.factor(&matrix).unwrap();
    let mut rise = [0.0; 3];
    lu.solve_factored_into(&rhs, &mut rise).unwrap();
    let temperatures = rise.map(|r| 300.0 + r);
    let faces = [0, 1].map(|i| {
        let ga = area() * K[0] / (LENGTHS[0] / 2.0);
        let gb = area() * K[i + 1] / (LENGTHS[i + 1] / 2.0);
        (heat[i] / DT + ga * temperatures[0] + gb * temperatures[i + 1]) / (ga + gb)
    });
    (temperatures, faces)
}

#[test]
fn shared_retained_microslip_and_external_work_close_the_independent_finite_heat_equations() {
    let f = Fixture::new();
    let old = thermal([299.0, 301.0, 303.0]);
    let laws = f.laws(&old);
    let states = laws.each_ref().map(TangentialLaw::initial_state);
    let (bodies, properties) = bodies();
    let routes = routes(&old);
    let contacts = contacts(&laws, &states, &routes);
    let outside = [0.001, -0.002, 0.0015];
    let external = drive();
    let native_contacts = contacts.map(|c| c.contact);
    let expected = PreparedRigidMindlinIsland::new_driven(
        &bodies,
        &properties,
        &native_contacts,
        &external,
        DT,
        controls().mechanical,
        &mut || false,
    )
    .unwrap()
    .advance(&mut || false)
    .unwrap();
    let result = propose_rigid_mindlin_island_driven_heat(
        &bodies,
        &properties,
        &contacts,
        &external,
        &old,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |t, ids, faces| f.source(t, ids, faces, outside),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(result.mechanical.bodies_after, expected.bodies_after);
    assert_eq!(result.mechanical.report, expected.report);
    let dissipation = [0, 1].map(|i| {
        expected
            .contacts
            .iter()
            .find(|c| c.contact_id == NAMES[i])
            .unwrap()
            .contact
            .dissipated_j
    });
    assert!(dissipation.iter().all(|d| *d > 0.0));
    let (temperatures, faces) = reference(&old, dissipation, outside);
    for i in 0..3 {
        close(
            result.thermal.temperature_k(IDS[i]).unwrap(),
            temperatures[i],
            2e-7,
        );
    }
    for i in 0..2 {
        let native = result
            .mechanical
            .contacts
            .iter()
            .find(|c| c.contact_id == NAMES[i])
            .unwrap();
        let baseline = expected
            .contacts
            .iter()
            .find(|c| c.contact_id == NAMES[i])
            .unwrap();
        assert_eq!(
            native.contact.state.constitutive_state(),
            baseline.contact.state.constitutive_state()
        );
        let heat = result
            .heat
            .interfaces
            .iter()
            .find(|r| r.link_id == LINKS[i])
            .unwrap();
        close(heat.interface_temperature_k, faces[i], 2e-7);
        close(heat.body_heat_j.iter().sum(), dissipation[i], 2e-11);
        assert_eq!(heat.dissipated_j, native.contact.dissipated_j);
        let source = result
            .sources
            .iter()
            .find(|s| s.contact_id == NAMES[i])
            .unwrap();
        assert_eq!(
            source.law.temperature_kelvin(),
            heat.interface_temperature_k
        );
        assert!(
            source
                .law
                .stored_energy_j(&native.contact.state)
                .is_finite()
        );
        assert!(laws[i].stored_energy_j(&native.contact.state).is_nan());
    }
    assert!(
        result
            .heat
            .interfaces
            .iter()
            .any(|h| h.body_heat_j.iter().any(|q| *q < 0.0))
    );
    close(
        result.heat.external_heat_j,
        DT * outside.iter().sum::<f64>(),
        1e-14,
    );
    let caloric: f64 = IDS
        .into_iter()
        .map(|id| {
            old.cells()[&id].mass_kg
                * (result.thermal.cells()[&id].specific_h_j_kg - old.cells()[&id].specific_h_j_kg)
        })
        .sum();
    close(
        expected.report.kinetic_change_j + expected.report.stored_change_j + caloric
            - expected.report.external_work_j
            - result.heat.external_heat_j,
        0.0,
        4e-11,
    );
    close(result.combined_energy_residual_j, 0.0, 4e-11);
    assert!(result.heat.entropy.is_some());
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(states, laws.each_ref().map(TangentialLaw::initial_state));
}

#[test]
fn repeated_temperature_updates_keep_loaded_reversal_memory_and_zero_drive_equivalence() {
    let f = Fixture::new();
    let mut thermal = thermal([300.0; 3]);
    let mut laws = f.laws(&thermal);
    let mut states = laws.each_ref().map(|law| {
        law.advance(
            &law.initial_state(),
            &prescribed(law, 0.7 * law.law().gross_slip_displacement_m()),
        )
        .unwrap()
        .state
    });
    let (mut bodies, properties) = bodies();
    for phase in 0..3 {
        let routes = routes(&thermal);
        let contacts = contacts(&laws, &states, &routes);
        let before = states.clone();
        let plain = propose_rigid_mindlin_island_heat(
            &bodies,
            &properties,
            &contacts,
            &thermal,
            thermal.time_s(),
            DT,
            controls(),
            &[0.0; 3],
            |t, ids, faces| f.source(t, ids, faces, [0.0; 3]),
            &CancelGate::new(),
        )
        .unwrap();
        let zero = propose_rigid_mindlin_island_driven_heat(
            &bodies,
            &properties,
            &contacts,
            &[Wrench::ZERO; 3],
            &thermal,
            thermal.time_s(),
            DT,
            controls(),
            &[0.0; 3],
            |t, ids, faces| f.source(t, ids, faces, [0.0; 3]),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(plain.mechanical, zero.mechanical);
        assert_eq!(plain.thermal.cells(), zero.thermal.cells());
        assert_eq!(states, before);
        if phase == 0 {
            assert!(plain.mechanical.report.stored_change_j < 0.0);
            assert!(plain.mechanical.report.kinetic_change_j > 0.0);
        }
        assert!(
            plain
                .mechanical
                .contacts
                .iter()
                .all(|c| c.contact.state.history_intervals() > 0)
        );
        for i in 0..2 {
            let source = plain
                .sources
                .iter()
                .find(|s| s.contact_id == NAMES[i])
                .unwrap();
            let contact = plain
                .mechanical
                .contacts
                .iter()
                .find(|c| c.contact_id == NAMES[i])
                .unwrap();
            laws[i] = source.law.clone();
            states[i] = contact.contact.state.clone();
        }
        bodies = plain.mechanical.bodies_after.try_into().unwrap();
        thermal = plain.thermal;
    }
    close(thermal.time_s(), 3.0 * DT, 1e-17);
}

#[test]
fn canonical_interfaces_common_acceleration_and_source_temperatures_preserve_physics() {
    let f = Fixture::new();
    let old = thermal([300.0, 301.0, 302.0]);
    let laws = f.laws(&old);
    let states = laws.each_ref().map(TangentialLaw::initial_state);
    let (bodies, properties) = bodies();
    let routes = routes(&old);
    let contacts = contacts(&laws, &states, &routes);
    let run = |c: &[RigidMindlinIslandHeatContact<'_>], w: &[Wrench]| {
        propose_rigid_mindlin_island_driven_heat(
            &bodies,
            &properties,
            c,
            w,
            &old,
            0.0,
            DT,
            controls(),
            &[0.0; 3],
            |t, ids, faces| f.source(t, ids, faces, [0.0; 3]),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let base = run(&contacts, &drive());
    let permuted = run(&[contacts[1], contacts[0]], &drive());
    assert_eq!(base.mechanical, permuted.mechanical);
    assert_eq!(base.thermal.cells(), permuted.thermal.cells());
    let acceleration = Vec3::new(3.0, -2.0, 1.0);
    let mut common = drive();
    for i in 0..3 {
        common[i].force_world = common[i]
            .force_world
            .add(acceleration.scale(properties[i].mass()));
    }
    let boosted = run(&contacts, &common);
    for i in 0..2 {
        close(
            boosted.heat.interfaces[i].dissipated_j,
            base.heat.interfaces[i].dissipated_j,
            1e-13,
        );
    }
    for i in 0..3 {
        close(
            boosted.thermal.temperature_k(IDS[i]).unwrap(),
            base.thermal.temperature_k(IDS[i]).unwrap(),
            2e-9,
        );
    }
    close(boosted.combined_energy_residual_j, 0.0, 4e-11);
    // Every endpoint source must use the receiving solid's own temperature.
    let failure = propose_rigid_mindlin_island_driven_heat(
        &bodies,
        &properties,
        &contacts,
        &drive(),
        &old,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |t, ids, faces| {
            let mut source = f.source(t, ids, faces, [0.0; 3])?;
            source.interfaces[0].law = f.law(1, [300.0, 300.0], faces[0]);
            Ok::<_, &'static str>(source)
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failure,
        Err(RigidMindlinIslandHeatError::Invalid(_))
    ));
}

#[test]
fn source_geometry_budget_and_final_cancellation_refuse_without_partial_publication() {
    let f = Fixture::new();
    let old = thermal([300.0; 3]);
    let laws = f.laws(&old);
    let states = laws.each_ref().map(TangentialLaw::initial_state);
    let (bodies, properties) = bodies();
    let saved = old.cells().clone();
    let checkpoint = states.clone();
    let routes = routes(&old);
    let contacts = contacts(&laws, &states, &routes);
    let calls = Cell::new(0);
    let result = propose_rigid_mindlin_island_driven_heat(
        &bodies,
        &properties,
        &contacts,
        &drive(),
        &old,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |t, ids, faces| {
            calls.set(calls.get() + 1);
            f.source(t, ids, faces, [0.0; 3])
        },
        &CancelGate::new(),
    )
    .unwrap();
    let total = calls.get();
    assert_eq!(total, result.trials);
    assert!(total > 1);
    let gate = CancelGate::new();
    let calls = Cell::new(0);
    let late = propose_rigid_mindlin_island_driven_heat(
        &bodies,
        &properties,
        &contacts,
        &drive(),
        &old,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |t, ids, faces| {
            calls.set(calls.get() + 1);
            let source = f.source(t, ids, faces, [0.0; 3])?;
            if calls.get() == total {
                gate.request();
            }
            Ok::<_, &'static str>(source)
        },
        &gate,
    );
    assert!(late.is_err());
    assert_eq!(calls.get(), total);
    let mut bounded = controls();
    bounded.maximum_trials = 1;
    assert!(matches!(
        propose_rigid_mindlin_island_driven_heat(
            &bodies,
            &properties,
            &contacts,
            &drive(),
            &old,
            0.0,
            DT,
            bounded,
            &[0.0; 3],
            |t, ids, faces| f.source(t, ids, faces, [0.0; 3]),
            &CancelGate::new()
        ),
        Err(RigidMindlinIslandHeatError::TrialBudget)
    ));
    let mut wrong = contacts;
    wrong[0].route = &routes[1];
    assert!(
        propose_rigid_mindlin_island_heat(
            &bodies,
            &properties,
            &wrong,
            &old,
            0.0,
            DT,
            controls(),
            &[0.0; 3],
            |t, ids, faces| f.source(t, ids, faces, [0.0; 3]),
            &CancelGate::new()
        )
        .is_err()
    );
    let foreign = Fixture::new();
    let mut changed = foreign.interfaces[0].context().clone();
    changed.history = "foreign heat source".into();
    let mut foreign = foreign;
    foreign.interfaces[0] = InterfaceSystemCard::assemble(
        foreign.interfaces[0].surface_a().clone(),
        foreign.interfaces[0].surface_b().clone(),
        changed,
        foreign.interfaces[0].claims().clone(),
        foreign.interfaces[0].models().to_vec(),
    )
    .unwrap();
    assert!(matches!(
        propose_rigid_mindlin_island_driven_heat(
            &bodies,
            &properties,
            &contacts,
            &drive(),
            &old,
            0.0,
            DT,
            controls(),
            &[0.0; 3],
            |t, ids, faces| foreign.source(t, ids, faces, [0.0; 3]),
            &CancelGate::new()
        ),
        Err(RigidMindlinIslandHeatError::MindlinSource(_))
    ));
    assert_eq!(old.cells(), &saved);
    assert_eq!(states, checkpoint);
    assert_eq!(old.time_s(), 0.0);
    let retry = propose_rigid_mindlin_island_driven_heat(
        &bodies,
        &properties,
        &contacts,
        &drive(),
        &old,
        0.0,
        DT,
        controls(),
        &[0.0; 3],
        |t, ids, faces| f.source(t, ids, faces, [0.0; 3]),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.mechanical, result.mechanical);
    assert_eq!(retry.thermal.cells(), result.thermal.cells());
}
