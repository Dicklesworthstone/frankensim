//! Shared native bodies retain each original source's partial-slip memory.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use fs_contact::interface_binding::NormalContactModelSelection;
use fs_contact::interface_binding::mindlin::{SourceMindlinContact, SourceMindlinState};
use fs_contact::normal_patch::source::{PointContactGeometry, state::PointContactMaterialCards};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::tangential_friction::mindlin::multiple::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_qty::Dims;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::{matrix::MatrixTangentialConfig, multiple::MultipleTangentialConfig};
use fs_tribo::{
    InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput,
};

const DT: f64 = 1e-3;
const MASS: [f64; 3] = [0.02, 0.03, 0.04];
const IZ: [f64; 3] = [0.0003, 0.0004, 0.0005];
const ARMS: [[f64; 2]; 2] = [[0.08, -0.06], [-0.12, 0.04]];

fn source(history: &str, young: f64, mu: f64) -> SourceMindlinContact {
    let a = common::solid("island-Mindlin-A", 1000.0, young);
    let b = common::solid("island-Mindlin-B", 2000.0, 2e7);
    let normal = common::bind(&a, &b, None, 0.0);
    let original = normal.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("local-friction", Dims::NONE),
            value: PropertyValue::Scalar {
                value: mu,
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic shared-body friction fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let mut context = original.context().clone();
    context.history = history.into();
    let interface = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        context,
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let normal = PointContactMaterialCards {
        surface_a: &a,
        surface_b: &b,
        interface: &interface,
        bulk_selection: [selection.clone(), selection.clone()],
        interface_selection: selection.clone(),
        model_selection: NormalContactModelSelection::SingleSupported,
        authority: InputAuthority::SyntheticFixture,
    }
    .resolve(
        &common::query(),
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.1,
        },
        common::envelope(),
    )
    .unwrap();
    SourceMindlinContact::resolve(
        normal,
        2e-5,
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

fn controls() -> RigidMindlinIslandControls {
    RigidMindlinIslandControls {
        solve: TangentialSolveConfig {
            max_iterations: 96,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-13,
        },
        multiple: MultipleTangentialConfig {
            maximum_contacts: 8,
            maximum_matrix_terms: 256,
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
        maximum_bodies: 8,
        maximum_step_s: DT,
        maximum_slip_speed_m_s: 1.0,
        geometry_tolerance_m: 1e-12,
        velocity_tolerance_m_s: 1e-10,
        angular_velocity_tolerance_rad_s: 1e-10,
        linear_momentum_tolerance_ns: 1e-12,
        angular_momentum_tolerance_nms: 1e-12,
        energy_tolerance_j: 1e-12,
    }
}
fn setup() -> ([RigidBodyState; 3], [MassProperties; 3]) {
    let properties = [0, 1, 2].map(|i| {
        MassProperties::new(
            MASS[i],
            Vec3::ZERO,
            Vec3::new(0.8 * IZ[i], 0.9 * IZ[i], IZ[i]),
        )
        .unwrap()
    });
    let y = [0.0, ARMS[0][0] - ARMS[0][1], ARMS[1][0] - ARMS[1][1]];
    let bodies = y.map(|y| {
        RigidBodyState::new(
            Pose::new(Vec3::new(0.0, y, 0.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap()
    });
    (bodies, properties)
}
fn sites<'a>(
    laws: &'a [SourceMindlinContact; 2],
    states: &'a [SourceMindlinState; 2],
) -> [RigidMindlinIslandContact<'a>; 2] {
    [0, 1].map(|i| RigidMindlinIslandContact {
        contact_id: ["upper", "lower"][i],
        law: &laws[i],
        history: &states[i],
        body_indices: [0, i + 1],
        arms_body: ARMS[i].map(|y| Vec3::new(0.0, y, 0.0)),
        frame: TangentFrame::new([0.0, if i == 0 { 1.0 } else { -1.0 }, 0.0], [1.0, 0.0, 0.0])
            .unwrap(),
    })
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn close_vec(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}
fn input(law: &SourceMindlinContact, du: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: law.normal_response().normal_force_n,
        slip_increment_m: [du, 0.0],
        dt_s: DT,
        contact_temperature_k: Some(law.temperature_kelvin()),
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: controls().solve,
    }
}
fn named<'a>(
    result: &'a RigidMindlinIslandCandidate,
    i: usize,
) -> &'a RigidMindlinIslandContactCandidate {
    result
        .contacts
        .iter()
        .find(|c| c.contact_id == ["upper", "lower"][i])
        .unwrap()
}
fn mobility() -> [[f64; 2]; 2] {
    // Independently expanded 3-body Jacobian: both contacts share central COM
    // recoil; the opposite central arms give NEGATIVE rotational cross recoil.
    let cross = 1.0 / MASS[0] + ARMS[0][0] * ARMS[1][0] / IZ[0];
    [0, 1].map(|i| {
        [0, 1].map(|j| {
            if i == j {
                1.0 / MASS[0]
                    + 1.0 / MASS[i + 1]
                    + ARMS[i][0].powi(2) / IZ[0]
                    + ARMS[i][1].powi(2) / IZ[i + 1]
            } else {
                cross
            }
        })
    })
}
fn actual_rate(bodies: &[RigidBodyState; 3], properties: &[MassProperties; 3], i: usize) -> f64 {
    let a = bodies[0]
        .point_kinematics(properties[0], Vec3::new(0.0, ARMS[i][0], 0.0))
        .unwrap();
    let b = bodies[i + 1]
        .point_kinematics(properties[i + 1], Vec3::new(0.0, ARMS[i][1], 0.0))
        .unwrap();
    a.point_velocity_world.x - b.point_velocity_world.x
}
fn reference_drive(
    bodies: &[RigidBodyState; 3],
    properties: &[MassProperties; 3],
    du: [f64; 2],
    expected: &[TangentialStep<SourceMindlinState>; 2],
) -> [Wrench; 3] {
    let w = mobility();
    let mut drive = [Wrench::ZERO; 3];
    for i in 0..2 {
        let released = du[i]
            + 0.5
                * DT
                * DT
                * (w[i][0] * expected[0].traction_n[0] + w[i][1] * expected[1].traction_n[0]);
        let acceleration = 2.0 * (released - DT * actual_rate(bodies, properties, i)) / DT.powi(2);
        drive[i + 1].force_world = Vec3::new(-MASS[i + 1] * acceleration, 0.0, 0.0);
    }
    drive
}
fn external_work(
    before: &[RigidBodyState; 3],
    after: &[RigidBodyState],
    properties: &[MassProperties; 3],
    drive: &[Wrench; 3],
) -> f64 {
    (0..3)
        .map(|i| {
            let v0 = before[i]
                .center_of_mass_velocity_world(properties[i])
                .unwrap();
            let v1 = after[i]
                .center_of_mass_velocity_world(properties[i])
                .unwrap();
            let w0 = properties[i].angular_velocity_body(before[i].angular_momentum_body());
            let w1 = properties[i].angular_velocity_body(after[i].angular_momentum_body());
            0.5 * DT * (drive[i].force_world.dot(v0.add(v1)) + drive[i].torque_body.dot(w0.add(w1)))
        })
        .sum()
}
fn virgin_work(capacity: f64, critical: f64, u: f64) -> f64 {
    capacity * (u - 0.4 * critical * (1.0 - (1.0 - u / critical).max(0.0).powf(2.5)))
}

#[test]
fn shared_bodies_follow_independent_elastic_partial_gross_and_reversal_paths() {
    let laws = [source("upper", 1e7, 0.6), source("lower", 1.5e7, 0.45)];
    let mut states = laws.each_ref().map(TangentialLaw::initial_state);
    let (mut bodies, properties) = setup();
    let critical = laws.each_ref().map(|l| l.law().gross_slip_displacement_m());
    let capacity = laws
        .each_ref()
        .map(|l| l.law().parameters().friction_coefficient * l.normal_response().normal_force_n);
    let schedule = [
        [1e-5, 1.6e-5],
        [0.35, 0.5],
        [3.0, 2.5],
        [-0.2, -0.1],
        [-4.0, -5.0],
        [0.08, 0.11],
    ];
    let mut accumulated = [0.0; 2];
    for (phase, increments) in schedule.into_iter().enumerate() {
        let checkpoint = states.clone();
        let before = bodies;
        let du = [0, 1].map(|i| increments[i] * critical[i]);
        // Constitutive references are original prescribed paths, with zero
        // body compliance. The hand-derived simultaneous load line supplies
        // the actuator command; neither native nor generic contact solver is
        // used to construct the reference body motion or interface work.
        let expected = [0, 1].map(|i| {
            laws[i]
                .advance(&states[i], &input(&laws[i], du[i]))
                .unwrap()
        });
        let drive = reference_drive(&bodies, &properties, du, &expected);
        let contacts = sites(&laws, &states);
        let prepared = PreparedRigidMindlinIsland::new_driven(
            &bodies,
            &properties,
            &contacts,
            &drive,
            DT,
            controls(),
            &mut || false,
        )
        .unwrap();
        let w = mobility();
        let full = prepared.tangential_mobility_kg_inverse();
        // IDs sort lower then upper; both lateral rows remain physical/nonzero.
        close(full[0], w[1][1], 5e-13);
        close(full[2], w[1][0], 5e-13);
        close(full[10], w[0][0], 5e-13);
        assert!(full[5] > 0.0 && full[15] > 0.0 && full[2] > 1.0);
        let result = prepared.advance(&mut || false).unwrap();
        assert_eq!(states, checkpoint);
        assert_eq!(bodies, before);
        let q = expected.each_ref().map(|s| s.traction_n[0]);
        for i in 0..2 {
            let contact = named(&result, i);
            close(contact.contact.traction_n[0], q[i], 3e-9);
            close(contact.report.relative_displacement_m[0], du[i], 3e-14);
            close(
                contact.contact.state.endpoint_force_n(),
                expected[i].state.endpoint_force_n(),
                3e-9,
            );
            close(
                contact.contact.dissipated_j,
                expected[i].dissipated_j,
                2e-12,
            );
            close(
                laws[i].stored_energy_j(&contact.contact.state),
                laws[i].stored_energy_j(&expected[i].state),
                2e-12,
            );
            assert_eq!(
                contact.contact.state.history_intervals(),
                expected[i].state.history_intervals()
            );
            close(contact.report.contact_energy_residual_j, 0.0, 1e-12);
            if phase == 0 {
                let k = laws[i].law().initial_stiffness_n_per_m();
                close(q[i], 0.5 * k * du[i], 1e-5 * k * du[i]);
            } else if phase < 3 {
                let work = virgin_work(capacity[i], critical[i], accumulated[i] + du[i])
                    - virgin_work(capacity[i], critical[i], accumulated[i]);
                close(q[i], work / du[i], 1e-10);
            } else if phase == 3 {
                let unload = -du[i];
                let mean = capacity[i]
                    - 4.0 * virgin_work(capacity[i], critical[i], 0.5 * unload) / unload;
                close(q[i], mean, 1e-10);
                assert!(contact.contact.stored_change_j < 0.0);
                assert!(contact.report.midpoint_impulse_work_j > 0.0);
            }
            if phase == 2 {
                close(contact.contact.state.endpoint_force_n(), capacity[i], 3e-9);
            }
            if phase == 4 {
                close(contact.contact.state.endpoint_force_n(), -capacity[i], 3e-9);
            }
            accumulated[i] += du[i];
            states[i] = contact.contact.state.clone();
        }
        // One summed native impulse per body, including BOTH original contacts.
        let force = [-(q[0] + q[1]), q[0], q[1]];
        let torque = [
            ARMS[0][0] * q[0] + ARMS[1][0] * q[1],
            -ARMS[0][1] * q[0],
            -ARMS[1][1] * q[1],
        ];
        for i in 0..3 {
            close_vec(
                result.bodies_after[i].linear_momentum_world(),
                before[i]
                    .linear_momentum_world()
                    .add(drive[i].force_world.scale(DT))
                    .add(Vec3::new(DT * force[i], 0.0, 0.0)),
                1e-11,
            );
            close_vec(
                result.bodies_after[i].angular_momentum_body(),
                before[i]
                    .angular_momentum_body()
                    .add(Vec3::new(0.0, 0.0, DT * torque[i])),
                1e-12,
            );
        }
        close(
            result.report.external_work_j,
            external_work(&before, &result.bodies_after, &properties, &drive),
            1e-14,
        );
        close(result.report.total_energy_residual_j, 0.0, 1e-12);
        bodies = result.bodies_after.try_into().unwrap();
    }
    assert!(states.iter().all(|s| s.history_intervals() > 1));
}

#[test]
fn full_body_inertia_drive_work_common_acceleration_and_world_frame_are_covariant() {
    let laws = [source("upper", 1e7, 0.6), source("lower", 1.5e7, 0.45)];
    let states = laws.each_ref().map(TangentialLaw::initial_state);
    let (bodies, properties) = setup();
    let contacts = sites(&laws, &states);
    let drive = [
        Wrench {
            force_world: Vec3::new(2.0, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.03),
        },
        Wrench {
            force_world: Vec3::new(-1.0, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.0, -0.02),
        },
        Wrench {
            force_world: Vec3::new(0.2, 0.0, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.01),
        },
    ];
    let run = |b: &[RigidBodyState],
               p: &[MassProperties],
               c: &[RigidMindlinIslandContact<'_>],
               d: &[Wrench]| {
        PreparedRigidMindlinIsland::new_driven(b, p, c, d, DT, controls(), &mut || false)
            .unwrap()
            .advance(&mut || false)
            .unwrap()
    };
    let base = run(&bodies, &properties, &contacts, &drive);
    assert!(base.report.external_work_j > 0.0 && base.report.dissipated_j > 0.0);
    close(
        base.report.external_work_j,
        external_work(&bodies, &base.bodies_after, &properties, &drive),
        1e-14,
    );
    // Re-express each SAME body in an oblique local basis. I_body now contains
    // actual products of inertia; its principal frame and native arms/torques
    // must all rotate together to recover the independently known response.
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.61).unwrap();
    let [w, x, y, z] = q.components();
    let inverse = UnitQuaternion::new(w, -x, -y, -z).unwrap();
    let rebody = [0, 1, 2].map(|i| {
        RigidBodyState::new(
            Pose::new(bodies[i].pose().position_world(), q).unwrap(),
            bodies[i].linear_momentum_world(),
            inverse.rotate_body_to_world(bodies[i].angular_momentum_body()),
        )
        .unwrap()
    });
    let reprops = [0, 1, 2].map(|i| {
        MassProperties::from_principal_frame(
            MASS[i],
            Vec3::ZERO,
            Vec3::new(0.8 * IZ[i], 0.9 * IZ[i], IZ[i]),
            inverse,
        )
        .unwrap()
    });
    assert!(reprops[0].inertia_tensor_body().m[1].abs() > 1e-6);
    let mut resites = contacts;
    for site in &mut resites {
        site.arms_body = site.arms_body.map(|v| q.rotate_world_to_body(v));
    }
    let redrive = drive.map(|d| Wrench {
        force_world: d.force_world,
        torque_body: q.rotate_world_to_body(d.torque_body),
    });
    let reexpressed = run(&rebody, &reprops, &resites, &redrive);
    for i in 0..3 {
        close_vec(
            reexpressed.bodies_after[i].linear_momentum_world(),
            base.bodies_after[i].linear_momentum_world(),
            1e-12,
        );
        close_vec(
            q.rotate_body_to_world(reexpressed.bodies_after[i].angular_momentum_body()),
            base.bodies_after[i].angular_momentum_body(),
            1e-12,
        );
    }
    close(
        reexpressed.report.external_work_j,
        base.report.external_work_j,
        1e-12,
    );
    let acceleration = Vec3::new(17.0, -5.0, 4.0);
    let common = [0, 1, 2].map(|i| Wrench {
        force_world: drive[i].force_world.add(acceleration.scale(MASS[i])),
        ..drive[i]
    });
    let accelerated = run(&bodies, &properties, &contacts, &common);
    for i in 0..3 {
        close_vec(
            accelerated.bodies_after[i].linear_momentum_world(),
            base.bodies_after[i]
                .linear_momentum_world()
                .add(acceleration.scale(MASS[i] * DT)),
            1e-12,
        );
    }
    for i in 0..2 {
        close(
            named(&accelerated, i).contact.dissipated_j,
            named(&base, i).contact.dissipated_j,
            1e-12,
        );
        close(
            named(&accelerated, i).contact.state.endpoint_force_n(),
            named(&base, i).contact.state.endpoint_force_n(),
            1e-10,
        );
    }
    close(
        accelerated.report.external_work_j,
        external_work(&bodies, &accelerated.bodies_after, &properties, &common),
        1e-14,
    );
    close(accelerated.report.total_energy_residual_j, 0.0, 1e-12);
    let rotated = bodies.map(|b| {
        RigidBodyState::new(
            Pose::new(q.rotate_body_to_world(b.pose().position_world()), q).unwrap(),
            q.rotate_body_to_world(b.linear_momentum_world()),
            b.angular_momentum_body(),
        )
        .unwrap()
    });
    let rotated_drive = drive.map(|d| Wrench {
        force_world: q.rotate_body_to_world(d.force_world),
        ..d
    });
    let mut rotated_contacts = contacts;
    for site in &mut rotated_contacts {
        let to_array = |v: Vec3| [v.x, v.y, v.z];
        let n = site.frame.normal();
        let t = site.frame.longitudinal();
        site.frame = TangentFrame::new(
            to_array(q.rotate_body_to_world(Vec3::new(n[0], n[1], n[2]))),
            to_array(q.rotate_body_to_world(Vec3::new(t[0], t[1], t[2]))),
        )
        .unwrap();
    }
    let world = run(&rotated, &properties, &rotated_contacts, &rotated_drive);
    for i in 0..3 {
        close_vec(
            world.bodies_after[i].linear_momentum_world(),
            q.rotate_body_to_world(base.bodies_after[i].linear_momentum_world()),
            1e-12,
        );
        close_vec(
            world.bodies_after[i].angular_momentum_body(),
            base.bodies_after[i].angular_momentum_body(),
            1e-12,
        );
    }
    close(
        world.report.external_work_j,
        base.report.external_work_j,
        1e-12,
    );
    close_vec(
        world.report.external_angular_impulse_nms,
        q.rotate_body_to_world(base.report.external_angular_impulse_nms),
        1e-12,
    );
}

#[test]
fn canonical_contacts_global_indices_and_zero_drive_preserve_native_results() {
    let laws = [source("upper", 1e7, 0.6), source("lower", 1.5e7, 0.45)];
    let states = laws.each_ref().map(|l| {
        l.advance(
            &l.initial_state(),
            &input(l, 0.4 * l.law().gross_slip_displacement_m()),
        )
        .unwrap()
        .state
    });
    let (bodies, properties) = setup();
    let contacts = sites(&laws, &states);
    let plain = PreparedRigidMindlinIsland::new(
        &bodies,
        &properties,
        &contacts,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap()
    .advance(&mut || false)
    .unwrap();
    let zero = PreparedRigidMindlinIsland::new_driven(
        &bodies,
        &properties,
        &contacts,
        &[Wrench::ZERO; 3],
        DT,
        controls(),
        &mut || false,
    )
    .unwrap()
    .advance(&mut || false)
    .unwrap();
    assert_eq!(plain, zero);
    let reversed = PreparedRigidMindlinIsland::new(
        &bodies,
        &properties,
        &[contacts[1], contacts[0]],
        DT,
        controls(),
        &mut || false,
    )
    .unwrap()
    .advance(&mut || false)
    .unwrap();
    assert_eq!(plain, reversed);
    // Existing unrelated globals surround a reordered shared island.
    let remap = [3, 1, 4];
    let whole_bodies = [bodies[0], bodies[1], bodies[2], bodies[0], bodies[2]];
    let whole_properties = [
        properties[0],
        properties[1],
        properties[2],
        properties[0],
        properties[2],
    ];
    let mut remapped = contacts;
    for contact in &mut remapped {
        contact.body_indices = contact.body_indices.map(|i| remap[i]);
    }
    let prepared = PreparedRigidMindlinIsland::new(
        &whole_bodies,
        &whole_properties,
        &remapped,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(prepared.body_indices(), &[1, 3, 4]);
    let global = prepared.advance(&mut || false).unwrap();
    for i in 0..3 {
        let local = global
            .body_indices
            .iter()
            .position(|&index| index == remap[i])
            .unwrap();
        assert_eq!(global.bodies_after[local], plain.bodies_after[i]);
    }
    assert_eq!(global.contacts, plain.contacts);
    let mut spectator_drive = [Wrench::ZERO; 5];
    spectator_drive[2].force_world.x = 1.0;
    assert!(
        PreparedRigidMindlinIsland::new_driven(
            &whole_bodies,
            &whole_properties,
            &remapped,
            &spectator_drive,
            DT,
            controls(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn sources_unsupported_physical_recoil_budgets_and_late_cancel_publish_nothing() {
    let laws = [source("upper", 1e7, 0.6), source("lower", 1.5e7, 0.45)];
    let states = laws.each_ref().map(TangentialLaw::initial_state);
    let checkpoint = states.clone();
    let (bodies, properties) = setup();
    let original_bodies = bodies;
    let contacts = sites(&laws, &states);
    let drive = [
        Wrench {
            force_world: Vec3::new(2.0, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
        Wrench::ZERO,
    ];
    let foreign = source("different ordered source history", 1e7, 0.6).initial_state();
    let mut wrong = contacts;
    wrong[0].history = &foreign;
    assert!(
        PreparedRigidMindlinIsland::new(&bodies, &properties, &wrong, DT, controls(), &mut || {
            false
        })
        .is_err()
    );
    let mut duplicate = contacts;
    duplicate[1].contact_id = duplicate[0].contact_id;
    assert!(
        PreparedRigidMindlinIsland::new(
            &bodies,
            &properties,
            &duplicate,
            DT,
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut overlapping = [contacts[0]; 2];
    overlapping[1].contact_id = "different-ID-same-Hertz-footprint";
    assert!(matches!(
        PreparedRigidMindlinIsland::new(
            &bodies,
            &properties,
            &overlapping,
            DT,
            controls(),
            &mut || false
        ),
        Err(RigidMindlinIslandError::Invalid(
            "independent shared-body Hertz footprints are not separated"
        ))
    ));
    // Admission must examine actual native velocities, and AFTER solving it
    // must still refuse unsupported recoil from a full product-of-inertia tensor.
    for kick in [Vec3::new(0.0, 1e-5, 0.0), Vec3::new(0.0, 0.0, 1e-5)] {
        let mut moving = bodies;
        moving[0] = RigidBodyState::new(bodies[0].pose(), kick, Vec3::ZERO).unwrap();
        assert!(matches!(
            PreparedRigidMindlinIsland::new(
                &moving,
                &properties,
                &contacts,
                DT,
                controls(),
                &mut || false
            ),
            Err(RigidMindlinIslandError::UnsupportedMotion { .. })
        ));
    }
    for axis in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)] {
        let mut coupled = properties;
        coupled[0] = MassProperties::from_principal_frame(
            MASS[0],
            Vec3::ZERO,
            Vec3::new(0.8 * IZ[0], 0.9 * IZ[0], IZ[0]),
            UnitQuaternion::from_axis_angle(axis, 0.4).unwrap(),
        )
        .unwrap();
        let prepared = PreparedRigidMindlinIsland::new_driven(
            &bodies,
            &coupled,
            &contacts,
            &drive,
            DT,
            controls(),
            &mut || false,
        )
        .unwrap();
        assert!(matches!(
            prepared.advance(&mut || false),
            Err(RigidMindlinIslandError::UnsupportedMotion { .. })
        ));
    }
    let mut normal_drive = drive;
    normal_drive[0].force_world.y = 1.0;
    let prepared = PreparedRigidMindlinIsland::new_driven(
        &bodies,
        &properties,
        &contacts,
        &normal_drive,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        prepared.advance(&mut || false),
        Err(RigidMindlinIslandError::UnsupportedMotion { .. })
    ));
    let mut source_time = controls();
    source_time.maximum_step_s = 2.0 * DT;
    assert!(matches!(
        PreparedRigidMindlinIsland::new(
            &bodies,
            &properties,
            &contacts,
            2.0 * DT,
            source_time,
            &mut || false
        ),
        Err(RigidMindlinIslandError::Normal(_))
    ));
    for budget in [0, 1, 2] {
        let mut bounded = controls();
        match budget {
            0 => bounded.multiple.maximum_law_evaluations = 1,
            1 => bounded.solve.maximum_force_n = 1e-8,
            _ => bounded.multiple.maximum_sweeps = 1,
        }
        let prepared = PreparedRigidMindlinIsland::new_driven(
            &bodies,
            &properties,
            &contacts,
            &drive,
            DT,
            bounded,
            &mut || false,
        )
        .unwrap();
        assert!(prepared.advance(&mut || false).is_err());
    }
    let prepared = PreparedRigidMindlinIsland::new_driven(
        &bodies,
        &properties,
        &contacts,
        &drive,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap();
    let mut calls = 0;
    let accepted = prepared
        .advance(&mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut late = 0;
    assert!(matches!(
        prepared.advance(&mut || {
            late += 1;
            late == calls
        }),
        Err(RigidMindlinIslandError::Cancelled)
    ));
    assert_eq!(prepared.advance(&mut || false).unwrap(), accepted);
    assert_eq!(states, checkpoint);
    assert_eq!(bodies, original_bodies);
}
