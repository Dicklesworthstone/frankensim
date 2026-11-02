//! G1/G3/G4: shared normal dynamics and retained shear use the same native work.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use fs_contact::interface_binding::mindlin::loading::{
    SourceHertzMindlinContact, SourceHertzMindlinState,
    path::{SourceMindlinPathInput, SourceMindlinPathStep},
};
use fs_contact::interface_binding::{NormalContactModelSelection, mindlin::SourceMindlinContact};
use fs_contact::normal_patch::source::{PointContactGeometry, state::PointContactMaterialCards};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::tangential_friction::mindlin::multiple::loading::*;
use fs_mbd::tangential_friction::mindlin::multiple::{
    RigidMindlinIslandControls, RigidMindlinIslandError,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_qty::Dims;
use fs_solver::{Globalization, LineSearchConfig, NewtonKrylovConfig};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::{matrix::MatrixTangentialConfig, multiple::MultipleTangentialConfig};
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig};

const DT: f64 = 1e-3;
const MASS: [f64; 3] = [0.02, 0.03, 0.04];
const IZ: [f64; 3] = [0.0003, 0.0004, 0.0005];
const ARMS: [[f64; 2]; 2] = [[0.08, -0.06], [-0.12, 0.04]];

fn source(history: &str, young: f64, mu: f64) -> SourceHertzMindlinContact {
    let a = common::solid("changing-island-A", 1000.0, young);
    let b = common::solid("changing-island-B", 2000.0, 2e7);
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
                source: "synthetic shared normal/shear fixture".into(),
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
    SourceHertzMindlinContact::from_fixed_source(
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
            128,
        )
        .unwrap(),
    )
    .unwrap()
}
fn controls() -> RigidHertzMindlinIslandControls {
    RigidHertzMindlinIslandControls {
        native: RigidMindlinIslandControls {
            solve: TangentialSolveConfig {
                max_iterations: 96,
                maximum_force_n: 100.0,
                force_absolute_tolerance_n: 1e-13,
                force_relative_tolerance: 1e-13,
            },
            multiple: MultipleTangentialConfig {
                maximum_contacts: 8,
                maximum_matrix_terms: 576,
                maximum_sweeps: 128,
                maximum_law_evaluations: 16384,
                local: MatrixTangentialConfig {
                    maximum_evaluations: 128,
                    displacement_tolerance_m: 1e-14,
                    energy_tolerance_j: 2e-12,
                    matrix_relative_tolerance: 1e-12,
                },
                total_energy_tolerance_j: 2e-12,
            },
            maximum_bodies: 8,
            maximum_step_s: DT,
            maximum_slip_speed_m_s: 1.0,
            geometry_tolerance_m: 1e-12,
            velocity_tolerance_m_s: 1e-10,
            angular_velocity_tolerance_rad_s: 1e-10,
            linear_momentum_tolerance_ns: 1e-12,
            angular_momentum_tolerance_nms: 1e-12,
            energy_tolerance_j: 2e-12,
        },
        normal_solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-16,
            relative_tolerance: 0.0,
            linear_restart: 8,
            max_linear_cycles: 4,
            forcing_minimum: 1e-12,
            forcing_maximum: 0.1,
            forcing_gamma: 0.9,
            forcing_exponent: 1.5,
            globalization: Globalization::LineSearch(LineSearchConfig {
                armijo: 1e-4,
                contraction: 0.5,
                minimum_step: 1.0 / 65536.0,
            }),
        },
        maximum_normal_iterations: 48,
        maximum_normal_evaluations: 4096,
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
    laws: &'a [SourceHertzMindlinContact; 2],
    states: &'a [SourceHertzMindlinState; 2],
) -> [RigidHertzMindlinIslandContact<'a>; 2] {
    [0, 1].map(|i| RigidHertzMindlinIslandContact {
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
fn input(indentation: f64, shear: f64) -> SourceMindlinPathInput {
    SourceMindlinPathInput {
        indentation_m: indentation,
        free_slip_increment_m: shear,
        dt_s: DT,
        compliance_m_per_n: 0.0,
        solve: controls().native.solve,
    }
}
fn named(
    result: &RigidHertzMindlinIslandCandidate,
    i: usize,
) -> &RigidHertzMindlinIslandContactCandidate {
    result
        .contacts
        .iter()
        .find(|c| c.contact_id == ["upper", "lower"][i])
        .unwrap()
}
fn normal_matrix() -> [[f64; 2]; 2] {
    [
        [1.0 / MASS[0] + 1.0 / MASS[1], -1.0 / MASS[0]],
        [-1.0 / MASS[0], 1.0 / MASS[0] + 1.0 / MASS[2]],
    ]
}
fn shear_matrix() -> [[f64; 2]; 2] {
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
fn hertz_mean(law: &SourceHertzMindlinContact, before: f64, after: f64) -> f64 {
    let k = law.fixed_source().normal_source().stiffness_n_per_m_3_2();
    if (after - before).abs() < 1e-8 * before.max(after) {
        // Independent close-endpoint limit; the midpoint error is O(delta^2)
        // and below floating-point resolution in this branch. Subtracting
        // rounded endpoint energies would destroy the reference at one ulp.
        k * (0.5 * before + 0.5 * after).powf(1.5)
    } else {
        0.4 * k * (after.powf(2.5) - before.powf(2.5)) / (after - before)
    }
}
fn reference_drive(
    bodies: &[RigidBodyState; 3],
    properties: &[MassProperties; 3],
    contacts: &[RigidHertzMindlinIslandContact<'_>; 2],
    endpoints: [f64; 2],
    du: [f64; 2],
    expected: &[SourceMindlinPathStep; 2],
) -> ([Wrench; 3], [f64; 2]) {
    let wn = normal_matrix();
    let wt = shear_matrix();
    let means = [0, 1].map(|i| {
        hertz_mean(
            contacts[i].law,
            contacts[i].history.indentation_m(),
            endpoints[i],
        )
    });
    let q = expected
        .each_ref()
        .map(|step| step.tangential.traction_n[0]);
    let mut drive = [Wrench::ZERO; 3];
    for i in 0..2 {
        let a = bodies[0]
            .point_kinematics(properties[0], contacts[i].arms_body[0])
            .unwrap();
        let b = bodies[i + 1]
            .point_kinematics(properties[i + 1], contacts[i].arms_body[1])
            .unwrap();
        let velocity = a.point_velocity_world.sub(b.point_velocity_world);
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let free_n = endpoints[i] - contacts[i].history.indentation_m()
            + 0.5 * DT * DT * (wn[i][0] * means[0] + wn[i][1] * means[1]);
        let free_t = du[i] + 0.5 * DT * DT * (wt[i][0] * q[0] + wt[i][1] * q[1]);
        let an = 2.0 * (free_n - DT * sign * velocity.y) / (DT * DT);
        let at = 2.0 * (free_t - DT * velocity.x) / (DT * DT);
        drive[i + 1].force_world = Vec3::new(-MASS[i + 1] * at, -MASS[i + 1] * sign * an, 0.0);
    }
    (drive, means)
}
fn external_work(
    before: &[RigidBodyState; 3],
    after: &[RigidBodyState],
    properties: &[MassProperties; 3],
    drive: &[Wrench; 3],
    dt: f64,
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
            0.5 * dt * (drive[i].force_world.dot(v0.add(v1)) + drive[i].torque_body.dot(w0.add(w1)))
        })
        .sum()
}

#[test]
fn independent_shared_normal_and_shear_loadlines_preserve_unloading_and_reversal_memory() {
    let laws = [source("upper", 1e7, 0.6), source("lower", 1.5e7, 0.45)];
    let mut states = laws.each_ref().map(TangentialLaw::initial_state);
    let (mut bodies, properties) = setup();
    let critical = laws
        .each_ref()
        .map(|l| l.fixed_source().law().gross_slip_displacement_m());
    let schedule = [
        ([2e-5, 2e-5], [0.3, 0.4]),
        ([4e-5, 3e-5], [1.2, 0.9]),
        ([8e-6, 1e-5], [0.0, 0.0]),
        ([3e-5, 2.5e-5], [-3.5, -2.5]),
        ([3e-5, 2.5e-5], [0.2, 0.3]),
    ];
    for (phase, (depth, fractions)) in schedule.into_iter().enumerate() {
        let old = states.clone();
        let before = bodies;
        let du = [0, 1].map(|i| fractions[i] * critical[i]);
        // Independent prescribed source histories; the hand-expanded native
        // mass/inertia equations supply drive, not either simultaneous solver.
        let expected = [0, 1].map(|i| {
            laws[i]
                .advance_coupled_path(&states[i], &input(depth[i], du[i]), &mut || false)
                .unwrap()
        });
        let contacts = sites(&laws, &states);
        let (drive, means) = reference_drive(&bodies, &properties, &contacts, depth, du, &expected);
        let prepared = PreparedRigidHertzMindlinIsland::new_driven(
            &bodies,
            &properties,
            &contacts,
            &drive,
            DT,
            controls(),
            &mut || false,
        )
        .unwrap();
        let physical = prepared.contact_mobility_kg_inverse();
        // Canonical lower/upper order, with N/T/L coordinates per patch.
        close(physical[0], normal_matrix()[1][1], 5e-13);
        close(physical[3], normal_matrix()[1][0], 5e-13);
        close(physical[7], shear_matrix()[1][1], 5e-13);
        close(physical[10], shear_matrix()[1][0], 5e-13);
        assert!(physical[14] > 0.0 && physical[35] > 0.0);
        let result = prepared.advance(&mut || false).unwrap();
        assert_eq!(states, old);
        assert_eq!(bodies, before);
        let q = expected.each_ref().map(|x| x.tangential.traction_n[0]);
        for i in 0..2 {
            let answer = named(&result, i);
            close(answer.normal_average_force_n, means[i], 2e-11);
            close(
                answer.contact.tangential.state.indentation_m(),
                depth[i],
                2e-14,
            );
            close(
                answer.relative_displacement_m[0],
                depth[i] - states[i].indentation_m(),
                2e-14,
            );
            close(answer.relative_displacement_m[1], du[i], 2e-14);
            close(answer.contact.tangential.traction_n[0], q[i], 3e-9);
            close(
                answer.contact.tangential.state.endpoint_force_n(),
                expected[i].tangential.state.endpoint_force_n(),
                3e-9,
            );
            close(
                answer.contact.tangential.dissipated_j,
                expected[i].tangential.dissipated_j,
                2e-12,
            );
            assert_eq!(
                answer.contact.tangential.state.history_intervals(),
                expected[i].tangential.state.history_intervals(),
            );
            close(
                answer.contact.normal_work_j,
                expected[i].normal_work_j,
                2e-12,
            );
            assert_eq!(
                answer.contact.normal_response.authority,
                InputAuthority::SyntheticFixture
            );
            close(answer.contact_energy_residual_j, 0.0, 2e-12);
            if phase == 0 {
                let fraction = fractions[i];
                let capacity =
                    laws[i].law().parameters().friction_coefficient * states[i].normal_force_n();
                let mean = capacity * (1.0 - 0.4 * (1.0 - (1.0 - fraction).powf(2.5)) / fraction);
                close(q[i], mean, 2e-12);
            }
            if phase == 2 {
                assert!(answer.contact.tangential.dissipated_j > 0.0);
                assert!(answer.contact.tangential.stored_change_j < 0.0);
                assert!(answer.contact.normal_work_j < 0.0);
                close(answer.relative_displacement_m[1], 0.0, 2e-14);
            }
            states[i] = answer.contact.tangential.state.clone();
        }
        let force = [
            Vec3::new(-q[0] - q[1], means[1] - means[0], 0.0),
            Vec3::new(q[0], means[0], 0.0),
            Vec3::new(q[1], -means[1], 0.0),
        ];
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
                    .add(drive[i].force_world.add(force[i]).scale(DT)),
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
            external_work(&before, &result.bodies_after, &properties, &drive, DT),
            2e-13,
        );
        close(result.report.total_energy_residual_j, 0.0, 2e-12);
        bodies = result.bodies_after.try_into().unwrap();
    }
    assert!(states.iter().all(|s| s.history_intervals() > 1));
}

// Separate classical coupled-Hertz ODE; no source mean, native phase or solver
// is used by this continuous-time reference.
fn ode_rhs(state: [f64; 4], stiffness: [f64; 2]) -> [f64; 4] {
    let force = [0, 1].map(|i| stiffness[i] * state[i].powf(1.5));
    let w = normal_matrix();
    [
        state[2],
        state[3],
        -w[0][0] * force[0] - w[0][1] * force[1],
        -w[1][0] * force[0] - w[1][1] * force[1],
    ]
}
fn ode(stiffness: [f64; 2], duration: f64) -> [f64; 4] {
    let mut x = [2e-5, 2e-5, 0.0, 0.0];
    let h = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = ode_rhs(x, stiffness);
        let k2 = ode_rhs([0, 1, 2, 3].map(|i| x[i] + 0.5 * h * k1[i]), stiffness);
        let k3 = ode_rhs([0, 1, 2, 3].map(|i| x[i] + 0.5 * h * k2[i]), stiffness);
        let k4 = ode_rhs([0, 1, 2, 3].map(|i| x[i] + h * k3[i]), stiffness);
        x = [0, 1, 2, 3].map(|i| x[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    x
}
#[test]
fn coupled_normal_dynamics_refine_against_independent_ode_and_close_energy() {
    let laws = [source("upper", 1e7, 0.0), source("lower", 1.5e7, 0.0)];
    let stiffness = laws
        .each_ref()
        .map(|l| l.fixed_source().normal_source().stiffness_n_per_m_3_2());
    let duration = 4e-4;
    let exact = ode(stiffness, duration);
    let errors = [8, 16, 32].map(|steps| {
        let (mut bodies, properties) = setup();
        let mut states = laws.each_ref().map(TangentialLaw::initial_state);
        let initial_energy: f64 = states.iter().map(|s| s.normal_energy_j()).sum();
        let dt = duration / steps as f64;
        for _ in 0..steps {
            let mut contacts = sites(&laws, &states);
            for (i, contact) in contacts.iter_mut().enumerate() {
                let a = bodies[0].pose().position_world();
                let b = bodies[i + 1].pose().position_world();
                let common = a.scale(0.5).add(b.scale(0.5));
                contact.arms_body = [common.sub(a), common.sub(b)];
            }
            let result = PreparedRigidHertzMindlinIsland::new(
                &bodies,
                &properties,
                &contacts,
                dt,
                controls(),
                &mut || false,
            )
            .unwrap()
            .advance(&mut || false)
            .unwrap();
            close(result.report.total_energy_residual_j, 0.0, 1e-13);
            assert_eq!(result.report.dissipated_j, 0.0);
            assert_eq!(result.report.external_work_j, 0.0);
            states = [0, 1].map(|i| named(&result, i).contact.tangential.state.clone());
            // Transport COM positions using actual midpoint velocities. Recenter
            // the normal-only attachment points on their evolving common line;
            // central Hertz normal mobility is independent of those arm lengths.
            bodies = [0, 1, 2].map(|i| {
                let a = bodies[i]
                    .center_of_mass_velocity_world(properties[i])
                    .unwrap();
                let b = result.bodies_after[i]
                    .center_of_mass_velocity_world(properties[i])
                    .unwrap();
                RigidBodyState::new(
                    Pose::new(
                        bodies[i]
                            .pose()
                            .position_world()
                            .add(a.add(b).scale(0.5 * dt)),
                        bodies[i].pose().orientation(),
                    )
                    .unwrap(),
                    result.bodies_after[i].linear_momentum_world(),
                    result.bodies_after[i].angular_momentum_body(),
                )
                .unwrap()
            });
        }
        let mut final_energy: f64 = states.iter().map(|s| s.normal_energy_j()).sum();
        for i in 0..3 {
            final_energy += bodies[i].kinetic_energy(properties[i]).unwrap();
        }
        close(final_energy, initial_energy, 2e-13);
        let maximum = [0, 1]
            .map(|i| (states[i].indentation_m() - exact[i]).abs())
            .into_iter()
            .fold(0.0, f64::max);
        assert!(states.iter().all(|s| s.indentation_m() > 0.0));
        maximum
    });
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 1e-10, "{errors:?}");
}

#[test]
fn source_driven_normal_shear_island_is_covariant_and_canonical() {
    let laws = [source("upper", 1e7, 0.6), source("lower", 1.5e7, 0.45)];
    let states = laws.each_ref().map(TangentialLaw::initial_state);
    let (bodies, properties) = setup();
    let contacts = sites(&laws, &states);
    let drive = [
        Wrench {
            force_world: Vec3::new(0.2, 0.01, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.001),
        },
        Wrench {
            force_world: Vec3::new(-0.3, -0.1, 0.0),
            torque_body: Vec3::new(0.0, 0.0, -0.001),
        },
        Wrench {
            force_world: Vec3::new(-0.2, 0.1, 0.0),
            torque_body: Vec3::ZERO,
        },
    ];
    let run = |b: &[RigidBodyState],
               p: &[MassProperties],
               c: &[RigidHertzMindlinIslandContact<'_>],
               d: &[Wrench]| {
        PreparedRigidHertzMindlinIsland::new_driven(b, p, c, d, DT, controls(), &mut || false)
            .unwrap()
            .advance(&mut || false)
            .unwrap()
    };
    let base = run(&bodies, &properties, &contacts, &drive);
    close(
        base.report.external_work_j,
        external_work(&bodies, &base.bodies_after, &properties, &drive, DT),
        1e-13,
    );
    assert_ne!(base.report.normal_stored_change_j, 0.0);
    assert!(base.report.dissipated_j > 0.0);
    assert_eq!(
        base,
        run(&bodies, &properties, &[contacts[1], contacts[0]], &drive)
    );
    let acceleration = Vec3::new(7.0, -3.0, 2.0);
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
            named(&accelerated, i)
                .contact
                .tangential
                .state
                .indentation_m(),
            named(&base, i).contact.tangential.state.indentation_m(),
            1e-14,
        );
        close(
            named(&accelerated, i)
                .contact
                .tangential
                .state
                .endpoint_force_n(),
            named(&base, i).contact.tangential.state.endpoint_force_n(),
            1e-10,
        );
    }
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.61).unwrap();
    let [w, x, y, z] = q.components();
    let inverse = UnitQuaternion::new(w, -x, -y, -z).unwrap();
    let rebody = bodies.map(|b| {
        RigidBodyState::new(
            Pose::new(b.pose().position_world(), q).unwrap(),
            b.linear_momentum_world(),
            q.rotate_world_to_body(b.angular_momentum_body()),
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
    for c in &mut resites {
        c.arms_body = c.arms_body.map(|v| q.rotate_world_to_body(v));
    }
    let redrive = drive.map(|d| Wrench {
        force_world: d.force_world,
        torque_body: q.rotate_world_to_body(d.torque_body),
    });
    let reexpressed = run(&rebody, &reprops, &resites, &redrive);
    let world_bodies = bodies.map(|b| {
        RigidBodyState::new(
            Pose::new(q.rotate_body_to_world(b.pose().position_world()), q).unwrap(),
            q.rotate_body_to_world(b.linear_momentum_world()),
            b.angular_momentum_body(),
        )
        .unwrap()
    });
    let world_drive = drive.map(|d| Wrench {
        force_world: q.rotate_body_to_world(d.force_world),
        ..d
    });
    let mut world_contacts = contacts;
    let array = |v: Vec3| [v.x, v.y, v.z];
    for c in &mut world_contacts {
        let n = c.frame.normal();
        let t = c.frame.longitudinal();
        c.frame = TangentFrame::new(
            array(q.rotate_body_to_world(Vec3::new(n[0], n[1], n[2]))),
            array(q.rotate_body_to_world(Vec3::new(t[0], t[1], t[2]))),
        )
        .unwrap();
    }
    let world = run(&world_bodies, &properties, &world_contacts, &world_drive);
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
        close_vec(
            world.bodies_after[i].linear_momentum_world(),
            q.rotate_body_to_world(base.bodies_after[i].linear_momentum_world()),
            1e-12,
        );
    }
    close(
        world.report.external_work_j,
        base.report.external_work_j,
        1e-12,
    );
    close(reexpressed.report.total_energy_residual_j, 0.0, 2e-12);
    let remap = [3, 1, 4];
    let all_bodies = [bodies[0], bodies[1], bodies[2], bodies[0], bodies[2]];
    let all_props = [
        properties[0],
        properties[1],
        properties[2],
        properties[0],
        properties[2],
    ];
    let mut all_drive = [Wrench::ZERO; 5];
    for i in 0..3 {
        all_drive[remap[i]] = drive[i];
    }
    let mut mapped = contacts;
    for c in &mut mapped {
        c.body_indices = c.body_indices.map(|i| remap[i]);
    }
    let global = run(&all_bodies, &all_props, &mapped, &all_drive);
    assert_eq!(global.body_indices, vec![1, 3, 4]);
    for i in 0..3 {
        let j = global
            .body_indices
            .iter()
            .position(|&index| index == remap[i])
            .unwrap();
        assert_eq!(global.bodies_after[j], base.bodies_after[i]);
    }
    assert_eq!(global.contacts, base.contacts);
    let plain = PreparedRigidHertzMindlinIsland::new(
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
    assert_eq!(
        plain,
        run(&bodies, &properties, &contacts, &[Wrench::ZERO; 3])
    );
}

#[test]
fn source_domains_physical_coupling_footprint_growth_and_late_failure_are_atomic() {
    let laws = [source("upper", 1e7, 0.6), source("lower", 1.5e7, 0.45)];
    let states = laws.each_ref().map(TangentialLaw::initial_state);
    let saved = states.clone();
    let (bodies, properties) = setup();
    let contacts = sites(&laws, &states);
    let drive = [
        Wrench {
            force_world: Vec3::new(0.2, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::new(-0.3, -0.1, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::new(-0.2, 0.1, 0.0),
            torque_body: Vec3::ZERO,
        },
    ];
    let original = bodies;
    let foreign = source("foreign original interface", 1e7, 0.6).initial_state();
    let mut wrong = contacts;
    wrong[0].history = &foreign;
    assert!(
        PreparedRigidHertzMindlinIsland::new(
            &bodies,
            &properties,
            &wrong,
            DT,
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut coupled = contacts;
    coupled[0].frame = TangentFrame::new([0.2, 1.0, 0.0], [1.0, -0.2, 0.0]).unwrap();
    assert!(
        PreparedRigidHertzMindlinIsland::new(
            &bodies,
            &properties,
            &coupled,
            DT,
            controls(),
            &mut || false
        )
        .is_err()
    );
    // Central normal arms preserve zero N/T coupling even when a shear impulse
    // excites excluded torsion/lateral velocity through a full inertia tensor.
    let mut inertias = properties;
    inertias[0] = MassProperties::from_principal_frame(
        MASS[0],
        Vec3::ZERO,
        Vec3::new(0.8 * IZ[0], 0.9 * IZ[0], IZ[0]),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.4).unwrap(),
    )
    .unwrap();
    let prepared = PreparedRigidHertzMindlinIsland::new_driven(
        &bodies,
        &inertias,
        &contacts,
        &drive,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        prepared.advance(&mut || false),
        Err(RigidHertzMindlinIslandError::Mechanics(
            RigidMindlinIslandError::UnsupportedMotion { .. }
        ))
    ));
    let mut rates = bodies;
    rates[0] =
        RigidBodyState::new(bodies[0].pose(), Vec3::new(0.0, 0.01, 0.0), Vec3::ZERO).unwrap();
    assert!(
        PreparedRigidHertzMindlinIsland::new(
            &rates,
            &properties,
            &contacts,
            DT,
            controls(),
            &mut || false
        )
        .is_err()
    );
    for kind in 0..3 {
        let mut budget = controls();
        match kind {
            0 => budget.maximum_normal_evaluations = 1,
            1 => budget.maximum_normal_iterations = 1,
            _ => budget.native.multiple.maximum_law_evaluations = 1,
        }
        let prepared = PreparedRigidHertzMindlinIsland::new_driven(
            &bodies,
            &properties,
            &contacts,
            &drive,
            DT,
            budget,
            &mut || false,
        )
        .unwrap();
        assert!(prepared.advance(&mut || false).is_err());
    }
    // Individual original radii begin disjoint; doubling both indentations
    // grows them into overlapping footprints and must not publish either state.
    let mut close_sites = contacts;
    let mut close_bodies = bodies;
    for i in 0..2 {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        close_sites[i].arms_body[0] = Vec3::new(0.0, sign * 0.0015, 0.0);
        close_bodies[i + 1] = RigidBodyState::new(
            Pose::new(
                Vec3::new(0.0, sign * 0.0015 - ARMS[i][1], 0.0),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap();
    }
    let end = [4e-5; 2];
    let expected = [0, 1].map(|i| {
        laws[i]
            .advance_coupled_path(&states[i], &input(end[i], 0.0), &mut || false)
            .unwrap()
    });
    let (grow, _) = reference_drive(
        &close_bodies,
        &properties,
        &close_sites,
        end,
        [0.0; 2],
        &expected,
    );
    let prepared = PreparedRigidHertzMindlinIsland::new_driven(
        &close_bodies,
        &properties,
        &close_sites,
        &grow,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert!(prepared.advance(&mut || false).is_err());
    let mut separating = drive;
    separating[1].force_world.y = 20.0;
    let prepared = PreparedRigidHertzMindlinIsland::new_driven(
        &bodies,
        &properties,
        &contacts,
        &separating,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap();
    assert!(prepared.advance(&mut || false).is_err());
    let prepared = PreparedRigidHertzMindlinIsland::new_driven(
        &bodies,
        &properties,
        &contacts,
        &drive,
        DT,
        controls(),
        &mut || false,
    )
    .unwrap();
    let mut polls = 0;
    let accepted = prepared
        .advance(&mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut last = 0;
    assert!(matches!(
        prepared.advance(&mut || {
            last += 1;
            last == polls
        }),
        Err(RigidHertzMindlinIslandError::Cancelled)
    ));
    assert_eq!(last, polls);
    assert_eq!(prepared.advance(&mut || false).unwrap(), accepted);
    assert_eq!(states, saved);
    assert_eq!(bodies, original);
}
