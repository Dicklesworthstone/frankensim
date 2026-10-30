//! Independent finite-body checks for simultaneous surface and material bulk loss.
//! Synthetic material data test the declared model, not experimental validation.

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::RigidHertzControls;
use fs_mbd::tangential_friction::hertz_rolling::multiple::*;
use fs_mbd::tangential_friction::hertz_rolling::*;
use fs_mbd::tangential_friction::rolling::RigidRollingHistory;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::rolling_loss::viscoelastic::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::sliding_patch::contact::*;
use fs_tribo::{
    ApplicabilityRange, HertzSpherePlane, InputAuthority, InterfaceMedium, InterfaceSystemRef,
};

const R: f64 = 0.1;
const L: f64 = 0.06;
const DT: f64 = 0.0002;

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e} tol {tol:e}");
}
fn close_vec(a: Vec3, b: Vec3, tol: f64) {
    close(a.x, b.x, tol);
    close(a.y, b.y, tol);
    close(a.z, b.z, tol);
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}

fn law(tau: f64) -> RollingLossLaw {
    RollingLossLaw::ViscoelasticSphere(
        ViscoelasticSphereRollingCard::new(
            VISCOELASTIC_SPHERE_LAW_ID,
            "synthetic shared sphere",
            InputAuthority::SyntheticFixture,
            ViscoelasticSphereParameters {
                sphere_surface: DeformingSphereSurface::SurfaceA,
                radius_m: R,
                young_modulus_pa: 1e6,
                poisson_ratio: 0.25,
                density_kg_m3: 1000.0,
                shear_viscosity_pa_s: tau * 1e6 / 2.5,
                bulk_viscosity_pa_s: tau * 1e6 / 1.5,
            },
            QuasistaticSphereLimits {
                maximum_indentation_radius_ratio: 0.1,
                maximum_speed_wave_ratio: 0.1,
                maximum_relaxation_ratio: 0.1,
            },
            RollingLossApplicability::new(
                ApplicabilityRange::new(280.0, 320.0).unwrap(),
                ApplicabilityRange::new(0.0, 100.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap(),
    )
}
fn controls() -> RigidHertzRollingIslandControls {
    RigidHertzRollingIslandControls {
        mechanical: RigidHertzRollingControls {
            mechanical: RigidHertzControls {
                maximum_bodies: 3,
                geometry_tolerance_m: 2e-12,
                velocity_tolerance_m_s: 3e-10,
                angular_velocity_tolerance_rad_s: 3e-9,
                linear_momentum_tolerance_ns: 2e-10,
                angular_momentum_tolerance_nms: 2e-10,
                energy_tolerance_j: 3e-11,
                contact: HertzWrenchConfig {
                    maximum_static_iterations: 96,
                    maximum_iterations: 256,
                    maximum_evaluations: 1024,
                    force_tolerance_n: 3e-9,
                    moment_tolerance_nm: 3e-10,
                    displacement_tolerance_m: 3e-13,
                    energy_tolerance_j: 3e-12,
                    relative_matrix_tolerance: 1e-13,
                    maximum_force_n: 20.0,
                    maximum_spin_moment_nm: 2.0,
                },
            },
            maximum_rolling_moment_nm: 1.0,
            rolling_moment_tolerance_nm: 3e-10,
        },
        maximum_contacts: 4,
        maximum_matrix_terms: 400,
        maximum_sweeps: 512,
        maximum_evaluations: 65536,
    }
}

struct Fixture {
    bodies: [RigidBodyState; 3],
    properties: [MassProperties; 3],
    bindings: [RigidContactBinding; 2],
    frames: [TangentFrame; 2],
    hertz: [HertzWrenchContact; 2],
    rolling: RollingLossLaw,
    patches: [RollingPatchReceipt; 2],
    interface: InterfaceSystemRef,
    work: [RollingWorkOwnership; 2],
    histories: [RigidRollingHistory; 2],
}
impl Fixture {
    fn contact(&self, i: usize) -> RigidHertzRollingContact<'_> {
        RigidHertzRollingContact {
            binding: self.bindings[i],
            frame: self.frames[i],
            hertz: &self.hertz[i],
            patch: &self.patches[i],
            interface: &self.interface,
            rolling: &self.rolling,
            temperature_k: 300.0,
            excitation_frequency_hz: 10.0,
            work: &self.work[i],
        }
    }
    fn contacts(&self) -> [RigidHertzRollingIslandContact<'_>; 2] {
        std::array::from_fn(|i| RigidHertzRollingIslandContact {
            patch_id: self.patches[i].patch_id(),
            contact: self.contact(i),
            history: &self.histories[i],
        })
    }
    fn laws(&self) -> [RigidHertzRollingIslandLaw<'_>; 2] {
        std::array::from_fn(|i| RigidHertzRollingIslandLaw {
            patch_id: self.patches[i].patch_id(),
            law: &self.hertz[i],
            temperature_k: 300.0,
        })
    }
    fn prepare(&self) -> PreparedRigidHertzRollingIsland<'_> {
        PreparedRigidHertzRollingIsland::new(
            &self.bodies,
            &self.properties,
            &self.contacts(),
            DT,
            controls(),
            &CancelGate::new(),
        )
        .unwrap()
    }
}
fn fixture(tau: f64, mu: [f64; 2], velocity: Vec3, omega: Vec3) -> Fixture {
    let interface = InterfaceSystemRef::new(
        "sphere/two supports",
        "independent material rolling and surface slip",
        "synthetic interface",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    let normals = [0.001, 0.0008].map(|delta| {
        HertzSpherePlane {
            effective_radius: R,
            reduced_modulus: 1e6 / (1.0 - 0.25_f64.powi(2)),
        }
        .response(&interface, delta)
        .unwrap()
    });
    let patches = std::array::from_fn(|i| {
        RollingPatchReceipt::new(
            if i == 0 { "lower patch" } else { "upper patch" },
            "original Hertz",
            "synthetic normal owner",
            InputAuthority::SyntheticFixture,
            normals[i].normal_force_n,
            PI * normals[i].contact_radius_m.powi(2),
            PatchCurvature::Principal {
                first_per_m: 1.0 / R,
                second_per_m: 1.0 / R,
            },
        )
        .unwrap()
    });
    let mass = 4.0 * PI * 1000.0 * R.powi(3) / 3.0;
    let inertia = 0.4 * mass * R * R;
    let properties = [
        MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap(),
        MassProperties::new(7.0, Vec3::ZERO, Vec3::new(0.08, 0.13, 0.18)).unwrap(),
        MassProperties::new(11.0, Vec3::ZERO, Vec3::new(0.17, 0.21, 0.25)).unwrap(),
    ];
    let bodies = std::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                Vec3::new(
                    0.0,
                    0.0,
                    if i == 0 {
                        0.0
                    } else if i == 1 {
                        -R - L
                    } else {
                        R + L
                    },
                ),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            if i == 0 {
                velocity.scale(mass)
            } else {
                Vec3::ZERO
            },
            if i == 0 {
                properties[0].angular_momentum_body_checked(omega).unwrap()
            } else {
                Vec3::ZERO
            },
        )
        .unwrap()
    });
    Fixture {
        bodies,
        properties,
        bindings: std::array::from_fn(|i| RigidContactBinding {
            body_a: 0,
            body_b: i + 1,
            arm_a_body: Vec3::new(0.0, 0.0, if i == 0 { -R } else { R }),
            arm_b_body: Vec3::new(0.0, 0.0, if i == 0 { L } else { -L }),
        }),
        frames: [
            TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap(),
        ],
        hertz: std::array::from_fn(|i| {
            HertzWrenchContact::new(HertzWrenchParameters {
                normal_force_n: normals[i].normal_force_n,
                contact_radius_m: normals[i].contact_radius_m,
                static_mu: mu[i],
                kinetic_mu: mu[i],
            })
            .unwrap()
        }),
        rolling: law(tau),
        patches,
        interface,
        work: std::array::from_fn(|i| {
            RollingWorkOwnership::new(
                if i == 0 { "lower patch" } else { "upper patch" },
                "same original interval",
                "relative tangent rolling vector",
                RollingLossChannel::RollingHysteresis,
            )
            .unwrap()
        }),
        histories: [RigidRollingHistory::new(), RigidRollingHistory::new()],
    }
}

fn audit(f: &Fixture, result: &RigidHertzRollingIslandCandidate) {
    assert_eq!(result.body_indices, vec![0, 1, 2]);
    let before = f
        .bodies
        .iter()
        .zip(&f.properties)
        .map(|(b, p)| b.kinetic_energy(*p).unwrap())
        .sum::<f64>();
    let after = result
        .bodies_after
        .iter()
        .zip(&f.properties)
        .map(|(b, p)| b.kinetic_energy(*p).unwrap())
        .sum::<f64>();
    close(before, result.report.kinetic_before_j, 2e-15);
    close(after, result.report.kinetic_after_j, 2e-15);
    close(
        after - before + result.report.surface_heat_j + result.report.bulk_heat_j,
        0.0,
        3e-11,
    );
    close_vec(result.report.linear_momentum_change_ns, Vec3::ZERO, 2e-12);
    close_vec(result.report.angular_momentum_change_nms, Vec3::ZERO, 2e-12);
    close(
        result.report.contacts.iter().map(|r| r.bulk_heat_j).sum(),
        result.report.bulk_heat_j,
        1e-20,
    );
    close(
        result
            .report
            .contacts
            .iter()
            .map(|r| r.surface_heat_j)
            .sum(),
        result.report.surface_heat_j,
        1e-20,
    );
    for (r, h) in result.report.contacts.iter().zip(&result.rolling_histories) {
        assert_eq!(r.patch_id, h.patch_id);
        close(h.history.state().cumulative_heat_j(), r.bulk_heat_j, 1e-20);
        close(r.rolling_work_residual_j, 0.0, 3e-12);
        close(r.surface_work_residual_j, 0.0, 3e-12);
    }
}

#[test]
fn two_original_patches_share_one_sphere_with_independent_dense_rolling_solution() {
    let omega = 0.025;
    let f = fixture(0.01, [0.0, 0.0], Vec3::ZERO, Vec3::new(0.0, omega, 0.0));
    let prepared = f.prepare();
    let result = prepared.advance(&f.laws(), &CancelGate::new()).unwrap();
    audit(&f, &result);
    assert_eq!(result.report.surface_heat_j, 0.0);
    assert!(result.report.bulk_heat_j > 0.0);
    let is = f.properties[0]
        .angular_momentum_body_checked(Vec3::new(0.0, 1.0, 0.0))
        .unwrap()
        .y;
    let ib = [0.13, 0.21];
    let c = f.patches.each_ref().map(|p| 0.01 * p.normal_load_n() * R);
    // Direct physical torque equations, with no reduced Hertz matrix or
    // production mobility/factor used in this reference.
    let a = 1.0 + 0.5 * DT * c[0] * (1.0 / is + 1.0 / ib[0]);
    let b = 0.5 * DT * c[0] / is;
    let cc = 0.5 * DT * c[1] / is;
    let d = 1.0 + 0.5 * DT * c[1] * (1.0 / is + 1.0 / ib[1]);
    let determinant = a * d - b * cc;
    let torques = [
        (-c[0] * omega * d + b * c[1] * omega) / determinant,
        (-a * c[1] * omega + cc * c[0] * omega) / determinant,
    ];
    close(
        result.bodies_after[0].angular_momentum_body().y,
        f.bodies[0].angular_momentum_body().y + DT * (torques[0] + torques[1]),
        2e-16,
    );
    for i in 0..2 {
        close(
            result.bodies_after[i + 1].angular_momentum_body().y,
            -DT * torques[i],
            2e-16,
        );
        let midpoint = -torques[i] / c[i];
        close(
            result.report.contacts[i].bulk_heat_j,
            DT * c[i] * midpoint.powi(2),
            2e-18,
        );
        close(
            result.report.contacts[i]
                .rolling_angular_impulse_world_nms
                .y,
            DT * torques[i],
            2e-16,
        );
    }
    let w = prepared.scaled_mobility_per_kg();
    assert!(
        w[4 * 10 + 9] != 0.0,
        "shared sphere must couple both rolling channels"
    );

    // One enclosing support can own the two opposite planar faces. Its center
    // of mass coincides with the sphere's; its geometry is not a second sphere.
    // This gives two distinct original pads but only five independent wrench
    // coordinates, so the full ten-coordinate response is dependent.
    let mut same_pair = fixture(0.01, [0.0, 0.0], Vec3::ZERO, Vec3::new(0.0, omega, 0.0));
    same_pair.bodies[1] = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap();
    for i in 0..2 {
        same_pair.bindings[i].body_b = 1;
        same_pair.bindings[i].arm_b_body = same_pair.bindings[i].arm_a_body;
    }
    let dependent = same_pair
        .prepare()
        .advance(&same_pair.laws(), &CancelGate::new())
        .unwrap();
    assert_eq!(dependent.body_indices, vec![0, 1]);
    let midpoint = omega / (1.0 + 0.5 * DT * (c[0] + c[1]) * (1.0 / is + 1.0 / ib[0]));
    let total_torque = -(c[0] + c[1]) * midpoint;
    close(
        dependent.bodies_after[0].angular_momentum_body().y,
        same_pair.bodies[0].angular_momentum_body().y + DT * total_torque,
        2e-16,
    );
    close(
        dependent.bodies_after[1].angular_momentum_body().y,
        -DT * total_torque,
        2e-16,
    );
    close(
        dependent.report.bulk_heat_j,
        DT * (c[0] + c[1]) * midpoint.powi(2),
        2e-18,
    );
    close(dependent.report.total_energy_residual_j, 0.0, 3e-13);
}

#[test]
fn slip_induces_rolling_loss_from_zero_initial_rolling_and_one_shared_recoil() {
    let f = fixture(0.01, [0.12, 0.2], Vec3::new(0.012, 0.003, 0.0), Vec3::ZERO);
    let result = f.prepare().advance(&f.laws(), &CancelGate::new()).unwrap();
    audit(&f, &result);
    assert!(result.report.surface_heat_j > 0.0);
    assert!(result.report.bulk_heat_j > 0.0);
    assert!(
        result
            .report
            .contacts
            .iter()
            .all(|r| r.scaled_rates_before_m_s[3..] == [0.0, 0.0])
    );
    assert!(
        result
            .report
            .contacts
            .iter()
            .any(|r| r.rolling_angular_impulse_world_nms != Vec3::ZERO)
    );
    let zero = fixture(0.0, [0.12, 0.2], Vec3::new(0.012, 0.003, 0.0), Vec3::ZERO);
    let unresisted = zero
        .prepare()
        .advance(&zero.laws(), &CancelGate::new())
        .unwrap();
    assert_eq!(unresisted.report.bulk_heat_j, 0.0);
    assert!(
        norm(
            result.bodies_after[0]
                .angular_momentum_body()
                .sub(unresisted.bodies_after[0].angular_momentum_body())
        ) > 1e-14
    );
}

#[test]
fn canonical_patch_order_coefficient_trials_and_untouched_spectators() {
    let f = fixture(
        0.01,
        [0.15, 0.2],
        Vec3::new(0.012, -0.003, 0.0),
        Vec3::new(0.015, 0.025, 0.1),
    );
    let prepared = f.prepare();
    let reference = prepared.advance(&f.laws(), &CancelGate::new()).unwrap();
    audit(&f, &reference);
    let mut contacts = f.contacts();
    contacts.reverse();
    let mut laws = f.laws();
    laws.reverse();
    let mut bodies = f.bodies.to_vec();
    let mut properties = f.properties.to_vec();
    bodies.push(
        RigidBodyState::new(
            Pose::new(Vec3::new(8.0, 3.0, -5.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(100.0, 200.0, 300.0),
            Vec3::ZERO,
        )
        .unwrap(),
    );
    properties.push(MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap());
    let reordered = PreparedRigidHertzRollingIsland::new(
        &bodies,
        &properties,
        &contacts,
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&laws, &CancelGate::new())
    .unwrap();
    assert_eq!(reference, reordered);
    let before = bodies.clone();
    let histories = f.histories.clone();
    let frictionless = f.hertz.each_ref().map(|h| {
        HertzWrenchContact::new(HertzWrenchParameters {
            static_mu: 0.0,
            kinetic_mu: 0.0,
            ..h.parameters()
        })
        .unwrap()
    });
    let other = std::array::from_fn::<_, 2, _>(|i| RigidHertzRollingIslandLaw {
        law: &frictionless[i],
        ..f.laws()[i]
    });
    let trial = prepared.advance(&other, &CancelGate::new()).unwrap();
    assert_eq!(trial.report.surface_heat_j, 0.0);
    assert_eq!(
        reference,
        prepared.advance(&f.laws(), &CancelGate::new()).unwrap()
    );
    assert_eq!(bodies, before);
    assert_eq!(f.histories, histories);
}

#[test]
fn singleton_matches_original_joint_owner_and_released_zero_ports() {
    let f = fixture(
        0.01,
        [0.15, 0.2],
        Vec3::new(0.012, 0.003, 0.0),
        Vec3::new(0.015, 0.025, 0.1),
    );
    let old = PreparedRigidHertzRollingPair::new(
        &f.bodies,
        &f.properties,
        &f.histories[0],
        f.contact(0),
        DT,
        controls().mechanical,
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.hertz[0], &CancelGate::new())
    .unwrap();
    let island = PreparedRigidHertzRollingIsland::new(
        &f.bodies,
        &f.properties,
        &f.contacts()[..1],
        DT,
        controls(),
        &CancelGate::new(),
    )
    .unwrap()
    .advance(&f.laws()[..1], &CancelGate::new())
    .unwrap();
    for i in 0..2 {
        close_vec(
            old.bodies_after[i].linear_momentum_world(),
            island.bodies_after[i].linear_momentum_world(),
            3e-12,
        );
        close_vec(
            old.bodies_after[i].angular_momentum_body(),
            island.bodies_after[i].angular_momentum_body(),
            3e-12,
        );
    }
    close(
        old.report.surface_heat_j,
        island.report.surface_heat_j,
        3e-12,
    );
    close(old.report.bulk_heat_j, island.report.bulk_heat_j, 3e-12);
    let mut released = f;
    for i in 0..2 {
        released.patches[i] = RollingPatchReceipt::new(
            if i == 0 { "lower patch" } else { "upper patch" },
            "original Hertz",
            "synthetic normal owner",
            InputAuthority::SyntheticFixture,
            0.0,
            0.0,
            PatchCurvature::Principal {
                first_per_m: 1.0 / R,
                second_per_m: 1.0 / R,
            },
        )
        .unwrap();
        released.hertz[i] = HertzWrenchContact::new(HertzWrenchParameters {
            normal_force_n: 0.0,
            ..released.hertz[i].parameters()
        })
        .unwrap();
    }
    let result = released
        .prepare()
        .advance(&released.laws(), &CancelGate::new())
        .unwrap();
    assert_eq!(result.bodies_after, released.bodies);
    assert_eq!(result.report.surface_heat_j, 0.0);
    assert_eq!(result.report.bulk_heat_j, 0.0);
}

#[test]
fn contradictory_materials_overlapping_patches_and_split_bulk_temperatures_refuse() {
    let f = fixture(
        0.01,
        [0.1, 0.2],
        Vec3::new(0.012, 0.0, 0.0),
        Vec3::new(0.0, 0.025, 0.0),
    );
    let other = law(0.02);
    let mut contacts = f.contacts();
    contacts[1].contact.rolling = &other;
    assert!(
        PreparedRigidHertzRollingIsland::new(
            &f.bodies,
            &f.properties,
            &contacts,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    contacts = f.contacts();
    contacts[1].contact.temperature_k = 301.0;
    assert!(
        PreparedRigidHertzRollingIsland::new(
            &f.bodies,
            &f.properties,
            &contacts,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    contacts = f.contacts();
    contacts[1].patch_id = "aliased original patch";
    assert!(
        PreparedRigidHertzRollingIsland::new(
            &f.bodies,
            &f.properties,
            &contacts,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    contacts = f.contacts();
    contacts[1].contact.binding = f.bindings[0];
    contacts[1].contact.frame = f.frames[0];
    assert!(
        PreparedRigidHertzRollingIsland::new(
            &f.bodies,
            &f.properties,
            &contacts,
            DT,
            controls(),
            &CancelGate::new()
        )
        .is_err()
    );
    let prepared = f.prepare();
    let reference = prepared.advance(&f.laws(), &CancelGate::new()).unwrap();
    let mut laws = f.laws();
    laws[1].temperature_k = 301.0;
    assert!(prepared.advance(&laws, &CancelGate::new()).is_err());
    laws = f.laws();
    laws[0].temperature_k = 500.0;
    laws[1].temperature_k = 500.0;
    assert!(prepared.advance(&laws, &CancelGate::new()).is_err());
    assert_eq!(
        reference,
        prepared.advance(&f.laws(), &CancelGate::new()).unwrap()
    );
}

#[test]
fn work_budgets_cancellation_and_late_physical_refusal_leave_inputs_unchanged() {
    let f = fixture(
        0.01,
        [0.1, 0.2],
        Vec3::new(0.012, 0.0, 0.0),
        Vec3::new(0.0, 0.025, 0.0),
    );
    let before = f.bodies;
    let histories = f.histories.clone();
    let mut c = controls();
    c.maximum_matrix_terms = 99;
    assert!(
        PreparedRigidHertzRollingIsland::new(
            &f.bodies,
            &f.properties,
            &f.contacts(),
            DT,
            c,
            &CancelGate::new()
        )
        .is_err()
    );
    c = controls();
    c.mechanical.maximum_rolling_moment_nm = 0.0;
    let limited = PreparedRigidHertzRollingIsland::new(
        &f.bodies,
        &f.properties,
        &f.contacts(),
        DT,
        c,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(matches!(
        limited.advance(&f.laws(), &CancelGate::new()),
        Err(RigidHertzRollingIslandError::RollingMomentBudget { .. })
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(
        PreparedRigidHertzRollingIsland::new(
            &f.bodies,
            &f.properties,
            &f.contacts(),
            DT,
            controls(),
            &gate
        )
        .is_err()
    );
    let prepared = f.prepare();
    let reference = prepared.advance(&f.laws(), &CancelGate::new()).unwrap();
    assert!(prepared.advance(&f.laws(), &gate).is_err());
    assert_eq!(
        reference,
        prepared.advance(&f.laws(), &CancelGate::new()).unwrap()
    );
    let coarse = PreparedRigidHertzRollingIsland::new(
        &f.bodies,
        &f.properties,
        &f.contacts(),
        0.1,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(coarse.advance(&f.laws(), &CancelGate::new()).is_err());
    assert_eq!(f.bodies, before);
    assert_eq!(f.histories, histories);
}

fn inverse(m: [f64; 9]) -> [f64; 9] {
    let [a, b, c, d, e, f, g, h, i] = m;
    let determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    [
        (e * i - f * h) / determinant,
        (c * h - b * i) / determinant,
        (b * f - c * e) / determinant,
        (f * g - d * i) / determinant,
        (a * i - c * g) / determinant,
        (c * d - a * f) / determinant,
        (d * h - e * g) / determinant,
        (b * g - a * h) / determinant,
        (a * e - b * d) / determinant,
    ]
}
fn multiply(m: [f64; 9], v: Vec3) -> Vec3 {
    Vec3::new(
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[3] * v.x + m[4] * v.y + m[5] * v.z,
        m[6] * v.x + m[7] * v.y + m[8] * v.z,
    )
}
fn vector(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn array(x: Vec3) -> [f64; 3] {
    [x.x, x.y, x.z]
}

#[test]
fn full_tensor_shared_response_matches_adjugate_inverse_and_world_frame_covariance() {
    let mut f = fixture(
        0.01,
        [0.15, 0.2],
        Vec3::new(0.012, -0.003, 0.0),
        Vec3::new(0.015, 0.025, 0.1),
    );
    for i in 1..3 {
        let old = f.properties[i];
        f.properties[i] = MassProperties::from_principal_frame(
            old.mass(),
            Vec3::ZERO,
            old.principal_inertia_body(),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.4 + 0.2 * i as f64)
                .unwrap(),
        )
        .unwrap();
    }
    let mut translation = vec![[Vec3::ZERO; 3]; 10];
    let mut rotation = vec![[Vec3::ZERO; 3]; 10];
    for i in 0..2 {
        let b = f.bindings[i];
        let axes = [
            vector(f.frames[i].longitudinal()),
            vector(f.frames[i].lateral()),
            vector(f.frames[i].normal()),
        ];
        for (body, arm, sign) in [
            (b.body_a, b.arm_a_body, 1.0),
            (b.body_b, b.arm_b_body, -1.0),
        ] {
            for axis in 0..5 {
                let direction = axes[if axis < 3 { axis } else { axis - 3 }];
                if axis < 2 {
                    translation[5 * i + axis][body] = direction.scale(sign);
                    rotation[5 * i + axis][body] = arm.cross(direction).scale(sign);
                } else {
                    rotation[5 * i + axis][body] = direction.scale(
                        sign * if axis == 2 {
                            f.hertz[i].parameters().contact_radius_m
                        } else {
                            R
                        },
                    );
                }
            }
        }
    }
    let prepared = f.prepare();
    for i in 0..10 {
        for j in 0..10 {
            let expected = (0..3)
                .map(|b| {
                    translation[i][b].dot(translation[j][b]) / f.properties[b].mass()
                        + rotation[i][b].dot(multiply(
                            inverse(f.properties[b].inertia_tensor_body().m),
                            rotation[j][b],
                        ))
                })
                .sum::<f64>();
            close(
                prepared.scaled_mobility_per_kg()[10 * i + j],
                expected,
                3e-14,
            );
        }
    }
    let reference = prepared.advance(&f.laws(), &CancelGate::new()).unwrap();
    audit(&f, &reference);
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.63).unwrap();
    let shift = Vec3::new(0.7, -0.3, 0.2);
    let boost = Vec3::new(0.002, -0.001, 0.0005);
    for i in 0..3 {
        let b = f.bodies[i];
        f.bodies[i] = RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(b.pose().position_world()).add(shift),
                q,
            )
            .unwrap(),
            q.rotate_body_to_world(b.linear_momentum_world())
                .add(boost.scale(f.properties[i].mass())),
            b.angular_momentum_body(),
        )
        .unwrap();
    }
    for i in 0..2 {
        let angle = [0.31_f64, -0.72][i];
        let n = Vec3::new(0.0, 0.0, if i == 0 { 1.0 } else { -1.0 });
        f.frames[i] = TangentFrame::new(
            array(q.rotate_body_to_world(n)),
            array(q.rotate_body_to_world(Vec3::new(angle.cos(), angle.sin(), 0.0))),
        )
        .unwrap();
    }
    let transformed = f.prepare().advance(&f.laws(), &CancelGate::new()).unwrap();
    audit(&f, &transformed);
    for i in 0..3 {
        close_vec(
            transformed.bodies_after[i].linear_momentum_world(),
            q.rotate_body_to_world(reference.bodies_after[i].linear_momentum_world())
                .add(boost.scale(f.properties[i].mass())),
            3e-11,
        );
        close_vec(
            transformed.bodies_after[i].angular_momentum_body(),
            reference.bodies_after[i].angular_momentum_body(),
            3e-12,
        );
    }
    close(
        transformed.report.surface_heat_j,
        reference.report.surface_heat_j,
        3e-13,
    );
    close(
        transformed.report.bulk_heat_j,
        reference.report.bulk_heat_j,
        3e-13,
    );
    for i in 0..2 {
        close_vec(
            transformed.report.contacts[i].force_impulse_world_ns,
            q.rotate_body_to_world(reference.report.contacts[i].force_impulse_world_ns),
            3e-11,
        );
        close_vec(
            transformed.report.contacts[i].rolling_angular_impulse_world_nms,
            q.rotate_body_to_world(reference.report.contacts[i].rolling_angular_impulse_world_nms),
            3e-12,
        );
    }
}
