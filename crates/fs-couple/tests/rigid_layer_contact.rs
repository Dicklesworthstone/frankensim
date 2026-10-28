//! G0/G1/G3/G4: native finite-thickness shear and its original finite heat store.
//! Synthetic material values exercise the declared linear director interpolation;
//! none is a calibrated layer specimen or a complete bending/inertial solid.
use fs_couple::friction_heat::third_body::rigid::*;
use fs_couple::friction_heat::third_body::{ThirdBodyHeatRoute, ThirdBodyHeatTolerances};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::{TangentialLaw, TangentialSolveConfig};
use std::cell::Cell;
use std::collections::BTreeMap;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance {tolerance}"
    );
}
fn vector_close(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn geometry() -> RigidLayerGeometry {
    RigidLayerGeometry {
        arms_body_m: [Vec3::new(0.5, 0.0, 0.0), Vec3::new(-0.5, 0.0, 0.0)],
        frame: TangentFrame::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap(),
        thickness_m: 1.0,
        area_m2: 1.0,
    }
}
fn properties(anisotropic: bool) -> [MassProperties; 2] {
    [
        MassProperties::new(
            2.0,
            Vec3::ZERO,
            if anisotropic {
                Vec3::new(0.4, 0.6, 0.8)
            } else {
                Vec3::new(1.0, 1.0, 1.0)
            },
        )
        .unwrap(),
        MassProperties::new(
            3.0,
            Vec3::ZERO,
            if anisotropic {
                Vec3::new(0.5, 0.7, 0.9)
            } else {
                Vec3::new(1.0, 1.0, 1.0)
            },
        )
        .unwrap(),
    ]
}
fn bodies(anisotropic: bool, velocity_b: Vec3) -> [RigidBodyState; 2] {
    [
        RigidBodyState::new(
            Pose::new(
                Vec3::new(-0.5, 0.0, 0.0),
                if anisotropic {
                    UnitQuaternion::from_axis_angle(
                        Vec3::new(1.0, 0.0, 0.0),
                        core::f64::consts::FRAC_PI_4,
                    )
                    .unwrap()
                } else {
                    UnitQuaternion::IDENTITY
                },
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap(),
        RigidBodyState::new(
            Pose::new(Vec3::new(1.5, 0.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
            velocity_b.scale(3.0),
            Vec3::ZERO,
        )
        .unwrap(),
    ]
}
fn law(k: f64, yield_n: f64, resistance: f64, temperature: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(k, yield_n, resistance, temperature, 1.0).unwrap()
}
fn law_solve() -> TangentialSolveConfig {
    TangentialSolveConfig {
        max_iterations: 32,
        maximum_force_n: 1000.0,
        force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-12,
    }
}
fn controls() -> RigidLayerSolve {
    RigidLayerSolve {
        maximum_trials: 1000,
        maximum_newton_iterations: 40,
        difference_step_n: 1e-4,
        geometry_tolerance_m: 1e-12,
        force_tolerance_n: 2e-9,
        linear_momentum_tolerance_ns: 2e-10,
        angular_momentum_tolerance_nms: 2e-10,
        energy_tolerance_j: 2e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..Default::default()
        },
    }
}

// Independent two-coordinate algebra: solve (A + beta I)F=t, where beta is
// the scalar Bingham flow factor. This is not the production native/Newton path.
fn reference_force(w: [[f64; 2]; 2], free: [f64; 2], k: f64, y: f64, c: f64, dt: f64) -> [f64; 2] {
    let a = [
        [
            1.0 + k * dt * dt * w[0][0] / 4.0,
            k * dt * dt * w[0][1] / 4.0,
        ],
        [
            k * dt * dt * w[1][0] / 4.0,
            1.0 + k * dt * dt * w[1][1] / 4.0,
        ],
    ];
    let t = free.map(|u| 0.5 * k * dt * u);
    let solve = |beta: f64| {
        let det = (a[0][0] + beta) * (a[1][1] + beta) - a[0][1] * a[1][0];
        [
            ((a[1][1] + beta) * t[0] - a[0][1] * t[1]) / det,
            ((a[0][0] + beta) * t[1] - a[1][0] * t[0]) / det,
        ]
    };
    let trial = solve(0.0);
    if trial[0].hypot(trial[1]) <= y {
        return trial;
    }
    let b = k * dt / (2.0 * c);
    let (mut lower, mut upper) = (0.0, b);
    for _ in 0..100 {
        let mid = 0.5 * (lower + upper);
        let f = solve(mid);
        let residual = mid - b * (1.0 - y / f[0].hypot(f[1])).max(0.0);
        if residual > 0.0 {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    solve(0.5 * (lower + upper))
}

#[test]
fn g1_actual_separated_face_wrenches_resolve_anisotropic_elastic_and_plastic_shear() {
    let p = properties(true);
    let old = bodies(true, Vec3::new(0.0, 1.0, 0.4));
    let dt = 0.05;
    // Face-to-COM arms are +/-0.5 x; the thickness is 1. The derived force
    // plus face-couple torque arm is 1, not the bare point arm 0.5.
    let avg = 0.5 * (1.0 / 0.6 + 1.0 / 0.8);
    let off = -0.5 * (1.0 / 0.6 - 1.0 / 0.8);
    let w = [
        [5.0 / 6.0 + avg + 1.0 / 0.9, off],
        [off, 5.0 / 6.0 + avg + 1.0 / 0.7],
    ];
    for y in [100.0, 0.2] {
        let model = law(100.0, y, 2.0, 300.0);
        let history = RigidLayerHistory::new(&model, geometry(), 0.0).unwrap();
        let candidate = step_rigid_layer_candidate(
            &old,
            p,
            &history,
            &model,
            dt,
            10.0,
            law_solve(),
            controls(),
            [0.0; 2],
            || false,
        )
        .unwrap();
        let report = candidate.report();
        let expected = reference_force(w, [1.0, 0.4], 100.0, y, 2.0, dt);
        close(report.applied_force_on_a_n.y, expected[0], 2e-9);
        close(report.applied_force_on_a_n.z, expected[1], 2e-9);
        for i in 0..2 {
            for j in 0..2 {
                close(report.inverse_shear_mass[i][j], w[i][j], 2e-14);
            }
        }
        assert!(report.inverse_shear_mass[0][1].abs() > 0.2);
        let f = Vec3::new(0.0, expected[0], expected[1]);
        let couple = Vec3::new(0.0, -0.5 * expected[1], 0.5 * expected[0]);
        for side in 0..2 {
            vector_close(report.face_couple_nm[side], couple, 2e-9);
            let after = candidate.bodies()[side];
            vector_close(
                after
                    .linear_momentum_world()
                    .sub(old[side].linear_momentum_world()),
                f.scale(if side == 0 { dt } else { -dt }),
                2e-10,
            );
            let spin = after
                .pose()
                .orientation()
                .rotate_body_to_world(after.angular_momentum_body());
            vector_close(spin, couple.scale(2.0 * dt), 2e-10);
        }
        vector_close(report.angular_momentum_change_nms, Vec3::ZERO, 2e-13);
        vector_close(report.linear_momentum_change_ns, Vec3::ZERO, 2e-13);
        close(report.total_energy_residual_j, 0.0, 2e-10);
        let r = expected[0].hypot(expected[1]);
        let plastic = dt / 2.0 * (r - y).max(0.0);
        let exact_heat = y * plastic + 2.0 / dt * plastic * plastic;
        close(report.constitutive.dissipated_j, exact_heat, 2e-10);
        for i in 0..2 {
            close(
                candidate.history().deformation().elastic_slip_m()[i],
                2.0 * expected[i] / 100.0,
                2e-10,
            );
        }
        if y > 1.0 {
            assert_eq!(exact_heat, 0.0);
        } else {
            assert!(exact_heat > 0.01);
        }
        // Removing the physical face couples would leave exactly this spurious
        // orbital impulse moment, which is materially nonzero in this fixture.
        assert!(couple.norm_squared().sqrt() * 2.0 * dt > 0.01);
        assert_eq!(old, bodies(true, Vec3::new(0.0, 1.0, 0.4)));
        assert_eq!(history.time_s(), 0.0);
    }
}

#[test]
fn g3_common_rigid_rotation_has_no_shear_and_rotated_frames_preserve_native_response() {
    let p = properties(false);
    let omega = Vec3::new(0.2, -0.3, 0.7);
    let translate = Vec3::new(2.0, -1.0, 0.6);
    let mut common = bodies(false, Vec3::ZERO);
    for side in 0..2 {
        let pose = common[side].pose();
        let velocity = translate.add(omega.cross(pose.position_world()));
        common[side] = RigidBodyState::new(pose, velocity.scale(p[side].mass()), omega).unwrap();
    }
    let model = law(100.0, 0.2, 2.0, 300.0);
    let history = RigidLayerHistory::new(&model, geometry(), 0.0).unwrap();
    let stopped = step_rigid_layer_candidate(
        &common,
        p,
        &history,
        &model,
        0.05,
        10.0,
        law_solve(),
        controls(),
        [0.0; 2],
        || false,
    )
    .unwrap();
    vector_close(stopped.report().applied_force_on_a_n, Vec3::ZERO, 2e-12);
    close(stopped.report().constitutive.dissipated_j, 0.0, 1e-20);
    for i in 0..2 {
        vector_close(
            stopped.bodies()[i].linear_momentum_world(),
            common[i].linear_momentum_world(),
            2e-12,
        );
    }

    let old = bodies(false, Vec3::new(0.0, 1.0, 0.4));
    let first = step_rigid_layer_candidate(
        &old,
        p,
        &history,
        &model,
        0.05,
        10.0,
        law_solve(),
        controls(),
        [0.0; 2],
        || false,
    )
    .unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.71).unwrap();
    let shift = Vec3::new(2.0, 3.0, -4.0);
    let rotated = old.map(|b| {
        RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(b.pose().position_world())
                    .add(shift),
                rotation,
            )
            .unwrap(),
            rotation.rotate_body_to_world(b.linear_momentum_world()),
            b.angular_momentum_body(),
        )
        .unwrap()
    });
    let mut geom = geometry();
    geom.frame = TangentFrame::new(
        array(rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
        array(rotation.rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0))),
    )
    .unwrap();
    let rhistory = RigidLayerHistory::new(&model, geom, 0.0).unwrap();
    let second = step_rigid_layer_candidate(
        &rotated,
        p,
        &rhistory,
        &model,
        0.05,
        10.0,
        law_solve(),
        controls(),
        [0.0; 2],
        || false,
    )
    .unwrap();
    vector_close(
        second.report().applied_force_on_a_n,
        rotation.rotate_body_to_world(first.report().applied_force_on_a_n),
        2e-9,
    );
    close(
        second.report().constitutive.dissipated_j,
        first.report().constitutive.dissipated_j,
        2e-10,
    );
    for i in 0..2 {
        vector_close(
            second.bodies()[i].linear_momentum_world(),
            rotation.rotate_body_to_world(first.bodies()[i].linear_momentum_world()),
            2e-10,
        );
        vector_close(
            second.bodies()[i].angular_momentum_body(),
            first.bodies()[i].angular_momentum_body(),
            2e-10,
        );
    }
}

#[test]
fn g1_retained_elastic_native_mode_refines_at_second_order_without_numerical_heat() {
    let model = law(100.0, 100.0, 2.0, 300.0);
    let time = 0.1_f64;
    let speed = 0.001_f64;
    let frequency = (100.0 * 17.0 / 6.0_f64).sqrt();
    let exact = speed / frequency * (frequency * time).sin();
    let mut errors = Vec::new();
    for n in [16, 32, 64] {
        let mut bodies = bodies(false, Vec3::new(0.0, speed, 0.0));
        let mut history = RigidLayerHistory::new(&model, geometry(), 0.0).unwrap();
        let mut maximum_balance = 0.0_f64;
        for i in 0..n {
            let candidate = step_rigid_layer_candidate(
                &bodies,
                properties(false),
                &history,
                &model,
                time * (i + 1) as f64 / n as f64,
                10.0,
                law_solve(),
                controls(),
                [0.0; 2],
                || false,
            )
            .unwrap();
            assert_eq!(candidate.report().constitutive.dissipated_j, 0.0);
            maximum_balance = maximum_balance.max(candidate.report().total_energy_residual_j.abs());
            (bodies, history, _) = candidate.into_parts();
        }
        errors.push((history.deformation().elastic_slip_m()[0] - exact).abs());
        assert!(maximum_balance < 1e-14);
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    assert!(errors[2] < 1e-9, "{errors:?}");
}

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && 300.0 + h > 0.0 {
            Ok(300.0 + h)
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
        Ok(Some((self.temperature_k(b)? / self.temperature_k(a)?).ln()))
    }
}
static CALORIC: Caloric = Caloric;
const LAYER: CellId = CellId(2);
const K: [f64; 3] = [2.0, 1.0, 0.5];
const CAP: [f64; 3] = [2.0, 3.0, 0.2];
fn thermal() -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..3 {
        let x = [-1.0, 1.0, 0.0][i];
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: CAP[i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [x, 0.0, 0.0],
                    end_m: [x + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
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
                id: 7 + i as u64,
                a: CellId(i as u64),
                b: LAYER,
                area_m2: 1.0,
                conductance_w_k: series_contact_conductance_w_k(1.0, [0.5; 2], [K[i], K[2]], 0.0)
                    .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 32,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn heat_controls() -> HeatedRigidLayerSolve {
    HeatedRigidLayerSolve {
        mechanical: controls(),
        maximum_cells: 3,
        energy_coordinate_j_per_n: 1.0,
        heat: ThirdBodyHeatTolerances {
            cell_j: 2e-10,
            thermal_j: 2e-10,
            layer_j: 2e-10,
            total_j: 2e-10,
        },
        combined_energy_tolerance_j: 5e-10,
        audit_entropy: true,
    }
}
fn source(temperature: f64) -> RigidLayerHeatSource<f64> {
    RigidLayerHeatSource {
        law: law(10.0, 0.2 - 0.5 * (temperature - 300.0), 0.3, temperature),
        conductivity_w_m_k: K,
        area_resistance_m2_k_w: [0.0; 2],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperature,
        },
    }
}

#[test]
fn g1_joint_native_layer_heat_matches_finite_storage_and_temperature_feedback_reference() {
    let thermal = thermal();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, [7, 8], K, [0.0; 2]).unwrap();
    let old = bodies(false, Vec3::new(0.0, 1.0, 0.0));
    let history = RigidLayerHistory::new(&source(300.0).law, geometry(), 0.0).unwrap();
    let result = step_heated_rigid_layer_candidate(
        &old,
        properties(false),
        &history,
        &thermal,
        &route,
        0.1,
        10.0,
        law_solve(),
        heat_controls(),
        [0.0; 2],
        &[0.0; 3],
        |_, t| Ok::<_, &str>(source(t)),
        || false,
    )
    .unwrap();
    let dt = 0.1;
    let g = [0.8, 2.0 / 3.0];
    let effective = CAP[2]
        + dt * g[0] * CAP[0] / (CAP[0] + dt * g[0])
        + dt * g[1] * CAP[1] / (CAP[1] + dt * g[1]);
    let reference = |theta: f64| {
        let y = 0.2 - 0.5 * theta;
        let b = 10.0 * dt / (2.0 * 0.3);
        let a = 1.0 + 10.0 * dt * dt * (17.0 / 6.0) / 4.0;
        let force = (10.0 * dt / 2.0 + b * y) / (a + b);
        let plastic = dt / 0.3 * (force - y);
        let heat = y * plastic + 0.3 / dt * plastic * plastic;
        (force, heat, plastic)
    };
    let (mut lo, mut hi) = (0.0, 0.3);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if effective * mid - reference(mid).1 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let theta = 0.5 * (lo + hi);
    let (force, heat, plastic) = reference(theta);
    let report = result.report();
    close(report.heat.layer_temperature_k, 300.0 + theta, 2e-9);
    close(report.mechanical.applied_force_on_a_n.y, force, 2e-9);
    close(report.heat.dissipated_j, heat, 2e-10);
    close(
        result.history().deformation().plastic_slip_m()[0],
        plastic,
        2e-9,
    );
    close(report.source, 300.0 + theta, 2e-9);
    assert!(
        (force - reference(0.0).0).abs() > 1e-3,
        "feedback axis must be live"
    );
    assert!(report.heat.layer_stored_change_j > 0.0);
    assert!(report.heat.layer_stored_change_j < heat);
    close(report.heat.layer_stored_change_j, CAP[2] * theta, 2e-10);
    for i in 0..2 {
        let theta_body = dt * g[i] * theta / (CAP[i] + dt * g[i]);
        close(
            result.thermal().temperature_k(CellId(i as u64)).unwrap(),
            300.0 + theta_body,
            2e-9,
        );
        close(report.heat.body_heat_j[i], CAP[i] * theta_body, 2e-10);
    }
    close(report.combined_energy_residual_j, 0.0, 5e-10);
    let entropy = report.heat.entropy.as_ref().unwrap();
    close(
        entropy.friction_generation_j_k,
        heat / (300.0 + theta),
        2e-12,
    );
    assert!(!entropy.has_unclassified_sources);
    assert_eq!(entropy.interface_conduction_generation_j_k, 0.0);
    assert!(entropy.network_conduction_generation_j_k > 0.0);
    // Only the existing three cells hold heat; accepted inputs are untouched.
    assert_eq!(result.thermal().cells().len(), 3);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(history.time_s(), 0.0);
    assert_eq!(old, bodies(false, Vec3::new(0.0, 1.0, 0.0)));
    // Continue with the SAME retained strain and freshly resolved accepted law.
    let (next_bodies, next_history, next_thermal, _) = result.into_parts();
    let next = step_heated_rigid_layer_candidate(
        &next_bodies,
        properties(false),
        &next_history,
        &next_thermal,
        &route,
        0.2,
        10.0,
        law_solve(),
        heat_controls(),
        [force, 0.0],
        &[0.0; 3],
        |_, t| Ok::<_, &str>(source(t)),
        || false,
    )
    .unwrap();
    assert!(
        next.history().deformation().plastic_slip_m()[0]
            > next_history.deformation().plastic_slip_m()[0]
    );
}

#[test]
fn g4_late_failures_preserve_original_bodies_history_and_caloric_state() {
    let thermal = thermal();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, [7, 8], K, [0.0; 2]).unwrap();
    let old = bodies(false, Vec3::new(0.0, 1.0, 0.0));
    let history = RigidLayerHistory::new(&source(300.0).law, geometry(), 0.0).unwrap();
    let snapshot = thermal
        .cells()
        .values()
        .map(|c| c.specific_h_j_kg.to_bits())
        .collect::<Vec<_>>();
    let calls = Cell::new(0);
    let fail = step_heated_rigid_layer_candidate(
        &old,
        properties(false),
        &history,
        &thermal,
        &route,
        0.1,
        10.0,
        law_solve(),
        heat_controls(),
        [0.0; 2],
        &[0.0; 3],
        |_, t| {
            calls.set(calls.get() + 1);
            if calls.get() == 3 {
                Err("late source domain")
            } else {
                Ok(source(t))
            }
        },
        || false,
    );
    assert!(matches!(
        fail,
        Err(HeatedRigidLayerError::Source("late source domain"))
    ));
    let mut budget = heat_controls();
    budget.mechanical.maximum_trials = 1;
    assert!(matches!(
        step_heated_rigid_layer_candidate(
            &old,
            properties(false),
            &history,
            &thermal,
            &route,
            0.1,
            10.0,
            law_solve(),
            budget,
            [0.0; 2],
            &[0.0; 3],
            |_, t| Ok::<_, &str>(source(t)),
            || false
        ),
        Err(HeatedRigidLayerError::Mechanical(
            RigidLayerError::TrialBudget
        ))
    ));
    let polls = Cell::new(0);
    assert!(
        step_heated_rigid_layer_candidate(
            &old,
            properties(false),
            &history,
            &thermal,
            &route,
            0.1,
            10.0,
            law_solve(),
            heat_controls(),
            [0.0; 2],
            &[0.0; 3],
            |_, t| Ok::<_, &str>(source(t)),
            || {
                polls.set(polls.get() + 1);
                polls.get() >= 12
            }
        )
        .is_err()
    );
    assert!(polls.get() >= 12);
    assert!(
        step_heated_rigid_layer_candidate(
            &old,
            properties(false),
            &history,
            &thermal,
            &route,
            0.1,
            10.0,
            law_solve(),
            heat_controls(),
            [0.0; 2],
            &[0.0; 3],
            |_, t| {
                let mut s = source(t);
                s.law = law(20.0, 0.2, 0.3, t);
                Ok::<_, &str>(s)
            },
            || false
        )
        .is_err()
    );
    let mut bad_geometry = geometry();
    bad_geometry.area_m2 = 0.5;
    let bad_history = RigidLayerHistory::new(&source(300.0).law, bad_geometry, 0.0).unwrap();
    assert!(
        step_heated_rigid_layer_candidate(
            &old,
            properties(false),
            &bad_history,
            &thermal,
            &route,
            0.1,
            10.0,
            law_solve(),
            heat_controls(),
            [0.0; 2],
            &[0.0; 3],
            |_, t| Ok::<_, &str>(source(t)),
            || false
        )
        .is_err()
    );
    assert_eq!(
        thermal
            .cells()
            .values()
            .map(|c| c.specific_h_j_kg.to_bits())
            .collect::<Vec<_>>(),
        snapshot
    );
    assert_eq!(
        history,
        RigidLayerHistory::new(&source(300.0).law, geometry(), 0.0).unwrap()
    );
    assert_eq!(old, bodies(false, Vec3::new(0.0, 1.0, 0.0)));
    let retry = || {
        step_heated_rigid_layer_candidate(
            &old,
            properties(false),
            &history,
            &thermal,
            &route,
            0.1,
            10.0,
            law_solve(),
            heat_controls(),
            [0.0; 2],
            &[0.0; 3],
            |_, t| Ok::<_, &str>(source(t)),
            || false,
        )
        .unwrap()
    };
    let a = retry();
    let b = retry();
    assert_eq!(a.bodies(), b.bodies());
    assert_eq!(a.history(), b.history());
    for (a, b) in a
        .thermal()
        .cells()
        .values()
        .zip(b.thermal().cells().values())
    {
        assert_eq!(a.specific_h_j_kg.to_bits(), b.specific_h_j_kg.to_bits());
    }
}
