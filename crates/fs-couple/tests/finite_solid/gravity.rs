//! Native dynamic gravity: independent trajectory, mass moments and atomic refusal.
//! These manufactured conservative fields are numerical oracles, not material data.
use super::*;
use fs_couple::finite_solid::contact::step::free::FreeFlightControls;
use fs_couple::finite_solid::gravity::StationaryGravity;
use fs_couple::finite_solid::impulse::FiniteSolidPoint;
use fs_ga::Mat3;
use fs_mbd::gravity::{GravityError, PointGravityField, PointGravitySample, UniformGravity};
use fs_mbd::{Gravity, Vec3};
use fs_mbd_material::density_field::{DensityVertexBinding, TetrahedralDensityBody};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use std::cell::Cell;
use std::sync::Arc;

const FRAME: u64 = 71;

fn vector(x: [f64; 3]) -> Vec3 {
    Vec3::new(x[0], x[1], x[2])
}
fn components(x: Vec3) -> [f64; 3] {
    [x.x, x.y, x.z]
}
fn uniform(g: [f64; 3]) -> UniformGravity {
    UniformGravity::new(FRAME, Gravity::new(vector(g)).unwrap(), Vec3::ZERO).unwrap()
}
fn field<'a>(source: &'a dyn PointGravityField) -> StationaryGravity<'a> {
    StationaryGravity::new(source, FRAME)
}
fn graded_density(rho: [f64; 4]) -> Arc<TetrahedralDensityBody> {
    let cards = rho.map(|value| source_card(value, 1.0, false, false));
    let context = point(300.0);
    let bindings: Vec<_> = (0..4)
        .map(|node| DensityVertexBinding {
            position_m: vector(NODES[node]),
            material: &cards[node],
            point: &context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    Arc::new(TetrahedralDensityBody::resolve(&bindings, &TETS, 4, || false).unwrap())
}

/// phi(x) = k x^2/2 - g0 x; all quantities follow this one scalar potential.
struct Tidal {
    g0: f64,
    k: f64,
}
impl PointGravityField for Tidal {
    fn sample_point(&self, frame: u64, x: Vec3) -> Result<PointGravitySample, GravityError> {
        if frame != FRAME {
            return Err(GravityError::FrameMismatch);
        }
        PointGravitySample {
            acceleration_world_m_s2: Vec3::new(self.g0 - self.k * x.x, 0.0, 0.0),
            specific_potential_j_kg: 0.5 * self.k * x.x * x.x - self.g0 * x.x,
            acceleration_gradient_s2: Mat3 {
                m: [-self.k, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            },
        }
        .validate()
    }
}

#[test]
fn bound_native_uniform_field_drives_free_fall_and_mixed_free_flight_account() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let g = [1.5, -2.0, 0.75];
        let source = uniform(g);
        let u0 = [0.2, -0.3, 0.4];
        let v0 = [-0.2, 0.1, 0.3];
        let mut body = FiniteSolidBody::new_with_gravity(
            problem(&material, &[], &[]),
            field(&source),
            &[u0; 4],
            &[v0; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        for (dof, actual) in body.native_state().a.iter().enumerate() {
            near(*actual, g[dof % 3], 2e-11);
        }
        let origin = [0.1, 0.4, -0.2];
        let h = 0.025;
        for step in 0..8 {
            let old = body.native_state().clone();
            let (next, rigid, report) = body
                .rigid_free_flight_candidate(
                    &[],
                    &[],
                    &[],
                    [0.0; 3],
                    &[],
                    Gravity::ZERO,
                    h,
                    origin,
                    FreeFlightControls {
                        solid: controls(1e-9),
                        energy_tolerance_j: 1e-9,
                        ..Default::default()
                    },
                    cx,
                )
                .unwrap();
            assert!(rigid.is_empty());
            assert_eq!(body.native_state(), &old);
            let t = (step + 1) as f64 * h;
            for node in 0..4 {
                for axis in 0..3 {
                    near(
                        next.displacement_m()[node][axis],
                        u0[axis] + t * v0[axis] + 0.5 * g[axis] * t * t,
                        2e-10,
                    );
                    near(
                        next.velocity_m_s()[node][axis],
                        v0[axis] + g[axis] * t,
                        2e-10,
                    );
                    near(next.native_state().a[3 * node + axis], g[axis], 2e-9);
                }
            }
            let gravity = report.solid.gravity.as_ref().unwrap();
            let center_mid: [f64; 3] = std::array::from_fn(|axis| {
                0.25 + u0[axis]
                    + v0[axis] * (t - 0.5 * h)
                    + 0.25 * g[axis] * (t * t + (t - h).powi(2))
            });
            let expected_impulse = vector(g).scale(h); // original total mass is 1 kg
            let expected_moment = vector(center_mid)
                .sub(vector(origin))
                .cross(expected_impulse);
            for axis in 0..3 {
                near(
                    components(report.external_impulse_kg_m_s)[axis],
                    components(expected_impulse)[axis],
                    2e-11,
                );
                near(
                    components(report.external_angular_impulse_kg_m2_s)[axis],
                    components(expected_moment)[axis],
                    2e-11,
                );
            }
            near(report.absolute_linear_residual_kg_m_s, 0.0, 2e-10);
            near(report.absolute_angular_residual_kg_m2_s, 0.0, 2e-10);
            near(report.solid.integration_defect_j, 0.0, 2e-10);
            near(gravity.quadrature_defect_j, 0.0, 2e-12);
            near(report.source_work_j, gravity.work_j, 0.0);
            assert_eq!(gravity.frame, FRAME);
            assert_eq!(next.stationary_gravity().unwrap().frame(), FRAME);
            body = next;
        }
        assert_eq!(body.native_state().steps, 8);
        assert_eq!(body.native_state().history.len(), 8);
    });
}

/// Exact P1 reference integrals for the original nodal density, independent of
/// sparse mass, native quadrature, force assembly and the production Newton solve.
fn moments(rho: [f64; 4]) -> (f64, [f64; 4], [f64; 4]) {
    let sum: f64 = rho.iter().sum();
    let row_sum = rho.map(|r| VOLUME * (sum + r) / 20.0);
    let column_one = std::array::from_fn(|node| {
        if node == 1 {
            VOLUME * (sum + 2.0 * rho[1]) / 60.0
        } else {
            VOLUME * (sum + rho[node] + rho[1]) / 120.0
        }
    });
    (VOLUME * sum / 4.0, row_sum, column_one)
}

fn tidal_oracle(
    old: f64,
    v: f64,
    h: f64,
    mass: f64,
    first_moment: f64,
    source: &Tidal,
    held: f64,
    pressure: f64,
    energy_momentum: bool,
) -> (f64, f64) {
    let residual = |next: f64| {
        let midpoint = 0.5 * (old + next);
        let material = if energy_momentum && (next - old).abs() > 1e-9 {
            (strain_energy(next, 1.0) - strain_energy(old, 1.0)) / (next - old)
        } else {
            VOLUME * nominal_stress(midpoint, 1.0)
        };
        2.0 * mass * (next - old - h * v) / (h * h) + material + pressure * VOLUME
            - held
            - source.g0 * first_moment
            + source.k * mass * midpoint
    };
    let (mut low, mut high) = (0.5, 2.5);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if residual(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let next = 0.5 * (low + high);
    (next, 2.0 * (next - old) / h - v)
}

#[test]
fn evolving_tidal_gravity_matches_scalar_motion_original_p1_moments_and_pressure() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|dof| *dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let source = Tidal {
            g0: 100.0,
            k: 200.0,
        };
        let held = [(3, 0.2)];
        let pressure = 0.7;
        let surface = Arc::new(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap(),
        );
        for graded in [false, true] {
            let rho = if graded {
                [3.0, 15.0, 3.0, 3.0]
            } else {
                [6.0; 4]
            };
            let (total_mass, row_sum, column_one) = moments(rho);
            let mass = column_one[1];
            for energy_momentum in [false, true] {
                let mut u = [[0.0; 3]; 4];
                let mut v = [[0.0; 3]; 4];
                let (mut stretch, mut velocity) = (1.35, 0.1);
                u[1][0] = stretch - 1.0;
                v[1][0] = velocity;
                let mut body = FiniteSolidBody::new_with_gravity_and_sources(
                    problem(&material, &fixed, &held),
                    graded.then(|| graded_density(rho)),
                    Some((Arc::clone(&surface), TetPressure::Uniform(pressure))),
                    field(&source),
                    &u,
                    &v,
                    0.0,
                    Default::default(),
                    cx,
                )
                .unwrap();
                let initial_force = source.g0 * row_sum[1] - source.k * mass * stretch;
                near(
                    body.native_state().a[0],
                    (held[0].1 + initial_force
                        - VOLUME * nominal_stress(stretch, 1.0)
                        - pressure * VOLUME)
                        / mass,
                    2e-9,
                );
                let potential = |lambda: f64| {
                    0.5 * source.k * mass * lambda * lambda - source.g0 * row_sum[1] * lambda
                };
                for _ in 0..6 {
                    let before = body.native_state().clone();
                    let h = 0.15;
                    let (expected_q, expected_v) = tidal_oracle(
                        stretch,
                        velocity,
                        h,
                        mass,
                        row_sum[1],
                        &source,
                        held[0].1,
                        pressure,
                        energy_momentum,
                    );
                    let mut limits = controls(if energy_momentum { 2e-9 } else { 1e-2 });
                    // A doubled/missing tidal tangent does not converge within
                    // this native Newton budget for the deliberately stiff field.
                    limits.solve.max_newton_iterations = 6;
                    if energy_momentum {
                        // Native trust-region acceptance evaluates an old-point
                        // Jacobian after the candidate residual. Receipt forces
                        // must remain those of the accepted residual.
                        limits.solve.newton.globalization = fs_solver::Globalization::TrustRegion(
                            fs_solver::TrustRegionConfig::default(),
                        );
                    }
                    let (next, report) = if energy_momentum {
                        body.step_energy_momentum_candidate_with_pressure(
                            h, &held, [0.0; 3], pressure, limits, cx,
                        )
                    } else {
                        body.step_candidate_with_pressure(h, &held, [0.0; 3], pressure, limits, cx)
                    }
                    .unwrap();
                    assert_eq!(body.native_state(), &before);
                    near(next.native_state().q[0] + 1.0, expected_q, 2e-9);
                    near(next.native_state().v[0], expected_v, 2e-8);
                    let gravity = report.gravity.as_ref().unwrap();
                    near(gravity.potential_before_j, potential(stretch), 2e-11);
                    near(gravity.potential_after_j, potential(expected_q), 2e-8);
                    near(
                        gravity.work_j,
                        potential(stretch) - potential(expected_q),
                        2e-8,
                    );
                    near(gravity.quadrature_defect_j, 0.0, 3e-11);
                    let midpoint = 0.5 * (stretch + expected_q);
                    for node in 0..4 {
                        near(
                            gravity.midpoint_force_n[3 * node],
                            source.g0 * row_sum[node] - source.k * column_one[node] * midpoint,
                            3e-8,
                        );
                    }
                    let dx = expected_q - stretch;
                    near(report.pressure_work_j, -pressure * VOLUME * dx, 2e-10);
                    near(
                        report.source_work_j,
                        held[0].1 * dx - pressure * VOLUME * dx + potential(stretch)
                            - potential(expected_q),
                        3e-8,
                    );
                    let support_x: f64 = report
                        .support_reactions_n
                        .iter()
                        .filter(|(dof, _)| dof % 3 == 0)
                        .map(|(_, force)| force)
                        .sum();
                    let total_gravity = source.g0 * total_mass - source.k * row_sum[1] * midpoint;
                    near(
                        support_x + held[0].1 + total_gravity,
                        row_sum[1] * (expected_v - velocity) / h,
                        8e-8,
                    );
                    body = next;
                    stretch = expected_q;
                    velocity = expected_v;
                }
            }
        }
    });
}

#[test]
fn impulse_retains_the_bound_field_and_native_acceleration_history() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let g = [0.0, -3.0, 0.0];
        let source = uniform(g);
        let body = FiniteSolidBody::new_with_gravity(
            problem(&material, &[], &[]),
            field(&source),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (body, _) = body
            .step_candidate(0.02, &[], [0.0; 3], controls(1e-9), cx)
            .unwrap();
        let before = body.native_state().clone();
        let impulse = [0.2, 0.1, -0.3];
        let (kicked, _) = body
            .point_impulses_candidate(
                &[(FiniteSolidPoint::new(0, [0.25; 4]).unwrap(), impulse)],
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        assert_eq!(body.native_state(), &before);
        assert_eq!(kicked.native_state().a, before.a);
        assert_eq!(kicked.native_state().history, before.history);
        assert_eq!(kicked.stationary_gravity().unwrap().frame(), FRAME);
        let h = 0.03;
        let (next, _) = kicked
            .step_candidate(h, &[], [0.0; 3], controls(1e-9), cx)
            .unwrap();
        for node in 0..4 {
            for axis in 0..3 {
                near(
                    next.velocity_m_s()[node][axis],
                    g[axis] * 0.05 + impulse[axis],
                    2e-10,
                );
                near(
                    next.displacement_m()[node][axis],
                    0.5 * g[axis] * 0.05_f64.powi(2) + impulse[axis] * h,
                    2e-10,
                );
            }
        }
        assert_eq!(next.native_state().steps, 2);
    });
}

struct Limited {
    inner: UniformGravity,
    maximum_x_m: f64,
}
impl PointGravityField for Limited {
    fn sample_point(&self, frame: u64, x: Vec3) -> Result<PointGravitySample, GravityError> {
        if x.x > self.maximum_x_m {
            return Err(GravityError::Invalid("outside authored field domain"));
        }
        self.inner.sample_point(frame, x)
    }
}

#[test]
fn inadmissible_endpoint_frame_and_unsupported_split_refuse_without_state_loss() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let source = Limited {
            inner: uniform([0.0; 3]),
            maximum_x_m: 0.74,
        };
        let p = problem(&material, &[], &[]);
        let zero = [[0.0; 3]; 4];
        assert!(
            FiniteSolidBody::new_with_gravity(
                p,
                StationaryGravity::new(&source, FRAME + 1),
                &zero,
                &zero,
                0.0,
                Default::default(),
                cx,
            )
            .is_err()
        );
        let body = FiniteSolidBody::new_with_gravity(
            p,
            field(&source),
            &zero,
            &[[1.0, 0.0, 0.0]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (body, _) = body
            .step_candidate(0.02, &[], [0.0; 3], controls(1e-9), cx)
            .unwrap();
        let before = body.native_state().clone();
        // Native midpoint samples are admitted; accepted endpoint samples leave
        // the field domain. Even after a prior accepted step, nothing publishes.
        assert!(
            body.step_candidate(0.2, &[], [0.0; 3], controls(1e-9), cx)
                .is_err()
        );
        assert_eq!(body.native_state(), &before);
        assert!(
            body.force_kick_candidate(0.01, &[], [0.0; 3], [0.0; 3], Default::default(), cx,)
                .is_err()
        );
        assert_eq!(body.native_state(), &before);
        // The retained field is still usable after the refused candidates.
        let (next, _) = body
            .step_candidate(0.01, &[], [0.0; 3], controls(1e-9), cx)
            .unwrap();
        near(next.native_state().t, 0.03, 1e-15);
        assert_eq!(next.native_state().history.len(), 2);
    });
}

struct Cancelling<'a> {
    inner: UniformGravity,
    gate: &'a CancelGate,
    armed: Cell<bool>,
}
impl PointGravityField for Cancelling<'_> {
    fn sample_point(&self, frame: u64, x: Vec3) -> Result<PointGravitySample, GravityError> {
        if self.armed.get() {
            self.gate.request();
        }
        self.inner.sample_point(frame, x)
    }
}

#[test]
fn cancellation_during_native_field_sampling_keeps_accepted_motion_and_history() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let source = Cancelling {
            inner: uniform([0.0, -1.0, 0.0]),
            gate,
            armed: Cell::new(false),
        };
        let body = FiniteSolidBody::new_with_gravity(
            problem(&material, &[], &[]),
            field(&source),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let (body, _) = body
            .step_candidate(0.02, &[], [0.0; 3], controls(1e-9), cx)
            .unwrap();
        let before = body.native_state().clone();
        source.armed.set(true);
        let error = body
            .step_candidate(0.02, &[], [0.0; 3], controls(1e-9), cx)
            .unwrap_err();
        assert!(
            error.to_string().to_ascii_lowercase().contains("cancel"),
            "{error}"
        );
        assert_eq!(body.native_state(), &before);
        assert_eq!(body.native_state().history.len(), 1);
    });
}

struct QuadraticAcceleration {
    coefficient: f64,
}
impl PointGravityField for QuadraticAcceleration {
    fn sample_point(&self, frame: u64, x: Vec3) -> Result<PointGravitySample, GravityError> {
        if frame != FRAME {
            return Err(GravityError::FrameMismatch);
        }
        PointGravitySample {
            acceleration_world_m_s2: Vec3::new(self.coefficient * x.x * x.x, 0.0, 0.0),
            specific_potential_j_kg: -self.coefficient * x.x * x.x * x.x / 3.0,
            acceleration_gradient_s2: Mat3 {
                m: [2.0 * self.coefficient * x.x, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            },
        }
        .validate()
    }
}

#[test]
fn nonlinear_field_retains_its_potential_work_defect_and_refuses_a_tight_energy_gate() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|dof| *dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        // Grade the fixed vertices while retaining a uniform x marginal:
        // sum(rho)=24 and rho_1=6. The original P1 rule then integrates this
        // cubic x potential exactly, and total/reference x mass stay unchanged.
        let density = graded_density([3.0, 6.0, 9.0, 6.0]);
        let source = QuadraticAcceleration { coefficient: 800.0 };
        let mut u = [[0.0; 3]; 4];
        let mut v = [[0.0; 3]; 4];
        let old = 1.1;
        u[1][0] = old - 1.0;
        v[1][0] = 0.5;
        let body = FiniteSolidBody::new_with_gravity_and_sources(
            problem(&material, &fixed, &[]),
            Some(density),
            None,
            field(&source),
            &u,
            &v,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let before = body.native_state().clone();
        let rejected = body
            .step_energy_momentum_candidate(0.02, &[], [0.0; 3], controls(1e-10), cx)
            .unwrap_err();
        assert!(matches!(rejected, FiniteSolidError::Energy { .. }));
        assert_eq!(body.native_state(), &before);
        let (next, report) = body
            .step_energy_momentum_candidate(0.02, &[], [0.0; 3], controls(1.0), cx)
            .unwrap();
        let new = next.native_state().q[0] + 1.0;
        // Exact reference integral of rho*N_1^3 = mass/20 = 0.05 kg.
        let cubic_mass = 0.05;
        let potential = |stretch: f64| {
            -source.coefficient * cubic_mass * stretch.powi(3) / 3.0
        };
        let gravity = report.gravity.as_ref().unwrap();
        near(gravity.potential_before_j, potential(old), 2e-10);
        near(gravity.potential_after_j, potential(new), 2e-10);
        near(gravity.work_j, potential(old) - potential(new), 2e-10);
        // Midpoint integration of this quadratic acceleration differs from
        // its cubic potential by -k*mass_3*(Delta stretch)^3/12.
        let defect = -source.coefficient * cubic_mass * (new - old).powi(3) / 12.0;
        assert!(defect.abs() > 1e-6);
        near(gravity.quadrature_defect_j, defect, 2e-9);
        near(report.integration_defect_j, defect, 2e-9);
        near(report.source_work_j, gravity.work_j, 0.0);
        assert_eq!(body.native_state(), &before);
    });
}
