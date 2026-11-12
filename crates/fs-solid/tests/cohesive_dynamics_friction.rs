//! Native deformable damaged-area reclosure: original inertia, damage and work.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    interface_stiffness_per_area_dims, resolve_interface_state_point,
};
use fs_qty::{Dims, Pressure, Temperature};
use fs_solid::cohesive::dynamics::{DynamicsConfig, DynamicsError};
use fs_solid::cohesive::temperature::TemperatureConfig;
use fs_solid::cohesive::{
    Budget, CohesivePatch, CohesiveState, CohesiveSystem, Config, Error, JumpTerm,
};
use fs_solver::{LinearOp, NewtonKrylovConfig};

const AREA: f64 = 0.01;
const KC: f64 = 2000.0;
const COMPRESSION: f64 = 0.001;

#[path = "cohesive_dynamics_friction/coulomb.rs"]
mod coulomb;

fn properties(viscous: bool) -> Vec<(&'static str, Dims, [f64; 2])> {
    let mut rows = vec![
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            [1000.0; 2],
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, [1.0; 2]),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            [0.002; 2],
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            [KC; 2],
        ),
        (RESIDUAL_KINETIC_FRICTION_PROPERTY, Dims::NONE, [0.0; 2]),
    ];
    if viscous {
        rows.extend([
            (
                COHESIVE_NORMAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                [2.0; 2],
            ),
            (
                COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                [4.0; 2],
            ),
        ]);
    }
    rows
}

fn card(mu: [f64; 2], viscous: bool) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, dims, mut values) in properties(viscous) {
        if name == RESIDUAL_KINETIC_FRICTION_PROPERTY {
            values = mu;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured damaged deformable interface; no specimen claim".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "original manufactured bonded member".into(),
            revision: 1,
        },
        texture_frame: "fixed normal/shear axes".into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit original bond and residual surfaces".into(),
            third_body: None,
            environment: "prescribed temperature".into(),
            history: "initially bonded".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(
    card: &InterfaceSystemCard,
    viscous: bool,
    include_mu: bool,
) -> ResolvedBilinearCohesive {
    let requirements: Vec<_> = properties(viscous).into_iter()
        .filter(|(name, _, _)| include_mu || *name != RESIDUAL_KINETIC_FRICTION_PROPERTY)
        // Finite-only extras test that the receiving friction owner enforces
        // its own nonnegative-mu domain, including after a temperature update.
        .map(|(name, dims, _)| ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap())
        .collect();
    let resolved = resolve_interface_state_point(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    if viscous {
        ResolvedBilinearCohesive::try_from_viscous_resolved(
            card,
            resolved,
            CohesiveEnvelope::CommonMode,
            0.02,
            1.0,
        )
        .unwrap()
    } else {
        ResolvedBilinearCohesive::try_from_resolved(card, resolved, 0.02).unwrap()
    }
}

struct Scalar(f64);
impl LinearOp for Scalar {
    fn n(&self) -> usize {
        1
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y[0] = self.0 * x[0];
    }
    // A 1x1 scalar, so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

fn system<'a>(
    elastic: &'a dyn LinearOp,
    source: &'a ResolvedBilinearCohesive,
    free_shear: bool,
) -> CohesiveSystem<'a> {
    let terms = if free_shear {
        vec![JumpTerm {
            dof: 0,
            local: [0.0, 1.0, 0.0],
        }]
    } else {
        vec![]
    };
    let patch = CohesivePatch::from_interface(
        7,
        AREA,
        source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
        terms,
    )
    .unwrap();
    CohesiveSystem::new(
        elastic,
        vec![patch],
        Budget {
            max_dofs: 1,
            max_patches: 1,
            max_jump_terms: 1,
            max_krylov_basis_entries: 128,
        },
    )
    .unwrap()
}

fn config() -> DynamicsConfig {
    DynamicsConfig {
        solve: Config {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 1e-12,
                linear_restart: 2,
                forcing_maximum: 0.01,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 40,
            max_operator_evaluations: 2000,
            force_tolerance_n: 1e-9,
            local_energy_tolerance_j_m2: 1e-10,
            total_cohesive_energy_tolerance_j: 1e-10,
            maximum_displacement_m: 0.02,
        },
        difference_step_m: 1e-8,
        maximum_step_s: 0.02,
        maximum_velocity_m_s: 2.0,
        maximum_steps: 50,
        absolute_energy_tolerance_j: 1e-10,
        relative_energy_tolerance: 1e-8,
    }
}

fn compressed<'m, 'a>(system: &'m CohesiveSystem<'a>, opening: f64) -> CohesiveState<'m, 'a> {
    system
        .bonded_state()
        .solve_over_time_controlled(
            &[0.0],
            &[[opening, 0.0, 0.0]],
            1.0,
            config().solve,
            &mut || false,
        )
        .unwrap()
        .state
        .solve_over_time_controlled(
            &[0.0],
            &[[-COMPRESSION, 0.0, 0.0]],
            1.0,
            config().solve,
            &mut || false,
        )
        .unwrap()
        .state
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn partial_and_complete_reclosure_drive_original_inertia_with_separate_losses() {
    let source = resolve(&card([0.3; 2], true), true, true);
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, true);
    for opening in [0.0, 0.002, 0.005] {
        let accepted = compressed(&system, opening);
        let damage = accepted.points()[0].damage();
        let expected_damage = if opening == 0.0 {
            0.0
        } else if opening < 0.004 {
            2.0 / 3.0
        } else {
            1.0
        };
        close(damage, expected_damage, 1e-13);
        let previous = accepted.points()[0];
        let motion = accepted
            .into_dynamics(&mass, vec![0.02], 2.0, config())
            .unwrap()
            .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
            .unwrap();
        let dt = 0.01;
        let resistance = AREA * 0.3 * damage * KC * COMPRESSION;
        let stiffness = 10.0 + AREA * 1000.0 * (1.0 - damage);
        let viscosity = AREA * 4.0 * (1.0 - damage);
        // Independent scalar momentum equation for the ORIGINAL mass, bulk
        // elasticity, surviving-bond stiffness/viscosity and damaged friction.
        let u = (2.0 * mass.0 * 0.02 / dt - resistance)
            / (2.0 * mass.0 / (dt * dt) + 0.5 * stiffness + viscosity / dt);
        let v = 2.0 * u / dt - 0.02;
        assert!(u / dt > 1e-4);
        let step = motion
            .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], dt, &mut || false)
            .unwrap();
        assert!(step.analytic_jacobian_evaluations > 0);
        close(step.motion.state().displacements_m()[0], u, 1e-12);
        close(step.motion.velocity_m_s()[0], v, 1e-10);
        close(step.friction_dissipated_j, resistance * u, 1e-13);
        close(step.viscous_dissipated_j, viscosity * u * u / dt, 1e-13);
        close(step.energy_residual_j, 0.0, 1e-12);
        assert_eq!(step.cohesive_energy_j.dissipated, 0.0);
        assert_eq!(step.motion.state().points()[0].damage(), previous.damage());
        assert_eq!(
            step.motion.state().points()[0].dissipated_j_m2(),
            previous.dissipated_j_m2()
        );
        assert_eq!(motion.state().points()[0], previous);
        assert_eq!(motion.time_s(), 2.0);
    }
    // A damaged but tensile bond still carries adhesive normal force. That
    // force must not be mistaken for compressive pressure available to friction.
    let tensile = system
        .bonded_state()
        .solve_over_time_controlled(
            &[0.0],
            &[[0.002, 0.0, 0.0]],
            1.0,
            config().solve,
            &mut || false,
        )
        .unwrap()
        .state
        .into_dynamics(&mass, vec![0.02], 0.0, config())
        .unwrap();
    assert!(tensile.state().points()[0].traction_pa()[0] > 0.0);
    let bound = tensile
        .clone()
        .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
        .unwrap();
    let baseline = tensile
        .advance_controlled(&[0.0], &[[0.002, 0.0, 0.0]], 0.01, &mut || false)
        .unwrap();
    let step = bound
        .advance_controlled(&[0.0], &[[0.002, 0.0, 0.0]], 0.01, &mut || false)
        .unwrap();
    assert_eq!(step.friction_dissipated_j, 0.0);
    assert_eq!(step.motion.velocity_m_s(), baseline.motion.velocity_m_s());
    assert_eq!(step.mean_traction_pa, baseline.mean_traction_pa);
}

#[test]
fn growing_damage_uses_integrated_compression_and_actual_prescribed_slip_work() {
    let source = resolve(&card([0.3; 2], false), false, true);
    let elastic = Scalar(1.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, false);
    for direction in [[1.0, 0.0], [0.6, 0.8], [-0.6, 0.8]] {
        let start = [-COMPRESSION, 0.0005 * direction[0], 0.0005 * direction[1]];
        let end = [-COMPRESSION, 0.003 * direction[0], 0.003 * direction[1]];
        let accepted = system
            .bonded_state()
            .solve_controlled(&[0.0], &[start], config().solve, &mut || false)
            .unwrap()
            .state;
        let motion = accepted
            .into_dynamics(&mass, vec![0.0], 0.0, config())
            .unwrap()
            .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
            .unwrap();
        let step = motion
            .advance_controlled(&[0.0], &[end], 0.01, &mut || false)
            .unwrap();
        // Integral over q from .001 to .003 of d(q)=4/3*(1-.001/q).
        // The interval starts at .0005 with zero damaged area.
        let integrated_damage_length = (4.0 / 3.0) * (0.002 - 0.001 * 3.0_f64.ln());
        let friction = AREA * 0.3 * KC * COMPRESSION * integrated_damage_length;
        let intrinsic_work = AREA
            * (0.5 * 1000.0 * (0.001_f64.powi(2) - 0.0005_f64.powi(2))
                + 0.5 * (1.0 + 1.0 / 3.0) * 0.002);
        close(step.friction_dissipated_j, friction, 1e-13);
        close(step.cohesive_energy_j.work, intrinsic_work, 1e-13);
        close(
            step.prescribed_jump_work_j,
            intrinsic_work + friction,
            1e-13,
        );
        close(step.energy_residual_j, 0.0, 1e-12);
        close(step.motion.state().points()[0].damage(), 8.0 / 9.0, 1e-12);
        assert_eq!(step.motion.velocity_m_s(), &[0.0]);
        assert_eq!(step.external_work_j, 0.0);
    }
}

#[test]
fn zero_mu_is_exact_original_motion_and_low_speed_regularization_is_physical_time_based() {
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    for mu in [0.0, 0.3] {
        let source = resolve(&card([mu; 2], false), false, true);
        let system = system(&elastic, &source, true);
        let original = compressed(&system, 0.005)
            .into_dynamics(&mass, vec![0.001], 0.0, config())
            .unwrap();
        let bound = original
            .clone()
            .with_damaged_area_friction_controlled(&[Some(0.01)], &mut || false)
            .unwrap();
        for dt in [0.005, 0.01] {
            let step = bound
                .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], dt, &mut || false)
                .unwrap();
            if mu == 0.0 {
                let baseline = original
                    .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], dt, &mut || false)
                    .unwrap();
                assert_eq!(
                    step.motion.state().displacements_m(),
                    baseline.motion.state().displacements_m()
                );
                assert_eq!(step.motion.velocity_m_s(), baseline.motion.velocity_m_s());
                assert_eq!(
                    step.motion.state().points(),
                    baseline.motion.state().points()
                );
                assert_eq!(step.mean_traction_pa, baseline.mean_traction_pa);
                assert_eq!(
                    step.energy_residual_j.to_bits(),
                    baseline.energy_residual_j.to_bits()
                );
                assert_eq!(step.friction_dissipated_j, 0.0);
            } else {
                let c = AREA * mu * KC * COMPRESSION / 0.01;
                let u = (2.0 * mass.0 * 0.001 / dt) / (2.0 * mass.0 / (dt * dt) + 5.0 + c / dt);
                assert!(u / dt < 0.01);
                close(step.motion.state().displacements_m()[0], u, 1e-13);
                close(step.friction_dissipated_j, c * u * u / dt, 1e-14);
            }
        }
    }
}

#[test]
fn smooth_friction_uses_original_path_jacobian_independently_of_difference_scale() {
    // G1: after fracture, compressive normal motion changes the available
    // friction DURING the step. All three coordinates are native inertial DOFs.
    // The independent triangular momentum equations expose the pressure/slip
    // cross derivative that a frozen-pressure shear tangent would omit.
    struct Diagonal([f64; 3]);
    impl LinearOp for Diagonal {
        fn n(&self) -> usize {
            3
        }
        fn apply(&self, x: &[f64], y: &mut [f64]) {
            for i in 0..3 {
                y[i] = self.0[i] * x[i];
            }
        }
        // Diagonal, so its transpose is itself.
        fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
            self.apply(x, y);
        }
    }
    let source = resolve(&card([0.3; 2], false), false, true);
    let elastic = Diagonal([10.0, 15.0, 15.0]);
    let mass = Diagonal([0.1, 0.2, 0.2]);
    let patch = CohesivePatch::from_interface(
        7,
        AREA,
        &source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
        (0..3)
            .map(|dof| JumpTerm {
                dof,
                local: core::array::from_fn(|i| if i == dof { 1.0 } else { 0.0 }),
            })
            .collect(),
    )
    .unwrap();
    let system = CohesiveSystem::new(
        &elastic,
        vec![patch],
        Budget {
            max_dofs: 3,
            max_patches: 1,
            max_jump_terms: 3,
            max_krylov_basis_entries: 128,
        },
    )
    .unwrap();
    let mut base_config = config();
    base_config.solve.newton.linear_restart = 3;
    let load = [-AREA * KC * COMPRESSION, 0.0, 0.0];
    let offsets = [[-COMPRESSION, 0.0, 0.0]];
    let broken = system
        .bonded_state()
        .solve_controlled(
            &[0.0; 3],
            &[[0.005, 0.0, 0.0]],
            base_config.solve,
            &mut || false,
        )
        .unwrap()
        .state;
    assert!(broken.points()[0].failed());
    let closed = broken
        .solve_controlled(&load, &offsets, base_config.solve, &mut || false)
        .unwrap()
        .state;
    assert_eq!(closed.displacements_m(), &[0.0; 3]);
    let velocity = [0.01, 0.03, 0.04];
    let dt = 0.01;
    let normal = (2.0 * mass.0[0] * velocity[0] / dt)
        / (2.0 * mass.0[0] / (dt * dt) + 0.5 * (elastic.0[0] + AREA * KC));
    let pressure = KC * (COMPRESSION - 0.5 * normal);
    for regularization in [1.0, 1e-4] {
        let expected: [f64; 3] = core::array::from_fn(|i| {
            if i == 0 {
                return normal;
            }
            let inertial = 2.0 * mass.0[i] / (dt * dt);
            if regularization == 1.0 {
                (2.0 * mass.0[i] * velocity[i] / dt)
                    / (inertial
                        + 0.5 * elastic.0[i]
                        + AREA * 0.3 * pressure / (regularization * dt))
            } else {
                (2.0 * mass.0[i] * velocity[i] / dt - AREA * 0.3 * pressure * velocity[i] / 0.05)
                    / (inertial + 0.5 * elastic.0[i])
            }
        });
        let mut reference = None;
        for difference_step in [1e-12, 0.5] {
            // The large value is outside the original adhesive applicability.
            // A smooth analytic path never probes there or needs that knob.
            let mut cfg = base_config;
            cfg.difference_step_m = difference_step;
            let motion = closed
                .clone()
                .into_dynamics(&mass, velocity.to_vec(), 2.0, cfg)
                .unwrap()
                .with_damaged_area_friction_controlled(&[Some(regularization)], &mut || false)
                .unwrap();
            let before = motion.state().points()[0];
            let step = motion
                .advance_controlled(&load, &offsets, dt, &mut || false)
                .unwrap();
            assert!(step.analytic_jacobian_evaluations > 0);
            for (i, &u) in step.motion.state().displacements_m().iter().enumerate() {
                close(u, expected[i], 1e-12);
                close(
                    step.motion.velocity_m_s()[i],
                    2.0 * expected[i] / dt - velocity[i],
                    1e-10,
                );
            }
            let slip = expected[1].hypot(expected[2]);
            let expected_work =
                AREA * 0.3 * pressure * slip * (slip / slip.max(regularization * dt));
            close(step.friction_dissipated_j, expected_work, 1e-13);
            close(step.energy_residual_j, 0.0, 1e-12);
            assert_eq!(step.viscous_dissipated_j, 0.0);
            assert_eq!(step.cohesive_energy_j.dissipated, 0.0);
            assert_eq!(motion.state().points()[0], before);
            let current = (
                step.motion.state().displacements_m().to_vec(),
                step.motion.velocity_m_s().to_vec(),
                step.residual_evaluations,
                step.analytic_jacobian_evaluations,
            );
            if let Some(original) = &reference {
                assert_eq!(original, &current);
            }
            reference = Some(current);
        }
    }
}

#[test]
fn temperature_updates_original_mu_and_bad_data_or_late_cancel_leave_motion_unchanged() {
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    for hot_mu in [0.6, -0.1] {
        let card = card([0.3, hot_mu], false);
        let cold = resolve(&card, false, true);
        let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
        let original = system(&elastic, &cold, true);
        let target = original
            .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
            .unwrap();
        let motion = compressed(&original, 0.005)
            .into_dynamics(&mass, vec![0.02], 2.0, config())
            .unwrap()
            .with_damaged_area_friction_controlled(&[Some(1e-4)], &mut || false)
            .unwrap();
        let point = motion.state().points()[0];
        let temperature_config = TemperatureConfig {
            local_energy_tolerance_j_m2: 1e-10,
            total_energy_tolerance_j: 1e-10,
        };
        let result =
            motion.transfer_temperature_controlled(&target, "T", temperature_config, &mut || false);
        if hot_mu < 0.0 {
            assert!(matches!(
                result,
                Err(DynamicsError::Model(Error::Invalid(_)))
            ));
        } else {
            let transferred = result.unwrap();
            assert_eq!(
                transferred.motion.damaged_area_friction().unwrap()[0]
                    .unwrap()
                    .kinetic_coefficient(),
                hot_mu
            );
            assert_eq!(transferred.motion.velocity_m_s(), motion.velocity_m_s());
            assert_eq!(transferred.motion.time_s(), motion.time_s());
            assert_eq!(transferred.motion.steps(), motion.steps());
            assert_eq!(transferred.energy.parameter_work_j, 0.0);
            assert_eq!(transferred.energy.damage_dissipated_j, 0.0);
            let dt = 0.01;
            let resistance = AREA * hot_mu * KC * COMPRESSION;
            let expected =
                (2.0 * mass.0 * 0.02 / dt - resistance) / (2.0 * mass.0 / (dt * dt) + 5.0);
            let mut polls = 0;
            let step = transferred
                .motion
                .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], dt, &mut || {
                    polls += 1;
                    false
                })
                .unwrap();
            close(step.motion.state().displacements_m()[0], expected, 1e-12);
            let mut count = 0;
            assert!(matches!(
                transferred.motion.advance_controlled(
                    &[0.0],
                    &[[-COMPRESSION, 0.0, 0.0]],
                    dt,
                    &mut || {
                        count += 1;
                        count == polls
                    }
                ),
                Err(DynamicsError::Model(Error::Cancelled))
            ));
            assert_eq!(transferred.motion.steps(), motion.steps());
        }
        assert_eq!(motion.state().points()[0], point);
        assert_eq!(
            motion.damaged_area_friction().unwrap()[0]
                .unwrap()
                .kinetic_coefficient(),
            0.3
        );
    }
    let missing = resolve(&card([0.3; 2], false), false, false);
    let original = system(&elastic, &missing, true);
    let motion = original
        .bonded_state()
        .into_dynamics(&mass, vec![0.0], 0.0, config())
        .unwrap();
    assert!(
        motion
            .clone()
            .with_damaged_area_friction_controlled(&[Some(0.01)], &mut || false)
            .is_err()
    );
    assert!(
        motion
            .clone()
            .with_damaged_area_friction_controlled(&[], &mut || false)
            .is_err()
    );
    assert!(
        motion
            .with_damaged_area_friction_controlled(&[Some(f64::NAN)], &mut || false)
            .is_err()
    );
}
