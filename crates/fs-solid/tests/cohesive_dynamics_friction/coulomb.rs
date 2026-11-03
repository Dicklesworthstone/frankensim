//! Native unregularized sticking, breakaway and retained source/history.
use super::*;

#[path = "impulse.rs"]
mod impulse;
use fs_solid::cohesive::dynamics::coulomb::{CoulombActivation, CoulombConfig, CoulombRegime};

#[path = "interval.rs"]
mod interval;

fn sourced(
    static_mu: [f64; 2],
    kinetic_mu: [f64; 2],
    viscous: bool,
) -> (InterfaceSystemCard, ResolvedBilinearCohesive) {
    let base = card(kinetic_mu, viscous);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, static_mu[0]), (400.0, static_mu[1])],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured Coulomb threshold; no specimen calibration".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    let requirements: Vec<_> = properties(viscous)
        .into_iter()
        .chain([(RESIDUAL_STATIC_FRICTION_PROPERTY, Dims::NONE, static_mu)])
        .map(|(name, dims, _)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
        })
        .collect();
    let point = resolve_interface_state_point(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = if viscous {
        ResolvedBilinearCohesive::try_from_viscous_resolved(
            &card,
            point,
            CohesiveEnvelope::CommonMode,
            0.02,
            1.0,
        )
        .unwrap()
    } else {
        ResolvedBilinearCohesive::try_from_resolved(&card, point, 0.02).unwrap()
    };
    (card, source)
}

fn controls() -> CoulombConfig {
    CoulombConfig {
        traction_scale_pa_per_m: 1000.0,
        slip_tolerance_m: 1e-10,
        velocity_tolerance_m_s: 1e-8,
    }
}
fn mechanical() -> DynamicsConfig {
    let mut c = config();
    c.solve.newton.linear_restart = 4;
    c
}

#[test]
fn damaged_bonds_hold_below_static_threshold_then_slide_at_kinetic_capacity() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], true);
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, true);
    for opening in [0.002, 0.005] {
        let accepted = compressed(&system, opening);
        let old = accepted.points()[0];
        let damage = old.damage();
        let motion = accepted
            .into_dynamics(&mass, vec![0.0], 2.0, mechanical())
            .unwrap()
            .with_damaged_area_coulomb_controlled(
                &[Some(CoulombActivation::AtRest)],
                controls(),
                &mut || false,
            )
            .unwrap();
        let static_force = AREA * 0.6 * damage * KC * COMPRESSION;
        let held = motion
            .advance_controlled(
                &[0.5 * static_force],
                &[[-COMPRESSION, 0.0, 0.0]],
                0.01,
                &mut || false,
            )
            .unwrap();
        close(held.motion.state().displacements_m()[0], 0.0, 1e-13);
        close(held.motion.velocity_m_s()[0], 0.0, 1e-11);
        close(
            held.mean_traction_pa[0][1] * AREA,
            0.5 * static_force,
            1e-10,
        );
        assert_eq!(held.friction_dissipated_j, 0.0);
        assert_eq!(held.patch_dissipation[0].friction_j, 0.0);
        assert_eq!(
            held.motion.damaged_area_coulomb_states().unwrap()[0].occupancy,
            CoulombRegime::Sticking
        );
        let force = 1.5 * static_force;
        let step = held
            .motion
            .advance_controlled(&[force], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || false)
            .unwrap();
        let friction = AREA * 0.3 * damage * KC * COMPRESSION;
        let stiffness = elastic.0 + AREA * 1000.0 * (1.0 - damage);
        let viscosity = AREA * 4.0 * (1.0 - damage);
        let expected = (force - friction)
            / (2.0 * mass.0 / 0.01_f64.powi(2) + 0.5 * stiffness + viscosity / 0.01);
        close(step.motion.state().displacements_m()[0], expected, 1e-12);
        close(step.friction_dissipated_j, friction * expected, 1e-13);
        close(
            step.viscous_dissipated_j,
            viscosity * expected * expected / 0.01,
            1e-13,
        );
        close(step.energy_residual_j, 0.0, 1e-12);
        assert_eq!(step.patch_dissipation[0].patch_id, 7);
        assert_eq!(
            step.patch_dissipation[0].fracture_j,
            step.cohesive_energy_j.dissipated
        );
        assert_eq!(
            step.patch_dissipation[0].viscous_j,
            step.viscous_dissipated_j
        );
        assert_eq!(
            step.patch_dissipation[0].friction_j,
            step.friction_dissipated_j
        );
        assert_eq!(
            step.motion.damaged_area_coulomb_states().unwrap()[0].occupancy,
            CoulombRegime::Sliding
        );
        assert_eq!(step.motion.state().points()[0].damage(), damage);
        assert_eq!(
            step.motion.state().points()[0].dissipated_j_m2(),
            old.dissipated_j_m2()
        );
        assert_eq!(motion.state().points()[0], old);
        assert_eq!(motion.steps(), 0);
        assert!(core::ptr::eq(step.motion.state().system(), &system));
    }
}

struct Pair(f64);
impl LinearOp for Pair {
    fn n(&self) -> usize {
        2
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for d in 0..2 {
            y[d] = self.0 * x[d];
        }
    }
}

#[test]
fn two_native_shear_modes_use_one_rotated_coulomb_disk() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], false);
    let elastic = Pair(10.0);
    let mass = Pair(0.1);
    let patch = CohesivePatch::from_interface(
        8,
        AREA,
        &source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
        vec![
            JumpTerm {
                dof: 0,
                local: [0.0, 0.6, 0.8],
            },
            JumpTerm {
                dof: 1,
                local: [0.0, -0.8, 0.6],
            },
        ],
    )
    .unwrap();
    let system = CohesiveSystem::new(
        &elastic,
        vec![patch],
        Budget {
            max_dofs: 2,
            max_patches: 1,
            max_jump_terms: 2,
            max_krylov_basis_entries: 128,
        },
    )
    .unwrap();
    let opened = system
        .bonded_state()
        .solve_controlled(
            &[0.0; 2],
            &[[0.005, 0.0, 0.0]],
            mechanical().solve,
            &mut || false,
        )
        .unwrap();
    let closed = opened
        .state
        .solve_controlled(
            &[0.0; 2],
            &[[-COMPRESSION, 0.0, 0.0]],
            mechanical().solve,
            &mut || false,
        )
        .unwrap();
    let motion = closed
        .state
        .into_dynamics(&mass, vec![0.0; 2], 0.0, mechanical())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::AtRest)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let held = motion
        .advance_controlled(
            &[0.004, -0.003],
            &[[-COMPRESSION, 0.0, 0.0]],
            0.01,
            &mut || false,
        )
        .unwrap();
    for value in held.motion.state().displacements_m() {
        close(*value, 0.0, 1e-13);
    }
    close(held.mean_traction_pa[0][1], 0.48, 1e-9);
    close(held.mean_traction_pa[0][2], 0.14, 1e-9);
    assert_eq!(held.friction_dissipated_j, 0.0);
    let force = [0.03, 0.04];
    let step = motion
        .advance_controlled(&force, &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || false)
        .unwrap();
    for (d, value) in force.iter().enumerate() {
        let expected =
            (value - 0.006 * value / 0.05) / (2.0 * mass.0 / 0.01_f64.powi(2) + 0.5 * elastic.0);
        close(step.motion.state().displacements_m()[d], expected, 1e-12);
    }
    close(step.mean_traction_pa[0][1], -0.168, 1e-8);
    close(step.mean_traction_pa[0][2], 0.576, 1e-8);
    close(step.energy_residual_j, 0.0, 1e-11);
}

#[test]
fn growing_damage_and_opening_use_actual_path_pressure_and_prescribed_work() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], false);
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let prescribed = system(&elastic, &source, false);
    for direction in [[0.6, 0.8], [-0.6, 0.8]] {
        let start = [-COMPRESSION, 0.0005 * direction[0], 0.0005 * direction[1]];
        let end = [-COMPRESSION, 0.003 * direction[0], 0.003 * direction[1]];
        let accepted = prescribed
            .bonded_state()
            .solve_controlled(&[0.0], &[start], mechanical().solve, &mut || false)
            .unwrap();
        let motion = accepted
            .state
            .into_dynamics(&mass, vec![0.0], 0.0, mechanical())
            .unwrap()
            .with_damaged_area_coulomb_controlled(
                &[Some(CoulombActivation::AtRest)],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = motion
            .advance_controlled(&[0.0], &[end], 0.01, &mut || false)
            .unwrap();
        let integrated_damage_length = (4.0 / 3.0) * (0.002 - 0.001 * 3.0_f64.ln());
        let expected = AREA * 0.3 * KC * COMPRESSION * integrated_damage_length;
        close(step.friction_dissipated_j, expected, 1e-12);
        close(
            step.prescribed_jump_work_j,
            step.cohesive_energy_j.work + expected,
            1e-12,
        );
        close(step.motion.state().points()[0].damage(), 8.0 / 9.0, 1e-12);
        close(step.energy_residual_j, 0.0, 1e-12);
    }
    let free = system(&elastic, &source, true);
    for velocity in [0.02, -0.02] {
        let motion = compressed(&free, 0.005)
            .into_dynamics(&mass, vec![velocity], 0.0, mechanical())
            .unwrap()
            .with_damaged_area_coulomb_controlled(
                &[Some(CoulombActivation::Sliding)],
                controls(),
                &mut || false,
            )
            .unwrap();
        let step = motion
            .advance_controlled(&[0.0], &[[COMPRESSION, 0.0, 0.0]], 0.01, &mut || false)
            .unwrap();
        // Damage is one; pressure decreases linearly to zero at half the
        // interval, then stays zero. Its path mean is Kc*compression/4.
        let friction = AREA * 0.3 * KC * COMPRESSION / 4.0;
        let expected = (2.0 * mass.0 * velocity / 0.01 - velocity.signum() * friction)
            / (2.0 * mass.0 / 0.01_f64.powi(2) + 0.5 * elastic.0);
        close(step.motion.state().displacements_m()[0], expected, 1e-12);
        close(step.friction_dissipated_j, friction * expected.abs(), 1e-12);
        assert_eq!(
            step.motion.damaged_area_coulomb_states().unwrap()[0].occupancy,
            CoulombRegime::Inactive
        );
        close(step.energy_residual_j, 0.0, 1e-12);
    }
}

#[test]
fn temperature_refreshes_static_and_kinetic_thresholds_without_rebonding() {
    let (card, source) = sourced([0.6, 0.4], [0.3, 0.2], false);
    let hot = source.at_temperature(&card, "T", 400.0).unwrap();
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let original = system(&elastic, &source, true);
    let target = original
        .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
        .unwrap();
    let motion = compressed(&original, 0.005)
        .into_dynamics(&mass, vec![0.0], 2.0, mechanical())
        .unwrap()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::AtRest)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let cold = motion
        .advance_controlled(&[0.01], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || false)
        .unwrap();
    assert_eq!(cold.friction_dissipated_j, 0.0);
    let transfer = motion
        .transfer_temperature_controlled(
            &target,
            "T",
            TemperatureConfig {
                local_energy_tolerance_j_m2: 1e-10,
                total_energy_tolerance_j: 1e-10,
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(transfer.energy.parameter_work_j, 0.0);
    assert_eq!(transfer.motion.time_s(), motion.time_s());
    assert_eq!(transfer.motion.steps(), motion.steps());
    assert_eq!(transfer.motion.velocity_m_s(), motion.velocity_m_s());
    assert_eq!(
        transfer.motion.damaged_area_coulomb_states(),
        motion.damaged_area_coulomb_states()
    );
    let step = transfer
        .motion
        .advance_controlled(&[0.01], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || false)
        .unwrap();
    let expected = (0.01 - 0.004) / 2005.0;
    close(step.motion.state().displacements_m()[0], expected, 1e-12);
    close(step.friction_dissipated_j, 0.004 * expected, 1e-13);
    assert_eq!(step.motion.state().points()[0].damage(), 1.0);
    assert_eq!(motion.state().points()[0].damage(), 1.0);
}

#[test]
fn resolved_motion_keeps_kinetic_heat_below_the_displacement_tolerance() {
    let (_, source) = sourced([0.6; 2], [0.3; 2], false);
    let elastic = Scalar(0.0);
    let mass = Scalar(0.1);
    let system = system(&elastic, &source, true);
    let accepted = compressed(&system, 0.005);
    for (dt, initial_velocity, force, activation) in [
        (1e-7, 4e-4, 0.006, CoulombActivation::Sliding),
        (1e-8, 4e-4, 0.006, CoulombActivation::Sliding),
        (1e-5, 0.0, 0.02, CoulombActivation::AtRest),
    ] {
        let mut mechanical = mechanical();
        mechanical.difference_step_m = 1e-14;
        let mut disk = controls();
        // Balance numerical traction coordinates against this physical mass
        // scale. This changes conditioning only, not either friction capacity.
        disk.traction_scale_pa_per_m = 2.0 * mass.0 / (AREA * dt * dt);
        let motion = accepted
            .clone()
            .into_dynamics(&mass, vec![initial_velocity], 0.0, mechanical)
            .unwrap()
            .with_damaged_area_coulomb_controlled(&[Some(activation)], disk, &mut || false)
            .unwrap();
        let step = motion
            .advance_controlled(&[force], &[[-COMPRESSION, 0.0, 0.0]], dt, &mut || false)
            .unwrap();
        let expected_slip = dt * initial_velocity + (force - 0.006) * dt * dt / (2.0 * mass.0);
        assert!(expected_slip > 0.0 && expected_slip < disk.slip_tolerance_m);
        let expected_velocity = 2.0 * expected_slip / dt - initial_velocity;
        assert!(expected_velocity > disk.velocity_tolerance_m_s);
        close(
            step.motion.state().displacements_m()[0],
            expected_slip,
            1e-16,
        );
        close(step.motion.velocity_m_s()[0], expected_velocity, 1e-10);
        close(step.friction_dissipated_j, 0.006 * expected_slip, 1e-19);
        assert!(step.friction_dissipated_j > 0.0);
        close(step.mean_traction_pa[0][1], 0.6, 1e-7);
        assert_eq!(
            step.motion.damaged_area_coulomb_states().unwrap()[0].occupancy,
            CoulombRegime::Sliding
        );
    }
}

#[test]
fn zero_capacity_matches_original_and_refusals_preserve_complete_motion() {
    let elastic = Scalar(10.0);
    let mass = Scalar(0.1);
    let (_, zero) = sourced([0.0; 2], [0.0; 2], false);
    let zero_system = system(&elastic, &zero, true);
    let original = compressed(&zero_system, 0.005)
        .into_dynamics(&mass, vec![0.02], 0.0, mechanical())
        .unwrap();
    let bound = original
        .clone()
        .with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::Sliding)],
            controls(),
            &mut || false,
        )
        .unwrap();
    let baseline = original
        .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || false)
        .unwrap();
    let zero_step = bound
        .advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || false)
        .unwrap();
    assert_eq!(
        baseline.motion.state().displacements_m(),
        zero_step.motion.state().displacements_m()
    );
    assert_eq!(
        baseline.motion.velocity_m_s(),
        zero_step.motion.velocity_m_s()
    );
    assert_eq!(
        baseline.motion.state().points(),
        zero_step.motion.state().points()
    );
    assert_eq!(
        baseline.energy_residual_j.to_bits(),
        zero_step.energy_residual_j.to_bits()
    );
    for (static_mu, kinetic_mu) in [(0.6, 0.3), (0.2, 0.3)] {
        let (_, source) = sourced([static_mu; 2], [kinetic_mu; 2], false);
        let system = system(&elastic, &source, true);
        let accepted = compressed(&system, 0.005);
        let motion = accepted
            .clone()
            .into_dynamics(&mass, vec![0.0], 0.0, mechanical())
            .unwrap();
        let bound = motion.clone().with_damaged_area_coulomb_controlled(
            &[Some(CoulombActivation::AtRest)],
            controls(),
            &mut || false,
        );
        if static_mu < kinetic_mu {
            assert!(bound.is_err());
            continue;
        }
        let bound = bound.unwrap();
        let mut polls = 0;
        bound
            .advance_controlled(&[0.02], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || {
                polls += 1;
                false
            })
            .unwrap();
        for stop in [polls / 2, polls - 1, polls] {
            let mut count = 0;
            assert!(matches!(
                bound.advance_controlled(&[0.02], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || {
                    count += 1;
                    count == stop
                }),
                Err(DynamicsError::Model(Error::Cancelled))
            ));
        }
        assert_eq!(bound.state().points(), motion.state().points());
        assert_eq!(bound.velocity_m_s(), motion.velocity_m_s());
        assert_eq!(bound.steps(), 0);
        let mut budget = mechanical();
        budget.solve.newton.linear_restart = 9; // original 110, augmented 150 > admitted 128.
        let small = accepted
            .clone()
            .into_dynamics(&mass, vec![0.0], 0.0, budget)
            .unwrap();
        assert!(matches!(
            small.with_damaged_area_coulomb_controlled(
                &[Some(CoulombActivation::AtRest)],
                controls(),
                &mut || false
            ),
            Err(DynamicsError::Model(Error::Budget(_)))
        ));
        assert!(
            bound
                .clone()
                .with_damaged_area_friction_controlled(&[Some(0.01)], &mut || false)
                .is_err()
        );
        // Both a zero-net-slip capture and a finite forward-slip midpoint
        // can reverse resolved velocity. Neither may manufacture path heat
        // across an unresolved physical stopping event.
        for initial_velocity in [1e-4, 4e-4] {
            if initial_velocity == 4e-4 {
                let predicted_slip = (2.0 * mass.0 * initial_velocity / 0.01 - 0.006) / 2005.0;
                assert!(predicted_slip > controls().slip_tolerance_m);
                assert!(2.0 * predicted_slip / 0.01 - initial_velocity < 0.0);
            }
            let stopping = accepted
                .clone()
                .into_dynamics(&mass, vec![initial_velocity], 0.0, mechanical())
                .unwrap()
                .with_damaged_area_coulomb_controlled(
                    &[Some(CoulombActivation::Sliding)],
                    controls(),
                    &mut || false,
                )
                .unwrap();
            assert!(
                matches!(stopping.advance_controlled(&[0.0], &[[-COMPRESSION, 0.0, 0.0]], 0.01, &mut || false),
                Err(DynamicsError::Model(Error::UnresolvedFrictionReversal { patch: 7 })))
            );
            assert_eq!(stopping.velocity_m_s(), &[initial_velocity]);
        }
    }
    let missing = resolve(&card([0.3; 2], false), false, true);
    let missing_system = system(&elastic, &missing, true);
    assert!(
        missing_system
            .bonded_state()
            .into_dynamics(&mass, vec![0.0], 0.0, mechanical())
            .unwrap()
            .with_damaged_area_coulomb_controlled(
                &[Some(CoulombActivation::AtRest)],
                controls(),
                &mut || false
            )
            .is_err()
    );
}
