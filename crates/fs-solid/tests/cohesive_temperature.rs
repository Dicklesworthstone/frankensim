//! G1/G3: temperature-dependent bonds retain the original deformable history.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::{Pressure, Temperature};
use fs_solid::cohesive::dynamics::{DynamicsConfig, DynamicsError};
use fs_solid::cohesive::temperature::TemperatureConfig;
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, Error, JumpTerm};
use fs_solver::{LinearOp, NewtonKrylovConfig};

fn card(provenance: &str) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            [10.0, 4.0],
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, [5.0, 1.0]),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            [5.0, 0.5],
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            [20.0; 2],
        ),
    ] {
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
                    source: provenance.into(),
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
            process: "manufactured thermal bond fixture".into(),
            revision: 1,
        },
        texture_frame: "fixed reference axes".into(),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit bonded interface".into(),
            third_body: None,
            environment: "explicit temperature; manufactured data only".into(),
            history: "initially intact".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &InterfaceSystemCard) -> ResolvedBilinearCohesive {
    resolve_bilinear_cohesive(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        3.0,
    )
    .unwrap()
}
struct Diagonal {
    n: usize,
    value: f64,
}
impl LinearOp for Diagonal {
    fn n(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (out, value) in y.iter_mut().zip(x) {
            *out = self.value * value;
        }
    }
    // A scalar multiple of the identity, so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}
fn budget() -> Budget {
    Budget {
        max_dofs: 2,
        max_patches: 2,
        max_jump_terms: 2,
        max_krylov_basis_entries: 128,
    }
}
fn patch(source: &ResolvedBilinearCohesive, id: u64, area: f64, dof: usize) -> CohesivePatch<'_> {
    CohesivePatch::from_interface(
        id,
        area,
        source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
        vec![JumpTerm {
            dof,
            local: [1.0, 0.0, 0.0],
        }],
    )
    .unwrap()
}
fn config() -> Config {
    Config {
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
        total_cohesive_energy_tolerance_j: 1e-9,
        maximum_displacement_m: 3.0,
    }
}
fn thermal_config() -> TemperatureConfig {
    TemperatureConfig {
        local_energy_tolerance_j_m2: 1e-10,
        total_energy_tolerance_j: 1e-9,
    }
}
fn dynamics_config() -> DynamicsConfig {
    DynamicsConfig {
        solve: config(),
        difference_step_m: 1e-6,
        maximum_step_s: 0.1,
        maximum_velocity_m_s: 2.0,
        maximum_steps: 2,
        absolute_energy_tolerance_j: 1e-9,
        relative_energy_tolerance: 1e-8,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn heating_native_loaded_patches_releases_damage_and_continues_without_rebonding() {
    let card = card("manufactured temperature curve, no measured claim");
    let cold = source(&card);
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    let elastic = Diagonal { n: 2, value: 20.0 };
    let original = CohesiveSystem::new(
        &elastic,
        vec![patch(&cold, 7, 2.0, 0), patch(&cold, 9, 3.0, 1)],
        budget(),
    )
    .unwrap();
    let heated = original
        .with_temperature_sources_controlled(&[&hot, &hot], "T", &mut || false)
        .unwrap();
    // At g=.5 each intact bond is exactly at initiation: traction=5 Pa.
    let loaded = original
        .bonded_state()
        .solve_controlled(&[20.0, 25.0], &[[0.0; 3]; 2], config(), &mut || false)
        .unwrap()
        .state;
    for &u in loaded.displacements_m() {
        close(u, 0.5, 1e-12);
    }
    let step = loaded
        .transfer_temperature_controlled(&heated, "T", thermal_config(), &mut || false)
        .unwrap();
    assert_eq!(step.state.displacements_m(), loaded.displacements_m());
    // Independent fixed-jump energy: psi_cold=1.25, psi_hot(d=0)=.5,
    // d_hot=2/3, psi_hot=1/6, release=1/3, all per original bonded m2.
    close(step.energy.parameter_work_j, -0.75 * 5.0, 1e-11);
    close(step.energy.damage_dissipated_j, 5.0 / 3.0, 1e-11);
    close(step.energy.stored_change_j, 5.0 * (1.0 / 6.0 - 1.25), 1e-11);
    close(step.energy.defect_j, 0.0, 1e-12);
    for point in step.state.points() {
        close(point.damage(), 2.0 / 3.0, 1e-12);
        close(point.traction_pa()[0], 2.0 / 3.0, 1e-12);
        close(point.dissipated_j_m2(), 1.0 / 3.0, 1e-12);
        assert!(core::ptr::eq(point.law(), hot.law()));
    }
    assert_eq!(
        heated.patches()[0].source().unwrap().resolved(),
        hot.resolved()
    );
    // Cooling strengthens the envelope but cannot restore damaged area or
    // rescale loss already accumulated at the hotter source state.
    let cooled = step
        .state
        .transfer_temperature_controlled(&original, "T", thermal_config(), &mut || false)
        .unwrap();
    for (before, after) in step.state.points().iter().zip(cooled.state.points()) {
        close(after.damage(), before.damage(), 1e-12);
        close(after.dissipated_j_m2(), before.dissipated_j_m2(), 1e-12);
    }
    close(cooled.energy.damage_dissipated_j, 0.0, 1e-12);
    // The same native structural solver now carries each bond to complete hot
    // failure. Past thermal release and remaining mechanical loss add once.
    let failed = step
        .state
        .solve_controlled(&[22.0; 2], &[[0.0; 3]; 2], config(), &mut || false)
        .unwrap();
    for point in failed.state.points() {
        assert!(point.failed());
        close(point.dissipated_j_m2(), 2.0 / 3.0, 1e-11);
    }
    close(failed.cohesive_energy_j.dissipated, 5.0 / 3.0, 1e-10);
    // The immutable cold trajectory was not mutated by either candidate.
    for point in loaded.points() {
        close(point.damage(), 0.0, 1e-12);
    }
}

#[test]
fn detached_source_bound_solid_can_change_temperature_beyond_the_adhesive_range() {
    let card = card("manufactured detached thermal continuation");
    let cold = source(&card);
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    let elastic = Diagonal { n: 1, value: 20.0 };
    let original = CohesiveSystem::new(&elastic, vec![patch(&cold, 7, 2.0, 0)], budget()).unwrap();
    let heated = original.with_temperature_sources_controlled(&[&hot], "T", &mut || false).unwrap();
    let failed = original.bonded_state()
        .solve_controlled(&[44.0], &[[0.0; 3]], config(), &mut || false).unwrap().state;
    assert!(failed.points()[0].failed());
    let mut detached_config = config();
    detached_config.maximum_displacement_m = 5.0;
    let far = failed.solve_controlled(&[80.0], &[[0.0; 3]], detached_config, &mut || false).unwrap();
    close(far.state.displacements_m()[0], 4.0, 1e-12);
    assert!(far.state.points()[0].jump_m()[0] > cold.law().maximum_separation_m());
    assert_eq!(far.cohesive_energy_j.dissipated, 0.0);
    let transfer = far.state.transfer_temperature_controlled(&heated, "T", thermal_config(), &mut || false).unwrap();
    assert_eq!(transfer.state.displacements_m(), far.state.displacements_m());
    assert_eq!(transfer.state.points()[0].damage(), 1.0);
    assert_eq!(transfer.state.points()[0].dissipated_j_m2(), failed.points()[0].dissipated_j_m2());
    assert_eq!(transfer.energy.parameter_work_j, 0.0);
    assert_eq!(transfer.energy.damage_dissipated_j, 0.0);
    assert_eq!(transfer.energy.stored_change_j, 0.0);
    assert_eq!(transfer.state.system().patches()[0].source().unwrap().resolved(), hot.resolved());
    let continued = transfer.state.solve_controlled(&[90.0], &[[0.0; 3]], detached_config, &mut || false).unwrap();
    close(continued.state.displacements_m()[0], 4.5, 1e-12);
    assert_eq!(continued.cohesive_energy_j.dissipated, 0.0);
    assert_eq!(continued.state.points()[0].traction_pa(), [0.0; 3]);
}

#[test]
fn thermal_transfer_keeps_actual_mass_velocity_clock_and_step_budget() {
    let card = card("manufactured inertial temperature fixture");
    let cold = source(&card);
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    let elastic = Diagonal { n: 1, value: 20.0 };
    let mass = Diagonal { n: 1, value: 2.0 };
    let original = CohesiveSystem::new(&elastic, vec![patch(&cold, 7, 2.0, 0)], budget()).unwrap();
    let heated = original
        .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
        .unwrap();
    let motion = original
        .bonded_state()
        .into_dynamics(&mass, vec![0.05], 3.0, dynamics_config())
        .unwrap();
    let moving = motion
        .advance_controlled(&[0.0], &[[0.0; 3]], 0.01, &mut || false)
        .unwrap()
        .motion;
    let mut polls = 0;
    let transferred = moving
        .transfer_temperature_controlled(&heated, "T", thermal_config(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(
        transferred.motion.state().displacements_m(),
        moving.state().displacements_m()
    );
    assert_eq!(transferred.motion.velocity_m_s(), moving.velocity_m_s());
    assert_eq!(transferred.motion.time_s(), moving.time_s());
    assert_eq!(transferred.motion.steps(), 1);
    let u0 = moving.state().displacements_m()[0];
    let v0 = moving.velocity_m_s()[0];
    close(transferred.energy.parameter_work_j, -6.0 * u0 * u0, 1e-13);
    assert_eq!(transferred.energy.damage_dissipated_j, 0.0);
    // Analytic midpoint motion with original M=2 and new K=20+2*4=28.
    let dt = 0.01;
    let expected_u = ((1.0 - dt * dt * 28.0 / 8.0) * u0 + dt * v0) / (1.0 + dt * dt * 28.0 / 8.0);
    let expected_v = 2.0 * (expected_u - u0) / dt - v0;
    let advanced = transferred
        .motion
        .advance_controlled(&[0.0], &[[0.0; 3]], dt, &mut || false)
        .unwrap();
    close(
        advanced.motion.state().displacements_m()[0],
        expected_u,
        1e-11,
    );
    close(advanced.motion.velocity_m_s()[0], expected_v, 1e-10);
    assert_eq!(advanced.motion.steps(), 2);
    assert!(matches!(
        advanced
            .motion
            .advance_controlled(&[0.0], &[[0.0; 3]], dt, &mut || false),
        Err(DynamicsError::Model(Error::Budget("trajectory steps")))
    ));
    let mut count = 0;
    assert!(matches!(
        moving.transfer_temperature_controlled(&heated, "T", thermal_config(), &mut || {
            count += 1;
            count == polls
        }),
        Err(DynamicsError::Model(Error::Cancelled))
    ));
    assert_eq!(moving.steps(), 1);
    assert_eq!(moving.velocity_m_s(), &[v0]);
    assert_eq!(moving.state().points()[0].law(), cold.law());
}

#[test]
fn changed_geometry_operator_source_and_late_cancellation_refuse_atomically() {
    let card = card("original thermal interface");
    let cold = source(&card);
    let hot = cold.at_temperature(&card, "T", 400.0).unwrap();
    let other = source(&self::card("different physical source, identical numbers"));
    let elastic = Diagonal { n: 2, value: 20.0 };
    let other_elastic = Diagonal { n: 2, value: 20.0 };
    let original = CohesiveSystem::new(
        &elastic,
        vec![patch(&cold, 7, 2.0, 0), patch(&cold, 9, 3.0, 1)],
        budget(),
    )
    .unwrap();
    assert!(matches!(
        original.with_temperature_sources_controlled(&[&hot, &other], "T", &mut || false),
        Err(Error::Source { .. })
    ));
    assert!(
        original
            .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
            .is_err()
    );
    let loaded = original
        .bonded_state()
        .solve_controlled(&[20.0, 25.0], &[[0.0; 3]; 2], config(), &mut || false)
        .unwrap()
        .state;
    let before = loaded.points().to_vec();
    let heated = original
        .with_temperature_sources_controlled(&[&hot, &hot], "T", &mut || false)
        .unwrap();
    let mut polls = 0;
    loaded
        .transfer_temperature_controlled(&heated, "T", thermal_config(), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        loaded.transfer_temperature_controlled(&heated, "T", thermal_config(), &mut || {
            count += 1;
            count == polls
        }),
        Err(Error::Cancelled)
    ));
    for (id, area, dof) in [(8, 2.0, 0), (7, 2.5, 0), (7, 2.0, 1)] {
        let changed = CohesiveSystem::new(
            &elastic,
            vec![patch(&hot, id, area, dof), patch(&hot, 9, 3.0, 1)],
            budget(),
        )
        .unwrap();
        assert!(matches!(
            loaded.transfer_temperature_controlled(&changed, "T", thermal_config(), &mut || false),
            Err(Error::Invalid(_))
        ));
    }
    let changed = CohesiveSystem::new(
        &other_elastic,
        vec![patch(&hot, 7, 2.0, 0), patch(&hot, 9, 3.0, 1)],
        budget(),
    )
    .unwrap();
    assert!(matches!(
        loaded.transfer_temperature_controlled(&changed, "T", thermal_config(), &mut || false),
        Err(Error::Invalid(_))
    ));
    assert_eq!(loaded.points(), before);
}
