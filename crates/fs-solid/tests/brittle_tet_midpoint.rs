//! Native sourced P1 midpoint dynamics: physical work, inertia and history.
//! Manufactured coefficients test mechanics; they do not qualify a material.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::source::dynamics::resolve_brittle_dynamic_material;
use fs_material::brittle::source::*;
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use fs_solid::brittle::solid::{BrittleTetError, BrittleTetSolid};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::IncrementSettings;

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const YOUNG: f64 = 1200.0;
const DENSITY: f64 = 6.0;
const STRENGTH: f64 = 1.2;
const GC: f64 = 0.0012;
const VOLUME: f64 = 1.0 / 6.0;

fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 79,
                kernel_id: 79,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate)
    })
}
fn budget() -> TetAssemblyBudget {
    TetAssemblyBudget {
        maximum_nodes: 4,
        maximum_tetrahedra: 1,
        maximum_free_dofs: 12,
        minimum_scaled_jacobian: 1e-10,
    }
}
fn settings() -> IncrementSettings {
    IncrementSettings {
        max_newton_iterations: 30,
        relative_tolerance: 1e-11,
        max_backtracks: 18,
        linear_tolerance: 1e-13,
        max_linear_iterations: 120,
    }
}
fn controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-10,
        equilibrium_work_tolerance_j: 1e-10,
        kinetic_work_tolerance_j: 1e-10,
        total_energy_tolerance_j: 1e-10,
        difference_step_strain: Some(1e-7),
    }
}
fn source(geometry: &TetLoadGeometry3<'_>) -> ResolvedBrittleCrackBand {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, YOUNG),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, DENSITY),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, STRENGTH),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 1.5),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            GC,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            2.0 * GC,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured P1 midpoint fixture; not measured data".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle solid".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1., 0., 0.]).unwrap();
    resolve_brittle_dynamic_material(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.02, 1e-9),
    )
    .unwrap()
}
fn prescribed(extension: f64) -> Vec<(usize, f64)> {
    (0..12)
        .map(|dof| (dof, if dof == 3 { extension } else { 0.0 }))
        .collect()
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:e} != {expected:e}"
    );
}
fn momentum(velocity: &[[f64; 3]]) -> [f64; 3] {
    core::array::from_fn(|c| (DENSITY * VOLUME / 4.0) * velocity.iter().map(|v| v[c]).sum::<f64>())
}
fn quadrature_kinetic(velocity: &[[f64; 3]]) -> f64 {
    // Independent positive degree-two quadrature, including velocity cross
    // terms. Row-sum lumping would give a different answer for these fields.
    let a = 0.5854101966249685;
    let b = 0.1381966011250105;
    let mut value = 0.0;
    for q in 0..4 {
        let v: [f64; 3] = core::array::from_fn(|c| {
            (0..4)
                .map(|node| (if node == q { a } else { b }) * velocity[node][c])
                .sum()
        });
        value += (DENSITY * VOLUME / 8.0) * v.iter().map(|v| v * v).sum::<f64>();
    }
    value
}

#[test]
fn free_elastic_vibration_preserves_consistent_energy_and_linear_momentum() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let velocity = [
            [0.003, -0.002, 0.001],
            [0.007, -0.002, 0.001],
            [0.003, -0.005, 0.001],
            [0.003, -0.002, 0.003],
        ];
        let mut motion = BrittleTetMotion::try_from_solid(&solid, &velocity, 0.0, 48, cx).unwrap();
        let expected_energy = quadrature_kinetic(&velocity);
        let expected_momentum = momentum(&velocity);
        near(motion.kinetic_energy_j(), expected_energy, 2e-18);
        let mut peak_storage = 0.0_f64;
        for index in 0..48 {
            let step = motion
                .advance_midpoint(0.01, &[], &[], settings(), controls(), cx)
                .unwrap();
            assert_eq!(step.energy.inertia.numerical_loss_j, 0.0);
            assert_eq!(step.energy.fracture_dissipation_j, 0.0);
            assert_eq!(step.energy.external_mean_work_j(), 0.0);
            near(
                step.motion.kinetic_energy_j(),
                quadrature_kinetic(step.motion.velocity_mps()),
                2e-17,
            );
            near(
                step.motion.kinetic_energy_j() + step.motion.solid().stored_energy_j(),
                expected_energy,
                2e-12,
            );
            for (actual, expected) in momentum(step.motion.velocity_mps())
                .into_iter()
                .zip(expected_momentum)
            {
                near(actual, expected, 2e-12);
            }
            peak_storage = peak_storage.max(step.motion.solid().stored_energy_j());
            assert_eq!(
                step.motion.solid().solution().element_states[0].damage(),
                0.0
            );
            assert_eq!(step.motion.accepted_steps(), index + 1);
            near(step.motion.reference_mass_kg(), DENSITY * VOLUME, 2e-15);
            assert!(core::ptr::eq(
                step.motion.solid().materials().as_ptr(),
                sources.as_ptr()
            ));
            motion = step.motion;
        }
        assert!(
            peak_storage > 0.01 * expected_energy,
            "the test must exercise strain, not just rigid travel"
        );
        near(motion.time_s(), 0.48, 2e-15);
    });
}

#[test]
fn sourced_rupture_unloading_and_reclosure_keep_gc_and_the_complete_work_balance() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let mut motion =
            BrittleTetMotion::try_from_solid(&solid, &[[0.; 3]; 4], 0.0, 8, cx).unwrap();
        let mut external = 0.0;
        let mut damage = 0.0;
        let mut loss = 0.0;
        for extension in [0.0003, 0.0025, 0.0007, 0.0045, 0.009, -0.001, 0.0004] {
            let step = motion
                .advance_midpoint(0.1, &prescribed(extension), &[], settings(), controls(), cx)
                .unwrap();
            external += step.energy.external_mean_work_j();
            loss += step.energy.fracture_dissipation_j;
            let point = &step.motion.solid().solution().element_states[0];
            assert!(point.damage() >= damage);
            damage = point.damage();
            near(step.motion.solid().fracture_dissipation_j(), loss, 2e-13);
            near(
                external,
                step.motion.kinetic_energy_j() + step.motion.solid().stored_energy_j() + loss,
                2e-11,
            );
            assert_eq!(step.energy.inertia.numerical_loss_j, 0.0);
            near(
                step.energy.mean_internal_work_j,
                step.energy.material_work_j,
                2e-13,
            );
            if extension == 0.009 {
                assert!(point.failed());
                near(
                    step.motion.solid().solution().element_stress[0][0],
                    0.0,
                    1e-13,
                );
            }
            if extension < 0.0 {
                assert!(point.failed());
                near(
                    step.motion.solid().solution().element_stress[0][0],
                    YOUNG * extension,
                    2e-12,
                );
                assert_eq!(step.energy.fracture_dissipation_j, 0.0);
            }
            motion = step.motion;
        }
        // Footprint is 1/2 m², width is 1/3 m: V*Gc/h = Gc*A.
        near(loss, 0.5 * GC, 2e-13);
        assert!(motion.solid().solution().element_states[0].failed());
        near(motion.solid().stored_energy_j(), 0.0, 2e-13);
        near(motion.solid().solution().element_stress[0][0], 0.0, 1e-13);
        assert_eq!(solid.fracture_dissipation_j(), 0.0);
    });
}

#[test]
fn prescribed_work_uses_mean_reactions_while_the_solid_retains_endpoint_stress() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let motion = BrittleTetMotion::try_from_solid(&solid, &[[0.; 3]; 4], 0.0, 2, cx).unwrap();
        let extension = 0.0002;
        let dt = 0.01;
        let step = motion
            .advance_midpoint(dt, &prescribed(extension), &[], settings(), controls(), cx)
            .unwrap();
        let elastic = 0.5 * YOUNG * extension * extension * VOLUME;
        let velocity = 2.0 * extension / dt;
        let kinetic = 0.5 * (DENSITY * VOLUME / 10.0) * velocity * velocity;
        near(
            step.mean_element_stress_pa[0][0],
            0.5 * YOUNG * extension,
            2e-13,
        );
        near(
            step.motion.solid().solution().element_stress[0][0],
            YOUNG * extension,
            2e-13,
        );
        near(step.motion.velocity_mps()[1][0], velocity, 2e-15);
        near(step.energy.stored_after_j, elastic, 2e-16);
        near(step.energy.inertia.kinetic_after_j, kinetic, 2e-16);
        near(step.energy.support_mean_work_j, elastic + kinetic, 2e-16);
        near(step.energy.mean_internal_work_j, elastic, 2e-16);
        near(step.energy.inertia.inertial_work_j, kinetic, 2e-16);
        assert_eq!(step.energy.dead_load_mean_work_j, 0.0);
        assert_eq!(step.energy.inertia.numerical_loss_j, 0.0);
    });
}

#[test]
fn analytic_only_and_explicit_fallback_keep_nonsmooth_and_domain_refusals_distinct() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let motion = BrittleTetMotion::try_from_solid(&solid, &[[0.; 3]; 4], 0.0, 4, cx).unwrap();
        let saved = motion.solid().solution().clone();
        let mut strict = controls();
        strict.difference_step_strain = None;
        assert!(
            motion
                .advance_midpoint(0.1, &prescribed(0.0002), &[], settings(), strict, cx)
                .is_err()
        );
        let mut oversized = controls();
        oversized.difference_step_strain = Some(0.1);
        assert!(
            motion
                .advance_midpoint(0.1, &prescribed(0.0002), &[], settings(), oversized, cx)
                .is_err()
        );
        assert_eq!(motion.solid().solution(), &saved);
        assert_eq!(motion.accepted_steps(), 0);
        let first = motion
            .advance_midpoint(0.1, &prescribed(0.0002), &[], settings(), controls(), cx)
            .unwrap();
        assert!(first.numerical_tangent_evaluations > 0);
        // A positive intact held jump has the native smooth mean derivative.
        let analytic = first
            .motion
            .advance_midpoint(0.1, &prescribed(0.0004), &[], settings(), strict, cx)
            .unwrap();
        let enabled = first
            .motion
            .advance_midpoint(0.1, &prescribed(0.0004), &[], settings(), controls(), cx)
            .unwrap();
        assert_eq!(analytic.numerical_tangent_evaluations, 0);
        assert_eq!(enabled.numerical_tangent_evaluations, 0);
        assert_eq!(
            analytic.motion.solid().solution(),
            enabled.motion.solid().solution()
        );
        assert_eq!(analytic.energy, enabled.energy);
    });
}

#[test]
fn cancellation_and_failed_newton_candidates_leave_motion_and_step_budget_retryable() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let sources = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &sources, budget(), cx).unwrap();
        let motion = BrittleTetMotion::try_from_solid(&solid, &[[0.; 3]; 4], 1.0, 1, cx).unwrap();
        let saved = motion.solid().solution().clone();
        let mut too_short = settings();
        too_short.max_newton_iterations = 1;
        assert!(
            motion
                .advance_midpoint(0.1, &prescribed(0.0002), &[], too_short, controls(), cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            motion.advance_midpoint(0.1, &prescribed(0.0002), &[], settings(), controls(), cx),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::Cancelled))
        ));
        assert_eq!(motion.solid().solution(), &saved);
        assert_eq!(motion.velocity_mps(), &[[0.; 3]; 4]);
        assert_eq!(motion.time_s(), 1.0);
        assert_eq!(motion.accepted_steps(), 0);
        let retry_gate = CancelGate::new_clock_free();
        let retry = Cx::new(
            &retry_gate,
            cx.arena(),
            cx.stream_key(),
            cx.budget(),
            cx.mode(),
        );
        let accepted = motion
            .advance_midpoint(
                0.1,
                &prescribed(0.0002),
                &[],
                settings(),
                controls(),
                &retry,
            )
            .unwrap();
        assert_eq!(accepted.motion.time_s(), 1.1);
        assert_eq!(accepted.motion.accepted_steps(), 1);
        assert_eq!(accepted.motion.remaining_steps(), 0);
        assert!(matches!(
            accepted.motion.advance_midpoint(
                0.1,
                &prescribed(0.0003),
                &[],
                settings(),
                controls(),
                &retry
            ),
            Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: 1,
                maximum_steps: 1
            })
        ));
        assert_eq!(motion.solid().solution(), &saved);
    });
}
