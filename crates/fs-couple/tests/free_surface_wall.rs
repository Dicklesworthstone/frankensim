//! Original population exchange, native rigid inertia and independent SI oracles.
//! Synthetic liquid claims establish numerical coupling, not measured validity.

use fs_couple::free_surface_wall::{GuidedFreeSurfaceWall, GuidedWallControls, GuidedWallError};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::{ContactModel3, PlanarWall3};
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialFreeSurfaceDomain3, plan_free_surface_from_material,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, resolve_liquid_surface_state,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection};
use fs_mbd::{MassProperties, Pose, RigidBodyState, Vec3};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

const DX: f64 = 1e-3;
const DT: f64 = 2e-4;
const CELLS: usize = 8 * 8 * 8;
const AREA: f64 = 64.0 * DX * DX;
const CELL_MASS: f64 = 1000.0 * DX * DX * DX;
// At rest, summing the original D3Q19 bounce-back weights over a full plane
// gives J_liquid,t = rho dx^3 N_face w / 3 for one native timestep.
const WALL_IMPEDANCE_KG: f64 = CELL_MASS * 64.0 / 3.0;

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn source(mu: f64) -> ResolvedLiquidSurfaceState {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (DENSITY_PROPERTY, 1000.0, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, mu, DYNAMIC_VISCOSITY_DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        290.0,
                        310.0,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        90_000.0,
                        110_000.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic finite-inertia wall fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall-driven liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap();
    resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [8, 8, 8],
        gravity_m_s2: [0.0; 3],
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: CELLS,
    }
}
fn liquid(source: &ResolvedLiquidSurfaceState, pressure: Option<f64>) -> MaterialFreeSurface3<'_> {
    let plan = plan_free_surface_from_material(source, DX, DT, 0.1).unwrap();
    let plan = match pressure {
        Some(value) => plan.with_pressure_reference(value).unwrap(),
        None => plan,
    };
    plan.free_surface3(domain(), |_, _, z| z < 4).unwrap()
}
fn properties(mass: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(mass, mass, mass)).unwrap()
}
fn body(mass: f64, speed: f64, tangent: Vec3) -> RigidBodyState {
    RigidBodyState::new(Pose::identity(), tangent.scale(mass * speed), Vec3::ZERO).unwrap()
}
fn controls() -> GuidedWallControls {
    GuidedWallControls {
        midpoint_speed_bracket_m_s: [-0.05, 0.05],
        maximum_native_trials: 24,
        maximum_trial_cells: 24 * CELLS,
        maximum_relative_mass_defect: 1e-12,
        momentum_absolute_kg_m_s: 1e-15,
        midpoint_absolute_m_s: 1e-10,
        work_absolute_j: 1e-16,
        relative_tolerance: 1e-11,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn speed(motion: &GuidedFreeSurfaceWall<'_>) -> f64 {
    motion
        .body()
        .center_of_mass_velocity_world(motion.properties())
        .unwrap()
        .dot(motion.tangent_world())
}
fn liquid_momentum(liquid: &MaterialFreeSurface3<'_>) -> Vec3 {
    let mut result = Vec3::ZERO;
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if let Some(s) = liquid.sample_at_si(x, y, z).unwrap() {
                    let m = s.fill_fraction * s.density_kg_m3 * DX * DX * DX;
                    result = result.add(
                        Vec3::new(s.velocity_m_s[0], s.velocity_m_s[1], s.velocity_m_s[2]).scale(m),
                    );
                }
            }
        }
    }
    result
}

#[test]
fn actual_finite_wall_and_liquid_match_independent_reflection_and_momentum_oracles() {
    let source = source(0.1);
    let gate = CancelGate::new_clock_free();
    for tangent in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.6, 0.8, 0.0)] {
        for mass in [
            0.5 * WALL_IMPEDANCE_KG,
            2.0 * WALL_IMPEDANCE_KG,
            1000.0 * WALL_IMPEDANCE_KG,
        ] {
            let initial_speed = 0.03;
            let fluid = liquid(&source, None);
            let before_p = liquid_momentum(&fluid);
            let original_bits = fluid.lattice().grid().state_bits();
            assert_eq!(fluid.lattice().domain_cells(), [8, 8, 8]);
            let old_body = body(mass, initial_speed, tangent);
            let motion = GuidedFreeSurfaceWall::new(
                fluid,
                properties(mass),
                old_body,
                PlanarWall3::ZMin,
                tangent,
                CELLS,
                &gate,
            )
            .unwrap();
            let step = motion.step_candidate(0.0, controls(), &gate).unwrap();
            let midpoint = 2.0 * mass * initial_speed / (2.0 * mass + WALL_IMPEDANCE_KG);
            let final_speed = 2.0 * midpoint - initial_speed;
            near(speed(&step.motion), final_speed, 2e-11);
            near(step.midpoint_velocity_world.dot(tangent), midpoint, 1e-11);
            let expected_impulse = WALL_IMPEDANCE_KG * midpoint;
            let after_p = liquid_momentum(step.motion.liquid());
            near(after_p.sub(before_p).dot(tangent), expected_impulse, 2e-13);
            near(
                after_p.sub(before_p).dot(tangent)
                    + step
                        .motion
                        .body()
                        .linear_momentum_world()
                        .sub(old_body.linear_momentum_world())
                        .dot(tangent),
                0.0,
                2e-13,
            );
            near(
                step.wall.fluid_work_j,
                WALL_IMPEDANCE_KG * midpoint * midpoint,
                1e-16,
            );
            near(
                step.momentum.work.kinetic_energy_change_j,
                0.5 * mass * (final_speed * final_speed - initial_speed * initial_speed),
                1e-16,
            );
            near(
                step.momentum.work.kinetic_energy_change_j + step.wall.fluid_work_j,
                0.0,
                1e-16,
            );
            near(
                step.motion.body().pose().position_world().dot(tangent),
                midpoint * DT,
                1e-14,
            );
            near(step.mechanical_work_residual_j, 0.0, 1e-16);
            assert!(
                step.native_trials <= 4,
                "affine native response should resolve promptly"
            );
            assert_eq!(step.trial_cells, step.native_trials * CELLS);
            assert_eq!(step.motion.liquid().lattice().grid().steps(), 1);
            assert_eq!(step.motion.accepted_steps(), 1);
            near(step.motion.time_s(), DT, 0.0);
            assert_eq!(step.motion.liquid().scaling().material(), &source);
            assert!(step.absolute_guide_impulse_world.is_none());
            assert_eq!(motion.body(), old_body);
            assert_eq!(motion.liquid().lattice().grid().state_bits(), original_bits);
        }
    }
}

#[test]
fn accepted_midpoint_replays_native_populations_and_continues_the_physical_history() {
    let source = source(0.1);
    let gate = CancelGate::new_clock_free();
    let tangent = Vec3::new(1.0, 0.0, 0.0);
    let mass = 2.0 * WALL_IMPEDANCE_KG;
    let mut motion = GuidedFreeSurfaceWall::new(
        liquid(&source, None),
        properties(mass),
        body(mass, 0.02, tangent),
        PlanarWall3::ZMin,
        tangent,
        CELLS,
        &gate,
    )
    .unwrap();
    let mut wall_work = 0.0;
    let mut applied_work = 0.0;
    let initial_energy = 0.5 * mass * 0.02_f64.powi(2);
    for i in 0..4 {
        let force = if i < 2 { 0.001 } else { -0.001 };
        let old_speed = speed(&motion);
        let step = motion.step_candidate(force, controls(), &gate).unwrap();
        let mut replay = motion.liquid().clone();
        let expected = replay
            .step_with_planar_wall_motion_si(
                PlanarWall3::ZMin,
                step.wall.velocity_m_s,
                controls().maximum_relative_mass_defect,
                &gate,
            )
            .unwrap();
        assert_eq!(step.wall, expected);
        assert_eq!(
            step.motion.liquid().lattice().grid().state_bits(),
            replay.lattice().grid().state_bits()
        );
        let mean = 0.5 * (old_speed + speed(&step.motion));
        near(step.external_work_j, force * DT * mean, 1e-17);
        near(
            step.momentum.work.kinetic_energy_change_j + step.wall.fluid_work_j
                - step.external_work_j
                - step.raw_guide_work_j,
            0.0,
            2e-16,
        );
        wall_work += step.wall.fluid_work_j;
        applied_work += step.external_work_j;
        motion = step.motion;
        assert_eq!(motion.accepted_steps(), i + 1);
        assert_eq!(motion.liquid().lattice().grid().steps(), i + 1);
    }
    near(
        0.5 * mass * speed(&motion).powi(2) - initial_energy + wall_work - applied_work,
        0.0,
        5e-16,
    );
    near(motion.time_s(), 4.0 * DT, 1e-18);
}

#[test]
fn numerical_and_absolute_pressure_guide_reactions_are_separate_and_do_no_tangent_work() {
    let source = source(0.1);
    let gate = CancelGate::new_clock_free();
    let tangent = Vec3::new(1.0, 0.0, 0.0);
    let mass = 2.0 * WALL_IMPEDANCE_KG;
    let step = |pressure| {
        GuidedFreeSurfaceWall::new(
            liquid(&source, Some(pressure)),
            properties(mass),
            body(mass, 0.03, tangent),
            PlanarWall3::ZMin,
            tangent,
            CELLS,
            &gate,
        )
        .unwrap()
        .step_candidate(0.0, controls(), &gate)
        .unwrap()
    };
    let first = step(100_000.0);
    let shifted = step(101_000.0);
    assert_eq!(first.motion.body(), shifted.motion.body());
    assert_eq!(
        first.motion.liquid().lattice().grid().state_bits(),
        shifted.motion.liquid().lattice().grid().state_bits()
    );
    assert_eq!(
        first.raw_guide_impulse_world,
        shifted.raw_guide_impulse_world
    );
    near(
        first.pressure_perturbation_guide_impulse_world.z,
        0.0,
        1e-15,
    );
    near(
        first.absolute_guide_impulse_world.unwrap().z,
        100_000.0 * AREA * DT,
        1e-15,
    );
    near(
        shifted.absolute_guide_impulse_world.unwrap().z
            - first.absolute_guide_impulse_world.unwrap().z,
        1000.0 * AREA * DT,
        1e-15,
    );
    assert!(first.raw_guide_impulse_world.z > 0.0);
    assert_ne!(
        first.raw_guide_impulse_world.z,
        first.absolute_guide_impulse_world.unwrap().z
    );
    near(first.raw_guide_work_j, 0.0, 1e-18);
}

#[test]
fn geometry_work_trials_and_cancellation_refuse_without_publishing_any_owner() {
    let source = source(0.1);
    let gate = CancelGate::new_clock_free();
    let tangent = Vec3::new(1.0, 0.0, 0.0);
    let mass = WALL_IMPEDANCE_KG;
    let fluid = liquid(&source, None);
    let old_body = body(mass, 0.03, tangent);
    assert!(
        GuidedFreeSurfaceWall::new(
            fluid.clone(),
            properties(mass),
            old_body,
            PlanarWall3::ZMin,
            tangent,
            CELLS - 1,
            &gate
        )
        .is_err()
    );
    assert!(
        GuidedFreeSurfaceWall::new(
            fluid.clone(),
            properties(mass),
            old_body,
            PlanarWall3::XMin,
            tangent,
            CELLS,
            &gate
        )
        .is_err()
    );
    let spinning = RigidBodyState::new(
        Pose::identity(),
        old_body.linear_momentum_world(),
        Vec3::new(0.0, 0.0, 1e-5),
    )
    .unwrap();
    assert!(
        GuidedFreeSurfaceWall::new(
            fluid.clone(),
            properties(mass),
            spinning,
            PlanarWall3::ZMin,
            tangent,
            CELLS,
            &gate
        )
        .is_err()
    );
    let motion = GuidedFreeSurfaceWall::new(
        fluid,
        properties(mass),
        old_body,
        PlanarWall3::ZMin,
        tangent,
        CELLS,
        &gate,
    )
    .unwrap();
    let bits = motion.liquid().lattice().grid().state_bits();
    let mut limited = controls();
    limited.maximum_native_trials = 1;
    assert!(matches!(
        motion.step_candidate(0.0, limited, &gate),
        Err(GuidedWallError::TrialBudget {
            completed: 1,
            maximum: 1
        })
    ));
    limited = controls();
    limited.maximum_trial_cells = CELLS;
    assert!(matches!(motion.step_candidate(0.0, limited, &gate),
        Err(GuidedWallError::WorkBudget { required_trial_cells, maximum })
            if required_trial_cells == 2 * CELLS && maximum == CELLS));
    limited = controls();
    limited.midpoint_speed_bracket_m_s = [0.04, 0.05];
    assert!(matches!(
        motion.step_candidate(0.0, limited, &gate),
        Err(GuidedWallError::Unbracketed { .. })
    ));
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        motion.step_candidate(0.0, controls(), &stopped),
        Err(GuidedWallError::Cancelled)
    ));
    assert_eq!(motion.body(), old_body);
    assert_eq!(motion.liquid().lattice().grid().state_bits(), bits);
    assert_eq!(motion.accepted_steps(), 0);
    let first = motion.step_candidate(0.0, controls(), &gate).unwrap();
    let retry = motion.step_candidate(0.0, controls(), &gate).unwrap();
    assert_eq!(first.motion.body(), retry.motion.body());
    assert_eq!(
        first.motion.liquid().lattice().grid().state_bits(),
        retry.motion.liquid().lattice().grid().state_bits()
    );
}

#[test]
fn native_late_liquid_speed_refusal_preserves_wall_and_liquid_together() {
    let source = source(0.1);
    let gate = CancelGate::new_clock_free();
    let fluid = plan_free_surface_from_material(&source, DX, DT, 0.02)
        .unwrap()
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                gravity_m_s2: [0.0, 0.0, -100.0],
                ..domain()
            },
            |_, _, z| z < 4,
        )
        .unwrap();
    let tangent = Vec3::new(1.0, 0.0, 0.0);
    let mass = WALL_IMPEDANCE_KG;
    let motion = GuidedFreeSurfaceWall::new(
        fluid,
        properties(mass),
        body(mass, 0.0, tangent),
        PlanarWall3::ZMin,
        tangent,
        CELLS,
        &gate,
    )
    .unwrap();
    let bits = motion.liquid().lattice().grid().state_bits();
    let old_body = motion.body();
    let mut c = controls();
    c.midpoint_speed_bracket_m_s = [-0.005, 0.005];
    assert!(matches!(
        motion.step_candidate(0.0, c, &gate),
        Err(GuidedWallError::Fluid(_))
    ));
    assert_eq!(motion.body(), old_body);
    assert_eq!(motion.liquid().lattice().grid().state_bits(), bits);
    assert_eq!(motion.liquid().lattice().grid().steps(), 0);
    assert_eq!(motion.time_s(), 0.0);
}
