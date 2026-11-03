//! Native bending work enters the original finite caloric cells exactly once.
use super::*;
use fs_couple::rod_motion::enthalpy::{
    RodAngularViscosity, RodBendingViscosity, RodTorsionalViscosity,
};

fn viscosity(value_pa_s: f64) -> ViscositySpec {
    ViscositySpec {
        value_pa_s,
        slope_pa_s_k: 0.0,
        lower_k: 200.0,
        upper_k: 500.0,
        upper_rad_s: f64::MAX,
    }
}
fn specimen(eta: f64) -> Spec {
    Spec {
        young_slope: 0.0,
        shear_slope: 0.0,
        alpha: 0.0,
        cp: 1e-4,
        bending_viscosity: Some(viscosity(eta)),
        ..Spec::default()
    }
}
fn pinned(
    source: &RodMaterialSource<'_>,
    temperatures: &[f64],
    angle: [f64; 3],
    speed: [f64; 3],
) -> RodMotion {
    let gate = CancelGate::new();
    let rod = source
        .resolve_at(temperatures, &gate)
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(LENGTH), Some(0.0), Some(0.0)],
                frame: None,
            }],
            &gate,
        )
        .unwrap();
    let model = Arc::new(
        rod.into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.frames[1] = So3::exp(So3Tangent::new(Vec3::new(angle[0], angle[1], angle[2]))).unwrap();
    RodMotion::new(model, pose, &[[0.0; 3]; 2], &[[0.0; 3], speed], 0.0, &gate).unwrap()
}
fn bending<'a>(selections: &'a [MaterialPropertySelection]) -> RodAngularViscosity<'a> {
    RodAngularViscosity {
        torsional: None,
        bending: Some(RodBendingViscosity {
            property: "bending_viscosity",
            selections,
        }),
    }
}
fn heat_solve() -> RodThermoelasticSolve {
    let mut config = solve();
    config.cell_energy_tolerance_j = 2e-14;
    config.thermal_energy_tolerance_j = 2e-14;
    config.coupling_energy_tolerance_j = 5e-13;
    config.newton.absolute_tolerance = 1e-15;
    config
}

// An independent pinned-section oracle: differentiate the finite-rotation
// shear/extension energy, then compose scalar dashpot and Verlet flows.
fn scalar_step(q: f64, w: f64, eta: f64, dt: f64) -> (f64, f64, f64) {
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let second = core::f64::consts::PI * RADIUS.powi(4) / 4.0;
    let inertia = 0.5 * RHO * LENGTH * second;
    let ea = YOUNG * area;
    let ga = 0.8 * SHEAR * area;
    let kb = YOUNG * second / LENGTH;
    let restoring = |q: f64| {
        let (s, c) = (0.5 * q).sin_cos();
        kb * q + 0.5 * LENGTH * s * (ga * c - ea * (c - 1.0))
    };
    let decay = (-0.5 * dt * eta * second / (LENGTH * inertia)).exp();
    let middle = decay * w - 0.5 * dt * restoring(q) / inertia;
    let next_q = q + dt * middle;
    let before_loss = middle - 0.5 * dt * restoring(next_q) / inertia;
    let heat = 0.5 * inertia * (1.0 - decay * decay) * (w * w + before_loss * before_loss);
    (next_q, decay * before_loss, heat)
}
fn angle(motion: &RodMotion, axis: usize) -> f64 {
    let q = motion.pose().frames[1]
        .body_minus(So3::identity())
        .unwrap()
        .angular;
    [q.x, q.y, q.z][axis]
}

#[test]
fn g1_both_bending_planes_deposit_annular_source_heat_and_close_entropy() {
    let specs = [specimen(60.0), specimen(4.0)];
    for axis in [1, 2] {
        fixture(
            &specs,
            1,
            &[300.0; 2],
            &[0.0; 2],
            false,
            |source, mut thermal, network| {
                let gate = CancelGate::new();
                let dt = 0.001;
                let (q0, w0) = (0.008, 0.6);
                let mut initial_angle = [0.0; 3];
                initial_angle[axis] = q0;
                let mut speed = [0.0; 3];
                speed[axis] = w0;
                let motion = pinned(source, &[300.0; 2], initial_angle, speed);
                let initial_caloric = total_heat(&thermal);
                let selections = vec![MaterialPropertySelection::SingleClaimOnly; 2];
                let mut owner = RodThermoelasticSystem::try_new(
                    motion,
                    source,
                    &thermal,
                    network,
                    &destinations(2, 1),
                    1e-11,
                    &gate,
                )
                .unwrap();
                let report = owner
                    .step_angular_viscous(
                        &mut thermal,
                        dt,
                        &[RodNodalWrench::default(); 2],
                        bending(&selections),
                        controls(dt),
                        heat_solve(),
                        &[0.0; 2],
                        &gate,
                    )
                    .unwrap();
                let (q1, w1, heat) = scalar_step(q0, w0, 7.5, dt);
                near(angle(owner.mechanical(), axis), q1, 3e-13);
                near(
                    owner.mechanical().angular_velocities_body_rad_s().unwrap()[1][axis],
                    w1,
                    3e-12,
                );
                let heating = report.bending.as_ref().unwrap();
                near(heating.mechanical_heat_j, heat, 3e-17);
                near(heating.deposited_heat_j, heat, 3e-17);
                near(
                    report
                        .mechanical
                        .bending_heat
                        .as_ref()
                        .unwrap()
                        .total_heat_j,
                    heat,
                    3e-17,
                );
                assert_eq!(heating.layers.len(), 2);
                assert_eq!(heating.source_sections.len(), 1);
                assert!(report.torsional.is_none());
                assert!(report.mechanical.torsional_heat.is_none());
                // The two annular second moments are 1:15. Eta=60:4 therefore
                // produces equal heat despite unequal 1:3 finite layer capacities.
                let mut entropy = 0.0;
                let mut generation = 0.0;
                for layer in 0..2 {
                    let (area, second) = geometry(layer, 2);
                    let capacity = RHO * LENGTH * area * specs[layer].cp;
                    let expected_t = 300.0 + 0.5 * heat / capacity;
                    near(
                        thermal.temperature_k(CellId(layer as u64)).unwrap(),
                        expected_t,
                        3e-8,
                    );
                    near(heating.layers[layer].second_moment_m4, second, 2e-20);
                    near(heating.layers[layer].heat_j, 0.5 * heat, 3e-17);
                    near(report.layers[layer].reversible_heat_j, 0.0, 3e-17);
                    entropy += capacity * (expected_t / 300.0).ln();
                    generation += 0.5 * heat / expected_t;
                }
                assert!(heating.deposited_heat_j > 0.0);
                near(total_heat(&thermal) - initial_caloric, heat, 3e-14);
                near(report.entropy.caloric_storage_change_j_k, entropy, 3e-15);
                near(report.entropy.bending_generation_j_k, generation, 3e-15);
                near(report.entropy.torsional_generation_j_k, 0.0, 0.0);
                near(report.entropy.accounting_residual_j_k, 0.0, 3e-15);
                near(report.coupling_energy_residual_j, 0.0, 4e-14);
                near(
                    report.total_energy_residual_j,
                    report.mechanical.integration_defect_j,
                    4e-14,
                );
            },
        );
    }
}

#[test]
fn g3_mixed_torsion_bending_sources_share_one_caloric_history_and_spatial_step() {
    let mut spec = specimen(80.0);
    spec.viscosity = Some(viscosity(12.0));
    spec.bending_viscosity.as_mut().unwrap().slope_pa_s_k = 10.0;
    spec.bending_viscosity.as_mut().unwrap().lower_k = 299.0;
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let dt = 0.001;
            let motion = pinned(source, &[300.0], [0.003, 0.005, -0.002], [0.6, 0.5, -0.3]);
            let mut spatial_thermal = thermal.clone();
            let mut owner = RodThermoelasticSystem::try_new(
                motion.clone(),
                source,
                &thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let mut spatial_owner = RodThermoelasticSystem::try_new(
                motion,
                source,
                &spatial_thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let selections = [MaterialPropertySelection::SingleClaimOnly];
            let losses = RodAngularViscosity {
                torsional: Some(RodTorsionalViscosity {
                    property: "shear_viscosity",
                    selections: &selections,
                }),
                ..bending(&selections)
            };
            let loads = [RodNodalWrench::default(); 2];
            let spatial = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
            let capacity = RHO * LENGTH * geometry(0, 1).0 * spec.cp;
            let mut temperature = 300.0;
            for step in 1..=3 {
                let before = total_heat(&thermal);
                let target = step as f64 * dt;
                let report = owner
                    .step_angular_viscous(
                        &mut thermal,
                        target,
                        &loads,
                        losses,
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let spatial_report = spatial_owner
                    .step_spatial_angular_viscous(
                        &mut spatial_thermal,
                        target,
                        &spatial,
                        losses,
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let t = report.torsional.as_ref().unwrap();
                let b = report.bending.as_ref().unwrap();
                near(
                    b.layers[0].viscosity_pa_s,
                    80.0 + 10.0 * (temperature - 300.0),
                    2e-10,
                );
                near(t.layers[0].viscosity_pa_s, 12.0, 0.0);
                assert!(t.deposited_heat_j > 0.0 && b.deposited_heat_j > 0.0);
                let heat = t.deposited_heat_j + b.deposited_heat_j;
                temperature += heat / capacity;
                near(thermal.temperature_k(CellId(0)).unwrap(), temperature, 3e-8);
                near(total_heat(&thermal) - before, heat, 3e-14);
                near(
                    report.entropy.bending_generation_j_k,
                    b.deposited_heat_j / temperature,
                    3e-15,
                );
                near(
                    report.entropy.torsional_generation_j_k,
                    t.deposited_heat_j / temperature,
                    3e-15,
                );
                near(report.entropy.accounting_residual_j_k, 0.0, 3e-15);
                assert_eq!(owner.mechanical().pose(), spatial_owner.mechanical().pose());
                assert_eq!(
                    owner.mechanical().body_spins_kg_m2_s(),
                    spatial_owner.mechanical().body_spins_kg_m2_s()
                );
                assert_eq!(thermal.cells(), spatial_thermal.cells());
                assert_eq!(report.bending, spatial_report.bending);
                assert_eq!(report.torsional, spatial_report.torsional);
                assert_eq!(owner.mechanical().accepted_steps(), step);
            }
        },
    );
}

#[test]
fn g0_explicit_zero_bending_preserves_original_torsional_history() {
    let mut spec = specimen(0.0);
    spec.viscosity = Some(viscosity(12.0));
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let dt = 0.001;
            let motion = pinned(source, &[300.0], [0.003, 0.005, -0.002], [0.6, 0.5, -0.3]);
            let mut legacy_thermal = thermal.clone();
            let mut owner = RodThermoelasticSystem::try_new(
                motion.clone(),
                source,
                &thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let mut legacy = RodThermoelasticSystem::try_new(
                motion,
                source,
                &legacy_thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let selections = [MaterialPropertySelection::SingleClaimOnly];
            let torsional = RodTorsionalViscosity {
                property: "shear_viscosity",
                selections: &selections,
            };
            let losses = RodAngularViscosity {
                torsional: Some(torsional),
                ..bending(&selections)
            };
            let loads = [RodNodalWrench::default(); 2];
            let actual = owner
                .step_angular_viscous(
                    &mut thermal,
                    dt,
                    &loads,
                    losses,
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            let expected = legacy
                .step_torsional(
                    &mut legacy_thermal,
                    dt,
                    &loads,
                    torsional,
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            assert_eq!(owner.mechanical().pose(), legacy.mechanical().pose());
            assert_eq!(
                owner.mechanical().body_spins_kg_m2_s(),
                legacy.mechanical().body_spins_kg_m2_s()
            );
            assert_eq!(thermal.cells(), legacy_thermal.cells());
            assert_eq!(actual.torsional, expected.torsional);
            assert_eq!(
                actual.total_energy_residual_j,
                expected.total_energy_residual_j
            );
            assert_eq!(actual.entropy, expected.entropy);
            assert_eq!(actual.bending.as_ref().unwrap().deposited_heat_j, 0.0);
            assert_eq!(actual.bending.as_ref().unwrap().mechanical_heat_j, 0.0);
            assert_eq!(actual.entropy.bending_generation_j_k, 0.0);
            assert!(actual.torsional.as_ref().unwrap().deposited_heat_j > 0.0);
        },
    );
}

#[test]
fn g4_bending_source_endpoint_refusal_and_cancellation_publish_neither_history() {
    let mut spec = specimen(100.0);
    spec.bending_viscosity.as_mut().unwrap().upper_k = 300.001;
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let motion = pinned(source, &[300.0], [0.0; 3], [0.0, 0.6, 0.0]);
            let before_thermal = thermal.clone();
            let mut owner = RodThermoelasticSystem::try_new(
                motion.clone(),
                source,
                &thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let before_layers = owner.layers().to_vec();
            let selections = [MaterialPropertySelection::SingleClaimOnly];
            let loads = [RodNodalWrench::default(); 2];
            let dt = 0.001;
            let heat = scalar_step(0.0, 0.6, 100.0, dt).2;
            let capacity = RHO * LENGTH * geometry(0, 1).0 * spec.cp;
            assert!(300.0 + heat / capacity > 300.001);
            assert!(source.resolve_at(&[300.0 + heat / capacity], &gate).is_ok());
            assert!(
                owner
                    .step_angular_viscous(
                        &mut thermal,
                        dt,
                        &loads,
                        bending(&selections),
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate
                    )
                    .is_err()
            );
            let cancelled = CancelGate::new();
            cancelled.request();
            assert!(
                owner
                    .step_angular_viscous(
                        &mut thermal,
                        1e-6,
                        &loads,
                        bending(&selections),
                        controls(1e-6),
                        heat_solve(),
                        &[0.0],
                        &cancelled
                    )
                    .is_err()
            );
            assert_eq!(thermal.cells(), before_thermal.cells());
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(owner.mechanical().pose(), motion.pose());
            assert_eq!(
                owner.mechanical().body_spins_kg_m2_s(),
                motion.body_spins_kg_m2_s()
            );
            assert_eq!(owner.layers(), before_layers);
            assert_eq!(owner.mechanical().accepted_steps(), 0);
            let mut fresh_thermal = before_thermal;
            let mut fresh = RodThermoelasticSystem::try_new(
                motion,
                source,
                &fresh_thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let dt = 1e-6;
            let retry = owner
                .step_angular_viscous(
                    &mut thermal,
                    dt,
                    &loads,
                    bending(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            let direct = fresh
                .step_angular_viscous(
                    &mut fresh_thermal,
                    dt,
                    &loads,
                    bending(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            assert!(retry.bending.as_ref().unwrap().deposited_heat_j > 0.0);
            assert!(thermal.temperature_k(CellId(0)).unwrap() < 300.001);
            assert_eq!(thermal.cells(), fresh_thermal.cells());
            assert_eq!(owner.mechanical().pose(), fresh.mechanical().pose());
            assert_eq!(
                owner.mechanical().body_spins_kg_m2_s(),
                fresh.mechanical().body_spins_kg_m2_s()
            );
            assert_eq!(
                retry.total_energy_residual_j,
                direct.total_energy_residual_j
            );
        },
    );
}
