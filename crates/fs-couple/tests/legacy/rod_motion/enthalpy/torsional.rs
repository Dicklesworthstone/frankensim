//! Actual finite-cell shaft heating, with independent scalar and annular oracles.
use super::*;
use fs_couple::rod_motion::enthalpy::RodTorsionalViscosity;

fn specimen(eta: f64) -> Spec {
    Spec {
        young_slope: 0.0,
        shear_slope: 0.0,
        alpha: 0.0,
        cp: 1e-4,
        viscosity: Some(ViscositySpec {
            value_pa_s: eta,
            slope_pa_s_k: 0.0,
            lower_k: 200.0,
            upper_k: 500.0,
            upper_rad_s: f64::MAX,
        }),
        ..Spec::default()
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
fn shaft(
    source: &RodMaterialSource<'_>,
    temperatures: &[f64],
    angle: f64,
    speed: f64,
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
    pose.frames[1] = So3::exp(So3Tangent::new(Vec3::new(angle, 0.0, 0.0))).unwrap();
    RodMotion::new(
        model,
        pose,
        &[[0.0; 3]; 2],
        &[[0.0; 3], [speed, 0.0, 0.0]],
        0.0,
        &gate,
    )
    .unwrap()
}
fn polar_moment() -> f64 {
    core::f64::consts::PI * RADIUS.powi(4) / 2.0
}
fn tip_inertia() -> f64 {
    0.5 * RHO * LENGTH * polar_moment()
}
fn caloric_capacity(spec: Spec, layer: usize, layers: usize) -> f64 {
    RHO * LENGTH * geometry(layer, layers).0 * spec.cp
}

// A pinned circular shaft has one exact twist coordinate. Independently compose
// its two scalar dashpot exponentials and elastic velocity-Verlet map. No native
// loss, section stiffness, heat receipt, or source resolver supplies this oracle.
fn scalar_step(q: f64, w: f64, shear: f64, eta: f64, dt: f64) -> (f64, f64, f64) {
    let inertia = tip_inertia();
    let stiffness = shear * polar_moment() / LENGTH;
    let resistance = eta * polar_moment() / LENGTH;
    let decay = (-0.5 * dt * resistance / inertia).exp();
    let middle = decay * w - 0.5 * dt * stiffness * q / inertia;
    let next_q = q + dt * middle;
    let before_final_loss = middle - 0.5 * dt * stiffness * next_q / inertia;
    let heat =
        0.5 * inertia * (1.0 - decay * decay) * (w * w + before_final_loss * before_final_loss);
    (next_q, decay * before_final_loss, heat)
}
fn angle(state: &RodMotion) -> f64 {
    state.pose().frames[1]
        .body_minus(So3::identity())
        .unwrap()
        .angular
        .x
}

#[test]
fn g1_equal_viscous_polar_weights_heat_unequal_layer_masses_and_close_entropy() {
    let specs = [specimen(60.0), specimen(4.0)];
    fixture(
        &specs,
        1,
        &[300.0; 2],
        &[0.0; 2],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let dt = 0.01;
            let initial_speed = 0.6;
            let motion = shaft(source, &[300.0; 2], 0.0, initial_speed);
            let initial_energy = 0.5 * tip_inertia() * initial_speed * initial_speed;
            let initial_caloric = total_heat(&thermal);
            let routes = destinations(2, 1);
            let selections = vec![MaterialPropertySelection::SingleClaimOnly; 2];
            let mut owner = RodThermoelasticSystem::try_new(
                motion, source, &thermal, network, &routes, 1e-11, &gate,
            )
            .unwrap();
            let report = owner
                .step_torsional(
                    &mut thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    RodTorsionalViscosity {
                        property: "shear_viscosity",
                        selections: &selections,
                    },
                    controls(dt),
                    heat_solve(),
                    &[0.0; 2],
                    &gate,
                )
                .unwrap();
            // At R/2, the annular polar moments are 1:15. With eta=60:4,
            // each layer generates half the heat, although their masses are 1:3.
            let (q, w, heat) = scalar_step(0.0, initial_speed, SHEAR, 7.5, dt);
            near(angle(owner.mechanical()), q, 2e-13);
            near(
                owner.mechanical().angular_velocities_body_rad_s().unwrap()[1][0],
                w,
                2e-12,
            );
            let heating = report.torsional.as_ref().unwrap();
            assert_eq!(heating.layers.len(), 2);
            assert_eq!(heating.source_sections.len(), 1);
            near(heating.mechanical_heat_j, heat, 2e-17);
            near(heating.deposited_heat_j, heat, 2e-17);
            near(heating.partition_residual_j, 0.0, 2e-20);
            near(
                report
                    .mechanical
                    .torsional_heat
                    .as_ref()
                    .unwrap()
                    .total_heat_j,
                heat,
                2e-17,
            );
            let mut expected_entropy = 0.0;
            let mut expected_generation = 0.0;
            let mut rises = [0.0; 2];
            for layer in 0..2 {
                let capacity = caloric_capacity(specs[layer], layer, 2);
                let expected_t = 300.0 + 0.5 * heat / capacity;
                let actual_t = thermal.temperature_k(routes[layer].cell).unwrap();
                near(actual_t, expected_t, 2e-8);
                rises[layer] = actual_t - 300.0;
                let receipt = &heating.layers[layer];
                assert_eq!(receipt.destination, routes[layer]);
                near(receipt.polar_moment_m4, 2.0 * geometry(layer, 2).1, 2e-20);
                assert_eq!(
                    receipt.viscosity_pa_s,
                    specs[layer].viscosity.unwrap().value_pa_s
                );
                near(receipt.heat_j, 0.5 * heat, 2e-17);
                near(capacity * (actual_t - 300.0), 0.5 * heat, 2e-14);
                expected_entropy += capacity * (expected_t / 300.0).ln();
                expected_generation += 0.5 * heat / expected_t;
            }
            near(rises[0], 3.0 * rises[1], 5e-8);
            near(report.reversible_heat_j, 0.0, 2e-18);
            near(report.external_thermal_heat_j, 0.0, 2e-14);
            near(report.thermal.external_heat_j, heat, 2e-14);
            near(total_heat(&thermal) - initial_caloric, heat, 2e-14);
            let independent_total_change = 0.5 * tip_inertia() * w * w
                + 0.5 * SHEAR * polar_moment() / LENGTH * q * q
                - initial_energy
                + total_heat(&thermal)
                - initial_caloric;
            near(
                independent_total_change,
                report.mechanical.integration_defect_j,
                3e-14,
            );
            assert!(report.mechanical.integration_defect_j.abs() > 1e-12);
            near(
                report.total_energy_residual_j,
                report.mechanical.integration_defect_j,
                3e-14,
            );
            near(report.coupling_energy_residual_j, 0.0, 3e-14);
            near(
                report.entropy.caloric_storage_change_j_k,
                expected_entropy,
                2e-15,
            );
            near(
                report.entropy.torsional_generation_j_k,
                expected_generation,
                2e-15,
            );
            assert!(report.entropy.torsional_generation_j_k > 0.0);
            near(report.entropy.thermoelastic_temporal_defect_j_k, 0.0, 2e-18);
            near(report.entropy.external_source_at_cell_j_k, 0.0, 2e-18);
            assert!(!report.entropy.has_unclassified_external_sources);
            near(report.entropy.accounting_residual_j_k, 0.0, 2e-15);
        },
    );
}

#[test]
fn g1_temperature_dependent_viscosity_refreshes_each_interval_without_reclassifying_reversible_heat()
 {
    let mut spec = specimen(40.0);
    spec.shear_slope = -0.4;
    spec.viscosity.as_mut().unwrap().slope_pa_s_k = 20.0;
    spec.viscosity.as_mut().unwrap().lower_k = 299.0;
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let dt = 0.005;
            let (mut q, mut w, mut t) = (0.005, 0.65, 300.0);
            let motion = shaft(source, &[t], q, w);
            let selections = [MaterialPropertySelection::SingleClaimOnly];
            let capacity = caloric_capacity(spec, 0, 1);
            // For affine G(T) and alpha=0, Psi+TS has this T-independent
            // coefficient; caloric energy is still measured in its original cell.
            let internal_stiffness = (SHEAR - 300.0 * spec.shear_slope) * polar_moment() / LENGTH;
            let initial_total = 0.5 * tip_inertia() * w * w + 0.5 * internal_stiffness * q * q;
            let initial_caloric = total_heat(&thermal);
            let mut owner = RodThermoelasticSystem::try_new(
                motion,
                source,
                &thermal,
                network,
                &destinations(1, 1),
                1e-11,
                &gate,
            )
            .unwrap();
            let mut defects = 0.0;
            let mut loss_heat = 0.0;
            let mut reversible_heat = 0.0;
            let mut largest_thermoelastic_defect = 0.0_f64;
            for _ in 0..8 {
                let eta = 40.0 + 20.0 * (t - 300.0);
                let shear = SHEAR + spec.shear_slope * (t - 300.0);
                let (next_q, next_w, heat) = scalar_step(q, w, shear, eta, dt);
                let entropy_change =
                    -0.5 * spec.shear_slope * polar_moment() / LENGTH * (next_q * next_q - q * q);
                let reversible = -t * entropy_change;
                let next_t = t + (heat + reversible) / capacity;
                let target = owner.mechanical().time_s() + dt;
                let report = owner
                    .step_torsional(
                        &mut thermal,
                        target,
                        &[RodNodalWrench::default(); 2],
                        RodTorsionalViscosity {
                            property: "shear_viscosity",
                            selections: &selections,
                        },
                        controls(dt * 1.001),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let heating = report.torsional.as_ref().unwrap();
                near(heating.layers[0].viscosity_pa_s, eta, 2e-6);
                let receipt = &heating.source_sections[0].shear().elastic().layers()[0].material;
                near(
                    receipt
                        .query_point()
                        .iter()
                        .find(|(axis, _)| axis == "T")
                        .unwrap()
                        .1,
                    t,
                    1e-7,
                );
                near(
                    receipt.property("shear_viscosity").unwrap().value_si(),
                    eta,
                    2e-6,
                );
                near(heating.deposited_heat_j, heat, 3e-14);
                near(angle(owner.mechanical()), next_q, 2e-10);
                near(
                    owner.mechanical().angular_velocities_body_rad_s().unwrap()[1][0],
                    next_w,
                    2e-8,
                );
                near(thermal.temperature_k(CellId(0)).unwrap(), next_t, 1e-7);
                near(report.reversible_heat_j, reversible, 3e-14);
                near(report.thermal.external_heat_j, heat + reversible, 3e-14);
                near(report.external_thermal_heat_j, 0.0, 3e-14);
                near(
                    report.entropy.torsional_generation_j_k,
                    heat / next_t,
                    2e-15,
                );
                near(
                    report.entropy.mechanical_storage_change_j_k,
                    entropy_change,
                    2e-15,
                );
                near(
                    report.entropy.thermoelastic_temporal_defect_j_k,
                    entropy_change + reversible / next_t,
                    2e-15,
                );
                near(report.entropy.accounting_residual_j_k, 0.0, 2e-15);
                assert!(!report.entropy.has_unclassified_external_sources);
                near(report.coupling_energy_residual_j, 0.0, 3e-14);
                largest_thermoelastic_defect = largest_thermoelastic_defect
                    .max(report.entropy.thermoelastic_temporal_defect_j_k.abs());
                defects += report.mechanical.integration_defect_j;
                loss_heat += heat;
                reversible_heat += reversible;
                (q, w, t) = (next_q, next_w, next_t);
            }
            assert!(
                t > 300.5,
                "the viscosity refresh must materially change the next trajectory"
            );
            assert!(largest_thermoelastic_defect > 1e-16);
            assert!(reversible_heat.abs() > 1e-9);
            near(
                total_heat(&thermal) - initial_caloric,
                loss_heat + reversible_heat,
                2e-13,
            );
            let final_internal = 0.5 * tip_inertia() * w * w + 0.5 * internal_stiffness * q * q;
            near(
                final_internal - initial_total + total_heat(&thermal) - initial_caloric,
                defects,
                2e-13,
            );
        },
    );
}

#[test]
fn g0_prepared_body_torque_composes_with_heating_and_explicit_zero_loss_preserves_native_step() {
    for eta in [0.0, 40.0] {
        fixture(
            &[specimen(eta)],
            1,
            &[300.0],
            &[0.0],
            false,
            |source, mut thermal, network| {
                let gate = CancelGate::new();
                let dt = 0.003;
                let motion = shaft(source, &[300.0], 0.005, 0.4);
                let original = motion.clone();
                let original_thermal = thermal.clone();
                let mut spatial_thermal = thermal.clone();
                let routes = destinations(1, 1);
                let mut world = RodThermoelasticSystem::try_new(
                    motion.clone(),
                    source,
                    &thermal,
                    network,
                    &routes,
                    1e-11,
                    &gate,
                )
                .unwrap();
                let mut spatial = RodThermoelasticSystem::try_new(
                    motion,
                    source,
                    &spatial_thermal,
                    network,
                    &routes,
                    1e-11,
                    &gate,
                )
                .unwrap();
                let selections = [MaterialPropertySelection::SingleClaimOnly];
                let viscosity = RodTorsionalViscosity {
                    property: "shear_viscosity",
                    selections: &selections,
                };
                let torque = 3e-5;
                let loads = [
                    RodNodalWrench::default(),
                    RodNodalWrench {
                        force_n: [0.0; 3],
                        torque_nm: [torque, 0.0, 0.0],
                    },
                ];
                // Rotation about x preserves the x axis, so this actual body-held
                // couple is independently equivalent to the world torque above.
                let prepared = PreparedRodLoads::new(
                    &[RodNodalWrench::default(); 2],
                    &[[0.0; 3], [torque, 0.0, 0.0]],
                    &[],
                    &gate,
                )
                .unwrap();
                let a = world
                    .step_torsional(
                        &mut thermal,
                        dt,
                        &loads,
                        viscosity,
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let b = spatial
                    .step_spatial_torsional(
                        &mut spatial_thermal,
                        dt,
                        &prepared,
                        viscosity,
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                assert_eq!(world.mechanical().pose(), spatial.mechanical().pose());
                assert_eq!(
                    world.mechanical().body_spins_kg_m2_s(),
                    spatial.mechanical().body_spins_kg_m2_s()
                );
                assert_eq!(thermal.cells(), spatial_thermal.cells());
                assert_eq!(a.mechanical.torsional_heat, b.mechanical.torsional_heat);
                assert_eq!(a.total_energy_residual_j, b.total_energy_residual_j);
                assert_eq!(
                    a.entropy.torsional_generation_j_k,
                    b.entropy.torsional_generation_j_k
                );
                assert!(a.mechanical.source_work_quadrature_j > 0.0);
                near(a.external_thermal_heat_j, 0.0, 3e-14);
                if eta == 0.0 {
                    let mut plain_thermal = original_thermal;
                    let mut plain = RodThermoelasticSystem::try_new(
                        original,
                        source,
                        &plain_thermal,
                        network,
                        &routes,
                        1e-11,
                        &gate,
                    )
                    .unwrap();
                    let native = plain
                        .step(
                            &mut plain_thermal,
                            dt,
                            &loads,
                            controls(dt),
                            heat_solve(),
                            &[0.0],
                            &gate,
                        )
                        .unwrap();
                    assert_eq!(world.mechanical().pose(), plain.mechanical().pose());
                    assert_eq!(
                        world.mechanical().body_spins_kg_m2_s(),
                        plain.mechanical().body_spins_kg_m2_s()
                    );
                    assert_eq!(thermal.cells(), plain_thermal.cells());
                    assert_eq!(a.total_energy_residual_j, native.total_energy_residual_j);
                    assert_eq!(a.torsional.as_ref().unwrap().deposited_heat_j, 0.0);
                    assert_eq!(a.entropy.torsional_generation_j_k, 0.0);
                    assert!(native.torsional.is_none());
                    assert!(native.mechanical.torsional_heat.is_none());
                } else {
                    assert!(a.torsional.as_ref().unwrap().deposited_heat_j > 0.0);
                }
            },
        );
    }
}

#[test]
fn g4_late_viscosity_domain_refusal_preserves_histories_and_allows_valid_retry() {
    let mut spec = specimen(100.0);
    spec.viscosity.as_mut().unwrap().upper_k = 300.01;
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let original = shaft(source, &[300.0], 0.0, 0.5);
            let original_thermal = thermal.clone();
            let routes = destinations(1, 1);
            let mut owner = RodThermoelasticSystem::try_new(
                original.clone(),
                source,
                &thermal,
                network,
                &routes,
                1e-11,
                &gate,
            )
            .unwrap();
            let layers = owner.layers().to_vec();
            let selections = [MaterialPropertySelection::SingleClaimOnly];
            let viscosity = RodTorsionalViscosity {
                property: "shear_viscosity",
                selections: &selections,
            };
            let attempted = 0.01;
            let (_, _, heat) = scalar_step(0.0, 0.5, SHEAR, 100.0, attempted);
            let forbidden_t = 300.0 + heat / caloric_capacity(spec, 0, 1);
            assert!(forbidden_t > 300.01);
            // Mechanical E/G/alpha and the caloric chart support this endpoint.
            // Only the original viscosity source ends at 300.01 K.
            assert!(source.resolve_at(&[forbidden_t], &gate).is_ok());
            assert!(
                owner
                    .step_torsional(
                        &mut thermal,
                        attempted,
                        &[RodNodalWrench::default(); 2],
                        viscosity,
                        controls(attempted),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .is_err()
            );
            assert_eq!(thermal.cells(), original_thermal.cells());
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(owner.mechanical().pose(), original.pose());
            assert_eq!(
                owner.mechanical().body_spins_kg_m2_s(),
                original.body_spins_kg_m2_s()
            );
            assert_eq!(
                owner.mechanical().linear_momenta_kg_m_s(),
                original.linear_momenta_kg_m_s()
            );
            assert_eq!(owner.mechanical().time_s(), 0.0);
            assert_eq!(owner.mechanical().accepted_steps(), 0);
            assert_eq!(owner.layers(), layers);
            let dt = 1e-5;
            let mut fresh_thermal = original_thermal;
            let mut fresh = RodThermoelasticSystem::try_new(
                original,
                source,
                &fresh_thermal,
                network,
                &routes,
                1e-11,
                &gate,
            )
            .unwrap();
            let retry = owner
                .step_torsional(
                    &mut thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    viscosity,
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            let direct = fresh
                .step_torsional(
                    &mut fresh_thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    viscosity,
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            assert!(thermal.temperature_k(CellId(0)).unwrap() < 300.01);
            assert!(retry.torsional.as_ref().unwrap().deposited_heat_j > 0.0);
            assert_eq!(thermal.cells(), fresh_thermal.cells());
            assert_eq!(owner.layers(), fresh.layers());
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
