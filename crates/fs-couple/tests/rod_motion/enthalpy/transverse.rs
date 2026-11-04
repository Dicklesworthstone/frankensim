//! Transverse material loss heats the original finite cells exactly once.
//! Manufactured viscosities and heat capacities, not measured calibrations.
use super::*;
use fs_couple::rod_motion::enthalpy::{
    RodAxialViscosity, RodBendingViscosity, RodMaterialViscosity, RodTorsionalViscosity,
    RodTransverseViscosity,
};

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

fn selection(selections: &[MaterialPropertySelection]) -> RodTransverseViscosity<'_> {
    RodTransverseViscosity {
        property: "shear_viscosity",
        selections,
    }
}

fn shear_motion(
    source: &RodMaterialSource<'_>,
    temperatures: &[f64],
    q: [f64; 2],
    v: [f64; 2],
) -> RodMotion {
    let gate = CancelGate::new();
    let rod = source
        .resolve_at(temperatures, &gate)
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(LENGTH), None, None],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
    let model = Arc::new(
        rod.into_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.positions_m[1][1..].copy_from_slice(&q);
    RodMotion::new(
        model,
        pose,
        &[[0.0; 3], [0.0, v[0], v[1]]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap()
}

fn tip_mass() -> f64 {
    0.5 * RHO * LENGTH * core::f64::consts::PI * RADIUS.powi(2)
}

// Two independent physical shear oscillators: exact viscous half steps around
// the elastic Verlet map. No native section, loss or thermal receipt supplies
// this reference; kappa enters both shear stiffness and viscosity exactly once.
fn scalar_step(q: [f64; 2], v: [f64; 2], eta: f64, dt: f64) -> ([f64; 2], [f64; 2], f64) {
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let mass = tip_mass();
    let stiffness = 0.8 * SHEAR * area / LENGTH;
    let resistance = 0.8 * eta * area / LENGTH;
    let decay = (-0.5 * dt * resistance / mass).exp();
    let mut next_q = [0.0; 2];
    let mut next_v = [0.0; 2];
    let mut heat = 0.0;
    for axis in 0..2 {
        let middle = decay * v[axis] - 0.5 * dt * stiffness * q[axis] / mass;
        next_q[axis] = q[axis] + dt * middle;
        let before_loss = middle - 0.5 * dt * stiffness * next_q[axis] / mass;
        next_v[axis] = decay * before_loss;
        heat +=
            0.5 * mass * (1.0 - decay * decay) * (v[axis] * v[axis] + before_loss * before_loss);
    }
    (next_q, next_v, heat)
}

#[test]
fn g1_equal_shear_weights_heat_unequal_layer_capacities_and_close_energy_entropy() {
    let specs = [
        specimen(60.0),
        Spec {
            cp: 2e-4,
            ..specimen(20.0)
        },
    ];
    fixture(
        &specs,
        1,
        &[300.0; 2],
        &[0.0; 2],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let (q0, v0, dt) = ([0.004, -0.003], [0.08, -0.04], 0.003);
            let (q, v, heat) = scalar_step(q0, v0, 30.0, dt);
            let motion = shear_motion(source, &[300.0; 2], q0, v0);
            let initial_caloric = total_heat(&thermal);
            let routes = destinations(2, 1);
            let selections = vec![MaterialPropertySelection::SingleClaimOnly; 2];
            let mut owner = RodThermoelasticSystem::try_new(
                motion, source, &thermal, network, &routes, 1e-11, &gate,
            )
            .unwrap();
            let initial_mechanical = owner.mechanical_internal_energy_j().unwrap();
            let report = owner
                .step_material_viscous_with_transverse(
                    &mut thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    RodMaterialViscosity::default(),
                    selection(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0; 2],
                    &gate,
                )
                .unwrap();
            let heating = report.transverse.as_ref().unwrap();
            assert!(heat > 0.0 && owner.mechanical().accepted_steps() == 1);
            assert!(
                report.axial.is_none() && report.bending.is_none() && report.torsional.is_none()
            );
            near(heating.mechanical_heat_j, heat, 3e-14);
            near(heating.deposited_heat_j, heat, 3e-14);
            near(heating.partition_residual_j, 0.0, 3e-20);
            near(
                report
                    .mechanical
                    .transverse_heat
                    .as_ref()
                    .unwrap()
                    .total_heat_j,
                heat,
                3e-14,
            );
            for axis in 0..2 {
                near(
                    owner.mechanical().pose().positions_m[1][axis + 1],
                    q[axis],
                    3e-13,
                );
                near(
                    owner.mechanical().velocities_m_s().unwrap()[1][axis + 1],
                    v[axis],
                    3e-12,
                );
            }
            let mut generation = 0.0;
            let mut entropy = 0.0;
            let mut rises = [0.0; 2];
            for layer in 0..2 {
                let capacity = RHO * LENGTH * geometry(layer, 2).0 * specs[layer].cp;
                let expected_t = 300.0 + 0.5 * heat / capacity;
                let actual_t = thermal.temperature_k(routes[layer].cell).unwrap();
                let receipt = &heating.layers[layer];
                assert_eq!(receipt.destination, routes[layer]);
                near(
                    receipt.viscosity_pa_s,
                    specs[layer].viscosity.unwrap().value_pa_s,
                    0.0,
                );
                near(receipt.area_m2, geometry(layer, 2).0, 2e-18);
                near(receipt.shear_correction, 0.8, 0.0);
                near(
                    heating.source_sections[0].layer_viscous_rigidity_n_s()[layer],
                    0.8 * receipt.viscosity_pa_s * receipt.area_m2,
                    2e-16,
                );
                near(receipt.heat_j, 0.5 * heat, 3e-14);
                near(actual_t, expected_t, 3e-8);
                near(capacity * (actual_t - 300.0), receipt.heat_j, 3e-14);
                rises[layer] = actual_t - 300.0;
                generation += 0.5 * heat / expected_t;
                entropy += capacity * (expected_t / 300.0).ln();
            }
            near(rises[0], 6.0 * rises[1], 3e-7);
            near(total_heat(&thermal) - initial_caloric, heat, 3e-14);
            near(report.thermal.external_heat_j, heat, 3e-14);
            near(report.external_thermal_heat_j, 0.0, 3e-14);
            near(report.reversible_heat_j, 0.0, 3e-17);
            near(
                owner.mechanical_internal_energy_j().unwrap() - initial_mechanical
                    + total_heat(&thermal)
                    - initial_caloric,
                report.mechanical.integration_defect_j,
                4e-14,
            );
            near(
                report.total_energy_residual_j,
                report.mechanical.integration_defect_j,
                4e-14,
            );
            near(report.coupling_energy_residual_j, 0.0, 4e-14);
            near(report.entropy.transverse_generation_j_k, generation, 3e-15);
            near(report.entropy.caloric_storage_change_j_k, entropy, 3e-15);
            near(report.entropy.axial_generation_j_k, 0.0, 0.0);
            near(report.entropy.bending_generation_j_k, 0.0, 0.0);
            near(report.entropy.torsional_generation_j_k, 0.0, 0.0);
            near(report.entropy.accounting_residual_j_k, 0.0, 3e-15);
        },
    );
}

#[test]
fn g1_transverse_viscosity_refreshes_from_actual_layer_temperature() {
    let mut spec = specimen(40.0);
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
            let (mut q, mut v, mut t) = ([0.004, -0.003], [0.08, -0.04], 300.0);
            let motion = shear_motion(source, &[t], q, v);
            let capacity = RHO * LENGTH * geometry(0, 1).0 * spec.cp;
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
            let selections = [MaterialPropertySelection::SingleClaimOnly];
            let mut total = 0.0;
            let initial_caloric = total_heat(&thermal);
            for _ in 0..5 {
                let eta = 40.0 + 20.0 * (t - 300.0);
                let dt = 0.001;
                let (next_q, next_v, heat) = scalar_step(q, v, eta, dt);
                let next_t = t + heat / capacity;
                let target = owner.mechanical().time_s() + dt;
                let report = owner
                    .step_material_viscous_with_transverse(
                        &mut thermal,
                        target,
                        &[RodNodalWrench::default(); 2],
                        RodMaterialViscosity::default(),
                        selection(&selections),
                        controls(dt * 1.001),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let heating = report.transverse.as_ref().unwrap();
                let receipt = &heating.source_sections[0].shear().elastic().layers()[0].material;
                near(heating.layers[0].viscosity_pa_s, eta, 3e-6);
                near(
                    receipt.property("shear_viscosity").unwrap().value_si(),
                    eta,
                    3e-6,
                );
                near(
                    receipt
                        .query_point()
                        .iter()
                        .find(|(axis, _)| axis == "T")
                        .unwrap()
                        .1,
                    t,
                    2e-7,
                );
                near(thermal.temperature_k(CellId(0)).unwrap(), next_t, 2e-7);
                for axis in 0..2 {
                    near(
                        owner.mechanical().pose().positions_m[1][axis + 1],
                        next_q[axis],
                        3e-11,
                    );
                    near(
                        owner.mechanical().velocities_m_s().unwrap()[1][axis + 1],
                        next_v[axis],
                        3e-9,
                    );
                }
                near(heating.deposited_heat_j, heat, 3e-14);
                near(report.external_thermal_heat_j, 0.0, 3e-14);
                near(
                    report.entropy.transverse_generation_j_k,
                    heat / next_t,
                    3e-15,
                );
                near(report.entropy.accounting_residual_j_k, 0.0, 3e-15);
                near(report.coupling_energy_residual_j, 0.0, 4e-14);
                total += heat;
                (q, v, t) = (next_q, next_v, next_t);
            }
            assert!(t > 300.5);
            near(total_heat(&thermal) - initial_caloric, total, 3e-13);
        },
    );
}

#[test]
fn g3_all_material_heat_shares_cells_and_zero_transverse_preserves_legacy_paths() {
    for shear_eta in [0.0, 12.0] {
        let template = specimen(shear_eta);
        let spec = Spec {
            axial_viscosity: Some(ViscositySpec {
                value_pa_s: 20.0,
                ..template.viscosity.unwrap()
            }),
            bending_viscosity: Some(ViscositySpec {
                value_pa_s: 80.0,
                ..template.viscosity.unwrap()
            }),
            ..template
        };
        fixture(
            &[spec],
            1,
            &[300.0],
            &[0.0],
            false,
            |source, mut thermal, network| {
                let gate = CancelGate::new();
                let model = Arc::new(
                    source
                        .resolve_at(&[300.0], &gate)
                        .unwrap()
                        .into_motion_model(RodRootBoundary::Clamped, &gate)
                        .unwrap(),
                );
                let mut pose = model.initial_pose();
                pose.positions_m[1][0] += 0.0008;
                pose.frames[1] =
                    So3::exp(So3Tangent::new(Vec3::new(0.003, 0.005, -0.002))).unwrap();
                let motion = RodMotion::new(
                    model,
                    pose,
                    &[[0.0; 3], [0.08, 0.02, -0.01]],
                    &[[0.0; 3], [0.6, 0.5, -0.3]],
                    0.0,
                    &gate,
                )
                .unwrap();
                let mut spatial_thermal = thermal.clone();
                let mut legacy_thermal = thermal.clone();
                let routes = destinations(1, 1);
                let mut owner = RodThermoelasticSystem::try_new(
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
                    motion.clone(),
                    source,
                    &spatial_thermal,
                    network,
                    &routes,
                    1e-11,
                    &gate,
                )
                .unwrap();
                let mut legacy = RodThermoelasticSystem::try_new(
                    motion,
                    source,
                    &legacy_thermal,
                    network,
                    &routes,
                    1e-11,
                    &gate,
                )
                .unwrap();
                let selections = [MaterialPropertySelection::SingleClaimOnly];
                let material = RodMaterialViscosity {
                    axial: Some(RodAxialViscosity {
                        property: "axial_viscosity",
                        selections: &selections,
                    }),
                    bending: Some(RodBendingViscosity {
                        property: "bending_viscosity",
                        selections: &selections,
                    }),
                    torsional: Some(RodTorsionalViscosity {
                        property: "shear_viscosity",
                        selections: &selections,
                    }),
                };
                let loads = [RodNodalWrench::default(); 2];
                let prepared = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
                let initial_caloric = total_heat(&thermal);
                let dt = 0.0005;
                let report = owner
                    .step_material_viscous_with_transverse(
                        &mut thermal,
                        dt,
                        &loads,
                        material,
                        selection(&selections),
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let spatial_report = spatial
                    .step_spatial_material_viscous_with_transverse(
                        &mut spatial_thermal,
                        dt,
                        &prepared,
                        material,
                        selection(&selections),
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let transverse = report.transverse.as_ref().unwrap().deposited_heat_j;
                let axial = report.axial.as_ref().unwrap().deposited_heat_j;
                let bending = report.bending.as_ref().unwrap().deposited_heat_j;
                let torsional = report.torsional.as_ref().unwrap().deposited_heat_j;
                assert!(axial > 0.0 && bending > 0.0);
                assert_eq!(transverse > 0.0, shear_eta > 0.0);
                assert_eq!(torsional > 0.0, shear_eta > 0.0);
                let total = axial + bending + torsional + transverse;
                let expected_t = 300.0 + total / (RHO * LENGTH * geometry(0, 1).0 * spec.cp);
                near(thermal.temperature_k(CellId(0)).unwrap(), expected_t, 3e-8);
                near(total_heat(&thermal) - initial_caloric, total, 3e-14);
                near(report.external_thermal_heat_j, 0.0, 3e-14);
                near(
                    report.entropy.transverse_generation_j_k,
                    transverse / expected_t,
                    3e-15,
                );
                near(report.entropy.accounting_residual_j_k, 0.0, 3e-15);
                near(report.coupling_energy_residual_j, 0.0, 4e-14);
                assert_eq!(owner.mechanical().pose(), spatial.mechanical().pose());
                assert_eq!(
                    owner.mechanical().linear_momenta_kg_m_s(),
                    spatial.mechanical().linear_momenta_kg_m_s()
                );
                assert_eq!(
                    owner.mechanical().body_spins_kg_m2_s(),
                    spatial.mechanical().body_spins_kg_m2_s()
                );
                assert_eq!(thermal.cells(), spatial_thermal.cells());
                assert_eq!(report.transverse, spatial_report.transverse);
                assert_eq!(report.axial, spatial_report.axial);
                assert_eq!(report.bending, spatial_report.bending);
                assert_eq!(report.torsional, spatial_report.torsional);
                if shear_eta == 0.0 {
                    let old = legacy
                        .step_material_viscous(
                            &mut legacy_thermal,
                            dt,
                            &loads,
                            material,
                            controls(dt),
                            heat_solve(),
                            &[0.0],
                            &gate,
                        )
                        .unwrap();
                    assert_eq!(owner.mechanical().pose(), legacy.mechanical().pose());
                    assert_eq!(
                        owner.mechanical().linear_momenta_kg_m_s(),
                        legacy.mechanical().linear_momenta_kg_m_s()
                    );
                    assert_eq!(
                        owner.mechanical().body_spins_kg_m2_s(),
                        legacy.mechanical().body_spins_kg_m2_s()
                    );
                    assert_eq!(thermal.cells(), legacy_thermal.cells());
                    assert!(old.transverse.is_none() && old.mechanical.transverse_heat.is_none());
                    assert_eq!(report.total_energy_residual_j, old.total_energy_residual_j);
                }
            },
        );
    }
}

#[test]
fn g4_transverse_endpoint_refusal_rolls_back_both_histories_and_allows_retry() {
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
            let original = shear_motion(source, &[300.0], [0.0; 2], [0.08, -0.04]);
            let original_thermal = thermal.clone();
            let routes = destinations(1, 1);
            let selections = [MaterialPropertySelection::SingleClaimOnly];
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
            let attempted = 0.01;
            let (_, _, heat) = scalar_step([0.0; 2], [0.08, -0.04], 100.0, attempted);
            let forbidden_t = 300.0 + heat / (RHO * LENGTH * geometry(0, 1).0 * spec.cp);
            assert!(forbidden_t > 300.01 && forbidden_t < 500.0);
            assert!(source.resolve_at(&[forbidden_t], &gate).is_ok());
            assert!(
                owner
                    .step_material_viscous_with_transverse(
                        &mut thermal,
                        attempted,
                        &[RodNodalWrench::default(); 2],
                        RodMaterialViscosity::default(),
                        selection(&selections),
                        controls(attempted),
                        heat_solve(),
                        &[0.0],
                        &gate
                    )
                    .is_err()
            );
            assert_eq!(thermal.cells(), original_thermal.cells());
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(owner.mechanical().pose(), original.pose());
            assert_eq!(
                owner.mechanical().linear_momenta_kg_m_s(),
                original.linear_momenta_kg_m_s()
            );
            assert_eq!(
                owner.mechanical().body_spins_kg_m2_s(),
                original.body_spins_kg_m2_s()
            );
            assert_eq!(owner.mechanical().time_s(), 0.0);
            assert_eq!(owner.mechanical().accepted_steps(), 0);
            assert_eq!(owner.layers(), layers);
            let dt = 1e-6;
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
                .step_material_viscous_with_transverse(
                    &mut thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    RodMaterialViscosity::default(),
                    selection(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            let direct = fresh
                .step_material_viscous_with_transverse(
                    &mut fresh_thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    RodMaterialViscosity::default(),
                    selection(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            assert!(thermal.temperature_k(CellId(0)).unwrap() < 300.01);
            assert!(retry.transverse.as_ref().unwrap().deposited_heat_j > 0.0);
            assert_eq!(thermal.cells(), fresh_thermal.cells());
            assert_eq!(owner.layers(), fresh.layers());
            assert_eq!(owner.mechanical().pose(), fresh.mechanical().pose());
            assert_eq!(
                owner.mechanical().linear_momenta_kg_m_s(),
                fresh.mechanical().linear_momenta_kg_m_s()
            );
            assert_eq!(retry.transverse, direct.transverse);
            assert_eq!(
                retry.total_energy_residual_j,
                direct.total_energy_residual_j
            );
        },
    );
}
