//! Axial viscous work heats the original finite material layers exactly once.
//! Manufactured section and caloric data, not experimental material calibration.
use super::enthalpy::{
    LENGTH, RADIUS, RHO, Spec, ViscositySpec, YOUNG, destinations, fixture, geometry, solve,
    total_heat,
};
use super::*;
use fs_couple::rod_motion::enthalpy::{
    RodAngularViscosity, RodAxialViscosity, RodBendingViscosity, RodMaterialViscosity,
    RodThermoelasticSolve, RodThermoelasticSystem, RodTorsionalViscosity,
};
use fs_couple::rod_motion::{PreparedRodLoads, RodMotion};
use fs_solid::rod::material_field::thermal::source::RodMaterialSource;
use fs_thermal::enthalpy::CellId;

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
        axial_viscosity: Some(viscosity(eta)),
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
fn selection(selections: &[MaterialPropertySelection]) -> RodMaterialViscosity<'_> {
    RodMaterialViscosity {
        axial: Some(RodAxialViscosity {
            property: "axial_viscosity",
            selections,
        }),
        bending: None,
        torsional: None,
    }
}
fn axial_motion(source: &RodMaterialSource<'_>, temperatures: &[f64], q: f64, v: f64) -> RodMotion {
    let gate = CancelGate::new();
    let rod = source
        .resolve_at(temperatures, &gate)
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [None, Some(0.0), Some(0.0)],
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
    pose.positions_m[1][0] += q;
    RodMotion::new(
        model,
        pose,
        &[[0.0; 3], [v, 0.0, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap()
}
fn tip_mass() -> f64 {
    0.5 * RHO * LENGTH * core::f64::consts::PI * RADIUS.powi(2)
}

// Independent one-coordinate split: two exact dashpot exponentials around
// elastic velocity Verlet. No native loss receipt, section integral, source
// resolver or thermal solver supplies this motion/physical-heat reference.
fn scalar_step(q: f64, v: f64, young: f64, eta: f64, dt: f64) -> (f64, f64, f64) {
    let area = core::f64::consts::PI * RADIUS.powi(2);
    let mass = tip_mass();
    let stiffness = young * area / LENGTH;
    let resistance = eta * area / LENGTH;
    let decay = (-0.5 * dt * resistance / mass).exp();
    let middle = decay * v - 0.5 * dt * stiffness * q / mass;
    let next_q = q + dt * middle;
    let before_loss = middle - 0.5 * dt * stiffness * next_q / mass;
    let heat = 0.5 * mass * (1.0 - decay * decay) * (v * v + before_loss * before_loss);
    (next_q, decay * before_loss, heat)
}

#[test]
fn g1_eta_area_partition_heats_actual_layer_capacities_and_closes_energy_entropy() {
    let specs = [
        Spec {
            cp: 2e-4,
            ..specimen(30.0)
        },
        specimen(10.0),
    ];
    fixture(
        &specs,
        1,
        &[300.0; 2],
        &[0.0; 2],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let (q0, v0, dt) = (0.0008, 0.08, 0.004);
            let original_caloric = total_heat(&thermal);
            let routes = destinations(2, 1);
            let selections = vec![MaterialPropertySelection::SingleClaimOnly; 2];
            let motion = axial_motion(source, &[300.0; 2], q0, v0);
            let initial_energy = 0.5 * tip_mass() * v0 * v0
                + 0.5 * YOUNG * core::f64::consts::PI * RADIUS.powi(2) / LENGTH * q0 * q0;
            let mut owner = RodThermoelasticSystem::try_new(
                motion, source, &thermal, network, &routes, 1e-11, &gate,
            )
            .unwrap();
            let report = owner
                .step_material_viscous(
                    &mut thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    selection(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0; 2],
                    &gate,
                )
                .unwrap();
            // Core/annulus areas are 1:3. Eta=30:10 gives equal heat, whereas
            // mass weighting gives 1:3 and bending-moment weighting gives 1:5.
            let (q, v, heat) = scalar_step(q0, v0, YOUNG, 15.0, dt);
            near(
                owner.mechanical().pose().positions_m[1][0],
                LENGTH + q,
                3e-13,
            );
            near(owner.mechanical().velocities_m_s().unwrap()[1][0], v, 3e-12);
            let heating = report.axial.as_ref().unwrap();
            assert_eq!(heating.layers.len(), 2);
            assert_eq!(heating.source_sections.len(), 1);
            assert!(report.bending.is_none() && report.torsional.is_none());
            near(heating.mechanical_heat_j, heat, 3e-15);
            near(heating.deposited_heat_j, heat, 3e-15);
            near(heating.partition_residual_j, 0.0, 3e-20);
            near(
                report.mechanical.axial_heat.as_ref().unwrap().total_heat_j,
                heat,
                3e-15,
            );
            let mut entropy = 0.0;
            let mut generation = 0.0;
            let mut rises = [0.0; 2];
            for layer in 0..2 {
                let area = geometry(layer, 2).0;
                let capacity = RHO * LENGTH * area * specs[layer].cp;
                let expected_t = 300.0 + 0.5 * heat / capacity;
                let actual_t = thermal.temperature_k(routes[layer].cell).unwrap();
                near(actual_t, expected_t, 3e-8);
                rises[layer] = actual_t - 300.0;
                let receipt = &heating.layers[layer];
                assert_eq!(receipt.destination, routes[layer]);
                near(receipt.area_m2, area, 2e-18);
                near(
                    receipt.viscosity_pa_s,
                    specs[layer].axial_viscosity.unwrap().value_pa_s,
                    0.0,
                );
                near(receipt.heat_j, 0.5 * heat, 3e-15);
                near(capacity * (actual_t - 300.0), receipt.heat_j, 3e-14);
                near(
                    heating.source_sections[0].layer_viscous_rigidity_n_s()[layer],
                    receipt.viscosity_pa_s * area,
                    2e-16,
                );
                entropy += capacity * (expected_t / 300.0).ln();
                generation += 0.5 * heat / expected_t;
            }
            near(rises[0], 1.5 * rises[1], 6e-8);
            near(report.reversible_heat_j, 0.0, 3e-17);
            near(report.external_thermal_heat_j, 0.0, 3e-14);
            near(report.thermal.external_heat_j, heat, 3e-14);
            near(total_heat(&thermal) - original_caloric, heat, 3e-14);
            let energy_change = 0.5 * tip_mass() * v * v
                + 0.5 * YOUNG * core::f64::consts::PI * RADIUS.powi(2) / LENGTH * q * q
                - initial_energy
                + total_heat(&thermal)
                - original_caloric;
            near(energy_change, report.mechanical.integration_defect_j, 4e-14);
            assert!(report.mechanical.integration_defect_j.abs() > 1e-12);
            near(
                report.total_energy_residual_j,
                report.mechanical.integration_defect_j,
                4e-14,
            );
            near(report.coupling_energy_residual_j, 0.0, 4e-14);
            near(report.entropy.caloric_storage_change_j_k, entropy, 3e-15);
            near(report.entropy.axial_generation_j_k, generation, 3e-15);
            near(report.entropy.bending_generation_j_k, 0.0, 0.0);
            near(report.entropy.torsional_generation_j_k, 0.0, 0.0);
            near(report.entropy.accounting_residual_j_k, 0.0, 3e-15);
        },
    );
}

#[test]
fn g1_local_temperature_refreshes_axial_claim_without_relabeling_reversible_exchange() {
    let mut spec = specimen(40.0);
    spec.young_slope = -1.5;
    spec.axial_viscosity.as_mut().unwrap().slope_pa_s_k = 20.0;
    spec.axial_viscosity.as_mut().unwrap().lower_k = 299.0;
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let (mut q, mut v, mut t) = (0.0008, 0.08, 300.0);
            let dt = 0.003;
            let area = geometry(0, 1).0;
            let capacity = RHO * LENGTH * area * spec.cp;
            let internal_stiffness = (YOUNG - 300.0 * spec.young_slope) * area / LENGTH;
            let initial_energy = 0.5 * tip_mass() * v * v + 0.5 * internal_stiffness * q * q;
            let initial_caloric = total_heat(&thermal);
            let motion = axial_motion(source, &[t], q, v);
            let selections = [MaterialPropertySelection::SingleClaimOnly];
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
            let (mut viscous_total, mut reversible_total, mut defects) = (0.0, 0.0, 0.0);
            for _ in 0..8 {
                let eta = 40.0 + 20.0 * (t - 300.0);
                let young = YOUNG + spec.young_slope * (t - 300.0);
                let (next_q, next_v, heat) = scalar_step(q, v, young, eta, dt);
                let entropy_change =
                    -0.5 * spec.young_slope * area / LENGTH * (next_q * next_q - q * q);
                let reversible = -t * entropy_change;
                let next_t = t + (heat + reversible) / capacity;
                let target = owner.mechanical().time_s() + dt;
                let report = owner
                    .step_material_viscous(
                        &mut thermal,
                        target,
                        &[RodNodalWrench::default(); 2],
                        selection(&selections),
                        controls(dt * 1.001),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let heating = report.axial.as_ref().unwrap();
                near(heating.layers[0].viscosity_pa_s, eta, 2e-6);
                let receipt = &heating.source_sections[0].elastic().layers()[0].material;
                near(
                    receipt.property("axial_viscosity").unwrap().value_si(),
                    eta,
                    2e-6,
                );
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
                    owner.mechanical().pose().positions_m[1][0] - LENGTH,
                    next_q,
                    3e-11,
                );
                near(
                    owner.mechanical().velocities_m_s().unwrap()[1][0],
                    next_v,
                    3e-9,
                );
                near(thermal.temperature_k(CellId(0)).unwrap(), next_t, 2e-7);
                near(heating.deposited_heat_j, heat, 3e-14);
                near(report.reversible_heat_j, reversible, 3e-14);
                near(report.thermal.external_heat_j, heat + reversible, 3e-14);
                near(report.external_thermal_heat_j, 0.0, 3e-14);
                near(report.entropy.axial_generation_j_k, heat / next_t, 3e-15);
                near(
                    report.entropy.mechanical_storage_change_j_k,
                    entropy_change,
                    3e-15,
                );
                near(
                    report.entropy.thermoelastic_temporal_defect_j_k,
                    entropy_change + reversible / next_t,
                    3e-15,
                );
                near(report.entropy.accounting_residual_j_k, 0.0, 3e-15);
                near(report.coupling_energy_residual_j, 0.0, 4e-14);
                viscous_total += heat;
                reversible_total += reversible;
                defects += report.mechanical.integration_defect_j;
                (q, v, t) = (next_q, next_v, next_t);
            }
            assert!(
                t > 300.5,
                "fixture must materially change subsequent viscosity"
            );
            assert!(reversible_total.abs() > 1e-9);
            near(
                total_heat(&thermal) - initial_caloric,
                viscous_total + reversible_total,
                3e-13,
            );
            near(
                0.5 * tip_mass() * v * v + 0.5 * internal_stiffness * q * q - initial_energy
                    + total_heat(&thermal)
                    - initial_caloric,
                defects,
                3e-13,
            );
        },
    );
}

#[test]
fn g3_axial_and_angular_heat_share_cells_and_world_spatial_zero_loss_paths_agree() {
    for axial_eta in [0.0, 20.0] {
        let spec = Spec {
            bending_viscosity: Some(viscosity(80.0)),
            viscosity: Some(viscosity(12.0)),
            ..specimen(axial_eta)
        };
        fixture(
            &[spec],
            1,
            &[300.0],
            &[0.0],
            false,
            |source, mut thermal, network| {
                let gate = CancelGate::new();
                let rod = source.resolve_at(&[300.0], &gate).unwrap();
                let model = Arc::new(
                    rod.into_motion_model(RodRootBoundary::Clamped, &gate)
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
                let angular = RodAngularViscosity {
                    bending: Some(RodBendingViscosity {
                        property: "bending_viscosity",
                        selections: &selections,
                    }),
                    torsional: Some(RodTorsionalViscosity {
                        property: "shear_viscosity",
                        selections: &selections,
                    }),
                };
                let losses = RodMaterialViscosity {
                    bending: angular.bending,
                    torsional: angular.torsional,
                    ..selection(&selections)
                };
                let loads = [RodNodalWrench::default(); 2];
                let prepared = PreparedRodLoads::new(&loads, &[[0.0; 3]; 2], &[], &gate).unwrap();
                let initial_caloric = total_heat(&thermal);
                let dt = 0.0005;
                let report = owner
                    .step_material_viscous(
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
                let spatial_report = spatial
                    .step_spatial_material_viscous(
                        &mut spatial_thermal,
                        dt,
                        &prepared,
                        losses,
                        controls(dt),
                        heat_solve(),
                        &[0.0],
                        &gate,
                    )
                    .unwrap();
                let axial = report.axial.as_ref().unwrap().deposited_heat_j;
                let bending = report.bending.as_ref().unwrap().deposited_heat_j;
                let torsional = report.torsional.as_ref().unwrap().deposited_heat_j;
                assert!(bending > 0.0 && torsional > 0.0);
                assert_eq!(axial > 0.0, axial_eta > 0.0);
                let total = axial + bending + torsional;
                let expected_t = 300.0 + total / (RHO * LENGTH * geometry(0, 1).0 * spec.cp);
                near(thermal.temperature_k(CellId(0)).unwrap(), expected_t, 3e-8);
                near(total_heat(&thermal) - initial_caloric, total, 3e-14);
                near(
                    report.entropy.axial_generation_j_k,
                    axial / expected_t,
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
                assert_eq!(report.axial, spatial_report.axial);
                assert_eq!(report.bending, spatial_report.bending);
                assert_eq!(report.torsional, spatial_report.torsional);
                if axial_eta == 0.0 {
                    let old = legacy
                        .step_angular_viscous(
                            &mut legacy_thermal,
                            dt,
                            &loads,
                            angular,
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
                    assert!(old.axial.is_none() && old.mechanical.axial_heat.is_none());
                    assert_eq!(report.total_energy_residual_j, old.total_energy_residual_j);
                }
            },
        );
    }
}

#[test]
fn g4_axial_endpoint_domain_refusal_keeps_both_histories_and_permits_valid_retry() {
    let mut spec = specimen(100.0);
    spec.axial_viscosity.as_mut().unwrap().upper_k = 300.01;
    fixture(
        &[spec],
        1,
        &[300.0],
        &[0.0],
        false,
        |source, mut thermal, network| {
            let gate = CancelGate::new();
            let original = axial_motion(source, &[300.0], 0.0, 0.08);
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
            let (_, _, heat) = scalar_step(0.0, 0.08, YOUNG, 100.0, attempted);
            let forbidden_t = 300.0 + heat / (RHO * LENGTH * geometry(0, 1).0 * spec.cp);
            assert!(forbidden_t > 300.01 && forbidden_t < 500.0);
            assert!(source.resolve_at(&[forbidden_t], &gate).is_ok());
            assert!(
                owner
                    .step_material_viscous(
                        &mut thermal,
                        attempted,
                        &[RodNodalWrench::default(); 2],
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
                .step_material_viscous(
                    &mut thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    selection(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            let direct = fresh
                .step_material_viscous(
                    &mut fresh_thermal,
                    dt,
                    &[RodNodalWrench::default(); 2],
                    selection(&selections),
                    controls(dt),
                    heat_solve(),
                    &[0.0],
                    &gate,
                )
                .unwrap();
            assert!(thermal.temperature_k(CellId(0)).unwrap() < 300.01);
            assert!(retry.axial.as_ref().unwrap().deposited_heat_j > 0.0);
            assert_eq!(thermal.cells(), fresh_thermal.cells());
            assert_eq!(owner.layers(), fresh.layers());
            assert_eq!(owner.mechanical().pose(), fresh.mechanical().pose());
            assert_eq!(
                owner.mechanical().linear_momenta_kg_m_s(),
                fresh.mechanical().linear_momenta_kg_m_s()
            );
            assert_eq!(retry.axial, direct.axial);
            assert_eq!(
                retry.total_energy_residual_j,
                direct.total_energy_residual_j
            );
        },
    );
}
