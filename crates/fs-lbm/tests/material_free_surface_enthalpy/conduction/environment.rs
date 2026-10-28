//! G1/G3/G4: environmental source data reaches moving-liquid heat and momentum.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceGasContact, FreeSurfaceGasFace,
};
use fs_material::thermochemical_gas::ThermochemicalGasModel;

// Reuse the existing authored NASA9/transport source fixture. No empirical-gas
// accuracy claim is made by these software/physics coupling checks.
#[path = "../../nsf3_thermochemical/source.rs"]
mod gas_source;

#[path = "atmosphere.rs"]
mod atmosphere;

#[path = "radiation.rs"]
mod radiation;

#[path = "surface.rs"]
mod surface;

const EXPOSED: [usize; 3] = [3, 3, 3];

fn gas_contact<'a>(
    model: &'a ThermochemicalGasModel,
    context: &'a QueryPoint,
) -> FreeSurfaceGasContact<'a> {
    FreeSurfaceGasContact {
        face: FreeSurfaceGasFace {
            liquid: EXPOSED,
            gas: [4, 3, 3],
            area_m2: 1e-6,
            distances_m: [0.0005; 2],
        },
        model,
        context,
    }
}

fn with_exposed_liquid<T>(
    card: &MaterialCard,
    moving: bool,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> T,
) -> T {
    with_exposed_caloric_liquid(card, moving, &caloric(card), 200000.0, run)
}

fn with_exposed_caloric_liquid<T>(
    card: &MaterialCard,
    moving: bool,
    chart: &EquilibriumEnthalpyPhaseCurve,
    initial_h: f64,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> T,
) -> T {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut cells = Vec::new();
    for z in 1..=6 {
        for y in 1..=6 {
            for x in 0..=3 {
                cells.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    velocity_m_s: if moving {
                        [0.0, 0.003 * (x as f64 - 1.5), 0.0]
                    } else {
                        [0.0; 3]
                    },
                });
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| initial_h,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

#[test]
fn g1_gas_film_matches_independent_layers_and_counts_bath_heat_once() {
    let liquid = conducting_card(1.0);
    let gas = gas_source::model(0.0, 250000.0, 100.0);
    let context = gas_source::point(450.0, 100000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let before = thermal.total_enthalpy_j().unwrap();
        let report = thermal
            .step_with_gas_environment(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[gas_contact(&gas, &context)],
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let environment = report.environment.unwrap();
        // Independent series resistance, source interpolation and finite liquid
        // storage. No production source or thermal helper computes this result.
        let kg = 0.026 * 100.0 * (0.5 + 99000.0 / 249000.0);
        let mut t = 400.0;
        let mut heat = 0.0;
        let mut last_g = 0.0;
        for _ in 0..4 {
            last_g = 1e-6 / (0.0005 / (t - 200.0) + 0.0005 / kg);
            let q = 0.00025 * last_g * (450.0 - t);
            t += q / 0.002;
            heat += q;
        }
        assert_eq!(report.thermal.substeps, 4);
        near(environment.last_applied[0].conductance_w_k, last_g, 2e-13);
        near(
            environment.thermal_endpoint[0].liquid_temperature_k,
            t,
            2e-13,
        );
        near(
            environment.thermal_endpoint[0].gas.state().temperature,
            0.5 * (t + 450.0),
            2e-13,
        );
        near(environment.heat_to_liquid_j, heat, 2e-12);
        near(report.flow.external_heat_j, heat, 2e-12);
        near(report.thermal.external_heat_j, heat, 2e-12);
        near(thermal.total_enthalpy_j().unwrap(), before + heat, 2e-12);
        assert!(report.flow.balance_defect_j.abs() < 1e-9);
        assert_eq!(
            environment.thermal_endpoint[0]
                .liquid
                .material()
                .card_identity(),
            liquid.content_hash()
        );
        assert_eq!(
            environment.thermal_endpoint[0].gas.source_packs(),
            gas.source_packs()
        );
    });
}

#[test]
fn g3_gas_material_data_changes_actual_liquid_viscosity_and_native_populations() {
    let run = |conductivity_scale| {
        let liquid = conducting_card(1.0);
        let gas = gas_source::model(0.0, 250000.0, conductivity_scale);
        let context = gas_source::point(450.0, 100000.0);
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let report = thermal
                .step_with_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[gas_contact(&gas, &context)],
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal
                    .material_at(EXPOSED)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                report.environment.unwrap().heat_to_liquid_j,
            )
        })
    };
    let slow = run(10.0);
    let fast = run(1000.0);
    assert!(fast.2 > slow.2);
    assert!(fast.1 < slow.1);
    assert_ne!(
        fast.0, slow.0,
        "gas data must reach native motion, not only a report"
    );
}

#[test]
fn g1_environment_wall_molecular_and_authored_heat_share_one_account() {
    let liquid = conducting_card(1.0);
    let gas = gas_source::model(0.0, 250000.0, 100.0);
    let gas_point = gas_source::point(350.0, 100000.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    let contact_source = body::interface_card(&liquid, &solid, 1e-5, 500.0, false);
    let contact_point = body::interface_point();
    let binding = body::interface(&liquid, &solid, &contact_source, &contact_point);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
        let heat = FreeSurfaceCellHeat {
            coordinates: EXPOSED,
            energy_j: 0.001,
        };
        let report = thermal
            .step_with_body_interfaces_and_viscous_heating(
                &[],
                &[heat],
                controls(),
                &mut bodies,
                &[body::contact()],
                std::slice::from_ref(&binding),
                &[gas_contact(&gas, &gas_point)],
                gate,
            )
            .unwrap()
            .conduction;
        let environment = report.environment.unwrap();
        let wall = report.bodies.unwrap();
        assert_eq!(
            wall.last_applied[0]
                .interface
                .as_ref()
                .unwrap()
                .card_identity(),
            contact_source.content_hash(),
        );
        assert!(environment.heat_to_liquid_j < 0.0);
        assert!(wall.heat_to_liquid_j > 0.0);
        assert!(report.flow.viscous_heat_j > 0.0);
        near(
            report.flow.external_heat_j,
            heat.energy_j + environment.heat_to_liquid_j,
            2e-12,
        );
        near(
            report.thermal.external_heat_j,
            report.flow.external_heat_j + report.flow.viscous_heat_j,
            2e-12,
        );
        near(
            thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j(),
            before + heat.energy_j + environment.heat_to_liquid_j + report.flow.viscous_heat_j,
            2e-12,
        );
        assert!(wall.combined_balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g1_combined_conduction_and_environment_limit_the_same_thermal_step() {
    let liquid = conducting_card(1.0);
    let gas = gas_source::model(0.0, 250000.0, 100.0);
    let context = gas_source::point(450.0, 100000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut contact = gas_contact(&gas, &context);
        contact.face.distances_m = [1e-6; 2];
        let face = FreeSurfaceConductionFace {
            cells: [[2, 3, 3], EXPOSED],
            area_m2: 1e-6,
            distances_m: [1e-4; 2],
        };
        let report = thermal
            .step_with_gas_environment(
                &[face],
                &[],
                FreeSurfaceConductionControls {
                    maximum_dt_s: 0.001,
                    ..controls()
                },
                &mut [],
                &[],
                &[contact],
                false,
                gate,
            )
            .unwrap()
            .conduction;
        assert!(
            report.thermal.substeps > 1,
            "gas-film stiffness must subcycle the owner"
        );
        let t = report.environment.unwrap().thermal_endpoint[0].liquid_temperature_k;
        assert!(t > 400.0 && t < 450.0);
        assert!(report.flow.balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g4_stale_faces_and_late_gas_source_refusal_preserve_liquid_and_wall() {
    let liquid = conducting_card(1.0);
    let gas = gas_source::model_with_temperature_limit(0.0, 250000.0, 100.0, 450.0);
    let gas_point = gas_source::point(450.0, 100000.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let bits = thermal.flow().lattice().grid().state_bits();
        let before = thermal.total_enthalpy_j().unwrap();
        let wall_before = bodies[0].specific_enthalpy_j_kg();
        let source_before = bodies[0].material().clone();
        let contact = gas_contact(&gas, &gas_point);
        let mut stale = contact;
        stale.face.gas = [2, 3, 3];
        for contacts in [vec![stale], vec![contact, contact]] {
            assert!(
                thermal
                    .step_with_gas_environment(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[body::contact()],
                        &contacts,
                        false,
                        gate,
                    )
                    .is_err()
            );
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
        }
        // First two substeps are supported. Applied heat then moves the gas
        // film above its 450 K source ceiling while the liquid stays in-domain.
        let heat = FreeSurfaceCellHeat {
            coordinates: EXPOSED,
            energy_j: 0.18,
        };
        let error = thermal
            .step_with_gas_environment(
                &[],
                &[heat],
                controls(),
                &mut bodies,
                &[body::contact()],
                &[contact],
                false,
                gate,
            )
            .unwrap_err();
        assert!(matches!(
            error,
            FreeSurfaceEnthalpyError::Flow(fs_lbm::material::MaterialLbmError::Transport(_))
        ));
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
        assert_eq!(bodies[0].material(), &source_before);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            thermal
                .step_with_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[body::contact()],
                    &[contact],
                    false,
                    &cancelled,
                )
                .is_err()
        );
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
        // A valid retry starts from exactly the previously accepted histories.
        thermal
            .step_with_gas_environment(
                &[],
                &[],
                controls(),
                &mut bodies,
                &[body::contact()],
                &[contact],
                false,
                gate,
            )
            .unwrap();
    });
}
