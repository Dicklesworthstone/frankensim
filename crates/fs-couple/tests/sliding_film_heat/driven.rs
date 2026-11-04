//! Pressure-fed wall friction heats the original stores with signed source work.
//! The equal-viscosity half-layers admit independent Stokes stress integrals.

use super::*;
use fs_tribo::resistive_film::shear::PressureDrivenShearedFilm;

const MU: f64 = 0.02;
const PRESSURE: f64 = 100.0;

fn sliding() -> ShearedFilm {
    let profile = GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: 0.5,
            dynamic_viscosity_pa_s: MU,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: MU,
        },
    ])
    .unwrap();
    let rows = [-1.0, 1.0]
        .into_iter()
        .map(|sign| WallVelocityPorts {
            lower: vec![sign / M1.sqrt(), 0.0],
            upper: vec![0.0, sign / M2.sqrt()],
        })
        .collect();
    ShearedFilm::from_stratified(
        StratifiedFilm::new(&graph(), vec![profile; 2]).unwrap(),
        rows,
    )
    .unwrap()
}

fn driven_mechanics(include_passive: bool) -> ImpactSystem {
    let driven = PressureDrivenShearedFilm::new(sliding(), &[Some(PRESSURE), Some(0.0)]).unwrap();
    let system = bare().with_driven_sliding_film(driven).unwrap();
    if include_passive {
        system.with_sliding_film(sliding()).unwrap()
    } else {
        system
    }
}

fn with_routes(
    c: &EquilibriumEnthalpyPhaseCurve,
    include_passive: bool,
    passive: Vec<SlidingFilmLayerHeatRoute>,
    driven: Vec<SlidingFilmLayerHeatRoute>,
) -> Result<HeatedImpactSystem<'_>, fs_couple::render::plate::impact::ImpactError> {
    HeatedImpactSystem::new_with_driven_sliding_losses(
        driven_mechanics(include_passive)
            .prepare_analytic()
            .unwrap(),
        stores(c),
        vec![],
        vec![],
        vec![],
        vec![],
        passive,
        driven,
        vec![],
        vec![],
        1e-12,
    )
}

fn heated_driven(c: &EquilibriumEnthalpyPhaseCurve) -> HeatedImpactSystem<'_> {
    with_routes(c, false, vec![], routes()).unwrap()
}

fn stress_integral(constant: f64, gradient: f64, lower: f64, upper: f64) -> f64 {
    // Exact integral of [constant + gradient*(z-h/2)]^2, evaluated
    // without subtracting nearby cubes or cancelling nonnegative terms.
    let width = upper - lower;
    let center = 0.5 * (lower + upper) - 0.5 * GAP;
    width * ((constant + gradient * center).powi(2) + gradient.powi(2) * width.powi(2) / 12.0)
}

#[test]
fn supplied_sliding_layer_heat_matches_stokes_work_and_closes_finite_enthalpy() {
    let c = curve(20.0);
    let conductance = WIDTH * GAP.powi(3) / (12.0 * MU * LENGTH);
    for analytic in [false, true] {
        let mut system = heated_driven(&c);
        system.set_analytic_newton(analytic);
        let mut native = driven_mechanics(false).prepare_analytic().unwrap();
        native.set_analytic_newton(analytic);
        let initial = system.mechanics().stored_energy_j();
        let mut supplied = 0.0;
        let mut expected_layers = [0.0; 2];
        for _ in 0..32 {
            let before = system.mechanics().state().to_vec();
            let frame = system
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let plain = native
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(system.mechanics().state(), native.state());
            assert_eq!(frame.mechanical.supplied_work_j, plain.supplied_work_j);
            assert_eq!(
                frame.mechanical.dissipated_energy_j,
                plain.dissipated_energy_j
            );
            let after = system.mechanics().state();
            let lower = 0.5 * (before[1] + after[1]) / M1.sqrt();
            let upper = 0.5 * (before[3] + after[3]) / M2.sqrt();
            let stress = MU * (lower - upper) / GAP;
            let gradient = PRESSURE / (2.0 * LENGTH);
            // Opposite channel orientations reverse the complete stress,
            // leaving identical layer heat. Cross terms differ by layer.
            let factor = DT * 2.0 * WIDTH * LENGTH / MU;
            let layer_heat = [
                factor * stress_integral(stress, gradient, 0.0, 0.5 * GAP),
                factor * stress_integral(stress, gradient, 0.5 * GAP, GAP),
            ];
            let expected_supply = DT
                * PRESSURE
                * (0.5 * conductance * PRESSURE + 0.5 * WIDTH * GAP * (lower + upper));
            let expected_heat = layer_heat[0] + layer_heat[1];
            near(frame.mechanical.supplied_work_j, expected_supply, 2e-12);
            near(frame.deposited_heat_j, expected_heat, 2e-12);
            near(frame.thermal_energy_change_j, expected_heat, 2e-12);
            near(frame.remaining_mechanical_loss_j, 0.0, 2e-12);
            near(frame.combined_energy_residual_j, 0.0, 2e-12);
            supplied += expected_supply;
            for i in 0..2 {
                expected_layers[i] += layer_heat[i];
                near(
                    system.reservoirs()[i].received_heat_j(),
                    expected_layers[i],
                    2e-12,
                );
            }
        }
        assert!(expected_layers[1] > expected_layers[0]);
        near(
            system.mechanics().stored_energy_j()
                + system
                    .reservoirs()
                    .iter()
                    .map(CaloricReservoir::storage_change_j)
                    .sum::<f64>(),
            initial + supplied,
            2e-11,
        );
    }
}

#[test]
fn passive_and_supplied_heat_indices_are_distinct_and_channel_coverage_is_exact() {
    let c = curve(20.0);
    let first = routes().remove(0);
    let mut mixed = with_routes(&c, true, vec![first.clone()], vec![first.clone()]).unwrap();
    let frame = mixed
        .step(&[0.0; 2], &CancelGate::new_clock_free())
        .unwrap();
    // Both islands contain two equal, oppositely oriented passages. Route
    // channel zero in BOTH namespaces and retain channel one's heat outside.
    near(
        frame.deposited_heat_j,
        0.5 * frame.mechanical.dissipated_energy_j,
        2e-12,
    );
    near(
        frame.remaining_mechanical_loss_j,
        frame.deposited_heat_j,
        2e-12,
    );
    for bad in [
        vec![first.clone(), first.clone()],
        vec![SlidingFilmLayerHeatRoute {
            film: 1,
            ..first.clone()
        }],
        vec![SlidingFilmLayerHeatRoute {
            channel: 2,
            ..first.clone()
        }],
        vec![SlidingFilmLayerHeatRoute {
            layer_destinations: vec!["lower".into()],
            ..first.clone()
        }],
        vec![SlidingFilmLayerHeatRoute {
            layer_destinations: vec!["lower".into(), "missing".into()],
            ..first.clone()
        }],
    ] {
        assert!(with_routes(&c, false, vec![], bad).is_err());
    }
}

#[test]
fn pressure_supply_heat_cancellation_and_late_thermal_refusal_are_atomic() {
    let c = curve(20.0);
    let mut system = heated_driven(&c);
    let before = system.mechanics().state().to_vec();
    let original_stores = system.reservoirs().to_vec();
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(system.step(&[0.0; 2], &stopped).is_err());
    assert_eq!(system.mechanics().state(), before);
    assert_eq!(system.mechanics().samples(), 0);
    assert_eq!(system.reservoirs(), original_stores);
    let retry = system
        .step(&[0.0; 2], &CancelGate::new_clock_free())
        .unwrap();
    let mut clean = heated_driven(&c);
    let expected = clean
        .step(&[0.0; 2], &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(system.mechanics().state(), clean.mechanics().state());
    assert_eq!(system.reservoirs(), clean.reservoirs());
    assert_eq!(
        retry.mechanical.supplied_work_j,
        expected.mechanical.supplied_work_j
    );
    assert_eq!(retry.deposited_heat_j, expected.deposited_heat_j);

    let narrow = curve(1e-12);
    let mut late = heated_driven(&narrow);
    let motion = late.mechanics().state().to_vec();
    let thermal_stores = late.reservoirs().to_vec();
    assert!(late.step(&[0.0; 2], &CancelGate::new_clock_free()).is_err());
    assert_eq!(late.mechanics().state(), motion);
    assert_eq!(late.mechanics().samples(), 0);
    assert_eq!(late.reservoirs(), thermal_stores);
}

#[test]
fn pressure_fed_layer_heat_enters_the_original_conducting_network_transaction() {
    let c = curve(20.0);
    let mut system = heated_driven(&c);
    let mut thermal = network(&c);
    let mut observer = ();
    let initial = system.mechanics().stored_energy_j();
    let mut supplied = 0.0;
    for tick in 1..=16 {
        let frame = system
            .step_with_thermal_network(
                &mut thermal,
                &mut observer,
                &bindings(),
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-12,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        supplied += frame.impact.thermal.impact.mechanical.supplied_work_j;
        near(frame.external_heat_j, 0.0, 2e-12);
        near(thermal.time_s(), tick as f64 * DT, 1e-14);
        for (i, binding) in bindings().iter().enumerate() {
            assert_eq!(
                system.reservoirs()[i].state().specific_enthalpy_j_kg(),
                thermal.cells()[&binding.cell].specific_h_j_kg,
            );
        }
    }
    assert!(
        thermal
            .cells()
            .values()
            .all(|cell| cell.boundary.source_w == 0.0)
    );
    near(
        system.mechanics().stored_energy_j()
            + system
                .reservoirs()
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>(),
        initial + supplied,
        2e-11,
    );
}
