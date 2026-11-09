//! Original-source local-temperature demixing, interfacial energy and atomic refusals.
//! Material tables and pair parameters are synthetic numerical fixtures.

use super::thermal::caloric;
use super::*;
use fs_flux::miscibility::adaptive::BinaryMiscibilityTimeControl;
use fs_flux::miscibility::material::nonisothermal::BinaryMiscibilityMaterialNonisothermalError as Error;
use fs_flux::miscibility::nonisothermal::adaptive::BinaryNonisothermalTimeControl;
use fs_flux::miscibility::nonisothermal::phase_field::{
    BinaryNonisothermalPhaseField, adaptive::BinaryNonisothermalPhaseFieldTimeControl,
};
use fs_flux::miscibility::nonisothermal::{BinaryNonisothermalFace, BinaryNonisothermalSolve};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn faces() -> Vec<BinaryNonisothermalFace> {
    vec![BinaryNonisothermalFace {
        thermal_conductivity_w_m_k: 0.03,
        heat_of_transport_j_per_mol: 0.0,
    }]
}

fn bracket() -> [Temperature; 2] {
    [Temperature::new(300.0), Temperature::new(600.0)]
}

fn time(step_s: f64) -> BinaryNonisothermalPhaseFieldTimeControl {
    BinaryNonisothermalPhaseFieldTimeControl {
        thermal: BinaryNonisothermalTimeControl {
            composition: BinaryMiscibilityTimeControl {
                initial_step_s: step_s,
                minimum_step_s: 1e-10,
                maximum_step_s: step_s,
                maximum_attempts: 1000,
                absolute_mole_fraction_tolerance: 1e-6,
                relative_mole_fraction_tolerance: 0.0,
            },
            absolute_temperature_tolerance_k: 1e-5,
            relative_temperature_tolerance: 0.0,
        },
        absolute_interfacial_energy_tolerance_j: 1e-7,
        relative_interfacial_energy_tolerance: 1e-6,
    }
}

// Integrate the original pure Cp curves independently. With this fixture's
// A/d=1 m and one physical edge, its complete interface energy is kappa*dx²/2.
fn energy(field: &BinaryNonisothermalPhaseField, kappa: f64) -> f64 {
    let bulk: f64 = field
        .component_moles()
        .iter()
        .zip(field.temperatures())
        .map(|(&[a, b], t)| {
            let delta = t.value() - 300.0;
            a * (40.0 * delta + delta * delta / 15.0)
                + b * (80.0 * delta + 2.0 * delta * delta / 15.0)
                + R * 900.0 * a * b / (a + b)
        })
        .sum();
    let difference = field.mole_fractions()[0][0] - field.mole_fractions()[1][0];
    bulk + 0.5 * kappa * difference * difference
}

#[test]
fn original_phase_field_sources_follow_local_states_and_conserve_total_energy() {
    let card = temperature_card(1.0, 600.0);
    let liquid = liquid(&card);
    let caloric = caloric(&liquid);
    let kappa = 200.0;
    let model = gradient(300.0, 10_000.0, kappa, 900.0);
    let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model.clone(), 1e-12).unwrap();
    let net = network(model, 1e-5);
    let initial = net
        .state_from_mole_fractions(vec![[0.4, 0.6], [0.6, 0.4]], || false)
        .unwrap();
    let resolved = net
        .resolve_liquid_field(&initial, &source, 2, || false)
        .unwrap();
    let local = resolved
        .bind_nonisothermal_phase_field(
            &caloric,
            vec![Temperature::new(320.0), Temperature::new(360.0)],
            faces(),
            bracket(),
            1 << 20,
            || false,
        )
        .unwrap();

    // Structurally identical laws are not the exact original liquid binding.
    let other_liquid = super::liquid(&card);
    let other_caloric = self::caloric(&other_liquid);
    assert!(matches!(
        resolved.bind_nonisothermal_phase_field(
            &other_caloric,
            vec![Temperature::new(320.0); 2],
            faces(),
            bracket(),
            1 << 20,
            || false,
        ),
        Err(Error::Invalid(_))
    ));

    let controls = BinaryNonisothermalSolve::default();
    let before = local.field().clone();
    let before_cells = local.cells().to_vec();
    let before_energy = energy(&before, kappa);
    let mut polls = 0;
    let candidate = local
        .step_candidate(0.001, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let end = candidate.field();
    assert_eq!(end.field().elapsed_s(), 0.001);
    assert!(std::ptr::eq(end.caloric_law(), &caloric));
    assert!(std::ptr::eq(end.source(), &source));
    assert_eq!(end.field().caloric_model(), caloric.model());
    assert_eq!(end.network().model(), source.gradient_model());
    assert_eq!(end.maximum_resolved_cells(), 2);
    assert_eq!(end.field().transport_faces(), faces());
    assert!(candidate.step().interfacial_energy_change_j.abs() > 1e-12);
    assert!((energy(end.field(), kappa) - before_energy).abs() <= 2e-8);
    assert!(candidate.step().energy_residual_j.abs() <= controls.energy_tolerance_j);
    assert!(candidate.step().entropy_change_j_per_k >= -controls.entropy_tolerance_j_per_k);
    for species in 0..2 {
        near(
            end.field()
                .component_moles()
                .iter()
                .map(|n| n[species])
                .sum(),
            1.0,
        );
    }
    for (i, cell) in end.cells().iter().enumerate() {
        let temperature = end.field().temperatures()[i];
        let amounts = end.field().component_moles()[i];
        let x = amounts[0] / (amounts[0] + amounts[1]);
        assert_eq!(cell.thermodynamics().temperature(), temperature);
        assert!(cell.thermodynamics().thermodynamic_factor() < 0.0);
        assert_eq!(cell.phase().component_moles(), amounts);
        near(cell.phase().volume_m3(), 1e-4);
        let transport = cell.phase().transport_state().unwrap();
        near(transport.density_kg_m3(), 1000.0 - 820.0 * x);
        near(
            transport.specific_heat_j_kg_k(),
            3500.0 + 500.0 * (temperature.value() - 300.0) / 300.0,
        );
        let material = cell.phase().liquid().material();
        assert_eq!(material.state_coordinate("T"), Some(temperature.value()));
        assert_eq!(material.state_coordinate("pressure"), Some(100_000.0));
        assert_eq!(
            material.card_identity(),
            liquid.reference_state().material().card_identity()
        );
    }
    let advanced = end
        .advance_candidate(0.002, time(0.001), controls, || false)
        .unwrap();
    assert_eq!(advanced.field().field().elapsed_s(), 0.001 + 0.002);
    assert!(advanced.advance().accepted_intervals > 0);
    assert!(advanced.advance().maximum_accepted_error_ratio <= 1.0);
    assert!(std::ptr::eq(advanced.field().caloric_law(), &caloric));
    assert!((energy(advanced.field().field(), kappa) - before_energy).abs() <= 4e-8);
    let (accepted, report) = advanced.into_parts();
    assert_eq!(report.duration_s, 0.002);
    for (i, cell) in accepted.cells().iter().enumerate() {
        assert_eq!(
            cell.thermodynamics().temperature(),
            accepted.field().temperatures()[i]
        );
    }
    for stop in [1, polls] {
        let mut calls = 0;
        assert!(
            local
                .step_candidate(0.001, controls, || {
                    calls += 1;
                    calls == stop
                })
                .is_err()
        );
        assert_eq!(local.field(), &before);
        assert_eq!(local.cells(), before_cells.as_slice());
    }
    let replay = local.step_candidate(0.001, controls, || false).unwrap();
    assert_eq!(replay.field().field(), candidate.field().field());
    assert_eq!(replay.field().cells(), candidate.field().cells());
    assert_eq!(replay.step(), candidate.step());
}

#[test]
fn phase_field_late_source_domain_or_volume_refusal_keeps_original_field() {
    for (density_scale, maximum_temperature) in [(0.9, 600.0), (1.0, 300.000001)] {
        let card = temperature_card(density_scale, maximum_temperature);
        let liquid = liquid(&card);
        let caloric = caloric(&liquid);
        let source = BinaryLiquidSquareGradientLaw::try_new(&liquid, model(), 1e-12).unwrap();
        let net = network(model(), 1e-5);
        let initial = net
            .state_from_mole_fractions(vec![[0.4, 0.6], [0.6, 0.4]], || false)
            .unwrap();
        let resolved = net
            .resolve_liquid_field(&initial, &source, 2, || false)
            .unwrap();
        let local = resolved
            .bind_nonisothermal_phase_field(
                &caloric,
                vec![Temperature::new(300.0); 2],
                faces(),
                bracket(),
                1 << 20,
                || false,
            )
            .unwrap();
        let controls = BinaryNonisothermalSolve::default();
        let numerical = local
            .field()
            .step_candidate(0.1, controls, || false)
            .unwrap();
        assert!(
            numerical
                .field()
                .temperatures()
                .iter()
                .any(|t| t.value() > 300.000001)
        );
        assert!(matches!(
            local.step_candidate(0.1, controls, || false),
            Err(Error::Material(
                BinaryMiscibilityMaterialError::Source { .. }
            ))
        ));
        assert!(matches!(
            local.advance_candidate(0.1, time(0.01), controls, || false),
            Err(Error::Material(
                BinaryMiscibilityMaterialError::Source { .. }
            ))
        ));
        assert_eq!(local.field().elapsed_s(), 0.0);
        assert_eq!(local.field().component_moles(), initial.component_moles());
        assert_eq!(local.cells(), resolved.cells());
    }
}
