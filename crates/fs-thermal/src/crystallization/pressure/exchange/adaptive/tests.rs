//! G1 finite-body cooling oracle, latent accuracy, and atomic retry behavior.
use super::super::super::super::{
    CrystallizationCaloric, CrystallizationConfig, CrystallizationLaw,
};
use super::*;

struct Branch {
    offset: f64,
    density: f64,
}
impl CrystallizationCaloric for Branch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut id = [0xE1; 32];
        id[..8].copy_from_slice(&self.offset.to_bits().to_le_bytes());
        id[8..16].copy_from_slice(&self.density.to_bits().to_le_bytes());
        Some(id)
    }
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, CrystallizationError> {
        Ok(1000.0 * (t - 300.0) + self.offset)
    }
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(1000.0)
    }
    fn density_kg_m3(&self, _: f64) -> Result<Option<f64>, CrystallizationError> {
        Ok(Some(self.density))
    }
}
static PARENT: Branch = Branch {
    offset: 0.0,
    density: 800.0,
};
static CRYSTAL: Branch = Branch {
    offset: -100_000.0,
    density: 1000.0,
};
fn zero_rate(_: f64) -> Result<f64, CrystallizationError> {
    Ok(0.0)
}
fn law(melting: f64) -> CrystallizationLaw<'static> {
    CrystallizationLaw::new(
        &PARENT,
        &CRYSTAL,
        &zero_rate,
        CrystallizationConfig {
            identity: [0xE2; 32],
            avrami_exponent: 1.0,
            maximum_crystalline_mass_fraction: 0.8,
            temperature_range_k: [250.0, 800.0],
            melting_temperature_k: melting,
        },
    )
    .unwrap()
}
fn side() -> IsobaricExchangeSideControls {
    let mut controls = IsobaricExchangeSideControls::default();
    controls.caloric.energy_tolerance_j = 1e-6;
    controls
}
fn exchange() -> IsobaricExchangeControls {
    IsobaricExchangeControls {
        heat_bracket_j: [0.0, 180_000.0],
        maximum_evaluations: 100,
        residual_tolerance_j: 1e-6,
        energy_tolerance_j: 1e-5,
    }
}
fn parcel<'a>(
    law: &'a CrystallizationLaw<'static>,
    t: f64,
    exposure: f64,
    mass: f64,
) -> IsobaricCrystallizationParcel<'a, 'static> {
    IsobaricCrystallizationParcel::new(
        law,
        law.state_at_temperature(t, exposure).unwrap(),
        mass,
        mass / 1000.0,
        2e6,
        side().caloric,
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn controls(tolerance: f64) -> IsobaricExchangeAdaptivity {
    IsobaricExchangeAdaptivity {
        initial_dt_s: 2.0,
        maximum_dt_s: 2.0,
        temperature_tolerance_k: tolerance,
        specific_enthalpy_tolerance_j_kg: 1000.0 * tolerance,
        relative_tolerance: 0.0,
        energy_tolerance_j: 1e-3,
        ..IsobaricExchangeAdaptivity::default()
    }
}

#[test]
fn tighter_contact_time_tolerance_converges_to_two_body_exponential() {
    let law = law(700.0);
    let left = parcel(&law, 500.0, 0.0, 2.0);
    let right = parcel(&law, 300.0, 0.0, 3.0);
    let (loose_left, _, loose) = advance_exchange_to(
        &left,
        &right,
        500.0,
        2.0,
        side(),
        side(),
        exchange(),
        controls(0.5),
        || false,
    )
    .unwrap();
    let (tight_left, tight_right, tight) = advance_exchange_to(
        &left,
        &right,
        500.0,
        2.0,
        side(),
        side(),
        exchange(),
        controls(0.01),
        || false,
    )
    .unwrap();
    let exact_left = 380.0 + 120.0 * (-500.0_f64 * (1.0 / 2000.0 + 1.0 / 3000.0) * 2.0).exp();
    let loose_error = (loose_left.state().temperature_k() - exact_left).abs();
    let tight_error = (tight_left.state().temperature_k() - exact_left).abs();
    assert!(
        tight_error < loose_error / 3.0,
        "loose {loose_error}, tight {tight_error}"
    );
    assert!(tight_error < 0.3);
    assert!(tight.accepted_substeps > loose.accepted_substeps);
    assert!(tight.rejected_trials > 0);
    assert!(tight.maximum_accepted_error.normalized <= 1.0);
    assert_eq!(tight_left.time_s(), 2.0);
    assert_eq!(tight_right.time_s(), 2.0);
    close(
        2000.0 * tight_left.state().temperature_k() + 3000.0 * tight_right.state().temperature_k(),
        1_900_000.0,
        1e-6,
    );
    close(
        tight.left.supplied_heat_j,
        -tight.right.supplied_heat_j,
        0.0,
    );
    close(
        tight.heat_left_to_right_j,
        -tight.left.enthalpy_change_j,
        1e-6,
    );
    assert!(tight.energy_residual_j.abs() <= tight.absolute_energy_residual_j);
    assert!(tight.absolute_energy_residual_j <= controls(0.01).energy_tolerance_j);
    assert!(tight.evaluated_candidates >= 3 * tight.attempted_trials);
    let (replay_left, replay_right, replay) = advance_exchange_to(
        &left,
        &right,
        500.0,
        2.0,
        side(),
        side(),
        exchange(),
        controls(0.01),
        || false,
    )
    .unwrap();
    assert_eq!(replay_left.state(), tight_left.state());
    assert_eq!(replay_right.state(), tight_right.state());
    assert_eq!(replay, tight);
    // Comparison steps and rejected heat never debit either accepted input.
    assert_eq!(left.state().temperature_k(), 500.0);
    assert_eq!(right.state().temperature_k(), 300.0);
    assert_eq!(left.time_s(), 0.0);
}

#[test]
fn latent_inventory_controls_refinement_even_with_loose_temperature_tolerance() {
    let hot_law = law(700.0);
    let melting_law = law(500.0);
    let hot = parcel(&hot_law, 560.0, 0.0, 2.0);
    let cold = parcel(&melting_law, 500.0, 2.0_f64.ln(), 2.0)
        .with_equilibrium_remelting(melting_law.equivalent_age_remelting().unwrap())
        .unwrap();
    close(cold.state().crystalline_mass_fraction(), 0.4, 1e-14);
    let contact = IsobaricExchangeControls {
        heat_bracket_j: [0.0, 40_000.0],
        ..exchange()
    };
    let accuracy = IsobaricExchangeAdaptivity {
        initial_dt_s: 0.5,
        maximum_dt_s: 0.5,
        temperature_tolerance_k: 1e6,
        specific_enthalpy_tolerance_j_kg: 1e6,
        crystalline_mass_fraction_tolerance: 1e-5,
        exposure_tolerance: 1.0,
        relative_tolerance: 0.0,
        energy_tolerance_j: 1e-3,
        ..IsobaricExchangeAdaptivity::default()
    };
    let (cooled, melted, report) = advance_exchange_to(
        &hot,
        &cold,
        400.0,
        0.5,
        side(),
        side(),
        contact,
        accuracy,
        || false,
    )
    .unwrap();
    assert!(report.rejected_trials > 0);
    assert!(report.maximum_accepted_error.crystalline_mass_fraction <= 1e-5);
    assert!(cooled.state().temperature_k() < hot.state().temperature_k());
    close(melted.state().temperature_k(), 500.0, 1e-10);
    assert!(melted.state().crystalline_mass_fraction() < 0.4);
    close(
        report.heat_left_to_right_j,
        200_000.0 * (0.4 - melted.state().crystalline_mass_fraction()),
        1e-5,
    );
    assert!(report.right.boundary_work_j > 0.0);
    close(
        report.right.internal_energy_change_j + report.right.boundary_work_j,
        report.heat_left_to_right_j,
        1e-5,
    );
    close(
        report.left.enthalpy_change_j + report.right.enthalpy_change_j,
        0.0,
        1e-5,
    );
}

#[test]
fn comparison_budget_and_last_checkpoint_refuse_without_partial_publication() {
    let law = law(700.0);
    let left = parcel(&law, 500.0, 0.0, 2.0);
    let right = parcel(&law, 300.0, 0.0, 3.0);
    let accuracy = controls(0.01);
    assert!(matches!(
        advance_exchange_to(
            &left,
            &right,
            500.0,
            2.0,
            side(),
            side(),
            exchange(),
            IsobaricExchangeAdaptivity {
                maximum_trials: 1,
                ..accuracy
            },
            || false,
        ),
        Err(CrystallizationError::NonConverged)
    ));
    let mut calls = 0;
    advance_exchange_to(
        &left,
        &right,
        500.0,
        2.0,
        side(),
        side(),
        exchange(),
        accuracy,
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut visited = 0;
    assert!(matches!(
        advance_exchange_to(
            &left,
            &right,
            500.0,
            2.0,
            side(),
            side(),
            exchange(),
            accuracy,
            || {
                visited += 1;
                visited == calls
            },
        ),
        Err(CrystallizationError::Cancelled)
    ));
    assert_eq!(visited, calls);
    assert_eq!(left.time_s(), 0.0);
    assert_eq!(right.time_s(), 0.0);
    assert_eq!(left.volume_m3(), 0.0025);
    assert_eq!(right.state().temperature_k(), 300.0);
}
