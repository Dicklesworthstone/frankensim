//! Continuing physical inventory, delegated kinetics and atomic refusal.
use super::super::super::{CrystallizationCaloric, CrystallizationConfig, CrystallizationSolve};
use super::*;

struct Branch {
    offset: f64,
    density: f64,
}
impl CrystallizationCaloric for Branch {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut identity = [0xC6; 32];
        identity[..8].copy_from_slice(&self.offset.to_bits().to_le_bytes());
        identity[8..16].copy_from_slice(&self.density.to_bits().to_le_bytes());
        Some(identity)
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
fn rate(t: f64) -> Result<f64, CrystallizationError> {
    Ok(0.002 * (500.0 - t).max(0.0))
}
fn law() -> CrystallizationLaw<'static> {
    CrystallizationLaw::new(
        &PARENT,
        &CRYSTAL,
        &rate,
        CrystallizationConfig {
            identity: [0xC7; 32],
            avrami_exponent: 1.0,
            maximum_crystalline_mass_fraction: 0.8,
            temperature_range_k: [250.0, 600.0],
            melting_temperature_k: 500.0,
        },
    )
    .unwrap()
}
fn controls() -> IsobaricCrystallizationControls {
    IsobaricCrystallizationControls {
        caloric: CrystallizationSolve::default(),
        energy_tolerance_j: 1e-6,
        volume_tolerance_m3: 1e-14,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

#[test]
fn continuing_adiabatic_parcel_accepts_contraction_work_and_later_remelting() {
    let law = law();
    let initial = law.state_at_temperature(380.0, 0.0).unwrap();
    let start =
        IsobaricCrystallizationParcel::new(&law, initial, 2.0, 0.002, 2e6, controls(), || false)
            .unwrap()
            .with_equilibrium_remelting(law.equivalent_age_remelting().unwrap())
            .unwrap();
    close(start.volume_m3(), 0.0025, 1e-18);
    assert_eq!(start.reference_volume_m3(), 0.002);
    let (grown, step) = start
        .advance_heat(
            0.0,
            2.0,
            controls(),
            CrystallizationAdaptivity::default(),
            || false,
        )
        .unwrap();
    let x = grown.state().crystalline_mass_fraction();
    assert!(x > 0.1);
    close(grown.state().temperature_k(), 380.0 + 100.0 * x, 1e-8);
    close(
        grown.volume_m3(),
        2.0 * ((1.0 - x) / 800.0 + x / 1000.0),
        2e-18,
    );
    assert_eq!(step.enthalpy_change_j, 0.0);
    close(step.boundary_work_j, -1000.0 * x, 2e-9);
    close(step.internal_energy_change_j, 1000.0 * x, 2e-9);
    close(
        grown.internal_energy_j().unwrap() - start.internal_energy_j().unwrap(),
        step.internal_energy_change_j,
        1e-8,
    );
    assert_eq!(grown.time_s(), 2.0);
    assert_eq!(grown.mass_kg(), 2.0);
    assert_eq!(grown.reference_volume_m3(), start.reference_volume_m3());
    assert!(step.adaptive.unwrap().accepted_substeps > 0);

    // The same original parcel heats to Tm and melts exactly half its crystals.
    let target_x = 0.5 * x;
    let target_h = 200_000.0 - 100_000.0 * target_x;
    let heat = 2.0 * (target_h - grown.state().specific_enthalpy_j_kg());
    let (melted, release) = grown.apply_heat(heat, controls(), || false).unwrap();
    close(melted.state().temperature_k(), 500.0, 1e-10);
    close(melted.state().crystalline_mass_fraction(), target_x, 2e-12);
    close(release.boundary_work_j, 500.0 * x, 2e-9);
    close(
        release.internal_energy_change_j + release.boundary_work_j,
        heat,
        1e-6,
    );
    assert_eq!(melted.time_s(), grown.time_s());
    assert!(melted.volume_m3() > grown.volume_m3());
    assert!(release.adaptive.is_none());
    let (held, stationary) = melted.apply_heat(0.0, controls(), || false).unwrap();
    assert_eq!(held.state(), melted.state());
    assert_eq!(held.volume_m3(), melted.volume_m3());
    assert_eq!(held.time_s(), melted.time_s());
    assert_eq!(stationary.boundary_work_j, 0.0);
    assert_eq!(start.state(), initial);
    assert_eq!(start.time_s(), 0.0);
}

#[test]
fn original_core_and_parcel_share_one_phase_path_and_replay() {
    let law = law();
    let initial = law.state_at_temperature(380.0, 0.0).unwrap();
    let start =
        IsobaricCrystallizationParcel::new(&law, initial, 2.0, 0.002, 2e6, controls(), || false)
            .unwrap();
    let heat = -100.0;
    let adaptive = CrystallizationAdaptivity::default();
    let direct = law
        .advance_specific_enthalpy_adaptive(
            initial,
            initial.specific_enthalpy_j_kg() + heat / 2.0,
            0.2,
            controls().caloric,
            adaptive,
            None,
            || false,
        )
        .unwrap();
    let (next, report) = start
        .advance_heat(heat, 0.2, controls(), adaptive, || false)
        .unwrap();
    let (replay, repeated) = start
        .advance_heat(heat, 0.2, controls(), adaptive, || false)
        .unwrap();
    assert_eq!(next.state(), direct.state);
    assert_eq!(report.adaptive, Some(direct));
    assert_eq!(next.state(), replay.state());
    assert_eq!(next.volume_m3(), replay.volume_m3());
    assert_eq!(next.time_s(), replay.time_s());
    assert_eq!(report, repeated);
    assert_eq!(next.remelting_model(), None);
    assert!(std::ptr::eq(next.law(), &law));
}

#[test]
fn late_cancellation_and_invalid_steps_preserve_all_accepted_inventory() {
    let law = law();
    let initial = law.state_at_temperature(380.0, 0.0).unwrap();
    let start =
        IsobaricCrystallizationParcel::new(&law, initial, 2.0, 0.002, 2e6, controls(), || false)
            .unwrap();
    let mut checkpoints = 0;
    start
        .advance_heat(
            0.0,
            0.2,
            controls(),
            CrystallizationAdaptivity::default(),
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut visited = 0;
    assert!(matches!(
        start.advance_heat(
            0.0,
            0.2,
            controls(),
            CrystallizationAdaptivity::default(),
            || {
                visited += 1;
                visited == checkpoints
            },
        ),
        Err(CrystallizationError::Cancelled)
    ));
    assert_eq!(visited, checkpoints);
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            start
                .advance_heat(
                    0.0,
                    dt,
                    controls(),
                    CrystallizationAdaptivity::default(),
                    || false
                )
                .is_err()
        );
    }
    for heat in [f64::INFINITY, f64::NAN, 1e12] {
        assert!(start.apply_heat(heat, controls(), || false).is_err());
    }
    let mut bad = controls();
    bad.energy_tolerance_j = -1.0;
    assert!(start.apply_heat(0.0, bad, || false).is_err());
    assert_eq!(start.state(), initial);
    assert_eq!(start.time_s(), 0.0);
    assert_eq!(start.mass_kg(), 2.0);
    assert_eq!(start.reference_volume_m3(), 0.002);
    close(start.volume_m3(), 0.0025, 1e-18);
    for reference in [0.0, -1.0, f64::INFINITY] {
        assert!(
            IsobaricCrystallizationParcel::new(
                &law,
                initial,
                2.0,
                reference,
                2e6,
                controls(),
                || false,
            )
            .is_err()
        );
    }
}
