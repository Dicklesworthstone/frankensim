use super::*;
use crate::healing::{HealingLaw, MobilityKnot, MobilityTable};

fn mobility() -> MobilityTable {
    MobilityTable::new(
        "declared-polymer-mobility".into(),
        vec![
            MobilityKnot {
                temperature_k: 300.0,
                rate_per_s: 0.0,
            },
            MobilityKnot {
                temperature_k: 400.0,
                rate_per_s: 0.0,
            },
            MobilityKnot {
                temperature_k: 500.0,
                rate_per_s: 1.0,
            },
        ],
    )
    .unwrap()
}

fn config() -> MixedModeStrengthConfig {
    MixedModeStrengthConfig {
        normal_stiffness_pa_m: 2e9,
        tangential_stiffness_pa_m: 5e8,
        bk_exponent: 2.0,
        compression_stiffness_pa_m: 4e9,
        maximum_separation_m: 0.01,
    }
}

fn knots() -> Vec<MixedModeStrengthKnot> {
    vec![
        MixedModeStrengthKnot {
            dose: 0.0,
            normal_strength_pa: 0.0,
            tangential_strength_pa: 0.0,
            mode_i_fracture_energy_j_m2: 0.0,
            mode_ii_fracture_energy_j_m2: 0.0,
        },
        MixedModeStrengthKnot {
            dose: 1.0,
            normal_strength_pa: 2e5,
            tangential_strength_pa: 3e5,
            mode_i_fracture_energy_j_m2: 100.0,
            mode_ii_fracture_energy_j_m2: 240.0,
        },
    ]
}

fn map(knots: Vec<MixedModeStrengthKnot>) -> Result<StrengthMap, HealingError> {
    StrengthMap::new_mixed_mode(
        "declared-mode-strengths".into(),
        "declared-polymer-mobility".into(),
        knots,
        config(),
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn exposure_selects_distinct_opening_and_shear_fracture_in_original_points() {
    let mobility = mobility();
    let map = map(knots()).unwrap();
    let virgin = ContactHistory::new(&mobility).unwrap();
    assert!(map.freeze(&virgin).unwrap().law.is_none());
    let exposed = virgin.advance(500.0, 500.0, 0.5).unwrap();
    let frozen = map.freeze(&exposed).unwrap();
    assert_eq!(frozen.map_identity, "declared-mode-strengths");
    assert_eq!(frozen.mobility_identity, mobility.identity());
    close(frozen.dose, 0.5);
    let law = frozen.law.unwrap();
    let p = law.mixed_mode_parameters().unwrap();
    close(p.normal_strength_pa, 1e5);
    close(p.tangential_strength_pa, 1.5e5);
    close(p.mode_i_fracture_energy_j_m2, 50.0);
    close(p.mode_ii_fracture_energy_j_m2, 120.0);
    for (axis, stiffness, strength, energy) in [
        (0, 2e9, 1e5, 50.0),
        (1, 5e8, 1.5e5, 120.0),
        (2, 5e8, 1.5e5, 120.0),
    ] {
        let mut at_peak = [0.0; 3];
        at_peak[axis] = strength / stiffness;
        let peak = law.bonded_point().trial(at_peak, 1e-9).unwrap();
        close(peak.traction_pa[axis], strength);
        let mut failed = [0.0; 3];
        failed[axis] = 1.01 * 2.0 * energy / strength;
        let failure = peak.point.trial(failed, 1e-9).unwrap();
        assert!(failure.point.failed());
        close(peak.energy.work + failure.energy.work, energy);
        close(failure.point.dissipated_j_m2(), energy);
        let reclosed = failure.point.trial([-1e-5, 0.0, 0.0], 1e-9).unwrap();
        assert!(reclosed.point.failed());
        close(reclosed.traction_pa[0], -4e4);
        close(reclosed.point.dissipated_j_m2(), energy);
    }
    assert_eq!(virgin.dose(), 0.0);
    assert_eq!(exposed.dose(), 0.5);
}

#[test]
fn frozen_bond_retains_original_bk_work_and_in_plane_rotation() {
    let mobility = mobility();
    let history = ContactHistory::new(&mobility)
        .unwrap()
        .advance(500.0, 500.0, 1.0)
        .unwrap();
    let law = map(knots()).unwrap().freeze(&history).unwrap().law.unwrap();
    // Independent proportional-loading calculation in physical separation.
    let direction = [0.6_f64, 0.8_f64, 0.0];
    let radial_stiffness = 2e9 * direction[0].powi(2) + 5e8 * direction[1].powi(2);
    let beta = 5e8 * direction[1].powi(2) / radial_stiffness;
    let expected_work = 100.0 + 140.0 * beta.powi(2);
    let onset =
        1.0 / ((2e9 * direction[0] / 2e5).powi(2) + (5e8 * direction[1] / 3e5).powi(2)).sqrt();
    let failure = 2.0 * expected_work / (radial_stiffness * onset);
    for e in [direction, [0.6, 0.0, -0.8], [0.6, 0.48, 0.64]] {
        let step = law
            .bonded_point()
            .trial(e.map(|v| 1.01 * failure * v), 1e-9)
            .unwrap();
        assert!(step.point.failed());
        close(step.energy.work, expected_work);
        close(
            step.point
                .mixed_mode_history()
                .unwrap()
                .energy_shear_fraction,
            beta,
        );
    }
}

#[test]
fn declared_induction_and_cold_pauses_preserve_the_original_exposure() {
    let mobility = mobility();
    let mut data = knots();
    data.insert(
        1,
        MixedModeStrengthKnot {
            dose: 0.25,
            ..data[0]
        },
    );
    let map = map(data).unwrap();
    let initial = ContactHistory::new(&mobility).unwrap();
    let induction = initial.advance(500.0, 500.0, 0.25).unwrap();
    assert!(map.freeze(&induction).unwrap().law.is_none());
    let hot = induction.advance(500.0, 500.0, 0.25).unwrap();
    let paused = hot.advance(350.0, 350.0, 10.0).unwrap();
    let law = map.freeze(&paused).unwrap().law.unwrap();
    let p = law.mixed_mode_parameters().unwrap();
    close(p.normal_strength_pa, 2e5 / 3.0);
    close(p.tangential_strength_pa, 1e5);
    close(p.mode_i_fracture_energy_j_m2, 100.0 / 3.0);
    close(p.mode_ii_fracture_energy_j_m2, 80.0);
    assert_eq!(map.freeze(&hot).unwrap().law.as_ref(), Some(&law));
    let reheated = paused.advance(500.0, 500.0, 0.5).unwrap();
    close(
        map.freeze(&reheated).unwrap().law.unwrap().strength_pa(),
        2e5,
    );
}

#[test]
fn incomplete_modes_and_inadmissible_maps_refuse() {
    for column in 0..4 {
        let mut data = knots();
        match column {
            0 => data[1].normal_strength_pa = 0.0,
            1 => data[1].tangential_strength_pa = 0.0,
            2 => data[1].mode_i_fracture_energy_j_m2 = 0.0,
            _ => data[1].mode_ii_fracture_energy_j_m2 = 0.0,
        }
        assert!(map(data).is_err());
    }
    let mut decreasing = knots();
    decreasing.push(MixedModeStrengthKnot {
        dose: 2.0,
        mode_ii_fracture_energy_j_m2: 230.0,
        ..decreasing[1]
    });
    assert!(map(decreasing).is_err());
    let mut impossible = knots();
    impossible[1].mode_i_fracture_energy_j_m2 = 1.0;
    assert!(map(impossible).is_err());
    for invalid in [
        MixedModeStrengthConfig {
            normal_stiffness_pa_m: f64::NAN,
            ..config()
        },
        MixedModeStrengthConfig {
            bk_exponent: 0.5,
            ..config()
        },
        MixedModeStrengthConfig {
            maximum_separation_m: 1e-6,
            ..config()
        },
    ] {
        assert!(
            StrengthMap::new_mixed_mode("map".into(), "mobility".into(), knots(), invalid).is_err()
        );
    }
}

#[test]
fn source_mismatch_and_overexposure_refuse_without_changing_history() {
    let mobility = mobility();
    let map = map(knots()).unwrap();
    let overexposed = ContactHistory::new(&mobility)
        .unwrap()
        .advance(500.0, 500.0, 2.0)
        .unwrap();
    assert_eq!(
        map.freeze(&overexposed).unwrap_err(),
        HealingError::Dose(2.0)
    );
    assert_eq!(overexposed.dose(), 2.0);
    let foreign = MobilityTable::new("other-polymer".into(), mobility.knots().to_vec()).unwrap();
    assert!(map.freeze(&ContactHistory::new(&foreign).unwrap()).is_err());
}
