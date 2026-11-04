//! G1/G3: original contact laws and line measures survive both moving-end clocks.
use super::*;
use crate::thin_plate::CompactBody;
use fs_scenario::{BendingRelaxationBranch, PronyBending};

fn obstacle() -> UnilateralObstacle {
    UnilateralObstacle {
        stations: vec![0.35],
        gaps_m: vec![0.0001],
        measures_m: vec![0.25],
        stiffness: 1e6,
        alpha: 2.0,
        mu_kinetic: 0.3,
        mu_static: 0.4,
        friction_regularization_m_s: Some(0.01),
        internal_loss: 2.0,
        provenance: "moving-end physical line contact".into(),
    }
}

fn render(
    obstacles: &[UnilateralObstacle],
    body: bool,
    prony: bool,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let string = PrestressedString {
        length_m: 0.5,
        tension_n: 64.0,
        lin_density_kg_m: 1.0,
        axial_stiffness_n: 0.0,
        width_m: 0.001,
        n_modes: 2,
        damping_ratio: 0.0,
        rayleigh: None,
        bending_stiffness_n_m2: 1e-3,
        kelvin_voigt_bending: None,
        relaxation_bending: prony.then_some(PronyBending {
            branches: vec![BendingRelaxationBranch {
                relaxing_stiffness_n_m2: 1e-3,
                relaxation_time_s: 0.1,
            }],
            omega_band_rad_s: (0.0, 1e5),
            material_state_identity: None,
            source_properties: None,
        }),
        polarization_detune: 0.0,
        moving_end: true,
        tension_loss_factor: None,
        bending_loss_factor: None,
    };
    let mut bank = PlateBank::default();
    if body {
        bank.linear.push(CompactBody::from_radiator(RadiatingPlate {
            area_m2: 0.01,
            mass_kg: 0.2,
            frequency_hz: 10.0,
            damping_ratio: 0.02,
        })?);
    }
    realize_dirac_join(
        &string,
        Some(Pluck {
            station_frac: 0.5,
            height_m: 0.001,
        }),
        None,
        obstacles,
        None,
        &mut bank,
        None,
        None,
        None,
        &gas_state(AmbientGas::sea_level())?,
        1.0,
        48_000,
        96,
    )
}

#[test]
fn split_line_measure_preserves_ode_and_descriptor_audio_with_material_history() {
    let whole = obstacle();
    let mut split = whole.clone();
    split.stations = vec![whole.stations[0]; 2];
    split.gaps_m = vec![whole.gaps_m[0]; 2];
    split.measures_m = vec![whole.measures_m[0] / 2.0; 2];
    for body in [false, true] {
        for prony in [false, true] {
            let a = render(&[whole.clone()], body, prony).unwrap();
            let b = render(&[split.clone()], body, prony).unwrap();
            let peak = a.iter().map(|p| p.abs()).fold(0.0_f64, f64::max);
            let error = a
                .iter()
                .zip(&b)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            assert!(peak > 0.0 && peak.is_finite());
            assert!(
                error <= 5e-10 * peak,
                "body={body} prony={prony}: {error} vs {peak}"
            );
            let mut frictionless = whole.clone();
            frictionless.mu_kinetic = 0.0;
            let dry = render(&[frictionless], body, prony).unwrap();
            let difference = a
                .iter()
                .zip(&dry)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            assert!(
                difference > 1e-9 * peak,
                "friction must reach the original force ports"
            );
        }
    }
}

#[test]
fn both_moving_end_clocks_reject_missing_measures_and_out_of_span_sites() {
    for body in [false, true] {
        let mut invalid = obstacle();
        invalid.measures_m.clear();
        assert!(
            render(&[invalid], body, false)
                .unwrap_err()
                .to_string()
                .contains("physical line measures")
        );
        let mut invalid = obstacle();
        invalid.stations[0] = 1.1;
        assert!(
            render(&[invalid], body, false)
                .unwrap_err()
                .to_string()
                .contains("strictly inside")
        );
    }
}
