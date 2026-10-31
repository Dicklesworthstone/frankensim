use super::*;
use crate::driving_point::regional::field::tests::{cylinder, field};

#[test]
fn gas_field_reaches_rendered_pressure_and_exact_inlet_source() {
    let geometry = cylinder();
    let uniform = field(300.0);
    let gradient = field(270.0);
    let inlet = *uniform.state_at(0.0).unwrap().state();
    let reed = BeatingReed {
        rest_opening_m: 4e-4,
        width_m: 0.013,
        closing_pressure_pa: 6000.0,
        blowing_pressure_pa: 2800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
        damping_ratio: 0.35,
    };
    let make = |field: &fs_material::gas_field::AxialGasField| {
        ReedBoreVoice::new_with_gas_field(
            &geometry,
            field,
            8,
            reed,
            Termination::UnflangedOpen,
            PlateBank::default(),
            1.0,
            8000,
            2048,
            None,
        )
        .unwrap()
    };
    let mut a = make(&uniform);
    let mut b = make(&gradient);
    let mut original = ReedBoreVoice::new(
        &geometry,
        &inlet,
        reed,
        Termination::UnflangedOpen,
        PlateBank::default(),
        1.0,
        8000,
        2048,
        None,
    )
    .unwrap();
    let area = core::f64::consts::PI * geometry.segments[0].inlet_radius().powi(2);
    for voice in [&a, &b] {
        assert_eq!(voice.rho.to_bits(), inlet.density.to_bits());
        assert_eq!(voice.inlet_gas.sound_speed.to_bits(), inlet.sound_speed.to_bits());
        assert_eq!(voice.zc, inlet.density * inlet.sound_speed / area);
    }
    let mut unchanged = [0.0; 2048];
    let mut baseline = [0.0; 2048];
    let mut changed = [0.0; 2048];
    a.step_block(&mut unchanged).unwrap();
    original.step_block(&mut baseline).unwrap();
    b.step_block(&mut changed).unwrap();
    assert_eq!(unchanged, baseline);
    assert!(changed.iter().all(|value| value.is_finite()));
    let energy: f64 = unchanged.iter().map(|value| value * value).sum();
    let difference: f64 = unchanged
        .iter()
        .zip(changed)
        .map(|(a, b)| (a - b).powi(2))
        .sum();
    assert!(energy > 0.0, "source must actually sound");
    assert!(
        difference > 1e-6 * energy,
        "the downstream thermal field must reach audio feedback"
    );
    assert_eq!(a.samples_rendered(), 2048);
    assert_eq!(b.samples_rendered(), 2048);
}

#[test]
fn continuous_smooth_and_fingered_bores_carry_material_gradients_to_rendered_pressure() {
    let smooth = fs_duct::Duct {
        segments: vec![fs_duct::Segment::Profile {
            profile: std::sync::Arc::new(
                fs_duct::RadiusProfile::exponential(0.0073, 0.6, 0.57).unwrap(),
            ),
            steps: fs_duct::RiccatiSteps(32),
        }],
    };
    let fingered = fs_duct::Duct {
        segments: vec![
            fs_duct::Segment::Cylinder {
                radius: 0.0073,
                length: 0.285,
            },
            fs_duct::Segment::ToneHole {
                hole_radius: 0.002,
                chimney_height: 0.005,
                bore_radius: 0.0073,
                state: fs_duct::HoleState::Open,
            },
            fs_duct::Segment::Cylinder {
                radius: 0.0073,
                length: 0.285,
            },
        ],
    };
    for geometry in [smooth, fingered] {
        check_continuous_bore_render(&geometry);
    }
}

fn check_continuous_bore_render(geometry: &fs_duct::Duct) {
    let uniform = field(300.0);
    let gradient = field(275.0);
    let inlet = *uniform.state_at(0.0).unwrap().state();
    let reed = BeatingReed {
        rest_opening_m: 4e-4,
        width_m: 0.013,
        closing_pressure_pa: 6000.0,
        blowing_pressure_pa: 2800.0,
        attack_s: 0.008,
        mass_kg: 0.0,
        stiffness_n_m: 0.0,
        damping_ratio: 0.35,
    };
    let make = |field: &fs_material::gas_field::AxialGasField, steps| {
        ReedBoreVoice::new_with_continuous_gas_field(
            geometry, field, fs_duct::RiccatiSteps(steps), reed,
            Termination::UnflangedOpen, PlateBank::default(), 1.0, 8000, 2048, None,
        ).unwrap()
    };
    let mut a = make(&uniform, 16);
    let mut b = make(&gradient, 16);
    let mut c = make(&gradient, 32);
    let mut original = ReedBoreVoice::new(
        geometry, &inlet, reed, Termination::UnflangedOpen,
        PlateBank::default(), 1.0, 8000, 2048, None,
    ).unwrap();
    for voice in [&a, &b, &c] {
        assert_eq!(voice.rho.to_bits(), inlet.density.to_bits());
        assert_eq!(voice.inlet_gas.sound_speed.to_bits(), inlet.sound_speed.to_bits());
    }
    let mut baseline = [0.0; 2048];
    let mut unchanged = [0.0; 2048];
    let mut changed = [0.0; 2048];
    let mut refined = [0.0; 2048];
    original.step_block(&mut baseline).unwrap();
    a.step_block(&mut unchanged).unwrap();
    b.step_block(&mut changed).unwrap();
    c.step_block(&mut refined).unwrap();
    assert_eq!(baseline, unchanged);
    assert!(changed.iter().chain(refined.iter()).all(|value| value.is_finite()));
    let energy: f64 = baseline.iter().map(|value| value * value).sum();
    let difference: f64 = baseline.iter().zip(changed.iter()).map(|(a, b)| (a - b).powi(2)).sum();
    let refinement: f64 = changed.iter().zip(refined.iter()).map(|(a, b)| (a - b).powi(2)).sum();
    assert!(energy > 0.0, "the material-bound source must actually sound");
    assert!(difference > 1e-6 * energy, "downstream gas must alter audible feedback");
    assert!(refinement < difference * 0.01, "material effect must exceed its discretization change");
}
