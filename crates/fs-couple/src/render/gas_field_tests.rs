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
        assert_eq!(voice.sound_speed.to_bits(), inlet.sound_speed.to_bits());
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
