use super::super::super::loss_amplitude;
use super::*;

#[test]
fn physical_heat_retains_source_admittance_and_loss_across_extreme_scales() {
    let cases = [
        (1e200, 0.5, C64::from_re(1e-200), 2.5e-201),
        (1e-200, 0.5, C64::from_re(1e200), 2.5e199),
        (1e-200, 1e-200, C64::from_re(1e100), 5e-201),
        (1e200, 0.5, C64::new(1e-200, -1e-200), 5e-201),
        (1e200, -0.5, C64::from_re(1e-200), -2.5e-201),
    ];
    for (y, a, source, expected) in cases {
        let actual = source_power(y, a, source).to_f64();
        assert!(
            (actual / expected - 1.0).abs() < 8e-15,
            "{actual:e} != {expected:e}"
        );
    }
    assert_eq!(source_power(1e200, 0.0, C64::from_re(1e200)).to_f64(), 0.0);
    assert_eq!(source_power(1e200, 0.5, C64::ZERO).to_f64(), 0.0);
}

#[test]
fn smallest_original_loss_survives_a_large_coherent_field() {
    let loss = f64::from_bits(1);
    let field = C64::from_re(2.0_f64.powi(500));
    let heat = loss_amplitude(field, loss).norm_sq();
    // (1/2)*2^-1074*(2^500)^2 = 2^-75, independently representable.
    assert!((heat / 2.0_f64.powi(-75) - 1.0).abs() < 5e-15);
    assert_eq!(loss_amplitude(field, 0.0), C64::ZERO);
}
