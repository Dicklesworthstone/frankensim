//! Coherent geometric limits for held expansion, independent of DC response.

use fs_material::visco::TensorProny;
use fs_math::c64::C64;

fn law(limit: f64) -> TensorProny {
    let c = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 2.0 } else { 0.0 }));
    TensorProny::new_mandel(c, vec![], limit)
        .unwrap()
        .with_band((1.0, 10.0))
        .unwrap()
}

#[test]
fn held_free_expansion_has_a_geometric_bound_without_authorizing_dc_stress() {
    let law = law(0.05);
    let spectrum = law.prepare_harmonic_spectrum(&[3.0, 3.0, 3.0]).unwrap();
    let mut fields = [[C64::ZERO; 6]; 3];
    fields[0][0] = C64::new(1e100, -1e100);
    fields[1][0] = -fields[0][0];
    fields[2][0] = C64::new(0.002, -0.001);
    let mut free = [0.0; 6];
    free[0] = 0.004;
    spectrum
        .validate_kinematic_envelope_with_bias(&fields, &free, 0.01)
        .unwrap();
    // Actual mechanical strain still has its independent source-band gate.
    assert!(spectrum.validate_strains_with_bias(&fields, &free).is_err());
    spectrum
        .validate_strains_with_bias(&fields, &[0.0; 6])
        .unwrap();
    free[0] = 0.009;
    assert!(
        spectrum
            .validate_kinematic_envelope_with_bias(&fields, &free, 0.01)
            .is_err()
    );
}

#[test]
fn component_and_original_norm_limits_cover_the_complete_real_waveform() {
    let material = law(0.05);
    let source = material.prepare_harmonic_spectrum(&[3.0]).unwrap();
    source
        .validate_kinematic_envelope_with_bias(&[[C64::ZERO; 6]], &[0.006; 6], 0.007)
        .unwrap();
    let mut field = [C64::ZERO; 6];
    field[0] = C64::new(0.006, 0.006);
    // Its starting value is small, but its true phase peak exceeds the ceiling.
    assert!(
        source
            .validate_kinematic_envelope_with_bias(&[field], &[0.0; 6], 0.008)
            .is_err()
    );
    assert!(
        source
            .validate_kinematic_envelope_with_bias(&[field], &[0.0; 6], f64::NAN)
            .is_err()
    );
    assert!(
        source
            .validate_kinematic_envelope_with_bias(&[field], &[0.0; 6], 0.0)
            .is_err()
    );
    let narrow = law(0.012);
    let source = narrow.prepare_harmonic_spectrum(&[3.0]).unwrap();
    // Every component fits; the original material's complete norm still refuses.
    assert!(
        source
            .validate_kinematic_envelope_with_bias(&[[C64::ZERO; 6]], &[0.006; 6], 0.007)
            .is_err()
    );
}
