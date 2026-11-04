//! Direct two-point observations retain both the envelope and its carrier.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumCoherence, AngularSpectrumCoherenceLimits, AngularSpectrumCorrelationLimits,
    AngularSpectrumObservation,
};

#[test]
fn g1_carrier_correlation_retains_small_envelope_bins_and_relative_optical_phase() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [2.0_f64.powi(60), 2.0_f64.powi(-60)],
    };
    let factors = [PlanarAngularSpectrum::prepare_with_carrier(
        &stack,
        grid,
        [0.4, 0.0],
        &[[C64::ZERO, C64::from_re(2.0)], [C64::ZERO; 2]],
        limits(),
        || false,
    )
    .unwrap()];
    let bank = AngularSpectrumCoherence::try_new(
        &factors[0],
        &factors,
        AngularSpectrumCoherenceLimits {
            maximum_modes: 1,
            maximum_point_mode_pairs: 2,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: limits().maximum_fft_work,
        },
        || false,
    )
    .unwrap();
    let query_limits = AngularSpectrumCorrelationLimits {
        maximum_mode_evaluations: bank.correlation_mode_evaluations().unwrap(),
        maximum_workspace_bytes: bank.correlation_workspace_bytes(),
    };
    let at = |x| AngularSpectrumObservation {
        side: Side::Incident,
        distance_from_face_m: 0.0,
        position_m: [x, 0.0],
    };
    let a = at(0.0);
    let dark = bank
        .mutual_coherence(a, at(grid.spacing_m[0]), query_limits, || false)
        .unwrap();
    // E_y = (1 + exp(+i*pi*x/dx))*exp(-i*0.4*x). The negative
    // Nyquist offset is below one ulp of the carrier, but its fringe remains.
    // The exact-source direct sum may retain the tiny binary64 pi residual.
    let sampled = factors[0].field(Side::Incident, 0.0, || false).unwrap();
    near(dark.electric_at_a_v2_m2()[1][1].re, 2.0, 2e-13);
    near(
        dark.electric_at_b_v2_m2()[1][1].re,
        0.5 * sampled.electric_v_per_m()[1][1].norm_sq(),
        1e-28,
    );
    complex_near(dark.electric_cross_v2_m2()[1][1], C64::ZERO, 2e-14);

    // The envelope variation over 0.5 m is below 2e-18 rad. The physical
    // cross moment still has the resolved carrier phase Gamma_ab=2*exp(+i*.2).
    let b = at(0.5);
    let mutual = bank.mutual_coherence(a, b, query_limits, || false).unwrap();
    complex_near(
        mutual.electric_cross_v2_m2()[1][1],
        phase(0.2).scale(2.0),
        2e-13,
    );
    let ey = [C64::ZERO, C64::ONE, C64::ZERO];
    let combined = bank
        .recombined_electric_mean_square_v2_m2(a, b, [ey; 2], 0.7, query_limits, || false)
        .unwrap();
    near(combined, 4.0 * (1.0 + det::cos(0.5)), 3e-13);
}
