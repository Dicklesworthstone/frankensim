//! Steinway-D-based string and hammer parameters which actually drive mechanics.
//!
//! Chabassier & Durufle, Physical parameters for piano modeling, INRIA RT-0425
//! (2012), section 3 and Appendix A (WRAPPED, not the six-metre virtual scale):
//! https://hal.inria.fr/hal-00688679
//! Author-uploaded readable copy:
//! https://www.researchgate.net/publication/265569903_Physical_parameters_for_piano_modeling
//!
//! These are the authors' fitted MODEL parameters based on six measured strings,
//! not 88 factory measurements. Effective diameter/density homogenize windings;
//! density is not the literal density of copper. Preserve their tabulated tension
//! and precision rather than quietly forcing equal temperament. A4 is near 441 Hz.
//! A0/Bb0/B0 and C8 are absent in Appendix A: the four extensions below are
//! explicit estimates. A0's 2.01 m length is from Steinway's Model D specifications.
//! Unison counts, duplex lengths, patch areas/thicknesses and relaxation spectrum
//! remain estimated. Do not relabel these as measured felt coupons.
use super::geometry::Course;
use fs_material::{WoolFelt, visco::GeneralizedMaxwell};
use fs_math::det;
use std::f64::consts::{PI, TAU};

const E: f64 = 2.02e11;
// MIDI 24..107: length [m], effective d [mm], density [kg/m3], T [N], strike x [m].
const WRAPPED: [[f64; 5]; 84] = [
    [2.007, 1.480, 57787.0, 1722.0, 0.241], // 24
    [1.997, 1.505, 52229.0, 1788.0, 0.240], // 25
    [1.981, 1.506, 47494.0, 1798.0, 0.238], // 26
    [1.965, 1.492, 43195.0, 1773.0, 0.236], // 27
    [1.938, 1.460, 39737.0, 1703.0, 0.233], // 28
    [1.911, 1.419, 36571.0, 1618.0, 0.229], // 29
    [1.879, 1.370, 33851.0, 1515.0, 0.225], // 30
    [1.842, 1.316, 31521.0, 1404.0, 0.221], // 31
    [1.805, 1.262, 29376.0, 1297.0, 0.217], // 32
    [1.762, 1.207, 27588.0, 1190.0, 0.211], // 33
    [1.709, 1.148, 26243.0, 1083.0, 0.205], // 34
    [1.655, 1.096, 25043.0, 991.0, 0.199], // 35
    [1.602, 1.051, 23919.0, 915.0, 0.192], // 36
    [1.548, 1.012, 22925.0, 853.0, 0.186], // 37
    [1.495, 0.982, 21997.0, 807.0, 0.179], // 38
    [1.442, 0.960, 21160.0, 774.0, 0.173], // 39
    [1.378, 0.937, 20738.0, 741.0, 0.165], // 40
    [1.837, 1.117, 7887.0, 799.0, 0.220], // 41
    [1.757, 1.108, 7734.0, 791.0, 0.211], // 42
    [1.660, 1.091, 7772.0, 773.0, 0.199], // 43
    [1.591, 1.095, 7589.0, 783.0, 0.191], // 44
    [1.482, 1.071, 7850.0, 754.0, 0.178], // 45
    [1.403, 1.067, 7850.0, 754.0, 0.168], // 46
    [1.329, 1.065, 7850.0, 756.0, 0.160], // 47
    [1.259, 1.063, 7850.0, 759.0, 0.151], // 48
    [1.192, 1.061, 7850.0, 762.0, 0.143], // 49
    [1.129, 1.059, 7850.0, 764.0, 0.136], // 50
    [1.070, 1.057, 7850.0, 766.0, 0.128], // 51
    [1.013, 1.053, 7850.0, 767.0, 0.122], // 52
    [0.960, 1.049, 7850.0, 766.0, 0.115], // 53
    [0.909, 1.045, 7850.0, 765.0, 0.109], // 54
    [0.861, 1.039, 7850.0, 763.0, 0.103], // 55
    [0.816, 1.033, 7850.0, 759.0, 0.098], // 56
    [0.773, 1.027, 7850.0, 755.0, 0.093], // 57
    [0.732, 1.020, 7850.0, 751.0, 0.088], // 58
    [0.694, 1.013, 7850.0, 746.0, 0.083], // 59
    [0.657, 1.006, 7850.0, 741.0, 0.079], // 60
    [0.622, 0.999, 7850.0, 735.0, 0.075], // 61
    [0.590, 0.991, 7850.0, 730.0, 0.071], // 62
    [0.559, 0.984, 7850.0, 725.0, 0.067], // 63
    [0.529, 0.977, 7850.0, 720.0, 0.064], // 64
    [0.501, 0.970, 7850.0, 715.0, 0.060], // 65
    [0.475, 0.964, 7850.0, 711.0, 0.057], // 66
    [0.450, 0.958, 7850.0, 707.0, 0.054], // 67
    [0.426, 0.952, 7850.0, 704.0, 0.051], // 68
    [0.404, 0.947, 7850.0, 701.0, 0.048], // 69
    [0.383, 0.941, 7850.0, 699.0, 0.046], // 70
    [0.363, 0.937, 7850.0, 697.0, 0.044], // 71
    [0.344, 0.932, 7850.0, 696.0, 0.041], // 72
    [0.326, 0.928, 7850.0, 695.0, 0.039], // 73
    [0.308, 0.924, 7850.0, 694.0, 0.037], // 74
    [0.292, 0.921, 7850.0, 694.0, 0.035], // 75
    [0.277, 0.917, 7850.0, 694.0, 0.033], // 76
    [0.262, 0.914, 7850.0, 694.0, 0.031], // 77
    [0.249, 0.910, 7850.0, 695.0, 0.030], // 78
    [0.236, 0.907, 7850.0, 695.0, 0.028], // 79
    [0.223, 0.904, 7850.0, 696.0, 0.027], // 80
    [0.211, 0.901, 7850.0, 697.0, 0.025], // 81
    [0.200, 0.898, 7850.0, 697.0, 0.024], // 82
    [0.190, 0.894, 7850.0, 697.0, 0.023], // 83
    [0.180, 0.891, 7850.0, 697.0, 0.022], // 84
    [0.171, 0.887, 7850.0, 697.0, 0.020], // 85
    [0.162, 0.883, 7850.0, 697.0, 0.019], // 86
    [0.153, 0.879, 7850.0, 696.0, 0.018], // 87
    [0.145, 0.875, 7850.0, 695.0, 0.017], // 88
    [0.138, 0.870, 7850.0, 693.0, 0.017], // 89
    [0.130, 0.866, 7850.0, 691.0, 0.016], // 90
    [0.124, 0.860, 7850.0, 689.0, 0.015], // 91
    [0.117, 0.855, 7850.0, 686.0, 0.014], // 92
    [0.111, 0.850, 7850.0, 683.0, 0.013], // 93
    [0.105, 0.844, 7850.0, 679.0, 0.013], // 94
    [0.100, 0.837, 7850.0, 674.0, 0.012], // 95
    [0.095, 0.831, 7850.0, 670.0, 0.011], // 96
    [0.090, 0.824, 7850.0, 664.0, 0.011], // 97
    [0.085, 0.817, 7850.0, 659.0, 0.010], // 98
    [0.081, 0.810, 7850.0, 652.0, 0.010], // 99
    [0.076, 0.802, 7850.0, 646.0, 0.009], // 100
    [0.072, 0.795, 7850.0, 639.0, 0.009], // 101
    [0.069, 0.787, 7850.0, 631.0, 0.008], // 102
    [0.065, 0.778, 7850.0, 623.0, 0.008], // 103
    [0.062, 0.770, 7850.0, 615.0, 0.007], // 104
    [0.058, 0.761, 7850.0, 606.0, 0.007], // 105
    [0.055, 0.752, 7850.0, 598.0, 0.007], // 106
    [0.052, 0.743, 7850.0, 588.0, 0.006], // 107
 ];

/// Section 3.1, equations (2)-(4), with piano key number i = MIDI - 20.
/// K is TOTAL hammer force coefficient [N/m^p], not stress and not per string.
pub fn hammer_parameters(midi: u8) -> Result<(f64, f64, f64), String> {
    if !(21..=108).contains(&midi) { return Err("hammer key outside A0..C8".into()); }
    let i = f64::from(midi - 20);
    Ok((0.0112 - 6.2348e-5 * i,
        2.4295e-4 * i * i - 0.007703 * i + 2.337,
        det::pow(10.0, 5.3097e-2 * i + 7.6425)))
}

/// Cold conversion of the published force envelope F=K*delta^p to the EXISTING
/// WoolFelt stress law. Dividing by TOTAL patch area ensures three unison patches
/// sum to K*delta^p rather than tripling the hammer's stiffness.
/// q, residual crush and Prony times are estimates; the spectrum is scaled to
/// each hammer's actual reference tangent, never the old uniform 5 MPa modulus.
pub fn hammer_material(c: &Course) -> Result<(WoolFelt, GeneralizedMaxwell), String> {
    c.validate()?;
    let (_, p, k) = hammer_parameters(c.midi)?;
    let reference = 0.2;
    let stress = k * det::pow(reference * c.felt_thickness_m, p) / c.felt_area_m2;
    let law = WoolFelt::new(stress, reference, p, p + 0.7, 0.25, 0.8)
        .map_err(|e| e.to_string())?;
    let tangent = p * stress / reference;
    let prony = GeneralizedMaxwell::new(0.5 * tangent,
        vec![(0.4 * tangent, 0.0002), (0.1 * tangent, 0.004)])
        .map_err(|e| e.to_string())?;
    Ok((law, prony))
}

pub fn courses() -> Result<Vec<Course>, String> {
    let mut result = Vec::with_capacity(88);
    for midi in 21..=108 {
        let (length, d_mm, rho, tension, strike) = if (24..=107).contains(&midi) {
            let [l, d, rho, t, x] = WRAPPED[usize::from(midi - 24)];
            (l, d, rho, t, x)
        } else {
            // Four explicitly ESTIMATED end extensions. Bass: preserve C1's
            // effective core and tension, solve wrapping mass at the desired f1.
            // Treble: continue diameter/length, solve tension with EI held fixed.
            let f = 441.0 * det::pow(2.0, (f64::from(midi) - 69.0) / 12.0);
            let (l, d): (f64, f64) = if midi < 24 {
                (2.010 - 0.001 * f64::from(midi - 21), 1.480e-3)
            } else { (0.052 / det::pow(2.0, 1.0 / 12.0), 0.734e-3) };
            let a = PI * d * d / 4.0;
            let ei = E * PI * d.powi(4) / 64.0;
            let wave = PI / l;
            let (rho, tension) = if midi < 24 {
                ((1722.0 * wave * wave + ei * wave.powi(4)) / (TAU * f).powi(2) / a, 1722.0)
            } else {
                (7850.0, (TAU * f).powi(2) * 7850.0 * a / (wave * wave) - ei * wave * wave)
            };
            (l, d * 1000.0, rho, tension, 0.12 * l)
        };
        let d = 0.001 * d_mm;
        let fraction = f64::from(midi - 21) / 87.0;
        let c = Course {
            midi, unison: if midi < 31 { 1 } else if midi <= 40 { 2 } else { 3 },
            length_m: length, linear_density_kg_m: rho * PI * d * d / 4.0,
            tension_n: tension, flexural_rigidity_nm2: E * PI * d.powi(4) / 64.0,
            strike_fraction: strike / length, hammer_mass_kg: hammer_parameters(midi)?.0,
            felt_thickness_m: 0.009 - 0.003 * fraction, felt_area_m2: 1.0e-4,
            duplex_length_m: if midi <= 40 { 0.0 } else { length / 3.03 },
            detune_cents: 0.8,
        };
        c.validate()?;
        result.push(c);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::Uniaxial;
    #[test]
    fn published_scale_break_and_tensions_are_not_smoothed_or_retuned() {
        let s = courses().unwrap();
        assert_eq!(s.len(), 88);
        assert_eq!(s[3].length_m, 2.007);
        assert_eq!(s[3].tension_n, 1722.0);
        assert_eq!(s[19].length_m, 1.378); // E2: wrapped bass
        assert_eq!(s[20].length_m, 1.837); // F2: long plain-string bridge
        assert_eq!(s[48].length_m, 0.404);
        assert_eq!(s[48].tension_n, 701.0);
        assert_eq!(s[86].length_m, 0.052);
        assert!((s[48].partial_hz(1, 701.0) - 441.0).abs() < 2.0);
        for c in &s {
            let nominal = 441.0 * 2.0f64.powf((f64::from(c.midi) - 69.0) / 12.0);
            // The source rounds millimetric lengths, especially in the treble.
            assert!((c.partial_hz(1, c.tension_n) / nominal - 1.0).abs() < 0.025);
        }
        assert_eq!(super::super::geometry::read_scale(&super::super::geometry::write_scale(&s)).unwrap(), s);
    }
    #[test]
    fn per_key_stress_reconstructs_total_hammer_force_without_unison_multiplication() {
        for c in courses().unwrap() {
            let (law, prony) = hammer_material(&c).unwrap();
            let (mass, p, k) = hammer_parameters(c.midi).unwrap();
            assert_eq!(mass, c.hammer_mass_kg);
            for strain in [0.01, 0.08, 0.2, 0.4] {
                let per_string = c.felt_area_m2 / c.unison as f64
                    * law.stress(strain, &law.initial_state());
                let force = k * (strain * c.felt_thickness_m).powf(p);
                assert!((per_string * c.unison as f64 / force - 1.0).abs() < 1e-10);
            }
            let instant = prony.e_inf + prony.terms.iter().map(|t| t.0).sum::<f64>();
            assert!((instant / law.tangent(0.2, &law.initial_state()) - 1.0).abs() < 1e-10);
        }
        assert!(hammer_parameters(20).is_err());
        assert!(hammer_parameters(109).is_err());
    }
    #[test]
    fn published_hammer_polynomials_match_table_examples() {
        let (m, p, k) = hammer_parameters(24).unwrap();
        assert!((m * 1000.0 - 10.95).abs() < 0.005);
        assert!((p - 2.310).abs() < 0.0005);
        assert!((k / 7.160e7 - 1.0).abs() < 1e-4);
        let (m, p, k) = hammer_parameters(69).unwrap();
        assert!((m * 1000.0 - 8.14).abs() < 0.005);
        assert!((p - 2.543).abs() < 0.0005);
        assert!((k / 1.755e10 - 1.0).abs() < 1e-4);
    }
}
