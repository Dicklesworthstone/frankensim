//! Phase-preserving reflection over the native conductor and coating data.
//!
//! The characteristic-matrix assembly is shared with scalar coated Fresnel.
//! Both bare and coated amplitudes use exp(-i omega t), n = eta + i k.
//! Conjugating n leaves scalar reflectance unchanged but reverses the U/V
//! phase coupling, so the scalar bare-metal convention cannot be copied as-is.

use super::{Complex, ComplexIor, forward_sqrt, propagation, reflection};
use crate::conductor::{CoherentCoating, ConductorError, ConductorOptics};
use crate::polarization::Mueller;
use fs_math::det;

impl ConductorOptics {
    /// Reflected Stokes power at a smooth conductor interface, including its
    /// explicitly authored transparent or absorbing coherent coating, when present.
    ///
    /// The incoming and reflected frames both have +Q along s. Their p axes
    /// make (s,p,propagation) right-handed. U=2 Re(E_s conj(E_p)),
    /// V=2 Im(E_s conj(E_p)); fields use exp(-i omega t). Thus normal reflection
    /// flips U/V, while oblique metal reflection can convert linear to elliptical
    /// polarization. Rotate ray frames before applying this matrix.
    ///
    /// Wavelength is vacuum nm; the surrounding dielectric index is real and
    /// positive. Optical tables, film thickness, dispersion and source identities
    /// retain their existing owners. This is an opaque half-space model: the
    /// non-reflected power enters the substrate, not a transmitted scene ray.
    /// It does not add roughness, depolarization or coherence between paths.
    /// No allocation or phase-angle trigonometry is needed for a bare interface.
    pub fn fresnel_mueller(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        incidence_cosine: f64,
    ) -> Result<Mueller, ConductorError> {
        let substrate = self.complex_ior(wavelength_nm)?;
        if !incident_eta.is_finite()
            || incident_eta <= 0.0
            || incident_eta > super::super::MAX_IOR_COMPONENT
        {
            return Err(ConductorError::InvalidIncidentMedium);
        }
        if !incidence_cosine.is_finite() || !(0.0..=1.0).contains(&incidence_cosine) {
            return Err(ConductorError::InvalidIncidence);
        }
        // Exact grazing limit: both electric-field amplitudes tend to -1.
        // Evaluating a film's transfer matrix at q0=0 can instead produce 0/0.
        if incidence_cosine == 0.0 {
            return Ok(Mueller::IDENTITY);
        }
        let amplitudes = match &self.coating {
            Some(CoherentCoating::Dielectric(film)) => {
                let film_eta = film.ior.eval(wavelength_nm)
                    .map_err(|_| ConductorError::NumericalFailure)?;
                if film.thickness_m == 0.0 || film_eta == incident_eta {
                    // A matched film contributes only a common phase, irrelevant
                    // to Stokes power. Keep the full bare matrix, not just M00.
                    bare_amplitudes(incident_eta, substrate, incidence_cosine)
                } else {
                    coated_amplitudes(wavelength_nm, incident_eta, substrate,
                        incidence_cosine, film_eta, film.thickness_m)
                }
            }
            Some(CoherentCoating::Absorbing(film)) => {
                let ior = film.complex_ior(wavelength_nm)?;
                if film.thickness_m() == 0.0
                    || (ior.k == 0.0 && ior.eta == incident_eta)
                    || (ior.eta == substrate.eta && ior.k == substrate.k) {
                    bare_amplitudes(incident_eta, substrate, incidence_cosine)
                } else {
                    super::super::absorbing_film::reflection_amplitudes(
                        wavelength_nm, incident_eta, substrate, incidence_cosine,
                        ior, film.thickness_m())
                }
            }
            None => bare_amplitudes(incident_eta, substrate, incidence_cosine),
        };
        let s = passive_amplitude(amplitudes[0])?;
        let p = if incidence_cosine == 1.0 {
            // Isotropy fixes this relation exactly, including an interference
            // null where separately rounded s/p formulas can lose their ratio.
            Complex::new(-s.re, -s.im)
        } else {
            passive_amplitude(amplitudes[1])?
        };
        let s_power = s.norm_squared();
        let p_power = p.norm_squared();
        let mean = 0.5 * s_power + 0.5 * p_power;
        let difference = 0.5 * s_power - 0.5 * p_power;
        let real = s.re * p.re + s.im * p.im;
        let imaginary = s.im * p.re - s.re * p.im;
        Ok(Mueller([
            [mean, difference, 0.0, 0.0],
            [difference, mean, 0.0, 0.0],
            [0.0, 0.0, real, -imaginary],
            [0.0, 0.0, imaginary, real],
        ]))
    }
}

// This is the original scalar coating assembly, retaining its operation order.
// Scalar and polarized callers differ only in what they retain from the same
// two complex amplitudes; there is no second coating interference solver.
pub(super) fn coated_amplitudes(
    wavelength_nm: f64,
    incident_eta: f64,
    substrate: ComplexIor,
    cosine: f64,
    film_eta: f64,
    thickness_m: f64,
) -> [Complex; 2] {
    let n0_squared = incident_eta * incident_eta;
    let n1_squared = film_eta * film_eta;
    let transverse_squared = n0_squared * (1.0 - cosine * cosine);
    let q0 = incident_eta * cosine;
    let q1_squared = n1_squared - transverse_squared;
    let n2_squared = Complex::new(substrate.eta, substrate.k).square();
    let q2 = forward_sqrt(n2_squared - Complex::new(transverse_squared, 0.0));
    let phase_scale = core::f64::consts::TAU * (thickness_m * 1.0e9) / wavelength_nm;
    let (cos_phase, sin_phase_over_q) = propagation(q1_squared, phase_scale);
    let r_s = reflection(q0, q2, cos_phase,
        sin_phase_over_q, q1_squared * sin_phase_over_q);
    let r_p = reflection(q0 / n0_squared, q2 / n2_squared, cos_phase,
        n1_squared * sin_phase_over_q, q1_squared * sin_phase_over_q / n1_squared);
    [r_s, r_p]
}

fn bare_amplitudes(incident_eta: f64, substrate: ComplexIor, cosine: f64) -> [Complex; 2] {
    let n = Complex::new(substrate.eta, substrate.k);
    if cosine == 1.0 {
        let incident = Complex::new(incident_eta, 0.0);
        let s = (incident - n) / (incident + n);
        return [s, Complex::new(-s.re, -s.im)];
    }
    let n0_squared = incident_eta * incident_eta;
    let n2_squared = n.square();
    let q2 = forward_sqrt(n2_squared - Complex::new(n0_squared * (1.0 - cosine * cosine), 0.0));
    let q0 = incident_eta * cosine;
    [
        reflection(q0, q2, 1.0, 0.0, 0.0),
        reflection(q0 / n0_squared, q2 / n2_squared, 1.0, 0.0, 0.0),
    ]
}

fn passive_amplitude(value: Complex) -> Result<Complex, ConductorError> {
    let norm = det::hypot(value.re, value.im);
    if !norm.is_finite() || norm > 1.0 + 256.0 * f64::EPSILON {
        return Err(ConductorError::NumericalFailure);
    }
    // Project only bounded roundoff, in amplitude space, preserving the
    // relative phase and the positive-semidefinite coherency transformation.
    Ok(if norm > 1.0 { value * (1.0 / norm) } else { value })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conductor::{CoherentDielectricCoating, ConductorIorSample,
        ConductorDataStatus, ConductorSource};
    use crate::dielectric::CauchyIor;
    use crate::polarization::Stokes;

    fn optics(eta: f64, k: f64, film: Option<(f64, f64)>) -> ConductorOptics {
        let source = ConductorSource::try_new(fs_blake3::ContentHash([83; 32]),
            ConductorDataStatus::Representative).unwrap();
        let bare = ConductorOptics::try_new(core::array::from_fn(|i|
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()), source).unwrap();
        match film {
            None => bare,
            Some((n, d)) => bare.with_coating(CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(n).unwrap(), d, source).unwrap()),
        }
    }
    fn close(a: f64, b: f64) {
        assert!((a-b).abs() < 4e-12, "{a:.17e} != {b:.17e}");
    }

    #[test]
    fn g0_metal_mueller_recovers_scalar_power_and_preserves_physical_states() {
        for film in [None, Some((1.45, 125e-9)), Some((1.0, 200e-9))] {
            for (eta, k) in [(0.2, 3.0), (2.7, 4.1), (1.5, 1e-8)] {
                let law = optics(eta, k, film);
                for n0 in [1.0, 1.33, 2.4] {
                    for j in 0..=50 {
                        let c = f64::from(j) / 50.0;
                        let matrix = law.fresnel_mueller(550.0, n0, c).unwrap();
                        close(matrix.0[0][0], law.fresnel(550.0, n0, c).unwrap());
                        for (q,u,v) in [(0.0,0.0,0.0), (1.0,0.0,0.0), (-1.0,0.0,0.0),
                            (0.0,1.0,0.0), (0.0,0.0,1.0), (0.3,-0.4,0.2)] {
                            let out = matrix.apply(Stokes::try_new(1.0,q,u,v).unwrap());
                            assert!(out.i >= -2e-14 && out.i <= 1.0 + 2e-14);
                            assert!(out.q.hypot(out.u).hypot(out.v) <= out.i + 2e-14);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn g1_normal_metal_reflection_has_opposite_sp_axes() {
        for film in [None, Some((1.45,125e-9))] {
            let law = optics(0.2,3.0,film);
            let matrix = law.fresnel_mueller(550.0,1.0,1.0).unwrap();
            let r = law.fresnel(550.0,1.0,1.0).unwrap();
            let out = matrix.apply(Stokes::try_new(1.0,0.0,0.6,0.8).unwrap());
            close(out.i,r); close(out.u,-0.6*r); close(out.v,-0.8*r);
            assert_eq!(matrix.0[3][2],0.0);
        }
    }

    #[test]
    fn g1_signed_phase_matches_independent_complex_airy_answers() {
        // Independent complex Airy evaluation, not scalar-power agreement.
        // In particular the authored coating reverses V generation here.
        for (n0, eta, k, c, film, w, expected) in [
            (1.0, 0.2, 3.0, 0.5, None, 550.0,
                [0.9184110846593692, 0.04400178561613738, -0.5749762549092698, 0.7148042176556463]),
            (1.0, 0.2, 3.0, 0.5, Some((1.45, 125e-9)), 550.0,
                [0.8925003326859151, 0.009894742965482706, -0.7934452898337169, -0.4085382600770766]),
            (1.5, 0.2, 3.0, 0.5, Some((1.0, 200e-9)), 550.0,
                [0.9920935748330739, 0.006283107851867198, 0.8205895346270934, 0.5575329581659549]),
            (1.33, 2.7, 4.1, 0.73, Some((2.1, 82e-9)), 430.0,
                [0.5367161366914662, 0.08477096576298163, -0.5233624604777231, -0.08348550599551866]),
        ] {
            let m = optics(eta, k, film).fresnel_mueller(w, n0, c).unwrap();
            for (actual, expected) in [m.0[0][0], m.0[0][1], m.0[2][2], m.0[3][2]]
                .into_iter().zip(expected) {
                close(actual, expected);
            }
        }
    }

    #[test]
    fn g3_zero_and_index_matched_films_preserve_the_entire_matrix() {
        let bare = optics(0.2,3.0,None);
        for n0 in [1.0,1.33,2.4] {
            for film in [Some((1.45,0.0)),Some((n0,200e-9))] {
                for c in [0.0,1e-8,0.3,0.8,1.0] {
                    assert_eq!(optics(0.2,3.0,film).fresnel_mueller(550.0,n0,c).unwrap(),
                        bare.fresnel_mueller(550.0,n0,c).unwrap());
                }
            }
        }
    }

    #[test]
    fn g1_coated_critical_and_evanescent_matrices_retain_phase() {
        let law = optics(0.2,3.0,Some((1.0,200e-9)));
        let critical = (1.0_f64 - (1.0_f64/1.5).powi(2)).sqrt();
        let reference = law.fresnel_mueller(550.0,1.5,critical).unwrap();
        for delta in [-1e-10,1e-10] {
            let matrix = law.fresnel_mueller(550.0,1.5,critical+delta).unwrap();
            for (a,b) in matrix.0.iter().flatten().zip(reference.0.iter().flatten()) {
                assert!((a-b).abs() < 2e-8);
            }
        }
        let out = law.fresnel_mueller(550.0,1.5,0.5).unwrap()
            .apply(Stokes::try_new(1.0,0.0,1.0,0.0).unwrap());
        assert!(out.i < 0.999 && out.v.abs() > 0.01);
        let dark = optics(2.25,1e-8,Some((1.5,100e-9)))
            .fresnel_mueller(600.0,1.0,1.0).unwrap();
        assert!(dark.0[0][0] > 0.0 && dark.0[0][0] < 1e-16);
    }

    #[test]
    fn g0_mueller_refuses_nonfinite_and_out_of_domain_interfaces() {
        let law = optics(0.2,3.0,None);
        for n0 in [0.0,-1.0,101.0,f64::NAN,f64::INFINITY] {
            assert!(law.fresnel_mueller(550.0,n0,0.5).is_err());
        }
        for c in [-0.01,1.01,f64::NAN,f64::INFINITY] {
            assert!(law.fresnel_mueller(550.0,1.0,c).is_err());
        }
        for w in [379.0,781.0,f64::NAN] { assert!(law.fresnel_mueller(w,1.0,0.5).is_err()); }
        assert_eq!(law.fresnel_mueller(550.0,1.0,0.0).unwrap(),Mueller::IDENTITY);
    }

    fn absorbing(eta: f64, k: f64, thickness: f64) -> crate::conductor::CoherentAbsorbingCoating {
        use crate::conductor::{CoatingIorSample, CoherentAbsorbingCoating};
        let source = ConductorSource::try_new(fs_blake3::ContentHash([89;32]),
            ConductorDataStatus::Representative).unwrap();
        CoherentAbsorbingCoating::try_new(core::array::from_fn(|i|
            CoatingIorSample::try_new(380.0+50.0*i as f64, eta, k).unwrap()),
            thickness, source).unwrap()
    }

    #[test]
    fn g1_absorbing_film_preserves_independent_signed_phase_answers() {
        for (n0,n1,k1,eta,k,d,w,c,expected) in [
            (1.33,2.2,0.35,3.0,4.0,77e-9,430.0,0.6,
                [0.13421992524904688,0.09353863303512977,-0.08247485506403528,-0.04963275879845326]),
            (1.5,1.0,0.01,0.2,3.0,80e-9,580.0,0.55,
                [0.7448873775190354,0.22411336846870045,0.2235208615118144,-0.6742913522579203]),
            (1.0,2.0,4.0,3.0,2.0,1e-4,500.0,0.37,
                [0.6428729360715091,0.22610297479708952,-0.37348035039731886,0.4718850332318299]),
        ] {
            let law = optics(eta,k,None).with_absorbing_coating(absorbing(n1,k1,d));
            let m = law.fresnel_mueller(w,n0,c).unwrap();
            close(m.0[0][0],law.fresnel(w,n0,c).unwrap());
            for (actual,expected) in [m.0[0][0],m.0[0][1],m.0[2][2],m.0[3][2]]
                .into_iter().zip(expected) { close(actual,expected); }
        }
    }

    #[test]
    fn g3_absorbing_zero_matched_and_opaque_limits_preserve_full_mueller_response() {
        let bare = optics(0.2,3.0,None);
        for coating in [absorbing(1.7,0.4,0.0),absorbing(1.33,0.0,1e-7),
            absorbing(0.2,3.0,1e-7)] {
            for c in [0.0,0.1,0.6,1.0] {
                assert_eq!(bare.with_absorbing_coating(coating).fresnel_mueller(550.0,1.33,c).unwrap(),
                    bare.fresnel_mueller(550.0,1.33,c).unwrap());
            }
        }
        let opaque = absorbing(1.8,0.6,1e-4);
        for c in [0.1,0.6,1.0] {
            let actual = bare.with_absorbing_coating(opaque).fresnel_mueller(550.0,1.33,c).unwrap();
            let expected = optics(1.8,0.6,None).fresnel_mueller(550.0,1.33,c).unwrap();
            for (a,b) in actual.0.iter().flatten().zip(expected.0.iter().flatten()) { close(*a,*b); }
        }
    }
}
