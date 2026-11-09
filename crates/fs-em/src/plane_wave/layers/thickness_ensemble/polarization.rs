//! Joint TE/TM statistics of the same authored thickness realizations.
//!
//! An ensemble of diagonal Jones responses acts as
//! J_out[a,b] = E[r_a conj(r_b)] J_in[a,b], not as the Jones matrix
//! diag(E[r_TE], E[r_TM]). The complex cross moment retains polarization
//! correlation and depolarization. See "The Complete Positivity of Classical
//! Polarization Maps", <https://arxiv.org/abs/1303.6990>.
//!
//! Axes are [TE, TM] with the parent's SAME tangential-E sign on incident
//! and reflected rays. Channel amplitudes are power normalized:
//! a_TE = sqrt(Y_TE/2) E_TE and a_TM = sqrt(Y_TM/2) E_TM, where the incident
//! lossless normal admittances are positive. The trace of J=E[a a*] is
//! incident normal power in the caller's units. Reflection uses the same
//! admittances, so each native tangential-E reflection is also its channel's
//! power-wave reflection. This frame is NOT a silently rotated optical s/p
//! frame. A consumer using right-handed propagation frames must apply its
//! explicit TM sign/basis change.
//!
//! The incident polarization statistics must be independent of the authored
//! thickness realization. The coherent/fluctuating split describes the
//! response operator at the common entrance plane, not inferred source
//! coherence. Transmission statistics require an explicit common downstream
//! plane; each realization propagates through its own remaining exit-medium
//! distance before any coherent average or covariance is formed.

use super::{Moments, PreparedLayerThicknessEnsemble, ThicknessEnsembleReflection, bad, poll};
use crate::plane_wave::{Polarization, ReflectionError};
use fs_math::{c64::C64, dd::Dd, det};

mod transmission;
use transmission::TransmissionMoments;
pub use transmission::{
    PolarizedThicknessEnsembleScattering, ThicknessEnsembleTransmission,
    ThicknessEnsembleTransmittedPower,
};

/// Joint polarization response, retaining the original scalar power statistics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedThicknessEnsembleReflection {
    channels: [ThicknessEnsembleReflection; 2],
    cross: C64,
    coherent_cross: C64,
    fluctuating_cross: C64,
}

/// Reflection of one admitted power-normalized incident coherency matrix.
/// All powers have the same units as the input coherency trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThicknessEnsemblePolarizedPower {
    /// Total reflected coherency, including the ensemble fluctuations.
    pub reflected_coherency: [[C64; 2]; 2],
    /// Action of the retained mean Jones response on the input coherency.
    pub coherent_reflected_coherency: [[C64; 2]; 2],
    /// Positive covariance contribution; this remains reflected light.
    pub fluctuating_reflected_coherency: [[C64; 2]; 2],
    /// Normal power entering each realization's own exit half-space.
    pub transmitted_power: f64,
    /// Original positive finite-film material loss, never polarization variance.
    pub absorbed_power: f64,
}

impl<'r, 's, 'm> PreparedLayerThicknessEnsemble<'r, 's, 'm> {
    /// Evaluate TE and TM on each SAME illuminated joint geometry row.
    ///
    /// The explicit work limit counts both polarizations:
    /// 2 * rows * max(original layer count, 1), including dark rows. It is
    /// checked before any solve. Materials are reused from the original cache,
    /// no allocation occurs, and cancellation is checked before either solve
    /// and before publication. The central covariance uses compensated online
    /// updates, never a subtraction of two almost equal cross moments.
    pub fn evaluate_polarized(
        &self,
        incidence_cosine: f64,
        maximum_layer_realizations: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedThicknessEnsembleReflection, ReflectionError> {
        self.evaluate_polarized_impl(
            incidence_cosine,
            None,
            maximum_layer_realizations,
            cancelled,
        )
        .map(|(reflected, _)| reflected)
    }

    fn evaluate_polarized_impl(
        &self,
        incidence_cosine: f64,
        reference_plane_m: Option<f64>,
        maximum_layer_realizations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(PolarizedThicknessEnsembleReflection, TransmissionMoments), ReflectionError> {
        poll(&mut cancelled)?;
        let work = self
            .layer_realization_work
            .checked_mul(2)
            .ok_or_else(|| bad("polarized thickness ensemble work overflow"))?;
        if work > maximum_layer_realizations {
            return Err(bad("polarized thickness ensemble exceeds its work budget"));
        }
        let transmitted_plane = reference_plane_m
            .map(|plane| transmission::admit_plane(self, incidence_cosine, plane, &mut cancelled))
            .transpose()?;
        let mut transmitted = TransmissionMoments::new();
        let mut te = Moments::new();
        let mut tm = Moments::new();
        let mut covariance = [Dd::ZERO; 2];
        for row in self.realizations {
            poll(&mut cancelled)?;
            if row.probability == 0.0 {
                continue;
            }
            let s = self.source.response_with_layer_thicknesses(
                incidence_cosine,
                Polarization::Te,
                Some(row.thicknesses_m),
            )?;
            poll(&mut cancelled)?;
            let p = self.source.response_with_layer_thicknesses(
                incidence_cosine,
                Polarization::Tm,
                Some(row.thicknesses_m),
            )?;
            // Both marginal means have the identical previous weight.
            let weight = Dd::from_f64(row.probability);
            let central_weight = weight * (te.weight / (te.weight + weight));
            let ds = [
                Dd::from_f64(s.reflected_tangential_e.re) - te.mean[0],
                Dd::from_f64(s.reflected_tangential_e.im) - te.mean[1],
            ];
            let dp = [
                Dd::from_f64(p.reflected_tangential_e.re) - tm.mean[0],
                Dd::from_f64(p.reflected_tangential_e.im) - tm.mean[1],
            ];
            covariance[0] = covariance[0] + central_weight * (ds[0] * dp[0] + ds[1] * dp[1]);
            covariance[1] = covariance[1] + central_weight * (ds[1] * dp[0] - ds[0] * dp[1]);
            if let Some(plane) = &transmitted_plane {
                transmitted.push(plane, row, [s, p])?;
            }
            te.push(row.probability, s);
            tm.push(row.probability, p);
        }
        poll(&mut cancelled)?;
        let mean_cross = [
            te.mean[0] * tm.mean[0] + te.mean[1] * tm.mean[1],
            te.mean[1] * tm.mean[0] - te.mean[0] * tm.mean[1],
        ];
        let covariance = covariance.map(|value| value / te.weight);
        let complex = |value: [Dd; 2]| C64::new(value[0].to_f64(), value[1].to_f64());
        let result = PolarizedThicknessEnsembleReflection {
            channels: [
                te.finish(self.source.layer_count())?,
                tm.finish(self.source.layer_count())?,
            ],
            cross: complex(core::array::from_fn(|i| mean_cross[i] + covariance[i])),
            coherent_cross: complex(mean_cross),
            fluctuating_cross: complex(covariance),
        };
        for moment in [
            result.reflected_response_moment(),
            result.coherent_response_moment(),
            result.fluctuating_response_moment(),
        ] {
            validate_coherency(moment)?;
        }
        Ok((result, transmitted))
    }
}

impl PolarizedThicknessEnsembleReflection {
    /// Scalar TE/TM statistics, in that order, from the same geometry rows.
    pub const fn channels(self) -> [ThicknessEnsembleReflection; 2] {
        self.channels
    }

    /// E[r_a conj(r_b)]. This is a response moment, not an illumination state.
    pub fn reflected_response_moment(self) -> [[C64; 2]; 2] {
        moment_matrix(
            self.channels.map(|x| x.mean_reflected_power_fraction),
            self.cross,
        )
    }

    /// E[r_a] conj(E[r_b]), the dyad of the mean Jones response.
    pub fn coherent_response_moment(self) -> [[C64; 2]; 2] {
        moment_matrix(
            self.channels.map(|x| x.coherent_reflected_power_fraction),
            self.coherent_cross,
        )
    }

    /// E[(r_a-E[r_a]) conj(r_b-E[r_b])], accumulated as a central moment.
    pub fn fluctuating_response_moment(self) -> [[C64; 2]; 2] {
        moment_matrix(
            self.channels
                .map(|x| x.fluctuating_reflected_power_fraction),
            self.fluctuating_cross,
        )
    }

    /// Total mean reflected Mueller matrix for S=[I,Q,U,V], where
    /// I=J00+J11, Q=J00-J11, U=2 Re J01 and V=2 Im J01 in the declared axes.
    ///
    /// A rounded I/Q block can erase a weak channel near Brewster incidence.
    /// Use reflect_coherency for physical propagation; this matrix is an
    /// interoperability/inspection representation, not a hidden Jones fit.
    pub fn reflected_mueller(self) -> [[f64; 4]; 4] {
        let [s, p] = self.channels.map(|x| x.mean_reflected_power_fraction);
        let mean = 0.5 * s + 0.5 * p;
        let difference = 0.5 * s - 0.5 * p;
        [
            [mean, difference, 0.0, 0.0],
            [difference, mean, 0.0, 0.0],
            [0.0, 0.0, self.cross.re, -self.cross.im],
            [0.0, 0.0, self.cross.im, self.cross.re],
        ]
    }

    /// Propagate any finite Hermitian positive incident polarization state.
    ///
    /// Diagonal powers are applied directly, retaining weak channels that
    /// would disappear when reconstructing them from rounded I/Q matrix
    /// entries. This constant-work query repeats no layer or material solve.
    /// Inputs are statistical power waves, not physical tangential E².
    pub fn reflect_coherency(
        self,
        incident: [[C64; 2]; 2],
    ) -> Result<ThicknessEnsemblePolarizedPower, ReflectionError> {
        validate_coherency(incident)?;
        let reflected = |moment: [[C64; 2]; 2]| {
            let diagonal = [
                incident[0][0].re * moment[0][0].re,
                incident[1][1].re * moment[1][1].re,
            ];
            moment_matrix(diagonal, incident[0][1] * moment[0][1])
        };
        let power = |coefficients: [f64; 2]| {
            (Dd::from_f64(incident[0][0].re) * Dd::from_f64(coefficients[0])
                + Dd::from_f64(incident[1][1].re) * Dd::from_f64(coefficients[1]))
            .to_f64()
        };
        let result = ThicknessEnsemblePolarizedPower {
            reflected_coherency: reflected(self.reflected_response_moment()),
            coherent_reflected_coherency: reflected(self.coherent_response_moment()),
            fluctuating_reflected_coherency: reflected(self.fluctuating_response_moment()),
            transmitted_power: power(self.channels.map(|x| x.mean_transmitted_power_fraction)),
            absorbed_power: power(self.channels.map(|x| x.mean_absorbed_power_fraction)),
        };
        for matrix in [
            result.reflected_coherency,
            result.coherent_reflected_coherency,
            result.fluctuating_reflected_coherency,
        ] {
            validate_coherency(matrix)?;
        }
        if !result.transmitted_power.is_finite() || !result.absorbed_power.is_finite() {
            return Err(bad("polarized thickness ensemble power is unrepresentable"));
        }
        Ok(result)
    }
}

fn moment_matrix(diagonal: [f64; 2], cross: C64) -> [[C64; 2]; 2] {
    [
        [C64::from_re(diagonal[0]), cross],
        [cross.conj(), C64::from_re(diagonal[1])],
    ]
}

fn validate_coherency(matrix: [[C64; 2]; 2]) -> Result<(), ReflectionError> {
    if matrix
        .iter()
        .flatten()
        .any(|z| !z.re.is_finite() || !z.im.is_finite())
        || matrix[0][0].im != 0.0
        || matrix[1][1].im != 0.0
        || matrix[0][0].re < 0.0
        || matrix[1][1].re < 0.0
        || matrix[1][0] != matrix[0][1].conj()
    {
        return Err(bad(
            "polarization coherency must be finite Hermitian and positive",
        ));
    }
    let large = matrix[0][0].re.max(matrix[1][1].re);
    let small = matrix[0][0].re.min(matrix[1][1].re);
    let cross = matrix[0][1].abs();
    // Divide by sqrt(large), not by large: extreme but representable cross
    // correlations survive even when small/large or large*small underflows.
    let valid = if small == 0.0 {
        cross == 0.0
    } else {
        cross / det::sqrt(large) <= det::sqrt(small) * (1.0 + 64.0 * f64::EPSILON)
    };
    if !valid {
        return Err(bad(
            "polarization coherency violates its positive covariance bound",
        ));
    }
    Ok(())
}
