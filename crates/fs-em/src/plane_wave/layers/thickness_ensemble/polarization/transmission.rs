//! Common-plane transmitted polarization and downstream material loss.
//!
//! The parent's exp(+i omega t) convention gives exp(-i k_z dz). Each row
//! shares z=0 and its incident phasor; only after propagating from its physical
//! exit to the authored common plane may its transmitted Jones response enter
//! a coherent mean. Normal-power waves retain the tangential-E phase and
//! sqrt(Re(Y_exit)/Y_inc) normalization, including passive absorbing exits.
//! Evanescent and exactly critical exits have zero normal-power channels;
//! this is not their nonzero near-field electric coherency.

use super::super::super::{LayeredReflection, admittance_state_from_normal};
use super::super::{LayerThicknessRealization, PreparedLayerThicknessEnsemble, bad, poll};
use super::{PolarizedThicknessEnsembleReflection, moment_matrix, validate_coherency};
use crate::plane_wave::{Polarization, ReflectionError};
use fs_math::{c64::C64, dd::Dd, det};

/// Both outgoing polarization responses from the same two native solves per row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedThicknessEnsembleScattering {
    /// Reflection is referenced to the original common entrance plane.
    pub reflected: PolarizedThicknessEnsembleReflection,
    /// Transmission is referenced to the explicitly authored downstream plane.
    pub transmitted: ThicknessEnsembleTransmission,
}

/// Power-normalized transmitted response moments at one physical plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThicknessEnsembleTransmission {
    reference_plane_m: f64,
    mean_power_wave: [C64; 2],
    total: [[C64; 2]; 2],
    coherent: [[C64; 2]; 2],
    fluctuating: [[C64; 2]; 2],
    exit_absorption: [f64; 2],
}

/// Propagated polarization and heat in the exit medium before the common plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThicknessEnsembleTransmittedPower {
    pub transmitted_coherency: [[C64; 2]; 2],
    pub coherent_transmitted_coherency: [[C64; 2]; 2],
    pub fluctuating_transmitted_coherency: [[C64; 2]; 2],
    /// Physical loss between each row's exit face and the common plane.
    /// Finite-film loss remains owned by the paired reflected result.
    pub exit_medium_absorbed_power: f64,
}

impl<'r, 's, 'm> PreparedLayerThicknessEnsemble<'r, 's, 'm> {
    /// Evaluate joint reflection and transmission at an authored common plane.
    ///
    /// The plane is a finite z coordinate [m] from the common entrance and
    /// must lie at or beyond EVERY row's summed film thickness, including
    /// zero-probability rows. No backwards extrapolation through an absorbing
    /// medium is allowed. Entrance illumination statistics must be independent
    /// of the row. This is an independent-column/time ensemble, not diffraction.
    ///
    /// The native layer-work bound is 2 * rows * max(layers, 1), as in
    /// evaluate_polarized, plus bounded geometry and moment contraction.
    /// No allocation or original material query occurs. Cancellation publishes
    /// neither partial channel. Flux at this plane plus downstream absorption
    /// equals the ordinary exit-face transmission; covariance is not heat.
    pub fn evaluate_polarized_at_plane(
        &self,
        incidence_cosine: f64,
        reference_plane_m: f64,
        maximum_layer_realizations: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedThicknessEnsembleScattering, ReflectionError> {
        let (reflected, moments) = self.evaluate_polarized_impl(
            incidence_cosine,
            Some(reference_plane_m),
            maximum_layer_realizations,
            cancelled,
        )?;
        let transmitted =
            moments.finish(reference_plane_m, reflected, self.source.layer_count())?;
        Ok(PolarizedThicknessEnsembleScattering {
            reflected,
            transmitted,
        })
    }
}

impl ThicknessEnsembleTransmission {
    pub const fn reference_plane_m(self) -> f64 {
        self.reference_plane_m
    }
    /// E[t_TE], E[t_TM] in normal-power coordinates at the common plane.
    pub const fn mean_transmitted_power_wave(self) -> [C64; 2] {
        self.mean_power_wave
    }
    /// E[t_a conj(t_b)], including all retained ensemble fluctuations.
    pub const fn transmitted_response_moment(self) -> [[C64; 2]; 2] {
        self.total
    }
    pub const fn coherent_response_moment(self) -> [[C64; 2]; 2] {
        self.coherent
    }
    pub const fn fluctuating_response_moment(self) -> [[C64; 2]; 2] {
        self.fluctuating
    }
    /// Mean exit-medium loss fraction for unit incident TE and TM power.
    pub const fn exit_medium_absorbed_power_fractions(self) -> [f64; 2] {
        self.exit_absorption
    }

    /// Mean transmitted Mueller matrix on the parent's [I,Q,U,V] axes,
    /// V=2 Im J01. Direct coherency propagation retains weak diagonal channels.
    pub fn transmitted_mueller(self) -> [[f64; 4]; 4] {
        let s = self.total[0][0].re;
        let p = self.total[1][1].re;
        let cross = self.total[0][1];
        let mean = 0.5 * s + 0.5 * p;
        let difference = 0.5 * s - 0.5 * p;
        [
            [mean, difference, 0.0, 0.0],
            [difference, mean, 0.0, 0.0],
            [0.0, 0.0, cross.re, -cross.im],
            [0.0, 0.0, cross.im, cross.re],
        ]
    }

    /// Apply any finite positive Hermitian incident power coherency.
    /// Constant work; no repeated geometry, material or layer solve.
    pub fn transmit_coherency(
        self,
        incident: [[C64; 2]; 2],
    ) -> Result<ThicknessEnsembleTransmittedPower, ReflectionError> {
        validate_coherency(incident)?;
        let apply = |moment: [[C64; 2]; 2]| {
            moment_matrix(
                [
                    incident[0][0].re * moment[0][0].re,
                    incident[1][1].re * moment[1][1].re,
                ],
                incident[0][1] * moment[0][1],
            )
        };
        let result = ThicknessEnsembleTransmittedPower {
            transmitted_coherency: apply(self.total),
            coherent_transmitted_coherency: apply(self.coherent),
            fluctuating_transmitted_coherency: apply(self.fluctuating),
            exit_medium_absorbed_power: (Dd::from_f64(incident[0][0].re)
                * Dd::from_f64(self.exit_absorption[0])
                + Dd::from_f64(incident[1][1].re) * Dd::from_f64(self.exit_absorption[1]))
            .to_f64(),
        };
        for matrix in [
            result.transmitted_coherency,
            result.coherent_transmitted_coherency,
            result.fluctuating_transmitted_coherency,
        ] {
            validate_coherency(matrix)?;
        }
        if !result.exit_medium_absorbed_power.is_finite() {
            return Err(bad("transmitted ensemble power is unrepresentable"));
        }
        Ok(result)
    }
}

pub(super) struct TransmissionPlane {
    z: f64,
    phase_wavenumber: Dd,
    decay_wavenumber: f64,
    exit_states: [[C64; 2]; 2],
}

pub(super) fn admit_plane(
    ensemble: &PreparedLayerThicknessEnsemble<'_, '_, '_>,
    cosine: f64,
    z: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransmissionPlane, ReflectionError> {
    if !z.is_finite() || z < 0.0 || !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
        return Err(bad(
            "transmitted ensemble requires a finite downstream plane and valid incidence",
        ));
    }
    for row in ensemble.realizations {
        poll(cancelled)?;
        gap_to_plane(z, row.thicknesses_m)?;
    }
    let source = ensemble.source;
    let q = source.exit_wave.normal(cosine)?;
    Ok(TransmissionPlane {
        z,
        phase_wavenumber: Dd::from_f64(q.re) * Dd::from_f64(source.wave_number),
        decay_wavenumber: q.im * source.wave_number,
        exit_states: [
            admittance_state_from_normal(source.exit_wave, cosine, Polarization::Te, q)?,
            admittance_state_from_normal(source.exit_wave, cosine, Polarization::Tm, q)?,
        ],
    })
}

fn gap_to_plane(z: f64, thicknesses: &[f64]) -> Result<Dd, ReflectionError> {
    let extent = thicknesses
        .iter()
        .fold(Dd::ZERO, |sum, &d| sum + Dd::from_f64(d));
    let gap = Dd::from_f64(z) - extent;
    if !gap.to_f64().is_finite() || gap.lt(Dd::ZERO) {
        return Err(bad(
            "common transmission plane lies before an authored exit face",
        ));
    }
    Ok(gap)
}

pub(super) struct TransmissionMoments {
    weight: Dd,
    // Complex TE/TM means, native plane powers, and positive exit-medium loss.
    mean: [Dd; 8],
    covariance: [Dd; 4],
}

impl TransmissionMoments {
    pub(super) fn new() -> Self {
        Self {
            weight: Dd::ZERO,
            mean: [Dd::ZERO; 8],
            covariance: [Dd::ZERO; 4],
        }
    }

    pub(super) fn push(
        &mut self,
        plane: &TransmissionPlane,
        row: &LayerThicknessRealization<'_>,
        responses: [LayeredReflection; 2],
    ) -> Result<(), ReflectionError> {
        let gap = gap_to_plane(plane.z, row.thicknesses_m)?;
        let decay = plane.decay_wavenumber * gap.to_f64();
        let attenuation = det::exp(decay);
        let power_attenuation = det::exp(2.0 * decay);
        let loss_fraction = -det::expm1(2.0 * decay);
        let amplitudes = if attenuation == 0.0
            || responses
                .iter()
                .all(|r| r.transmitted_power_fraction == 0.0)
        {
            [C64::ZERO; 2]
        } else {
            let phase = plane.phase_wavenumber * gap;
            if !phase.to_f64().is_finite() {
                return Err(bad("common-plane transmission phase is unrepresentable"));
            }
            let (s, c) = det::sin_cos_dd(phase);
            let propagation = C64::new(c.to_f64(), -s.to_f64()).scale(attenuation);
            let mut amplitudes = [C64::ZERO; 2];
            for i in 0..2 {
                let response = responses[i];
                if response.transmitted_power_fraction == 0.0 {
                    continue;
                }
                // Retain the native positive normal power instead of forming a
                // possibly overflowing admittance ratio. The better E/H chart
                // supplies tangential-E phase; H is rotated by the SAME exit
                // admittance pair, including its absorbing-medium phase.
                let e = response.transmitted_tangential_e;
                let h = response.transmitted_conjugate_h;
                let phase = if maximum(e) >= maximum(h) {
                    unit_phase(e)?
                } else {
                    unit_phase(h)?
                        * unit_phase(plane.exit_states[i][0])?
                        * unit_phase(plane.exit_states[i][1])?.conj()
                };
                amplitudes[i] =
                    phase.scale(det::sqrt(response.transmitted_power_fraction)) * propagation;
            }
            amplitudes
        };
        let values = [
            amplitudes[0].re,
            amplitudes[0].im,
            amplitudes[1].re,
            amplitudes[1].im,
            responses[0].transmitted_power_fraction * power_attenuation,
            responses[1].transmitted_power_fraction * power_attenuation,
            responses[0].transmitted_power_fraction * loss_fraction,
            responses[1].transmitted_power_fraction * loss_fraction,
        ];
        if values.iter().any(|x| !x.is_finite()) {
            return Err(bad("common-plane transmitted response is unrepresentable"));
        }
        let weight = Dd::from_f64(row.probability);
        let total = self.weight + weight;
        let delta: [Dd; 8] = core::array::from_fn(|i| Dd::from_f64(values[i]) - self.mean[i]);
        let central_weight = weight * (self.weight / total);
        let [sr, si, pr, pi, ..] = delta;
        let terms = [
            sr * sr + si * si,
            pr * pr + pi * pi,
            sr * pr + si * pi,
            si * pr - sr * pi,
        ];
        for (sum, term) in self.covariance.iter_mut().zip(terms) {
            *sum = *sum + central_weight * term;
        }
        let ratio = weight / total;
        for (mean, difference) in self.mean.iter_mut().zip(delta) {
            *mean = *mean + ratio * difference;
        }
        self.weight = total;
        Ok(())
    }

    fn finish(
        self,
        reference_plane_m: f64,
        reflected: PolarizedThicknessEnsembleReflection,
        layers: usize,
    ) -> Result<ThicknessEnsembleTransmission, ReflectionError> {
        let [sr, si, pr, pi, ts, tp, as_, ap] = self.mean;
        let covariance = self.covariance.map(|x| x / self.weight);
        let cross = [sr * pr + si * pi, si * pr - sr * pi];
        let coherent_diagonal = [(sr * sr + si * si).to_f64(), (pr * pr + pi * pi).to_f64()];
        let result = ThicknessEnsembleTransmission {
            reference_plane_m,
            mean_power_wave: [
                C64::new(sr.to_f64(), si.to_f64()),
                C64::new(pr.to_f64(), pi.to_f64()),
            ],
            total: moment_matrix(
                [ts.to_f64(), tp.to_f64()],
                C64::new(
                    (cross[0] + covariance[2]).to_f64(),
                    (cross[1] + covariance[3]).to_f64(),
                ),
            ),
            coherent: moment_matrix(
                coherent_diagonal,
                C64::new(cross[0].to_f64(), cross[1].to_f64()),
            ),
            fluctuating: moment_matrix(
                [covariance[0].to_f64(), covariance[1].to_f64()],
                C64::new(covariance[2].to_f64(), covariance[3].to_f64()),
            ),
            exit_absorption: [as_.to_f64(), ap.to_f64()],
        };
        for moment in [result.total, result.coherent, result.fluctuating] {
            validate_coherency(moment)?;
        }
        let tolerance = 1024.0 * f64::EPSILON * (layers + 1) as f64;
        for (i, channel) in reflected.channels().into_iter().enumerate() {
            let power = result.total[i][i].re;
            let absorbed = result.exit_absorption[i];
            if !absorbed.is_finite()
                || absorbed < 0.0
                || (power + absorbed - channel.mean_transmitted_power_fraction).abs() > tolerance
                || (result.coherent[i][i].re + result.fluctuating[i][i].re - power).abs()
                    > tolerance
            {
                return Err(bad("common-plane transmission power balance failed"));
            }
        }
        Ok(result)
    }
}

fn maximum(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}

fn unit_phase(z: C64) -> Result<C64, ReflectionError> {
    let magnitude = maximum(z);
    if magnitude == 0.0 || !magnitude.is_finite() {
        return Err(bad(
            "nonzero transmitted power has no representable field phase",
        ));
    }
    let scaled = C64::new(z.re / magnitude, z.im / magnitude);
    Ok(scaled.scale(1.0 / det::sqrt(scaled.norm_sq())))
}
