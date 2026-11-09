//! Nonthermal emission on the original spectral multiple-scattering walk.
//!
//! The integral source is `j = P/(4*pi*V)`, independently of absorptivity or
//! temperature. Each native flight integrates that source over its complete
//! ray with total extinction. A sampled scattering vertex supplies only the
//! continuation term. This is the same source/continuation split as native
//! LTE transport, with line/bin-integrated radiance instead of spectral density.

use super::super::{
    VolumeTransportError, cells::validate_ray, checkpoint, mixture::PositiveSum, positive_product,
    scattering::ScatteringGrid,
};
use super::{LuminescenceGrid, segment_emission};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;
use fs_rand::StreamKey;

/// Work allowance for one complete camera/receiver history.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LuminescenceScatteringBudget {
    /// Maximum scattering vertices. A surviving continuation at this cap
    /// refuses the entire sample; it is never silently replaced by darkness.
    pub max_scattering_events: u32,
    /// Actual DDA visits summed over the entire history, including transparent
    /// emitting cells and complete source integrals beyond sampled vertices.
    /// Each segment also retains the original medium's independent limit.
    pub max_cell_visits: u32,
}

/// One completed nonthermal radiance history, with its actual work counts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteredLuminescenceSample<const M: usize> {
    /// Line/bin-integrated W/(m² sr), in the original emission-node order.
    /// Includes ballistic and multiply scattered emission, with absorption.
    pub radiance_w_m2_sr: [f64; M],
    pub scattering_events: u32,
    pub ray_segments: u32,
    pub cell_visits: u32,
}

/// Borrow accepted cell emission and the same immutable absorption snapshot,
/// extended by its original elastic scattering and reciprocal angular laws.
/// Binding allocates nothing and does not copy the emitting-cell map.
///
/// The medium is open and index one. Emission remains uniform and isotropic
/// within each supplied original cell. These are stationary or quasistatic
/// observations of accepted endpoint/interval-mean powers, not a new kinetics,
/// heating, photon-recycling, surface or time-of-flight solve. LTE radiation is
/// separate and is not evaluated by this observer.
pub struct ScatteringLuminescenceGrid<'s, 'm, 'g, const M: usize> {
    field: &'s LuminescenceGrid<'m, 'g, M>,
    scattering: &'s ScatteringGrid<'g>,
    identity: ContentHash,
}

impl<'s, 'm, 'g, const M: usize> ScatteringLuminescenceGrid<'s, 'm, 'g, M> {
    /// Bind the complete original absorption identity, including density,
    /// geometry and cell-selected absorption. A matching cell count alone is
    /// insufficient. Scattering may use independent cell-selected spectra and
    /// phase functions, but its wavelength domain must cover every source node.
    /// Additional active opacity on the emission field is refused: the native
    /// scattering flight sampler must own that same extinction before it can
    /// sample the corresponding physical path measure.
    pub fn try_new(
        cx: &Cx<'_>,
        field: &'s LuminescenceGrid<'m, 'g, M>,
        scattering: &'s ScatteringGrid<'g>,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if field.medium.content_identity() != scattering.absorption_grid().content_identity()
            || !field.additional_absorption.is_empty()
        {
            return Err(VolumeTransportError::SourceIdentity);
        }
        let (minimum, maximum) = scattering.wavelength_bounds_nm();
        if field
            .wavelengths_nm()
            .iter()
            .any(|value| !(minimum..=maximum).contains(value))
        {
            return Err(VolumeTransportError::Wavelength);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.scattered-luminescence.v1");
        hash.update(field.content_identity().as_bytes());
        hash.update(scattering.content_identity().as_bytes());
        Ok(Self {
            field,
            scattering,
            identity: hash.finalize(),
        })
    }

    /// Numerical emission, absorption, scattering and phase snapshot identity.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.field.source_identity()
    }

    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; M] {
        self.field.wavelengths_nm()
    }

    /// Sample one complete path from an observer into the medium. The seed and
    /// logical path ID select an independent deterministic Philox stream; no
    /// thread identity or mutable shared random state enters this observation.
    /// Average independent completed histories to estimate radiance. There is
    /// no sampling-error or experimental-validation claim from one sample.
    ///
    /// The zero-scattering limit retains the exact original cell integral and
    /// consumes no variates. A collision cap, crossing cap, cancellation or
    /// numerical refusal returns no partial radiance. To retry a refused image,
    /// increase its work allowance and retain all original sample identities;
    /// do not discard individual refused histories or replace them with zero.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        budget: LuminescenceScatteringBudget,
        seed: u64,
        path_id: u64,
    ) -> Result<ScatteredLuminescenceSample<M>, VolumeTransportError> {
        let mut hash = DomainHasher::new("org.frankensim.render.luminescence-camera-stream.v1");
        hash.update(&seed.to_le_bytes());
        hash.update(&path_id.to_le_bytes());
        let seed = u64::from_le_bytes(
            hash.finalize().as_bytes()[..8]
                .try_into()
                .expect("eight hash bytes"),
        );
        let mut stream = StreamKey {
            seed,
            kernel: 0x6c75_6d73,
            tile: 0,
        }
        .stream();
        self.sample_ray_with_samples(cx, origin, direction, budget, |_| {
            std::array::from_fn(|_| stream.next_f64())
        })
    }

    /// Expose the original four flight/phase variates for deterministic
    /// quadrature and independent checks. Each visited segment requests
    /// `[spectral_lane, distance_survival, phase_cosine, phase_azimuth]` in
    /// `[0,1)`. The callback's index is the number of previous scattering
    /// vertices. General random observations require independent uniform draws.
    /// No cutoff, Russian roulette, absorptivity multiplier or wavelength width
    /// is applied to the source radiance.
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        mut origin: Point3,
        mut direction: Vec3,
        budget: LuminescenceScatteringBudget,
        mut samples: impl FnMut(u32) -> [f64; 4],
    ) -> Result<ScatteredLuminescenceSample<M>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, direction, f64::INFINITY)?;
        if !self.scattering.has_scattering() {
            let sample = self.field.transport_counted(
                cx,
                origin,
                direction,
                f64::INFINITY,
                u64::from(budget.max_cell_visits),
            )?;
            return Ok(ScatteredLuminescenceSample {
                radiance_w_m2_sr: sample.transport.emitted_radiance_w_m2_sr,
                scattering_events: 0,
                ray_segments: 1,
                cell_visits: sample.cell_visits,
            });
        }
        let mut radiance: [PositiveSum; M] = std::array::from_fn(|_| PositiveSum::default());
        // Keep weights in logarithmic form so a bright emitter can remain
        // visible after a product of small scattering/survival probabilities.
        let mut log_weight = [0.0; M];
        let mut scattering_events = 0_u32;
        let mut ray_segments = 0_u32;
        let mut cell_visits = 0_u32;
        loop {
            checkpoint(cx)?;
            let draw = samples(scattering_events);
            if draw
                .iter()
                .any(|value| !value.is_finite() || !(0.0..1.0).contains(value))
            {
                return Err(VolumeTransportError::Sample);
            }
            let (flight, visits) = self.scattering.sample_segment_observed_source(
                cx,
                origin,
                direction,
                self.field.wavelengths_nm(),
                [draw[0], draw[1]],
                budget.max_cell_visits - cell_visits,
                |cell, length, rho, extinction, depth| {
                    if length == 0.0 {
                        return Ok(());
                    }
                    if let Ok(index) = self.field.emissivity.binary_search_by_key(&cell, |e| e.0) {
                        let source = &self.field.emissivity[index].1;
                        for lane in 0..M {
                            if log_weight[lane] == f64::NEG_INFINITY {
                                continue;
                            }
                            let tau = if rho == 0.0 || extinction[lane] == 0.0 {
                                0.0
                            } else {
                                positive_product(rho, extinction[lane], length)
                            };
                            radiance[lane].add(segment_emission(
                                source[lane],
                                rho,
                                extinction[lane],
                                length,
                                tau,
                                depth[lane] - log_weight[lane],
                            ))?;
                        }
                    }
                    Ok(())
                },
            )?;
            cell_visits = cell_visits
                .checked_add(visits)
                .ok_or(VolumeTransportError::Budget)?;
            ray_segments = ray_segments
                .checked_add(1)
                .ok_or(VolumeTransportError::Budget)?;
            let Some(event) = flight.event else {
                break;
            };
            multiply_weights(&mut log_weight, &flight.weight)?;
            if log_weight.iter().all(|&value| value == f64::NEG_INFINITY) {
                break;
            }
            if scattering_events == budget.max_scattering_events {
                return Err(VolumeTransportError::Budget);
            }
            let phase = self.scattering.phase_at_cell(event.cell_index)?;
            let angular =
                phase.sample_packet(self.field.wavelengths_nm(), direction, draw[2], draw[3])?;
            multiply_weights(&mut log_weight, &angular.weight)?;
            scattering_events += 1;
            if log_weight.iter().all(|&value| value == f64::NEG_INFINITY) {
                break;
            }
            origin = event.position;
            direction = angular.direction;
        }
        let mut result = ScatteredLuminescenceSample {
            radiance_w_m2_sr: [0.0; M],
            scattering_events,
            ray_segments,
            cell_visits,
        };
        for (lane, sum) in radiance.into_iter().enumerate() {
            result.radiance_w_m2_sr[lane] = sum.finish()?;
        }
        checkpoint(cx)?;
        Ok(result)
    }
}

fn multiply_weights<const M: usize>(
    logarithms: &mut [f64; M],
    factors: &[f64; M],
) -> Result<(), VolumeTransportError> {
    for (value, &factor) in logarithms.iter_mut().zip(factors) {
        if !factor.is_finite() || factor < 0.0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        if factor == 0.0 {
            *value = f64::NEG_INFINITY;
        } else if factor != 1.0 {
            *value += det::ln(factor);
            if value.is_nan() || *value == f64::INFINITY {
                return Err(VolumeTransportError::Unrepresentable);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
