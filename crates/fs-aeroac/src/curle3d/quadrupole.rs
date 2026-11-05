//! Volume-stress radiation from resolved three-dimensional flow.
//!
//! Lighthill's analogy gives `p = integral T_ij(y) d_i d_j G(x-y) dV`,
//! with `G = exp(ikr)/(4 pi r)` for peak phasors exp(-i omega t). Each input
//! is an actual symmetric stress moment `Q_ij = integral_cell T_ij dV` [N m].
//! Distributed cells retain their own propagation phases; only each cell,
//! not the whole source region, is represented as compact. No source-size
//! error or spatial-convergence claim is inferred from a point quadrature.
//!
//! The complete tensor, including its trace, is propagated. The resolved
//! low-Mach term is `rho u_i u_j`; omitted pressure/density, viscous and SGS
//! stresses are not reconstructed from an eddy viscosity or a noise table.
//! Transform the actual stress history: multiplying velocity phasors and
//! keeping the same frequency is wrong because the source is quadratic.
//!
//! Reference: Lighthill (1952), *On sound generated aerodynamically I*,
//! <https://doi.org/10.1098/rspa.1952.0060>; resolved-flow application:
//! Sarkar & Hussaini, NASA-CR-191543 (1993),
//! <https://ntrs.nasa.gov/citations/19940017407>.

use super::statistics::{ComplexSum, StatisticsSum, finite, invalid, nonfinite};
use super::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, det};

mod kernel;
/// Physical-clock acquisition and propagation of resolved volume histories.
pub mod recorded;
mod welch;

/// Physical scope retained by every quadrupole observation.
pub const QUADRUPOLE_SCOPE: &str = "Prescribed symmetric 3D Lighthill stress moments [N m]; \
homogeneous stationary lossless free-space point quadrature with coherent near and far fields; \
no source-size validation, unresolved SGS sound, boundaries, moving-source correction, \
mean-flow propagation, absorption, source feedback or measured-SPL qualification.";

/// Six independent entries of one physical symmetric stress moment.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct QuadrupolePhasor {
    /// Qxx, Qyy, Qzz, Qxy, Qxz, Qyz [N m], as complex peak amplitudes.
    /// Shear entries are physical Qxy etc.; do not multiply them by two.
    /// Values share one emission clock. Propagation phase is applied once.
    pub stress_moment_n_m: [C64; 6],
}

/// Extract the instantaneous resolved Reynolds-stress moment of one cell.
///
/// Returns xx, yy, zz, xy, xz, yz of `rho u_i u_j volume` in N m. Density is
/// the cell's actual value; velocity is expressed in the propagation frame.
/// Integrate or transform these instantaneous moments using the original
/// sampling clock. No mean subtraction, detrending, SGS completion, force
/// term or acoustic calibration is inferred. Full Lighthill stresses may
/// instead be supplied directly through [`QuadrupolePhasor`].
pub fn resolved_reynolds_moment_n_m(
    density_kg_m3: f64,
    cell_volume_m3: f64,
    velocity_m_s: [f64; 3],
) -> Result<[f64; 6], AeroacError> {
    kernel::reynolds_moment(density_kg_m3, cell_volume_m3, velocity_m_s)
        .ok_or_else(|| invalid("resolved stress requires positive finite density/volume and representable rho u_i u_j dV"))
}

#[derive(Clone, Debug)]
struct Transfer {
    source_m: [f64; 3],
    delay_s: f64,
    field: [[C64; 6]; 4],
}

/// Prepared pressure and full vector velocity from physical volume stresses.
///
/// Preparation stores 24 complex coefficients per source. Evaluation is O(N),
/// allocation-free, and contains no transcendental evaluation. Fixed-order
/// double-double reductions retain cancellation between strong sources.
/// Existing monopole/dipole observers keep their smaller original transfers.
#[derive(Clone, Debug)]
pub struct PreparedQuadrupoleObserver {
    acoustic: FluidAcousticState,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer>,
}

impl PreparedQuadrupoleObserver {
    /// Prepare stationary source/receiver paths in one original acoustic state.
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
    ) -> Result<Self, AeroacError> {
        Self::new_cancellable(
            acoustic,
            frequency_hz,
            observer_m,
            source_positions_m,
            max_sources,
            || false,
        )
    }

    /// Bounded preparation, polling before every source and before publication.
    /// Wrong counts, coincident paths, invalid scales or allocation failure
    /// refuse the entire observer. Empty source arrays represent silence.
    pub fn new_cancellable(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if source_positions_m.len() > max_sources {
            return Err(invalid(
                "quadrupole source count exceeds preparation budget",
            ));
        }
        let omega = core::f64::consts::TAU * frequency_hz;
        if [
            frequency_hz,
            omega,
            omega / acoustic.sound_speed_m_s(),
            acoustic.density_kg_m3() * omega,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "quadrupole frequency and wave scales must be positive finite",
            ));
        }
        if observer_m.iter().any(|x| !x.is_finite()) {
            return Err(nonfinite("quadrupole observer position"));
        }
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("quadrupole transfer allocation exceeds available memory"))?;
        for source in source_positions_m {
            poll(&mut cancelled)?;
            if source.iter().any(|x| !x.is_finite()) {
                return Err(nonfinite("quadrupole source position"));
            }
            let displacement = core::array::from_fn(|i| observer_m[i] - source[i]);
            let field = kernel::transfer(
                acoustic.density_kg_m3(),
                acoustic.sound_speed_m_s(),
                frequency_hz,
                displacement,
            )
            .ok_or_else(|| {
                invalid("quadrupole separation or field transfer is not representable")
            })?;
            let range = det::hypot(
                det::hypot(displacement[0], displacement[1]),
                displacement[2],
            );
            let delay_s = range / acoustic.sound_speed_m_s();
            if !delay_s.is_finite() || delay_s <= 0.0 {
                return Err(invalid("quadrupole propagation delay is not representable"));
            }
            transfers.push(Transfer {
                source_m: *source,
                delay_s,
                field,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            frequency_hz,
            observer_m,
            transfers,
        })
    }

    /// Original material/EOS acoustic snapshot used by every path.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    /// Common source frequency [Hz].
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    /// Receiver position in the original Cartesian frame [m].
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    /// Number of source moments required by every evaluation.
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }
    /// Original physical source position [m].
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|t| t.source_m)
    }
    /// Travel time r/c [s], already included in the harmonic transfer.
    pub fn source_delay_s(&self, index: usize) -> Option<f64> {
        self.transfers.get(index).map(|t| t.delay_s)
    }

    /// Sum coherent stress sources before forming signed energy transport.
    pub fn evaluate(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    /// Allocation-free field evaluation with a poll at most every 256 sources.
    pub fn evaluate_cancellable(
        &self,
        sources: &[QuadrupolePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(sources, &mut cancelled)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("quadrupole field intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: QUADRUPOLE_SCOPE,
        })
    }

    /// Feed the original finite-window observer, retaining phase and near fields.
    /// Equal-frequency force/volume/quadrupole lines are coherently coalesced by
    /// [`crate::observation::PreparedHarmonicField`] before computing intensity.
    /// Each physical source must have one owner: do not add a volume mechanism
    /// already represented by a permeable acoustic surface or a noise model.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<crate::observation::HarmonicFieldLine<3>, AeroacError> {
        let [pressure_pa, ux, uy, uz] = self.amplitudes(sources, &mut || false)?;
        Ok(crate::observation::HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }

    /// Explicit low-rank source coherence factors, with no inferred independence.
    /// Factors contain peak stress moments. A source PSD in a bin of width df
    /// must be factored as 2 df S; outputs are line moments, not PSD per hertz.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[QuadrupolePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        let pairs = modes.len().checked_mul(self.source_count());
        if modes.len() > budget.max_modes
            || pairs.is_none_or(|n| n > budget.max_source_mode_pairs)
            || modes.iter().any(|mode| mode.len() != self.source_count())
        {
            return Err(invalid(
                "quadrupole coherence shape or work budget is invalid",
            ));
        }
        let mut statistics = StatisticsSum::default();
        for mode in modes {
            statistics.add(self.amplitudes(mode, &mut || false)?)?;
        }
        statistics.finish(QUADRUPOLE_SCOPE)
    }

    fn amplitudes(
        &self,
        sources: &[QuadrupolePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[C64; 4], AeroacError> {
        poll(cancelled)?;
        if sources.len() != self.source_count() {
            return Err(invalid(
                "one symmetric moment is required per prepared quadrupole source",
            ));
        }
        let mut sums = [ComplexSum::default(); 4];
        for (index, (transfer, source)) in self.transfers.iter().zip(sources).enumerate() {
            if index.is_multiple_of(256) {
                poll(cancelled)?;
            }
            for (entry, value) in source.stress_moment_n_m.into_iter().enumerate() {
                if !finite(value) {
                    return Err(nonfinite("quadrupole stress amplitude"));
                }
                if value != C64::ZERO {
                    for (row, sum) in transfer.field.iter().zip(&mut sums) {
                        sum.add(row[entry] * value)?;
                    }
                }
            }
        }
        poll(cancelled)?;
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("accumulated quadrupole field"));
        }
        Ok(values)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "quadrupole observation",
        })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
