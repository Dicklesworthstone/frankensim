//! Physical observations of the shared contour-discovered descriptor modes.

use super::{C64, HarmonicModel, HarmonicPole, HarmonicPoleError, PolePencil, poll};
use fs_la::nonlinear_eig::beyn::{BeynError, BeynReport, beyn_with_cancel};

pub use fs_la::nonlinear_eig::beyn::{BeynOptions, Ellipse as ResonanceRegion};

/// Resolved finite modes in a declared region of the Laplace s-plane [1/s].
/// Numerical quadrature/rank agreement does not constitute an interval proof
/// of spectral completeness. Every reported mode passes the original physical
/// descriptor residual and carries independent storage/loss observations.
#[derive(Debug, Clone)]
pub struct HarmonicSpectrum {
    pub poles: Vec<HarmonicPole>,
    pub unrefined_laplace_poles_per_s: Vec<C64>,
    pub argument_count: usize,
    pub resolved_rank: usize,
    pub quadrature_nodes: usize,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum HarmonicSpectrumError {
    Invalid { what: &'static str },
    Contour(BeynError),
    PhysicalPole(HarmonicPoleError),
}
impl core::fmt::Display for HarmonicSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical resonance spectrum refused: {self:?}")
    }
}
impl std::error::Error for HarmonicSpectrumError {}

impl HarmonicModel {
    /// Discover finite free modes in a Laplace-plane ellipse without supplying
    /// individual mode seeds. The original complete sE-(J-R)W pencil is linear
    /// and entire, including when E is singular. Infinite descriptor roots are
    /// outside every finite region; no E inverse or eliminated constraint is
    /// introduced. Internal material memory and zero-output modes remain.
    ///
    /// Uses exactly the same pencil, refinement and physical observations as
    /// `resonance_pole`. The contour owner checks nested quadrature, numerical
    /// rank and an independent argument-principle count. Insufficient probes,
    /// unresolved clusters, poor contour conditioning or limits refuse the
    /// whole result. These checks supply numerical evidence, not a rigorous
    /// complete-spectrum certificate. Refusal leaves AC factors unchanged.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicSpectrum, HarmonicSpectrumError> {
        poll(&mut cancelled).map_err(|e| HarmonicSpectrumError::PhysicalPole(e.into()))?;
        if options.probe_columns == 0 || options.probe_columns > self.n {
            return Err(HarmonicSpectrumError::Invalid {
                what: "physical spectrum probes must lie between one and the complete dimension",
            });
        }
        let n = self.n as u128;
        let m = self.m as u128;
        let modes = options.probe_columns as u128;
        // Reserve physical result vectors and original-row/energy observations
        // in addition to the numerical owner's conservative bound.
        let bytes = 16 * modes * (3 * n + m) + 256 * (n + modes);
        let work = modes * (8 * n * n + 4 * n * m);
        if bytes > options.max_dense_bytes as u128 || work > options.max_work as u128 {
            return Err(HarmonicSpectrumError::Invalid {
                what: "insufficient physical spectrum observation budget",
            });
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let mut pencil = PolePencil::new(self, options.refinement.spectral_scale);
        let BeynReport {
            poles,
            unrefined_values,
            argument_integral,
            coarse_argument_integral,
            argument_count,
            rank,
            nodes,
            moment_relative_change,
            minimum_reciprocal_condition,
        } = beyn_with_cancel(
            self.n,
            region,
            &numerical,
            |s, cancel| pencil.evaluate(s, cancel),
            &mut cancelled,
        )
        .map_err(HarmonicSpectrumError::Contour)?;
        let mut physical = Vec::with_capacity(poles.len());
        for pole in poles {
            physical.push(
                self.observe_pole(pole, &mut cancelled)
                    .map_err(HarmonicSpectrumError::PhysicalPole)?,
            );
        }
        poll(&mut cancelled).map_err(|e| HarmonicSpectrumError::PhysicalPole(e.into()))?;
        Ok(HarmonicSpectrum {
            poles: physical,
            unrefined_laplace_poles_per_s: unrefined_values,
            argument_count,
            resolved_rank: rank,
            quadrature_nodes: nodes,
            argument_integral,
            coarse_argument_integral,
            moment_relative_change,
            minimum_reciprocal_condition,
        })
    }
}
