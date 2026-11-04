//! Seed-free contour discovery and the original coupled physical observations.

use super::{C64, CancelGate, DgLoadedModalPole, DgLoadedPoleError, PreparedModalDgPoles};
use fs_la::nonlinear_eig::beyn::{BeynError, BeynReport};

pub use fs_la::nonlinear_eig::beyn::{BeynOptions, Ellipse as ResonanceRegion};

/// Finite modes resolved in a declared Laplace-plane ellipse [1/s].
/// Argument-count, rank and nested-quadrature agreement are numerical evidence,
/// not interval enclosures or a rigorous certificate of spectral completeness.
#[derive(Debug, Clone, PartialEq)]
pub struct DgLoadedModalSpectrum {
    /// Original mixed-unit state coordinates and independently checked rows.
    /// The order is the numerical owner's deterministic real/imaginary ordering.
    pub poles: Vec<DgLoadedModalPole>,
    /// Independently sorted pre-refinement estimates, not mode identifiers.
    pub unrefined_laplace_poles_per_s: Vec<C64>,
    pub argument_count: usize,
    pub resolved_rank: usize,
    pub quadrature_nodes: usize,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub moment_relative_change: f64,
    /// Coordinate-dependent sampled inverse margin, not an interval bound.
    pub minimum_reciprocal_condition: f64,
}

#[derive(Debug)]
pub enum DgLoadedSpectrumError {
    Contour(BeynError),
    PhysicalPole(DgLoadedPoleError),
}
impl core::fmt::Display for DgLoadedSpectrumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete modal/DG spectrum refused: {self:?}")
    }
}
impl std::error::Error for DgLoadedSpectrumError {}

impl PreparedModalDgPoles<'_> {
    /// Discover coupled resonances and relaxations without individual mode seeds.
    ///
    /// The same entire A+s*B+s^2*M pencil retains every original Maxwell, duct,
    /// branch, terminal and finite-channel state and algebraic reaction. It has
    /// no eliminated material poles or isolated-component impedance inverses.
    /// The native Beyn owner reuses nested contour nodes, checks an independent
    /// argument-principle count against moment rank, and refines its discovered
    /// vectors with the native Newton owner. Each result then passes the same
    /// original-row residual and physical observer as `resonance_pole`.
    ///
    /// Two moments resolve at most the available independent probe directions.
    /// A quadratic pencil can have more finite roots than its state dimension;
    /// choose smaller regions when the count/rank exceeds this method's
    /// capacity. Insufficient probes, unresolved clusters, contour singularities,
    /// or an inadmissible physical pole refuse the whole result. No mode is
    /// silently filtered. Empty accepted regions carry the same count and
    /// quadrature evidence. Neither finite nor infinite spectrum completeness,
    /// discretization convergence, or transfer observability is certified.
    ///
    /// Quadrature nodes evaluate the mathematical continuation of the declared
    /// finite model. They need not lie in an empirical real-frequency band.
    /// Every returned non-real pole, including a numerically small imaginary
    /// part, must still pass each original duct/terminal band at abs(Im(s)).
    /// Exactly real relaxation has no oscillation frequency; general poles are
    /// never snapped to the real axis. Band admission does not establish the
    /// empirical accuracy of complex-frequency continuation.
    ///
    /// In addition to the numerical owner's bounds, the query reserves
    /// 128*n*l+256*(n+l) bytes and
    /// (256*l+16*(nodes+l*refinement_iterations))*n^2 work for callbacks and
    /// physical observations, where l is the probe count. The already admitted
    /// immutable preparation cache and source storage are excluded. Source
    /// history and clocks are unchanged, including on cancellation or refusal.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalSpectrum, DgLoadedSpectrumError> {
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
        } = self
            .pencil
            .discover(region, options, || gate.is_requested())
            .map_err(DgLoadedSpectrumError::Contour)?;
        let mut physical = Vec::with_capacity(poles.len());
        for pole in poles {
            physical.push(
                self.observe_pole(pole, options.refinement.residual_tolerance, gate)
                    .map_err(DgLoadedSpectrumError::PhysicalPole)?,
            );
        }
        super::poll(Some(gate))
            .map_err(|error| DgLoadedSpectrumError::PhysicalPole(error.into()))?;
        Ok(DgLoadedModalSpectrum {
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
