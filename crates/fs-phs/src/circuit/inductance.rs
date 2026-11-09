//! Reversible nonlinear inductors on the native Kirchhoff flux coordinates.
//!
//! The original branch voltage remains dλ/dt, while current is dH/dλ from one
//! bounded, strictly monotone I(λ) law. Saturation therefore changes differential
//! inductance and stored energy together. H is the integral of current over flux
//! linkage; λ I / 2 is generally wrong. No hysteresis, remanence, extrapolation,
//! or nonlinear mutual coupling is inferred from a scalar curve.
//!
//! Every branch owns its exact finite-step work and analytic discrete-current
//! derivative, using the same scalar integral kernel as nonlinear capacitors.
//! This keeps unrelated winding currents and capacitor voltages independent.

use super::{
    CircuitDae, CircuitError, CircuitGraph, MutualInductance,
    capacitance::{self, CapacitanceBlock, CapacitorCurve},
    scalar::{CurveError, ScalarCurve},
};

/// Source-bounded common-core saturation with independent winding leakage.
pub mod coupled;
pub use coupled::{MagneticStorageBlock, SharedCoreInductance};

/// Maximum retained knots per inductor, bounding each crossed-segment sum.
pub const MAX_INDUCTOR_KNOTS: usize = super::scalar::MAX_SCALAR_KNOTS;

fn curve_error(error: CurveError) -> CircuitError {
    CircuitError::Invalid {
        what: match error {
            CurveError::Knots => "bounded finite inductor flux/current knots required",
            CurveError::ZeroKnot => "inductor curve requires an explicit unenergized (0,0) knot",
            CurveError::Polarities => "inductor curve must support both flux polarities",
            CurveError::PositiveSlope => {
                "inductor knots require positive representable differential inductance"
            }
            CurveError::Arithmetic => "inductor constitutive arithmetic is not representable",
            CurveError::Energy => "inductor energy is not representable",
            CurveError::StoredEnergy => "inductor stored energy is not representable",
            CurveError::CoordinateSupport => "inductor flux is outside the supplied curve support",
            CurveError::EffortSupport => "inductor current is outside the supplied curve support",
            CurveError::Effort => "inductor current is not representable",
            CurveError::Coordinate => "inductor flux is not representable",
        },
    }
}

/// Immutable piecewise-affine current I(λ) [A] versus flux linkage λ [Wb].
///
/// Knots increase strictly in both coordinates. An interior (0,0) knot fixes
/// the zero-energy datum; both flux polarities are supported, and need not be
/// symmetric. The complete source support is closed and never extrapolated.
/// Immutable knots, slopes and energy integrals are shared across clones.
///
/// Source authority and geometry belong to the caller. A secant inductance
/// λ/I is not the incremental inductance dλ/dI and cannot replace the integral.
#[derive(Clone, Debug, PartialEq)]
pub struct InductorFluxCurve {
    pub(super) scalar: ScalarCurve,
}

impl InductorFluxCurve {
    /// Admit 3..=MAX_INDUCTOR_KNOTS finite ordered (flux linkage [Wb], current
    /// [A]) knots, including (0,0) strictly inside the support. Every segment
    /// requires positive finite differential inductance and representable energy.
    pub fn try_new(knots: Vec<(f64, f64)>) -> Result<Self, CircuitError> {
        Ok(Self {
            scalar: ScalarCurve::try_new(knots).map_err(curve_error)?,
        })
    }

    /// Original (flux linkage [Wb], current [A]) knots, without resampling.
    pub fn knots(&self) -> &[(f64, f64)] {
        self.scalar.knots()
    }

    /// Closed flux-linkage support [Wb]; Newton probes also need support.
    pub fn flux_range_wb(&self) -> (f64, f64) {
        self.scalar.range()
    }

    /// Actual branch current dH/dλ [A] at the stored flux linkage.
    pub fn current_a(&self, flux_wb: f64) -> Result<f64, CircuitError> {
        self.scalar.effort(flux_wb).map_err(curve_error)
    }

    /// Initialize flux linkage [Wb] from current using the SAME monotone law.
    pub fn flux_wb(&self, current_a: f64) -> Result<f64, CircuitError> {
        self.scalar.coordinate(current_a).map_err(curve_error)
    }

    /// dλ/dI [H]. An interior knot selects its right segment; the final knot
    /// selects its left segment. A unique derivative at a kink is not claimed.
    pub fn differential_inductance_h(&self, flux_wb: f64) -> Result<f64, CircuitError> {
        self.scalar.inverse_slope(flux_wb).map_err(curve_error)
    }

    /// Integral from zero to λ of I(s) ds [J], with the unenergized datum fixed.
    pub fn stored_energy_j(&self, flux_wb: f64) -> Result<f64, CircuitError> {
        self.scalar.energy(flux_wb).map_err(curve_error)
    }

    /// Stable energy difference [J] across every retained source segment.
    pub fn energy_change_j(&self, before_wb: f64, after_wb: f64) -> Result<f64, CircuitError> {
        self.scalar
            .energy_change(before_wb, after_wb)
            .map_err(curve_error)
    }

    /// Work-conjugate mean current [A] over this branch's OWN flux change.
    /// Satisfies ΔH = I_bar Δλ across knots without subtracting biased energies
    /// or applying a global correction to other branches. Coincidence returns
    /// the original current, including the same admitted endpoint support.
    pub fn discrete_current_a(&self, before_wb: f64, after_wb: f64) -> Result<f64, CircuitError> {
        self.scalar
            .discrete_effort(before_wb, after_wb)
            .map_err(curve_error)
    }

    /// Derivative of [`Self::discrete_current_a`] with respect to end flux
    /// linkage. Differentiates the same integral across every crossed knot.
    /// Coincident knots use the same one-sided segment convention as
    /// [`Self::differential_inductance_h`].
    pub fn discrete_current_end_derivative(
        &self,
        before_wb: f64,
        after_wb: f64,
    ) -> Result<f64, CircuitError> {
        self.scalar
            .discrete_effort_end_derivative(before_wb, after_wb)
            .map_err(curve_error)
    }

    /// Integral of I(λ)-I(reference_wb) from before_wb to after_wb [J].
    /// Retains small work beside large stored bias without subtracting the
    /// endpoint energies. The reference must lie in the same source support.
    pub fn energy_change_remainder_j(
        &self,
        before_wb: f64,
        after_wb: f64,
        reference_wb: f64,
    ) -> Result<f64, CircuitError> {
        self.scalar
            .midpoint_remainder(before_wb, after_wb, reference_wb)
            .map_err(curve_error)
    }
}

/// Replacement storage for one ORIGINAL graph branch, not its compressed slot.
#[derive(Clone, Debug, PartialEq)]
pub struct InductorCurve {
    /// Existing inductor index in [`CircuitGraph::branches`].
    pub branch: usize,
    /// Complete constitutive storage for this branch, including its self term.
    pub curve: InductorFluxCurve,
}

/// Bind nonlinear inductors beside any disjoint linear mutual-inductance block.
/// Empty bindings preserve the original exact-quadratic native circuit.
///
/// Bound graph inductances remain positive admitted placeholders. Flux states,
/// incidence, KCL/KVL, ports, resistors and ideal transformers are unchanged.
/// A bound branch must not occur in any mutual link: a scalar law cannot replace
/// a coupled magnetic potential. Unsupported links and out-of-support states
/// refuse; no additional leakage, clamp, or extrapolation is manufactured.
pub fn assemble_circuit_with_inductor_curves(
    graph: &CircuitGraph,
    curves: &[InductorCurve],
    mutual: &[MutualInductance],
) -> Result<CircuitDae, CircuitError> {
    assemble_circuit_with_storage_curves(graph, &[], &[], curves, mutual)
}

/// Compose scalar I(λ), scalar V(q), reciprocal C blocks, and disjoint linear
/// mutual L in ONE native descriptor with the original charge and flux states.
///
/// Nonlinear inductors must be disjoint from every supplied mutual link.
/// Capacitance blocks must be disjoint from scalar capacitor curves and each
/// other. Each nonlinear coordinate uses its own exact constitutive work;
/// all remaining quadratic cross terms retain their symmetric mean efforts.
/// Branch-domain and Newton refusals publish no partial state.
pub fn assemble_circuit_with_storage_curves(
    graph: &CircuitGraph,
    capacitance_blocks: &[CapacitanceBlock],
    capacitor_curves: &[CapacitorCurve],
    inductor_curves: &[InductorCurve],
    mutual: &[MutualInductance],
) -> Result<CircuitDae, CircuitError> {
    assemble_circuit_with_storage_blocks(
        graph,
        capacitance_blocks,
        capacitor_curves,
        inductor_curves,
        &[],
        mutual,
    )
}

/// Compose common-core nonlinear magnetic blocks, scalar flux/charge laws,
/// reciprocal capacitance blocks and disjoint linear mutual inductance on the
/// original Kirchhoff descriptor. Each block owns its complete magnetic energy;
/// overlapping branches or external mutual links refuse rather than erasing a
/// physical cross term. Positive graph coefficients are topology placeholders
/// on bound branches, never extra magnetic storage.
pub fn assemble_circuit_with_storage_blocks(
    graph: &CircuitGraph,
    capacitance_blocks: &[CapacitanceBlock],
    capacitor_curves: &[CapacitorCurve],
    inductor_curves: &[InductorCurve],
    magnetic_blocks: &[MagneticStorageBlock],
    mutual: &[MutualInductance],
) -> Result<CircuitDae, CircuitError> {
    capacitance::assemble_storage_blocks(
        graph,
        capacitance_blocks,
        capacitor_curves,
        inductor_curves,
        magnetic_blocks,
        mutual,
    )
}
