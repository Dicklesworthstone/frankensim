//! Voltage-dependent capacitance on the original Kirchhoff DAE.
//!
//! A reversible capacitor is specified by its charge/voltage relation, not by
//! replacing a constant C with C(V) in q=C*V. Differential capacitance is dq/dV;
//! its energy is integral V(q) dq, generally neither q*V/2 nor C(V)*V²/2.
//! See Analog Devices, "How to Use LTspice Simulations to Account for the Effect
//! of Voltage Dependence" (RAQ 192) for the charge-retention distinction.
//!
//! Piecewise-linear V(q) retains every supplied knot, positive differential
//! capacitance and exact finite support. No fitting, extrapolation, clipping,
//! dielectric memory, leakage, thermal dependence or measured-data authority
//! is inferred. Source selection and interpolation accuracy belong to the
//! caller. The existing bounded descriptor solve supplies the time update;
//! this module supplies constitutive energy, componentwise mean effort and
//! its analytic end-state tangent to that same solver.

use super::{
    Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance,
    assemble_circuit_with_mutual_inductance,
};
use super::{
    inductance::InductorCurve,
    scalar::{CurveError, ScalarCurve},
};
use crate::{ConstitutiveDiscreteGradient, QuadraticStorage, Storage, linear_exact::QuadraticPart};

/// Prepared reciprocal multi-terminal storage on the original charge states.
pub mod block;
pub use block::CapacitanceBlock;

/// Maximum retained knots per capacitor, bounding each crossed-segment sum.
pub const MAX_CAPACITOR_KNOTS: usize = super::scalar::MAX_SCALAR_KNOTS;

fn invalid(what: &'static str) -> CircuitError {
    CircuitError::Invalid { what }
}

fn curve_error(error: CurveError) -> CircuitError {
    invalid(match error {
        CurveError::Knots => "bounded finite capacitor charge/voltage knots required",
        CurveError::ZeroKnot => "capacitor curve requires an explicit uncharged (0,0) knot",
        CurveError::Polarities => "capacitor curve must support both charge polarities",
        CurveError::PositiveSlope => {
            "capacitor knots require positive representable differential capacitance"
        }
        CurveError::Arithmetic => "capacitor constitutive arithmetic is not representable",
        CurveError::Energy => "capacitor energy is not representable",
        CurveError::StoredEnergy => "capacitor stored energy is not representable",
        CurveError::CoordinateSupport => "capacitor charge is outside the supplied curve support",
        CurveError::EffortSupport => "capacitor voltage is outside the supplied curve support",
        CurveError::Effort => "capacitor voltage is not representable",
        CurveError::Coordinate => "capacitor charge is not representable",
    })
}

fn finite(value: f64) -> Result<f64, CircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(curve_error(CurveError::Arithmetic))
    }
}

/// Immutable V(q), with strictly increasing charge [C] and voltage [V].
///
/// An explicit (0,0) knot sets the uncharged minimum of this reversible energy.
/// Both polarities need support so an initially uncharged element admits the
/// original centered Newton probes. Asymmetry between polarities is allowed.
/// Prepared knots, slopes and energies are immutable and shared across clones,
/// so material/thermal circuit trials do not copy the constitutive tables.
#[derive(Clone, Debug, PartialEq)]
pub struct CapacitorChargeCurve {
    pub(super) scalar: ScalarCurve,
}

impl CapacitorChargeCurve {
    /// Admit 3..=MAX_CAPACITOR_KNOTS finite ordered knots, including (0,0)
    /// strictly inside the support. Every segment must have a positive,
    /// representable differential capacitance and stored energy.
    pub fn try_new(knots: Vec<(f64, f64)>) -> Result<Self, CircuitError> {
        Ok(Self {
            scalar: ScalarCurve::try_new(knots).map_err(curve_error)?,
        })
    }

    /// Original (charge [C], voltage [V]) knots, without resampling.
    pub fn knots(&self) -> &[(f64, f64)] {
        self.scalar.knots()
    }

    /// Closed source support [C]; stepping may need room for Newton probes.
    pub fn charge_range_c(&self) -> (f64, f64) {
        self.scalar.range()
    }

    /// Actual terminal voltage dH/dq at the stored charge.
    pub fn voltage_v(&self, charge_c: f64) -> Result<f64, CircuitError> {
        self.scalar.effort(charge_c).map_err(curve_error)
    }

    /// Initialize charge from voltage using the SAME monotone relation.
    pub fn charge_c(&self, voltage_v: f64) -> Result<f64, CircuitError> {
        self.scalar.coordinate(voltage_v).map_err(curve_error)
    }

    /// dq/dV [F]. At an interior knot the right segment is returned; the last
    /// knot uses its left segment. No derivative continuity is claimed there.
    pub fn differential_capacitance_f(&self, charge_c: f64) -> Result<f64, CircuitError> {
        self.scalar.inverse_slope(charge_c).map_err(curve_error)
    }

    /// Integral from zero to q of V(s) ds [J], with the uncharged datum fixed.
    pub fn stored_energy_j(&self, charge_c: f64) -> Result<f64, CircuitError> {
        self.scalar.energy(charge_c).map_err(curve_error)
    }

    /// Stable energy difference across any retained segments [J].
    pub fn energy_change_j(&self, before_c: f64, after_c: f64) -> Result<f64, CircuitError> {
        self.scalar
            .energy_change(before_c, after_c)
            .map_err(curve_error)
    }

    /// Work-conjugate mean voltage over this capacitor's OWN charge change.
    /// Uses the original integral, not a global energy correction involving
    /// other charges or fluxes. At coincidence it is the original voltage.
    /// Same-segment arithmetic retains a rounded midpoint's missing offset;
    /// no squared charge increment can underflow and erase a knot correction.
    pub fn discrete_voltage_v(&self, before_c: f64, after_c: f64) -> Result<f64, CircuitError> {
        self.scalar
            .discrete_effort(before_c, after_c)
            .map_err(curve_error)
    }

    /// Derivative of [`Self::discrete_voltage_v`] with respect to its end
    /// charge. This differentiates the actual piecewise integral across every
    /// crossed knot. At coincident knots the right segment supplies a one-sided
    /// Newton representative (the last knot uses its left segment); a unique
    /// sensitivity at a kink is not claimed.
    pub fn discrete_voltage_end_derivative(
        &self,
        before_c: f64,
        after_c: f64,
    ) -> Result<f64, CircuitError> {
        self.scalar
            .discrete_effort_end_derivative(before_c, after_c)
            .map_err(curve_error)
    }

    /// Integral of `V(q)-V(reference_c)` from `before_c` to `after_c` [J].
    /// Evaluates the original piecewise law without subtracting large endpoint
    /// energies. A composite storage owner may use its represented midpoint
    /// as the reference while retaining every other physical energy term.
    pub fn energy_change_remainder_j(
        &self,
        before_c: f64,
        after_c: f64,
        reference_c: f64,
    ) -> Result<f64, CircuitError> {
        self.scalar
            .midpoint_remainder(before_c, after_c, reference_c)
            .map_err(curve_error)
    }
}

/// A replacement law for one ORIGINAL capacitor branch, not its compressed slot.
#[derive(Clone, Debug, PartialEq)]
pub struct CapacitorCurve {
    pub branch: usize,
    pub curve: CapacitorChargeCurve,
}

#[derive(Clone)]
struct CurveStorage {
    linear: QuadraticStorage,
    curves: Vec<(usize, ScalarCurve)>,
}

impl Storage for CurveStorage {
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }

    fn coordinate_scale(&self, coordinate: usize) -> Option<f64> {
        self.curves
            .iter()
            .find(|(i, _)| *i == coordinate)
            .map(|(_, curve)| {
                let (lo, hi) = curve.range();
                lo.abs().max(hi.abs())
            })
    }

    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(self.clone()))
    }

    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.curves
            .iter()
            .fold(self.linear.hamiltonian(x), |sum, (i, curve)| {
                sum + curve.energy(x[*i]).unwrap_or(f64::NAN)
            })
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.linear.gradient(x, out);
        for (i, curve) in &self.curves {
            out[*i] = curve.effort(x[*i]).unwrap_or(f64::NAN);
        }
    }

    fn quadratic_part(&self) -> Option<QuadraticPart> {
        Some(QuadraticPart {
            matrix: self.linear.matrix().to_vec(),
            exact: false,
        })
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let mut remainder = self.linear.midpoint_energy_remainder(a, b, mid)?;
        for (i, curve) in &self.curves {
            remainder += curve
                .midpoint_remainder(a[*i], b[*i], mid[*i])
                .unwrap_or(f64::NAN);
        }
        Some(remainder)
    }
}

impl ConstitutiveDiscreteGradient for CurveStorage {
    fn evaluate(&self, a: &[f64], b: &[f64], mid: &[f64], out: &mut [f64]) {
        self.linear.mean_gradient(a, b, mid, out);
        for (i, curve) in &self.curves {
            out[*i] = curve.discrete_effort(a[*i], b[*i]).unwrap_or(f64::NAN);
        }
    }

    fn end_jacobian_action(
        &self,
        a: &[f64],
        b: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> bool {
        self.linear.mean_end_jacobian_action(direction, out);
        for (i, curve) in &self.curves {
            let (lo, hi) = curve.range();
            if !a[*i].is_finite()
                || !b[*i].is_finite()
                || a[*i] < lo
                || a[*i] > hi
                || b[*i] < lo
                || b[*i] > hi
            {
                return false;
            }
            out[*i] = if direction[*i] == 0.0 {
                0.0
            } else {
                curve
                    .discrete_effort_end_derivative(a[*i], b[*i])
                    .unwrap_or(f64::NAN)
                    * direction[*i]
            };
        }
        out.iter().all(|value| value.is_finite())
    }
}

/// Bind nonlinear capacitors into the SAME descriptor, including optional
/// reciprocal magnetic coupling. Unbound capacitors retain their original C.
/// Graph coefficients remain ordinary positive placeholders for bound branches.
///
/// Charges, fluxes, graph incidence, source ports, resistors, transformers and
/// state ordering are unchanged. The native discrete gradient uses the exact
/// piecewise integral and analytic voltage, preserving its energy/work identity
/// even across knots. Native bounded Newton/domain refusals publish no state;
/// timestep convergence and source interpolation error remain separate checks.
/// An empty curve list returns the original assembled circuit unchanged.
pub fn assemble_circuit_with_capacitor_curves(
    graph: &CircuitGraph,
    curves: &[CapacitorCurve],
    mutual: &[MutualInductance],
) -> Result<CircuitDae, CircuitError> {
    assemble_circuit_with_capacitance_blocks(graph, &[], curves, mutual)
}

/// Assemble reciprocal coupled C, disjoint scalar V(q), and mutual L into
/// the SAME Kirchhoff descriptor and original conserved charge/flux states.
///
/// Each prepared block owns its complete capacitance matrix, including its
/// diagonal. The graph's capacitor coefficients remain positive placeholders
/// on bound branches. Other capacitors retain their ordinary scalar storage.
/// Block/curve overlap refuses: no individual law may erase a cross term.
/// All supplied block inverses are reused without refactorization; exact
/// quadratic storage stays available when no nonlinear curve is supplied.
///
/// References, branch orientations, and physical source authority belong to
/// the caller. Singular unreferenced Maxwell matrices are refused when the
/// block is prepared, never grounded or replaced by pair-capacitor fits.
pub fn assemble_circuit_with_capacitance_blocks(
    graph: &CircuitGraph,
    blocks: &[CapacitanceBlock],
    curves: &[CapacitorCurve],
    mutual: &[MutualInductance],
) -> Result<CircuitDae, CircuitError> {
    assemble_storage_curves(graph, blocks, curves, &[], mutual)
}

// One storage owner composes both scalar laws with the remaining quadratic
// block. A scalar branch must never erase a reciprocal cross term.
pub(super) fn assemble_storage_curves(
    graph: &CircuitGraph,
    blocks: &[CapacitanceBlock],
    capacitors: &[CapacitorCurve],
    inductors: &[InductorCurve],
    mutual: &[MutualInductance],
) -> Result<CircuitDae, CircuitError> {
    let mut ordered: Vec<_> = capacitors.iter().collect();
    ordered.sort_by_key(|entry| entry.branch);
    if ordered
        .windows(2)
        .any(|pair| pair[0].branch == pair[1].branch)
        || ordered.iter().any(|entry| {
            !graph
                .branches
                .get(entry.branch)
                .is_some_and(|b| matches!(b.2, Branch::Capacitor { .. }))
        })
    {
        return Err(invalid(
            "capacitor curves require distinct existing capacitor branch indices",
        ));
    }
    let mut ordered_inductors: Vec<_> = inductors.iter().collect();
    ordered_inductors.sort_by_key(|entry| entry.branch);
    if ordered_inductors
        .windows(2)
        .any(|pair| pair[0].branch == pair[1].branch)
        || ordered_inductors.iter().any(|entry| {
            !graph
                .branches
                .get(entry.branch)
                .is_some_and(|b| matches!(b.2, Branch::Inductor { .. }))
        })
    {
        return Err(invalid(
            "inductor curves require distinct existing inductor branch indices",
        ));
    }
    if mutual.iter().any(|link| {
        [link.branches.0, link.branches.1].iter().any(|branch| {
            ordered_inductors
                .binary_search_by_key(branch, |entry| entry.branch)
                .is_ok()
        })
    }) {
        return Err(invalid(
            "scalar inductor curves cannot share a mutual-inductance link",
        ));
    }
    let mut dae = assemble_circuit_with_mutual_inductance(graph, mutual)?;
    if ordered.is_empty() && ordered_inductors.is_empty() && blocks.is_empty() {
        return Ok(dae);
    }
    let n = dae.system.state_dim();
    let mut matrix = dae
        .system
        .storage()
        .exact_quadratic_matrix()
        .expect("native circuit storage is quadratic")
        .to_vec();
    block::stamp(graph, &dae, blocks, &ordered, &mut matrix)?;
    if ordered.is_empty() && ordered_inductors.is_empty() {
        dae.system.storage = Box::new(QuadraticStorage::new(matrix, n)?);
        return Ok(dae);
    }
    let mut bound = Vec::with_capacity(ordered.len() + ordered_inductors.len());
    for entry in ordered {
        let slot = graph.branches[..entry.branch]
            .iter()
            .filter(|b| matches!(b.2, Branch::Capacitor { .. }))
            .count();
        let coordinate = dae.charge_index[slot];
        matrix[coordinate * n + coordinate] = 0.0;
        bound.push((coordinate, entry.curve.scalar.clone()));
    }
    for entry in ordered_inductors {
        let slot = graph.branches[..entry.branch]
            .iter()
            .filter(|b| matches!(b.2, Branch::Inductor { .. }))
            .count();
        let coordinate = dae.flux_index[slot];
        matrix[coordinate * n + coordinate] = 0.0;
        bound.push((coordinate, entry.curve.scalar.clone()));
    }
    dae.system.storage = Box::new(CurveStorage {
        linear: QuadraticStorage::new(matrix, n)?,
        curves: bound,
    });
    Ok(dae)
}

impl CircuitDae {
    /// Actual capacitor voltages [V] at stored charges, in capacitor branch
    /// order. These are endpoint constitutive efforts, not the discrete-gradient
    /// mean voltages used by an energy-preserving finite step.
    pub fn capacitor_voltages(&self, state: &[f64]) -> Result<Vec<f64>, CircuitError> {
        if state.len() != self.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "capacitor voltages require a complete finite circuit state",
            ));
        }
        let mut effort = vec![0.0; state.len()];
        self.system.storage().gradient(state, &mut effort);
        self.charge_index
            .iter()
            .map(|&i| finite(effort[i]))
            .collect()
    }
}
