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
//! this module only supplies its constitutive energy and analytic gradient.

use super::{
    Branch, CircuitDae, CircuitError, CircuitGraph, MutualInductance,
    assemble_circuit_with_mutual_inductance,
};
use crate::{QuadraticStorage, Storage, linear_exact::QuadraticPart};

/// Maximum retained knots per capacitor, bounding each crossed-segment sum.
pub const MAX_CAPACITOR_KNOTS: usize = 4096;

fn invalid(what: &'static str) -> CircuitError {
    CircuitError::Invalid { what }
}

fn finite(value: f64) -> Result<f64, CircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "capacitor constitutive arithmetic is not representable",
        ))
    }
}

/// Immutable V(q), with strictly increasing charge [C] and voltage [V].
///
/// An explicit (0,0) knot sets the uncharged minimum of this reversible energy.
/// Both polarities need support so an initially uncharged element admits the
/// original centered Newton probes. Asymmetry between polarities is allowed.
#[derive(Clone, Debug, PartialEq)]
pub struct CapacitorChargeCurve {
    knots: Vec<(f64, f64)>,
    inverse_capacitance: Vec<f64>,
    energy_at_knots: Vec<f64>,
}

impl CapacitorChargeCurve {
    /// Admit 3..=MAX_CAPACITOR_KNOTS finite ordered knots, including (0,0)
    /// strictly inside the support. Every segment must have a positive,
    /// representable differential capacitance and stored energy.
    pub fn try_new(knots: Vec<(f64, f64)>) -> Result<Self, CircuitError> {
        if !(3..=MAX_CAPACITOR_KNOTS).contains(&knots.len())
            || knots.iter().any(|(q, v)| !q.is_finite() || !v.is_finite())
        {
            return Err(invalid(
                "bounded finite capacitor charge/voltage knots required",
            ));
        }
        let zero = knots
            .iter()
            .position(|&(q, v)| q == 0.0 && v == 0.0)
            .ok_or_else(|| invalid("capacitor curve requires an explicit uncharged (0,0) knot"))?;
        if zero == 0 || zero + 1 == knots.len() {
            return Err(invalid(
                "capacitor curve must support both charge polarities",
            ));
        }
        let mut inverse_capacitance = Vec::with_capacity(knots.len() - 1);
        for pair in knots.windows(2) {
            let dq = finite(pair[1].0 - pair[0].0)?;
            let dv = finite(pair[1].1 - pair[0].1)?;
            let slope = dv / dq;
            if dq <= 0.0
                || dv <= 0.0
                || !slope.is_finite()
                || slope <= 0.0
                || !slope.recip().is_finite()
                || slope.recip() <= 0.0
            {
                return Err(invalid(
                    "capacitor knots require positive representable differential capacitance",
                ));
            }
            inverse_capacitance.push(slope);
        }
        let mut energy_at_knots = vec![0.0; knots.len()];
        for i in zero..knots.len() - 1 {
            let area = f64::midpoint(knots[i].1, knots[i + 1].1) * (knots[i + 1].0 - knots[i].0);
            energy_at_knots[i + 1] = finite(energy_at_knots[i] + area)?;
            if area <= 0.0 || energy_at_knots[i + 1] <= 0.0 {
                return Err(invalid("capacitor energy is not representable"));
            }
        }
        for i in (0..zero).rev() {
            let area = f64::midpoint(knots[i].1, knots[i + 1].1) * (knots[i + 1].0 - knots[i].0);
            energy_at_knots[i] = finite(energy_at_knots[i + 1] - area)?;
            if area >= 0.0 || energy_at_knots[i] <= 0.0 {
                return Err(invalid("capacitor energy is not representable"));
            }
        }
        Ok(Self {
            knots,
            inverse_capacitance,
            energy_at_knots,
        })
    }

    /// Original (charge [C], voltage [V]) knots, without resampling.
    pub fn knots(&self) -> &[(f64, f64)] {
        &self.knots
    }

    /// Closed source support [C]; stepping may need room for Newton probes.
    pub fn charge_range_c(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.knots.len() - 1].0)
    }

    fn segment(&self, charge: f64) -> Result<usize, CircuitError> {
        let (lo, hi) = self.charge_range_c();
        if !charge.is_finite() || charge < lo || charge > hi {
            return Err(invalid(
                "capacitor charge is outside the supplied curve support",
            ));
        }
        Ok(self
            .knots
            .partition_point(|&(q, _)| q <= charge)
            .saturating_sub(1)
            .min(self.knots.len() - 2))
    }

    /// Actual terminal voltage dH/dq at the stored charge.
    pub fn voltage_v(&self, charge_c: f64) -> Result<f64, CircuitError> {
        let i = self.segment(charge_c)?;
        if charge_c == self.knots[i + 1].0 {
            return Ok(self.knots[i + 1].1);
        }
        let anchor = if charge_c < 0.0 { i + 1 } else { i };
        let value = finite(
            self.inverse_capacitance[i]
                .mul_add(charge_c - self.knots[anchor].0, self.knots[anchor].1),
        )?;
        if charge_c != 0.0 && value == 0.0 {
            return Err(invalid("capacitor voltage is not representable"));
        }
        Ok(value)
    }

    /// Initialize charge from voltage using the SAME monotone relation.
    pub fn charge_c(&self, voltage_v: f64) -> Result<f64, CircuitError> {
        if !voltage_v.is_finite()
            || voltage_v < self.knots[0].1
            || voltage_v > self.knots[self.knots.len() - 1].1
        {
            return Err(invalid(
                "capacitor voltage is outside the supplied curve support",
            ));
        }
        let i = self
            .knots
            .partition_point(|&(_, v)| v <= voltage_v)
            .saturating_sub(1)
            .min(self.knots.len() - 2);
        if voltage_v == self.knots[i + 1].1 {
            return Ok(self.knots[i + 1].0);
        }
        let anchor = if voltage_v < 0.0 { i + 1 } else { i };
        let charge = finite(
            (voltage_v - self.knots[anchor].1) / self.inverse_capacitance[i] + self.knots[anchor].0,
        )?;
        if voltage_v != 0.0 && charge == 0.0 {
            return Err(invalid("capacitor charge is not representable"));
        }
        self.segment(charge)?;
        Ok(charge)
    }

    /// dq/dV [F]. At an interior knot the right segment is returned; the last
    /// knot uses its left segment. No derivative continuity is claimed there.
    pub fn differential_capacitance_f(&self, charge_c: f64) -> Result<f64, CircuitError> {
        Ok(self.inverse_capacitance[self.segment(charge_c)?].recip())
    }

    /// Integral from zero to q of V(s) ds [J], with the uncharged datum fixed.
    pub fn stored_energy_j(&self, charge_c: f64) -> Result<f64, CircuitError> {
        let i = self.segment(charge_c)?;
        // Evaluate from the knot nearer zero: subtracting almost all the
        // energy of a negative knot would erase small negative-charge energy.
        let anchor = if charge_c < 0.0 { i + 1 } else { i };
        let dq = charge_c - self.knots[anchor].0;
        let voltage = self.voltage_v(charge_c)?;
        let energy = finite(
            self.energy_at_knots[anchor] + dq * f64::midpoint(self.knots[anchor].1, voltage),
        )?;
        if energy < 0.0 || charge_c != 0.0 && energy == 0.0 {
            return Err(invalid("capacitor stored energy is not representable"));
        }
        Ok(energy)
    }

    // Integral_base^end (V(q)-V(base)) dq. Walking outwards from base adds
    // positive segment areas without subtracting large biased endpoint energies.
    fn relative_energy(&self, base: f64, end: f64) -> Result<f64, CircuitError> {
        let mut i = self.segment(base)?;
        self.segment(end)?;
        let mut position = base;
        let mut voltage_change = 0.0;
        let mut area = 0.0;
        while position != end {
            if end < position && position == self.knots[i].0 {
                i -= 1;
            }
            let next = if end > position {
                end.min(self.knots[i + 1].0)
            } else {
                end.max(self.knots[i].0)
            };
            let dq = next - position;
            let dv = finite(self.inverse_capacitance[i] * dq)?;
            area = finite(area + dq * (voltage_change + 0.5 * dv))?;
            voltage_change = finite(voltage_change + dv)?;
            position = next;
            if position != end && end > position {
                i += 1;
            }
        }
        Ok(area)
    }

    fn midpoint_remainder(
        &self,
        before: f64,
        after: f64,
        midpoint: f64,
    ) -> Result<f64, CircuitError> {
        let a = self.segment(before)?;
        let b = self.segment(after)?;
        let m = self.segment(midpoint)?;
        if a == b && a == m {
            // Includes the supplied represented midpoint's rounding error.
            return finite(
                (after - before)
                    * (0.5 * self.inverse_capacitance[m])
                    * ((before - midpoint) + (after - midpoint)),
            );
        }
        finite(self.relative_energy(midpoint, after)? - self.relative_energy(midpoint, before)?)
    }

    /// Stable energy difference across any retained segments [J].
    pub fn energy_change_j(&self, before_c: f64, after_c: f64) -> Result<f64, CircuitError> {
        let midpoint = f64::midpoint(before_c, after_c);
        finite(
            self.voltage_v(midpoint)? * (after_c - before_c)
                + self.midpoint_remainder(before_c, after_c, midpoint)?,
        )
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
    curves: Vec<(usize, CapacitorChargeCurve)>,
}

impl Storage for CurveStorage {
    fn coordinate_scale(&self, coordinate: usize) -> Option<f64> {
        self.curves
            .iter()
            .find(|(i, _)| *i == coordinate)
            .map(|(_, curve)| {
                let (lo, hi) = curve.charge_range_c();
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
                sum + curve.stored_energy_j(x[*i]).unwrap_or(f64::NAN)
            })
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.linear.gradient(x, out);
        for (i, curve) in &self.curves {
            out[*i] = curve.voltage_v(x[*i]).unwrap_or(f64::NAN);
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
    let mut ordered: Vec<_> = curves.iter().collect();
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
    let mut dae = assemble_circuit_with_mutual_inductance(graph, mutual)?;
    if ordered.is_empty() {
        return Ok(dae);
    }
    let n = dae.system.state_dim();
    let mut matrix = dae
        .system
        .storage()
        .exact_quadratic_matrix()
        .expect("native circuit storage is quadratic")
        .to_vec();
    let mut bound = Vec::with_capacity(ordered.len());
    for entry in ordered {
        let slot = graph.branches[..entry.branch]
            .iter()
            .filter(|b| matches!(b.2, Branch::Capacitor { .. }))
            .count();
        let coordinate = dae.charge_index[slot];
        matrix[coordinate * n + coordinate] = 0.0;
        bound.push((coordinate, entry.curve.clone()));
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
