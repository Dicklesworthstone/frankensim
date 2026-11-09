//! Factored physical load tangents in original conjugate port coordinates.

use super::{
    HarmonicPortParameterResponse, HarmonicResponse, HarmonicResponseError, PowerSum,
    PreparedHarmonicResponse, Term, finite, invalid, poll,
};
use crate::quadratic_fraction;
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd};

mod mixed;
pub use mixed::HarmonicPortFrequencyParameterDirection;
pub(super) use mixed::frequency_parameter_action;

/// Authored load derivatives along the same physical scalar parameter.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPortParameterDirection<'a> {
    /// Complete real row-major dB/dtheta, n*ports entries in the original
    /// force/flow attachment units. The same map supplies B*effort and B^T*v.
    pub force_map: &'a [f64],
    /// Complete complex row-major dZ/dtheta, ports*ports entries in the
    /// original conjugate port units. Signed, mutual and nonreciprocal
    /// directions are allowed; this tangent is not itself a passive load.
    pub impedance: &'a [C64],
}

impl Term {
    const ZERO: Self = Self {
        fraction: Dd::ZERO,
        exponent: 0,
    };

    fn abs(self) -> Self {
        Self {
            fraction: self.fraction.abs(),
            exponent: self.exponent,
        }
    }

    fn times(self, factor: f64, shift: i32) -> Self {
        if factor == 0.0 || self.fraction == Dd::ZERO {
            return Self::ZERO;
        }
        let (fraction, exponent) = quadratic_fraction(factor);
        Self {
            fraction: self.fraction * Dd::from_f64(fraction.copysign(factor)),
            exponent: self.exponent + exponent + shift,
        }
    }
}

impl PowerSum {
    // Keep a contracted action normalized even when its physical magnitude
    // exceeds f64. Later multiplication by B or cancellation can be finite.
    fn term(self) -> Result<Term, HarmonicResponseError> {
        if self.fraction == Dd::ZERO {
            return Ok(Term::ZERO);
        }
        if !self.fraction.hi.is_finite() || !self.fraction.lo.is_finite() {
            return Err(invalid("normalized port action is not representable"));
        }
        let exponent = self.exponent + quadratic_fraction(self.fraction.hi).1;
        let original = Term {
            fraction: self.fraction,
            exponent: self.exponent,
        };
        Ok(Term {
            fraction: original.scaled(exponent)?,
            exponent,
        })
    }
}

#[derive(Default)]
struct ComplexSum([PowerSum; 2]);

impl ComplexSum {
    fn real_product(&mut self, factor: f64, value: C64) -> Result<(), HarmonicResponseError> {
        self.0[0].add(Term::product([factor, value.re, 1.0, 1.0], 0))?;
        self.0[1].add(Term::product([factor, value.im, 1.0, 1.0], 0))
    }

    fn complex_product(&mut self, a: C64, b: C64) -> Result<(), HarmonicResponseError> {
        self.0[0].add(Term::product([a.re, b.re, 1.0, 1.0], 0))?;
        self.0[0].add(Term::product([-a.im, b.im, 1.0, 1.0], 0))?;
        self.0[1].add(Term::product([a.re, b.im, 1.0, 1.0], 0))?;
        self.0[1].add(Term::product([a.im, b.re, 1.0, 1.0], 0))
    }

    fn scaled_product(&mut self, a: C64, b: [Term; 2]) -> Result<(), HarmonicResponseError> {
        self.0[0].add(b[0].times(a.re, 0))?;
        self.0[0].add(b[1].times(-a.im, 0))?;
        self.0[1].add(b[1].times(a.re, 0))?;
        self.0[1].add(b[0].times(a.im, 0))
    }

    fn terms(self) -> Result<[Term; 2], HarmonicResponseError> {
        let [real, imaginary] = self.0;
        Ok([real.term()?, imaginary.term()?])
    }

    fn value(self) -> Result<C64, HarmonicResponseError> {
        let [real, imaginary] = self.0;
        Ok(C64::new(real.value()?, imaginary.value()?))
    }
}

pub(super) fn admit_direction(
    owner: &PreparedHarmonicResponse,
    direction: Option<HarmonicPortParameterDirection<'_>>,
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    let Some(direction) = direction else {
        return Ok(());
    };
    let load = owner
        .port_load
        .as_ref()
        .ok_or_else(|| invalid("port parameter direction requires an attached physical load"))?;
    if direction.force_map.len() != load.force_map.len()
        || direction.impedance.len() != load.impedance.len()
    {
        return Err(invalid(
            "complete n*ports attachment and ports*ports impedance directions required",
        ));
    }
    for row in direction.force_map.chunks_exact(load.ports) {
        poll(cx)?;
        if row.iter().any(|value| !value.is_finite()) {
            return Err(invalid("finite physical attachment direction required"));
        }
    }
    for row in direction.impedance.chunks_exact(load.ports) {
        poll(cx)?;
        if row.iter().any(|&value| !finite(value)) {
            return Err(invalid("finite physical impedance direction required"));
        }
    }
    Ok(())
}

// Only the direct effort and its uncancelled envelope survive preparation.
// Intermediate flow/effort components retain binary exponents throughout;
// physical conversion is reserved for the final published port derivatives.
pub(super) struct PortAction<'a> {
    map_direction: Option<&'a [f64]>,
    effort: Vec<[Term; 2]>,
    envelope: Vec<Term>,
}

impl PortAction<'_> {
    pub(super) fn ports(&self) -> usize {
        self.effort.len()
    }

    pub(super) fn visit_rhs(
        &self,
        owner: &PreparedHarmonicResponse,
        response: &HarmonicResponse,
        cx: &Cx<'_>,
        visit: impl FnMut(usize, usize, Term, Term) -> Result<(), HarmonicResponseError>,
    ) -> Result<(), HarmonicResponseError> {
        let Some(primal) = &response.port_load else {
            return Err(invalid(
                "prepared port action lost its original physical load",
            ));
        };
        self.visit_projected_rhs(owner, &primal.effort, cx, visit)
    }

    // Share the original factored force projection between parameter and
    // mixed-frequency directions. The caller chooses the matching physical
    // effort (e or e_omega), without constructing a fictitious response.
    fn visit_projected_rhs(
        &self,
        owner: &PreparedHarmonicResponse,
        effort: &[C64],
        cx: &Cx<'_>,
        mut visit: impl FnMut(usize, usize, Term, Term) -> Result<(), HarmonicResponseError>,
    ) -> Result<(), HarmonicResponseError> {
        let Some(load) = &owner.port_load else {
            return Err(invalid(
                "prepared port action lost its original physical load",
            ));
        };
        let p = load.ports;
        if effort.len() != p {
            return Err(invalid("one original effort per physical port required"));
        }
        for i in 0..owner.n {
            poll(cx)?;
            let shift = -owner.coordinate_exponents[i];
            for a in 0..p {
                if let Some(map) = self.map_direction {
                    let db = map[i * p + a];
                    for (axis, e) in [effort[a].re, effort[a].im].into_iter().enumerate() {
                        let term = Term::product([-db, e, 1.0, 1.0], shift);
                        visit(i, axis, term, term.abs())?;
                    }
                }
                let b = load.force_map[i * p + a];
                for axis in 0..2 {
                    let bound = if axis == 0 {
                        self.envelope[a].times(b.abs(), shift)
                    } else {
                        Term::ZERO
                    };
                    visit(i, axis, self.effort[a][axis].times(-b, shift), bound)?;
                }
            }
        }
        Ok(())
    }
}

pub(super) fn direct_action<'a>(
    owner: &PreparedHarmonicResponse,
    response: &HarmonicResponse,
    direction: Option<HarmonicPortParameterDirection<'a>>,
    frequency_derivative: f64,
    cx: &Cx<'_>,
) -> Result<Option<PortAction<'a>>, HarmonicResponseError> {
    if direction.is_none() && frequency_derivative == 0.0 {
        return Ok(None);
    }
    let (Some(load), Some(primal)) = (&owner.port_load, &response.port_load) else {
        return Ok(None);
    };
    let p = load.ports;
    let mut flow = Vec::with_capacity(p);
    let mut flow_envelope = Vec::with_capacity(p);
    for a in 0..p {
        poll(cx)?;
        let mut value = ComplexSum::default();
        let mut bound = PowerSum::default();
        for i in 0..owner.n {
            // Direct flow at fixed q includes both attachment motion and the
            // explicit velocity derivative -i*omega_theta*q, even at DC.
            let db = direction.map_or(0.0, |d| d.force_map[i * p + a]);
            let b = load.force_map[i * p + a];
            let v = response.velocity[i];
            let q = response.displacement[i];
            for (axis, terms) in [
                [
                    Term::product([db, v.re, 1.0, 1.0], 0),
                    Term::product([frequency_derivative, b, q.im, 1.0], 0),
                ],
                [
                    Term::product([db, v.im, 1.0, 1.0], 0),
                    Term::product([-frequency_derivative, b, q.re, 1.0], 0),
                ],
            ]
            .into_iter()
            .enumerate()
            {
                for term in terms {
                    value.0[axis].add(term)?;
                    bound.add(term.abs())?;
                }
            }
        }
        flow.push(value.terms()?);
        flow_envelope.push(bound.term()?);
    }
    let mut effort = Vec::with_capacity(p);
    let mut envelope = Vec::with_capacity(p);
    for a in 0..p {
        poll(cx)?;
        let mut value = ComplexSum::default();
        let mut bound = PowerSum::default();
        for b in 0..p {
            let z = load.impedance[a * p + b];
            value.scaled_product(z, flow[b])?;
            for part in [z.re, z.im] {
                bound.add(flow_envelope[b].times(part.abs(), 0))?;
            }
            if let Some(direction) = direction {
                let dz = direction.impedance[a * p + b];
                let y = primal.flow[b];
                value.complex_product(dz, y)?;
                for z in [dz.re, dz.im] {
                    for y in [y.re, y.im] {
                        bound.add(Term::product([z, y, 1.0, 1.0], 0).abs())?;
                    }
                }
            }
        }
        effort.push(value.terms()?);
        envelope.push(bound.term()?);
    }
    Ok(Some(PortAction {
        map_direction: direction.map(|d| d.force_map),
        effort,
        envelope,
    }))
}

pub(super) fn observe(
    owner: &PreparedHarmonicResponse,
    response: &HarmonicResponse,
    velocity: &[C64],
    direction: Option<HarmonicPortParameterDirection<'_>>,
    cx: &Cx<'_>,
) -> Result<Option<HarmonicPortParameterResponse>, HarmonicResponseError> {
    let (Some(load), Some(primal)) = (&owner.port_load, &response.port_load) else {
        return Ok(None);
    };
    let p = load.ports;
    let mut flow = Vec::with_capacity(p);
    for a in 0..p {
        poll(cx)?;
        let mut value = ComplexSum::default();
        for (i, &dv) in velocity.iter().enumerate() {
            value.real_product(load.force_map[i * p + a], dv)?;
            if let Some(direction) = direction {
                value.real_product(direction.force_map[i * p + a], response.velocity[i])?;
            }
        }
        flow.push(value.value()?);
    }
    let mut effort = Vec::with_capacity(p);
    let mut power = PowerSum::default();
    for a in 0..p {
        poll(cx)?;
        let mut value = ComplexSum::default();
        for b in 0..p {
            value.complex_product(load.impedance[a * p + b], flow[b])?;
            if let Some(direction) = direction {
                value.complex_product(direction.impedance[a * p + b], primal.flow[b])?;
            }
        }
        let derivative = value.value()?;
        power.half_inner(flow[a], primal.effort[a], 1.0)?;
        power.half_inner(primal.flow[a], derivative, 1.0)?;
        effort.push(derivative);
    }
    let power = power.value()?;
    if !power.is_finite() || flow.iter().chain(&effort).any(|&value| !finite(value)) {
        return Err(invalid(
            "physical parameter port tangent or work is not representable",
        ));
    }
    Ok(Some(HarmonicPortParameterResponse {
        flow_derivative: flow,
        effort_derivative: effort,
        average_absorbed_power_derivative: power,
    }))
}
