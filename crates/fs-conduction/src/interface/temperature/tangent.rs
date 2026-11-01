//! Exact linearization of admitted per-face mean or two-sided contact laws.
//! The original matching producer supplies its vertex correspondence and mass
//! rule. This adds dK/dT * T; it does not replace K or finite-difference a solve.

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_matdb::{InterpolationPolicy, PreparedSurfaceResponse, PropertyValue};
use fs_sparse::{Coo, Csr, ops};
use crate::{ConductionError, ThermalBoundary};
use super::{InterfaceResistance, TemperatureContact, TemperatureContactCoordinates,
    TemperatureContactEvaluation, TemperatureContacts, finite, invalid, poll, positive};
use crate::transient::backward_euler::StepLinearization;

/// A rechecked capacity endpoint with complete R(T) contact derivatives.
/// Existing source, history, capacity and Robin operations use the SAME
/// complete nonsymmetric response. No frozen-contact solve is exposed here.
pub struct ContactCapacityLinearization<'a> {
    step: StepLinearization<'a>,
    contact: ContactTemperatureLinearization,
}

impl<'a> std::ops::Deref for ContactCapacityLinearization<'a> {
    type Target = StepLinearization<'a>;
    fn deref(&self) -> &Self::Target { &self.step }
}

/// One complete physical transpose and its named contact-scale contractions.
#[derive(Debug)]
pub struct ContactCapacityGradient {
    /// Existing watt-coordinate source, history and boundary multiplier.
    pub transport: crate::adjoint::robin::RobinGradient,
    /// Independent multipliers of each interface's entire resistance law at one.
    pub resistance_multipliers: Vec<crate::interface::ContactResistanceGradient>,
}

impl<'a> ContactCapacityLinearization<'a> {
    pub(crate) fn new(step: StepLinearization<'a>, contact: ContactTemperatureLinearization) -> Self {
        Self { step, contact }
    }
    /// Physical local source states and the complete contact-only Jacobian.
    #[must_use]
    pub fn contact(&self) -> &ContactTemperatureLinearization { &self.contact }
    /// Differentiate nodal temperatures and selected Robin observations once,
    /// then contract each whole resistance-law multiplier with that same dual.
    /// Geometry, pressure and source selection stay fixed. No extra dt belongs
    /// in these watt-coordinate contractions; named controls remain independent
    /// even when their source cards match.
    pub fn pullback_with_resistance_controls(
        &self, cx: &Cx<'_>, temperature_weights: &[f64], mean_weights: &[f64],
        heat_weights: &[f64], max_face_pairs: usize,
    ) -> Result<ContactCapacityGradient, ConductionError> {
        poll(cx, 0)?;
        if self.contact.face_pairs > max_face_pairs {
            return Err(invalid("<contact-adjoint>", "contact controls exceed the paired-face allowance"));
        }
        let transport = self.step.pullback(cx, temperature_weights, mean_weights, heat_weights)?;
        let resistance_multipliers = self.contact.physical_contacts().frozen_interfaces()
            .matching_resistance_scale_pullback(
                cx, self.contact.temperature(), &transport.nodal_load, max_face_pairs,
            )?;
        Ok(ContactCapacityGradient { transport, resistance_multipliers })
    }
}

enum ContactDerivative<'a> {
    Constant,
    Curve { axis: &'a str, knots: &'a [(f64, f64)] },
    Surface(PreparedSurfaceResponse),
}

impl TemperatureContact<'_> {
    fn axis_weights(&self, axis: &str) -> [f64; 2] {
        match &self.coordinates {
            TemperatureContactCoordinates::Mean(temperature) => {
                if axis == temperature { [0.5; 2] } else { [0.0; 2] }
            }
            TemperatureContactCoordinates::TwoSided(axes) =>
                axes.each_ref().map(|temperature| if axis == temperature { 1.0 } else { 0.0 }),
        }
    }

    // Prepare each original source law once per complete linearization, not
    // once per face. Fixed external axes require no artificial cell tangent:
    // even an isolated exact source sample can be constant in temperature.
    fn derivative_response(&self) -> Result<ContactDerivative<'_>, ConductionError> {
        let claim = self.card.claims().claim(self.resistance.receipt().selected)
            .ok_or_else(|| invalid(self.name(), "the selected contact claim is missing"))?;
        match (&claim.value, &claim.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) =>
                Ok(ContactDerivative::Constant),
            (PropertyValue::Curve { abscissa, knots, .. }, InterpolationPolicy::LinearInside) => {
                if self.axis_weights(abscissa) == [0.0; 2] { return Ok(ContactDerivative::Constant); }
                Ok(ContactDerivative::Curve { axis: abscissa, knots })
            }
            (PropertyValue::Surface2D { axes, .. }, InterpolationPolicy::LinearInside) => {
                if axes.iter().all(|axis| self.axis_weights(axis) == [0.0; 2]) {
                    return Ok(ContactDerivative::Constant);
                }
                claim.value.prepare_surface_response().map(ContactDerivative::Surface)
                    .map_err(|error| invalid(self.name(), format!("contact surface derivative refused: {error}")))
            }
            _ => Err(invalid(self.name(), "a contact derivative requires a supported constant scalar, linear curve or original linear surface")),
        }
    }

    fn temperature_partials(
        &self, cx: &Cx<'_>, temperatures: [f64; 2], response: &ContactDerivative<'_>,
    ) -> Result<[f64; 2], ConductionError> {
        let point = self.query_point_at(temperatures)?;
        let axes = match &self.coordinates {
            TemperatureContactCoordinates::Mean(axis) => std::slice::from_ref(axis),
            TemperatureContactCoordinates::TwoSided(axes) => axes.as_slice(),
        };
        // The actual point already passed the physical pinned query. Its
        // neighbours admit all declared thermal validity coordinates; they do
        // not estimate a derivative. Perturb the source coordinate directly:
        // half-weighting a trace ULP could otherwise round back to the mean.
        for axis in axes {
            let t = *point.axes().get(axis)
                .ok_or_else(|| invalid(self.name(), "contact derivative lost a thermal query coordinate"))?;
            for temperature in [t.next_down(), t.next_up()] {
                poll(cx, 0)?;
                let neighbour = point.clone().with(axis, positive(temperature)?)
                    .map_err(|error| invalid(self.name(), format!("contact derivative query refused: {error}")))?;
                InterfaceResistance::from_card_pinned(
                    self.name(), self.card, &neighbour, self.resistance.receipt().selected,
                )?;
            }
        }
        match response {
            ContactDerivative::Constant => Ok([0.0; 2]),
            ContactDerivative::Curve { axis, knots } => {
                let t = *point.axes().get(*axis)
                    .ok_or_else(|| invalid(self.name(), "contact curve derivative lost its source coordinate"))?;
                let mut slope = None;
                for (index, pair) in knots.windows(2).enumerate() {
                    if index % 256 == 0 { poll(cx, index)?; }
                    if pair[0].0 <= t && t <= pair[1].0 {
                        let value = finite((pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0))?;
                        if slope.is_some_and(|previous| previous != value) {
                            return Err(invalid(self.name(), "contact resistance has a slope corner at the evaluated temperature"));
                        }
                        slope = Some(value);
                    }
                }
                let slope = slope.ok_or_else(|| invalid(self.name(), "contact temperature is outside the selected curve segments"))?;
                Ok(self.axis_weights(axis).map(|weight| weight * slope))
            }
            ContactDerivative::Surface(surface) => {
                let mut coordinates = [0.0; 2];
                let mut varying = [false; 2];
                for (i, axis) in surface.axes().iter().enumerate() {
                    coordinates[i] = *point.axes().get(axis)
                        .ok_or_else(|| invalid(self.name(), "contact surface derivative lost its source coordinate"))?;
                    varying[i] = self.axis_weights(axis) != [0.0; 2];
                }
                let gradient = surface.gradient_for_axes(coordinates, varying)
                    .map_err(|error| invalid(self.name(), format!("contact surface derivative refused: {error}")))?;
                let mut partials = [0.0; 2];
                for (axis, derivative) in surface.axes().iter().zip(gradient) {
                    for (partial, weight) in partials.iter_mut().zip(self.axis_weights(axis)) {
                        *partial = finite(*partial + finite(weight * derivative)?)?;
                    }
                }
                Ok(partials)
            }
        }
    }
}

/// Physical contact state and both actions of its complete temperature tangent.
/// Units: loads are W and tangent entries W/K. Geometry and every nonthermal
/// coordinate are fixed. This is a constitutive derivative, not a PDE solve,
/// invertibility certificate or derivative of material-source selection.
#[derive(Debug)]
pub struct ContactTemperatureLinearization {
    face_pairs: usize,
    physical: TemperatureContactEvaluation,
    temperature: Vec<f64>,
    load: Vec<f64>,
    frozen: Csr,
    feedback: Csr,
    feedback_transpose: Csr,
}

impl ContactTemperatureLinearization {
    /// Actual source-backed state at the unchanged linearization point.
    #[must_use]
    pub fn physical_contacts(&self) -> &TemperatureContactEvaluation { &self.physical }
    /// Temperature field at which all local queries and derivatives were made.
    #[must_use]
    pub fn temperature(&self) -> &[f64] { &self.temperature }
    /// Nodal outward contact loads, internally conservative; not external heat.
    #[must_use]
    pub fn nodal_outflow_w(&self) -> &[f64] { &self.load }
    /// Complete J_contact*dT, including both trace-temperature partials at every face.
    pub fn apply(&self, cx: &Cx<'_>, direction: &[f64]) -> Result<Vec<f64>, ConductionError> {
        self.action(cx, direction, false)
    }
    /// Complete J_contact^T*z, not the symmetric frozen-resistance operator.
    pub fn apply_transpose(&self, cx: &Cx<'_>, weights: &[f64]) -> Result<Vec<f64>, ConductionError> {
        self.action(cx, weights, true)
    }
    fn action(&self, cx: &Cx<'_>, input: &[f64], transpose: bool) -> Result<Vec<f64>, ConductionError> {
        poll(cx, 0)?;
        if input.len() != self.temperature.len() {
            return Err(invalid("<contact-tangent>", "contact tangent field shape differs"));
        }
        for (i, &value) in input.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            finite(value)?;
        }
        let mut result = vec![0.0; input.len()];
        let mut extra = vec![0.0; input.len()];
        // The original fixed-contact operator is symmetric by construction.
        self.frozen.spmv(input, &mut result);
        if transpose { self.feedback_transpose.spmv(input, &mut extra); }
        else { self.feedback().spmv(input, &mut extra); }
        for (i, (value, extra)) in result.iter_mut().zip(extra).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            *value = finite(*value + extra)?;
        }
        poll(cx, input.len())?;
        Ok(result)
    }
    /// Only checked in-crate transient binders may attach this physical term.
    pub(crate) fn feedback(&self) -> &Csr { &self.feedback }
}

impl TemperatureContacts<'_, '_> {
    /// Linearize local source laws at the supplied physical field without
    /// advancing it. Every source must admit a two-sided temperature
    /// neighbourhood and unique partials, including for an eventual zero seed.
    /// Equal-slope knots and equal-gradient internal source edges are admitted;
    /// a genuine source corner or thermal domain end refuses. Fixed external
    /// coordinates do not acquire artificial derivative requirements.
    ///
    /// `max_face_pairs` bounds ALL fixed and variable faces. At most 36 extra
    /// entries per variable pair are staged before sparse assembly; both
    /// orientations and nodal action workspace have checked size products.
    /// This bound is not a measurement of allocator/source-card storage.
    pub fn linearize_temperature(
        &self, cx: &Cx<'_>, boundary: &ThermalBoundary, temperature: &[f64], max_face_pairs: usize,
    ) -> Result<ContactTemperatureLinearization, ConductionError> {
        poll(cx, 0)?;
        let count = self.fixed.iter().try_fold(self.variable_pairs, |n, surface|
            n.checked_add(surface.face_pairs().len()))
            .ok_or_else(|| invalid("<contact-tangent>", "contact face allowance overflows"))?;
        let n = self.mesh.vertex_count();
        if count > max_face_pairs || count.checked_mul(144 * 24)
            .and_then(|bytes| n.checked_mul(8 * 8).and_then(|nodes| bytes.checked_add(nodes)))
            .is_none_or(|bytes| bytes > isize::MAX as usize)
        {
            return Err(invalid("<contact-tangent>", "contact tangent exceeds its paired-face or addressable storage allowance"));
        }
        let physical = self.evaluate(cx, boundary, temperature)?;
        let mut responses = BTreeMap::new();
        for (i, law) in self.variable.iter().enumerate() {
            poll(cx, i)?;
            responses.insert(law.name(), (law, law.derivative_response()?));
        }
        let mut partials = BTreeMap::new();
        for (i, state) in physical.states().iter().enumerate() {
            poll(cx, i)?;
            let (law, response) = responses.get(state.interface.as_str())
                .ok_or_else(|| invalid("<contact-tangent>", "local contact state lost its source owner"))?;
            let derivative = law.temperature_partials(cx,
                [state.side_a_temperature_k, state.side_b_temperature_k], response)?;
            if partials.insert(state.pair, derivative).is_some() {
                return Err(invalid(law.name(), "duplicate local contact derivative"));
            }
        }
        let feedback = physical.interfaces.matching.temperature_feedback(
            cx, temperature, &partials, max_face_pairs,
        )?;
        let mut coo = Coo::new(n, n);
        physical.interfaces.assemble_into(cx, &mut coo)?;
        let frozen = coo.assemble();
        let feedback_transpose = ops::transpose(&feedback);
        let mut load = vec![0.0; n];
        frozen.spmv(temperature, &mut load);
        for (i, &value) in load.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            finite(value)?;
        }
        poll(cx, n)?;
        Ok(ContactTemperatureLinearization {
            face_pairs: count, physical, temperature: temperature.to_vec(), load, frozen, feedback, feedback_transpose,
        })
    }
}
