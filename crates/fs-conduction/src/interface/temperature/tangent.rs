//! Exact linearization of the admitted per-face R(T_interface) contact law.
//! The original matching producer supplies its vertex correspondence and mass
//! rule. This adds dK/dT * T; it does not replace K or finite-difference a solve.

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_matdb::{InterpolationPolicy, PropertyValue};
use fs_sparse::{Coo, Csr, ops};
use crate::{ConductionError, ThermalBoundary};
use super::{TemperatureContact, TemperatureContactEvaluation, TemperatureContacts, finite, invalid, poll};

impl TemperatureContact<'_> {
    fn temperature_slope(&self, cx: &Cx<'_>, t: f64) -> Result<f64, ConductionError> {
        // All three evaluations keep the original pin and nonthermal context.
        // Neighbours admit the source domain; they do NOT estimate the slope.
        for temperature in [t.next_down(), t, t.next_up()] {
            poll(cx, 0)?;
            self.at(temperature)?;
        }
        let claim = self.card.claims().claim(self.resistance.receipt().selected)
            .ok_or_else(|| invalid(self.name(), "the selected contact claim is missing"))?;
        match (&claim.value, &claim.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => Ok(0.0),
            (PropertyValue::Curve { abscissa, knots, .. }, InterpolationPolicy::LinearInside) => {
                // A curve in a fixed nonthermal coordinate is constant in T.
                if abscissa != &self.temperature_axis { return Ok(0.0); }
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
                slope.ok_or_else(|| invalid(self.name(), "contact temperature is outside the selected curve segments"))
            }
            _ => Err(invalid(self.name(), "a contact derivative requires a continuously supported constant scalar or linear curve")),
        }
    }
}

/// Physical contact state and both actions of its complete temperature tangent.
/// Units: loads are W and tangent entries W/K. Geometry and every nonthermal
/// coordinate are fixed. This is a constitutive derivative, not a PDE solve,
/// invertibility certificate or derivative of material-source selection.
#[derive(Debug)]
pub struct ContactTemperatureLinearization {
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
    /// Complete J_contact*dT, including R'(T_interface) at every face.
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
    /// neighbourhood and a unique slope, including for an eventual zero seed.
    /// A redundant equal-slope knot is admitted; a corner or domain end refuses.
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
        let mut slopes = BTreeMap::new();
        for (i, state) in physical.states().iter().enumerate() {
            poll(cx, i)?;
            let law = self.variable.iter().find(|law| law.name() == state.interface)
                .ok_or_else(|| invalid("<contact-tangent>", "local contact state lost its source owner"))?;
            let slope = law.temperature_slope(cx, state.interface_temperature_k)?;
            if slopes.insert(state.pair, slope).is_some() {
                return Err(invalid(law.name(), "duplicate local contact derivative"));
            }
        }
        let feedback = physical.interfaces.matching.temperature_feedback(
            cx, temperature, &slopes, max_face_pairs,
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
            physical, temperature: temperature.to_vec(), load, frozen, feedback, feedback_transpose,
        })
    }
}
