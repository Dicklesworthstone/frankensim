//! Exact linearization of admitted face-mean or pointwise contact quadrature.
//! Matching faces retain their declared mass or pointwise rule; nonmatching
//! faces retain common-refinement samples. This adds dK/dT * T; it does not
//! replace K, change the quadrature, or finite-difference a solve.

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_matdb::{InterpolationPolicy, PreparedSurfaceResponse, PropertyValue};
use fs_sparse::{Coo, Csr, ops};
use crate::{ConductionError, ThermalBoundary};
use super::{InterfaceResistance, TemperatureContact, TemperatureContactCoordinates,
    TemperatureContactEvaluation, TemperatureContacts, finite, invalid, poll};
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
            .resistance_scale_pullback(
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

#[derive(Clone, Copy)]
enum DerivativePolicy {
    TwoSided,
    PrimalSemismooth,
}

impl TemperatureContact<'_> {
    // Ordered partial coordinates are T_A, T_B, f_A, f_B. Mean T, when
    // explicitly present, contributes half to each independent thermal trace.
    fn axis_weights(&self, axis: &str) -> [f64; 4] {
        if let Some(material) = self.material_states {
            return std::array::from_fn(|i| if i < 2 {
                if axis == material.point.temperature_axes[i] { 1.0 }
                else if material.point.mean_temperature_axis == Some(axis) { 0.5 } else { 0.0 }
            } else if material.point.liquid_fraction_axes[i-2] == Some(axis) { 1.0 } else { 0.0 });
        }
        let thermal = match &self.coordinates {
            TemperatureContactCoordinates::Mean(temperature) =>
                if axis == temperature { [0.5; 2] } else { [0.0; 2] },
            TemperatureContactCoordinates::TwoSided(axes) =>
                axes.each_ref().map(|temperature| if axis == temperature { 1.0 } else { 0.0 }),
        };
        [thermal[0], thermal[1], 0.0, 0.0]
    }

    // Each distinct pinned response is prepared once, not once per face.
    fn derivative_response(&self, selected: fs_matdb::ClaimId) -> Result<ContactDerivative<'_>, ConductionError> {
        let claim = self.card.claims().claim(selected)
            .ok_or_else(|| invalid(self.name(), "the selected contact claim is missing"))?;
        match (&claim.value, &claim.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) =>
                Ok(ContactDerivative::Constant),
            (PropertyValue::Curve { abscissa, knots, .. }, InterpolationPolicy::LinearInside) => {
                if self.axis_weights(abscissa) == [0.0; 4] { return Ok(ContactDerivative::Constant); }
                Ok(ContactDerivative::Curve { axis: abscissa, knots })
            }
            (PropertyValue::Surface2D { axes, .. }, InterpolationPolicy::LinearInside) => {
                if axes.iter().all(|axis| self.axis_weights(axis) == [0.0; 4]) {
                    return Ok(ContactDerivative::Constant);
                }
                claim.value.prepare_surface_response().map(ContactDerivative::Surface)
                    .map_err(|error| invalid(self.name(), format!("contact surface derivative refused: {error}")))
            }
            _ => Err(invalid(self.name(), "a contact derivative requires a supported constant scalar, linear curve or original linear surface")),
        }
    }

    fn state_neighbourhood(&self, cx: &Cx<'_>, side_a: Option<usize>, temperatures: [f64; 2],
        fractions: Option<[f64; 2]>, active: [bool; 4]) -> Result<(), ConductionError>
    {
        let point = self.query_point_at_state(side_a, temperatures, fractions)?;
        let selected = self.source_context(side_a)?.1.receipt().selected;
        // These are support checks, not numerical derivatives. A phase or
        // temperature coordinate that is locally constant in h needs no
        // invented independent perturbation outside its physical chart.
        for (axis, &value) in point.axes() {
            if !self.axis_weights(axis).iter().zip(active).any(|(&w, a)| w != 0.0 && a) { continue; }
            for coordinate in [value.next_down(), value.next_up()] {
                poll(cx, 0)?;
                let weights = self.axis_weights(axis);
                if weights[0] != 0.0 || weights[1] != 0.0 { super::positive(coordinate)?; }
                let neighbour = point.clone().with(axis, coordinate)
                    .map_err(|error| invalid(self.name(), format!("contact derivative query refused: {error}")))?;
                InterfaceResistance::from_card_pinned(self.name(), self.card, &neighbour, selected)?;
            }
        }
        Ok(())
    }

    fn state_partials(&self, cx: &Cx<'_>, side_a: Option<usize>, temperatures: [f64; 2],
        fractions: Option<[f64; 2]>, active: [bool; 4], response: &ContactDerivative<'_>,
        policy: DerivativePolicy) -> Result<[f64; 4], ConductionError>
    {
        let point = self.query_point_at_state(side_a, temperatures, fractions)?;
        let varying = |axis: &str| self.axis_weights(axis).iter().zip(active)
            .any(|(&weight, active)| weight != 0.0 && active);
        let primal = matches!(policy, DerivativePolicy::PrimalSemismooth)
            && (!matches!(response, ContactDerivative::Surface(_)) || self.material_states.is_some());
        if !primal { self.state_neighbourhood(cx, side_a, temperatures, fractions, active)?; }
        match response {
            ContactDerivative::Constant => Ok([0.0; 4]),
            ContactDerivative::Curve { axis, knots } => {
                if !varying(axis) { return Ok([0.0; 4]); }
                let value = *point.axes().get(*axis)
                    .ok_or_else(|| invalid(self.name(), "contact curve derivative lost its source coordinate"))?;
                if primal {
                    if knots.len() < 2 || value < knots[0].0 || value > knots[knots.len()-1].0 {
                        return Err(invalid(self.name(), "contact coordinate is outside the selected curve segments"));
                    }
                    let upper = knots.partition_point(|knot| knot.0 <= value).min(knots.len()-1);
                    let a = knots[upper-1]; let b = knots[upper];
                    let slope = finite((b.1-a.1)/(b.0-a.0))?;
                    return Ok(self.axis_weights(axis).map(|weight| weight*slope));
                }
                let mut slope = None;
                for (index, pair) in knots.windows(2).enumerate() {
                    if index % 256 == 0 { poll(cx, index)?; }
                    if pair[0].0 <= value && value <= pair[1].0 {
                        let next = finite((pair[1].1-pair[0].1)/(pair[1].0-pair[0].0))?;
                        if slope.is_some_and(|old| old != next) {
                            return Err(invalid(self.name(), "contact resistance has a source-coordinate slope corner"));
                        }
                        slope = Some(next);
                    }
                }
                let slope = slope.ok_or_else(|| invalid(self.name(), "contact coordinate is outside the selected curve segments"))?;
                Ok(self.axis_weights(axis).map(|weight| weight*slope))
            }
            ContactDerivative::Surface(surface) => {
                let mut coordinates = [0.0; 2]; let mut axes_vary = [false; 2];
                for (i, axis) in surface.axes().iter().enumerate() {
                    coordinates[i] = *point.axes().get(axis)
                        .ok_or_else(|| invalid(self.name(), "contact surface derivative lost its source coordinate"))?;
                    axes_vary[i] = varying(axis);
                }
                let gradient = if primal {
                    // Original selected source-cell tangent permits inward
                    // phase motion from f=0/1. Actual residuals enforce support;
                    // this generalized derivative is never an accepted adjoint.
                    let response = surface.evaluate(coordinates).map_err(|error|
                        invalid(self.name(), format!("contact source-cell derivative refused: {error}")))?;
                    std::array::from_fn(|i| if axes_vary[i] { response.gradient[i] } else { 0.0 })
                } else {
                    surface.gradient_for_axes(coordinates, axes_vary).map_err(|error|
                        invalid(self.name(), format!("contact surface derivative refused: {error}")))?
                };
                let mut partials = [0.0; 4];
                for (axis, derivative) in surface.axes().iter().zip(gradient) {
                    for (partial, weight) in partials.iter_mut().zip(self.axis_weights(axis)) {
                        *partial = finite(*partial + finite(weight*derivative)?)?;
                    }
                }
                Ok(partials)
            }
        }
    }
}

/// Physical contact state and its original constitutive derivatives.
/// Temperature-only laws expose complete W/K actions. An explicit material
/// binding instead exposes complete enthalpy actions in W/(J/kg), retaining
/// independent phase feedback when T is constant. Geometry, external conditions
/// and source selection stay fixed; this is not a PDE solve or invertibility proof.
/// The explicitly primal constructor may select a generalized derivative at
/// a source-curve corner; that selection is not a unique physical adjoint.
#[derive(Debug)]
pub struct ContactTemperatureLinearization {
    face_pairs: usize,
    physical: TemperatureContactEvaluation,
    temperature: Vec<f64>,
    load: Vec<f64>,
    frozen: Csr,
    feedback: Csr,
    feedback_transpose: Csr,
    specific_enthalpy: Option<Vec<f64>>,
    enthalpy_slopes: Option<Vec<f64>>,
    phase_feedback: Option<Csr>,
    phase_feedback_transpose: Option<Csr>,
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
    /// Complete temperature-coordinate action for a temperature-only law.
    /// Material-state responses require `apply_enthalpy`: independent thermal
    /// perturbations need not exist on the original latent chart.
    pub fn apply(&self, cx: &Cx<'_>, direction: &[f64]) -> Result<Vec<f64>, ConductionError> {
        if self.phase_feedback.is_some() {
            return Err(invalid("<contact-tangent>", "material contact action requires its original enthalpy coordinates"));
        }
        self.action(cx, direction, false)
    }
    /// Complete temperature transpose for a temperature-only source law.
    /// Material-state responses require `apply_enthalpy_transpose`.
    pub fn apply_transpose(&self, cx: &Cx<'_>, weights: &[f64]) -> Result<Vec<f64>, ConductionError> {
        if self.phase_feedback.is_some() {
            return Err(invalid("<contact-tangent>", "material contact transpose requires its original enthalpy coordinates"));
        }
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
    /// Original h field when this response was bound through the native h owner.
    #[must_use]
    pub fn specific_enthalpy_j_kg(&self) -> Option<&[f64]> { self.specific_enthalpy.as_deref() }

    /// Complete dQ/dh, including phase feedback on latent plateaus.
    pub fn apply_enthalpy(&self, cx: &Cx<'_>, direction: &[f64]) -> Result<Vec<f64>, ConductionError> {
        let slopes = self.enthalpy_slopes.as_ref().ok_or_else(|| invalid("<contact-tangent>",
            "enthalpy action needs a response bound at original enthalpy states"))?;
        if direction.len() != slopes.len() { return Err(invalid("<contact-tangent>", "enthalpy direction shape differs")); }
        let mut scaled = Vec::with_capacity(direction.len());
        for (i, (&value, &slope)) in direction.iter().zip(slopes).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; } scaled.push(finite(value*slope)?);
        }
        let mut result = self.action(cx, &scaled, false)?;
        if let Some(feedback) = &self.phase_feedback {
            let mut extra = vec![0.0; result.len()]; feedback.spmv(direction, &mut extra);
            for (value, extra) in result.iter_mut().zip(extra) { *value = finite(*value+extra)?; }
        }
        poll(cx, result.len())?; Ok(result)
    }

    /// Transpose of the same complete enthalpy-coordinate contact tangent.
    pub fn apply_enthalpy_transpose(&self, cx: &Cx<'_>, weights: &[f64]) -> Result<Vec<f64>, ConductionError> {
        let slopes = self.enthalpy_slopes.as_ref().ok_or_else(|| invalid("<contact-tangent>",
            "enthalpy transpose needs a response bound at original enthalpy states"))?;
        let mut result = self.action(cx, weights, true)?;
        for (value, &slope) in result.iter_mut().zip(slopes) { *value = finite(*value*slope)?; }
        if let Some(feedback) = &self.phase_feedback_transpose {
            let mut extra = vec![0.0; result.len()]; feedback.spmv(weights, &mut extra);
            for (value, extra) in result.iter_mut().zip(extra) { *value = finite(*value+extra)?; }
        }
        poll(cx, result.len())?; Ok(result)
    }

    pub(crate) fn phase_feedback(&self) -> Option<&Csr> { self.phase_feedback.as_ref() }

    /// Matching pairs plus retained overlap triangles in this physical binding.
    pub(crate) fn trace_elements(&self) -> usize { self.face_pairs }

    /// Only checked in-crate endpoint binders may attach this physical term.
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
    /// `max_face_pairs` bounds matching pairs plus overlap triangles. At most
    /// 36 extra entries per matching pair or 108 per overlap triangle are staged;
    /// both
    /// orientations and nodal action workspace have checked size products.
    /// This bound is not a measurement of allocator/source-card storage.
    pub fn linearize_temperature(
        &self, cx: &Cx<'_>, boundary: &ThermalBoundary, temperature: &[f64], max_face_pairs: usize,
    ) -> Result<ContactTemperatureLinearization, ConductionError> {
        self.linearize_temperature_with_policy(
            cx, boundary, temperature, max_face_pairs, DerivativePolicy::TwoSided,
        )
    }

    /// Build a primal Newton tangent for the same physical contact residual.
    ///
    /// Continuous piecewise-linear source curves use their right segment at
    /// interior knots and their inward segment at either tabulated endpoint.
    /// Scalar laws have zero temperature partials. These generalized tangents
    /// let a physical step start at an equilibrium source corner or move inward
    /// from a source-domain boundary without smoothing or extrapolating R(T).
    /// Every actual state still passes its original pinned source query; an
    /// out-of-domain Newton trial must be rejected by the residual owner.
    /// Two-dimensional surfaces retain the strict unique-gradient requirement.
    ///
    /// This method is for forward nonlinear solves only. Accepted-endpoint
    /// sensitivities and adjoints must use [`Self::linearize_temperature`],
    /// which still requires a unique two-sided derivative and neighbourhood.
    pub fn linearize_temperature_for_newton(
        &self, cx: &Cx<'_>, boundary: &ThermalBoundary, temperature: &[f64], max_face_pairs: usize,
    ) -> Result<ContactTemperatureLinearization, ConductionError> {
        self.linearize_temperature_with_policy(
            cx, boundary, temperature, max_face_pairs, DerivativePolicy::PrimalSemismooth,
        )
    }

    fn linearize_temperature_with_policy(
        &self, cx: &Cx<'_>, boundary: &ThermalBoundary, temperature: &[f64],
        max_face_pairs: usize, policy: DerivativePolicy,
    ) -> Result<ContactTemperatureLinearization, ConductionError> {
        self.linearize_state_with_policy(cx, boundary, temperature, None, None, None, max_face_pairs, policy)
    }

    /// Bind the complete contact h derivative to the original native charts.
    /// Strict derivatives require unique two-sided chart and source branches.
    pub fn linearize_enthalpy(&self, cx: &Cx<'_>, boundary: &ThermalBoundary,
        engine: &crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        h: &[f64], max_face_pairs: usize)
        -> Result<ContactTemperatureLinearization, crate::transient::enthalpy::EnthalpyError>
    {
        let (temperature, slopes) = self.enthalpy_fields(cx, engine, h, true)?;
        Ok(self.linearize_enthalpy_at_fields(cx, boundary, h, &temperature, &slopes,
            max_face_pairs, |vertex| engine.phase_curve_at(vertex))?)
    }

    // Shared strict binder for native BE endpoints and genuine RK stages.
    // The caller owns the original native fields; source chart identity and
    // two-sided constitutive admission remain enforced here.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn linearize_enthalpy_at_fields<'c>(&self, cx: &Cx<'_>, boundary: &ThermalBoundary,
        h: &[f64], temperature: &[f64], slopes: &[f64], max_face_pairs: usize,
        mut curve_at: impl FnMut(usize) -> Option<&'c fs_material::phase::EquilibriumEnthalpyPhaseCurve>)
        -> Result<ContactTemperatureLinearization, ConductionError>
    {
        self.check_material_curves(cx, &mut curve_at)?;
        let phase = self.phase_field(cx, h, temperature, true, true)?;
        let response = self.linearize_state_with_policy(cx, boundary, temperature, phase.as_ref(),
            Some(h), Some(slopes), max_face_pairs, DerivativePolicy::TwoSided)?;
        // Fixed and temperature-only contacts also inherit T(h). The retained
        // frozen operator identifies their actual participating vertices,
        // including fixed nonmatching faces, without another geometry owner.
        for vertex in 0..h.len() {
            if vertex % 256 == 0 { poll(cx, vertex)?; }
            if !response.frozen.row(vertex).0.is_empty()
                && !self.phase_curves.get(vertex).is_some_and(|binding| binding.is_some()) {
                super::phase::admit_chart_derivative(
                    curve_at(vertex).ok_or_else(|| invalid("<contact-tangent>",
                        "original native chart missing at a participating vertex"))?,
                    h[vertex], false,
                )?;
            }
        }
        Ok(response)
    }

    pub(crate) fn linearize_enthalpy_for_newton(&self, cx: &Cx<'_>, boundary: &ThermalBoundary,
        h: &[f64], temperature: &[f64], slopes: &[f64], max_face_pairs: usize)
        -> Result<ContactTemperatureLinearization, ConductionError>
    {
        let phase = self.phase_field(cx, h, temperature, true, false)?;
        self.linearize_state_with_policy(cx, boundary, temperature, phase.as_ref(),
            Some(h), Some(slopes), max_face_pairs, DerivativePolicy::PrimalSemismooth)
    }

    #[allow(clippy::too_many_arguments)]
    fn linearize_state_with_policy(&self, cx: &Cx<'_>, boundary: &ThermalBoundary,
        temperature: &[f64], phase: Option<&super::phase::PhaseField>, h: Option<&[f64]>,
        h_slopes: Option<&[f64]>, max_face_pairs: usize, policy: DerivativePolicy)
        -> Result<ContactTemperatureLinearization, ConductionError>
    {
        poll(cx, 0)?;
        let base = self.variable_pairs.checked_add(self.nonmatching_triangles)
            .ok_or_else(|| invalid("<contact-tangent>", "contact triangle count overflows"))?;
        let count = self.fixed.iter().try_fold(base, |n, surface|
            n.checked_add(surface.face_pairs().len()))
            .ok_or_else(|| invalid("<contact-tangent>", "contact face allowance overflows"))?;
        let n = self.mesh.vertex_count();
        let entry_budget = if phase.is_some() { 864 } else { 432 };
        if h.is_some_and(|values| values.len() != n) || h_slopes.is_some_and(|values| values.len() != n) {
            return Err(invalid("<contact-tangent>", "enthalpy contact fields must match the original nodal layout"));
        }
        if count > max_face_pairs || count.checked_mul(entry_budget * 24)
            .and_then(|bytes| n.checked_mul(8 * 8).and_then(|nodes| bytes.checked_add(nodes)))
            .is_none_or(|bytes| bytes > isize::MAX as usize)
        {
            return Err(invalid("<contact-tangent>", "contact tangent exceeds its paired-face or addressable storage allowance"));
        }
        let physical = self.evaluate_state(cx, boundary, temperature, phase, &mut |_| {})?;
        let phase_coordinates: BTreeMap<_, _> = physical.material_trace_states().iter().map(|state|
            ((state.interface.as_str(), state.pair, state.quadrature_index), state.liquid_mass_fractions)).collect();
        let active_at_vertices = |law: &TemperatureContact<'_>, vertices: [usize; 2]| -> [bool; 4] {
            if law.material_states.is_none() { return [true, true, false, false]; }
            let phase = phase.expect("material field admitted");
            [phase.temperature_slopes[vertices[0]] != 0.0, phase.temperature_slopes[vertices[1]] != 0.0,
                phase.fraction_slopes[vertices[0]] != 0.0, phase.fraction_slopes[vertices[1]] != 0.0]
        };
        let active_at_pair = |law: &TemperatureContact<'_>, pair: super::InterfaceFacePair| -> [bool; 4] {
            let a = self.mesh.boundary()[pair.side_a].vertices;
            let b = self.mesh.boundary()[pair.side_b].vertices;
            let mut active = [false; 4];
            for i in 0..3 {
                for (result, value) in active.iter_mut().zip(active_at_vertices(law, [a[i] as usize, b[i] as usize])) { *result |= value; }
            }
            active
        };
        let mut responses = BTreeMap::new();
        for (i, law) in self.variable.iter().enumerate() {
            poll(cx, i)?;
            // Prepare each distinct pinned source chart once, even when many
            // mapped faces use different fixed contexts on that same chart.
            for (index, source) in std::iter::once(&law.resistance)
                .chain(law.face_references.iter().map(|entry| &entry.resistance)).enumerate() {
                poll(cx, index)?;
                let key = (law.name(), source.receipt().selected);
                if let std::collections::btree_map::Entry::Vacant(entry) = responses.entry(key) {
                    entry.insert((law, law.derivative_response(key.1)?));
                }
            }
        }
        for surface in &physical.interfaces.nonmatching {
            if let Some(law) = self.variable.iter().find(|law| law.name() == surface.name) {
                let (_, response) = responses.get(&(law.name(), law.resistance.receipt().selected))
                    .ok_or_else(|| invalid(law.name(), "contact support check lost its original source"))?;
                if matches!(policy, DerivativePolicy::TwoSided)
                    || matches!(response, ContactDerivative::Surface(_)) && law.material_states.is_none() {
                    let fractions = if law.material_states.is_some() {
                        Some(phase.expect("material field admitted").fractions.as_slice())
                    } else { None };
                    surface.matching_state_support(cx, temperature, fractions, |vertices, traces, fractions|
                        law.state_neighbourhood(cx, None, traces, fractions, active_at_vertices(law, vertices)))?;
                }
            }
        }
        let mut thermal_partials = BTreeMap::new(); let mut phase_partials = BTreeMap::new();
        for (i, state) in physical.states().iter().enumerate() {
            poll(cx, i)?;
            let (law, response) = responses.get(&(state.interface.as_str(), state.resistance.receipt().selected))
                .ok_or_else(|| invalid("<contact-tangent>", "local contact state lost its source owner"))?;
            let fractions = phase_coordinates.get(&(state.interface.as_str(), state.pair, None)).copied();
            let d = law.state_partials(cx, Some(state.pair.side_a),
                [state.side_a_temperature_k, state.side_b_temperature_k], fractions,
                active_at_pair(law, state.pair), response, policy)?;
            if thermal_partials.insert(state.pair, [d[0], d[1]]).is_some() {
                return Err(invalid(law.name(), "duplicate local contact derivative"));
            }
            phase_partials.insert(state.pair, [d[2], d[3]]);
        }
        let matching = physical.interfaces.matching_integration();
        let matching_thermal = matching.temperature_feedback(cx, temperature, &thermal_partials, max_face_pairs)?;
        let matching_phase = if phase.is_some() {
            Some(matching.temperature_feedback(cx, temperature, &phase_partials, max_face_pairs)?)
        } else { None };
        let mut local = BTreeMap::<&str, (Vec<[f64; 2]>, Vec<[f64; 2]>)>::new();
        for (i, state) in physical.quadrature_states().iter().enumerate() {
            poll(cx, i)?;
            let (law, response) = responses.get(&(state.interface.as_str(), state.resistance.receipt().selected))
                .ok_or_else(|| invalid("<contact-tangent>", "quadrature state lost its source owner"))?;
            let entries = local.entry(state.interface.as_str()).or_default();
            if state.quadrature_index != entries.0.len() {
                return Err(invalid(law.name(), "quadrature derivative order differs from the physical rule"));
            }
            let fractions = phase_coordinates.get(&(state.interface.as_str(), state.pair, Some(state.quadrature_index))).copied();
            let d = law.state_partials(cx, Some(state.pair.side_a),
                [state.side_a_temperature_k, state.side_b_temperature_k], fractions,
                active_at_pair(law, state.pair), response, policy)?;
            entries.0.push([d[0], d[1]]); entries.1.push([d[2], d[3]]);
        }
        let mut thermal_coo = Coo::new(n, n); let mut phase_coo = Coo::new(n, n);
        for row in 0..n {
            if row % 256 == 0 { poll(cx, row)?; }
            let (columns, values) = matching_thermal.row(row);
            for (&column, &value) in columns.iter().zip(values) { thermal_coo.push(row, column, value); }
            if let Some(matrix) = &matching_phase {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) { phase_coo.push(row, column, value); }
            }
        }
        for surface in &physical.interfaces.nonmatching {
            if let Some((thermal, fraction)) = local.get(surface.name.as_str()) {
                surface.temperature_feedback(cx, temperature, thermal, &mut thermal_coo)?;
                if phase.is_some() { surface.temperature_feedback(cx, temperature, fraction, &mut phase_coo)?; }
            }
        }
        let feedback = thermal_coo.assemble();
        let phase_feedback = if let Some(phase) = phase {
            let fraction = phase_coo.assemble(); let mut enthalpy = Coo::new(n, n);
            for row in 0..n {
                if row % 256 == 0 { poll(cx, row)?; }
                let (columns, values) = fraction.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    let slope = phase.fraction_slopes[column];
                    let derivative = finite(value*slope)?;
                    if value != 0.0 && slope != 0.0 && derivative == 0.0 {
                        return Err(invalid("<contact-tangent>", "nonzero contact phase derivative underflows"));
                    }
                    enthalpy.push(row, column, derivative);
                }
            }
            Some(enthalpy.assemble())
        } else { None };
        let phase_feedback_transpose = phase_feedback.as_ref().map(ops::transpose);
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
            specific_enthalpy: h.map(<[f64]>::to_vec), enthalpy_slopes: h_slopes.map(<[f64]>::to_vec),
            phase_feedback, phase_feedback_transpose,
        })
    }
}
