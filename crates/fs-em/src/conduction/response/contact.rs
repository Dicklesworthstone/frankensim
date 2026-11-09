//! Source closure at original interface coordinates, with optional bounded quadrature.
//! Uniform-trace admission remains the default. Explicit quadrature closes each
//! local electrical state on the same points used by native current and heat.

use crate::conduction::{ConductionError, admit_contact_source};
use fs_material::electrical_interface::{
    ResolvedElectricalInterface, response::ScalarElectricalInterfaceResponse,
};
use fs_material::state_point::MaterialStatePointError;
use std::sync::Arc;

/// Separate dimensional contact-law gates. Voltage and density tolerances also
/// bound variation across electrically varying traces; static traces retain
/// the native nonuniform P1 response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactResponseControls {
    pub absolute_resistance_tolerance_ohm_m2: f64,
    pub relative_resistance_tolerance: f64,
    pub voltage_tolerance_v: f64,
    pub current_density_tolerance_a_m2: f64,
}
impl ContactResponseControls {
    pub fn validate(self) -> Result<(), ContactResponseError> {
        if !self.absolute_resistance_tolerance_ohm_m2.is_finite()
            || self.absolute_resistance_tolerance_ohm_m2 < 0.0
            || !self.relative_resistance_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_resistance_tolerance)
            || (self.absolute_resistance_tolerance_ohm_m2 == 0.0
                && self.relative_resistance_tolerance == 0.0)
            || ![
                self.voltage_tolerance_v,
                self.current_density_tolerance_a_m2,
            ]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(ContactResponseError::Invalid(
                "invalid dimensional contact response gates",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum ContactResponseError {
    Cancelled,
    Material(MaterialStatePointError),
    Invalid(&'static str),
    NonuniformTrace {
        contact: usize,
        voltage_spread_v: f64,
        current_density_spread_a_m2: f64,
    },
}
impl From<MaterialStatePointError> for ContactResponseError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}

/// Explicit composite-Simpson resolution for nonuniform electrical contact laws.
/// Each edge uses `2*subdivisions+1` points, including both original vertices.
/// Triangular faces use `3*subdivisions²+3*subdivisions+1` shared points.
/// The point budget covers all original contacts, independently of iteration
/// and field-solve bounds. Resolution is caller-owned, not an error certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ContactResponseQuadrature {
    pub subdivisions: usize,
    pub maximum_points: usize,
}
impl ContactResponseQuadrature {
    pub(crate) fn edge_points(self, contacts: usize) -> Result<usize, ContactResponseError> {
        self.admit_points(
            contacts,
            self.subdivisions
                .checked_mul(2)
                .and_then(|n| n.checked_add(1)),
        )
    }
    pub(crate) fn triangle_points(self, contacts: usize) -> Result<usize, ContactResponseError> {
        let n = self.subdivisions;
        self.admit_points(
            contacts,
            n.checked_mul(n)
                .and_then(|v| v.checked_mul(3))
                .and_then(|v| v.checked_add(n.checked_mul(3)?))
                .and_then(|v| v.checked_add(1)),
        )
    }
    fn admit_points(
        self,
        contacts: usize,
        points: Option<usize>,
    ) -> Result<usize, ContactResponseError> {
        if self.subdivisions == 0
            || contacts == 0
            || points
                .and_then(|n| n.checked_mul(contacts))
                .is_none_or(|n| n > self.maximum_points)
        {
            return Err(ContactResponseError::Invalid(
                "positive bounded contact quadrature resolution and complete point budget required",
            ));
        }
        Ok(points.unwrap())
    }
}

/// Complete source bindings in the native owner's canonical contact order.
#[derive(Clone, Debug)]
pub struct ContactResponseSources {
    laws: Arc<[ScalarElectricalInterfaceResponse]>,
    controls: ContactResponseControls,
    quadrature: Option<ContactResponseQuadrature>,
}
impl ContactResponseSources {
    pub fn new(
        laws: Vec<ScalarElectricalInterfaceResponse>,
        controls: ContactResponseControls,
    ) -> Result<Self, ContactResponseError> {
        controls.validate()?;
        Ok(Self {
            laws: laws.into(),
            controls,
            quadrature: None,
        })
    }
    pub fn laws(&self) -> &[ScalarElectricalInterfaceResponse] {
        &self.laws
    }
    pub fn controls(&self) -> ContactResponseControls {
        self.controls
    }
    pub fn quadrature(&self) -> Option<ContactResponseQuadrature> {
        self.quadrature
    }
    pub(crate) fn with_quadrature(
        mut self,
        quadrature: ContactResponseQuadrature,
    ) -> Result<Self, ContactResponseError> {
        quadrature.edge_points(self.laws.len())?;
        self.quadrature = Some(quadrature);
        Ok(self)
    }
}

pub(crate) fn admit_source(
    index: usize,
    source: &ResolvedElectricalInterface,
    axes: &[Option<String>; 2],
    bindings: Option<&ContactResponseSources>,
) -> Result<f64, ConductionError> {
    let Some(bindings) = bindings else {
        return admit_contact_source(source, axes);
    };
    let law = bindings.laws().get(index).ok_or(ConductionError::Invalid(
        "missing original contact response source",
    ))?;
    if axes != law.liquid_mass_fraction_axes() {
        return Err(ConductionError::Invalid(
            "contact response changed original phase ownership",
        ));
    }
    law.admit(source).map_err(|_| {
        ConductionError::Invalid("contact response changed source identity or domain")
    })?;
    Ok(source.resistance_ohm_m2())
}

/// Original source evaluated at every vertex or admitted quadrature point, with independent local
/// resistance, jump and current defects. These residuals are never heat credits.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactResponseSample {
    pub contact: usize,
    /// Complete source states that supplied the native operator, in the same
    /// vertex/quadrature ordering as `evaluated`.
    pub applied: Vec<ResolvedElectricalInterface>,
    pub evaluated: Vec<ResolvedElectricalInterface>,
    pub maximum_resistance_defect_ohm_m2: f64,
    pub maximum_voltage_defect_v: f64,
    pub maximum_current_density_defect_a_m2: f64,
    /// Conservative coefficient-work defect summed over admitted quadrature,
    /// or the original uniform-trace bound. This is not physical heat.
    pub absolute_power_defect_w: f64,
}

pub(crate) fn evaluate<const N: usize>(
    contact: usize,
    law: &ScalarElectricalInterfaceResponse,
    seed: &ResolvedElectricalInterface,
    applied: &ResolvedElectricalInterface,
    voltage: [f64; N],
    density: [f64; N],
    native_power_w: f64,
    controls: ContactResponseControls,
) -> Result<(ContactResponseSample, bool), ContactResponseError> {
    if N == 0 || !native_power_w.is_finite() || native_power_w < 0.0 {
        return Err(ContactResponseError::Invalid(
            "finite passive native contact response required",
        ));
    }
    law.admit(seed)?;
    law.admit(applied)?;
    let mut voltage_range = [voltage[0]; 2];
    let mut density_range = [density[0]; 2];
    for (&v, &j) in voltage.iter().zip(&density) {
        if !v.is_finite() || !j.is_finite() {
            return Err(ContactResponseError::Invalid(
                "finite local contact fields required",
            ));
        }
        voltage_range = [voltage_range[0].min(v), voltage_range[1].max(v)];
        density_range = [density_range[0].min(j), density_range[1].max(j)];
    }
    let voltage_spread = voltage_range[1] - voltage_range[0];
    let density_spread = density_range[1] - density_range[0];
    if !law.axes().is_empty()
        && (voltage_spread > controls.voltage_tolerance_v
            || density_spread > controls.current_density_tolerance_a_m2)
    {
        return Err(ContactResponseError::NonuniformTrace {
            contact,
            voltage_spread_v: voltage_spread,
            current_density_spread_a_m2: density_spread,
        });
    }
    let mut sample = ContactResponseSample {
        contact,
        applied: Vec::with_capacity(N),
        evaluated: Vec::with_capacity(N),
        maximum_resistance_defect_ohm_m2: 0.0,
        maximum_voltage_defect_v: 0.0,
        maximum_current_density_defect_a_m2: 0.0,
        absolute_power_defect_w: 0.0,
    };
    let used = applied.resistance_ohm_m2();
    let mut converged = true;
    for (&v, &j) in voltage.iter().zip(&density) {
        let evaluated = law.at_operating_point(seed, v, j)?;
        let resistance = evaluated.resistance_ohm_m2();
        let dr = (resistance - used).abs();
        let dv = (v - resistance * j).abs();
        let dj = (j - v / resistance).abs();
        let power = native_power_w * (dr / used);
        let tolerance = controls.absolute_resistance_tolerance_ohm_m2
            + controls.relative_resistance_tolerance * resistance.max(used);
        if ![dr, dv, dj, power, tolerance].iter().all(|v| v.is_finite()) {
            return Err(ContactResponseError::Invalid(
                "unrepresentable local contact constitutive residual",
            ));
        }
        converged &= dr <= tolerance
            && dv <= controls.voltage_tolerance_v
            && dj <= controls.current_density_tolerance_a_m2;
        sample.maximum_resistance_defect_ohm_m2 = sample.maximum_resistance_defect_ohm_m2.max(dr);
        sample.maximum_voltage_defect_v = sample.maximum_voltage_defect_v.max(dv);
        sample.maximum_current_density_defect_a_m2 =
            sample.maximum_current_density_defect_a_m2.max(dj);
        sample.absolute_power_defect_w = sample.absolute_power_defect_w.max(power);
        sample.applied.push(applied.clone());
        sample.evaluated.push(evaluated);
    }
    Ok((sample, converged))
}

/// Close each sampled source against its own native voltage/current state.
/// The applied resistance and positive quadrature power come from the SAME
/// native operator; numerical source defects never enter physical heat.
#[allow(clippy::too_many_arguments)]
pub(crate) fn evaluate_sampled<const N: usize>(
    contact: usize,
    law: &ScalarElectricalInterfaceResponse,
    seed: &ResolvedElectricalInterface,
    applied: &[ResolvedElectricalInterface],
    observations: &[crate::conduction::ConductionContactQuadraturePoint<N>],
    controls: ContactResponseControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ContactResponseSample, bool), ContactResponseError> {
    if applied.is_empty() || applied.len() != observations.len() {
        return Err(ContactResponseError::Invalid(
            "complete applied contact quadrature states required",
        ));
    }
    law.admit(seed)?;
    let mut sample = ContactResponseSample {
        contact,
        applied: Vec::with_capacity(applied.len()),
        evaluated: Vec::with_capacity(applied.len()),
        maximum_resistance_defect_ohm_m2: 0.0,
        maximum_voltage_defect_v: 0.0,
        maximum_current_density_defect_a_m2: 0.0,
        absolute_power_defect_w: 0.0,
    };
    let mut converged = true;
    for (applied, point) in applied.iter().zip(observations) {
        if cancelled() {
            return Err(ContactResponseError::Cancelled);
        }
        law.admit(applied)?;
        let used = applied.resistance_ohm_m2();
        let v = point.voltage_jump_v;
        let j = point.current_density_a_m2;
        if point.resistance_ohm_m2 != used
            || ![v, j, point.power_w].iter().all(|value| value.is_finite())
            || point.power_w < 0.0
        {
            return Err(ContactResponseError::Invalid(
                "native contact quadrature lost its applied source state",
            ));
        }
        let evaluated = law.at_operating_point(seed, v, j)?;
        let resistance = evaluated.resistance_ohm_m2();
        let dr = (resistance - used).abs();
        let dv = (v - resistance * j).abs();
        let dj = (j - v / resistance).abs();
        let power = point.power_w * (dr / used);
        let tolerance = controls.absolute_resistance_tolerance_ohm_m2
            + controls.relative_resistance_tolerance * resistance.max(used);
        if ![dr, dv, dj, power, tolerance]
            .iter()
            .all(|value| value.is_finite())
        {
            return Err(ContactResponseError::Invalid(
                "unrepresentable sampled contact constitutive residual",
            ));
        }
        converged &= dr <= tolerance
            && dv <= controls.voltage_tolerance_v
            && dj <= controls.current_density_tolerance_a_m2;
        sample.maximum_resistance_defect_ohm_m2 = sample.maximum_resistance_defect_ohm_m2.max(dr);
        sample.maximum_voltage_defect_v = sample.maximum_voltage_defect_v.max(dv);
        sample.maximum_current_density_defect_a_m2 =
            sample.maximum_current_density_defect_a_m2.max(dj);
        sample.absolute_power_defect_w += power;
        if !sample.absolute_power_defect_w.is_finite() {
            return Err(ContactResponseError::Invalid(
                "unrepresentable sampled contact power defect",
            ));
        }
        sample.applied.push(applied.clone());
        sample.evaluated.push(evaluated);
    }
    Ok((sample, converged))
}
