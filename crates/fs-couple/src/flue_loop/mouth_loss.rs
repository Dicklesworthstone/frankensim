//! Optional separated-flow pressure drop inside the existing mouth/bore solve.
//! The fixed-window Bernoulli contraction owner supplies the constitutive law;
//! the original jet, propagation, radiation and time-discretization remain.

use super::{FlueError, FlueVoice, MouthStep};
use crate::bernoulli_aperture::contraction::{ContractionLoss, ContractionSpec};
use crate::wind_line::WindBoundary;

impl FlueVoice {
    /// Attach a quasi-steady separation loss to this voice's ACTUAL mouth window
    /// and gas, before any sample. It acts on acoustic mouth flow, not jet flow.
    /// The coefficient/source are supplied data, not inferred from a flute name.
    ///
    /// Compatible with the original TMM and material/regional/fractional tube
    /// images. Separate linear radiation resistance remains unchanged. All
    /// performers that own this voice automatically advance the same loss.
    /// No claim of resolved shedding, turbulent noise, high-Mach flow, measured
    /// flute timbre, or transient validity follows from this reduced closure.
    ///
    /// # Errors
    /// Repeated/late attachment, poisoned voice, or original loss admission.
    pub fn with_mouth_contraction(mut self, spec: ContractionSpec) -> Result<Self, FlueError> {
        if self.poisoned || self.samples != 0 || self.mouth_contraction.is_some() {
            return Err(FlueError::Invalid { what: "mouth contraction requires one cold attachment" });
        }
        self.mouth_contraction = Some(ContractionLoss::new(self.geometry.flue_width_m,
            self.geometry.cut_up_m, self.gas, spec).map_err(FlueError::Physics)?);
        Ok(self)
    }
    /// Original supplied loss and its actual gas/geometry; None is the unchanged
    /// linear mouth image. Unity area ratio does NOT mean zero separation loss.
    #[must_use]
    pub fn mouth_contraction(&self) -> Option<&ContractionLoss> { self.mouth_contraction.as_ref() }
    /// Actual accepted acoustic mouth flow [m^3/s], not injected jet volume flow.
    #[must_use]
    pub const fn mouth_flow_m3_s(&self) -> f64 { self.q_prev }
    /// Complete fixed material gas used by the original jet and mouth.
    #[must_use]
    pub const fn gas(&self) -> &fs_material::gas::GasState { &self.gas }
}

// Eliminate the linear characteristic port, not the nonlinear constitutive law.
// B>0: (L/dt+R+Z(1+r)/B) Q + dp_loss(Q) = source+L/dt*Qold-2h/B.
// B=0: the wave boundary itself constrains Q=-h/Z; pressure remains solvable.
#[allow(clippy::too_many_arguments)]
pub(super) fn solve(boundary: WindBoundary, source: f64, previous: f64,
    inertance: f64, resistance: f64, zc: f64, dt: f64, loss: &ContractionLoss)
    -> Result<(MouthStep, f64), FlueError>
{
    let bad = || FlueError::Invalid { what: "nonlinear mouth needs a finite passive same-sample response" };
    let r = boundary.direct_reflection;
    let history = boundary.history_pressure_pa;
    if ![r,history,source,previous,inertance,resistance,zc,dt].iter().all(|v|v.is_finite())
        || r.abs()>1.0 || inertance<0.0 || resistance<0.0 || zc<=0.0 || dt<=0.0 {
        return Err(bad());
    }
    let alpha = inertance/dt;
    let series = alpha+resistance;
    let drive = source+alpha*previous;
    if ![alpha,series,drive].iter().all(|x|x.is_finite()) { return Err(bad()); }
    let frame = if r == 1.0 {
        loss.at_flow(-history/zc)
    } else {
        loss.solve(series+zc*((1.0+r)/(1.0-r)), drive-2.0*(history/(1.0-r)))
    }.map_err(FlueError::Physics)?;
    let q = frame.volume_flow_m3_s;
    // Choose a reconstruction whose boundary denominator is >=1. Keep the
    // incoming arithmetic identical to WindLineBank's actual commitment.
    let outgoing = if r <= 0.0 {
        (zc*q+history)/(1.0-r)
    } else {
        (drive-series*q-frame.pressure_drop_pa-history)/(1.0+r)
    };
    let incoming = r*outgoing+history;
    let pressure = outgoing+incoming;
    let momentum = source-alpha*(q-previous)-resistance*q-frame.pressure_drop_pa-pressure;
    let continuity = outgoing-incoming-zc*q;
    let momentum_scale = source.abs().max((alpha*q).abs()).max((alpha*previous).abs())
        .max((resistance*q).abs()).max(frame.pressure_drop_pa.abs()).max(pressure.abs());
    let flow_scale = outgoing.abs().max(incoming.abs()).max((zc*q).abs());
    if ![q,outgoing,incoming,pressure,momentum,continuity,momentum_scale,flow_scale]
        .iter().all(|x|x.is_finite())
        || momentum.abs()>2048.0*f64::EPSILON*momentum_scale
        || continuity.abs()>2048.0*f64::EPSILON*flow_scale {
        return Err(bad());
    }
    Ok((MouthStep { flow_m3_s:q, outgoing_pa:outgoing, incoming_pa:incoming,
        pressure_pa:pressure },frame.dissipated_power_w))
}

#[cfg(test)]
mod tests;
