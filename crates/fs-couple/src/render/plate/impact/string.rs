//! Fixed-end filament stretching in the existing nonlinear impact composition.
//!
//! fs-nlmodal owns the Kirchhoff--Carrier stress channel and prestressed-beam
//! frequencies. This adapter preserves the wire's existing mass-normalized sine
//! coordinates, bending, contacts and damping. It supplies no new integrator.
//! A coiled wire needs its effective axial rigidity supplied independently;
//! neither winding mass nor transverse bending rigidity identifies E*A.
use std::rc::Rc;
use fs_math::det;
use fs_nlmodal::{KcStringParams, SosModalStorage, kirchhoff_carrier_string};
use fs_phs::Storage;
use super::{BodyPotential, ImpactBody, ImpactError, MAX_IMPACT_MODES, invalid};
use super::linear::wire::WireSpan;
use crate::modal_acoustic_time::ModalAcousticState;

/// An explicit constitutive input and moderate-slope validity bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringStretching {
    /// Effective axial rigidity E*A [N], not tension or spring rate E*A/L.
    /// Zero retains the linear potential for controlled comparisons.
    pub axial_rigidity_n: f64,
    /// Upper bound on |dw/dx|; finite in (0,0.3]. No clipping is performed.
    pub maximum_slope: f64,
}
impl StringStretching {
    pub fn validate(self) -> Result<(), ImpactError> {
        if !self.axial_rigidity_n.is_finite() || self.axial_rigidity_n < 0.0
            || !self.maximum_slope.is_finite() || self.maximum_slope <= 0.0
            || self.maximum_slope > 0.3 {
            return Err(invalid("string stretching needs finite nonnegative axial rigidity and a slope limit in (0,0.3]"));
        }
        Ok(())
    }
}

/// Read-only physical diagnostics. Stretching storage is ALREADY in total H.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StringObservation {
    /// Continuous-span bound sum_n |Q_n|*n*pi/L, not a sampled maximum.
    pub slope_bound: f64,
    /// Averaged geometric extension / reference length (dimensionless).
    pub additional_strain: f64,
    /// Initial tension plus E*A times additional_strain [N].
    pub tension_n: f64,
    pub stretching_energy_j: f64,
}

/// Immutable, shared storage compiled by the existing nonlinear-string owner.
#[derive(Clone)]
pub struct StringPotential {
    storage: Rc<SosModalStorage>,
    stretching: StringStretching,
    rest_tension_n: f64,
}
impl core::fmt::Debug for StringPotential {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("StringPotential").field("omegas", &self.storage.omegas)
            .field("stretching", &self.stretching).field("rest_tension_n", &self.rest_tension_n).finish()
    }
}
impl WireSpan {
    /// Keep all original modes, initial motion, damping and contact projections;
    /// add the existing averaged-tension potential using supplied axial rigidity.
    /// The initial continuous-span slope bound must pass before publication.
    pub fn stretching_body(&self, initial: Vec<ModalAcousticState>, stretching: StringStretching)
        -> Result<ImpactBody, ImpactError>
    {
        stretching.validate()?;
        if initial.is_empty() || initial.len() > MAX_IMPACT_MODES {
            return Err(invalid("stretching wire exceeds the nonlinear mechanical mode budget"));
        }
        let mut body = self.body(initial)?; // original geometry/linear-law admission
        let n = body.initial.len();
        let mut storage = kirchhoff_carrier_string(&KcStringParams {
            length: self.length_m(), tension: self.tension_n,
            lin_density: self.linear_density_kg_m, ea: stretching.axial_rigidity_n,
        }, n).map_err(|e| ImpactError::Owner(e.to_string()))?;
        let BodyPotential::Linear(omegas) = body.potential else {
            return Err(invalid("wire body changed its expected linear preparation"));
        };
        // The stress channel adds geometric extension, not a second copy of T0.
        // Preserve EI and the exact original prestressed-beam frequencies.
        storage.omegas = omegas;
        if storage.channels.len() != 1 || !storage.channels[0].coefficient.is_finite()
            || storage.channels[0].coefficient < 0.0 || storage.channels[0].coupling.len() != n*n
            || storage.channels[0].coupling.iter().any(|v| !v.is_finite())
            || (0..n).any(|i| storage.channels[0].coupling[i*n+i] <= 0.0) {
            return Err(invalid("derived string stretching channel is unrepresentable"));
        }
        let potential = StringPotential { storage: Rc::new(storage), stretching,
            rest_tension_n: self.tension_n };
        let mut q = [0.0; MAX_IMPACT_MODES];
        for (i, s) in body.initial.iter().enumerate() { q[i] = s.displacement_m_sqrt_kg; }
        potential.observe(&q[..n])?;
        body.potential = BodyPotential::String(potential);
        Ok(body)
    }
}
impl StringPotential {
    pub fn mode_count(&self) -> usize { self.storage.omegas.len() }
    pub fn omegas(&self) -> &[f64] { &self.storage.omegas }
    pub fn slope_limit(&self) -> f64 { self.stretching.maximum_slope }

    fn valid_coordinates(&self, q: &[f64]) -> bool {
        q.len() == self.mode_count() && q.iter().all(|v| v.is_finite())
    }
    /// Complete linear plus stretching potential. Momentum is owned by the host.
    pub fn potential(&self, q: &[f64]) -> f64 {
        if !self.valid_coordinates(q) { return f64::NAN; }
        let mut x = [0.0; 2*MAX_IMPACT_MODES];
        for (i, q) in q.iter().enumerate() { x[2*i] = *q; }
        self.storage.hamiltonian(&x[..2*q.len()])
    }
    /// Exact gradient of the same fs-nlmodal storage, without allocation.
    pub fn gradient(&self, q: &[f64], out: &mut [f64]) {
        if !self.valid_coordinates(q) || out.len() != q.len() { out.fill(f64::NAN); return; }
        let (mut x, mut g) = ([0.0; 2*MAX_IMPACT_MODES], [0.0; 2*MAX_IMPACT_MODES]);
        for (i, q) in q.iter().enumerate() { x[2*i] = *q; }
        self.storage.gradient(&x[..2*q.len()], &mut g[..2*q.len()]);
        for (i, v) in out.iter_mut().enumerate() { *v = g[2*i]; }
    }
    /// Exact derivative of the owner's admitted diagonal stress channel.
    /// Coefficients come from that SAME storage, never a fitted Duffing constant.
    pub fn hessian_vector(&self, q: &[f64], d: &[f64], out: &mut [f64]) {
        let n = self.mode_count();
        if !self.valid_coordinates(q) || d.len() != n || out.len() != n
            || d.iter().any(|v| !v.is_finite()) { out.fill(f64::NAN); return; }
        let ch = &self.storage.channels[0];
        let mut s = 0.0; let mut cross = 0.0;
        for i in 0..n {
            let e = ch.coupling[i*n+i]; s += e*q[i]*q[i]; cross += e*q[i]*d[i];
        }
        for i in 0..n {
            out[i] = self.omegas()[i].powi(2)*d[i]
                + ch.coefficient*ch.coupling[i*n+i]*(s*d[i] + 2.0*cross*q[i]);
        }
    }
    /// Observe and enforce the admitted continuous-span slope bound. It can
    /// conservatively refuse cancelling modal shapes; it never misses a peak
    /// between quadrature stations. This is a validity limit, not error control.
    pub fn observe(&self, q: &[f64]) -> Result<StringObservation, ImpactError> {
        if !self.valid_coordinates(q) { return Err(invalid("string observation needs every finite modal coordinate")); }
        let n = self.mode_count(); let ch = &self.storage.channels[0];
        let mut s = 0.0; let mut slope_bound = 0.0;
        for i in 0..n {
            let e = ch.coupling[i*n+i]; s += e*q[i]*q[i]; slope_bound += det::sqrt(e)*q[i].abs();
        }
        let additional_strain = 0.25*s;
        let tension_n = self.rest_tension_n + self.stretching.axial_rigidity_n*additional_strain;
        let stretching_energy_j = 0.25*ch.coefficient*s*s;
        if ![slope_bound,additional_strain,tension_n,stretching_energy_j].iter().all(|x| x.is_finite())
            || slope_bound > self.stretching.maximum_slope {
            return Err(invalid("string trial exceeds finite/moderate-slope validity"));
        }
        Ok(StringObservation { slope_bound, additional_strain, tension_n, stretching_energy_j })
    }
    pub(super) fn observe_interleaved(&self, state: &[f64], first: usize) -> Result<StringObservation, ImpactError> {
        let end = first.checked_add(self.mode_count()).and_then(|n| n.checked_mul(2))
            .ok_or_else(|| invalid("string state address overflow"))?;
        if state.len() < end { return Err(invalid("string observation needs the complete mechanical state")); }
        let mut q = [0.0; MAX_IMPACT_MODES];
        for i in 0..self.mode_count() { q[i] = state[2*(first+i)]; }
        self.observe(&q[..self.mode_count()])
    }
}

#[cfg(test)]
#[path = "string_tests.rs"]
mod tests;
