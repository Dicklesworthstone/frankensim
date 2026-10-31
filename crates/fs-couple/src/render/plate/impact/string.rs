//! Fixed-end filament stretching in the existing nonlinear impact composition.
//!
//! fs-nlmodal owns the Kirchhoff--Carrier stress channel and prestressed-beam
//! frequencies. This adapter preserves the wire's existing mass-normalized sine
//! coordinates, bending, contacts and damping. It supplies no new integrator.
//! A coiled wire needs its effective axial rigidity supplied independently;
//! neither winding mass nor transverse bending rigidity identifies E*A.
use std::sync::Arc;
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
    /// Continuous-span bound sum_n ||Q_n||*n*pi/L, not a sampled maximum.
    /// Q_n is a scalar for one plane or the two-component transverse amplitude.
    pub slope_bound: f64,
    /// Mean additional axial strain (dimensionless). Moving axial endpoints
    /// contribute their actual extension when longitudinal motion is retained.
    pub additional_strain: f64,
    /// Initial tension plus equilibrium E*A times additional_strain [N]. A
    /// material-bound string bank also includes its accepted hereditary stress.
    pub tension_n: f64,
    /// Axial geometric storage [J], plus accepted axial material storage when
    /// observed through a string bank carrying the original material history.
    pub stretching_energy_j: f64,
}

/// Immutable, shared storage compiled by the existing nonlinear-string owner.
#[derive(Clone)]
pub struct StringPotential {
    storage: Arc<SosModalStorage>,
    stretching: StringStretching,
    rest_tension_n: f64,
    // One physical span, with one or two mass-normalized sine-coordinate blocks.
    // Both transverse components contribute to the SAME stress channel.
    transverse_planes: usize,
    // Private potential-only lift coordinate in metres, never a new unit-mass
    // oscillator. SupportedStrings applies its physical kinetic transformation.
    affine_endpoint: bool,
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
        let potential = StringPotential { storage: Arc::new(storage), stretching,
            rest_tension_n: self.tension_n, transverse_planes: 1, affine_endpoint: false };
        let mut q = [0.0; MAX_IMPACT_MODES];
        for (i, s) in body.initial.iter().enumerate() { q[i] = s.displacement_m_sqrt_kg; }
        potential.observe(&q[..n])?;
        body.potential = BodyPotential::String(potential);
        Ok(body)
    }

    /// Two transverse components of the SAME isotropic span, ordered as the
    /// complete first sine block followed by the complete second sine block.
    ///
    /// The existing fs-nlmodal channel is lifted from E to diag(E,E), leaving
    /// ONE coefficient and ONE squared geometric extension. Its cross term
    /// transfers stretching work between the planes; two independent scalar
    /// stretching bodies would omit that term. The original mass, EI, tension
    /// and damping apply in both directions without detuning or rescaling.
    pub fn stretching_body_two_planes(
        &self, initial: Vec<ModalAcousticState>, stretching: StringStretching,
    ) -> Result<ImpactBody, ImpactError> {
        let n = self.damping_per_s.len();
        let total = n.checked_mul(2)
            .ok_or_else(|| invalid("two-plane string coordinate count overflow"))?;
        if n == 0 || total > MAX_IMPACT_MODES || initial.len() != total
            || initial.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite()
                || !s.velocity_m_sqrt_kg_per_s.is_finite()) {
            return Err(invalid("two-plane stretching needs two complete finite transverse blocks within the mode budget"));
        }
        let mut body = self.stretching_body(initial[..n].to_vec(), stretching)?;
        let BodyPotential::String(original) = body.potential else {
            return Err(invalid("stretching wire changed its expected nonlinear preparation"));
        };
        let mut storage = SosModalStorage {
            omegas: original.storage.omegas.clone(),
            channels: original.storage.channels.clone(),
        };
        let mut coupling = vec![0.0; total * total];
        for plane in 0..2 {
            for i in 0..n {
                let start = (plane * n + i) * total + plane * n;
                coupling[start..start+n]
                    .copy_from_slice(&storage.channels[0].coupling[i*n..(i+1)*n]);
            }
        }
        storage.channels[0].coupling = coupling;
        storage.omegas.extend_from_within(..);
        let potential = StringPotential {
            storage: Arc::new(storage), transverse_planes: 2, ..original
        };
        let mut q = [0.0; MAX_IMPACT_MODES];
        for (i, state) in initial.iter().enumerate() {
            q[i] = state.displacement_m_sqrt_kg;
        }
        potential.observe(&q[..total])?;
        body.potential = BodyPotential::String(potential);
        body.initial = initial;
        body.damping_per_s.extend_from_within(..);
        Ok(body)
    }
}
impl StringPotential {
    /// The SAME source-owned stress channels without quadratic bending/tension.
    /// Used only when the original exact linear owner already carries that
    /// storage and inertia. Zeroing the separate frequency list avoids both
    /// double stiffness and cancellation from subtracting two large energies.
    pub(crate) fn geometric_remainder(&self) -> Self {
        let storage = SosModalStorage {
            omegas: vec![0.0; self.mode_count()],
            channels: self.storage.channels.clone(),
        };
        Self { storage: Arc::new(storage), ..self.clone() }
    }

    pub fn mode_count(&self) -> usize { self.storage.omegas.len() }
    pub fn omegas(&self) -> &[f64] { &self.storage.omegas }
    pub fn slope_limit(&self) -> f64 { self.stretching.maximum_slope }
    /// Transverse components sharing this one physical elongation and tension.
    pub const fn transverse_planes(&self) -> usize { self.transverse_planes }

    /// Work-conjugate derivative of this SAME mean geometric strain. Material
    /// histories use the original modal normalization and both transverse planes.
    pub(crate) fn mean_strain_gradient(&self, q: &[f64], out: &mut [f64]) -> Result<(), ImpactError> {
        self.observe(q)?;
        if out.len() != self.mode_count() { return Err(invalid("mean strain gradient requires every original coordinate")); }
        let n = self.mode_count();
        let channel = &self.storage.channels[0];
        for i in 0..n { out[i] = 0.5*channel.coupling[i*n+i]*q[i]; }
        if out.iter().any(|v| !v.is_finite()) { return Err(invalid("mean strain gradient is unrepresentable")); }
        Ok(())
    }

    /// Exact norm of the original diagonal positive mean-strain Hessian.
    pub(crate) fn mean_strain_hessian_norm(&self) -> f64 {
        let n = self.mode_count();
        (0..n).map(|i| 0.5*self.storage.channels[0].coupling[i*n+i]).fold(0.0, f64::max)
    }

    /// Extend the SAME stress channel by the affine lift w(s)=d*s*X/L+u(s)
    /// (or d*(1-s/L)*X+u), with unit transverse direction d. This also applies
    /// to vector u: since it vanishes at both ends, integral(u')=0:
    /// strain = sum(e_n*q_n^2)/4 + X^2/(2 L^2). The linear lift stores
    /// T0*X^2/(2L); its curvature is zero, so EI is not counted again.
    /// Only potential/gradient/Hessian/observation are consumed by the support.
    pub(super) fn with_affine_endpoint(self, length_m: f64) -> Result<Self, ImpactError> {
        let n = self.mode_count();
        let stiffness = self.rest_tension_n / length_m;
        let strain_weight = 2.0 / length_m / length_m;
        if self.affine_endpoint || n >= MAX_IMPACT_MODES
            || !length_m.is_finite() || length_m <= 0.0
            || !stiffness.is_finite() || stiffness <= 0.0
            || !strain_weight.is_finite() || strain_weight <= 0.0 {
            return Err(invalid("affine string endpoint exceeds finite geometry or coordinate budget"));
        }
        let mut storage = SosModalStorage {
            omegas: self.storage.omegas.clone(),
            channels: self.storage.channels.clone(),
        };
        let old = &storage.channels[0].coupling;
        let mut coupling = vec![0.0; (n + 1) * (n + 1)];
        for i in 0..n { coupling[i*(n+1)..i*(n+1)+n].copy_from_slice(&old[i*n..(i+1)*n]); }
        coupling[n*(n+1)+n] = strain_weight;
        storage.channels[0].coupling = coupling;
        storage.omegas.push(det::sqrt(stiffness));
        Ok(Self { storage: Arc::new(storage), affine_endpoint: true, ..self })
    }

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
            let e = ch.coupling[i*n+i]; s += e*q[i]*q[i];
            if self.transverse_planes == 1 {
                // Sine modes have mean-square slope = peak^2/2; a linear lift
                // has constant slope. Do not give its bound an artificial sqrt(2).
                let weight = if self.affine_endpoint && i + 1 == n { 0.5 * e } else { e };
                slope_bound += det::sqrt(weight)*q[i].abs();
            }
        }
        if self.transverse_planes == 2 {
            // Triangle inequality for the continuous vector slope. Pair modes
            // before taking their norm so rotating the transverse basis changes
            // neither the validity envelope nor the physical strain.
            let per_plane = (n - usize::from(self.affine_endpoint)) / 2;
            for i in 0..per_plane {
                slope_bound += det::sqrt(ch.coupling[i*n+i])
                    * det::hypot(q[i], q[per_plane+i]);
            }
            if self.affine_endpoint {
                // The support moves along a unit vector in this transverse
                // frame. Its constant vector slope has norm |X|/L, whatever
                // that direction; the lift is one physical coordinate.
                slope_bound += det::sqrt(0.5 * ch.coupling[n*n-1]) * q[n-1].abs();
            }
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
