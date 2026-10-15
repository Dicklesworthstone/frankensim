//! Source-bound normal contacts in the existing joint impact time equation.
//!
//! Material resolution and Hertz geometry stay in fs-contact; the same
//! fs-dcontact potential and passive port feed the original fs-phs stepper.
//! Source envelopes survive reference, prepared, analytic and substepped use.
/// Whole homogeneous specimens using the same bulk states for inertia and contact.
pub mod pair;
/// Atomic deposition of actual contact work into finite material enthalpy.
pub mod heating;

use fs_contact::normal_patch::{NormalPatchError, PointNormalPatchReceipt};
use fs_contact::normal_patch::source::SourcePointContact;
use fs_dcontact::Obstacle;
use super::{ImpactBody, ImpactConfig, ImpactError, ImpactSystem, MAX_IMPACT_MODES,
    VolumeSpring, felt::FeltPad, invalid};

/// One physical point contact in the concatenated mass-normalized body basis.
/// Positive closure is weights dot q - gap. Weights must include both moving
/// participants with their physical signs; unrelated coordinates have zeros.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialContact {
    key: String,
    source: SourcePointContact,
    weights: Vec<f64>,
    gap_m: f64,
}
impl MaterialContact {
    /// Bind immutable source/geometry data to a fixed work-conjugate contact map.
    /// There is no duplicate stiffness, exponent or damping argument to go stale.
    pub fn new(key: String, source: SourcePointContact, weights: Vec<f64>, gap_m: f64)
        -> Result<Self, ImpactError> {
        if key.trim().is_empty() || key.len()>128 || weights.is_empty()
            || weights.len()>MAX_IMPACT_MODES || weights.iter().any(|w|!w.is_finite())
            || weights.iter().all(|w|*w==0.0) || !gap_m.is_finite() {
            return Err(invalid("material contact needs a bounded key, finite nonzero map and finite gap"));
        }
        Ok(Self { key, source, weights, gap_m })
    }
    /// Stable physical contact key, unique within its assembled system.
    pub fn key(&self) -> &str { &self.key }
    /// Exact source-bound interface and finite-patch operating box.
    pub const fn source(&self) -> &SourcePointContact { &self.source }
    /// Fixed physical closure map [1/sqrt(kg)].
    pub fn weights(&self) -> &[f64] { &self.weights }
    /// Signed reference gap [m]; a negative value declares precompression.
    pub const fn gap_m(&self) -> f64 { self.gap_m }

    fn lower(&self, modes: usize) -> Result<Obstacle, ImpactError> {
        if self.weights.len()!=modes { return Err(invalid("material contact map must cover the complete body basis")); }
        Obstacle::new(self.weights.clone(),1,modes,vec![self.gap_m],vec![1.0],
            self.source.stiffness_n_per_m_3_2(),1.5,
            format!("source-point-contact/{}/{}",self.key,self.source.model().identity()))
            .and_then(|o|o.with_internal_loss(self.source.dissipation_s_per_m()))
            .map_err(|e|ImpactError::Owner(e.to_string()))
    }

    // For dissipation e is the numerical owner's ACTUAL discrete gradient,
    // not a guessed endpoint or midpoint velocity. For observations it is x.
    fn kinematics(&self, x: &[f64], e: &[f64]) -> Result<(f64,f64), NormalPatchError> {
        if x.len()<2*self.weights.len() || e.len()<2*self.weights.len() {
            return Err(NormalPatchError::InvalidInput { field:"material contact mechanical prefix" });
        }
        let mut closure=-self.gap_m; let mut rate=0.0;
        for (i,w) in self.weights.iter().enumerate() {
            closure+=w*x[2*i];rate+=w*e[2*i+1];
        }
        // Check BEFORE max: f64::max would otherwise hide NaN as a dry contact.
        if !closure.is_finite() || !rate.is_finite() {
            return Err(NormalPatchError::InvalidInput { field:"material contact kinematics" });
        }
        Ok((closure.max(0.0),rate))
    }
}

impl ImpactSystem {
    /// Compose source-defined point contacts with the unchanged body, felt and
    /// cavity storage. Every contact has a complete map in original mode order.
    /// No authored penalty or restitution coefficient can override its source.
    ///
    /// The selected material state and relative local curvature are frozen.
    /// Validity is checked at initialization, both accepted endpoints and the
    /// actual discrete-gradient rate used inside the nonlinear solve. Newton
    /// trials outside the envelope may refuse and require a smaller timestep.
    /// This is a compliant fixed-map model, not collision detection or CCD.
    pub fn new_with_material_contacts(bodies: Vec<ImpactBody>, contacts: Vec<MaterialContact>,
        pads: Vec<FeltPad>, volumes: Vec<VolumeSpring>, config: ImpactConfig)
        -> Result<Self, ImpactError> {
        if contacts.len()>32 { return Err(invalid("at most 32 material contacts are supported")); }
        let modes=bodies.iter().try_fold(0usize,|n,b|n.checked_add(b.potential.count()))
            .ok_or_else(||invalid("material contact mode count overflow"))?;
        let mut raw=Vec::with_capacity(contacts.len());
        for (i,c) in contacts.iter().enumerate() {
            if contacts[..i].iter().any(|previous|previous.key==c.key) {
                return Err(invalid("material contact keys must be unique"));
            }
            c.source.validate_kinematics(0.0,0.0,config.dt_s)
                .map_err(|source|ImpactError::MaterialContact { contact:i, source })?;
            raw.push(c.lower(modes)?);
        }
        let mut system=Self::new(bodies,raw,pads,volumes,config)?;
        system.material_contacts=contacts;
        system.validate_material_contacts(&system.x,&system.x)?;
        Ok(system)
    }

    /// Retained source states in physical contact order, also accessible through
    /// prepared execution's read-only view. Observing never advances time.
    pub fn material_contacts(&self) -> &[MaterialContact] { &self.material_contacts }

    /// Original-owner endpoint force, pressure, footprint and instantaneous loss.
    /// Its power-times-dt is NOT the integrated step's dissipated energy; use
    /// ImpactFrame for the time owner's actual energy balance. This observation
    /// constructs source receipts and is outside allocation-free stepping.
    pub fn material_contact_observation(&self, contact: usize)
        -> Result<PointNormalPatchReceipt, ImpactError> {
        let c=self.material_contacts.get(contact)
            .ok_or_else(||invalid("unknown material contact"))?;
        let (d,v)=c.kinematics(&self.x,&self.x)
            .map_err(|source|ImpactError::MaterialContact { contact,source })?;
        c.source.evaluate(&format!("{}/sample/{}",c.key,self.sample),d,v,self.config.dt_s)
            .map_err(|source|ImpactError::MaterialContact { contact,source })
    }

    pub(super) fn validate_material_contacts(&self, x: &[f64], e: &[f64]) -> Result<(),ImpactError> {
        for (contact,c) in self.material_contacts.iter().enumerate() {
            let (d,v)=c.kinematics(x,e)
                .map_err(|source|ImpactError::MaterialContact { contact,source })?;
            c.source.validate_kinematics(d,v,self.config.dt_s)
                .map_err(|source|ImpactError::MaterialContact { contact,source })?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path="material_tests.rs"]
mod tests;
