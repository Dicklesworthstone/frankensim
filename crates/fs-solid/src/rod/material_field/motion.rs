//! Frozen source-bound Cosserat configurations for native dynamic consumers.
//!
//! This is a force/inertia adapter, not another rod law or time integrator.
//! Inertia uses endpoint row lumping of the ORIGINAL reference inventory and
//! physical cross-section rotary inertia. Material frames have the axial x
//! principal axis. No deformed-length mass update or equivalent density is used.
//! Offset-link inertia would couple translation and rotation; this diagonal
//! lane refuses such links instead of discarding their eccentric inertia.

use std::sync::Arc;
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_material::section::ConcentricShearSection;
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};
use crate::rod::{Rod, derivatives::RodEnergy};
use super::{SpatialMaterialRod, SpatialRodError, checkpoint};

/// Explicit root policy; releasing the static root never happens implicitly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RodRootBoundary { Clamped, Free }

/// Caller-owned candidate; every evaluation validates dimensions and supports.
#[derive(Clone, Debug, PartialEq)]
pub struct RodMotionPose {
    pub positions_m: Vec<[f64; 3]>,
    pub frames: Vec<So3>,
}

/// Row-lumped translational mass and cross-section principal rotary inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodNodalInertia {
    mass_kg: f64,
    principal_kg_m2: [f64; 3],
}
impl RodNodalInertia {
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    pub const fn principal_kg_m2(&self) -> [f64; 3] { self.principal_kg_m2 }
}

/// Complete native force, including fixed nodes, for impulse/reaction accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct RodMotionResponse {
    /// Elastic/link storage PLUS the original gravitational potential [J].
    pub potential_energy_j: f64,
    /// Force ON each node, world frame [N].
    pub force_n: Vec<[f64; 3]>,
    /// Torque ON each node, its current body frame [N m].
    pub torque_body_nm: Vec<[f64; 3]>,
}

/// Immutable sources, supports and reference inertia shared across a trajectory.
///
/// Geometry and thermal/material states are frozen declarations. Dynamics may
/// change the pose, not the sources. Existing 64-segment and link budgets apply.
/// Cross-section inertia is not the rod's orbital inertia, which the consumer
/// obtains separately from nodal positions and linear momenta.
#[derive(Clone, Debug)]
pub struct PreparedRodMotion {
    source: Arc<SpatialMaterialRod>,
    inertia: Vec<RodNodalInertia>,
    fixed_position: Vec<[bool; 3]>,
    fixed_frame: Vec<bool>,
}

impl SpatialMaterialRod {
    /// Freeze this specimen for a time-domain consumer without copying its laws.
    /// Only centerline-attached links are admitted by this diagonal inertia lane.
    /// Explicit finite frequency validity bands refuse: transient forcing is not
    /// certified band-limited merely because its timestep is small.
    pub fn into_motion_model(self, root: RodRootBoundary, gate: &CancelGate)
        -> Result<PreparedRodMotion, SpatialRodError> {
        checkpoint(gate)?;
        for section in &self.sections { admit_section(section)?; checkpoint(gate)?; }
        let reference = self.reference_mass();
        let mut inertia: Vec<_> = reference.nodal_mass_kg().iter().map(|&mass_kg|
            RodNodalInertia { mass_kg, principal_kg_m2: [0.0; 3] }).collect();
        let half_length = 0.5 * reference.segment_length_m();
        for (segment, section) in reference.sections().iter().enumerate() {
            checkpoint(gate)?;
            admit_section(section)?;
            let mut diameter = 0.0;
            for (layer, geometry) in section.elastic().layers().iter()
                .zip(section.elastic().layer_integrals()) {
                let rho = layer.material.property(DENSITY_PROPERTY)
                    .ok_or(SpatialRodError::Invalid("reference inertia lost its original density"))?;
                let value = half_length * (rho.value_si() * geometry.second_moment_m4());
                if !value.is_finite() || value <= 0.0 {
                    return Err(SpatialRodError::Invalid("reference rotary inertia is unrepresentable"));
                }
                diameter += value;
            }
            for node in [segment, segment + 1] {
                for (sum, value) in inertia[node].principal_kg_m2.iter_mut()
                    .zip([2.0 * diameter, diameter, diameter]) { *sum += value; }
            }
        }
        for link in &self.tension_links {
            checkpoint(gate)?;
            if link.attachments().iter().flat_map(|a| a.offset_body_m).any(|v| v != 0.0) {
                return Err(SpatialRodError::Invalid(
                    "dynamic offset links require coupled spatial inertia, not diagonal node inertia"));
            }
            for layer in link.section().layers() {
                for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY] {
                    admit_property(&layer.material, name)?;
                }
            }
            // Declared endpoint lumping of the strand inventory. This is not
            // distributed cable inertia or a dynamic sagging-string model.
            for a in link.attachments() { inertia[a.node].mass_kg += 0.5 * link.reference_mass_kg(); }
        }
        if inertia.iter().any(|i| !i.mass_kg.is_finite() || i.mass_kg <= 0.0
            || i.principal_kg_m2.iter().any(|v| !v.is_finite() || *v <= 0.0)) {
            return Err(SpatialRodError::Invalid("positive finite nodal mass and rotary inertia required"));
        }
        let mut fixed_position = vec![[false; 3]; inertia.len()];
        let mut fixed_frame = vec![false; inertia.len()];
        if root == RodRootBoundary::Clamped { fixed_position[0] = [true; 3]; fixed_frame[0] = true; }
        for support in &self.supports {
            fixed_position[support.node] = support.position_m.map(|v| v.is_some());
            fixed_frame[support.node] = support.frame.is_some();
        }
        let model = PreparedRodMotion { source: Arc::new(self), inertia, fixed_position, fixed_frame };
        model.evaluate(&model.initial_pose(), gate)?;
        checkpoint(gate)?;
        Ok(model)
    }
}

fn admit_section(section: &ConcentricShearSection) -> Result<(), SpatialRodError> {
    for layer in section.elastic().layers() {
        for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, section.shear_property()] {
            admit_property(&layer.material, name)?;
        }
    }
    Ok(())
}
fn admit_property(material: &fs_material::state_point::ResolvedMaterialStatePoint, name: &str)
    -> Result<(), SpatialRodError> {
    let property = material.property(name)
        .ok_or(SpatialRodError::Invalid("dynamic rod lost an original section property"))?;
    if !property.is_constant_at_fixed_temperature()
        || property.angular_frequency_band_rad_s()?.is_some_and(|(lo, hi)| lo > 0.0 || hi < f64::MAX) {
        return Err(SpatialRodError::Invalid(
            "transient rod needs frozen elastic coefficients without a finite frequency validity band"));
    }
    Ok(())
}

impl PreparedRodMotion {
    pub fn source(&self) -> &SpatialMaterialRod { &self.source }
    pub fn inertia(&self) -> &[RodNodalInertia] { &self.inertia }
    pub fn fixed_positions(&self) -> &[[bool; 3]] { &self.fixed_position }
    pub fn fixed_frames(&self) -> &[bool] { &self.fixed_frame }
    pub fn initial_pose(&self) -> RodMotionPose {
        RodMotionPose { positions_m: self.source.rod.positions.clone(), frames: self.source.rod.frames.clone() }
    }

    /// Evaluate full analytical native forces; no tangent matrix is allocated.
    /// Original source bounds, relative-rotation branch cuts, exact stationary
    /// supports and cancellation remain active. No candidate is retained here.
    pub fn evaluate(&self, pose: &RodMotionPose, gate: &CancelGate)
        -> Result<RodMotionResponse, SpatialRodError> {
        checkpoint(gate)?;
        let count = self.inertia.len();
        if pose.positions_m.len() != count || pose.frames.len() != count
            || pose.positions_m.iter().flatten().any(|v| !v.is_finite()) {
            return Err(SpatialRodError::Invalid("dynamic rod pose must match the finite reference mesh"));
        }
        for node in 0..count {
            checkpoint(gate)?;
            for c in 0..3 {
                if self.fixed_position[node][c]
                    && pose.positions_m[node][c] != self.source.rod.positions[node][c] {
                    return Err(SpatialRodError::Invalid("dynamic rod moved a stationary position support"));
                }
            }
            if self.fixed_frame[node] && pose.frames[node] != self.source.rod.frames[node] {
                return Err(SpatialRodError::Invalid("dynamic rod moved a stationary frame support"));
            }
        }
        let rod = Rod { positions: pose.positions_m.clone(), frames: pose.frames.clone(),
            l0: self.source.rod.l0, section: self.source.rod.section };
        let energy = super::derivatives::SpatialEnergy::new(&self.source);
        let potential_energy_j = energy.value(&rod)?;
        let response = energy.linearize(&rod, &[], false, &|| checkpoint(gate))?
            .ok_or(SpatialRodError::Invalid("native dynamic rod requires analytical forces"))?;
        let force_n = response.gradient.iter().map(|g| [-g[0], -g[1], -g[2]]).collect();
        let torque_body_nm = response.gradient.iter().map(|g| [-g[3], -g[4], -g[5]]).collect();
        checkpoint(gate)?;
        Ok(RodMotionResponse { potential_energy_j, force_n, torque_body_nm })
    }
}
