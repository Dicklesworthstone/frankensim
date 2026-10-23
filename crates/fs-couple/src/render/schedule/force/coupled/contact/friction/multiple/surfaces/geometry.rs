//! Held, modulated tangential velocity/force ports on the original mechanics.
//!
//! The supplied map B is held over THIS interval: surface velocity is B^T qdot
//! and reaction is B F. This is a power-conjugate velocity transformer, NOT the
//! derivative of B(t)^T q and not an extra elastic contact coordinate. Changing
//! B creates no stored energy, impulse or artificial Bdot*q sliding velocity.
//! Normal-contact geometry retains its separate product-rule work accounting.
use super::*;
use super::obstacles::geometry::ObstacleGeometry;

/// One held tangential map in the original component bases. Only geometry may
/// change: component identities, material/interface law and limits stay fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct FrictionGeometry {
    pub left: ModalAttachment,
    pub right: Option<ModalAttachment>,
}

/// Explicit per-update work allowance, including mixed normal/tangent responses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrictionGeometryConfig {
    pub max_update_terms: usize,
}

impl MultiFrictionModalSystem {
    /// Check the complete response rebuild before allocating any new maps.
    /// Includes exact-operator material and radiation memory, not only modes.
    pub fn validate_friction_geometry_config(&self, config: FrictionGeometryConfig)
        -> Result<(), ModalFrictionError>
    {
        let n = self.mode_count();
        let k = self.core.network.columns.len();
        let s = self.port_count().checked_add(self.obstacle_count())
            .ok_or_else(|| invalid("moving friction port count overflow"))?;
        let d = n.checked_mul(2).and_then(|v| v.checked_add(self.material_memory().len()))
            .and_then(|v| v.checked_add(self.core.network.radiation_state().len()))
            .ok_or_else(|| invalid("moving friction state dimension overflow"))?;
        let dense = d.checked_mul(d).and_then(|v| v.checked_mul(3));
        let linked = k.checked_add(1).and_then(|v| v.checked_mul(v))
            .and_then(|v| n.checked_mul(k+2).and_then(|w| v.checked_add(w)));
        let terms = dense.zip(linked).map(|(a,b)| a.max(b))
            .and_then(|v| v.checked_mul(s))
            .and_then(|v| n.checked_mul(s).and_then(|w| w.checked_mul(s))
                .and_then(|w| v.checked_add(w)))
            // Allow a complete normal-geometry rebuild by the existing owner.
            .and_then(|v| v.checked_mul(2));
        if config.max_update_terms == 0 || terms.is_none_or(|v| v > config.max_update_terms) {
            return Err(invalid("moving friction exceeds max_update_terms").into());
        }
        Ok(())
    }

    /// Solve with new held friction maps and optional original normal controls.
    ///
    /// Geometry, histories, observations and clock publish together. Every
    /// inner normal/friction/DG trial uses the same maps and full coupled force
    /// response. Refusal restores the previous friction maps and invalidates
    /// only scratch normal matrices; material memory and vibration are never
    /// reconstructed. An unchanged map takes the original path exactly.
    ///
    /// A friction port whose load comes from a normal contact must stay fixed:
    /// relocating the complete normal/tangent patch needs a joint spatial model.
    /// Other bows may move alongside those fixed contacts and independent stops.
    /// Spatial/time refinement remains necessary; this is neither continuous
    /// contact kinematics nor an alias-free or allocation-free rendering claim.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_friction_geometry(
        &mut self, external: &[f64], drives: &[FrictionDrive], geometry: &[FrictionGeometry],
        normal_geometry: Option<&[ObstacleGeometry]>, normal_offsets_m: Option<&[f64]>,
        config: FrictionGeometryConfig, gate: &CancelGate,
    ) -> Result<&MultiFrictionFrame, ModalFrictionError> {
        poll(Some(gate))?;
        self.validate_friction_geometry_config(config)?;
        if geometry.len() != self.points.len() {
            return Err(invalid("moving friction needs one map per original port").into());
        }
        let mut changed = false;
        for (i, (next, old)) in geometry.iter().zip(&self.points).enumerate() {
            if next.left.component != old.port.left.component
                || next.right.as_ref().map(|r| r.component)
                    != old.port.right.as_ref().map(|r| r.component) {
                return Err(invalid("moving friction must preserve both original component identities").into());
            }
            let moving = next.left != old.port.left || next.right != old.port.right;
            if moving && self.contact_loads[i].is_some() {
                return Err(invalid("contact-loaded friction cannot migrate without its normal patch").into());
            }
            if moving && self.shear_layers[i].is_some() {
                return Err(invalid("material shear cannot migrate without layer-history transport").into());
            }
            changed |= moving;
        }
        if !changed {
            return self.step_inner(external, drives, normal_offsets_m, normal_geometry, Some(gate));
        }
        let mut ports = Vec::with_capacity(self.points.len());
        for (next, old) in geometry.iter().zip(&self.points) {
            poll(Some(gate))?;
            let port = ModalFrictionPort {
                left: next.left.clone(), right: next.right.clone(), law: old.port.law.clone(),
            };
            ports.push(port);
        }
        let (points, compliance) = self.prepare_surface_maps(ports, gate)?;
        poll(Some(gate))?;
        let normal_snapshot = match &mut self.obstacles {
            Some(set) => set.stage_friction_geometry(&self.core.network, &points, gate)?,
            None => None,
        };
        // All fallible geometry preparation precedes these swaps. The existing
        // mechanical owner is atomic on failure; restore its port caches too.
        let core_port = std::mem::replace(&mut self.core.port, points[0].port.clone());
        let core_column = std::mem::replace(&mut self.core.column, points[0].column.clone());
        let core_compliance = std::mem::replace(&mut self.core.compliance, compliance[0]);
        let points = std::mem::replace(&mut self.points, points);
        let compliance = std::mem::replace(&mut self.compliance, compliance);
        let result = self.step_inner(external, drives, normal_offsets_m, normal_geometry, Some(gate))
            .map(|_| ());
        if let Err(error) = result {
            self.core.port = core_port;
            self.core.column = core_column;
            self.core.compliance = core_compliance;
            self.points = points;
            self.compliance = compliance;
            if let Some(set) = &mut self.obstacles { set.restore_friction_geometry(normal_snapshot); }
            return Err(error);
        }
        Ok(&self.frame)
    }
}
