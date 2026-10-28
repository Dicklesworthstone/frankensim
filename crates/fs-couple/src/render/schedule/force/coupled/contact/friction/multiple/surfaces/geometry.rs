//! Held, modulated tangential velocity/force ports on the original mechanics.
//!
//! The supplied map B is held over THIS interval: surface velocity is B^T qdot
//! and reaction is B F. This is a power-conjugate velocity transformer, NOT the
//! derivative of B(t)^T q and not an extra elastic contact coordinate. Changing
//! B creates no stored energy, impulse or artificial Bdot*q sliding velocity.
//! Normal-contact geometry retains its separate product-rule work accounting.
use super::*;
use super::obstacles::geometry::ObstacleGeometry;
use super::super::super::super::contact_attachment_column;

/// One held tangential map in the original component bases. Only geometry may
/// change: component identities, material/interface law and limits stay fixed.
#[derive(Clone, Debug, PartialEq)]
pub struct FrictionGeometry {
    pub left: ModalAttachment,
    pub right: Option<ModalAttachment>,
}

/// Actual start/end velocity maps for source-domain admission. The held map
/// still supplies the work-conjugate integration and interval-average check.
/// Normal load comes from the original mechanical drive/normal-contact owner;
/// it cannot be overridden or set to zero to bypass endpoint applicability.
#[derive(Clone, Debug, PartialEq)]
pub struct FrictionEndpointGeometry {
    pub start: FrictionGeometry,
    pub end: FrictionGeometry,
    pub start_speed_m_s: f64,
    pub end_speed_m_s: f64,
}

pub(crate) struct EndpointControls {
    start_column: Vec<f64>,
    end_column: Vec<f64>,
    start_speed_m_s: f64,
    end_speed_m_s: f64,
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

    /// Same held-map transaction, checking original material applicability at
    /// the ACTUAL start/end velocities and attachments before publication.
    /// Each column comes from the original signed port admission. The original
    /// solver still checks its work-conjugate average slip and owns every force,
    /// loss and state update; no endpoint force or extra integration is added.
    /// Endpoint column work reserves an additional `2 * modes * ports` terms.
    /// These sampled checks do not certify unsampled continuous-time extrema.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_friction_geometry_and_endpoints(
        &mut self, external: &[f64], drives: &[FrictionDrive], geometry: &[FrictionGeometry],
        endpoints: &[FrictionEndpointGeometry], normal_geometry: Option<&[ObstacleGeometry]>,
        normal_offsets_m: Option<&[f64]>, config: FrictionGeometryConfig, gate: &CancelGate,
    ) -> Result<&MultiFrictionFrame, ModalFrictionError> {
        poll(Some(gate))?;
        if endpoints.len() != self.points.len() || self.endpoint_controls.is_some() {
            return Err(invalid("endpoint applicability needs one complete pair per original port").into());
        }
        let endpoint_terms = self.mode_count().checked_mul(self.points.len())
            .and_then(|terms| terms.checked_mul(2))
            .ok_or_else(|| invalid("friction endpoint work count overflow"))?;
        let config = FrictionGeometryConfig { max_update_terms: config.max_update_terms
            .checked_sub(endpoint_terms).ok_or_else(|| invalid("friction endpoints exceed max_update_terms"))? };
        self.validate_friction_geometry_config(config)?;
        let mut controls = Vec::with_capacity(endpoints.len());
        for (endpoint, old) in endpoints.iter().zip(&self.points) {
            poll(Some(gate))?;
            if !endpoint.start_speed_m_s.is_finite() || !endpoint.end_speed_m_s.is_finite() {
                return Err(invalid("friction endpoint speeds must be finite").into());
            }
            let column = |map: &FrictionGeometry| {
                if map.left.component != old.port.left.component
                    || map.right.as_ref().map(|r| r.component)
                        != old.port.right.as_ref().map(|r| r.component) {
                    return Err(invalid("endpoint maps must preserve original contact participants").into());
                }
                scalar::admit_port(&self.core.network, &ModalFrictionPort {
                    left: map.left.clone(), right: map.right.clone(), law: old.port.law.clone(),
                }, old.config)
            };
            controls.push(EndpointControls {
                start_column: column(&endpoint.start)?, end_column: column(&endpoint.end)?,
                start_speed_m_s: endpoint.start_speed_m_s, end_speed_m_s: endpoint.end_speed_m_s,
            });
        }
        self.endpoint_controls = Some(controls);
        let result = self.step_with_friction_geometry(external, drives, geometry,
            normal_geometry, normal_offsets_m, config, gate).map(|_| ());
        // Clear per-call data after success OR failure. The original geometry
        // transaction has already restored held maps on any refused trial.
        self.endpoint_controls = None;
        result?;
        Ok(&self.frame)
    }

    pub(crate) fn check_interface_endpoint(&self, port: usize,
        models: &[ModalAcousticTimeModel], drive: FrictionDrive, final_endpoint: bool,
    ) -> Result<(), ModalFrictionError> {
        let (column, speed_m_s) = match &self.endpoint_controls {
            Some(controls) if final_endpoint => (&controls[port].end_column,
                controls[port].end_speed_m_s),
            Some(controls) => (&controls[port].start_column, controls[port].start_speed_m_s),
            None => (&self.points[port].column, drive.speed_m_s),
        };
        interface::check_endpoint(self.interface(port), column, models,
            FrictionDrive { speed_m_s, normal_force_n: drive.normal_force_n })
    }

    /// Solve with new held friction maps and optional original normal controls.
    ///
    /// Geometry, histories, observations and clock publish together. Every
    /// inner normal/friction/DG trial uses the same maps and full coupled force
    /// response. Refusal restores the previous friction maps and invalidates
    /// only scratch normal matrices; material memory and vibration are never
    /// reconstructed. An unchanged map takes the original path exactly.
    ///
    /// A friction port whose load comes from a normal contact may move only
    /// with an explicit complete `normal_geometry` endpoint. Its normal uses
    /// the existing product-rule work account and the same solved reaction
    /// supplies tangential friction. The geometry producer remains responsible
    /// for physical co-location and compatible normal/tangent directions, as
    /// with the original `ContactNormalLoad` binding: modal maps alone cannot
    /// prove these. Retained material shear layers still require transport of
    /// their patch histories and refuse migration.
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
            if moving && self.contact_loads[i].is_some() && normal_geometry.is_none() {
                return Err(invalid("contact-loaded friction cannot migrate without its normal patch").into());
            }
            if moving && self.shear_layers[i].is_some() {
                return Err(invalid("material shear cannot migrate without layer-history transport").into());
            }
            if self.shear_layers[i].is_some() {
                if let (Some(contact), Some(endpoints)) = (self.contact_loads[i], normal_geometry) {
                    let next = endpoints.get(contact).ok_or_else(||
                        invalid("moving friction needs the complete normal geometry"))?;
                    let column = contact_attachment_column(&self.core.network, &next.left, &next.right)?;
                    let accepted = self.obstacle_normal_column(contact)
                        .expect("bound normal belongs to the original obstacle set");
                    if column.iter().zip(accepted).any(|(a, b)| a.to_bits() != b.to_bits()) {
                        return Err(invalid("material shear cannot migrate without layer-history transport").into());
                    }
                }
            }
            changed |= moving;
        }
        if !changed {
            return self.step_with_patch_geometry_inner(external, drives, normal_offsets_m,
                normal_geometry, normal_geometry.is_some(), Some(gate));
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
        let result = self.step_with_patch_geometry_inner(external, drives, normal_offsets_m,
            normal_geometry, normal_geometry.is_some(), Some(gate))
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
