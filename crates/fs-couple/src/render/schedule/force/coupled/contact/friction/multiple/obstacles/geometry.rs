//! Prescribed contact-map motion on the original joint contact owner.
//!
//! For x = N^T q - s, the endpoint product rule is exact:
//! delta_x = mean_N^T delta_q + delta_N^T mean_q - delta_s.
//! The held reaction therefore applies -mean_N R to the original mechanics,
//! while R*(delta_N^T mean_q-delta_s) is prescribed geometry work. The original
//! contact secant and loss still use the complete delta_x. No constitutive law,
//! mechanical state or nonlinear storage is replaced by this adapter.

use super::*;
use super::super::super::super::contact_attachment_column;

/// Endpoint displacement maps in the original mass-normalized component bases.
/// The producer owns physical location, directions and motion limits. This is
/// geometry only: the original contact's material law and clearance are retained.
#[derive(Clone, Debug, PartialEq)]
pub struct ObstacleGeometry {
    /// Updated map on the original left component.
    pub left: ModalAttachment,
    /// Updated map on the original right component; zeros denote a fixed support.
    pub right: ModalAttachment,
}

/// Bound on geometry projection and rebuilding the original mixed responses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ObstacleGeometryConfig {
    /// Per endpoint update: s*max(n*(k+2)+(k+1)^2,9*n^2)+n*s^2+3*n*q, with
    /// n modes, k bilateral links, q normal contacts and s total contact ports.
    /// The bound is checked before allocating geometry buffers. It is reused
    /// once per mechanical sample, not once per outer nonlinear-storage trial.
    pub max_update_terms: usize,
}

pub(super) struct GeometryState {
    pub(super) end_columns: Vec<Vec<f64>>,
    pub(super) force_columns: Vec<Vec<f64>>,
    pub(super) change_columns: Vec<Vec<f64>>,
    pub(super) moving: Vec<bool>,
    pub(super) start_transport: Vec<f64>,
    // Mixed matrices are scratch. After moving publication they describe the
    // previous midpoint, so a subsequent HOLD must rebuild at the new station.
    matrices_at_accepted: bool,
    config: ObstacleGeometryConfig,
}

impl MultiFrictionModalSystem {
    /// Admit prescribed changing attachment maps before the first sample.
    ///
    /// Component identities and the complete original modal basis stay fixed.
    /// The maps may move over that basis; this is not moving-basis truncation.
    /// Endpoint geometry defines a discrete product-rule coupling, with its
    /// actuator work, original contact potential and nonnegative loss all in
    /// the same acceptance balance. Spatial and temporal refinement remain
    /// necessary; no rigid-impact or arbitrary moving-patch limit is promised.
    ///
    /// A normal mapped to a tangential friction port moves only through
    /// `step_with_friction_geometry`, which declares both parts of the patch.
    /// Normal-only stepping still refuses to detach the two contact maps.
    pub fn with_obstacle_geometry_motion(mut self, config: ObstacleGeometryConfig,
        gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 {
            return Err(invalid("attachment motion binds before the first physical sample").into());
        }
        let set = self.obstacles.as_mut()
            .ok_or_else(|| invalid("attachment motion requires existing compliant obstacles"))?;
        if set.geometry.is_some() {
            return Err(invalid("attachment motion cannot replace an existing binding").into());
        }
        if set.has_residual() {
            // Moving stations change the omitted modes' compliance every step.
            return Err(invalid("attachment motion does not carry the obstacles' residual compliance").into());
        }
        let q = set.points.len(); let n = self.core.network.mode_count();
        let s = q+self.points.len(); let k = self.core.network.columns.len();
        let dense = n.checked_mul(n).and_then(|v| v.checked_mul(9))
            .ok_or_else(|| invalid("contact geometry dense response work overflows"))?;
        let links = k.checked_add(1).and_then(|v| v.checked_mul(v))
            .ok_or_else(|| invalid("contact geometry connection work overflows"))?;
        let terms = k.checked_add(2).and_then(|v| n.checked_mul(v)).and_then(|v| v.checked_add(links))
            .map(|v| v.max(dense))
            .and_then(|v| v.checked_mul(s))
            .and_then(|v| n.checked_mul(s*s).and_then(|w| v.checked_add(w)))
            .and_then(|v| n.checked_mul(q).and_then(|w| w.checked_mul(3))
                .and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("contact geometry update work overflows"))?;
        if terms > config.max_update_terms {
            return Err(invalid("contact geometry exceeds max_update_terms").into());
        }
        let columns: Vec<_> = set.points.iter().map(|p| p.column.clone()).collect();
        set.geometry = Some(GeometryState { end_columns: columns.clone(), force_columns: columns,
            change_columns: vec![vec![0.0; n]; q], moving: vec![false; q],
            start_transport: vec![0.0; q], matrices_at_accepted: true, config });
        Ok(self)
    }

    /// Bound update work retained with this realization.
    #[must_use]
    pub fn obstacle_geometry_config(&self) -> Option<ObstacleGeometryConfig> {
        Some(self.obstacles.as_ref()?.geometry.as_ref()?.config)
    }

    /// Last accepted complete normal map; the original source remains available
    /// through `obstacle`. Rejected endpoints never alter this map.
    #[must_use]
    pub fn obstacle_normal_column(&self, index: usize) -> Option<&[f64]> {
        Some(&self.obstacles.as_ref()?.points.get(index)?.column)
    }

    /// Advance with one endpoint geometry per normal contact and optional normal
    /// translations. Omitted offsets hold accepted offsets. Ordinary `step`
    /// holds both accepted maps and offsets, including after a refused command.
    /// Positive diagonal endpoint/midpoint compliance is required by the
    /// original scalar normal root; an excessive geometry step can refuse.
    /// Acceleration pressure uses the same held scalar reactions projected
    /// through endpoint maps, plus any existing endpoint nonlinear gradient.
    pub fn step_with_obstacle_geometry(&mut self, external: &[f64], drives: &[FrictionDrive],
        endpoints: &[ObstacleGeometry], offsets_m: Option<&[f64]>)
        -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        self.step_inner(external, drives, offsets_m, Some(endpoints), None)
    }

    /// Same complete transaction with cancellation polling and atomic refusal.
    pub fn step_with_obstacle_geometry_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive],
        endpoints: &[ObstacleGeometry], offsets_m: Option<&[f64]>, gate: &CancelGate)
        -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        self.step_inner(external, drives, offsets_m, Some(endpoints), Some(gate))
    }
}
impl ObstacleSet {
    pub(in super::super) fn prepare_geometry(&mut self, network: &CoupledModalSystem,
        friction: &[Point], contact_loads: &[Option<usize>], endpoints: Option<&[ObstacleGeometry]>,
        contact_patch_motion: bool, gate: Option<&CancelGate>) -> Result<(), ModalCouplingError>
    {
        let Some(geometry) = &mut self.geometry else {
            return if endpoints.is_some() { Err(invalid("attachment motion has not been admitted")) }
                else { Ok(()) };
        };
        if endpoints.is_some_and(|e| e.len() != self.points.len()) {
            return Err(invalid("attachment motion needs one geometry per normal contact"));
        }
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            if let Some(endpoints) = endpoints {
                let end = &endpoints[i];
                if end.left.component != point.contact.left.component
                    || end.right.component != point.contact.right.component {
                    return Err(invalid("moving contacts must retain their original component identities"));
                }
                let column = contact_attachment_column(network, &end.left, &end.right)?;
                geometry.end_columns[i].copy_from_slice(&column);
            } else { geometry.end_columns[i].copy_from_slice(&point.column); }
            geometry.moving[i] = geometry.end_columns[i].iter().zip(&point.column)
                .any(|(end, start)| end.to_bits() != start.to_bits());
            if geometry.moving[i] && contact_loads.contains(&Some(i)) && !contact_patch_motion {
                return Err(invalid("moving a contact-loaded normal requires its tangential migration and slip map"));
            }
            for j in 0..point.column.len() {
                if geometry.moving[i] {
                    let change = finite(geometry.end_columns[i][j]-point.column[j])?;
                    geometry.change_columns[i][j] = change;
                    geometry.force_columns[i][j] = finite(point.column[j]+0.5*change)?;
                } else {
                    geometry.change_columns[i][j] = 0.0;
                    geometry.force_columns[i][j] = point.column[j];
                }
            }
        }
        let moving = geometry.moving.iter().any(|x| *x);
        if !moving && geometry.matrices_at_accepted { return Ok(()); }
        // Invalidate BEFORE writing scratch: any partial failure or cancellation
        // must cause a valid retry/hold to reconstruct the accepted geometry.
        geometry.matrices_at_accepted = false;
        let quiet = CancelGate::new_clock_free(); let gate = gate.unwrap_or(&quiet);
        let q = self.points.len(); let p = friction.len();
        for j in 0..q {
            let response = network_response(network, &geometry.force_columns[j], gate)?;
            for i in 0..q {
                poll(Some(gate))?;
                self.nn[i*q+j] = dot(&geometry.end_columns[i], &response)?;
            }
            for i in 0..p {
                poll(Some(gate))?;
                self.tn[i*q+j] = dot(&friction[i].column, &response)?;
            }
            if self.nn[j*q+j] <= 0.0 {
                return Err(invalid("moving normal needs positive endpoint/midpoint network compliance"));
            }
        }
        for j in 0..p {
            let response = network_response(network, &friction[j].column, gate)?;
            for i in 0..q {
                poll(Some(gate))?;
                self.nt[i*p+j] = dot(&geometry.end_columns[i], &response)?;
            }
        }
        poll(Some(gate))?;
        geometry.matrices_at_accepted = !moving;
        Ok(())
    }

    pub(in super::super) fn endpoint_pressure(&self, mut pressure: f64, weights: &[f64],
        gate: Option<&CancelGate>) -> Result<f64, ModalCouplingError>
    {
        if let Some(geometry) = &self.geometry {
            for (i, &reaction) in self.reactions.iter().enumerate() {
                poll(gate)?;
                if geometry.moving[i] && reaction != 0.0 {
                    // The same jointly solved held scalar reaction, projected at
                    // the physical endpoint rather than the work-conjugate mean.
                    // Nonlinear storage adds its own endpoint gradient separately.
                    for (j, &weight) in weights.iter().enumerate() {
                        if weight != 0.0 {
                            let correction = finite((geometry.end_columns[i][j]
                                - geometry.force_columns[i][j])*reaction)?;
                            pressure = finite(pressure-weight*correction)?;
                        }
                    }
                }
            }
        }
        Ok(pressure)
    }
}

impl ObstacleSet {
    // A moved tangent changes BOTH mixed force responses even when every
    // normal stays fixed. Normal geometry, storage and source laws are unchanged.
    pub(in super::super) fn stage_friction_geometry(&mut self, network: &CoupledModalSystem,
        friction: &[Point], gate: &CancelGate)
        -> Result<Option<(Vec<f64>, Vec<f64>)>, ModalCouplingError>
    {
        poll(Some(gate))?;
        if let Some(geometry) = &mut self.geometry {
            // prepare_geometry rebuilds all mixed maps once for this step,
            // including any simultaneous normal endpoint motion.
            geometry.matrices_at_accepted = false;
            return Ok(None);
        }
        let q = self.points.len(); let p = friction.len();
        let mut nt = vec![0.0; q*p]; let mut tn = vec![0.0; p*q];
        for j in 0..q {
            let response = network_response(network, &self.points[j].column, gate)?;
            for i in 0..p { poll(Some(gate))?; tn[i*q+j] = dot(&friction[i].column, &response)?; }
        }
        for j in 0..p {
            let response = network_response(network, &friction[j].column, gate)?;
            for i in 0..q { poll(Some(gate))?; nt[i*p+j] = dot(&self.points[i].column, &response)?; }
        }
        poll(Some(gate))?;
        Ok(Some((std::mem::replace(&mut self.nt, nt), std::mem::replace(&mut self.tn, tn))))
    }

    pub(in super::super) fn restore_friction_geometry(&mut self, old: Option<(Vec<f64>, Vec<f64>)>) {
        if let Some((nt, tn)) = old { self.nt = nt; self.tn = tn; }
        if let Some(geometry) = &mut self.geometry {
            // A refused combined update may have overwritten trial normal maps.
            // The next HOLD/retry must reconstruct from accepted geometry.
            geometry.matrices_at_accepted = false;
        }
    }
}
