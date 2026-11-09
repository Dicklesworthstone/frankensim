//! Original resolved fluid fields loaded by full BEM aperture reaction.
//!
//! The original mixed velocity/pressure/temperature factor retains every fluid,
//! molecular, finite thermal-wall and moving-wall mechanism. A small port solve
//! connects its total-traction/net-flow coordinates to the native exterior's
//! complete self/mutual matrix. No pressure-only substitute or added attenuation
//! changes the interior. All phasors here use exp(+i omega t); the outgoing
//! exp(-i omega t) BEM impedance is conjugated exactly once.
//!
//! Radiators must be original unit, planar, exclusively owned traction ports.
//! Their scalar inward flow Q lifts to exterior outward velocity -Q/A on the
//! SAME physical triangles. Unit-traction admission is not a uniform-velocity
//! theorem: unresolved transverse opening motion is outside this port model.
//! Signed/nonuniform traction modes require an explicit dual map and refuse.
//! The caller declares the lossless homogeneous exterior Medium independently;
//! interior molecular absorption does not silently become exterior absorption.
//! The BEM is rigid outside the selected apertures. Moving structural walls may
//! excite the interior but do not acquire direct exterior radiation here.

use super::{ThermoviscousInterior, ThermoviscousSourceError};
use fs_bem::helmholtz::aperture::{
    Aperture, ApertureError, ApertureImpedance, ApertureRoute, aperture_impedance,
};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_policy::GeometryPolicy;
use fs_exec::CancelGate;
use fs_feec::acoustics::thermoviscous::body_force::ThermoviscousBodyForceSource;
use fs_feec::acoustics::thermoviscous::heat_source::ThermoviscousHeatSource;
use fs_feec::acoustics::thermoviscous::{
    PreparedThermoviscousScattering, ThermoviscousError, ThermoviscousField,
};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use std::collections::BTreeMap;

/// Ceiling for the original physical pressure/flow interface residual.
pub const DEFAULT_PORT_RELATIVE_TOLERANCE: f64 = 1e-8;

/// Explicit cold-work bounds in addition to the native fluid and BEM limits.
#[derive(Clone, Copy, Debug)]
pub struct ThermoviscousRadiationBudget {
    pub maximum_mixed_dofs: usize,
    pub maximum_panels: usize,
    pub maximum_trace_entries: usize,
    pub maximum_projection_terms: u64,
}

/// Source, geometry, original solver or physical-interface refusal.
#[derive(Debug)]
pub enum ThermoviscousRadiationError {
    Invalid(&'static str),
    Source(ThermoviscousSourceError),
    Interior(ThermoviscousError),
    Geometry(HelmholtzError),
    Exterior(ApertureError),
    Closure(LuComplexError),
    Cancelled,
}
impl core::fmt::Display for ThermoviscousRadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved interior/exterior radiation: {self:?}")
    }
}
impl std::error::Error for ThermoviscousRadiationError {}
impl From<ThermoviscousSourceError> for ThermoviscousRadiationError {
    fn from(value: ThermoviscousSourceError) -> Self { Self::Source(value) }
}
impl From<ThermoviscousError> for ThermoviscousRadiationError {
    fn from(value: ThermoviscousError) -> Self { Self::Interior(value) }
}
fn bad(what: &'static str) -> ThermoviscousRadiationError {
    ThermoviscousRadiationError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), ThermoviscousRadiationError> {
    if gate.is_requested() { Err(ThermoviscousRadiationError::Cancelled) } else { Ok(()) }
}
fn finite(z: C64) -> bool { z.re.is_finite() && z.im.is_finite() }
fn triangle_key(points: [[f64; 3]; 3]) -> [[u64; 3]; 3] {
    let mut key = points.map(|p| p.map(|x| if x == 0.0 { 0 } else { x.to_bits() }));
    key.sort_unstable();
    key
}

/// Original source-bound LNS owner and explicitly declared exterior openings.
pub struct RadiatingThermoviscousInterior<'a, 'c, 'surface, S = GasState> {
    interior: &'a ThermoviscousInterior<'c, S>,
    surface: &'surface SpherePanels,
    exterior: Medium,
    radiating_modes: Vec<usize>,
    driven_modes: Vec<usize>,
    apertures: Vec<Aperture>,
}

/// One loaded original fluid LU, full BEM artifact and small port factor.
/// Metadata is owned; the original material/fluid remains borrowed immutably.
pub struct PreparedRadiatingThermoviscousInterior<'a, 'c> {
    interior: PreparedThermoviscousScattering<'a, 'c>,
    source: ApertureImpedance,
    reaction_impedance: Vec<C64>,
    radiating_modes: Vec<usize>,
    driven_modes: Vec<usize>,
    roots: Vec<f64>,
    closure: LuComplex,
}

/// Complete original fields and separately evaluated physical work [W].
#[derive(Clone, Debug)]
pub struct ThermoviscousRadiationResponse {
    pub interior: ThermoviscousField,
    /// Outward volume flows, in the selected radiating-mode order.
    pub outward_flow_m3_s: Vec<C64>,
    /// Exterior reaction Z U in exp(+i omega t), in the same order.
    pub exterior_pressure_pa: Vec<C64>,
    /// Signed work at the prescribed-flow interior feeds only.
    pub source_power_w: f64,
    /// Re(U^H Z U)/2, not an independently added interior damping channel.
    pub exterior_power_w: f64,
    pub port_relative_residual: f64,
    /// All physical source work minus exterior work and native dissipation.
    /// Neither this residual nor conductive availability is caloric heat.
    pub power_balance_residual_w: f64,
}
impl ThermoviscousRadiationResponse {
    pub fn total_source_power_w(&self) -> f64 {
        self.source_power_w + self.interior.moving_wall.external_input_w
            + self.interior.power.body_force_work_w + self.interior.power.heat_source_work_w
    }
}

impl<'a, 'c, 'surface, S> RadiatingThermoviscousInterior<'a, 'c, 'surface, S> {
    /// Bind exact native plane-port triangles to a closed exterior surface.
    /// The independently authored exterior has its own density and sound speed;
    /// no molecular, thermal or viscous exterior law is inferred from them.
    /// All nonradiating native traction modes are prescribed-flow feeds.
    pub fn new(
        interior: &'a ThermoviscousInterior<'c, S>,
        radiating_modes: &[usize],
        surface: &'surface SpherePanels,
        exterior: Medium,
        budget: ThermoviscousRadiationBudget,
        gate: &CancelGate,
    ) -> Result<Self, ThermoviscousRadiationError> {
        poll(gate)?;
        let fe = interior.inner();
        let n = fe.port_mode_count();
        let count = radiating_modes.len();
        let panels = surface.areas().len();
        if n == 0 || count == 0 || count > n || panels == 0
            || fe.ndof() > budget.maximum_mixed_dofs || panels > budget.maximum_panels
            || panels.checked_mul(count).is_none_or(|v| v > budget.maximum_trace_entries)
            || panels as u128 * count as u128 * count as u128 > u128::from(budget.maximum_projection_terms)
        {
            return Err(bad("resolved opening bank exceeds its declared work limits"));
        }
        if !exterior.density.is_finite() || exterior.density <= 0.0
            || !exterior.sound_speed.is_finite() || exterior.sound_speed <= 0.0
        {
            return Err(bad("exterior needs finite positive density and sound speed"));
        }
        let mut selected = vec![false; n];
        for &mode in radiating_modes {
            poll(gate)?;
            if mode >= n || selected[mode] || fe.plane_wave_port_area_m2(mode).is_none() {
                return Err(bad("radiators need unique native unit planar exclusive pressure ports"));
            }
            selected[mode] = true;
        }
        GeometryPolicy::new(surface).map_err(ThermoviscousRadiationError::Geometry)?;
        let triangles = surface.triangles()
            .ok_or_else(|| bad("exterior requires original retained triangles"))?;
        let mut panels_by_triangle = BTreeMap::new();
        for (i, &triangle) in triangles.iter().enumerate() {
            poll(gate)?;
            if panels_by_triangle.insert(triangle_key(triangle), i).is_some() {
                return Err(bad("duplicate physical exterior triangle"));
            }
        }
        let mut apertures = Vec::with_capacity(count);
        for &mode in radiating_modes {
            let mut shape = vec![0.0; panels];
            let mut reference_point = [0.0; 3];
            let mut normal = [0.0; 3];
            let mut area = 0.0;
            for (index, triangle) in fe.plane_wave_port_triangles(mode)
                .expect("native plane port admitted").enumerate()
            {
                poll(gate)?;
                let panel = *panels_by_triangle.get(&triangle_key(triangle.vertices_m))
                    .ok_or_else(|| bad("native LNS opening triangle is absent from exterior"))?;
                if surface.normals()[panel].iter().zip(triangle.outward_normal)
                    .any(|(a, b)| (*a - b).abs() > 1e-12)
                {
                    return Err(bad("native LNS and exterior opening normals disagree"));
                }
                if index == 0 {
                    reference_point = triangle.vertices_m[0];
                    normal = triangle.outward_normal;
                }
                shape[panel] = 1.0;
                area += surface.areas()[panel];
            }
            let native_area = fe.port_areas_m2()[mode];
            if !area.is_finite() || (area - native_area).abs() > 128.0 * f64::EPSILON * native_area {
                return Err(bad("exterior area differs from the original LNS port"));
            }
            apertures.push(Aperture { shape, reference_point, normal });
        }
        poll(gate)?;
        Ok(Self {
            interior, surface, exterior, radiating_modes: radiating_modes.to_vec(),
            driven_modes: (0..n).filter(|&i| !selected[i]).collect(), apertures,
        })
    }

    pub fn interior(&self) -> &'a ThermoviscousInterior<'c, S> { self.interior }
    pub fn exterior_medium(&self) -> Medium { self.exterior }
    pub fn radiating_modes(&self) -> &[usize] { &self.radiating_modes }
    pub fn driven_modes(&self) -> &[usize] { &self.driven_modes }
    pub fn apertures(&self) -> &[Aperture] { &self.apertures }

    /// Factor native loaded fields once, then close P+ZQ=0 at every radiator.
    /// Reference resistances are numerical power coordinates, never added loss.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        reference_resistance: &[f64],
        route: &ApertureRoute,
        gate: &CancelGate,
    ) -> Result<PreparedRadiatingThermoviscousInterior<'a, 'c>, ThermoviscousRadiationError> {
        poll(gate)?;
        let mut cancelled = || gate.is_requested();
        let interior = self.interior.prepare_scattering(frequency_hz, reference_resistance, &mut cancelled)?;
        let source = aperture_impedance(
            self.surface, &self.apertures, &[core::f64::consts::TAU * frequency_hz],
            self.exterior, route, gate,
        ).map_err(ThermoviscousRadiationError::Exterior)?;
        let row = &source.rows[0];
        if !row.passivity_margin.is_finite() || row.passivity_margin < -1e-12 {
            return Err(bad("original BEM aperture matrix has a resolved passivity deficit"));
        }
        let nr = self.radiating_modes.len();
        let areas = self.interior.inner().port_areas_m2();
        let reaction_impedance: Vec<_> = row.impedance.iter().enumerate().map(|(entry, z)| {
            z.conj().scale(areas[self.radiating_modes[entry / nr]].recip())
                .scale(areas[self.radiating_modes[entry % nr]].recip())
        }).collect();
        if reaction_impedance.iter().any(|&z| !finite(z)) {
            return Err(bad("opening pressure/flow normalization is unrepresentable"));
        }
        let scattering = interior.scattering(&mut cancelled)?;
        let roots: Vec<_> = reference_resistance.iter().map(|&r| det::sqrt(r)).collect();
        let n = roots.len();
        let mut matrix = vec![C64::ZERO; n * n];
        for &i in &self.driven_modes {
            poll(gate)?;
            for j in 0..n {
                matrix[i * n + j] = C64::from_re(if i == j { 1.0 } else { 0.0 }) - scattering[i * n + j];
            }
        }
        for (r, &i) in self.radiating_modes.iter().enumerate() {
            poll(gate)?;
            for j in 0..n {
                let mut entry = C64::from_re(if i == j { 1.0 } else { 0.0 }) + scattering[i * n + j];
                for (s, &k) in self.radiating_modes.iter().enumerate() {
                    let flow = C64::from_re(if k == j { 1.0 } else { 0.0 }) - scattering[k * n + j];
                    entry = entry + reaction_impedance[r * nr + s]
                        .scale(roots[i].recip() * roots[k].recip()) * flow;
                }
                matrix[i * n + j] = entry;
            }
        }
        if matrix.iter().any(|&z| !finite(z)) {
            return Err(bad("resolved radiation port closure is unrepresentable"));
        }
        let closure = lu_complex_with_cancel(&matrix, n, &cancelled)
            .map_err(ThermoviscousRadiationError::Closure)?;
        poll(gate)?;
        Ok(PreparedRadiatingThermoviscousInterior {
            interior, source, reaction_impedance, roots, closure,
            radiating_modes: self.radiating_modes.clone(), driven_modes: self.driven_modes.clone(),
        })
    }
}

impl PreparedRadiatingThermoviscousInterior<'_, '_> {
    pub fn interior(&self) -> &PreparedThermoviscousScattering<'_, '_> { &self.interior }
    pub fn source(&self) -> &ApertureImpedance { &self.source }
    pub fn reaction_impedance(&self) -> &[C64] { &self.reaction_impedance }
    pub fn radiating_modes(&self) -> &[usize] { &self.radiating_modes }
    pub fn driven_modes(&self) -> &[usize] { &self.driven_modes }

    pub fn response(
        &self,
        driven_inward_flow_m3_s: &[C64],
        gate: &CancelGate,
    ) -> Result<ThermoviscousRadiationResponse, ThermoviscousRadiationError> {
        self.response_with_volume_sources(driven_inward_flow_m3_s, &[], &[], &[], gate)
    }

    /// Coherent feed flows, original structural forces [N/sqrt(kg)], body
    /// forces [N] and heat release [W], all in the original exp(+i omega t).
    /// Source identity/profile admission belongs to the unchanged native owner.
    pub fn response_with_volume_sources(
        &self,
        driven_inward_flow_m3_s: &[C64],
        wall_forces: &[C64],
        force_sources: &[(&ThermoviscousBodyForceSource<'_, '_>, C64)],
        heat_sources: &[(&ThermoviscousHeatSource<'_, '_>, C64)],
        gate: &CancelGate,
    ) -> Result<ThermoviscousRadiationResponse, ThermoviscousRadiationError> {
        self.response_with_volume_sources_and_tolerance(
            driven_inward_flow_m3_s, wall_forces, force_sources, heat_sources,
            DEFAULT_PORT_RELATIVE_TOLERANCE, gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn response_with_volume_sources_and_tolerance(
        &self,
        driven_inward_flow_m3_s: &[C64],
        wall_forces: &[C64],
        force_sources: &[(&ThermoviscousBodyForceSource<'_, '_>, C64)],
        heat_sources: &[(&ThermoviscousHeatSource<'_, '_>, C64)],
        maximum_port_relative_residual: f64,
        gate: &CancelGate,
    ) -> Result<ThermoviscousRadiationResponse, ThermoviscousRadiationError> {
        poll(gate)?;
        if driven_inward_flow_m3_s.len() != self.driven_modes.len()
            || driven_inward_flow_m3_s.iter().any(|&q| !finite(q))
            || !maximum_port_relative_residual.is_finite()
            || maximum_port_relative_residual <= 0.0 || maximum_port_relative_residual >= 1.0
        {
            return Err(bad("finite feed flows and positive port residual tolerance below one required"));
        }
        let mut cancelled = || gate.is_requested();
        let mut incident = vec![C64::ZERO; self.roots.len()];
        for (&i, &q) in self.driven_modes.iter().zip(driven_inward_flow_m3_s) {
            incident[i] = q.scale(self.roots[i]);
        }
        let has_sources = !wall_forces.is_empty() || !force_sources.is_empty() || !heat_sources.is_empty();
        if has_sources {
            let outgoing = self.interior.source_outgoing_with_volume_sources(
                wall_forces, force_sources, heat_sources, &mut cancelled,
            )?;
            // b=S a+e. Feed rows are (I-S)a=sqrt(R)Q+e; radiators
            // are [(I+S)+R^-1/2 Z R^-1/2(I-S)]a=-e+R^-1/2 Z R^-1/2 e.
            for &i in &self.driven_modes { incident[i] = incident[i] + outgoing[i]; }
            let nr = self.radiating_modes.len();
            for (r, &i) in self.radiating_modes.iter().enumerate() {
                poll(gate)?;
                let mut value = -outgoing[i];
                for (s, &k) in self.radiating_modes.iter().enumerate() {
                    value = value + self.reaction_impedance[r * nr + s]
                        .scale(self.roots[i].recip() * self.roots[k].recip()) * outgoing[k];
                }
                incident[i] = value;
            }
        }
        if incident.iter().any(|&z| !finite(z)) { return Err(bad("source power coordinates overflow")); }
        self.closure.solve_with_cancel(&mut incident, &cancelled)
            .map_err(ThermoviscousRadiationError::Closure)?;
        let interior = if has_sources {
            self.interior.response_with_volume_sources(
                &incident, wall_forces, force_sources, heat_sources, &mut cancelled,
            )?
        } else {
            self.interior.response(&incident, &mut cancelled)?
        };
        let outward_flow_m3_s: Vec<_> = self.radiating_modes.iter()
            .map(|&i| -interior.inward_flow_m3_s[i]).collect();
        let exterior_pressure_pa: Vec<_> = self.reaction_impedance.chunks(self.radiating_modes.len())
            .map(|row| row.iter().zip(&outward_flow_m3_s)
                .fold(C64::ZERO, |sum, (&z, &q)| sum + z * q)).collect();
        let exterior_power_w = 0.5 * exterior_pressure_pa.iter().zip(&outward_flow_m3_s)
            .map(|(&p, &q)| (p.conj() * q).re).sum::<f64>();
        let source_power_w = 0.5 * self.driven_modes.iter().zip(driven_inward_flow_m3_s)
            .map(|(&i, &q)| (interior.pressure_drive_pa[i].conj() * q).re).sum::<f64>();
        let mut scale = 0.0_f64;
        let mut residual = 0.0_f64;
        for (i, (&p, &q)) in interior.pressure_drive_pa.iter().zip(&interior.inward_flow_m3_s).enumerate() {
            let pressure = p.abs() / self.roots[i];
            let flow = q.abs() * self.roots[i];
            if !pressure.is_finite() || !flow.is_finite() { return Err(bad("port residual scale overflow")); }
            scale = scale.max(pressure).max(flow);
        }
        for (&i, &q) in self.driven_modes.iter().zip(driven_inward_flow_m3_s) {
            let error = (interior.inward_flow_m3_s[i] - q).abs() * self.roots[i];
            let input = q.abs() * self.roots[i];
            if !error.is_finite() || !input.is_finite() { return Err(bad("feed residual overflow")); }
            residual = residual.max(error);
            scale = scale.max(input);
        }
        for (&i, &p) in self.radiating_modes.iter().zip(&exterior_pressure_pa) {
            let error = (interior.pressure_drive_pa[i] - p).abs() / self.roots[i];
            let exterior = p.abs() / self.roots[i];
            if !error.is_finite() || !exterior.is_finite() { return Err(bad("exterior residual overflow")); }
            residual = residual.max(error);
            scale = scale.max(exterior);
        }
        let port_relative_residual = if scale == 0.0 { 0.0 } else { residual / scale };
        if port_relative_residual > maximum_port_relative_residual {
            return Err(bad("original LNS/BEM port constraints exceed residual tolerance"));
        }
        let total_source = source_power_w + interior.moving_wall.external_input_w
            + interior.power.body_force_work_w + interior.power.heat_source_work_w;
        let power_balance_residual_w = total_source - exterior_power_w - interior.power.dissipated_w();
        if !source_power_w.is_finite() || !exterior_power_w.is_finite()
            || !total_source.is_finite() || !power_balance_residual_w.is_finite()
        {
            return Err(bad("resolved interior/exterior work overflow"));
        }
        poll(gate)?;
        Ok(ThermoviscousRadiationResponse {
            interior, outward_flow_m3_s, exterior_pressure_pa, source_power_w,
            exterior_power_w, port_relative_residual, power_balance_residual_w,
        })
    }
}

#[cfg(test)]
mod tests;
