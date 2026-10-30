//! Spatial viscothermal wall memory in the original DG state owner.
//!
//! A positive Gauss-point channel has projection v = sqrt(w*m)*M^(-1/2)*phi,
//! with m=C0 for pressure or L0 for flow. Each ZK arm adds the negative outer
//! product [sqrt(g)*v; -sqrt(r)] and stores its own physical memory energy.
//! The same exact exponential advances all fields and histories together.
//!
//! This is the local-cylinder Webster wall approximation on tapered ducts,
//! not the spherical-cone or full 3-D wall problem. The exact Bessel owner
//! supplies the material coefficients (DLMF 10.10.1); the variational
//! auxiliary-state formulation is described by Thibault & Chabassier,
//! JASA 150 (2021), DOI 10.1121/10.0005537.

use super::*;
use fs_phs::{MAX_ZK_DIFFUSIVE_POLES, ZkDiffusive, ZkDiffusiveTerm};

mod update;

/// Independent wall diffusion and spatial-integration resolution axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DgWallResolution {
    /// Positive continued-fraction arms for EACH of viscosity and heat.
    pub poles_per_field: usize,
    /// Gauss points per axial element, from order+1 through 64.
    /// Refine this separately for variable geometry or material fields.
    pub quadrature_points: usize,
    /// Physical angular-frequency band admitted by the exact Bessel owner.
    pub band_rad_s: [f64; 2],
    /// Ceiling for sampled relative complex AND real-loss errors.
    pub maximum_relative_error: f64,
}

/// Constitutive admission over every actual wall quadrature point.
/// This is not a bound on spatial discretization error or on an input's
/// out-of-band content; the caller still owns the spatial and band ladders.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DgWallAdmission {
    /// Supplied independent resolution and band.
    pub resolution: DgWallResolution,
    /// Number of additional retained physical wall-memory coordinates.
    pub memory_dofs: usize,
    /// Largest sampled complex excess-wall error over all material points.
    pub max_complex_error: f64,
    /// Largest sampled real resistance/conductance error.
    pub max_real_loss_error: f64,
}

#[derive(Debug, Clone)]
struct WallChannel {
    field_start: usize,
    /// Original physical quadrature point and radius identify a fixed section.
    point_m: f64,
    radius_m: f64,
    /// sqrt(w*m) for viscous inertia; sqrt(w*m*(gamma-1)) for thermal
    /// compliance. The retained dimensionless Jacobi weight is gas-independent.
    storage_root: f64,
    /// sqrt(w*m) M^(-1/2) phi; m is C0 or L0, not a damping fit.
    projection: Vec<f64>,
    drag_per_s: f64,
    terms: Vec<ZkDiffusiveTerm>,
    memory_start: usize,
}

/// A retained positive quadrature functional of the physical pressure field.
/// Molecular relaxation consumes this original thermal spectral measure; it
/// does not resample geometry or replace the local radius by an average.
#[derive(Debug)]
pub(in crate::dg) struct ThermalQuadrature {
    pub(in crate::dg) field_start: usize,
    pub(in crate::dg) point_m: f64,
    pub(in crate::dg) projection: Vec<f64>,
    pub(in crate::dg) terms: Vec<ZkDiffusiveTerm>,
}

#[derive(Debug, Clone)]
pub(in crate::dg) struct SpatialViscothermal {
    pub(super) admission: DgWallAdmission,
    channels: Vec<WallChannel>,
}

impl SpatialViscothermal {
    /// Original pointwise thermal measure, before molecular replacement.
    pub(in crate::dg) fn thermal_quadrature(&self, field_dofs: usize) -> Vec<ThermalQuadrature> {
        self.channels.iter().filter(|channel| channel.field_start < field_dofs)
            .map(|channel| ThermalQuadrature {
                field_start: channel.field_start,
                point_m: channel.point_m,
                projection: channel.projection.clone(),
                terms: channel.terms.clone(),
            }).collect()
    }
    /// Reindex only the surviving viscous coordinates. Thermal coordinates are
    /// replaced by the joint heat/molecular law after the returned prefix.
    /// Original thermal channels remain available as constitutive metadata.
    pub(in crate::dg) fn split_for_molecular(
        &mut self,
        field_dofs: usize,
    ) -> (Vec<usize>, Vec<ThermalQuadrature>) {
        let mut kept: Vec<_> = (0..2 * field_dofs).collect();
        let mut thermal = Vec::new();
        for channel in &mut self.channels {
            if channel.field_start < field_dofs {
                thermal.push(ThermalQuadrature {
                    field_start: channel.field_start,
                    point_m: channel.point_m,
                    projection: channel.projection.clone(),
                    terms: channel.terms.clone(),
                });
            } else {
                let original = channel.memory_start;
                channel.memory_start = kept.len();
                kept.extend(original..original + channel.terms.len());
            }
        }
        (kept, thermal)
    }
}

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

impl DgDuct {
    /// Attach causal viscous and thermal memory over the original variable
    /// geometry and primitive material-gas field.
    ///
    /// Every positive Gauss point reevaluates the material owner at its own
    /// temperature/composition and samples the original radius function.
    /// Neither inlet properties nor averaged gas coefficients are substituted.
    /// Refine wall quadrature independently of the DG mesh and diffusion poles.
    /// For uniform gas/radius, with_viscothermal remains the exact-mass shortcut;
    /// this explicit spatial route uses quadrature even in that limit.
    ///
    /// Between prescribed material updates the gas field is frozen; walls are
    /// rigid and isothermal. rebind_spatial_viscothermal_gas_field continues
    /// these native histories at a new admitted background state. The local
    /// cylindrical ZK law does not provide mean flow, finite wall temperature,
    /// spherical-cone corrections or a 3-D solution. For uniform material gas,
    /// with_molecular_relaxation can subsequently couple internal molecular
    /// reservoirs to these original thermal quadrature channels. A nonuniform
    /// source field instead uses with_molecular_relaxation_field, evaluating
    /// the same material law at each original point without inlet averaging.
    /// All memories start at zero when an acoustic initial state is projected.
    /// The original ledger retains total energy and still infers cumulative
    /// loss; wall_dissipated_power_w supplies an independent instantaneous rate.
    ///
    /// # Errors
    /// Repeated attachment, source/phase failure, invalid quadrature/pole band,
    /// nonfinite coefficients, or total state above the existing dense cap.
    pub fn with_spatial_viscothermal(
        mut self,
        resolution: DgWallResolution,
    ) -> Result<Self, GasFieldDuctError> {
        if self.compliant_wall.is_some() || self.molecular.is_some() {
            return Err(invalid("attach DG viscothermal memory before compliant wall mechanics").into());
        }
        if self.viscothermal.is_some() {
            return Err(invalid("DG viscothermal memory is already attached").into());
        }
        let geometry = self.geometry.clone();
        let field = self.gas_field.clone();
        let uniform = self.gas;
        self.attach_spatial_viscothermal::<GasFieldDuctError>(&geometry, resolution, &mut |x| {
            match &field {
                Some(field) => Ok(*field.state_at(x)?.state()),
                None => Ok(uniform),
            }
        })?;
        Ok(self)
    }

    fn attach_spatial_viscothermal<E: From<DuctError>>(
        &mut self,
        duct: &Duct,
        resolution: DgWallResolution,
        gas_at: &mut impl FnMut(f64) -> Result<GasState, E>,
    ) -> Result<(), E> {
        let local = self.mesh.order + 1;
        if !(1..=MAX_ZK_DIFFUSIVE_POLES).contains(&resolution.poles_per_field)
            || !(local..=64).contains(&resolution.quadrature_points)
        {
            return Err(invalid("DG wall pole or quadrature count is outside its domain").into());
        }
        // GasState has gamma>1, so both physical fields have exactly this
        // many arms. Check the full state before any dense allocation.
        let memory_dofs =
            2 * self.mesh.elements * resolution.quadrature_points * resolution.poles_per_field;
        let acoustic_dofs = 2 * self.field_dofs;
        let n = acoustic_dofs + memory_dofs;
        if n > MAX_DG_PROPAGATOR_DOFS {
            return Err(
                invalid("DG acoustic plus wall-memory DOFs exceed MAX_DG_PROPAGATOR_DOFS").into(),
            );
        }
        let mut spans = Vec::with_capacity(duct.segments.len());
        let mut x_start = 0.0;
        for segment in &duct.segments {
            spans.push(AxialSpan {
                x_start,
                length: segment.axial_length(),
                geometry_start: 0.0,
                segment: segment.clone(),
            });
            x_start += segment.axial_length();
        }
        let (nodes, weights) = gauss_legendre_rule(resolution.quadrature_points);
        let mut channels = Vec::new();
        let mut cursor = acoustic_dofs;
        let mut admission = DgWallAdmission {
            resolution,
            memory_dofs,
            max_complex_error: 0.0,
            max_real_loss_error: 0.0,
        };
        for (e, element) in self.elements.iter().enumerate() {
            let jac = 0.5 * element.length;
            for (&xi, &weight) in nodes.iter().zip(&weights) {
                let x = element.x_left + jac * (xi + 1.0);
                let span = spans
                    .iter()
                    .find(|span| x >= span.x_start && x <= span.x_start + span.length)
                    .ok_or_else(|| invalid("DG wall quadrature escaped the original geometry"))?;
                let radius = span.radius_at_local(x - span.x_start)?;
                let area = core::f64::consts::PI * radius * radius;
                let gas = gas_at(x)?;
                crate::validate_gas(&gas)?;
                let law = ZkDiffusive::new(
                    radius,
                    gas.density,
                    gas.dynamic_viscosity,
                    gas.gamma,
                    gas.prandtl,
                    resolution.poles_per_field,
                    resolution.band_rad_s,
                    resolution.maximum_relative_error,
                )
                .map_err(|reason| DuctError::WallDiffusionUnresolved {
                    reason,
                    radius_m: radius,
                    x_m: x,
                })?;
                admission.max_complex_error =
                    admission.max_complex_error.max(law.max_complex_error());
                admission.max_real_loss_error =
                    admission.max_real_loss_error.max(law.max_real_loss_error());
                let phi = self.lgl.eval_basis(xi);
                let compliance = area / (gas.density * gas.sound_speed * gas.sound_speed);
                let inertance = gas.density / area;
                for (is_pressure, coefficient, drag, terms) in [
                    (true, compliance, 0.0, law.shunt_terms()),
                    (
                        false,
                        inertance,
                        law.series_drag_per_s(),
                        law.series_terms(),
                    ),
                ] {
                    let mass = jac * weight * coefficient;
                    if !mass.is_finite() || mass <= 0.0 {
                        return Err(
                            invalid("DG wall quadrature mass is not positive and finite").into(),
                        );
                    }
                    let storage_root = det::sqrt(mass)
                        * if is_pressure { det::sqrt(gas.gamma - 1.0) } else { 1.0 };
                    if !storage_root.is_finite() || storage_root <= 0.0 {
                        return Err(invalid("DG wall physical arm storage is unrepresentable").into());
                    }
                    let inverse = if is_pressure {
                        &element.mp_inv_sqrt
                    } else {
                        &element.mq_inv_sqrt
                    };
                    let projection: Vec<f64> = (0..local)
                        .map(|i| {
                            det::sqrt(mass)
                                * (0..local)
                                    .map(|j| inverse[i * local + j] * phi[j])
                                    .sum::<f64>()
                        })
                        .collect();
                    if projection.iter().any(|v| !v.is_finite()) {
                        return Err(invalid("DG wall energy projection is unrepresentable").into());
                    }
                    channels.push(WallChannel {
                        field_start: e * local + if is_pressure { 0 } else { self.field_dofs },
                        point_m: x,
                        radius_m: radius,
                        storage_root,
                        projection,
                        drag_per_s: drag,
                        terms: terms.to_vec(),
                        memory_start: cursor,
                    });
                    cursor += terms.len();
                }
            }
        }
        if cursor != n {
            return Err(invalid("DG wall memory size disagrees with admitted material").into());
        }
        let mut operator = vec![0.0; n * n];
        for i in 0..acoustic_dofs {
            for j in 0..acoustic_dofs {
                operator[i * n + j] = self.l_tilde[i * acoustic_dofs + j];
            }
        }
        for channel in &channels {
            let total_drag = channel.drag_per_s
                + channel
                    .terms
                    .iter()
                    .map(|term| term.gain_per_s)
                    .sum::<f64>();
            for (i, &left) in channel.projection.iter().enumerate() {
                let row = channel.field_start + i;
                for (j, &right) in channel.projection.iter().enumerate() {
                    let col = channel.field_start + j;
                    operator[row * n + col] -= total_drag * left * right;
                }
                for (arm, term) in channel.terms.iter().enumerate() {
                    let memory = channel.memory_start + arm;
                    let coupling = left * det::sqrt(term.gain_per_s) * det::sqrt(term.rate_per_s);
                    operator[row * n + memory] += coupling;
                    operator[memory * n + row] += coupling;
                }
            }
            for (arm, term) in channel.terms.iter().enumerate() {
                let memory = channel.memory_start + arm;
                operator[memory * n + memory] -= term.rate_per_s;
            }
        }
        if operator.iter().any(|v| !v.is_finite()) {
            return Err(invalid("DG wall augmented operator is unrepresentable").into());
        }
        for (_, port) in &mut self.ports {
            port.resize(n, 0.0);
        }
        self.l_tilde = operator;
        self.g_block = None;
        self.memory_dofs = memory_dofs;
        self.viscothermal = Some(ViscothermalMemory::Spatial(SpatialViscothermal {
            admission,
            channels,
        }));
        self.operator_identity = hash_operator(
            self.mesh,
            self.inlet,
            self.outlet,
            &self.elements,
            &self.l_tilde,
            &self.ports,
        );
        Ok(())
    }

}

impl<M: DgMedium> DgDuct<M> {
    /// Admission over every spatial wall quadrature point. A uniform
    /// exact-mass attachment is reported by viscothermal_model instead.
    #[must_use]
    pub fn spatial_viscothermal_admission(&self) -> Option<&DgWallAdmission> {
        match &self.viscothermal {
            Some(ViscothermalMemory::Spatial(wall)) => Some(&wall.admission),
            _ => None,
        }
    }

    /// Whether any physical viscothermal memory is attached, independent of
    /// the uniform exact-mass or spatial quadrature representation.
    #[must_use]
    pub fn has_viscothermal_memory(&self) -> bool {
        self.viscothermal.is_some()
    }

    /// Instantaneous nonnegative wall power loss [W] from physical drag and
    /// relaxation squares, independent of the stored-energy deficit.
    /// Includes an attached compliant wall's material viscosity once.
    /// Matched-boundary absorption is excluded; this is not an integrated
    /// finite-wall heat transfer or temperature-feedback model.
    ///
    /// # Errors
    /// Foreign/nonfinite state or power overflow.
    pub fn wall_dissipated_power_w(&self, state: &DgState) -> Result<f64, DuctError> {
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != self.total_dofs()
            || state.tilde.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "DG wall loss needs a finite state of this exact operator",
            ));
        }
        let mut power = self.compliant_wall_loss(state);
        match &self.viscothermal {
            Some(ViscothermalMemory::Spatial(wall)) => {
                for channel in &wall.channels {
                    if self.molecular_replaces_thermal() && channel.field_start < self.field_dofs {
                        continue;
                    }
                    let value: f64 = channel
                        .projection
                        .iter()
                        .enumerate()
                        .map(|(i, weight)| weight * state.tilde[channel.field_start + i])
                        .sum();
                    power += channel.drag_per_s * value * value;
                    for (arm, term) in channel.terms.iter().enumerate() {
                        let slip = det::sqrt(term.gain_per_s) * value
                            - det::sqrt(term.rate_per_s) * state.tilde[channel.memory_start + arm];
                        power += slip * slip;
                    }
                }
            }
            Some(ViscothermalMemory::Uniform(model)) => {
                let m = self.field_dofs;
                power += model.series_drag_per_s()
                    * state.tilde[m..2 * m].iter().map(|v| v * v).sum::<f64>();
                for (branch, (offset, term)) in model
                    .series_terms()
                    .iter()
                    .map(|t| (m, t))
                    .chain(model.shunt_terms().iter().filter(|_| !self.molecular_replaces_thermal())
                        .map(|t| (0, t)))
                    .enumerate()
                {
                    for i in 0..m {
                        let slip = det::sqrt(term.gain_per_s) * state.tilde[offset + i]
                            - det::sqrt(term.rate_per_s) * state.tilde[2 * m + branch * m + i];
                        power += slip * slip;
                    }
                }
            }
            None => {}
        }
        power += self.molecular_thermal_loss(state);
        if !power.is_finite() {
            return Err(invalid("DG wall loss is unrepresentable"));
        }
        Ok(power)
    }
}

impl DgState {
    /// Retained auxiliary energy [J], excluding the primary gas p/q fields.
    /// This legacy-named getter includes attached molecular/thermal gas memory.
    #[must_use]
    pub fn wall_stored_energy_j(&self) -> f64 {
        0.5 * self.tilde[2 * self.p.len()..]
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
    }
}

#[cfg(test)]
mod tests;
