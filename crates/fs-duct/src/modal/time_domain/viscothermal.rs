//! Work-conjugate cylindrical wall diffusion on the full retained modal field.
//!
//! The wall basis is psi_n(R)=sgn J0(gamma_n), NOT one independent resistor
//! for every radial mode. For a cylinder, pressure and tangential velocity at
//! the wall are sum s_n P_n and -(1/rho) sum s_n q'_n. Let T and V be their
//! integrated, ZK-normalized Gram factors: T^T T = (A/(rho*c²)) int N^T ss^T N,
//! V^T V = (A/rho) int N'^T ss^T N'. The physical per-area shear impedance is
//! rho*R/2 times the ZK excess series response; thermal admittance is R/(2rho c²)
//! times its excess shunt response. Thus these factors retain real wall work,
//! signed cross-mode coupling, and exactly the original plane-wave ZK law.
//!
//! Shear acts on the existing velocity subspace: K q_dot = K P - V^T z_v.
//! Its effort port is B_v=K^+ V^T, since B_v^T Kq=Vq. Only the known constant
//! plane-potential gauge is fixed in that solve; every acoustic coordinate
//! remains. Thermal ports use B_t=[0;T^T], so B_t^T effort=TP. Each original
//! positive arm g*s/(s+r) appends unit storage y and R's source square
//! (sqrt(g)*B^T effort - sqrt(r)*y)^2. The direct Poiseuille term also retains
//! its original positive work. No negative surface stiffness is introduced.
//!
//! This is a uniform-cylinder, thin-layer modal closure. Its first-order wall
//! work is the shear/tangential-gradient and thermal/pressure law derived in
//! Berggren, Bernland & Noreland, JCP 371 (2018), 633–650,
//! doi:10.1016/j.jcp.2018.06.005. A direct time-domain Wentzell truncation is
//! NOT passive (Hägg & Berggren, arXiv:2207.03152); here the uneliminated
//! positive velocity-port impedance is retained instead. The higher-mode
//! model is asymptotic, not a full radial linearized Navier–Stokes solution.
//! Thin-layer bounds are domain controls, not continuum error certificates.

use super::*;
use crate::modal::bessel_j0_certified;
use fs_math::{c64::C64, det};
use fs_phs::{ResistanceWitness, ZkDiffusive, ZkDiffusiveTerm};

/// Independently refinable wall diffusion and explicit thin-layer domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalViscothermalResolution {
    /// Positive poles in each original viscous/thermal Bessel spectrum.
    pub poles: usize,
    /// Physical angular-frequency interval [rad/s].
    pub band_rad_s: [f64; 2],
    /// Original sampled complex-wall and real-loss relative error limit.
    pub maximum_diffusive_relative_error: f64,
    /// At most 0.1 for both retained-mode/curvature and acoustic length ratios.
    /// This is an authored asymptotic domain, not an error certificate.
    pub maximum_boundary_layer_ratio: f64,
}

/// Complete-state preparation resources. No resolution is lowered to fit.
#[derive(Debug, Clone, Copy)]
pub struct ModalViscothermalBudget {
    /// Every acoustic and retained wall-history coordinate.
    pub maximum_states: usize,
    /// Conservative dense matrices, source factors and preparation scratch.
    pub maximum_payload_bytes: usize,
    /// Conservative scalar preparation work, including structural admission.
    pub maximum_preparation_work: u64,
}

/// Physical and numerical wall domain retained with the original operator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalViscothermalAdmission {
    /// Original authored diffusion resolution and physical domain.
    pub resolution: ModalViscothermalResolution,
    /// max(delta_v,delta_t)*(1+gamma_max)/R at the lowest frequency.
    pub radial_boundary_layer_ratio: f64,
    /// max(delta_v,delta_t)*omega/c at the highest frequency.
    pub acoustic_boundary_layer_ratio: f64,
    /// Complete acoustic and wall-history state count.
    pub states: usize,
}

/// Independently evaluated wall source work [W].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalViscothermalHeat {
    /// Tangential shear diffusion and the original direct viscous drag.
    pub viscous_w: f64,
    /// Thermal exchange with the isothermal cylindrical wall.
    pub thermal_w: f64,
}

#[derive(Debug, Clone)]
pub(super) struct ViscothermalSource {
    model: ZkDiffusive,
    admission: ModalViscothermalAdmission,
    /// Channel rows on the original potential / physical pressure coordinates.
    velocity: Vec<Vec<f64>>,
    pressure: Vec<Vec<f64>>,
    pub(super) states: usize,
}

fn invalid(what: &'static str) -> ModalTimeError {
    ModalTimeError::Invalid(what)
}

/// Exact wall Gram factors on the admitted cylinder. A pressure channel is one
/// column of the original plane-mode consistent-mass Cholesky; a shear channel
/// is one original element derivative. Factoring the common spatial Gram
/// removes no acoustic mode and is algebraically identical to wall quadrature.
#[allow(clippy::too_many_arguments)]
fn wall_channels(
    mass: &[f64],
    modes: usize,
    nodes: usize,
    elements: usize,
    length: f64,
    radius: f64,
    density: f64,
    signs: &[f64],
) -> Result<(Vec<Vec<f64>>, Vec<Vec<f64>>), ModalTimeError> {
    let n = nodes * modes;
    let mut axial_mass = vec![0.0; nodes * nodes];
    for i in 0..nodes {
        for j in 0..nodes {
            axial_mass[i * nodes + j] = mass[(i * modes) * n + j * modes];
        }
    }
    let lower = cholesky(&axial_mass, nodes)
        .map_err(|_| invalid("modal wall pressure Gram matrix is not positive definite"))?;
    let pressure: Vec<_> = (0..nodes)
        .map(|channel| {
            (0..n)
                .map(|i| lower.l(i / modes, channel) * signs[i % modes])
                .collect()
        })
        .collect();
    let scale =
        det::sqrt(core::f64::consts::PI * radius * radius / (density * (length / elements as f64)));
    if !scale.is_finite() || scale <= 0.0 {
        return Err(invalid("modal wall velocity scale is unrepresentable"));
    }
    let velocity: Vec<_> = (0..elements)
        .map(|element| {
            let mut row = vec![0.0; n];
            for mode in 0..modes {
                row[element * modes + mode] = -scale * signs[mode];
                if element + 1 < nodes {
                    row[(element + 1) * modes + mode] = scale * signs[mode];
                }
            }
            row
        })
        .collect();
    Ok((velocity, pressure))
}

/// Lift physical wall velocity ports through the original kinetic storage.
/// For a free endpoint the sole null vector is the constant plane potential.
/// Its known rank-one gauge is admitted and checked, never eigenvalue-clipped.
fn lift_velocity(
    stiffness: &[f64],
    modes: usize,
    end: ModalEnd,
    channels: &[Vec<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Vec<f64>>, ModalTimeError> {
    let n = channels[0].len();
    let nodes = n / modes;
    let mut positive = stiffness.to_vec();
    let gauge = end != ModalEnd::PressureRelease;
    if gauge {
        let scale = (0..n).map(|i| stiffness[i * n + i]).fold(0.0_f64, f64::max);
        for i in 0..n {
            let residual: f64 = (0..nodes).map(|j| stiffness[i * n + j * modes]).sum();
            let magnitude: f64 = (0..nodes).map(|j| stiffness[i * n + j * modes].abs()).sum();
            if residual.abs() > 64.0 * f64::EPSILON * magnitude.max(scale * f64::EPSILON) {
                return Err(invalid(
                    "modal wall stiffness does not retain the known plane gauge",
                ));
            }
        }
        for i in 0..nodes {
            for j in 0..nodes {
                positive[(i * modes) * n + j * modes] += scale / nodes as f64;
            }
        }
    }
    let factor = cholesky(&positive, n)
        .map_err(|_| invalid("modal wall kinetic port lift is not positive off the plane gauge"))?;
    let mut result = Vec::with_capacity(channels.len());
    for (channel, original) in channels.iter().enumerate() {
        cancelled_at(cancelled, channel)?;
        let mut lifted = original.clone();
        factor.solve(&mut lifted);
        if gauge {
            // Exact gauge choice, not removal of an acoustic eigenmode.
            let mean = (0..nodes).map(|i| lifted[i * modes]).sum::<f64>() / nodes as f64;
            for i in 0..nodes {
                lifted[i * modes] -= mean;
            }
        }
        for i in 0..n {
            let actual: f64 = (0..n).map(|j| stiffness[i * n + j] * lifted[j]).sum();
            let magnitude: f64 = (0..n)
                .map(|j| (stiffness[i * n + j] * lifted[j]).abs())
                .sum();
            if !actual.is_finite()
                || (actual - original[i]).abs()
                    > 2e-11 * magnitude.max(original[i].abs()).max(f64::MIN_POSITIVE)
            {
                return Err(invalid(
                    "modal wall velocity lift does not reproduce its physical trace",
                ));
            }
        }
        result.push(lifted);
    }
    Ok(result)
}

fn stamp_square(resistance: &mut [f64], rows: &mut Vec<f64>, row: &[f64]) {
    let states = row.len();
    for i in 0..states {
        for j in i..states {
            let value = row[i] * row[j];
            resistance[i * states + j] += value;
            if i != j {
                resistance[j * states + i] += value;
            }
        }
    }
    rows.extend_from_slice(row);
}

/// Stamp the original positive Bessel terms at the physical effort channels.
fn wall_resistance(
    n: usize,
    states: usize,
    velocity_ports: &[Vec<f64>],
    pressure_ports: &[Vec<f64>],
    model: &ZkDiffusive,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, Vec<f64>), ModalTimeError> {
    let mut r = vec![0.0; states * states];
    let mut rows = Vec::new();
    let mut row = vec![0.0; states];
    let drag = det::sqrt(model.series_drag_per_s());
    for channel in velocity_ports {
        cancelled_at(cancelled, 0)?;
        row.fill(0.0);
        for (i, &coefficient) in channel.iter().enumerate() {
            row[i] = drag * coefficient;
        }
        stamp_square(&mut r, &mut rows, &row);
    }
    let mut memory = 2 * n;
    for (offset, terms, channels) in [
        (0, model.series_terms(), velocity_ports),
        (n, model.shunt_terms(), pressure_ports),
    ] {
        for term in terms {
            let gain = det::sqrt(term.gain_per_s);
            let rate = det::sqrt(term.rate_per_s);
            for channel in channels {
                cancelled_at(cancelled, memory)?;
                row.fill(0.0);
                for (i, &coefficient) in channel.iter().enumerate() {
                    row[offset + i] = gain * coefficient;
                }
                row[memory] = -rate;
                stamp_square(&mut r, &mut rows, &row);
                memory += 1;
            }
        }
    }
    if memory != states || r.iter().chain(&rows).any(|v| !v.is_finite()) {
        return Err(invalid("complete modal wall resistance is unrepresentable"));
    }
    Ok((r, rows))
}

impl ModalGalerkin {
    /// Retain viscous and thermal wall histories on every original m=0 mode.
    ///
    /// Uniform cylindrical geometry and gas are required. Signed wall traces
    /// couple the modes; higher-mode shear depends on tangential velocity, not
    /// a plane-wave damping scalar. Original positive Bessel memories recover
    /// the plane-mode wall law and its finite-band admission. For higher modes
    /// this is a thin-layer closure with an explicit retained-scale domain.
    /// No finite-mode radial Navier–Stokes or continuum error claim is made.
    ///
    /// All original q,p coordinates and physical face ports remain, followed
    /// by series histories (pole, axial element), then thermal histories
    /// (pole, consistent axial mass channel). q is the velocity potential;
    /// physical pressure stays M^-1 p and q_dot=P only in the lossless limit.
    /// Source heat is available independently of a time integrator's balance.
    /// Molecular attachment and wall attachment currently exclude each other:
    /// their shared temperature physics cannot be double counted.
    ///
    /// # Errors
    /// Nonuniform shape/material, repeated or molecular attachment, invalid
    /// thin-layer domain, original Bessel/dense-capacity refusals, unresolved
    /// physical port lift, source-square admission or early/late cancellation.
    pub fn with_viscothermal(
        mut self,
        resolution: ModalViscothermalResolution,
        budget: ModalViscothermalBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        cancelled_at(&mut cancelled, 0)?;
        if self.source.viscothermal.is_some() || self.source.molecular.is_some() {
            return Err(invalid(
                "modal viscothermal walls need the original molecular-free lossless source",
            ));
        }
        let source = &self.source;
        let radius = source.profile.constant_radius_m().ok_or_else(|| {
            invalid("modal viscothermal walls require an exactly cylindrical profile")
        })?;
        let gas = source.gas();
        if source
            .regional_gases()
            .is_some_and(|gases| gases.iter().any(|g| g != gas))
            || source
                .gas_field()
                .is_some_and(|field| field.uniform_state().is_none())
        {
            return Err(invalid(
                "modal viscothermal walls require uniform original gas",
            ));
        }
        let limit = resolution.maximum_boundary_layer_ratio;
        let band = resolution.band_rad_s;
        if !limit.is_finite()
            || limit <= 0.0
            || limit > 0.1
            || band.iter().any(|f| !f.is_finite() || *f <= 0.0)
            || band[0] >= band[1]
        {
            return Err(invalid(
                "modal viscothermal wall band or thin-layer domain is invalid",
            ));
        }
        let diffusivity = gas.dynamic_viscosity / gas.density * 1.0_f64.max(1.0 / gas.prandtl);
        let radial_ratio = det::sqrt(2.0 * diffusivity / band[0])
            * (1.0 + source.radial_eigenvalues.last().copied().unwrap_or(0.0))
            / radius;
        let acoustic_ratio = det::sqrt(2.0 * diffusivity * band[1]) / gas.sound_speed;
        if !radial_ratio.is_finite()
            || !acoustic_ratio.is_finite()
            || radial_ratio > limit
            || acoustic_ratio > limit
        {
            return Err(invalid(
                "retained modal wall scales leave the declared thin-boundary-layer domain",
            ));
        }
        let n = source.coordinates();
        let nodes = source.nodes.len();
        let elements = source.elements;
        let model = ZkDiffusive::new(
            radius,
            gas.density,
            gas.dynamic_viscosity,
            gas.gamma,
            gas.prandtl,
            resolution.poles,
            band,
            resolution.maximum_diffusive_relative_error,
        )?;
        cancelled_at(&mut cancelled, 0)?;
        let memories = elements
            .checked_mul(model.series_terms().len())
            .and_then(|s| {
                nodes
                    .checked_mul(model.shunt_terms().len())
                    .and_then(|t| s.checked_add(t))
            })
            .ok_or_else(|| invalid("modal wall memory extent overflow"))?;
        let states = (2 * n)
            .checked_add(memories)
            .ok_or_else(|| invalid("complete modal wall state extent overflow"))?;
        let ports = source.port_count();
        let block = states
            .checked_add(ports)
            .and_then(|v| v.checked_add(1))
            .and_then(|v| v.checked_mul(2))
            .ok_or_else(|| invalid("complete modal wall exact-flow extent overflow"))?;
        let d = states as u128;
        if states > budget.maximum_states
            || block > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || 160 * d * d + 64 * d * (ports as u128 + 1) > budget.maximum_payload_bytes as u128
            || 80 * d * d * d + 32 * d * d * (ports as u128 + 1)
                > budget.maximum_preparation_work as u128
        {
            return Err(invalid(
                "complete modal wall states exceed declared preparation resources",
            ));
        }
        let signs: Vec<_> = source
            .radial_eigenvalues
            .iter()
            .map(|&gamma| {
                if bessel_j0_certified(gamma) > 0.0 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();
        let (velocity, pressure) = wall_channels(
            &source.mass,
            source.modes.0,
            nodes,
            elements,
            source.profile.length(),
            radius,
            gas.density,
            &signs,
        )?;
        let velocity_ports = lift_velocity(
            &source.stiffness,
            source.modes.0,
            source.end,
            &velocity,
            &mut cancelled,
        )?;
        let (r, rows) = wall_resistance(
            n,
            states,
            &velocity_ports,
            &pressure,
            &model,
            &mut cancelled,
        )?;
        let mut q = vec![0.0; states * states];
        let mut j = vec![0.0; states * states];
        let mut g = vec![0.0; states * ports];
        let (old_j, old_r, old_g) = self.system.structure();
        if self.system.state_dim() != 2 * n || old_r.iter().any(|v| *v != 0.0) {
            return Err(invalid(
                "modal wall attachment requires the complete original lossless operator",
            ));
        }
        for i in 0..n {
            for k in 0..n {
                q[i * states + k] = source.stiffness[i * n + k];
                q[(n + i) * states + n + k] = source.inverse_mass[i * n + k];
            }
        }
        for i in 0..2 * n {
            j[i * states..i * states + 2 * n].copy_from_slice(&old_j[i * 2 * n..(i + 1) * 2 * n]);
            g[i * ports..(i + 1) * ports].copy_from_slice(&old_g[i * ports..(i + 1) * ports]);
        }
        for i in 2 * n..states {
            q[i * states + i] = 1.0;
        }
        let witness = ResistanceWitness::try_from_source_rows(states, r.clone(), rows)?;
        let system = PortHamiltonian::new(
            states,
            ports,
            j,
            r,
            g,
            Box::new(QuadraticStorage::new(q, states)?),
        )?
        .with_resistance_witness(witness)?;
        cancelled_at(&mut cancelled, elements)?;
        self.system = system;
        self.source.viscothermal = Some(ViscothermalSource {
            model,
            admission: ModalViscothermalAdmission {
                resolution,
                radial_boundary_layer_ratio: radial_ratio,
                acoustic_boundary_layer_ratio: acoustic_ratio,
                states,
            },
            velocity,
            pressure,
            states,
        });
        Ok(self)
    }
}

fn channel_heat(
    values: &[C64],
    states: &[C64],
    cursor: &mut usize,
    terms: &[ZkDiffusiveTerm],
    factor: f64,
) -> f64 {
    let mut heat = 0.0;
    for term in terms {
        let gain = det::sqrt(term.gain_per_s);
        let rate = det::sqrt(term.rate_per_s);
        for value in values {
            heat += factor * (value.scale(gain) - states[*cursor].scale(rate)).norm_sq();
            *cursor += 1;
        }
    }
    heat
}

impl ModalGalerkinSource {
    /// Original wall resolution, admitted thin-layer ratios and full state count.
    pub fn viscothermal_admission(&self) -> Option<&ModalViscothermalAdmission> {
        self.viscothermal.as_ref().map(|wall| &wall.admission)
    }
    /// Original positive cylindrical diffusion spectrum, without refitting.
    pub fn viscothermal_model(&self) -> Option<&ZkDiffusive> {
        self.viscothermal.as_ref().map(|wall| &wall.model)
    }
    /// All wall-history coordinates following the original acoustic prefix.
    pub fn viscothermal_state_range(&self) -> Option<core::ops::Range<usize>> {
        self.viscothermal
            .as_ref()
            .map(|wall| 2 * self.coordinates()..wall.states)
    }
    /// Independent instantaneous original viscous/thermal resistor power.
    pub fn viscothermal_heat_w(
        &self,
        state: &[f64],
    ) -> Result<ModalViscothermalHeat, ModalTimeError> {
        if state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("modal wall heat requires finite complete states"));
        }
        self.wall_heat(
            &state.iter().copied().map(C64::from_re).collect::<Vec<_>>(),
            1.0,
        )
    }
    /// Independent cycle-mean heat for peak physical-state phasors.
    pub fn viscothermal_mean_heat_w(
        &self,
        state: &[C64],
    ) -> Result<ModalViscothermalHeat, ModalTimeError> {
        self.wall_heat(state, 0.5)
    }
    fn wall_heat(
        &self,
        state: &[C64],
        factor: f64,
    ) -> Result<ModalViscothermalHeat, ModalTimeError> {
        if state.len() != self.state_dimension()
            || state.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(invalid("modal wall heat requires finite complete states"));
        }
        let Some(wall) = &self.viscothermal else {
            return Ok(ModalViscothermalHeat {
                viscous_w: 0.0,
                thermal_w: 0.0,
            });
        };
        let n = self.coordinates();
        let pressure: Vec<_> = (0..n)
            .map(|i| {
                (0..n).fold(C64::ZERO, |sum, k| {
                    sum + state[n + k].scale(self.inverse_mass[i * n + k])
                })
            })
            .collect();
        let project = |rows: &[Vec<f64>], field: &[C64]| -> Vec<C64> {
            rows.iter()
                .map(|row| {
                    row.iter()
                        .zip(field)
                        .fold(C64::ZERO, |sum, (&a, &x)| sum + x.scale(a))
                })
                .collect()
        };
        let v = project(&wall.velocity, &state[..n]);
        let t = project(&wall.pressure, &pressure);
        let mut cursor = 2 * n;
        let viscous_w =
            factor * wall.model.series_drag_per_s() * v.iter().map(|v| v.norm_sq()).sum::<f64>()
                + channel_heat(&v, state, &mut cursor, wall.model.series_terms(), factor);
        let thermal_w = channel_heat(&t, state, &mut cursor, wall.model.shunt_terms(), factor);
        if !viscous_w.is_finite() || !thermal_w.is_finite() {
            return Err(invalid("modal wall heat is unrepresentable"));
        }
        Ok(ModalViscothermalHeat {
            viscous_w,
            thermal_w,
        })
    }
}

#[cfg(test)]
mod tests;
