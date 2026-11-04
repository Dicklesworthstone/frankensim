//! Original molecular reservoirs in the complete modal pressure field.
//!
//! The original material's cp(s)=cp_inf+sum C_j/(1+s*tau_j) supplies the
//! molecular compliance 1/K(s)=(1-R/cp(s))/(rho*R*T). No wall attenuation,
//! fitted modal damping, merged relaxation time or second bulk viscosity enters.
//! With M_inf=L L^T and P=M_inf^-1*p, x=L^T P is the exact consistent-mass
//! pressure coordinate. Append t_i=y_i+sqrt(R/cv_inf)*x_i and each original
//! z_ji. Their loss is sum (z_ji-sqrt(C_j/cp_inf)*t_i)^2/tau_j.
//! Negative symmetric outer products realize that law in the original pHS
//! effort coordinates. All added states have unit storage; q and every radial
//! and axial coordinate remain. The source equilibrium gas and face references
//! stay fixed. These are frozen uniform material properties, not heat/species
//! evolution, Stokes--Kirchhoff absorption or wall boundary-layer physics.

use super::*;
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, det};

mod field;

/// Additional complete-source resources, checked before dense allocation.
#[derive(Clone, Copy, Debug)]
pub struct ModalMolecularBudget {
    /// Acoustic coordinates, temperature coordinates and all active reservoirs.
    pub maximum_states: usize,
    /// Conservative dense preparation and observation payload, excluding allocator overhead.
    pub maximum_payload_bytes: usize,
    /// Conservative scalar work including shared structural admission.
    pub maximum_preparation_work: u64,
}

#[derive(Debug, Clone)]
pub(super) enum MolecularSource {
    Uniform(UniformMolecularSource),
    Field(field::FieldMolecularSource),
}

impl MolecularSource {
    fn states(&self) -> usize {
        match self {
            Self::Uniform(source) => source.states,
            Self::Field(source) => source.states,
        }
    }

    fn uniform(&self) -> Option<&UniformMolecularSource> {
        match self {
            Self::Uniform(source) => Some(source),
            Self::Field(_) => None,
        }
    }
}

#[derive(Debug, Clone)]
pub(super) struct UniformMolecularSource {
    original: RelaxationState,
    /// Original process indices: equal rates stay distinct, zero capacities have no state.
    processes: Vec<usize>,
    rates: Vec<f64>,
    reservoir_scales: Vec<f64>,
    temperature_scale: f64,
    /// Original frozen consistent-mass Cholesky; x=L^T*P.
    mass_lower: Vec<f64>,
    states: usize,
}

fn invalid(what: &'static str) -> ModalTimeError {
    ModalTimeError::Invalid(what)
}

impl ModalGalerkin {
    /// Couple the original molecular heat-capacity reservoirs to every retained
    /// pressure coordinate, keeping the full original spatial basis.
    ///
    /// The canonical acoustic state remains [q; p], now p=M_inf*P. It is
    /// followed by one unit-storage temperature block and one unit-storage
    /// block for EACH positive original reservoir, in original process order.
    /// Pressure is the frozen inverse mass times p, and face work remains P*U.
    /// The original source's mass/inverse_mass accessors return this actual
    /// frozen canonical map. Stiffness and equilibrium face references remain.
    ///
    /// Uniform material is required, including regionally authored copies of
    /// the same exact state. Empty reservoirs leave every matrix bitwise intact.
    /// A material with a different gas, repeated attachment or insufficient
    /// complete-state/work/payload capacity refuses without truncating states.
    /// This adiabatic interior model adds molecular absorption only; wall
    /// conduction/shear and evolving thermodynamic fields are not inferred.
    pub fn with_molecular_relaxation(
        mut self,
        original: RelaxationState,
        budget: ModalMolecularBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        cancelled_at(&mut cancelled, 0)?;
        if self.source.molecular.is_some() {
            return Err(invalid("molecular relaxation is already attached"));
        }
        if self.source.viscothermal.is_some()
            || self.system.state_dim() != 2 * self.source.coordinates()
            || self.system.structure().1.iter().any(|r| *r != 0.0)
        {
            return Err(invalid(
                "attach molecular relaxation to the original lossless source",
            ));
        }
        if original.gas() != self.source.gas()
            || self
                .source
                .regional_gases()
                .is_some_and(|gases| gases.iter().any(|gas| gas != original.gas()))
            || self
                .source
                .gas_field()
                .is_some_and(|field| field.uniform_state().is_none())
        {
            return Err(invalid(
                "molecular modal source needs its exact uniform gas state",
            ));
        }
        let n = self.source.coordinates();
        let count = original
            .processes()
            .iter()
            .filter(|p| p.heat_capacity_j_kg_k > 0.0)
            .count();
        let states = n
            .checked_mul(if count == 0 { 2 } else { 3 + count })
            .ok_or_else(|| invalid("molecular modal state extent overflow"))?;
        let ports = self.source.port_count();
        let block = states
            .checked_add(ports)
            .and_then(|v| v.checked_add(1))
            .and_then(|v| v.checked_mul(2))
            .ok_or_else(|| invalid("molecular exact-flow extent overflow"))?;
        let d = states as u128;
        let payload = 128 * d * d + 64 * d * (ports as u128 + 1);
        let work = 64 * d * d * d + 32 * d * d * (ports as u128 + 1);
        if states > budget.maximum_states
            || block > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || payload > budget.maximum_payload_bytes as u128
            || work > budget.maximum_preparation_work as u128
        {
            return Err(invalid(
                "complete molecular modal source exceeds declared resources",
            ));
        }
        if count == 0 {
            self.source.molecular = Some(MolecularSource::Uniform(UniformMolecularSource {
                original,
                processes: vec![],
                rates: vec![],
                reservoir_scales: vec![],
                temperature_scale: 0.0,
                mass_lower: vec![],
                states,
            }));
            cancelled_at(&mut cancelled, n)?;
            return Ok(self);
        }
        let cv = original.frozen_cv_j_kg_k();
        let cp = cv + original.gas().specific_gas_constant;
        let ratio = original.frozen_sound_speed_m_s() / self.source.gas().sound_speed;
        let scale = ratio * ratio;
        let temperature_scale = det::sqrt(original.gas().specific_gas_constant / cv);
        if [cv, cp, scale, temperature_scale]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "molecular modal thermodynamic scales are unrepresentable",
            ));
        }
        let mut processes = Vec::with_capacity(count);
        let mut rates = Vec::with_capacity(count);
        let mut reservoir_scales = Vec::with_capacity(count);
        for (index, process) in original.processes().iter().enumerate() {
            if process.heat_capacity_j_kg_k == 0.0 {
                continue;
            }
            let rate = process.relaxation_time_s.recip();
            let b = det::sqrt(process.heat_capacity_j_kg_k / cp);
            if !rate.is_finite() || rate <= 0.0 || !b.is_finite() || b <= 0.0 {
                return Err(invalid("molecular modal reservoir is unrepresentable"));
            }
            processes.push(index);
            rates.push(rate);
            reservoir_scales.push(b);
        }
        let mass: Vec<_> = self.source.mass.iter().map(|v| v / scale).collect();
        let inverse_mass: Vec<_> = self.source.inverse_mass.iter().map(|v| v * scale).collect();
        if mass.iter().chain(&inverse_mass).any(|v| !v.is_finite()) {
            return Err(invalid("frozen modal pressure storage is unrepresentable"));
        }
        let factor = cholesky(&mass, n)
            .map_err(|_| invalid("frozen consistent mass is not positive definite"))?;
        let mut mass_lower = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                mass_lower[i * n + j] = factor.l(i, j);
            }
        }
        let (old_j, old_r, old_g) = self.system.structure();
        if old_r.iter().any(|r| *r != 0.0) {
            return Err(invalid(
                "attach molecular relaxation to the original lossless source",
            ));
        }
        let old = 2 * n;
        let mut q = vec![0.0; states * states];
        let mut j = vec![0.0; states * states];
        let mut r = vec![0.0; states * states];
        let mut g = vec![0.0; states * ports];
        for i in 0..old {
            cancelled_at(&mut cancelled, i)?;
            j[i * states..i * states + old].copy_from_slice(&old_j[i * old..(i + 1) * old]);
            g[i * ports..(i + 1) * ports].copy_from_slice(&old_g[i * ports..(i + 1) * ports]);
        }
        for i in 0..n {
            for k in 0..n {
                q[i * states + k] = self.source.stiffness[i * n + k];
                q[(n + i) * states + n + k] = inverse_mass[i * n + k];
            }
        }
        for i in old..states {
            q[i * states + i] = 1.0;
        }
        let mut slip = vec![0.0; states];
        for (&b, (process, &rate)) in reservoir_scales.iter().zip(rates.iter().enumerate()) {
            for coordinate in 0..n {
                cancelled_at(&mut cancelled, coordinate)?;
                slip.fill(0.0);
                for row in coordinate..n {
                    slip[n + row] = -b * temperature_scale * mass_lower[row * n + coordinate];
                }
                slip[2 * n + coordinate] = -b;
                slip[(3 + process) * n + coordinate] = 1.0;
                for row in n..states {
                    if slip[row] == 0.0 {
                        continue;
                    }
                    for col in n..=row {
                        r[row * states + col] += rate * slip[row] * slip[col];
                    }
                }
            }
        }
        for i in 0..states {
            for k in 0..i {
                r[k * states + i] = r[i * states + k];
            }
        }
        if q.iter().chain(&r).any(|v| !v.is_finite()) {
            return Err(invalid(
                "molecular modal storage or loss is unrepresentable",
            ));
        }
        let system = PortHamiltonian::new(
            states,
            ports,
            j,
            r,
            g,
            Box::new(QuadraticStorage::new(q, states)?),
        )?;
        cancelled_at(&mut cancelled, states)?;
        self.system = system;
        self.source.mass = mass;
        self.source.inverse_mass = inverse_mass;
        self.source.molecular = Some(MolecularSource::Uniform(UniformMolecularSource {
            original,
            processes,
            rates,
            reservoir_scales,
            temperature_scale,
            mass_lower,
            states,
        }));
        Ok(self)
    }
}

impl ModalGalerkinSource {
    /// Exact uniform material evaluation, retaining every species and process.
    /// A spatial attachment instead exposes each original evaluation through
    /// [`Self::molecular_relaxation_at_quadrature`].
    pub fn molecular_relaxation(&self) -> Option<&RelaxationState> {
        self.molecular
            .as_ref()
            .and_then(MolecularSource::uniform)
            .map(|m| &m.original)
    }
    /// Original uniform process indices for memory blocks, in source order.
    /// Equal rates are never merged; zero capacities remain in the source.
    /// Spatial attachments expose [`Self::molecular_process_indices_at_quadrature`].
    pub fn molecular_process_indices(&self) -> &[usize] {
        self.molecular
            .as_ref()
            .and_then(MolecularSource::uniform)
            .map_or(&[], |m| m.processes.as_slice())
    }
    /// Full pressure-temperature-reservoir state extent.
    pub fn state_dimension(&self) -> usize {
        self.molecular.as_ref().map_or_else(
            || {
                self.viscothermal
                    .as_ref()
                    .map_or(2 * self.coordinates(), |wall| wall.states)
            },
            MolecularSource::states,
        )
    }
    /// All unit-storage thermal and reservoir states. The uniform attachment
    /// uses full-coordinate blocks; a spatial attachment groups its states by
    /// Gauss point, then temperature/process, then radial mode.
    pub fn molecular_state_range(&self) -> Option<core::ops::Range<usize>> {
        self.molecular
            .as_ref()
            .map(|m| 2 * self.coordinates()..m.states())
    }

    /// Independent cycle-mean heat [W] for each retained original process.
    /// Spatial attachments return Gauss-point/process order, summed over modes;
    /// [`Self::molecular_heat_range_at_quadrature`] identifies each point's slice.
    /// Peak phasors use exp(-i omega t); no attenuation is converted into heat.
    pub fn molecular_mean_heat_per_process_w(
        &self,
        state: &[C64],
    ) -> Result<Vec<f64>, ModalTimeError> {
        self.molecular_heat(state, 0.5)
    }

    /// Independent instantaneous nonnegative heat [W] by original process.
    pub fn molecular_heat_per_process_w(&self, state: &[f64]) -> Result<Vec<f64>, ModalTimeError> {
        if state.len() != self.state_dimension() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "molecular heat needs the complete finite original state",
            ));
        }
        let complex: Vec<_> = state.iter().map(|&v| C64::from_re(v)).collect();
        self.molecular_heat(&complex, 1.0)
    }

    fn molecular_heat(&self, state: &[C64], factor: f64) -> Result<Vec<f64>, ModalTimeError> {
        if state.len() != self.state_dimension()
            || state.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(invalid(
                "molecular heat needs the complete finite original state",
            ));
        }
        let Some(memory) = &self.molecular else {
            return Ok(vec![]);
        };
        let memory = match memory {
            MolecularSource::Uniform(memory) => memory,
            MolecularSource::Field(memory) => return memory.heat(self, state, factor),
        };
        if memory.processes.is_empty() {
            return Ok(vec![]);
        }
        let n = self.coordinates();
        let pressure: Vec<_> = (0..n)
            .map(|i| {
                (0..n).fold(C64::ZERO, |sum, j| {
                    sum + state[n + j].scale(self.inverse_mass[i * n + j])
                })
            })
            .collect();
        let mut temperature = vec![C64::ZERO; n];
        for i in 0..n {
            let x = (i..n).fold(C64::ZERO, |sum, j| {
                sum + pressure[j].scale(memory.mass_lower[j * n + i])
            });
            temperature[i] = state[2 * n + i] + x.scale(memory.temperature_scale);
        }
        let result: Vec<_> = memory
            .rates
            .iter()
            .zip(&memory.reservoir_scales)
            .enumerate()
            .map(|(j, (&rate, &b))| {
                factor
                    * rate
                    * (0..n)
                        .map(|i| (state[(3 + j) * n + i] - temperature[i].scale(b)).norm_sq())
                        .sum::<f64>()
            })
            .collect();
        if result.iter().any(|v| !v.is_finite()) {
            return Err(invalid("molecular modal heat is unrepresentable"));
        }
        Ok(result)
    }
}
