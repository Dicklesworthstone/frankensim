//! Original local caloric reservoirs on the actual spatial quadrature.
//!
//! At point l, w_l=A_l h/(2 rho_l c_inf,l²), x_l=sqrt(w_l) N_l P,
//! a_l=sqrt(R_l/cv_inf,l), b_lj=sqrt(C_lj/cp_inf,l). Unit-storage
//! coordinates y_l,z_lj dissipate (z_lj-b_lj(y_l+a_l x_l))²/tau_lj.
//! Eliminating them gives w_l [1+a_l² S_l/(1+S_l)] = A_l h/(2 K_l(s)),
//! where S_l=sum b_lj²/(1+s tau_lj). This is the original local caloric
//! law, including every distinct process, before spatial integration.

use super::*;
use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium};

#[derive(Debug, Clone)]
struct Point {
    original: RelaxationState,
    processes: Vec<usize>,
    rates: Vec<f64>,
    scales: Vec<f64>,
    temperature_scale: f64,
    sqrt_weight: f64,
    shape: [f64; 2],
    state_start: usize,
    heat_start: usize,
}

#[derive(Debug, Clone)]
pub(in crate::modal::time_domain) struct FieldMolecularSource {
    law: RelaxationAbsorption,
    points: Vec<Point>,
    pub(super) states: usize,
}

impl ModalGalerkin {
    /// Attach the authored molecular law to each ORIGINAL retained gas-field
    /// quadrature state. Temperature, composition, rates and strengths remain
    /// local; no interpolation of relaxation times or averaged gas is used.
    ///
    /// The acoustic prefix is [q; p=M_inf P]. Each Gauss point with at least
    /// one positive reservoir then owns a temperature block and one block per
    /// positive original process, each containing all retained radial modes.
    /// Zero-capacity processes stay in the original material evaluation without
    /// dynamical states. Equal rates are never merged. Empty laws preserve all
    /// original numerical matrices bitwise, while retaining local evaluations.
    ///
    /// This frozen-background, adiabatic molecular model requires the original
    /// lossless field source. Wall thermal diffusion needs a shared temperature
    /// realization and cannot be independently attached to this operator.
    ///
    /// # Errors
    /// Missing field, repeated/mixed attachment, original local material refusal,
    /// complete-state/work/payload capacity, cancellation or structural refusal.
    pub fn with_field_molecular_relaxation(
        mut self,
        law: RelaxationAbsorption,
        budget: ModalMolecularBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ModalTimeError> {
        cancelled_at(&mut cancelled, 0)?;
        let n = self.source.coordinates();
        let old = 2 * n;
        let modes = self.source.modes.0;
        let ports = self.source.port_count();
        if self.source.gas_field().is_none()
            || self.source.molecular.is_some()
            || self.source.viscothermal.is_some()
            || self.system.state_dim() != old
            || self.system.structure().1.iter().any(|r| *r != 0.0)
        {
            return Err(invalid(
                "spatial molecular relaxation needs the original lossless gas-field source",
            ));
        }
        admit(old, ports, budget)?;
        let h = self.source.profile.length() / self.source.elements as f64;
        let offset = 0.5 / det::sqrt(3.0);
        let mut points = Vec::with_capacity(2 * self.source.elements);
        let mut states = old;
        let mut heat_start = 0;
        for element in 0..self.source.elements {
            for (point, fraction) in [0.5 - offset, 0.5 + offset].into_iter().enumerate() {
                cancelled_at(&mut cancelled, 2 * element + point)?;
                let resolved = self
                    .source
                    .gas_at_quadrature(element, point)
                    .ok_or_else(|| {
                        invalid("spatial molecular source is missing its original gas evaluation")
                    })?;
                let original = law.evaluate(&RelaxingMedium::from_thermochemical(resolved)?)?;
                let mut processes = Vec::new();
                let mut rates = Vec::new();
                let mut scales = Vec::new();
                let cv = original.frozen_cv_j_kg_k();
                let gas = original.gas();
                let cp = cv + gas.specific_gas_constant;
                let temperature_scale = det::sqrt(gas.specific_gas_constant / cv);
                let speed = original.frozen_sound_speed_m_s();
                let radius = self
                    .source
                    .profile
                    .radius((element as f64 + fraction) * h)?;
                let weight = 0.5 * h * core::f64::consts::PI * radius * radius
                    / gas.density
                    / (speed * speed);
                let sqrt_weight = det::sqrt(weight);
                if [cv, cp, temperature_scale, weight, sqrt_weight]
                    .iter()
                    .any(|v| !v.is_finite() || *v <= 0.0)
                {
                    return Err(invalid(
                        "local molecular thermodynamic weight is unrepresentable",
                    ));
                }
                for (index, process) in original.processes().iter().enumerate() {
                    if process.heat_capacity_j_kg_k == 0.0 {
                        continue;
                    }
                    let rate = process.relaxation_time_s.recip();
                    let scale = det::sqrt(process.heat_capacity_j_kg_k / cp);
                    if !rate.is_finite() || rate <= 0.0 || !scale.is_finite() || scale <= 0.0 {
                        return Err(invalid("local molecular reservoir is unrepresentable"));
                    }
                    processes.push(index);
                    rates.push(rate);
                    scales.push(scale);
                }
                let state_start = states;
                if !processes.is_empty() {
                    states = states
                        .checked_add(modes * (1 + processes.len()))
                        .ok_or_else(|| invalid("spatial molecular state extent overflow"))?;
                    admit(states, ports, budget)?;
                }
                let count = processes.len();
                points.push(Point {
                    original,
                    processes,
                    rates,
                    scales,
                    temperature_scale,
                    sqrt_weight,
                    shape: [1.0 - fraction, fraction],
                    state_start,
                    heat_start,
                });
                heat_start += count;
            }
        }
        let memory = FieldMolecularSource {
            law,
            points,
            states,
        };
        if states == old {
            cancelled_at(&mut cancelled, memory.points.len())?;
            self.source.molecular = Some(MolecularSource::Field(memory));
            return Ok(self);
        }
        let mut mass = vec![0.0; n * n];
        for (index, point) in memory.points.iter().enumerate() {
            cancelled_at(&mut cancelled, index)?;
            let element = index / 2;
            let weight = point.sqrt_weight * point.sqrt_weight;
            for mode in 0..modes {
                for a in 0..2 {
                    let row = (element + a) * modes + mode;
                    if row >= n {
                        continue;
                    }
                    for b in 0..=a {
                        let column = (element + b) * modes + mode;
                        mass[row * n + column] += weight * point.shape[a] * point.shape[b];
                    }
                }
            }
        }
        for i in 0..n {
            for j in 0..i {
                mass[j * n + i] = mass[i * n + j];
            }
        }
        // Same positive Cholesky/Gram inverse as the original spatial owner.
        let factor = cholesky(&mass, n)
            .map_err(|_| invalid("local frozen consistent mass is not positive definite"))?;
        let mut inverse_l = vec![0.0; n * n];
        for column in 0..n {
            cancelled_at(&mut cancelled, column)?;
            for row in 0..n {
                let mut value = if row == column { 1.0 } else { 0.0 };
                for k in 0..row {
                    value -= factor.l(row, k) * inverse_l[k * n + column];
                }
                inverse_l[row * n + column] = value / factor.l(row, row);
            }
        }
        let mut inverse_mass = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let value = (i..n)
                    .map(|k| inverse_l[k * n + i] * inverse_l[k * n + j])
                    .sum();
                inverse_mass[i * n + j] = value;
                inverse_mass[j * n + i] = value;
            }
        }
        let (old_j, _, old_g) = self.system.structure();
        let mut q = vec![0.0; states * states];
        let mut j = vec![0.0; states * states];
        let mut r = vec![0.0; states * states];
        let mut g = vec![0.0; states * ports];
        let source_count = memory
            .points
            .iter()
            .map(|point| point.rates.len())
            .sum::<usize>()
            * modes;
        let mut source_rows = vec![0.0; source_count * states];
        let mut source_row = 0;
        for row in 0..old {
            j[row * states..row * states + old].copy_from_slice(&old_j[row * old..(row + 1) * old]);
            g[row * ports..(row + 1) * ports]
                .copy_from_slice(&old_g[row * ports..(row + 1) * ports]);
        }
        for row in 0..n {
            for column in 0..n {
                q[row * states + column] = self.source.stiffness[row * n + column];
                q[(n + row) * states + n + column] = inverse_mass[row * n + column];
            }
        }
        for row in old..states {
            q[row * states + row] = 1.0;
        }
        for (index, point) in memory.points.iter().enumerate() {
            cancelled_at(&mut cancelled, index)?;
            for (process, (&rate, &b)) in point.rates.iter().zip(&point.scales).enumerate() {
                for mode in 0..modes {
                    // At most two pressure coefficients and two local memory
                    // coefficients. Accumulate the exact sparse outer product.
                    let mut slip = [(0, 0.0); 4];
                    for side in 0..2 {
                        let coordinate = (index / 2 + side) * modes + mode;
                        if coordinate < n {
                            slip[side] = (
                                n + coordinate,
                                -b * point.temperature_scale
                                    * point.sqrt_weight
                                    * point.shape[side],
                            );
                        }
                    }
                    slip[2] = (point.state_start + mode, -b);
                    slip[3] = (point.state_start + (process + 1) * modes + mode, 1.0);
                    let sqrt_rate = det::sqrt(rate);
                    for &(column, value) in &slip {
                        if value != 0.0 {
                            source_rows[source_row * states + column] = sqrt_rate * value;
                        }
                    }
                    source_row += 1;
                    for &(row, a) in &slip {
                        if a == 0.0 {
                            continue;
                        }
                        for &(column, b) in &slip {
                            if b != 0.0 && column <= row {
                                r[row * states + column] += rate * a * b;
                            }
                        }
                    }
                }
            }
        }
        for row in 0..states {
            for column in 0..row {
                r[column * states + row] = r[row * states + column];
            }
        }
        if q.iter().chain(&r).any(|v| !v.is_finite()) {
            return Err(invalid(
                "local molecular storage or dissipation is unrepresentable",
            ));
        }
        // Retain the ORIGINAL positive loss squares for exact-flow heat.
        // This keeps weak reservoirs and nullspaces through composition; the
        // shared componentwise reconstruction and normal PSD gates both apply.
        let witness =
            fs_phs::ResistanceWitness::try_from_source_rows(states, r.clone(), source_rows)?;
        let system = PortHamiltonian::new(
            states,
            ports,
            j,
            r,
            g,
            Box::new(QuadraticStorage::new(q, states)?),
        )?
        .with_resistance_witness(witness)?;
        cancelled_at(&mut cancelled, states)?;
        self.system = system;
        self.source.mass = mass;
        self.source.inverse_mass = inverse_mass;
        self.source.molecular = Some(MolecularSource::Field(memory));
        Ok(self)
    }
}

fn admit(states: usize, ports: usize, budget: ModalMolecularBudget) -> Result<(), ModalTimeError> {
    let block = states
        .checked_add(ports)
        .and_then(|v| v.checked_add(1))
        .and_then(|v| v.checked_mul(2))
        .ok_or_else(|| invalid("spatial molecular exact-flow extent overflow"))?;
    let d = states as u128;
    if states > budget.maximum_states
        || block > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
        || 128 * d * d + 64 * d * (ports as u128 + 1) > budget.maximum_payload_bytes as u128
        || 64 * d * d * d + 32 * d * d * (ports as u128 + 1)
            > budget.maximum_preparation_work as u128
    {
        return Err(invalid(
            "complete spatial molecular source exceeds declared resources",
        ));
    }
    Ok(())
}

impl ModalGalerkinSource {
    /// Authored species/reservoir/rate law for a spatial molecular attachment.
    pub fn molecular_field_law(&self) -> Option<&RelaxationAbsorption> {
        self.field_memory().map(|m| &m.law)
    }

    /// Exact evaluation on the retained original material at this Gauss point.
    pub fn molecular_relaxation_at_quadrature(
        &self,
        element: usize,
        point: usize,
    ) -> Option<&RelaxationState> {
        self.molecular_point(element, point).map(|p| &p.original)
    }

    /// Indices of positive-capacity processes into this point's original
    /// evaluation. Equal rates remain separate; zero capacities remain in it.
    pub fn molecular_process_indices_at_quadrature(
        &self,
        element: usize,
        point: usize,
    ) -> Option<&[usize]> {
        self.molecular_point(element, point)
            .map(|p| p.processes.as_slice())
    }

    /// This point's thermal block followed by all positive process blocks,
    /// each in radial-mode order. Empty ranges add no fictitious temperature.
    pub fn molecular_state_range_at_quadrature(
        &self,
        element: usize,
        point: usize,
    ) -> Option<core::ops::Range<usize>> {
        self.molecular_point(element, point).map(|p| {
            let count = if p.processes.is_empty() {
                0
            } else {
                (p.processes.len() + 1) * self.modes.0
            };
            p.state_start..p.state_start + count
        })
    }

    /// This point's slice of the per-process heat observations, in the same
    /// source order as [`Self::molecular_process_indices_at_quadrature`].
    pub fn molecular_heat_range_at_quadrature(
        &self,
        element: usize,
        point: usize,
    ) -> Option<core::ops::Range<usize>> {
        self.molecular_point(element, point)
            .map(|p| p.heat_start..p.heat_start + p.processes.len())
    }

    fn field_memory(&self) -> Option<&FieldMolecularSource> {
        match self.molecular.as_ref()? {
            MolecularSource::Field(memory) => Some(memory),
            MolecularSource::Uniform(_) => None,
        }
    }

    fn molecular_point(&self, element: usize, point: usize) -> Option<&Point> {
        if element >= self.elements || point >= 2 {
            return None;
        }
        self.field_memory()?.points.get(2 * element + point)
    }
}

impl FieldMolecularSource {
    pub(super) fn heat(
        &self,
        source: &ModalGalerkinSource,
        state: &[C64],
        factor: f64,
    ) -> Result<Vec<f64>, ModalTimeError> {
        let n = source.coordinates();
        let modes = source.modes.0;
        let pressure: Vec<_> = (0..n)
            .map(|row| {
                (0..n).fold(C64::ZERO, |sum, column| {
                    sum + state[n + column].scale(source.inverse_mass[row * n + column])
                })
            })
            .collect();
        let mut result = Vec::new();
        for (index, point) in self.points.iter().enumerate() {
            for (process, (&rate, &b)) in point.rates.iter().zip(&point.scales).enumerate() {
                let mut heat = 0.0;
                for mode in 0..modes {
                    let p = (0..2).fold(C64::ZERO, |sum, side| {
                        let coordinate = (index / 2 + side) * modes + mode;
                        sum + pressure
                            .get(coordinate)
                            .copied()
                            .unwrap_or(C64::ZERO)
                            .scale(point.shape[side])
                    });
                    let temperature = state[point.state_start + mode]
                        + p.scale(point.temperature_scale * point.sqrt_weight);
                    let slip = state[point.state_start + (process + 1) * modes + mode]
                        - temperature.scale(b);
                    heat += factor * rate * slip.norm_sq();
                }
                result.push(heat);
            }
        }
        if result.iter().any(|v| !v.is_finite()) {
            return Err(invalid("spatial molecular heat is unrepresentable"));
        }
        Ok(result)
    }
}
