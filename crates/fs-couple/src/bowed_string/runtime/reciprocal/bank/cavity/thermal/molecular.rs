//! Joint FE heat diffusion and source molecular caloric storage.
//!
//! The original heat eigenvectors define H from pressure space into heat space.
//! Each retained heat mode carries its own gas temperature and the original
//! local internal reservoirs. Its positive bulk spectrum is computed by the
//! same caloric solver as cylindrical material walls. The unresolved spatial
//! complement R=I-H' H retains the SOURCE adiabatic molecular law. Consequently
//! K(s)=sum_l H_l' K_l(s) H_l + R K_adiabatic(s): heat diffusion is truncated,
//! but molecular storage is neither dropped from R nor counted twice on H.

use super::*;
use crate::cavity_modal::shared::thermal::spectrum;
use fs_material::relaxation::RelaxationState;

/// Additional complete work allowances for molecular/FE caloric coupling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellCavityMolecularOptions {
    /// Gas temperature plus active local reservoirs, at most 128 per heat mode.
    pub maximum_caloric_coordinates: usize,
    /// Complete joint heat and unresolved-complement memory count, at most 4096.
    /// Original native setup/state budgets remain independently binding.
    pub maximum_memory_states: usize,
}

impl ShellCavityThermalLoss {
    pub(in crate::bowed_string::runtime::reciprocal::bank::cavity) fn with_molecular(
        mut self,
        gas: &GasState,
        coupling: &CavityCoupling,
        source: RelaxationState,
        options: ShellCavityMolecularOptions,
        gate: &CancelGate,
    ) -> Result<Self, ShellCavityError> {
        check_cancel(gate)?;
        if source.gas() != gas
            || !(1..=spectrum::MAX_COORDINATES).contains(&options.maximum_caloric_coordinates)
            || !(1..=4096).contains(&options.maximum_memory_states)
        {
            return Err(bad(
                "joint FE caloric law requires the exact source gas and bounded complete memory",
            ));
        }
        let processes: Vec<_> = source
            .processes()
            .iter()
            .filter(|process| process.heat_capacity_j_kg_k > 0.0)
            .collect();
        let width = processes
            .len()
            .checked_add(1)
            .filter(|n| *n <= options.maximum_caloric_coordinates)
            .ok_or_else(|| bad("local gas reservoirs exceed the caloric coordinate allowance"))?;
        if processes.is_empty() {
            if self.memory_states() > options.maximum_memory_states {
                return Err(bad(
                    "original heat memory exceeds the complete caloric allowance",
                ));
            }
            self.molecular = Some(source);
            return Ok(self);
        }
        let bulk = source
            .bulk_relaxation()
            .map_err(|_| bad("original molecular bulk source refused its caloric realization"))?;
        let pressure_count = self.pressure_norms.len();
        let heat_count = self.overlaps.len();
        let memories = heat_count.checked_mul(width).and_then(|n| {
            pressure_count
                .checked_mul(bulk.arms().len())
                .and_then(|m| n.checked_add(m))
        });
        let work = width
            .checked_mul(width)
            .and_then(|n| n.checked_mul(width))
            .and_then(|n| n.checked_mul(heat_count))
            .and_then(|n| {
                pressure_count
                    .checked_mul(pressure_count)
                    .and_then(|p| p.checked_mul(pressure_count))
                    .and_then(|p| n.checked_add(p))
            })
            .and_then(|n| {
                memories.and_then(|m| {
                    pressure_count
                        .checked_mul(pressure_count)
                        .and_then(|p| {
                            p.checked_add(coupling.total_modes().checked_mul(pressure_count)?)
                        })
                        .and_then(|p| p.checked_mul(m))
                        .and_then(|p| n.checked_add(p))
                })
            });
        if memories.is_none_or(|n| n > options.maximum_memory_states)
            || work.is_none_or(|n| n > self.options.maximum_projection_terms)
        {
            return Err(bad(
                "complete joint caloric projection exceeds its memory or work allowance",
            ));
        }
        let frozen = bulk.frozen_bulk_modulus_pa();
        let previous = self.instantaneous_bulk_modulus_pa;
        let pressure_scale = frozen / previous;
        if !pressure_scale.is_finite() || pressure_scale <= 0.0 {
            return Err(bad("frozen molecular pressure scale is unrepresentable"));
        }
        let mut rows = Rows::new(
            coupling,
            &self.pressure_norms,
            frozen,
            options.maximum_memory_states,
        );
        let mut capacities = vec![0.0; width];
        capacities[0] = gas.density * (source.frozen_cv_j_kg_k() + gas.specific_gas_constant);
        for (i, process) in processes.iter().enumerate() {
            capacities[i + 1] = gas.density * process.heat_capacity_j_kg_k;
        }
        if capacities.iter().any(|c| !c.is_finite() || *c <= 0.0) {
            return Err(bad("source caloric capacities are unrepresentable"));
        }
        let roots: Vec<_> = capacities.iter().map(|c| det::sqrt(*c)).collect();
        let mut forcing = vec![0.0; width];
        forcing[0] = roots[0].recip();
        let mut local = vec![0.0; width * width];
        for (i, process) in processes.iter().enumerate() {
            let j = i + 1;
            let conductance = capacities[j] / process.relaxation_time_s;
            local[0] += conductance / capacities[0];
            local[j * width + j] = conductance / capacities[j];
            let coupling = -conductance / (roots[0] * roots[j]);
            local[j] = coupling;
            local[j * width] = coupling;
        }
        for (mode, overlap) in self.overlaps.iter().enumerate() {
            check_cancel(gate)?;
            if overlap.iter().all(|value| *value == 0.0) {
                continue;
            }
            let mut a = local.clone();
            a[0] +=
                gas.thermal_conductivity * self.diffusion.eigenvalues_per_m2[mode] / capacities[0];
            let realized = spectrum::realize(
                gas.pressure,
                gas.pressure / gas.temperature,
                a,
                &forcing,
                &mut || gate.is_requested(),
            )
            .map_err(|reason| {
                if gate.is_requested() {
                    ShellCavityError::Cancelled
                } else {
                    bad(reason)
                }
            })?;
            if (realized.frozen_bulk_modulus_pa - frozen).abs() > 3e-11 * frozen {
                return Err(bad(
                    "joint FE caloric law lost the original frozen molecular storage",
                ));
            }
            for arm in realized.arms {
                rows.push(
                    overlap,
                    arm.increment_pa,
                    arm.relaxation_time_s,
                    Some(mode),
                    gate,
                )?;
            }
        }

        // Do not confuse unresolved heat with equilibrated molecular gas. The
        // unresolved pressure-space form has no retained Fourier heat mode but
        // still owns every local molecular reservoir and its original rate.
        let n = pressure_count;
        let mut remainder = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                remainder[i * n + j] = if i == j { 1.0 } else { 0.0 }
                    - self.overlaps.iter().map(|h| h[i] * h[j]).sum::<f64>();
            }
        }
        let (values, vectors) = fs_la::eigen::jacobi_eigh(&remainder, n);
        check_cancel(gate)?;
        let tolerance = 128.0 * f64::EPSILON * n as f64;
        for k in 0..n {
            check_cancel(gate)?;
            let mut norm = 0.0;
            let mut residual = 0.0_f64;
            for i in 0..n {
                let v = vectors[i * n + k];
                norm += v * v;
                let product: f64 = (0..n)
                    .map(|j| remainder[i * n + j] * vectors[j * n + k])
                    .sum();
                residual = residual.max((product - values[k] * v).abs());
            }
            if !values[k].is_finite()
                || values[k] < -tolerance
                || (norm - 1.0).abs() > 2e-11
                || residual > 2e-11
            {
                return Err(bad(
                    "unresolved heat complement is not a positive resolved pressure form",
                ));
            }
            if values[k] <= 0.0 {
                continue;
            }
            let root = det::sqrt(values[k]);
            let projection: Vec<_> = (0..n).map(|i| root * vectors[i * n + k]).collect();
            for arm in bulk.arms() {
                rows.push(
                    &projection,
                    arm.delta_bulk_modulus_pa,
                    arm.relaxation_time_s,
                    None,
                    gate,
                )?;
            }
        }
        let equilibrium = fs_la::factor::cholesky(&rows.equilibrium, n).map_err(|_| {
            bad("joint molecular and wall compression storage is not positive definite")
        })?;
        let mut equilibrium_rows = vec![vec![0.0; coupling.total_modes()]; n];
        for (i, target) in equilibrium_rows.iter_mut().enumerate() {
            check_cancel(gate)?;
            for (j, spring) in coupling.compression_storage().iter().enumerate() {
                let weight = equilibrium.l(j, i) / self.pressure_norms[j];
                for (v, a) in target.iter_mut().zip(&spring.areas) {
                    *v += weight * a;
                }
            }
        }
        if equilibrium_rows.iter().flatten().any(|v| !v.is_finite()) {
            return Err(bad("joint molecular compression projection overflow"));
        }
        self.instantaneous_bulk_modulus_pa = frozen;
        self.instantaneous_pressure_scale = pressure_scale;
        self.rates = rows.rates;
        self.mode_indices = rows.mode_indices;
        self.memory_rows = rows.memory_rows;
        self.pressure_rows = rows.pressure_rows;
        self.equilibrium_rows = equilibrium_rows;
        self.molecular = Some(source);
        check_cancel(gate)?;
        Ok(self)
    }
}

struct Rows<'a> {
    coupling: &'a CavityCoupling,
    norms: &'a [f64],
    maximum: usize,
    equilibrium: Vec<f64>,
    rates: Vec<f64>,
    mode_indices: Vec<usize>,
    memory_rows: Vec<Vec<f64>>,
    pressure_rows: Vec<Vec<f64>>,
}
impl<'a> Rows<'a> {
    fn new(coupling: &'a CavityCoupling, norms: &'a [f64], bulk: f64, maximum: usize) -> Self {
        let mut equilibrium = vec![0.0; norms.len() * norms.len()];
        for i in 0..norms.len() {
            equilibrium[i * norms.len() + i] = bulk;
        }
        Self {
            coupling,
            norms,
            maximum,
            equilibrium,
            rates: Vec::new(),
            mode_indices: Vec::new(),
            memory_rows: Vec::new(),
            pressure_rows: Vec::new(),
        }
    }
    fn push(
        &mut self,
        overlap: &[f64],
        increment: f64,
        tau: f64,
        mode: Option<usize>,
        gate: &CancelGate,
    ) -> Result<(), ShellCavityError> {
        check_cancel(gate)?;
        if self.rates.len() >= self.maximum
            || !increment.is_finite()
            || increment <= 0.0
            || !tau.is_finite()
            || tau <= 0.0
            || !tau.recip().is_finite()
            || tau.recip() <= 0.0
        {
            return Err(bad(
                "complete caloric memory or positive source rate is unrepresentable",
            ));
        }
        let root = det::sqrt(increment);
        let n = self.norms.len();
        let mut row = vec![0.0; self.coupling.total_modes()];
        let mut pressure = vec![0.0; n];
        for (i, spring) in self.coupling.compression_storage().iter().enumerate() {
            let normalized = root * overlap[i];
            pressure[i] = normalized / self.norms[i];
            for (value, a) in row.iter_mut().zip(&spring.areas) {
                *value += pressure[i] * a;
            }
            for j in 0..n {
                self.equilibrium[i * n + j] -= normalized * root * overlap[j];
            }
        }
        if row
            .iter()
            .chain(&pressure)
            .chain(&self.equilibrium)
            .any(|v| !v.is_finite())
        {
            return Err(bad("caloric mechanical or pressure projection overflow"));
        }
        // Even a mechanical null projection can carry a physical pressure
        // direction. Retain its thermal coordinate so other admissible traces
        // and observers see the original complete constitutive matrix.
        self.rates.push(tau.recip());
        if let Some(mode) = mode {
            self.mode_indices.push(mode);
        }
        self.memory_rows.push(row);
        self.pressure_rows.push(pressure);
        Ok(())
    }
}
