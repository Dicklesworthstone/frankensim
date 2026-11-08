//! Cold-prepared observation of vibration about the complete held-load equilibrium.
//!
//! K = diag(omega^2) + sum k b b^T. For each component's quadrature row w,
//! solve K z = w once with the existing static owner. Then w^T q_eq is
//! z^T (F + sum k r b): a sparse dot product, never a sample-time static solve.
//! The mechanical force, state, energy and reaction equations are not projected.

use super::*;
use fs_math::{c64::C64, dd::Dd};
use std::collections::BTreeSet;

enum ComponentPressure {
    Unshifted,
    Diagonal(Vec<C64>),
    Coupled {
        coefficients: Vec<(usize, f64)>,
        rest: Dd,
    },
}

#[derive(Default)]
pub(crate) struct EquilibriumPressure {
    components: Vec<ComponentPressure>,
    retained_bytes: usize,
    peak_bytes: usize,
    setup_terms: u64,
}

impl EquilibriumPressure {
    #[allow(clippy::too_many_lines)] // one cold layout/admission/static-owner preparation
    pub(crate) fn prepare(
        models: &[ModalAcousticTimeModel],
        columns: &[Vec<f64>],
        links: &[ConnectionLaw],
        config: ModalCouplingConfig,
        budget: Option<PreparedModalBudget>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        if models
            .iter()
            .flat_map(|m| m.modes())
            .all(|m| !has_quadrature(m))
        {
            return Ok(Self::default());
        }
        poll(Some(gate))?;
        let modes: Vec<_> = models.iter().flat_map(|m| m.modes()).copied().collect();
        let n = modes.len();
        let mut labels: Vec<_> = (0..n).collect();
        let mut connected = vec![false; n];
        for (column, law) in columns.iter().zip(links) {
            poll(Some(gate))?;
            if law.stiffness_n_m == 0.0 {
                continue;
            }
            let mut first = None;
            for (i, &value) in column.iter().enumerate() {
                if value == 0.0 {
                    continue;
                }
                connected[i] = true;
                if let Some(j) = first {
                    let (a, b) = (root(&labels, j), root(&labels, i));
                    labels[a.max(b)] = a.min(b);
                } else {
                    first = Some(i);
                }
            }
        }
        for i in 0..n {
            labels[i] = root(&labels, i);
        }
        let mut sizes = vec![0usize; n];
        for &label in &labels {
            sizes[label] += 1;
        }
        let mut required = vec![false; n];
        let mut row_bounds = vec![0usize; models.len()];
        let mut offsets = Vec::with_capacity(models.len() + 1);
        let (mut offset, mut diagonal_entries, mut rows) = (0, 0usize, 0usize);
        for (component, model) in models.iter().enumerate() {
            poll(Some(gate))?;
            offsets.push(offset);
            let range = offset..offset + model.modes().len();
            let quadrature = range.clone().filter(|&i| has_quadrature(&modes[i]));
            if quadrature.clone().any(|i| connected[i]) {
                let roots: BTreeSet<_> = quadrature.map(|i| labels[i]).collect();
                for label in roots {
                    required[label] = true;
                    row_bounds[component] += sizes[label];
                }
                rows += 1;
            } else if quadrature.count() != 0 {
                diagonal_entries += model.modes().len();
            }
            offset = range.end;
        }
        offsets.push(offset);
        let active: Vec<_> = labels.iter().map(|&label| required[label]).collect();
        // Keep only springs in observed stiffness islands. In particular a
        // disconnected floating assembly cannot make an unrelated readout singular.
        let selected: Vec<_> = columns
            .iter()
            .zip(links)
            .enumerate()
            .filter(|(_, (column, law))| {
                law.stiffness_n_m > 0.0
                    && column
                        .iter()
                        .zip(&active)
                        .any(|(&b, &used)| used && b != 0.0)
            })
            .map(|(i, _)| i)
            .collect();
        let free = modes
            .iter()
            .zip(&active)
            .filter(|(mode, used)| **used && mode.angular_frequency_rad_s == 0.0)
            .count();
        let (nn, mm, ff, rr) = (
            n as u128,
            selected.len() as u128,
            free as u128,
            rows as u128,
        );
        let terms = if rows == 0 {
            0
        } else {
            nn*(mm+1)*(mm+1) + mm*mm*mm + ff*mm*mm + ff*ff*mm + ff*ff*ff
                // Each row includes the static solve, original-K residual,
                // and signed rest-load projection, including free balance.
                + rr*(6*nn*(mm+1) + 4*mm*mm + 2*ff*ff + 4*ff*mm)
        };
        let original = nn * (columns.len() as u128 + 1).pow(2);
        if original + terms > config.max_setup_terms as u128 || terms > u64::MAX as u128 {
            return Err(invalid(
                "held-load pressure preparation exceeds max_setup_terms",
            ));
        }
        let retained = row_bounds.iter().map(|&v| v as u128).sum::<u128>()
            * core::mem::size_of::<(usize, f64)>() as u128
            + diagonal_entries as u128 * core::mem::size_of::<C64>() as u128
            + models.len() as u128 * core::mem::size_of::<ComponentPressure>() as u128;
        let peak = retained + 192 * nn + 16 * nn * mm + 32 * mm * mm + 32 * ff * ff + 32 * ff * mm;
        if retained > usize::MAX as u128 || peak > usize::MAX as u128 {
            return Err(invalid("held-load pressure preparation storage overflows"));
        }
        let mut result = Self {
            components: Vec::new(),
            retained_bytes: retained as usize,
            peak_bytes: peak as usize,
            setup_terms: terms as u64,
        };
        if let Some(budget) = budget {
            result.reserve(budget)?;
        }
        let local_columns: Vec<_> = selected.iter().map(|&i| columns[i].clone()).collect();
        let local_links: Vec<_> = selected
            .iter()
            .map(|&i| ConnectionLaw {
                stiffness_n_m: links[i].stiffness_n_m,
                damping_n_s_m: 0.0,
                rest_extension_m: links[i].rest_extension_m,
                constant_loss_reference_omega_rad_s: None,
            })
            .collect();
        let system = StaticSystem {
            models,
            columns: &local_columns,
            connections: &local_links,
            config,
            count: n,
        };
        let response = if rows == 0 {
            None
        } else {
            // Free directions required by an observed island must be resolved
            // by the existing support owner. No pseudoinverse, tether or new
            // rank tolerance is introduced for an unresolved singular island.
            Some(StaticResponse::for_system(
                &system,
                &[],
                Some(&active),
                gate,
            )?)
        };
        result.components = Vec::with_capacity(models.len());
        for (component, model) in models.iter().enumerate() {
            poll(Some(gate))?;
            if row_bounds[component] == 0 {
                result
                    .components
                    .push(if model.modes().iter().any(has_quadrature) {
                        let mut transfers = Vec::with_capacity(model.modes().len());
                        transfers
                            .extend(model.modes().iter().map(|m| m.pressure_per_modal_velocity));
                        ComponentPressure::Diagonal(transfers)
                    } else {
                        ComponentPressure::Unshifted
                    });
                continue;
            }
            let mut weights = vec![0.0; n];
            for (i, mode) in (offsets[component]..offsets[component + 1]).zip(model.modes()) {
                weights[i] =
                    finite(mode.pressure_per_modal_velocity.im * mode.angular_frequency_rad_s)?;
            }
            let row = response
                .as_ref()
                .expect("coupled rows have a static owner")
                .solve_system(&system, &weights, false, false, None, gate)?;
            check_adjoint(&system, &row, &weights, gate)?;
            let count = row.iter().filter(|&&value| value != 0.0).count();
            if count > row_bounds[component] {
                return Err(invalid(
                    "static observation escaped its admitted stiffness islands",
                ));
            }
            let mut coefficients = Vec::with_capacity(count);
            coefficients.extend(
                row.iter()
                    .copied()
                    .enumerate()
                    .filter(|(_, value)| *value != 0.0),
            );
            let mut rest = Dd::ZERO;
            for (column, law) in local_columns.iter().zip(&local_links) {
                let mut projection = Dd::ZERO;
                for (&z, &b) in row.iter().zip(column) {
                    projection = projection + Dd::from_f64(z) * Dd::from_f64(b);
                }
                rest = rest
                    + projection
                        * Dd::from_f64(law.stiffness_n_m)
                        * Dd::from_f64(law.rest_extension_m);
            }
            finite_dd(rest)?;
            result
                .components
                .push(ComponentPressure::Coupled { coefficients, rest });
        }
        Ok(result)
    }

    // Reserve retained observation storage beside the exact operator, and
    // charge both cold preparations to the caller's existing dense-work cap.
    pub(crate) fn reserve(
        &self,
        mut budget: PreparedModalBudget,
    ) -> Result<PreparedModalBudget, ModalCouplingError> {
        if self.peak_bytes > budget.max_dense_bytes || self.setup_terms > budget.max_multiply_terms
        {
            return Err(invalid(
                "held-load pressure preparation exceeds the explicit dense byte/work budget",
            ));
        }
        budget.max_dense_bytes = budget
            .max_dense_bytes
            .checked_sub(self.retained_bytes)
            .ok_or_else(|| invalid("held-load pressure storage exceeds the dense byte budget"))?;
        budget.max_multiply_terms -= self.setup_terms;
        Ok(budget)
    }

    pub(crate) fn needs_correction(&self, component: usize) -> bool {
        matches!(
            self.components.get(component),
            Some(ComponentPressure::Diagonal(_) | ComponentPressure::Coupled { .. })
        )
    }

    pub(crate) fn observe(
        &self,
        component: usize,
        model: &ModalAcousticTimeModel,
        local_forces: &[f64],
        forces: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<f64, ModalCouplingError> {
        match &self.components[component] {
            ComponentPressure::Diagonal(transfers) => model
                .observer_pressure_with_transfers_about_static_equilibrium(transfers, local_forces)
                .map_err(|source| match source {
                    ModalAcousticTimeError::BudgetExceeded { value, limit, .. } => {
                        ModalCouplingError::Budget {
                            what: "component observer pressure",
                            value,
                            limit,
                        }
                    }
                    source => ModalCouplingError::Component { component, source },
                }),
            ComponentPressure::Coupled { coefficients, rest } => {
                let mut pressure = -*rest;
                for (i, (mode, state)) in model.modes().iter().zip(model.states()).enumerate() {
                    if i % 64 == 0 {
                        poll(gate)?;
                    }
                    pressure = pressure
                        + Dd::from_f64(mode.pressure_per_modal_velocity.re)
                            * Dd::from_f64(state.velocity_m_sqrt_kg_per_s);
                    let weight = mode.pressure_per_modal_velocity.im * mode.angular_frequency_rad_s;
                    pressure = pressure
                        + Dd::from_f64(weight) * Dd::from_f64(state.displacement_m_sqrt_kg);
                }
                for (i, &(mode, weight)) in coefficients.iter().enumerate() {
                    if i % 64 == 0 {
                        poll(gate)?;
                    }
                    pressure = pressure - Dd::from_f64(weight) * Dd::from_f64(forces[mode]);
                }
                finite_dd(pressure)
            }
            ComponentPressure::Unshifted => {
                Err(invalid("unshifted pressure has no equilibrium correction"))
            }
        }
    }
}

fn root(labels: &[usize], mut index: usize) -> usize {
    while labels[index] != index {
        index = labels[index];
    }
    index
}

fn has_quadrature(mode: &crate::modal_acoustic_time::ModalAcousticMode) -> bool {
    mode.angular_frequency_rad_s != 0.0 && mode.pressure_per_modal_velocity.im != 0.0
}

fn finite_dd(value: Dd) -> Result<f64, ModalCouplingError> {
    if !value.lo.is_finite() {
        return Err(invalid("nonfinite equilibrium pressure accumulation"));
    }
    finite(value.to_f64())
}

// Check the ORIGINAL K z = w, including the supported free coordinates,
// after the Schur solve. A small reduced residual alone does not admit a map.
fn check_adjoint(
    system: &StaticSystem<'_>,
    row: &[f64],
    weights: &[f64],
    gate: &CancelGate,
) -> Result<(), ModalCouplingError> {
    let mut applied = vec![Dd::ZERO; row.len()];
    let mut scales = vec![0.0; row.len()];
    for (i, (mode, &z)) in system
        .models
        .iter()
        .flat_map(|m| m.modes())
        .zip(row)
        .enumerate()
    {
        let value = Dd::from_f64(mode.angular_frequency_rad_s)
            * Dd::from_f64(mode.angular_frequency_rad_s)
            * Dd::from_f64(z);
        applied[i] = value;
        scales[i] = finite(finite_dd(value)?.abs() + weights[i].abs())?;
    }
    for (column, law) in system.columns.iter().zip(system.connections) {
        poll(Some(gate))?;
        let mut strain = Dd::ZERO;
        for (&z, &b) in row.iter().zip(column) {
            strain = strain + Dd::from_f64(z) * Dd::from_f64(b);
        }
        let force = strain * Dd::from_f64(law.stiffness_n_m);
        for (i, &b) in column.iter().enumerate() {
            let value = force * Dd::from_f64(b);
            applied[i] = applied[i] + value;
            scales[i] = finite(scales[i] + finite_dd(value)?.abs())?;
        }
    }
    let mut relative = 0.0_f64;
    for ((value, &weight), scale) in applied.into_iter().zip(weights).zip(scales) {
        let residual = finite_dd(value - Dd::from_f64(weight))?.abs();
        relative = relative.max(if scale == 0.0 {
            residual
        } else {
            residual / scale
        });
    }
    if relative > system.config.solve_relative_tolerance {
        return Err(ModalCouplingError::SolveResidual {
            relative,
            tolerance: system.config.solve_relative_tolerance,
        });
    }
    Ok(())
}
