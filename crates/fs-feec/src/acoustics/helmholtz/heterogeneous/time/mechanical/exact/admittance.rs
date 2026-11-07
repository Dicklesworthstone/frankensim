//! Attach complete original admittance memories to the fluid/mechanical operator.

use super::*;

pub(super) const BASE_FORMS: usize = 4;
pub(super) const HEAT: usize = 4;
pub(super) const WORK: usize = 5;
pub(super) const FORM_COUNT: usize = 6;

pub(super) struct ExactAdmittance {
    pub(super) trace: TraceMemory,
    pub(super) offset: usize,
}

pub(super) struct Extension {
    pub(super) generator: Vec<f64>,
    pub(super) scales: Vec<f64>,
    pub(super) forms: Vec<Vec<f64>>,
    pub(super) observations: Vec<Vec<f64>>,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn extend(
    model: &HeterogeneousHelmholtzFe<'_, '_>,
    duration_s: f64,
    release: &PressureReleaseConstraints,
    fluid_mass: &MassMap,
    pressure_map: &[f64],
    arms: &[ExactAdmittance],
    base: usize,
    states: usize,
    q: &[f64],
    j: &[f64],
    r: &[f64],
    physical_forms: &[Vec<f64>],
    source_observations: &[Vec<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Extension, HeterogeneousError> {
    let n = release.free_dofs.len();
    let mut generator = vec![0.0; states * states];
    let mut storage = vec![0.0; states * states];
    let mut forms = vec![vec![0.0; states * states]; FORM_COUNT];
    // Lift the complete mechanical backbone without changing its coordinates,
    // Hamiltonian, mechanical reaction, Maxwell spectrum or source columns.
    for row in 0..base {
        poll(cancelled)?;
        for col in 0..base {
            storage[row * states + col] = q[row * base + col];
            let mut value = Dd::from_f64(0.0);
            for k in 0..base {
                value = value
                    + (Dd::from_f64(j[row * base + k]) - Dd::from_f64(r[row * base + k]))
                        * Dd::from_f64(q[k * base + col]);
            }
            generator[row * states + col] = value.to_f64();
            for (form, original) in forms[..BASE_FORMS].iter_mut().zip(physical_forms) {
                form[row * states + col] = original[row * base + col];
            }
        }
    }
    let mut scales = Vec::with_capacity(states - base);
    for arm in arms {
        let trace = &arm.trace;
        let m = trace.dofs.len();
        let mut load = vec![0.0; n * m];
        for (row, &global) in trace.dofs.iter().enumerate() {
            poll(cancelled)?;
            let free = release.original_to_free[global];
            // This row is retained even if pressure is constrained or Cj is zero.
            generator[(arm.offset + row) * states + arm.offset + row] = -1.0 / trace.tau;
            if free != usize::MAX {
                for col in 0..n {
                    generator[(arm.offset + row) * states + n + col] =
                        pressure_map[free * n + col] / trace.tau;
                }
            }
            let (columns, values) = trace.conductance.row(row);
            for (&col, &value) in columns.iter().zip(values) {
                forms[HEAT][(arm.offset + row) * states + arm.offset + col] = value;
                storage[(arm.offset + row) * states + arm.offset + col] = trace.tau * value;
                if free != usize::MAX {
                    generator[(n + free) * states + arm.offset + col] =
                        -value / fluid_mass.scales[free];
                    load[free * m + col] = value;
                }
            }
            // Scaling is numerical only. Zero residue does not confer energy,
            // remove a memory, or permit discarding its pressure-driven dynamics.
            let diagonal = columns.binary_search(&row).ok().map_or(0.0, |i| values[i]);
            let scale = if diagonal > 0.0 {
                det::sqrt(trace.tau * diagonal)
            } else {
                let (mass_columns, mass_values) = model.mass.row(global);
                let at = mass_columns
                    .binary_search(&global)
                    .map_err(|_| invalid("original admittance trace lacks fluid mass"))?;
                let scale = trace.tau.max(duration_s);
                let ratio = (trace.tau / scale) / (trace.tau / scale + duration_s / scale);
                det::sqrt(mass_values[at]) * det::sqrt(ratio)
            };
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("exact composed admittance scale is unrepresentable").into());
            }
            scales.push(scale);
        }
        // Independent original pressure-work pairing. The symmetric half blocks
        // integrate signed p^T Cj q, not a heat or an endpoint energy deficit.
        for row in 0..n {
            poll(cancelled)?;
            for col in 0..m {
                let mut value = Dd::from_f64(0.0);
                for k in 0..n {
                    value = value
                        + Dd::from_f64(pressure_map[k * n + row]) * Dd::from_f64(load[k * m + col]);
                }
                let half = 0.5 * value.to_f64();
                forms[WORK][(n + row) * states + arm.offset + col] = half;
                forms[WORK][(arm.offset + col) * states + n + row] = half;
            }
        }
    }
    admit_physical_identity(&storage, &generator, &forms, states, cancelled)?;
    let mut observations = Vec::with_capacity(source_observations.len());
    for original in source_observations {
        poll(cancelled)?;
        let mut observation = vec![0.0; states];
        observation[..base].copy_from_slice(original);
        observations.push(observation);
    }
    Ok(Extension {
        generator,
        scales,
        forms,
        observations,
    })
}

// Audit original storage and all physical heat against the complete generator.
// Both reversible exchanges are internal and therefore excluded from this sum.
fn admit_physical_identity(
    storage: &[f64],
    generator: &[f64],
    forms: &[Vec<f64>],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), HeterogeneousError> {
    for row in 0..n {
        poll(cancelled)?;
        for col in 0..n {
            let mut sum = Dd::from_f64(0.0);
            let mut scale = 0.0;
            for k in 0..n {
                for term in [
                    Dd::from_f64(storage[row * n + k]) * Dd::from_f64(generator[k * n + col]),
                    Dd::from_f64(generator[k * n + row]) * Dd::from_f64(storage[k * n + col]),
                ] {
                    sum = sum + term;
                    scale += term.to_f64().abs();
                }
            }
            for heat in [0, 1, 2, HEAT] {
                let term = 2.0 * forms[heat][row * n + col];
                sum = sum + Dd::from_f64(term);
                scale += term.abs();
            }
            if !sum.to_f64().is_finite()
                || !scale.is_finite()
                || sum.to_f64().abs() > 2048.0 * f64::EPSILON * n as f64 * scale
            {
                return Err(invalid("complete mechanical/admittance generator violates original storage and heat identity").into());
            }
        }
    }
    Ok(())
}
