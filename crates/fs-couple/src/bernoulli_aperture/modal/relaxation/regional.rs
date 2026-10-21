//! Spatially heterogeneous Maxwell spectra, without diagonal modal averaging.
use super::*;
use crate::bernoulli_aperture::modal::exact::{CoupledArm, ExactModes};
use fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM;

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

/// Positive factors of the complete regional form. Exact zero rows carry no
/// strain energy; an unresolved singular nonzero subspace refuses. No eigenvalue
/// clipping, diagonal replacement or rank cutoff changes the supplied material.
fn factor(
    form: &[f64],
    n: usize,
    gate: &CancelGate,
) -> Result<Vec<Vec<f64>>, AcousticRealizeError> {
    let active: Vec<_> = (0..n)
        .filter(|&i| form[n * i..n * (i + 1)].iter().any(|v| *v != 0.0))
        .collect();
    let size = active.len();
    if size == 0 {
        return Ok(Vec::new());
    }
    let scales: Vec<_> = active.iter().map(|&i| det::sqrt(form[n * i + i])).collect();
    if scales.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err(invalid(
            "regional bending memory needs a resolved positive strain subspace",
        ));
    }
    let normalized: Vec<_> = (0..size * size)
        .map(|i| {
            form[n * active[i / size] + active[i % size]] / scales[i / size] / scales[i % size]
        })
        .collect();
    let root = fs_la::factor::cholesky(&normalized, size).map_err(|_| {
        invalid("regional bending memory cannot discard singular nonzero strain channels")
    })?;
    let mut rows = Vec::with_capacity(size);
    for channel in 0..size {
        poll(gate)?;
        let mut row = vec![0.0; n];
        for (i, &mode) in active.iter().enumerate() {
            row[mode] = scales[i] * root.l(i, channel);
        }
        if row.iter().any(|v| !v.is_finite()) || row[active[channel]] == 0.0 {
            return Err(invalid(
                "regional material strain factor is unrepresentable",
            ));
        }
        rows.push(row);
    }
    for i in 0..n {
        for j in 0..n {
            poll(gate)?;
            let actual: f64 = rows.iter().map(|row| row[i] * row[j]).sum();
            let scale = det::sqrt(form[n * i + i]) * det::sqrt(form[n * j + j]);
            if !actual.is_finite()
                || !scale.is_finite()
                || (actual - form[n * i + j]).abs() > 1e-10 * scale
            {
                return Err(invalid(
                    "regional material channels do not reconstruct the original bending form",
                ));
            }
        }
    }
    Ok(rows)
}

impl ModalAperture {
    /// Attach independent source spectra to every ORIGINAL material region.
    ///
    /// Each selected (modulus [Pa], time [s]) pair is resolved on each region.
    /// Within a region its excess modulus scales that region's original rotated
    /// plane-stress tensor. Regions may have DIFFERENT normalized spectra and
    /// relaxation times. Full DKT rotations and transverse modes from the same
    /// original eigensolve project the resulting off-diagonal bending forms.
    /// This is not a branchwise reduction of a three-dimensional Prony tensor.
    ///
    /// Positive channels H=(c^T q-z)^2/2 retain their own internal strain, exact
    /// held-force propagation and independently integrated dissipation. Pressure,
    /// contact, swept volume and microphone output use the same accepted state.
    /// Equilibrium elastic stiffness and geometric pretension remain unrelaxed.
    ///
    /// # Errors
    /// Cold ExactZoh material plates, 1..=8 selected source terms, no separate
    /// material damping, original source bands and angular caps are required.
    /// Positive relaxation rates require dt/tau<=2. Every state must fit the
    /// existing 256-row exact block; no modes or material channels are dropped.
    /// Singular nonzero regional forms, projection budgets or cancellation refuse.
    pub fn with_regional_material_relaxation(
        mut self,
        properties: &[(String, String)],
        initial: InitialModalMemory,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if self.samples != 0
            || self.relaxation.is_some()
            || self.viscosity.is_some()
            || self.blade_propagator() != BladePropagator::ExactZoh
            || self.material_plates.is_empty()
            || properties.is_empty()
            || properties.len() > 8
            || self.spec.modes.iter().any(|m| m.damping_n_s_m != 0.0)
        {
            return Err(invalid(
                "regional material memory requires one cold exact attachment with original plates, 1..=8 terms and no separate damping",
            ));
        }
        let count = self.state.len();
        let maximum_arms = (MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2)
            .checked_sub(3 * count)
            .ok_or_else(|| {
                invalid("complete regional memory exceeds the exact propagation capacity")
            })?;
        let mut spectra = Vec::with_capacity(self.material_plates.len());
        let mut arms = Vec::<CoupledArm>::new();
        for plate in &self.material_plates {
            poll(gate)?;
            if plate.source().bending_viscosity().is_some() {
                return Err(invalid(
                    "regional material memory cannot duplicate selected bending viscosity",
                ));
            }
            let range = plate.mode_range();
            let n = range.len();
            let low = self.spec.modes[range.clone()]
                .iter()
                .map(|m| det::sqrt(m.stiffness_n_m / m.mass_kg))
                .fold(f64::INFINITY, f64::min);
            let high = self.spec.modes[range.clone()]
                .iter()
                .map(|m| det::sqrt(m.stiffness_n_m / m.mass_kg))
                .fold(0.0, f64::max);
            let mut regional = Vec::with_capacity(plate.source().regions().len());
            for region in plate.source().regions() {
                poll(gate)?;
                let material = &region.input.material;
                let key = match region.input.model {
                    PlateMaterialModel::Isotropic => "young_modulus",
                    _ => "young_modulus_1",
                };
                let equilibrium = property(material, key, Pressure::DIMS, [low, high])?;
                let terms = properties
                    .iter()
                    .map(|(modulus, time)| {
                        Ok((
                            property(material, modulus, Pressure::DIMS, [low, high])?,
                            property(material, time, Dims([0, 0, 1, 0, 0, 0]), [low, high])?,
                        ))
                    })
                    .collect::<Result<Vec<_>, AcousticRealizeError>>()?;
                let law =
                    GeneralizedMaxwell::new(equilibrium, terms).map_err(super::super::nonlinear)?;
                if law.terms.iter().any(|&(e, tau)| {
                    e > 0.0 && (e / equilibrium == 0.0 || self.spec.time_step_s / tau > 2.0)
                }) {
                    return Err(invalid(
                        "regional material relaxation is unrepresentable or unresolved by the original clock",
                    ));
                }
                regional.push(law);
            }
            let (basis, scale) = plate.constitutive_basis();
            let forms = basis.regional_bending_forms(
                plate.source(),
                scale,
                plate.options().max_projection_terms,
                gate,
            )?;
            let first_arm = arms.len();
            for (form, law) in forms.iter().zip(&regional) {
                if law.terms.iter().all(|(e, _)| *e == 0.0) {
                    continue;
                }
                let rows = factor(form, n, gate)?;
                let positive = law.terms.iter().filter(|(e, _)| *e > 0.0).count();
                if rows
                    .len()
                    .checked_mul(positive)
                    .and_then(|a| arms.len().checked_add(a))
                    .is_none_or(|a| a > maximum_arms)
                {
                    return Err(invalid(
                        "complete regional material history exceeds the exact block capacity",
                    ));
                }
                for &(modulus, relaxation_time_s) in &law.terms {
                    if modulus == 0.0 {
                        continue;
                    }
                    let scale = det::sqrt(modulus / law.e_inf);
                    for row in &rows {
                        poll(gate)?;
                        let mut coupling = vec![0.0; count];
                        for (i, &value) in row.iter().enumerate() {
                            coupling[range.start + i] = scale * value;
                        }
                        if coupling.iter().any(|v| !v.is_finite())
                            || coupling.iter().all(|v| *v == 0.0)
                        {
                            return Err(invalid(
                                "regional projection lost a positive material channel",
                            ));
                        }
                        arms.push(CoupledArm {
                            coupling,
                            relaxation_time_s,
                        });
                    }
                }
            }
            // Symmetric mass-normalized instantaneous tangent. Its maximum
            // absolute row sum bounds every retained frequency, including all
            // cross-mode material reactions; no mean material is substituted.
            let mut tangent = vec![0.0; n * n];
            for i in 0..n {
                let mode = &self.spec.modes[range.start + i];
                tangent[n * i + i] = mode.stiffness_n_m / mode.mass_kg;
            }
            for arm in &arms[first_arm..] {
                for i in 0..n {
                    for j in 0..n {
                        tangent[n * i + j] += arm.coupling[range.start + i]
                            / det::sqrt(self.spec.modes[range.start + i].mass_kg)
                            * arm.coupling[range.start + j]
                            / det::sqrt(self.spec.modes[range.start + j].mass_kg);
                    }
                }
            }
            let instantaneous = det::sqrt(
                tangent
                    .chunks_exact(n)
                    .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
                    .fold(0.0, f64::max),
            );
            if tangent.iter().any(|v| !v.is_finite())
                || !instantaneous.is_finite()
                || instantaneous * self.spec.time_step_s > plate.options().maximum_angular_step
            {
                return Err(invalid(
                    "regional instantaneous material stiffness exceeds its original clock cap",
                ));
            }
            for region in plate.source().regions() {
                poll(gate)?;
                let material = &region.input.material;
                for (modulus, time) in properties {
                    property(material, modulus, Pressure::DIMS, [low, instantaneous])?;
                    property(
                        material,
                        time,
                        Dims([0, 0, 1, 0, 0, 0]),
                        [low, instantaneous],
                    )?;
                }
                let elastic: &[(&str, Dims)] = match region.input.model {
                    PlateMaterialModel::Isotropic => &[
                        ("density", Density::DIMS),
                        ("young_modulus", Pressure::DIMS),
                        ("poisson_ratio", Dims::NONE),
                    ],
                    _ => &[
                        ("density", Density::DIMS),
                        ("young_modulus_1", Pressure::DIMS),
                        ("young_modulus_2", Pressure::DIMS),
                        ("shear_modulus_12", Pressure::DIMS),
                        ("poisson_ratio_12", Dims::NONE),
                    ],
                };
                for &(key, dims) in elastic {
                    property(material, key, dims, [low, instantaneous])?;
                }
            }
            spectra.push(regional);
        }
        let states = vec![
            arms.iter()
                .map(|arm| match initial {
                    InitialModalMemory::Relaxed => arm
                        .coupling
                        .iter()
                        .zip(&self.state)
                        .map(|(c, state)| c * state.displacement_m)
                        .sum(),
                    InitialModalMemory::Unrelaxed => 0.0,
                })
                .collect::<Vec<f64>>(),
        ];
        if states[0].iter().any(|z| !z.is_finite()) {
            return Err(invalid(
                "initial regional material history is unrepresentable",
            ));
        }
        let exact = ExactModes::with_coupled_arms(&self.spec, arms.clone())?;
        let relaxation = ModalRelaxation {
            spectra,
            initial,
            arms: Vec::new(),
            states,
            coupled_arms: Some(arms),
        };
        let energy = self.energy(&self.state) + relaxation.energy(&self.state, &relaxation.states);
        if !energy.is_finite() || energy > self.spec.limits.maximum_energy_j {
            return Err(invalid(
                "initial regional material history exceeds the original energy budget",
            ));
        }
        poll(gate)?;
        self.exact = Some(exact);
        self.relaxation = Some(relaxation);
        Ok(self)
    }
}
