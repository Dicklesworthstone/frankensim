//! Physical joint norm and accepted-only accounting; no surrogate field update.

use super::*;
use fs_flux::bdm::eval_basis;
use fs_math::dd::Dd;

#[allow(clippy::too_many_arguments)]
pub(super) fn compare<G, T, E>(
    cx: &Cx<'_>,
    system: &FluxSystem<'_>,
    old: &DepositionSimulation<'_>,
    old_velocity: &[f64],
    coarse: &LaminarEnthalpyCandidate<'_, G, T>,
    first: &LaminarEnthalpyCandidate<'_, G, T>,
    fine: &LaminarEnthalpyCandidate<'_, G, T>,
    correction: f64,
    controls: LaminarAdaptiveControls,
) -> Result<LaminarTimeEstimate, LaminarEnthalpyError<E>> {
    let mut result = LaminarTimeEstimate {
        start_s: old.time_s(),
        end_s: fine.thermal().time_s(),
        energy_error_ratio: 0.0,
        temperature_error_ratio: 0.0,
        velocity_error_ratio: 0.0,
        transfer_error_ratio: 0.0,
        maximum_energy_difference_j: 0.0,
        maximum_temperature_difference_k: 0.0,
        maximum_velocity_difference_m_s: 0.0,
        maximum_transfer_difference_j: 0.0,
    };
    if old.cells().len() != coarse.thermal().cells().len()
        || old.cells().len() != fine.thermal().cells().len()
    {
        return Err(LaminarEnthalpyError::Invalid(
            "temporal comparison changed original inventory",
        ));
    }
    for (&id, previous) in old.cells() {
        checkpoint(cx)?;
        let a = coarse
            .thermal()
            .cells()
            .get(&id)
            .ok_or(LaminarEnthalpyError::Invalid("missing coarse thermal cell"))?;
        let b = fine
            .thermal()
            .cells()
            .get(&id)
            .ok_or(LaminarEnthalpyError::Invalid("missing fine thermal cell"))?;
        if a.mass_kg != previous.mass_kg
            || b.mass_kg != previous.mass_kg
            || a.material_id != previous.material_id
            || b.material_id != previous.material_id
        {
            return Err(LaminarEnthalpyError::Invalid(
                "temporal comparison changed original mass or law",
            ));
        }
        // Subtract original specific energies BEFORE multiplying by mass:
        // formation-energy offsets never become temporal-error allowances.
        let mass = Dd::from_f64(previous.mass_kg);
        let h0 = Dd::from_f64(previous.specific_h_j_kg);
        let ha = Dd::from_f64(a.specific_h_j_kg);
        let hb = Dd::from_f64(b.specific_h_j_kg);
        let difference = ((hb - ha) * mass).to_f64().abs();
        let change = ((ha - h0) * mass)
            .to_f64()
            .abs()
            .max(((hb - h0) * mass).to_f64().abs());
        result.maximum_energy_difference_j = result.maximum_energy_difference_j.max(difference);
        result.energy_error_ratio = result.energy_error_ratio.max(ratio(
            difference,
            change,
            controls.energy_absolute_tolerance_j,
            controls.relative_tolerance,
            correction,
        )?);
        let t0 = old.temperature_k(id)?;
        let ta = coarse.thermal().temperature_k(id)?;
        let tb = fine.thermal().temperature_k(id)?;
        let difference = (tb - ta).abs();
        let change = (ta - t0).abs().max((tb - t0).abs());
        result.maximum_temperature_difference_k =
            result.maximum_temperature_difference_k.max(difference);
        result.temperature_error_ratio = result.temperature_error_ratio.max(ratio(
            difference,
            change,
            controls.temperature_absolute_tolerance_k,
            controls.relative_tolerance,
            correction,
        )?);
    }
    let a = coarse.flow().step.solution.velocity_dofs();
    let b = fine.flow().step.solution.velocity_dofs();
    if [old_velocity.len(), a.len(), b.len()]
        .iter()
        .any(|&n| n != system.n_u)
    {
        return Err(LaminarEnthalpyError::Invalid(
            "incomplete original BDM velocity",
        ));
    }
    // BDM1 is affine on each original triangle. Its Euclidean error norm
    // reaches its maximum at a vertex; centroid-only comparison would miss
    // shear modes. Original global basis orientations supply all edge signs.
    for (cell, vertices) in system.mesh.tris.iter().enumerate() {
        checkpoint(cx)?;
        for &vertex in vertices {
            let point = system.mesh.verts[vertex];
            let difference = physical_difference(system, cell, point, b, a)?;
            let change = physical_difference(system, cell, point, a, old_velocity)?
                .max(physical_difference(system, cell, point, b, old_velocity)?);
            result.maximum_velocity_difference_m_s =
                result.maximum_velocity_difference_m_s.max(difference);
            result.velocity_error_ratio = result.velocity_error_ratio.max(ratio(
                difference,
                change,
                controls.velocity_absolute_tolerance_m_s,
                controls.relative_tolerance,
                correction,
            )?);
        }
    }
    for (a, (b, c)) in terms(coarse.report())
        .into_iter()
        .zip(terms(first.report()).into_iter().zip(terms(fine.report())))
    {
        let fine_total = Dd::from_f64(b) + Dd::from_f64(c);
        let difference = (fine_total - Dd::from_f64(a)).to_f64().abs();
        let change = a.abs().max(fine_total.to_f64().abs());
        result.maximum_transfer_difference_j = result.maximum_transfer_difference_j.max(difference);
        result.transfer_error_ratio = result.transfer_error_ratio.max(ratio(
            difference,
            change,
            controls.transfer_absolute_tolerance_j,
            controls.relative_tolerance,
            correction,
        )?);
    }
    Ok(result)
}

fn physical_difference<E>(
    system: &FluxSystem<'_>,
    cell: usize,
    point: [f64; 2],
    a: &[f64],
    b: &[f64],
) -> Result<f64, LaminarEnthalpyError<E>> {
    let mut difference = [Dd::from_f64(0.0); 2];
    for i in 0..6 {
        let (edge, _) = system.mesh.tri_edges[cell][i / 2];
        let index = 2 * edge + i % 2;
        let value = eval_basis(&system.bases[cell], i, point);
        let dof = Dd::from_f64(a[index]) - Dd::from_f64(b[index]);
        for axis in 0..2 {
            difference[axis] = difference[axis] + dof * Dd::from_f64(value[axis]);
        }
    }
    let norm = difference[0].to_f64().hypot(difference[1].to_f64());
    if norm.is_finite() {
        Ok(norm)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}

fn ratio<E>(
    difference: f64,
    change: f64,
    absolute: f64,
    relative: f64,
    correction: f64,
) -> Result<f64, LaminarEnthalpyError<E>> {
    if !difference.is_finite() || !change.is_finite() {
        return Err(DepositionError::Arithmetic.into());
    }
    let scale = absolute + relative * change;
    if !scale.is_finite() || scale <= 0.0 {
        return Err(DepositionError::Arithmetic.into());
    }
    let value = (difference / scale) * correction;
    if value.is_nan() {
        Err(DepositionError::Arithmetic.into())
    } else {
        // Positive infinity requests refinement, never apparent convergence.
        Ok(value)
    }
}

fn terms<T>(report: &LaminarEnthalpyReport<T>) -> [f64; 5] {
    [
        report.bulk_viscous_heat_j,
        report.slip_heat_j,
        report.external_stream_enthalpy_j,
        report.other_external_heat_j,
        report.buoyancy.as_ref().map_or(0.0, |b| b.endpoint_work_j),
    ]
}

pub(super) struct Account {
    terms: [Dd; 5],
}
impl Default for Account {
    fn default() -> Self {
        Self {
            terms: [Dd::from_f64(0.0); 5],
        }
    }
}
impl Account {
    pub(super) fn add<T, E>(
        &mut self,
        report: &LaminarEnthalpyReport<T>,
    ) -> Result<(), LaminarEnthalpyError<E>> {
        for (total, value) in self.terms.iter_mut().zip(terms(report)) {
            if !value.is_finite() {
                return Err(DepositionError::Arithmetic.into());
            }
            *total = *total + Dd::from_f64(value);
            if !total.to_f64().is_finite() {
                return Err(DepositionError::Arithmetic.into());
            }
        }
        Ok(())
    }

    pub(super) fn finish<E>(
        self,
        cx: &Cx<'_>,
        initial: &DepositionSimulation<'_>,
        final_state: &DepositionSimulation<'_>,
        allowance_j: f64,
    ) -> Result<[f64; 7], LaminarEnthalpyError<E>> {
        if initial.cells().len() != final_state.cells().len() {
            return Err(LaminarEnthalpyError::Invalid(
                "window changed original inventory",
            ));
        }
        let mut stored = Dd::from_f64(0.0);
        for (&id, old) in initial.cells() {
            checkpoint(cx)?;
            let new = final_state
                .cells()
                .get(&id)
                .ok_or(LaminarEnthalpyError::Invalid("missing final original cell"))?;
            if new.mass_kg != old.mass_kg || new.material_id != old.material_id {
                return Err(LaminarEnthalpyError::Invalid(
                    "window changed original mass or law",
                ));
            }
            stored = stored
                + Dd::from_f64(old.mass_kg)
                    * (Dd::from_f64(new.specific_h_j_kg) - Dd::from_f64(old.specific_h_j_kg));
        }
        // Measure original storage independently; never reconstruct it by
        // summing already-accepted residuals or repairing external sources.
        let heat = self.terms[..4]
            .iter()
            .fold(Dd::from_f64(0.0), |sum, &x| sum + x);
        let residual = (stored - heat).to_f64();
        let scale = stored.to_f64().abs()
            + self.terms[..4]
                .iter()
                .map(|x| x.to_f64().abs())
                .sum::<f64>();
        let tolerance = allowance_j + 64.0 * f64::EPSILON * scale;
        if !residual.is_finite() || !tolerance.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        if residual.abs() > tolerance {
            return Err(DepositionError::EnergyBalance(residual).into());
        }
        Ok([
            stored.to_f64(),
            self.terms[0].to_f64(),
            self.terms[1].to_f64(),
            self.terms[2].to_f64(),
            self.terms[3].to_f64(),
            self.terms[4].to_f64(),
            residual,
        ])
    }
}
