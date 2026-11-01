//! Monolithic vibrations and reciprocal physical ports about an accepted body.
//!
//! The loaded tangent and reference mass enter the original fs-modal slicer
//! unchanged. No unloaded geometry or element stiffness is reconstructed here.
//! Mode shapes and forces use the same junction transformation as equilibrium.

use super::PreloadedShellAssembly;
use fs_exec::Cx;
use fs_modal::{ModePair, MultiSliceReport, SliceOptions, WindowPlan, slice_windows};

/// Output/projection ceilings in addition to the existing eigensolver budgets.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyVibrationBudget {
    /// Largest accepted total number of modes across all requested windows.
    pub maximum_modes: usize,
    /// Largest dense glue-damping projection, in scalar entries (modes squared).
    /// This bounds the actual complete projection, including its off-diagonals.
    pub maximum_damping_entries: usize,
}

/// Certified windows tied by a borrow to their unchanged physical equilibrium.
#[derive(Debug)]
pub struct AssemblyVibration<'a> {
    equilibrium: &'a PreloadedShellAssembly,
    spectrum: MultiSliceReport,
    glue_damping: Vec<f64>,
    modes: usize,
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), String> {
    cx.checkpoint()
        .map_err(|_| "assembled shell vibration cancelled".to_string())
}

impl AssemblyVibration<'_> {
    /// Original per-window inertia counts, residual intervals and work receipts.
    pub fn spectrum(&self) -> &MultiSliceReport {
        &self.spectrum
    }

    /// Modes in ascending window/eigenvalue order, without copying their vectors.
    pub fn modes(&self) -> impl Iterator<Item = &ModePair> {
        self.spectrum.reports.iter().flat_map(|r| r.modes.iter())
    }

    /// Total number of accepted modes in that same order.
    pub fn mode_count(&self) -> usize {
        self.modes
    }

    /// Full row-major Phi^T C_glue Phi [1/s], including all cross-mode terms.
    /// This is a projection onto the declared windows, not a certificate that
    /// omitted modes decouple or that a diagonal damping approximation is valid.
    /// Material, radiation and other loss owners are not inferred or duplicated.
    pub fn glue_damping(&self) -> &[f64] {
        &self.glue_damping
    }

    /// Expand one mass-normalized mode into full physical component coordinates.
    /// Essential supports stay zero and eliminated slave motion comes from T.
    /// The original junction expansion is an indivisible cancellation unit.
    pub fn expand_mode(&self, mode: usize, cx: &Cx<'_>) -> Result<Vec<Vec<f64>>, String> {
        checkpoint(cx)?;
        let mode = self
            .modes()
            .nth(mode)
            .ok_or("unknown assembled vibration mode")?;
        self.expand(&mode.phi, cx)
    }

    /// Project incremental full-component forces/moments by Phi^T T^T.
    /// These are perturbing loads: the static dead forces and connection preload
    /// are already represented by equilibrium and its tangent, not applied again.
    /// Forces on eliminated slave nodes retain their exact work-conjugate image.
    /// Polls per mode; the original junction load projection is indivisible.
    pub fn project_loads(
        &self,
        component_loads: &[&[f64]],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, String> {
        checkpoint(cx)?;
        let retained = self
            .equilibrium
            .assembly
            .project_loads(component_loads)
            .map_err(|e| e.to_string())?;
        let mut loads = Vec::with_capacity(self.modes);
        for mode in self.modes() {
            checkpoint(cx)?;
            let load: f64 = mode.phi.iter().zip(&retained).map(|(a, b)| a * b).sum();
            if !load.is_finite() {
                return Err("assembled modal force projection overflow".into());
            }
            loads.push(load);
        }
        checkpoint(cx)?;
        Ok(loads)
    }

    /// Full component displacement/rotation PERTURBATIONS T Phi q.
    /// The same linear map gives velocity from modal velocity and is conjugate
    /// to `project_loads`. No constant equilibrium offset enters a velocity.
    pub fn perturbation_displacements(
        &self,
        q: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, String> {
        let retained = self.retained_perturbation(q, cx)?;
        self.expand(&retained, cx)
    }

    /// Total physical displacement from the original unloaded component meshes:
    /// T (u_equilibrium + Phi q). Zero q therefore retains the loaded shape.
    /// This is a linear perturbation image; its validity requires the original
    /// moderate-deformation range and an unchanged active connection branch.
    /// It is not a nonlinear update, time integrator or truncation certificate.
    pub fn displacements(&self, q: &[f64], cx: &Cx<'_>) -> Result<Vec<Vec<f64>>, String> {
        let mut retained = self.retained_perturbation(q, cx)?;
        for (value, base) in retained.iter_mut().zip(&self.equilibrium.coordinates) {
            *value += base;
        }
        self.expand(&retained, cx)
    }

    fn expand(&self, retained: &[f64], cx: &Cx<'_>) -> Result<Vec<Vec<f64>>, String> {
        checkpoint(cx)?;
        let expanded = self
            .equilibrium
            .assembly
            .expand_components(retained)
            .map_err(|e| e.to_string())?;
        checkpoint(cx)?;
        Ok(expanded)
    }

    fn retained_perturbation(&self, q: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, String> {
        checkpoint(cx)?;
        if q.len() != self.modes || q.iter().any(|v| !v.is_finite()) {
            return Err(
                "assembled modal coordinates must be finite and match the window basis".into(),
            );
        }
        let mut retained = vec![0.; self.equilibrium.assembly.model.free];
        for (mode, amplitude) in self.modes().zip(q) {
            checkpoint(cx)?;
            for (value, shape) in retained.iter_mut().zip(&mode.phi) {
                *value += amplitude * shape;
            }
        }
        checkpoint(cx)?;
        Ok(retained)
    }
}

/// Slice the exact accepted assembled tangent, then prepare its physical ports.
///
/// Uses fs-modal's shared symbolic pencil, explicit ordering, scoped parallel
/// windows, mass normalization, inertia counts and original-pencil residuals.
/// The equilibrium is borrowed so the mean shape and accepted tangent cannot be
/// changed while these modes are in use. Tether/strut/glue stiffness is already
/// in that tangent and is never added a second time. Mass is unchanged.
///
/// Output caps are checked after the slicer's own bounded computation and before
/// allocating the dense glue projection. Cancellation is polled around slicing,
/// between projected columns and immediately before publication. Empty windows
/// remain valid empty reports. A nonpositive or zero-crossing mode interval is
/// refused as an unstable/unresolved vibration, never clipped to a frequency.
///
/// This consumes an accepted equilibrium, not a separately mutated pencil, and
/// adds no fresh static convergence claim. Window projection retains every glue
/// damping entry; omitted spectral response and loss coupling need the existing
/// owners' separate convergence/certification work before a truncated consumer
/// can claim full-body accuracy.
pub fn eigen_windows<'a>(
    equilibrium: &'a PreloadedShellAssembly,
    plan: &WindowPlan,
    options: &SliceOptions,
    budget: AssemblyVibrationBudget,
    cx: &Cx<'_>,
) -> Result<AssemblyVibration<'a>, String> {
    checkpoint(cx)?;
    let model = &equilibrium.assembly.model;
    if budget.maximum_modes == 0 || budget.maximum_damping_entries == 0 {
        return Err("assembled vibration requires positive mode and projection budgets".into());
    }
    if equilibrium.coordinates.len() != model.free
        || equilibrium.coordinates.iter().any(|v| !v.is_finite())
    {
        return Err("assembled vibration equilibrium coordinates disagree with the tangent".into());
    }
    let spectrum =
        slice_windows(&model.k, &model.m, plan, options, cx).map_err(|e| e.to_string())?;
    checkpoint(cx)?;
    let count = spectrum.reports.iter().try_fold(0usize, |count, report| {
        count
            .checked_add(report.modes.len())
            .ok_or("assembled vibration mode count overflow")
    })?;
    if count > budget.maximum_modes {
        return Err(format!(
            "assembled vibration mode cap exceeded: {count} > {}",
            budget.maximum_modes
        ));
    }
    let entries = count
        .checked_mul(count)
        .ok_or("assembled damping projection size overflow")?;
    if entries > budget.maximum_damping_entries {
        return Err("assembled vibration damping projection budget exceeded".into());
    }
    let modes: Vec<_> = spectrum.reports.iter().flat_map(|r| &r.modes).collect();
    if modes
        .iter()
        .any(|m| !m.lambda.is_finite() || m.lambda <= 0. || m.interval.0 <= 0.)
    {
        return Err(
            "assembled vibration has an unstable or unresolved nonpositive mode interval".into(),
        );
    }
    let mut glue_damping = vec![0.; entries];
    if let Some(glue) = &equilibrium.assembly.glue_damping {
        if glue.nrows() != model.free || glue.ncols() != model.free {
            return Err("assembled glue damping and tangent dimensions disagree".into());
        }
        let mut product = vec![0.; model.free];
        // Same Phi^T C Phi contraction as the junction owner, using borrowed
        // mode vectors so no second full eigenbasis allocation is needed.
        for (j, mode) in modes.iter().enumerate() {
            checkpoint(cx)?;
            glue.spmv(&mode.phi, &mut product);
            for i in 0..=j {
                let value: f64 = modes[i].phi.iter().zip(&product).map(|(a, b)| a * b).sum();
                if !value.is_finite() {
                    return Err("assembled glue damping projection overflow".into());
                }
                glue_damping[i * count + j] = value;
                glue_damping[j * count + i] = value;
            }
        }
    }
    checkpoint(cx)?;
    Ok(AssemblyVibration {
        equilibrium,
        spectrum,
        glue_damping,
        modes: count,
    })
}

#[cfg(test)]
mod tests;
