//! Exact material memory in the original enriched mechanical coordinates.

use super::*;
use crate::shell::mitc::relaxation::ShearShellRelaxationBranch;
use fs_material::visco::PronyRamp;

/// Bounds on the additional original-facet constitutive payload and assembly.
#[derive(Clone, Copy, Debug)]
pub struct MitcRelaxationBudget {
    /// Total supplied facet branches, including explicitly zero increments.
    pub max_element_branches: usize,
    /// Original MITC assembly work, charged once per supplied facet branch.
    pub max_assembly_terms: usize,
    /// Branch operators, prepared ramps and accepted/proposed local histories.
    /// Existing shell geometry, mechanical states and factors have their own bounds.
    pub max_payload_bytes: usize,
}

#[derive(Clone, Debug)]
struct Branch {
    element: usize,
    source: ShearShellRelaxationBranch,
    stiffness: [[f64; 20]; 20],
}

#[derive(Clone, Debug)]
pub(super) struct EnrichedRelaxation {
    branches: Vec<Branch>,
}

#[derive(Clone, Debug)]
pub(super) struct History {
    source: Arc<EnrichedRelaxation>,
    pub(super) coordinates: Vec<[f64; 20]>,
    pub(super) stored_energy_j: f64,
}

impl EnrichedMitcShell {
    /// Attach physical excess stiffness spectra above the EQUILIBRIUM shell
    /// sections. Every supplied A/B/D/shear tensor retains its own time and
    /// original facet axes; neither drilling nor inertia acquires memory.
    /// The original nodal and two interior director coordinates remain fixed.
    ///
    /// Histories use twenty local coordinates per facet branch, independently
    /// of the global mesh size. They must be explicitly initialized through
    /// [`Self::with_relaxation_history`] before advancing this owner. Material
    /// refresh keeps this spectrum frozen; changing a spectrum is not an
    /// implicit transfer of accepted material memory or stored energy.
    pub fn with_relaxation(
        &self,
        sections: &[Vec<ShearShellRelaxationBranch>],
        budget: MitcRelaxationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        check(&mut cancelled)?;
        if sections.len() != self.reference_fields.len() {
            return Err(invalid(
                "one explicit relaxation spectrum per original shell facet required",
            ));
        }
        let mut count = 0usize;
        for section in sections {
            check(&mut cancelled)?;
            count = count
                .checked_add(section.len())
                .ok_or_else(|| invalid("shell relaxation branch extent overflows"))?;
        }
        let bytes_per_branch = core::mem::size_of::<Branch>()
            + core::mem::size_of::<PronyRamp>()
            + 2 * core::mem::size_of::<[f64; 20]>();
        if count > budget.max_element_branches
            || count
                .checked_mul(bytes_per_branch)
                .is_none_or(|n| n > budget.max_payload_bytes)
            || crate::shell::mitc::viscosity::projection_terms(0, count, 0)
                .is_none_or(|n| n > budget.max_assembly_terms)
        {
            return Err(invalid(
                "enriched shell relaxation exceeds its branch, work or payload allowance",
            ));
        }
        let mut branches = Vec::with_capacity(count);
        for (element, section) in sections.iter().enumerate() {
            for source in section {
                check(&mut cancelled)?;
                let stiffness = self.reference_fields[element]
                    .complete_viscosity(source.excess(), &mut cancelled)?;
                branches.push(Branch {
                    element,
                    source: source.clone(),
                    stiffness,
                });
            }
        }
        let mut out = self.clone();
        out.relaxation = if branches.is_empty() {
            None
        } else {
            Some(Arc::new(EnrichedRelaxation { branches }))
        };
        check(&mut cancelled)?;
        Ok(out)
    }

    /// Number of original facet branches, in facet-major/source order.
    pub fn relaxation_branch_count(&self) -> usize {
        self.relaxation
            .as_ref()
            .map_or(0, |source| source.branches.len())
    }

    /// Initialize explicit elastic branch coordinates on a compatible mechanical
    /// state. Each entry contains the facet's eighteen GLOBAL nodal coordinates
    /// followed by its two residual director coordinates, in original triangle
    /// order. Zero means relaxed history, even at a nonzero mechanical shape;
    /// copying that shape means an instantaneous change from the zero reference.
    /// This initial-condition operation does not infer unobserved prehistory.
    pub fn with_relaxation_history(
        &self,
        initial: &MitcTransientState,
        coordinates: &[[f64; 20]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MitcTransientState, PlateError> {
        check(&mut cancelled)?;
        self.admit_state(initial)?;
        if coordinates.len() != self.relaxation_branch_count() {
            return Err(invalid(
                "one finite local elastic history per original shell branch required",
            ));
        }
        for q in coordinates {
            check(&mut cancelled)?;
            if q.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "one finite local elastic history per original shell branch required",
                ));
            }
        }
        let mut out = initial.clone();
        out.relaxation = if let Some(source) = &self.relaxation {
            let mut stored = Dd::ZERO;
            for (branch, q) in source.branches.iter().zip(coordinates) {
                check(&mut cancelled)?;
                let nodal = core::array::from_fn(|i| q[i]);
                let energy = 0.5
                    * self.reference_fields[branch.element].viscous_power(
                        &nodal,
                        [q[18], q[19]],
                        branch.source.factor(),
                        &mut cancelled,
                    )?;
                stored = stored + Dd::from_f64(energy);
            }
            Some(History {
                source: Arc::clone(source),
                coordinates: coordinates.to_vec(),
                stored_energy_j: finite(stored.to_f64())?,
            })
        } else {
            None
        };
        check(&mut cancelled)?;
        Ok(out)
    }

    pub(super) fn admit_relaxation(
        &self,
        before: &Self,
        initial: &MitcTransientState,
    ) -> Result<(), PlateError> {
        match (&self.relaxation, &before.relaxation, &initial.relaxation) {
            (None, None, None) => Ok(()),
            (Some(a), Some(b), Some(history))
                if Arc::ptr_eq(a, b)
                    && Arc::ptr_eq(a, &history.source)
                    && history.coordinates.len() == a.branches.len() =>
            {
                Ok(())
            }
            _ => Err(invalid(
                "midpoint interval requires the same frozen shell spectrum and its initialized histories",
            )),
        }
    }

    pub(super) fn relaxation_tangent(
        &self,
        dt_s: f64,
        effective: &mut Coo,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<PronyRamp>, PlateError> {
        let Some(source) = &self.relaxation else {
            return Ok(Vec::new());
        };
        let mut ramps = Vec::with_capacity(source.branches.len());
        for branch in &source.branches {
            check(cancelled)?;
            let ramp = PronyRamp::new(dt_s / branch.source.relaxation_time_s());
            let half_dt = 0.5 * dt_s;
            let scale = finite((half_dt * half_dt) * (2.0 * ramp.apply(0.0, 1.0).1))?;
            material::scatter(
                effective,
                &self.model,
                branch.element,
                &branch.stiffness,
                scale,
                cancelled,
            )?;
            ramps.push(ramp);
        }
        Ok(ramps)
    }

    pub(super) fn relaxation_resistance(
        &self,
        initial: &MitcTransientState,
        ramps: &[PronyRamp],
        half_dt: f64,
        rhs: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PlateError> {
        let Some(history) = &initial.relaxation else {
            return Ok(());
        };
        for ((branch, q), ramp) in history
            .source
            .branches
            .iter()
            .zip(&history.coordinates)
            .zip(ramps)
        {
            check(cancelled)?;
            let mut relative = *q;
            let nodal = core::array::from_fn(|i| q[i]);
            relative[..18]
                .copy_from_slice(&crate::shell::mitc::element::relative_translations(&nodal));
            let mut local = [0.0; 20];
            for (force, row) in local.iter_mut().zip(&branch.stiffness) {
                *force = product(row, &relative)?;
            }
            let scale = -half_dt * ramp.apply(1.0, 0.0).1;
            let triangle = self.model.mesh.tris[branch.element];
            for i in 0..18 {
                let node = triangle[i / 6];
                let Some(row) = self.model.pencil.dof_map[6 * node + i % 6] else {
                    continue;
                };
                let block = 6 * (i / 6) + 3 * ((i % 6) / 3);
                let force: f64 = (0..3)
                    .map(|axis| self.model.frames[node][i % 3][axis] * local[block + axis])
                    .sum();
                rhs[row] = finite(rhs[row] + scale * force)?;
            }
            for axis in 0..2 {
                let row = self.model.pencil.free + 2 * branch.element + axis;
                rhs[row] = finite(rhs[row] + scale * local[18 + axis])?;
            }
        }
        Ok(())
    }

    pub(super) fn advance_relaxation(
        &self,
        initial: &MitcTransientState,
        velocity: &[f64],
        dt_s: f64,
        ramps: &[PronyRamp],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Option<History>, f64, f64, Vec<f64>), PlateError> {
        let Some(history) = &initial.relaxation else {
            return Ok((None, 0.0, 0.0, Vec::new()));
        };
        let nodal_count = self.model.pencil.free;
        let velocity_full = self.model.expand(&velocity[..nodal_count])?;
        let mut coordinates = Vec::with_capacity(history.coordinates.len());
        let mut facets = vec![0.0; self.reference_fields.len()];
        let (mut change, mut stored, mut heat) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for ((branch, q), ramp) in history
            .source
            .branches
            .iter()
            .zip(&history.coordinates)
            .zip(ramps)
        {
            check(cancelled)?;
            let triangle = self.model.mesh.tris[branch.element];
            let du: [f64; 20] = core::array::from_fn(|i| {
                dt_s * if i < 18 {
                    velocity_full[6 * triangle[i / 6] + i % 6]
                } else {
                    velocity[nodal_count + 2 * branch.element + i - 18]
                }
            });
            let next = core::array::from_fn(|i| ramp.apply(q[i], du[i]).0);
            if du.iter().chain(&next).any(|x| !x.is_finite()) {
                return Err(invalid(
                    "enriched shell relaxation history is unrepresentable",
                ));
            }
            let (delta, energy, loss) = self.reference_fields[branch.element].relaxation_interval(
                q,
                &du,
                &next,
                branch.source.factor(),
                ramp,
                cancelled,
            )?;
            change = change + Dd::from_f64(delta);
            stored = stored + Dd::from_f64(energy);
            heat = heat + Dd::from_f64(loss);
            facets[branch.element] = finite(facets[branch.element] + loss)?;
            coordinates.push(next);
        }
        Ok((
            Some(History {
                source: Arc::clone(&history.source),
                coordinates,
                stored_energy_j: finite(stored.to_f64())?,
            }),
            finite(change.to_f64())?,
            finite(heat.to_f64())?,
            facets,
        ))
    }
}
