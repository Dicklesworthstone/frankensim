//! Causal bulk absorption and dispersion in the complete pressure FE space.
//!
//! Original frozen compressibility M_inf and velocity potential obey
//! M_inf*p_dot + sum(delta_beta*M_region*q_dot) + K*psi + C*p = b,
//! psi_dot=p, tau*q_dot+q=region(p). Each region retains the full original
//! consistent volume mass and every original memory coefficient. Pressure
//! memory stores delta_beta*q^T*M_region*q/2; heat is the independent positive
//! delta_beta*(p-q)^T*M_region*(p-q)/tau. No frequency sample determines a time law.
//!
//! Midpoint condenses memories into ONE native sparse positive pressure factor.
//! Region mass forms are shared between their arms; no dense mass inverse,
//! pressure-mode projection or material query occurs during a step. Molecular
//! fields are allowed to jump across material interfaces. Source material
//! adapters must supply M_inf, rather than equilibrium compressibility, once.

use super::super::super::AcousticFeError;
use super::super::{BoundaryVelocitySource, VolumeVelocitySource};
use super::{
    HeterogeneousError, HeterogeneousHelmholtzFe, HeterogeneousTimeState, RealVelocityTrace,
    ResistiveBoundary, SpatialResistiveBoundary,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use crate::highorder::simplex::duffy_quadrature;
use fs_math::dd::Dd;
use fs_sparse::Coo;
use std::collections::BTreeSet;

mod physics;
pub use physics::{
    BulkMemoryArm, BulkMemoryWork, BulkPressureState, BulkTimeEnergy, BulkTimeError,
    BulkTimeLimits, BulkTimeStep, BulkTimeWork,
};
use physics::{PreparedBulkDynamics, Region, bound, poll};

/// Constant original material coefficients on a conforming set of tetrahedra.
/// Different regions are disjoint; reservoir order remains authored order.
#[derive(Clone, Debug)]
pub struct BulkMemoryRegion {
    pub tetrahedra: Vec<usize>,
    pub arms: Vec<BulkMemoryArm>,
}
#[derive(Debug)]
pub enum BulkInteriorError {
    Interior(HeterogeneousError),
    Dynamics(BulkTimeError),
}
impl core::fmt::Display for BulkInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Interior(e) => e.fmt(f),
            Self::Dynamics(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for BulkInteriorError {}
impl From<HeterogeneousError> for BulkInteriorError {
    fn from(e: HeterogeneousError) -> Self {
        Self::Interior(e)
    }
}
impl From<AcousticFeError> for BulkInteriorError {
    fn from(e: AcousticFeError) -> Self {
        Self::Interior(e.into())
    }
}
impl From<BulkTimeError> for BulkInteriorError {
    fn from(e: BulkTimeError) -> Self {
        Self::Dynamics(e)
    }
}

/// Source-bound original spatial operator plus complete volume memories.
pub struct PreparedBulkTime<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    dynamics: PreparedBulkDynamics,
    release_faces: BTreeSet<[u32; 3]>,
}
impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare frozen-bulk pressure fields with positive original material arms.
    ///
    /// This model's mass MUST contain beta_inf. The original molecular adapter
    /// constructs that coherent mass from its source. Assigning equilibrium
    /// compressibility here would double-count the relaxed storage. Density
    /// remains inside the unchanged stiffness divergence. Regions must conform
    /// to original element faces; smooth varying rates require another owner.
    ///
    /// Boundary resistance and exact pressure release retain their existing
    /// source semantics. A release fixes psi=p=0; volume memories remain live,
    /// including coefficients on that trace. Wall memory, mean flow and thermal
    /// diffusion are not silently inferred or added by this owner.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_bulk_time<'m>(
        &'m self,
        duration_s: f64,
        regions: &[BulkMemoryRegion],
        boundaries: &[ResistiveBoundary],
        spatial: &[&SpatialResistiveBoundary<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: BulkTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedBulkTime<'m, 'f, 'c>, BulkInteriorError> {
        poll(cancelled)?;
        bound("pressure coordinates", self.ndof(), limits.maximum_dofs)?;
        bound("memory regions", regions.len(), limits.maximum_regions)?;
        let constraints = self.geometry().pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let conductance = self.resistive_matrix(boundaries, spatial, &constraints, cancelled)?;
        let geometry = self.geometry();
        let n = self.ndof();
        let quad = duffy_quadrature(self.quadrature_order());
        let mut seen = BTreeSet::new();
        let mut prepared = Vec::with_capacity(regions.len());
        let (mut total_arms, mut memory_dofs, mut total_entries) = (0usize, 0usize, 0usize);
        for region in regions {
            poll(cancelled)?;
            if region.tetrahedra.is_empty() || region.arms.is_empty() {
                return Err(BulkTimeError::Invalid(
                    "nonempty original material region and memory arms required",
                )
                .into());
            }
            total_arms = total_arms
                .checked_add(region.arms.len())
                .ok_or(BulkTimeError::Invalid("material arm extent overflow"))?;
            bound("memory arms", total_arms, limits.maximum_arms)?;
            let mut indices = BTreeSet::new();
            for &element in &region.tetrahedra {
                if element >= geometry.tetrahedra().len() || !seen.insert(element) {
                    return Err(BulkTimeError::Invalid(
                        "material memory regions must be disjoint original tetrahedra",
                    )
                    .into());
                }
                indices.extend(geometry.space.element_dofs(element).iter().map(|(i, _)| *i));
            }
            let dofs: Vec<_> = indices.into_iter().collect();
            memory_dofs = memory_dofs
                .checked_add(
                    dofs.len()
                        .checked_mul(region.arms.len())
                        .ok_or(BulkTimeError::Invalid("material memory extent overflow"))?,
                )
                .ok_or(BulkTimeError::Invalid("material memory extent overflow"))?;
            bound(
                "memory coordinates",
                memory_dofs,
                limits.maximum_memory_dofs,
            )?;
            let mut map = vec![usize::MAX; n];
            for (i, &dof) in dofs.iter().enumerate() {
                map[dof] = i;
            }
            let mut mass = Coo::new(dofs.len(), dofs.len());
            for &element in &region.tetrahedra {
                poll(cancelled)?;
                let basis = geometry.space.element_dofs(element);
                let count = basis.len();
                let mut local = vec![0.0; count * count];
                let mut values = vec![0.0; count];
                for &(barycentric, weight) in &quad {
                    poll(cancelled)?;
                    let jacobian = if geometry.is_curved() {
                        geometry
                            .geometry_point(element, barycentric)?
                            .determinant_m3
                    } else {
                        6.0 * geometry.affine_geometry.vol_signed[element].abs()
                    };
                    let w = weight * jacobian;
                    if !w.is_finite() || w <= 0.0 {
                        return Err(BulkTimeError::Invalid(
                            "positive finite original volume quadrature required",
                        )
                        .into());
                    }
                    for (i, (_, shape)) in basis.iter().enumerate() {
                        values[i] = shape.eval(barycentric, geometry.order());
                    }
                    for i in 0..count {
                        for j in i..count {
                            local[i * count + j] += w * values[i] * values[j];
                        }
                    }
                }
                for (i, &(row, _)) in basis.iter().enumerate() {
                    for (j, &(col, _)) in basis.iter().enumerate().skip(i) {
                        let v = local[i * count + j];
                        if !v.is_finite() {
                            return Err(
                                BulkTimeError::Invalid("original volume mass overflow").into()
                            );
                        }
                        bound(
                            "region staging entries",
                            mass.staged() + if row == col { 1 } else { 2 },
                            limits.maximum_matrix_entries,
                        )?;
                        mass.push(map[row], map[col], v);
                        if row != col {
                            mass.push(map[col], map[row], v);
                        }
                    }
                }
            }
            let mass = mass.assemble();
            total_entries = total_entries
                .checked_add(mass.nnz())
                .ok_or(BulkTimeError::Invalid("material matrix extent overflow"))?;
            bound(
                "region matrix entries",
                total_entries,
                limits.maximum_matrix_entries,
            )?;
            prepared.push(Region {
                dofs,
                mass,
                arms: region.arms.clone(),
            });
        }
        let dynamics = PreparedBulkDynamics::new(
            self.mass().clone(),
            self.stiffness().clone(),
            conductance,
            prepared,
            constraints.free_dofs,
            duration_s,
            limits,
            cancelled,
        )?;
        Ok(PreparedBulkTime {
            model: self,
            dynamics,
            release_faces: constraints.faces,
        })
    }
}
impl<'m, 'f, 'c> PreparedBulkTime<'m, 'f, 'c> {
    pub const fn model(&self) -> &'m HeterogeneousHelmholtzFe<'f, 'c> {
        self.model
    }
    pub fn duration_s(&self) -> f64 {
        self.dynamics.duration()
    }
    pub fn arm_count(&self) -> usize {
        self.dynamics.arm_count()
    }
    pub fn arm_dofs(&self, arm: usize) -> Option<&[usize]> {
        self.dynamics.arm_dofs(arm)
    }
    pub fn zero_state(&self) -> BulkPressureState<'_> {
        self.dynamics.zero_state()
    }
    /// Import every original history explicitly. Neither relaxed nor frozen
    /// histories are inferred from the pressure endpoints.
    pub fn state(
        &self,
        potential: &[f64],
        pressure: &[f64],
        memories: &[Vec<f64>],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BulkPressureState<'_>, BulkTimeError> {
        self.dynamics
            .state(potential, pressure, memories, cancelled)
    }
    pub fn energy(
        &self,
        state: &BulkPressureState<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BulkTimeEnergy, BulkTimeError> {
        self.dynamics.energy(state, cancelled)
    }
    /// Copy physical pressure/potential once for the existing spatial probes.
    /// The lossless state's energy omits molecular storage; use energy above
    /// for the complete physical account. This export does not rebind memories.
    pub fn acoustic_state(
        &self,
        state: &BulkPressureState<'_>,
    ) -> Result<HeterogeneousTimeState<'m, 'f, 'c>, BulkTimeError> {
        self.dynamics.check(state)?;
        Ok(HeterogeneousTimeState {
            model: self.model,
            potential: state.pressure_potential_pa_s().to_vec(),
            pressure: state.pressure_pa().to_vec(),
        })
    }
    /// Already integrated original load [m^3/s] in the same pressure basis.
    /// Its conjugate work is h*p_mid^T*load; no density or frequency multiplier.
    pub fn step_projected<'p>(
        &'p self,
        state: &BulkPressureState<'p>,
        load: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BulkTimeStep<'p>, BulkTimeError> {
        self.dynamics.step(state, load, cancelled)
    }
    /// Simultaneous linear boundary velocity and original prepared source maps.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_sources<'p>(
        &'p self,
        initial: &BulkPressureState<'p>,
        inward: &[RealVelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, f64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BulkTimeStep<'p>, BulkInteriorError> {
        poll(cancelled)?;
        self.dynamics.check(initial)?;
        super::super::boundary_source::admit_count(boundary_sources.len())?;
        super::super::source::admit_count(volume_sources.len())?;
        for &(face, _) in inward {
            let mut key = face;
            key.sort_unstable();
            if self.release_faces.contains(&key) {
                return Err(BulkTimeError::Invalid(
                    "release and imposed-velocity faces must not overlap",
                )
                .into());
            }
        }
        for &(source, amplitude) in boundary_sources {
            source.validate(self.model)?;
            if !amplitude.is_finite() || source.faces.iter().any(|f| self.release_faces.contains(f))
            {
                return Err(BulkTimeError::Invalid(
                    "finite source amplitude and disjoint pressure-release support required",
                )
                .into());
            }
        }
        for &(source, amplitude) in volume_sources {
            source.validate(self.model)?;
            if !amplitude.is_finite() {
                return Err(
                    BulkTimeError::Invalid("finite volume-source amplitude required").into(),
                );
            }
        }
        let mut load = self.model.geometry().boundary_linear_load(inward)?;
        for (i, value) in load.iter_mut().enumerate() {
            poll(cancelled)?;
            let mut s = Dd::from_f64(*value);
            for &(source, amplitude) in boundary_sources {
                s = s + Dd::from_f64(source.load[i]) * Dd::from_f64(amplitude);
            }
            for &(source, amplitude) in volume_sources {
                s = s + Dd::from_f64(source.load[i]) * Dd::from_f64(amplitude);
            }
            *value = s.to_f64();
        }
        let mut step = self.step_projected(initial, &load, cancelled)?;
        let mut boundary = Dd::ZERO;
        let mut volume = Dd::ZERO;
        for &(s, a) in boundary_sources {
            boundary = boundary
                + Dd::from_f64(self.duration_s())
                    * Dd::from_f64(a)
                    * Dd::from_f64(
                        s.midpoint_pressure(initial.pressure_pa(), step.state.pressure_pa()),
                    );
        }
        for &(s, a) in volume_sources {
            volume = volume
                + Dd::from_f64(self.duration_s())
                    * Dd::from_f64(a)
                    * Dd::from_f64(
                        s.midpoint_pressure(initial.pressure_pa(), step.state.pressure_pa()),
                    );
        }
        step.work.boundary_source_work_j = boundary.to_f64();
        step.work.volume_source_work_j = volume.to_f64();
        if !step.work.boundary_source_work_j.is_finite()
            || !step.work.volume_source_work_j.is_finite()
        {
            return Err(BulkTimeError::Invalid("original source-work overflow").into());
        }
        poll(cancelled)?;
        Ok(step)
    }
}

#[cfg(test)]
mod tests;
