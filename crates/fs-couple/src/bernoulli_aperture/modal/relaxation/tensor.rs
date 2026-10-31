//! Complete 3D Maxwell materials on the original retained DKT reed basis.
use super::super::exact::{CoupledArm, ExactModes};
use super::*;
use crate::tensor_material::{PlaneStressPoint, PlaneStressRelaxation, TensorMaterialError};
use core::ops::Range;
use fs_material::visco_tensor::{TensorProjectionBudget, TensorProny};
use fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM;

/// Explicit total material-quadrature/history/work ceilings. The original
/// exact structural block also imposes its unchanged complete-state limit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TensorModalRelaxationLimits {
    pub max_material_points: usize,
    pub max_memory_states: usize,
    pub max_projection_terms: usize,
}

/// One original area quadrature point, retaining its exact normal constraint.
/// Homogeneous linear bending integrates thickness analytically. All six
/// histories of every source branch survive; no weak loss floor is added.
#[derive(Clone, Debug)]
pub struct TensorModalMaterialPoint {
    blade: usize,
    region: usize,
    modes: Range<usize>,
    memories: Range<usize>,
    image: PlaneStressRelaxation,
}
impl TensorModalMaterialPoint {
    pub const fn blade(&self) -> usize {
        self.blade
    }
    pub const fn region(&self) -> usize {
        self.region
    }
    pub fn mode_range(&self) -> Range<usize> {
        self.modes.clone()
    }
    pub fn memory_range(&self) -> Range<usize> {
        self.memories.clone()
    }
    pub fn image(&self) -> &PlaneStressRelaxation {
        &self.image
    }
    /// Full tensor-shear strain at the positive physical outer face. The other
    /// face has its negative. This enforces the 3D strain limit at h/2, not just
    /// at the h/sqrt(12) energy-equivalent thickness point.
    pub fn surface_strain(
        &self,
        modes: &[BladeState],
        memory: &[f64],
    ) -> Result<fs_material::Voigt, AcousticRealizeError> {
        let q = modes
            .get(self.modes.clone())
            .ok_or_else(|| invalid("tensor material mode layout changed"))?;
        let z = memory
            .get(self.memories.clone())
            .ok_or_else(|| invalid("tensor material history layout changed"))?;
        let surface = det::sqrt(3.0);
        self.image
            .strain(
                &q.iter()
                    .map(|m| surface * m.displacement_m)
                    .collect::<Vec<_>>(),
                &z.iter().map(|z| surface * z).collect::<Vec<_>>(),
            )
            .map_err(super::super::nonlinear)
    }
}

fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

fn material_error(error: TensorMaterialError) -> AcousticRealizeError {
    match error {
        TensorMaterialError::Cancelled => AcousticRealizeError::Cancelled,
        error => super::super::nonlinear(error),
    }
}

impl ModalAperture {
    /// Attach full 3D Maxwell laws in original blade/region order.
    ///
    /// Inputs use the material's original 1/2/3 axes. The original section's
    /// ordered plane and angle rotate the COMPLETE law, including all branches.
    /// Its equilibrium plane-stress section must match the original elastic
    /// chart; that stiffness is already in the modes and is never added twice.
    /// Every DKT area point uses the shared full-energy normal elimination.
    /// The resulting exact passive channels enter the existing pressure/contact
    /// solve, swept volume, source work and independently integrated dissipation.
    ///
    /// This is fixed-state, homogeneous-through-thickness linear bending. It
    /// does not infer wet-cane data, mix independent Kelvin normal-rate terms,
    /// approximate a general tensor by a scalar spectrum, or claim measured
    /// instrument accuracy. Source bands, all histories and outer-face strain
    /// limits are admitted explicitly. No modes or material points are dropped.
    pub fn with_tensor_material_relaxation(
        mut self,
        laws: Vec<Vec<TensorProny>>,
        initial: InitialModalMemory,
        limits: TensorModalRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if self.samples != 0
            || self.relaxation.is_some()
            || self.viscosity.is_some()
            || self.blade_propagator() != BladePropagator::ExactZoh
            || self.material_plates.is_empty()
            || laws.len() != self.material_plates.len()
            || self.spec.modes.iter().any(|m| m.damping_n_s_m != 0.0)
            || !(1..=4096).contains(&limits.max_material_points)
            || limits.max_memory_states > 256
            || limits.max_projection_terms == 0
        {
            return Err(invalid(
                "tensor material memory needs cold exact original plates, complete laws and explicit budgets without duplicate damping",
            ));
        }
        let count = self.state.len();
        let maximum = (MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2)
            .checked_sub(3 * count)
            .ok_or_else(|| {
                invalid("complete tensor material memory exceeds exact propagation capacity")
            })?
            .min(limits.max_memory_states);
        let mut point_count = 0_usize;
        let mut memory_count = 0_usize;
        for (plate, laws) in self.material_plates.iter().zip(&laws) {
            if laws.len() != plate.source().regions().len()
                || plate.source().bending_viscosity().is_some()
            {
                return Err(invalid(
                    "tensor laws must cover every original material region without duplicate viscosity",
                ));
            }
            point_count = point_count
                .checked_add(3 * plate.source().chart().mesh.tris.len())
                .ok_or_else(|| invalid("tensor point extent overflow"))?;
            for (region, law) in plate.source().regions().iter().zip(laws) {
                region.input.model.require_canonical_relaxation()?;
                memory_count = law
                    .branches()
                    .len()
                    .checked_mul(18)
                    .and_then(|v| v.checked_mul(region.input.region.triangle_indices.len()))
                    .and_then(|v| memory_count.checked_add(v))
                    .ok_or_else(|| invalid("complete tensor material history extent overflow"))?;
            }
        }
        if point_count == 0 || point_count > limits.max_material_points || memory_count > maximum {
            return Err(invalid(
                "complete tensor points and histories exceed original capacity; no rows are removed",
            ));
        }
        // Disjoint thirds bound original DKT projection, point condensation and
        // exact history transformation. All limits are totals across all blades.
        let geometry_work = limits.max_projection_terms / 3 / self.material_plates.len();
        let point_work = limits.max_projection_terms / 3 / point_count;
        let mut points = Vec::with_capacity(point_count);
        let mut arms = Vec::<CoupledArm>::with_capacity(memory_count);
        for (blade, (plate, laws)) in self.material_plates.iter().zip(laws).enumerate() {
            poll(gate)?;
            let range = plate.mode_range();
            let n = range.len();
            let mut rotated = Vec::with_capacity(laws.len());
            for (region, law) in plate.source().regions().iter().zip(laws) {
                poll(gate)?;
                let law = if region.input.model == PlateMaterialModel::Orthotropic21 {
                    law.rotated([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
                        .map_err(super::super::nonlinear)?
                } else {
                    law
                };
                let law = if region.input.material_angle_rad == 0.0 {
                    law
                } else {
                    law.rotated(fs_material::tensor::rotation(
                        [0.0, 0.0, 1.0],
                        region.input.material_angle_rad,
                    ))
                    .map_err(super::super::nonlinear)?
                };
                // Compare at equilibrium, even if the source's operating band
                // begins above zero; this is constitutive admission, not a
                // claim to operate out of band. The source law is left intact.
                let section = law
                    .plane_stress(0.0)
                    .map_err(super::super::nonlinear)?
                    .section(region.section.thickness)
                    .map_err(super::super::nonlinear)?;
                for i in 0..3 {
                    for j in 0..3 {
                        let scale = det::sqrt(region.section.d[3 * i + i])
                            * det::sqrt(region.section.d[3 * j + j]);
                        let actual = section.bending[i][j];
                        if !scale.is_finite()
                            || scale <= 0.0
                            || actual.im != 0.0
                            || (actual.re - region.section.d[3 * i + j]).abs() > 1e-9 * scale
                        {
                            return Err(invalid(
                                "tensor equilibrium must reproduce the original material plate section",
                            ));
                        }
                    }
                }
                rotated.push(law);
            }
            let (basis, scale) = plate.constitutive_basis();
            let quadrature =
                basis.bending_strain_points(plate.source(), scale, geometry_work, gate)?;
            let first = arms.len();
            for point in quadrature {
                poll(gate)?;
                let image = PlaneStressPoint::new(
                    &rotated[point.region],
                    &point.columns,
                    point.volume_m3,
                    TensorProjectionBudget {
                        max_coordinates: n + 3,
                        max_scalar_terms: point_work,
                    },
                    gate,
                )
                .map_err(material_error)?
                .relaxation_image(point_work, gate)
                .map_err(material_error)?;
                let start = arms.len();
                for (row, rate) in image
                    .coupling_rows()
                    .iter()
                    .zip(image.relaxation_rates_per_s())
                {
                    if self.spec.time_step_s * rate > 2.0 {
                        return Err(invalid(
                            "complete constrained tensor relaxation is unresolved by the original clock",
                        ));
                    }
                    let mut coupling = vec![0.0; count];
                    coupling[range.clone()].copy_from_slice(row);
                    arms.push(CoupledArm {
                        coupling,
                        relaxation_time_s: rate.recip(),
                    });
                }
                points.push(TensorModalMaterialPoint {
                    blade,
                    region: point.region,
                    modes: range.clone(),
                    memories: start..arms.len(),
                    image,
                });
            }
            let mut tangent = vec![0.0; n * n];
            let low = self.spec.modes[range.clone()]
                .iter()
                .map(|m| det::sqrt(m.stiffness_n_m / m.mass_kg))
                .fold(f64::INFINITY, f64::min);
            for i in 0..n {
                let mode = &self.spec.modes[range.start + i];
                tangent[n * i + i] = mode.stiffness_n_m / mode.mass_kg;
            }
            for arm in &arms[first..] {
                for i in 0..n {
                    for j in 0..n {
                        tangent[n * i + j] += arm.coupling[range.start + i]
                            / det::sqrt(self.spec.modes[range.start + i].mass_kg)
                            * arm.coupling[range.start + j]
                            / det::sqrt(self.spec.modes[range.start + j].mass_kg);
                    }
                }
            }
            let high = det::sqrt(
                tangent
                    .chunks_exact(n)
                    .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
                    .fold(0.0, f64::max),
            );
            if !high.is_finite()
                || high * self.spec.time_step_s > plate.options().maximum_angular_step
            {
                return Err(invalid(
                    "complete instantaneous tensor stiffness exceeds the original structural clock cap",
                ));
            }
            for (region, law) in plate.source().regions().iter().zip(&rotated) {
                let band = law.band_rad_s();
                if low < band.0 || high > band.1 {
                    return Err(invalid(
                        "original tensor material band does not cover all retained equilibrium and instantaneous frequencies",
                    ));
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
                    property(&region.input.material, key, dims, [low, high])?;
                }
            }
        }
        let states = vec![
            arms.iter()
                .map(|arm| match initial {
                    InitialModalMemory::Relaxed => arm
                        .coupling
                        .iter()
                        .zip(&self.state)
                        .map(|(c, s)| c * s.displacement_m)
                        .sum(),
                    InitialModalMemory::Unrelaxed => 0.0,
                })
                .collect::<Vec<f64>>(),
        ];
        let relaxation = ModalRelaxation {
            spectra: Vec::new(),
            initial,
            arms: Vec::new(),
            states,
            coupled_arms: Some(arms.clone()),
            tensor_points: points,
        };
        relaxation.validate_material_state(&self.state, &relaxation.states)?;
        let energy = self.energy(&self.state) + relaxation.energy(&self.state, &relaxation.states);
        if !energy.is_finite() || energy > self.spec.limits.maximum_energy_j {
            return Err(invalid(
                "initial tensor material history exceeds original energy budget",
            ));
        }
        let exact = ExactModes::with_coupled_arms(&self.spec, arms)?;
        poll(gate)?;
        self.exact = Some(exact);
        self.relaxation = Some(relaxation);
        Ok(self)
    }
}
