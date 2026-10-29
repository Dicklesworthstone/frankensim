//! Complete source-resolved tensor memory on the original DKT body basis.
//! The shared material/plane-stress owners retain every normal and shear
//! history; only their exact change of coordinates is used by the modal owner.
use super::*;
use crate::tensor_material::{PlaneStressPoint, PlaneStressRelaxation, TensorMaterialError};
use core::ops::Range;
use fs_material::visco_tensor::source::{TensorElasticity, TensorPronyProperties};
use fs_material::visco_tensor::{TensorProjectionBudget, TensorProny};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TensorPlateRelaxationLimits {
    pub max_modes: usize,
    pub max_material_points: usize,
    pub max_memory_states: usize,
    pub max_projection_terms: usize,
}

/// One original area point; thickness uses exact homogeneous bending energy.
#[derive(Clone, Debug)]
pub struct TensorPlateMaterialPoint {
    region: usize,
    memory: Range<usize>,
    image: PlaneStressRelaxation,
}
impl TensorPlateMaterialPoint {
    pub const fn region(&self) -> usize {
        self.region
    }
    pub fn memory_range(&self) -> Range<usize> {
        self.memory.clone()
    }
    pub fn image(&self) -> &PlaneStressRelaxation {
        &self.image
    }
}

/// Full rotated source laws, original geometry and all passive history rows.
#[derive(Clone, Debug)]
pub struct TensorPlateRelaxation {
    chart_identity: ContentHash,
    properties: Vec<TensorPronyProperties>,
    laws: Vec<TensorProny>,
    bands: Vec<(f64, f64)>,
    points: Vec<TensorPlateMaterialPoint>,
    rows: Vec<Vec<f64>>,
    times: Vec<f64>,
    modes: usize,
    limits: TensorPlateRelaxationLimits,
}
fn tensor_error(error: TensorMaterialError) -> AcousticRealizeError {
    match error {
        TensorMaterialError::Cancelled => AcousticRealizeError::Cancelled,
        error => owner(error),
    }
}
impl TensorPlateRelaxation {
    pub fn chart_identity(&self) -> ContentHash {
        self.chart_identity
    }
    pub fn properties(&self) -> &[TensorPronyProperties] {
        &self.properties
    }
    pub fn material_laws(&self) -> &[TensorProny] {
        &self.laws
    }
    pub fn material_points(&self) -> &[TensorPlateMaterialPoint] {
        &self.points
    }
    pub fn coupling_rows_per_s(&self) -> &[Vec<f64>] {
        &self.rows
    }
    pub fn relaxation_times_s(&self) -> &[f64] {
        &self.times
    }
    pub const fn mode_count(&self) -> usize {
        self.modes
    }
    pub const fn limits(&self) -> TensorPlateRelaxationLimits {
        self.limits
    }

    /// Bind original source properties, never authored matrix substitutes.
    /// Equilibrium elasticity remains counted once in the existing plate modes.
    pub(crate) fn from_chart(
        chart: &ResolvedPlateChart,
        basis: &RetainedPlateBasis,
        properties: &[TensorPronyProperties],
        limits: TensorPlateRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let n = basis.mass_normalized_modes.len();
        if n == 0
            || n > limits.max_modes
            || !(1..=64).contains(&limits.max_modes)
            || !(1..=4096).contains(&limits.max_material_points)
            || limits.max_memory_states > 256
            || limits.max_projection_terms == 0
            || chart.bending_viscosity().is_some()
            || properties.len() != chart.regions().len()
            || properties.is_empty()
        {
            return Err(invalid(
                "complete tensor plate memory needs original regions, modes and explicit point/history/work budgets without parallel viscosity",
            ));
        }
        let points = chart
            .chart()
            .mesh
            .tris
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("tensor plate quadrature count overflow"))?;
        let mut states = 0_usize;
        for (region, selection) in chart.regions().iter().zip(properties) {
            let expected = if region.input.model == PlateMaterialModel::Isotropic {
                TensorElasticity::Isotropic
            } else {
                TensorElasticity::Orthotropic
            };
            if selection.elasticity != expected {
                return Err(invalid(
                    "tensor equilibrium family must match the original material chart",
                ));
            }
            states = selection
                .branches
                .len()
                .checked_mul(18)
                .and_then(|v| v.checked_mul(region.input.region.triangle_indices.len()))
                .and_then(|v| states.checked_add(v))
                .ok_or_else(|| invalid("complete tensor plate history count overflow"))?;
        }
        if points == 0 || points > limits.max_material_points || states > limits.max_memory_states {
            return Err(invalid(
                "complete tensor plate quadrature or histories exceed their explicit budget",
            ));
        }
        let mut laws = Vec::with_capacity(properties.len());
        let mut bands = Vec::with_capacity(properties.len());
        for (region, selection) in chart.regions().iter().zip(properties) {
            poll(gate)?;
            let mut law = selection.bind(&region.input.material).map_err(owner)?;
            if region.input.model == PlateMaterialModel::Orthotropic21 {
                law = law
                    .rotated([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]])
                    .map_err(owner)?;
            }
            if region.input.material_angle_rad != 0.0 {
                law = law
                    .rotated(fs_material::tensor::rotation(
                        [0.0, 0.0, 1.0],
                        region.input.material_angle_rad,
                    ))
                    .map_err(owner)?;
            }
            let section = law
                .plane_stress(0.0)
                .map_err(owner)?
                .section(region.section.thickness)
                .map_err(owner)?;
            for i in 0..3 {
                for j in 0..3 {
                    let scale = det::sqrt(region.section.d[3 * i + i])
                        * det::sqrt(region.section.d[3 * j + j]);
                    if !scale.is_finite()
                        || scale <= 0.0
                        || section.bending[i][j].im != 0.0
                        || (section.bending[i][j].re - region.section.d[3 * i + j]).abs()
                            > 1e-9 * scale
                    {
                        return Err(invalid(
                            "complete tensor equilibrium must reproduce the original plate section",
                        ));
                    }
                }
            }
            let mut band = law.band_rad_s();
            coefficient(&region.input.material, "density", Density::DIMS, &mut band)?;
            if band.1 < band.0 {
                return Err(invalid(
                    "tensor solid and original mass source bands do not intersect",
                ));
            }
            laws.push(law);
            bands.push(band);
        }
        let point_work = limits.max_projection_terms / 3 / points;
        let quadrature = basis.bending_strain_points(
            chart,
            &vec![1.0; n],
            limits.max_projection_terms / 3,
            gate,
        )?;
        let mut images = Vec::with_capacity(points);
        let mut rows = Vec::with_capacity(states);
        let mut times = Vec::with_capacity(states);
        for point in quadrature {
            poll(gate)?;
            let image = PlaneStressPoint::new(
                &laws[point.region],
                &point.columns,
                point.volume_m3,
                TensorProjectionBudget {
                    max_coordinates: n + 3,
                    max_scalar_terms: point_work,
                },
                gate,
            )
            .map_err(tensor_error)?
            .relaxation_image(point_work, gate)
            .map_err(tensor_error)?;
            let start = rows.len();
            rows.extend_from_slice(image.coupling_rows());
            times.extend(
                image
                    .relaxation_rates_per_s()
                    .iter()
                    .map(|rate| rate.recip()),
            );
            images.push(TensorPlateMaterialPoint {
                region: point.region,
                memory: start..rows.len(),
                image,
            });
        }
        if rows.len() != states {
            return Err(invalid(
                "tensor projection must retain every original material history",
            ));
        }
        poll(gate)?;
        Ok(Self {
            chart_identity: chart.specimen_identity(),
            properties: properties.to_vec(),
            laws,
            bands,
            points: images,
            rows,
            times,
            modes: n,
            limits,
        })
    }

    pub(crate) fn validate_band(&self, high: f64, dt: f64) -> Result<(), AcousticRealizeError> {
        if !high.is_finite()
            || high < 0.0
            || !dt.is_finite()
            || dt <= 0.0
            || self
                .bands
                .iter()
                .any(|&(low, upper)| low > 0.0 || high > upper)
            || self.times.iter().any(|time| dt / time > 2.0)
        {
            return Err(invalid(
                "complete coupled spectrum or clock exceeds original tensor material support",
            ));
        }
        Ok(())
    }

    pub(crate) fn bindings(
        &self,
        mode_start: usize,
        memory_start: usize,
    ) -> Vec<crate::render::schedule::force::coupled::ModalPlaneStressPoint> {
        self.points
            .iter()
            .map(
                |point| crate::render::schedule::force::coupled::ModalPlaneStressPoint {
                    modes: mode_start..mode_start + self.modes,
                    memories: memory_start + point.memory.start..memory_start + point.memory.end,
                    strain_scale: det::sqrt(3.0),
                    image: point.image.clone(),
                },
            )
            .collect()
    }
}
