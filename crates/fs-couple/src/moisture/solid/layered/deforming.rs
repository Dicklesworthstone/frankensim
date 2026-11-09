//! Spatial Fick diffusion pulled back through the actual compatible swelling.
//!
//! For dry-basis M, j = -rho_d,current D grad_x(M). Conservation of the dry
//! skeleton gives the reference flux Q = -rho_d,reference F^-1 D F^-T grad_X(M).
//! An axis-aligned, transversely uniform prism has axial metric 1/lambda_x².
//! Its changing area therefore cancels from bulk transport, but still multiplies
//! an explicitly area-specific surface/contact transfer law. This distinction
//! prevents replacing conserved dry mass by a new density/volume calculation.
//!
//! The selected scalar D is explicitly interpreted in the CURRENT spatial frame
//! by this opt-in API. The source must support that constitutive interpretation.
//! Reference-frame diffusivity belongs to the ordinary `advance` method. Native
//! mechanics remains small-strain quasistatic; this geometric correction does
//! not introduce stress-driven chemical potentials, strain-dependent D, sorption
//! heat, general ALE transport, a finite-strain solid, or free-liquid transport.

use super::*;
use crate::moisture::layered::MoistureLayeredMetric;

#[path = "deforming/atmosphere.rs"]
mod atmosphere;
#[path = "deforming/solution.rs"]
mod solution;
pub use atmosphere::{
    MoistureSpatialAtmosphereCandidate, MoistureSpatialAtmosphereControls,
    MoistureSpatialAtmosphereReport,
};
pub use solution::{
    MoistureSpatialPairedSolutionCandidate, MoistureSpatialPairedSolutionControls,
    MoistureSpatialPairedSolutionReport,
};

// One counter follows every RH trial and nonlinear geometry evaluation.
// A native mechanical solve is reserved BEFORE entering the original owner.
struct SpatialSolveBudget {
    used: usize,
    maximum: usize,
}
impl SpatialSolveBudget {
    fn reserve(&mut self, cx: &Cx<'_>) -> Result<(), MoistureError> {
        poll(cx)?;
        if self.used >= self.maximum {
            return Err(MoistureError::Invalid(
                "cumulative spatial solid-solve budget exhausted",
            ));
        }
        self.used += 1;
        Ok(())
    }
}

/// Original conservative transport, together with the measured spatial metric.
#[derive(Clone, Debug)]
pub struct MoistureSpatialDiffusionReport {
    /// End-state fluxes and residuals use the actual end-state geometry.
    pub transport: MoistureLayeredReport,
    /// Actual axial length/reference length for each original moisture cell.
    pub axial_stretch: Vec<f64>,
    /// Common current yz area [m²], also used by finite transfer laws.
    pub current_cross_section_m2: f64,
    /// Largest nodal departure from the admitted separable prism map [m].
    pub maximum_geometry_residual_m: f64,
    /// Original native geometry/source identity of the returned solid state.
    pub geometry_identity: ContentHash,
    /// Initial plus proposed/relaxed endpoint native mechanical evaluations.
    pub native_solid_solves: usize,
}

struct AxialMetric {
    pullback: MoistureLayeredMetric,
    axial_stretch: Vec<f64>,
    area_m2: f64,
    residual_m: f64,
}

impl AxialMetric {
    fn from_geometry(
        reference: &[[f64; 3]],
        current: &[[f64; 3]],
        tolerance_m: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        if reference.len() != current.len() || current.len() < 8 || current.len() % 4 != 0 {
            return Err(MoistureError::Invalid(
                "spatial moisture requires the original prism mesh",
            ));
        }
        let reference_area =
            (reference[1][1] - reference[0][1]) * (reference[3][2] - reference[0][2]);
        let width = current[1][1] - current[0][1];
        let thickness = current[3][2] - current[0][2];
        let area_m2 = width * thickness;
        let area_ratio = area_m2 / reference_area;
        if !positive(width) || !positive(thickness) || !positive(area_ratio) {
            return Err(MoistureError::Invalid(
                "spatial moisture section is unrepresentable",
            ));
        }
        let section = [
            [current[0][1], current[0][2]],
            [current[1][1], current[0][2]],
            [current[1][1], current[3][2]],
            [current[0][1], current[3][2]],
        ];
        let mut residual_m: f64 = 0.0;
        // A spatially isotropic law and a field depending only on X admit this
        // exact 1-D reduction only when F is diagonal. Do not silently project
        // bending, shear, taper or transverse warping into an axial coefficient.
        for face in current.chunks_exact(4) {
            poll(cx)?;
            for (node, yz) in face.iter().zip(section) {
                residual_m = residual_m.max((node[0] - face[0][0]).abs());
                residual_m = residual_m.max((node[1] - yz[0]).abs());
                residual_m = residual_m.max((node[2] - yz[1]).abs());
            }
        }
        if !residual_m.is_finite() || residual_m > tolerance_m {
            return Err(MoistureError::Invalid(
                "spatial moisture requires axis-aligned sections without shear, bending or taper",
            ));
        }
        let cells = current.len() / 4 - 1;
        let mut axial_stretch = Vec::with_capacity(cells);
        let mut inverse_axial_stretch_squared = Vec::with_capacity(cells);
        for i in 0..cells {
            poll(cx)?;
            let reference_dx = reference[4 * (i + 1)][0] - reference[4 * i][0];
            let current_dx = current[4 * (i + 1)][0] - current[4 * i][0];
            let stretch = current_dx / reference_dx;
            let inverse = reference_dx / current_dx;
            let metric = inverse * inverse;
            if !positive(reference_dx) || !positive(stretch) || !positive(metric) {
                return Err(MoistureError::Invalid(
                    "spatial moisture cell stretch is unrepresentable",
                ));
            }
            axial_stretch.push(stretch);
            inverse_axial_stretch_squared.push(metric);
        }
        Ok(Self {
            pullback: MoistureLayeredMetric {
                inverse_axial_stretch_squared,
                face_area_ratio: area_ratio,
            },
            axial_stretch,
            area_m2,
            residual_m,
        })
    }
}

impl MoistureLayeredSolid {
    /// Solve moisture and its compatible swelling together, with the selected
    /// scalar diffusivity explicitly interpreted as current-frame isotropic D.
    ///
    /// Every nonlinear candidate invokes the original source-bound native
    /// solid solve, then evaluates the conservative water residual using that
    /// candidate's actual axial lengths and exposed/contact area. The original
    /// nonlinear iteration limit also bounds these geometry evaluations (at
    /// most twice that limit plus the initial state when relaxation is used).
    /// The same source support, maximum principle, true linear/nonlinear
    /// residuals and actual water balance remain hard acceptance conditions.
    ///
    /// The geometry must remain a separable axis-aligned prism: x may stretch
    /// differently per cell, while the yz section stays common and rectangular.
    /// `maximum_geometry_residual_m` bounds numerical departure from that map.
    /// A shear/bending/tapered candidate outside this tolerance is refused;
    /// it needs a spatial transport formulation. No initial state is mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_spatial_diffusion(
        &self,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        faces: &[MoistureLayeredFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        maximum_geometry_residual_m: f64,
        cx: &Cx<'_>,
    ) -> Result<MoistureLayeredSolidStep<MoistureSpatialDiffusionReport>, MoistureError> {
        let mut work = SpatialSolveBudget {
            used: 0,
            maximum: usize::MAX,
        };
        self.advance_spatial_diffusion_budgeted(
            moisture,
            dt_s,
            faces,
            transport,
            nonlinear,
            maximum_geometry_residual_m,
            &mut work,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_spatial_diffusion_budgeted(
        &self,
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        faces: &[MoistureLayeredFace; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        maximum_geometry_residual_m: f64,
        work: &mut SpatialSolveBudget,
        cx: &Cx<'_>,
    ) -> Result<MoistureLayeredSolidStep<MoistureSpatialDiffusionReport>, MoistureError> {
        poll(cx)?;
        if !positive(maximum_geometry_residual_m) {
            return Err(MoistureError::Invalid(
                "positive spatial geometry tolerance required",
            ));
        }
        let mut final_state = None;
        let initial_solid_solves = work.used;
        let mut resolve_metric = |candidate: &MoistureLayeredSlab| {
            work.reserve(cx)?;
            let solid = self.resolve(candidate, cx)?;
            let geometry = AxialMetric::from_geometry(
                &self.nodes_m,
                solid.geometry().nodes_m(),
                maximum_geometry_residual_m,
                cx,
            )?;
            let pullback = geometry.pullback;
            final_state = Some((
                solid,
                geometry.axial_stretch,
                geometry.area_m2,
                geometry.residual_m,
            ));
            Ok(pullback)
        };
        let (moisture, transport) = moisture.advance_with_metric(
            dt_s,
            faces,
            transport,
            nonlinear,
            Some(&mut resolve_metric),
            cx,
        )?;
        let (solid, axial_stretch, current_cross_section_m2, maximum_geometry_residual_m) =
            final_state.expect("successful metric endpoint retains its native solid");
        let report = MoistureSpatialDiffusionReport {
            transport,
            axial_stretch,
            current_cross_section_m2,
            maximum_geometry_residual_m,
            geometry_identity: solid.geometry().identity(),
            native_solid_solves: work.used - initial_solid_solves,
        };
        poll(cx)?;
        Ok(MoistureLayeredSolidStep {
            moisture,
            transport: report,
            solid,
        })
    }
}
