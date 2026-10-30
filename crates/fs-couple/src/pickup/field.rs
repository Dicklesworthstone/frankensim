//! Source material and moving geometry -> native magnetic fields -> Faraday EMF.
//!
//! The cold compiler solves the existing planar fixed-remanence winding problem
//! at declared gaps. It keeps the original cell materials, coil turns, magnet
//! orientation and mesh lineage. Permeability and remanence changes therefore
//! change the field and voltage without any metal-name or density heuristic.
//! Central differences and independent midpoint solves check the interpolation
//! at explicit SI tolerances; neither is a continuum error certificate.
//!
//! This is an open-circuit, magnetostatic, one-displacement observation. It does
//! not add static magnetic attraction, current-dependent inductance, electrical
//! loading, eddy currents, hysteresis or a string time integrator. The supplied
//! planar geometry and the pickup's aperture reduction still need physical
//! validation. A changed material state requires a new field map.

use super::{Pickup, PickupError, PickupFluxKnot, PickupFluxMap, PickupPose};
use crate::acoustic_realize::LinearMaterialStringRuntime;
use fs_blake3::ContentHash;
use fs_em::{
    planar::{MagneticError, PlanarMagnetostatics, permanent::PlanarMagneticMaterial},
    winding::{PlanarWindingBasis, WindingConfig, WindingError},
};
use fs_exec::Cx;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_math::det;
use fs_rep_mesh::TriComplex2;
use fs_solver::{CheckedCgConfig, SolveReport};
use std::sync::Arc;

/// One frozen physical field specimen, before any gap sampling.
#[derive(Clone, Debug)]
pub struct PlanarPickupSpec {
    /// Original planar mesh, thickness, material-cell order and feature keys.
    pub mesh: TriComplex2,
    /// Positive physical gap represented by the original vertices [m].
    pub reference_gap_m: f64,
    /// Vertex displacement per metre of increasing gap, in the mesh frame.
    /// All vertices of a magnet or winding cell must be stationary. Adjacent
    /// exterior cells may deform with the mesh; inverted triangles refuse.
    pub vertex_gap_directions: Vec<[f64; 2]>,
    /// Original constitutive model for every triangle, including the exterior.
    pub materials: Vec<PlanarMagneticMaterial>,
    /// Signed winding density J_z/I [1/m²], including turns and return path.
    /// Nonzero cells stay fixed, so their original turn counts do not drift.
    pub winding_profile_per_m2: Vec<f64>,
    /// Complete exterior edges held at A_z=0, as admitted by the field owner.
    pub prescribed_vertices: Vec<usize>,
    /// Explicit series end/leakage inductance [H]; zero neglects the ends.
    /// It is used only by native winding admission, not as an output filter.
    pub end_inductance_h: f64,
    /// Original induced-material cells representing the sensed moving body.
    /// A material-string binding checks their card, state and rigid translation.
    /// The field compiler also supports deformation, without that binding.
    pub sensed_cells: Vec<usize>,
}

/// Independent field, finite-difference and interpolation work/accuracy controls.
#[derive(Clone, Copy, Debug)]
pub struct PickupFieldConfig {
    /// Maximum complete field evaluations. Each evaluation uses one bias and
    /// one unit-current solve on the same native magnetic problem.
    pub max_field_evaluations: usize,
    /// Central-difference step h [m]. The accepted derivative uses h/2 and is
    /// compared with h, including at independently checked interval midpoints.
    pub derivative_step_m: f64,
    /// Maximum difference between the two central derivatives [Wb/m].
    pub derivative_tolerance_wb_m: f64,
    /// Maximum midpoint interpolation error in permanent linkage [Wb].
    pub interpolation_tolerance_wb: f64,
    /// Maximum midpoint interpolation error in linkage derivative [Wb/m].
    pub interpolation_slope_tolerance_wb_m: f64,
    /// Original FEEC assembly dimensions and temporary-storage bounds.
    pub assembly: WeightedAssemblyLimits,
    /// Original winding storage, source-pairing and inverse-residual admission.
    pub winding: WindingConfig,
    /// Original fixed-remanence stationarity allowance [Wb].
    pub bias_stationarity_tolerance_wb: f64,
    /// Native checked-CG residual and total iteration budget per field solve.
    pub linear: CheckedCgConfig,
}

/// A refusal from the original physics owner or the explicit map resolution.
#[derive(Debug)]
pub enum PickupFieldError {
    /// Invalid geometry correspondence, source binding or numerical controls.
    Invalid(&'static str),
    /// No partially compiled map or changed string state is published.
    Cancelled,
    /// Existing magnetic material, assembly or residual refusal.
    Magnetic(MagneticError),
    /// Existing winding-current, energy-pairing or field refusal.
    Winding(WindingError),
    /// Existing bounded Faraday observer refusal.
    Pickup(PickupError),
    /// h and h/2 derivatives disagree beyond the declared allowance.
    Derivative {
        /// Physical sample gap [m].
        gap_m: f64,
        /// Measured absolute derivative disagreement [Wb/m].
        disagreement_wb_m: f64,
    },
    /// A fresh midpoint field solve disagrees with the Hermite map.
    Interpolation {
        /// Physical midpoint gap [m].
        gap_m: f64,
        /// Absolute permanent-linkage difference [Wb].
        linkage_error_wb: f64,
        /// Absolute derivative difference [Wb/m].
        slope_error_wb_m: f64,
    },
}
impl core::fmt::Display for PickupFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material pickup field: {self:?}")
    }
}
impl std::error::Error for PickupFieldError {}
impl From<MagneticError> for PickupFieldError {
    fn from(error: MagneticError) -> Self {
        match error {
            MagneticError::Cancelled => Self::Cancelled,
            error => Self::Magnetic(error),
        }
    }
}
impl From<WindingError> for PickupFieldError {
    fn from(error: WindingError) -> Self {
        match error {
            WindingError::Magnetic(MagneticError::Cancelled) => Self::Cancelled,
            error => Self::Winding(error),
        }
    }
}
impl From<PickupError> for PickupFieldError {
    fn from(error: PickupError) -> Self {
        Self::Pickup(error)
    }
}
fn invalid(message: &'static str) -> PickupFieldError {
    PickupFieldError::Invalid(message)
}
fn poll(cx: &Cx<'_>) -> Result<(), PickupFieldError> {
    cx.checkpoint().map_err(|_| PickupFieldError::Cancelled)
}
fn finite(value: f64) -> Result<f64, PickupFieldError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable pickup field response"))
    }
}

/// One actual native field evaluation, retained in evaluation order.
#[derive(Clone, Debug)]
pub struct PickupFieldSample {
    /// Gap used to move the original mesh [m].
    pub gap_m: f64,
    /// Signed zero-current winding linkage including turns and thickness [Wb].
    pub permanent_linkage_wb: f64,
    /// Native self-inductance, retained as field evidence; not applied to EMF.
    pub inductance_h: f64,
    /// True-residual report of the bias solve, if a magnet is present.
    pub bias_report: Option<SolveReport>,
    /// True-residual report of the 1-A response used by winding admission.
    pub unit_current_report: Option<SolveReport>,
}

/// An immutable field-derived linkage map with all original material sources.
#[derive(Clone, Debug)]
pub struct MaterialPickupMap {
    source: PlanarPickupSpec,
    map: PickupFluxMap,
    samples: Vec<PickupFieldSample>,
    maximum_derivative_disagreement_wb_m: f64,
    maximum_midpoint_linkage_error_wb: f64,
    maximum_midpoint_slope_error_wb_m: f64,
}
impl MaterialPickupMap {
    /// Solve source-defined fields at every knot and derivative stencil, then
    /// compare the resulting Hermite law with independent midpoint fields.
    /// The source, boundary, materials and coil profile remain identical across
    /// the sweep; only declared vertex positions move. No coefficient fitting,
    /// magnetic-name inference, extrapolation or new field solver is used.
    ///
    /// Every knot and midpoint consumes five field evaluations. All counts are
    /// checked before assembly. Derivatives require positive gaps at +/-h.
    /// Numerical comparisons are estimates at sampled positions; refine mesh,
    /// h and knots separately before claiming a physically resolved response.
    pub fn compile(
        cx: &Cx<'_>,
        source: PlanarPickupSpec,
        gaps_m: &[f64],
        config: PickupFieldConfig,
    ) -> Result<Self, PickupFieldError> {
        poll(cx)?;
        validate(&source, gaps_m, config)?;
        let count = 5 * (2 * gaps_m.len() - 1);
        let mut samples = Vec::with_capacity(count);
        let mut maximum_derivative_disagreement_wb_m = 0.0_f64;
        let mut knots = Vec::with_capacity(gaps_m.len());
        for &gap in gaps_m {
            let (knot, disagreement) = knot(cx, &source, gap, config, &mut samples)?;
            maximum_derivative_disagreement_wb_m =
                maximum_derivative_disagreement_wb_m.max(disagreement);
            knots.push(knot);
        }
        let map = PickupFluxMap::new(knots)?;
        let mut maximum_midpoint_linkage_error_wb = 0.0_f64;
        let mut maximum_midpoint_slope_error_wb_m = 0.0_f64;
        for pair in gaps_m.windows(2) {
            let gap = pair[0] + 0.5 * (pair[1] - pair[0]);
            let (reference, disagreement) = knot(cx, &source, gap, config, &mut samples)?;
            maximum_derivative_disagreement_wb_m =
                maximum_derivative_disagreement_wb_m.max(disagreement);
            let (linkage, slope) = map.sample(gap)?;
            let linkage_error_wb = finite(linkage - reference.linkage_wb)?.abs();
            let slope_error_wb_m = finite(slope - reference.slope_wb_per_m)?.abs();
            if linkage_error_wb > config.interpolation_tolerance_wb
                || slope_error_wb_m > config.interpolation_slope_tolerance_wb_m
            {
                return Err(PickupFieldError::Interpolation {
                    gap_m: gap,
                    linkage_error_wb,
                    slope_error_wb_m,
                });
            }
            maximum_midpoint_linkage_error_wb =
                maximum_midpoint_linkage_error_wb.max(linkage_error_wb);
            maximum_midpoint_slope_error_wb_m =
                maximum_midpoint_slope_error_wb_m.max(slope_error_wb_m);
        }
        poll(cx)?;
        Ok(Self {
            source,
            map,
            samples,
            maximum_derivative_disagreement_wb_m,
            maximum_midpoint_linkage_error_wb,
            maximum_midpoint_slope_error_wb_m,
        })
    }

    /// Exact frozen geometry, moving-cell assignment, material states and turns.
    #[must_use]
    pub const fn source(&self) -> &PlanarPickupSpec {
        &self.source
    }
    /// Every knot, stencil and independent midpoint field sample.
    #[must_use]
    pub fn samples(&self) -> &[PickupFieldSample] {
        &self.samples
    }
    /// Largest measured h versus h/2 derivative discrepancy [Wb/m].
    #[must_use]
    pub const fn maximum_derivative_disagreement_wb_m(&self) -> f64 {
        self.maximum_derivative_disagreement_wb_m
    }
    /// Largest fresh-midpoint linkage discrepancy [Wb].
    #[must_use]
    pub const fn maximum_midpoint_linkage_error_wb(&self) -> f64 {
        self.maximum_midpoint_linkage_error_wb
    }
    /// Largest fresh-midpoint slope discrepancy [Wb/m].
    #[must_use]
    pub const fn maximum_midpoint_slope_error_wb_m(&self) -> f64 {
        self.maximum_midpoint_slope_error_wb_m
    }
    /// Attach the original Faraday aperture operator in physical sinusoidal
    /// amplitude coordinates [m]. For a source-bound mass-normalized runtime,
    /// use `bind_material_string` so the source check and mass scale survive.
    pub fn bind(&self, pose: PickupPose, mode_count: usize) -> Result<Pickup, PickupFieldError> {
        if pose.height_m != self.source.reference_gap_m {
            return Err(invalid(
                "pickup rest height must equal the original field reference gap",
            ));
        }
        Ok(Pickup::bind_flux_map(pose, mode_count, self.map.clone())?)
    }

    /// Observe the existing homogeneous string without duplicating its state.
    /// The sensed magnetic region must carry the SAME immutable material card
    /// and query coordinates as the string and translate rigidly at unit speed
    /// per gap metre. Geometry remains the caller's admitted planar reduction,
    /// not an automatically certified reconstruction of the circular specimen.
    pub fn bind_material_string(
        self: &Arc<Self>,
        pose: PickupPose,
        string: &LinearMaterialStringRuntime,
    ) -> Result<MaterialStringPickup, PickupFieldError> {
        let source = string.specimen().material();
        let mut direction = None;
        for &cell in &self.source.sensed_cells {
            let PlanarMagneticMaterial::Linear(material) = &self.source.materials[cell] else {
                return Err(invalid(
                    "sensed string cells require induced magnetic material",
                ));
            };
            if !same_material_state(material, source) {
                return Err(invalid(
                    "pickup and string must share one material card and physical query state",
                ));
            }
            for &vertex in &self.source.mesh.faces()[cell] {
                let next = self.source.vertex_gap_directions[vertex as usize];
                if (det::hypot(next[0], next[1]) - 1.0).abs() > 1e-12
                    || direction.is_some_and(|prior| prior != next)
                {
                    return Err(invalid(
                        "material-string pickup needs one rigid unit gap translation of its sensed cells",
                    ));
                }
                direction = Some(next);
            }
        }
        let input = string.specimen().string();
        let mass = finite(0.5 * input.lin_density_kg_m * input.length_m)?;
        if mass <= 0.0 {
            return Err(invalid("material-string modal mass must be positive"));
        }
        let pickup = self.bind(pose, string.states().len())?;
        Ok(MaterialStringPickup {
            map: Arc::clone(self),
            pickup,
            specimen_identity: string.specimen().specimen_identity(),
            inverse_modal_mass_root: finite(1.0 / det::sqrt(mass))?,
        })
    }
}

fn same_material_state(a: &ResolvedMaterialStatePoint, b: &ResolvedMaterialStatePoint) -> bool {
    a.card_identity() == b.card_identity()
        && a.query_point() == b.query_point()
        && a.axis_quantities() == b.axis_quantities()
}

/// Open-circuit electrical observation of the original accepted string state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialStringPickupObservation {
    /// The string's actual accepted material/mechanical/gas epoch.
    pub epoch: u64,
    /// Window-weighted physical displacement toward the pole [m].
    pub displacement_m: f64,
    /// Window-weighted physical velocity toward the pole [m/s].
    pub velocity_m_s: f64,
    /// Field-derived signed linkage [Wb turns].
    pub linkage_wb: f64,
    /// Faraday voltage from this same motion [V].
    pub emf_v: f64,
}

/// Read-only material-string binding. It owns no vibration or thermal history.
#[derive(Clone, Debug)]
pub struct MaterialStringPickup {
    map: Arc<MaterialPickupMap>,
    pickup: Pickup,
    specimen_identity: ContentHash,
    inverse_modal_mass_root: f64,
}
impl MaterialStringPickup {
    /// Original material/field map, shared rather than rebuilt during playback.
    #[must_use]
    pub fn map(&self) -> &MaterialPickupMap {
        &self.map
    }
    /// Borrow the current owner's states, including its actual mass-normalized
    /// basis. Material/geometry replacement or basis refinement refuses until
    /// the map is recompiled/rebound; a gas-only update preserves this binding.
    /// No allocation or field solve is performed on this observation path.
    pub fn observe(
        &self,
        cx: &Cx<'_>,
        string: &LinearMaterialStringRuntime,
    ) -> Result<MaterialStringPickupObservation, PickupFieldError> {
        poll(cx)?;
        if string.specimen().specimen_identity() != self.specimen_identity
            || string.states().len() != self.pickup.gains.len()
        {
            return Err(invalid(
                "pickup field map is stale after string material, geometry or basis replacement",
            ));
        }
        let mut displacement = 0.0;
        let mut velocity = 0.0;
        for (index, (gain, state)) in self.pickup.gains.iter().zip(string.states()).enumerate() {
            if index % 64 == 0 {
                poll(cx)?;
            }
            let shape = gain * self.inverse_modal_mass_root;
            displacement += shape * state.displacement_m_sqrt_kg;
            velocity += shape * state.velocity_m_sqrt_kg_per_s;
        }
        let linkage_wb = self
            .pickup
            .linkage_and_gain_at_displacement(displacement)?
            .0;
        let emf_v = self.pickup.emf_at_displacement(displacement, velocity)?;
        poll(cx)?;
        Ok(MaterialStringPickupObservation {
            epoch: string.epoch(),
            displacement_m: displacement,
            velocity_m_s: velocity,
            linkage_wb,
            emf_v,
        })
    }
}

fn validate(
    source: &PlanarPickupSpec,
    gaps: &[f64],
    config: PickupFieldConfig,
) -> Result<(), PickupFieldError> {
    let count = gaps
        .len()
        .checked_mul(2)
        .and_then(|n| n.checked_sub(1))
        .and_then(|n| n.checked_mul(5));
    if gaps.len() < 2
        || count.is_none_or(|n| n > config.max_field_evaluations)
        || !source.reference_gap_m.is_finite()
        || source.reference_gap_m <= 0.0
        || gaps
            .first()
            .is_some_and(|gap| source.reference_gap_m < *gap)
        || gaps.last().is_some_and(|gap| source.reference_gap_m > *gap)
        || !config.derivative_step_m.is_finite()
        || config.derivative_step_m <= 0.0
        || gaps.iter().any(|gap| {
            !gap.is_finite()
                || *gap <= config.derivative_step_m
                || !(*gap + config.derivative_step_m).is_finite()
                || *gap - 0.5 * config.derivative_step_m == *gap
                || *gap + 0.5 * config.derivative_step_m == *gap
        })
        || gaps.windows(2).any(|p| p[0] >= p[1])
        || [
            config.derivative_tolerance_wb_m,
            config.interpolation_tolerance_wb,
            config.interpolation_slope_tolerance_wb_m,
        ]
        .iter()
        .any(|t| !t.is_finite() || *t < 0.0)
    {
        return Err(invalid(
            "pickup field knots, derivative resolution, tolerances or work budget are invalid",
        ));
    }
    if source.vertex_gap_directions.len() != source.mesh.vertices().len()
        || source.mesh.vertices().len() > config.assembly.max_dofs
        || source.mesh.edges().len() > config.assembly.max_dofs
        || source.mesh.faces().len() > config.assembly.max_cells
        || source
            .mesh
            .faces()
            .len()
            .checked_mul(9)
            .is_none_or(|n| n > config.assembly.max_triplets)
        || source
            .vertex_gap_directions
            .iter()
            .flatten()
            .any(|v| !v.is_finite())
        || source.materials.len() != source.mesh.faces().len()
        || source.winding_profile_per_m2.len() != source.mesh.faces().len()
        || source.winding_profile_per_m2.iter().any(|v| !v.is_finite())
        || source.sensed_cells.is_empty()
        || source.sensed_cells.iter().enumerate().any(|(i, &cell)| {
            cell >= source.materials.len() || source.sensed_cells[..i].contains(&cell)
        })
    {
        return Err(invalid(
            "pickup field geometry, material, winding and sensed-cell correspondences are incomplete",
        ));
    }
    for (cell, face) in source.mesh.faces().iter().enumerate() {
        let stationary = source.winding_profile_per_m2[cell] != 0.0
            || matches!(source.materials[cell], PlanarMagneticMaterial::Recoil(_));
        if stationary
            && face
                .iter()
                .any(|&v| source.vertex_gap_directions[v as usize] != [0.0; 2])
        {
            return Err(invalid(
                "pickup gap sweep must retain stationary magnet geometry and winding turns",
            ));
        }
    }
    if source
        .prescribed_vertices
        .iter()
        .any(|&v| v >= source.mesh.vertices().len() || source.vertex_gap_directions[v] != [0.0; 2])
    {
        return Err(invalid(
            "pickup prescribed field boundary must stay fixed throughout the gap sweep",
        ));
    }
    Ok(())
}

fn knot(
    cx: &Cx<'_>,
    source: &PlanarPickupSpec,
    gap: f64,
    config: PickupFieldConfig,
    samples: &mut Vec<PickupFieldSample>,
) -> Result<(PickupFluxKnot, f64), PickupFieldError> {
    let mut values = [0.0; 5];
    let h = config.derivative_step_m;
    for (value, delta) in values.iter_mut().zip([0.0, -h, h, -0.5 * h, 0.5 * h]) {
        let sample = evaluate(cx, source, gap + delta, config)?;
        *value = sample.permanent_linkage_wb;
        samples.push(sample);
    }
    let coarse = finite((values[2] - values[1]) / ((gap + h) - (gap - h)))?;
    let fine = finite((values[4] - values[3]) / ((gap + 0.5 * h) - (gap - 0.5 * h)))?;
    let disagreement = finite(coarse - fine)?.abs();
    if disagreement > config.derivative_tolerance_wb_m {
        return Err(PickupFieldError::Derivative {
            gap_m: gap,
            disagreement_wb_m: disagreement,
        });
    }
    Ok((
        PickupFluxKnot {
            gap_m: gap,
            linkage_wb: values[0],
            slope_wb_per_m: fine,
        },
        disagreement,
    ))
}

fn evaluate(
    cx: &Cx<'_>,
    source: &PlanarPickupSpec,
    gap: f64,
    config: PickupFieldConfig,
) -> Result<PickupFieldSample, PickupFieldError> {
    poll(cx)?;
    let displacement = finite(gap - source.reference_gap_m)?;
    let mut vertices = Vec::with_capacity(source.mesh.vertices().len());
    for (index, (point, direction)) in source
        .mesh
        .vertices()
        .iter()
        .zip(&source.vertex_gap_directions)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(cx)?;
        }
        vertices.push([
            finite(point[0] + displacement * direction[0])?,
            finite(point[1] + displacement * direction[1])?,
        ]);
    }
    let orientation = |vertices: &[[f64; 2]], face: [u32; 3]| {
        let [a, b, c] = face.map(|i| vertices[i as usize]);
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    };
    for (index, face) in source.mesh.faces().iter().enumerate() {
        if index % 64 == 0 {
            poll(cx)?;
        }
        let before = orientation(source.mesh.vertices(), *face);
        let after = finite(orientation(&vertices, *face))?;
        if after == 0.0 || before.is_sign_positive() != after.is_sign_positive() {
            return Err(invalid(
                "pickup gap sweep inverts or collapses an original triangle",
            ));
        }
    }
    let mesh = TriComplex2::from_triangles(
        source.mesh.lineage(),
        vertices,
        source.mesh.vertex_keys().to_vec(),
        source.mesh.faces().to_vec(),
        source.mesh.metric(),
    )
    .map_err(|_| invalid("pickup gap sweep leaves the original planar mesh domain"))?;
    let mut cancelled = || cx.checkpoint().is_err();
    let problem = PlanarMagnetostatics::new_with_material_models(
        mesh,
        source.materials.clone(),
        config.assembly,
        &mut cancelled,
    )?;
    let basis = PlanarWindingBasis::new_with_fixed_remanence(
        &problem,
        vec![source.winding_profile_per_m2.clone()],
        &source.prescribed_vertices,
        vec![source.end_inductance_h],
        config.winding,
        config.bias_stationarity_tolerance_wb,
        config.linear,
        &mut cancelled,
    )?;
    let sample = PickupFieldSample {
        gap_m: gap,
        permanent_linkage_wb: basis.bias_flux_linkages_wb()[0],
        inductance_h: basis.inductance_h()[0],
        bias_report: basis.bias_field().and_then(|field| field.report().cloned()),
        unit_current_report: basis.unit_fields()[0].report().cloned(),
    };
    poll(cx)?;
    Ok(sample)
}
