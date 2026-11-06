//! Scalar source closure at the actual tetrahedral electric field and current.
//!
//! Each bounded Picard iteration uses the original native bulk/contact operator.
//! Only potentials are combined inside a frozen terminal response; the entire
//! response is rebuilt after a conductivity change. Nonlinear fields cannot be
//! scaled to another drive. Optional coordinate relaxation requeries the same
//! original source; no coefficient blending or numerical heat is introduced.
//! Positive source values do not guarantee convergence or uniqueness.

use super::*;
use fs_material::conductor::{ConductorError, response::ScalarConductorResponse};

/// Constitutive and total native-solve limits, separate from field tolerances.
#[derive(Clone, Copy, Debug)]
pub struct SpatialConductionResponseControls {
    pub maximum_iterations: usize,
    /// Boundary drives cost one solve per iteration; N electrodes cost N-1.
    pub maximum_field_solves: usize,
    pub absolute_conductivity_tolerance_s_m: f64,
    pub relative_conductivity_tolerance: f64,
    pub current_density_tolerance_a_m2: f64,
}

#[derive(Debug)]
pub enum SpatialConductionResponseError {
    Field(SpatialConductionError),
    Material(ConductorError),
    Invalid(&'static str),
    Cancelled,
    FieldSolveBudget {
        completed: usize,
        required: usize,
        limit: usize,
    },
    NotConverged {
        iterations: usize,
        maximum_conductivity_defect_s_m: f64,
        maximum_current_density_defect_a_m2: f64,
        absolute_power_defect_w: f64,
    },
}
impl From<SpatialConductionError> for SpatialConductionResponseError {
    fn from(error: SpatialConductionError) -> Self {
        match error {
            SpatialConductionError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<ConductorError> for SpatialConductionResponseError {
    fn from(error: ConductorError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for SpatialConductionResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral conductor source response: {self:?}")
    }
}
impl std::error::Error for SpatialConductionResponseError {}
type Result<T> = core::result::Result<T, SpatialConductionResponseError>;

/// Actual-amplitude drive held throughout this stationary constitutive closure.
#[derive(Clone, Copy, Debug)]
pub enum SpatialConductionResponseDrive<'a> {
    Boundary(&'a SpatialConductionBoundary),
    /// Signed current INTO the selected electrodes, including zero for floating
    /// electrodes. Other electrodes prescribe voltage in original terminal order;
    /// currents follow the supplied current-terminal order.
    MixedTerminals {
        terminals: &'a [&'a [usize]],
        current_terminals: &'a [usize],
        voltages_v: &'a [f64],
        currents_a: &'a [f64],
        response: SpatialConductanceControls,
        mixed: SpatialMixedConductanceControls,
    },
}
impl<'a> From<&'a SpatialConductionBoundary> for SpatialConductionResponseDrive<'a> {
    fn from(boundary: &'a SpatialConductionBoundary) -> Self {
        Self::Boundary(boundary)
    }
}

/// Same immutable geometry and original source laws; no second current model.
/// The frozen operator stays private so a nonlinear operating point cannot be
/// accidentally exposed as a reusable linear resistance.
#[derive(Clone)]
pub struct SpatialConductionResponse {
    native: SpatialConduction,
    sources: Arc<[ScalarConductorResponse]>,
    relaxation: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpatialConductionResponseSample {
    pub tetrahedron: usize,
    /// Complete original source evaluated at the actual reported E and J.
    pub conductor: ResolvedConductorState,
    /// Evaluated minus applied conductivity.
    pub conductivity_defect_s_m: f64,
    /// Native J minus source-evaluated sigma(E,J) E, in original mesh axes.
    pub current_density_defect_a_m2: [f64; 3],
    pub power_defect_w: f64,
}

#[derive(Clone, Debug)]
pub struct SpatialConductionResponseReport {
    /// Original native field, including separately retained contact Joule heat.
    pub field: SpatialConductionField,
    pub evaluated: Vec<SpatialConductionResponseSample>,
    pub iterations: usize,
    /// Includes every basis solve in discarded constitutive iterations.
    pub field_solves: usize,
    /// Actual-field and damped-coordinate source queries, excluding initial
    /// binding. Full steps and fixed sources need no extra damping queries.
    pub source_evaluations: usize,
}

pub struct SpatialConductionResponseCandidate {
    problem: SpatialConductionResponse,
    report: SpatialConductionResponseReport,
}
impl SpatialConductionResponseCandidate {
    pub fn problem(&self) -> &SpatialConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &SpatialConductionResponseReport {
        &self.report
    }
    pub fn into_parts(self) -> (SpatialConductionResponse, SpatialConductionResponseReport) {
        (self.problem, self.report)
    }
}

impl SpatialConductionResponse {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        sources: Vec<ScalarConductorResponse>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_contacts(
            positions_m,
            tetrahedra,
            sources,
            Vec::new(),
            limits,
            cancelled,
        )
    }

    /// One complete scalar source per tetrahedron. The original contact laws
    /// retain fixed external states; this owner closes bulk E/J coordinates only.
    /// Temperature/phase contexts are prescribed here, not advanced thermally.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contacts(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        sources: Vec<ScalarConductorResponse>,
        contacts: Vec<SpatialConductionContact>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if sources.len() != tetrahedra.len()
            || sources.len() > limits.max_cells
            || contacts.len() > limits.max_triplets / 36
        {
            return Err(SpatialConductionResponseError::Invalid(
                "complete bounded source/contact coverage required",
            ));
        }
        let mut conductors = Vec::with_capacity(sources.len());
        for source in &sources {
            check(cancelled)?;
            source.admit(source.original())?;
            conductors.push(source.original().clone());
        }
        let geometry = Arc::new(
            prepare_geometry(positions_m, tetrahedra, limits, cancelled)
                .map_err(SpatialConductionError::from)?,
        );
        // The original laws above own EVERY E/J coordinate. This private
        // assembly seam never relaxes admission for ordinary Ohmic snapshots.
        let contact_axes = vec![[None, None]; contacts.len()];
        let native =
            SpatialConduction::assemble_admitted(geometry, conductors, limits, cancelled, |_| {
                true
            })?
            .bind_contacts(contacts, contact_axes, limits, cancelled)?;
        Ok(Self {
            native,
            sources: sources.into(),
            relaxation: 1.0,
        })
    }

    /// Select a fixed source-coordinate relaxation weight in `(0, 1]`.
    /// One preserves the original undamped Picard path. Smaller weights can
    /// stabilize alternating resistance responses by requerying original laws
    /// at interpolated owned E/J magnitudes. Temperature/phase, source pins and
    /// all actual-field current/energy acceptance gates remain unchanged.
    /// No global convergence or uniqueness claim follows from relaxation.
    pub fn with_relaxation(mut self, weight: f64) -> Result<Self> {
        if !(super::terminals::positive(weight) && weight <= 1.0) {
            return Err(SpatialConductionResponseError::Invalid(
                "source-coordinate relaxation must be in (0, 1]",
            ));
        }
        self.relaxation = weight;
        Ok(self)
    }

    /// Numerical weight retained across accepted source-response candidates.
    pub fn relaxation(&self) -> f64 {
        self.relaxation
    }

    pub fn mesh(&self) -> &TetComplex {
        self.native.mesh()
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        self.native.positions_m()
    }
    pub fn exterior_faces(&self) -> impl Iterator<Item = usize> + '_ {
        self.native.exterior_faces()
    }
    pub fn conductors(&self) -> &[ResolvedConductorState] {
        self.native.conductors()
    }
    pub fn sources(&self) -> &[ScalarConductorResponse] {
        &self.sources
    }
    pub fn contacts(&self) -> &[SpatialConductionContact] {
        self.native.contacts()
    }

    /// Close the original sources against actual E/J at the authored drive.
    /// Native current, Hodge and supplied-work gates remain active. Separate
    /// sigma, current-density and sum-ABSOLUTE-local-power gates admit the
    /// constitutive response; errors in different cells cannot cancel.
    ///
    /// Geometry is reused. Every iteration rebuilds native conductivity/contact
    /// assembly and terminal response, charging all N-1 basis solves before work.
    /// Source bounds, nonconvergence and cancellation publish no candidate.
    /// This is a bounded Picard closure, not a claim of global convergence,
    /// uniqueness, field-dependent contact physics or a transient thermal solve.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_candidate<'d>(
        &self,
        drive: impl Into<SpatialConductionResponseDrive<'d>>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: SpatialConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionResponseCandidate> {
        check(cancelled)?;
        if !(1..=crate::conduction::response::MAX_RESPONSE_ITERATIONS)
            .contains(&controls.maximum_iterations)
            || controls.maximum_field_solves == 0
            || !controls.absolute_conductivity_tolerance_s_m.is_finite()
            || controls.absolute_conductivity_tolerance_s_m < 0.0
            || !controls.relative_conductivity_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_conductivity_tolerance)
            || (controls.absolute_conductivity_tolerance_s_m == 0.0
                && controls.relative_conductivity_tolerance == 0.0)
            || !super::terminals::positive(controls.current_density_tolerance_a_m2)
        {
            return Err(SpatialConductionResponseError::Invalid(
                "invalid bounded source response controls",
            ));
        }
        let drive = drive.into();
        let solves_per_iteration = match drive {
            SpatialConductionResponseDrive::Boundary(_) => 1,
            SpatialConductionResponseDrive::MixedTerminals { terminals, .. } => {
                terminals.len().checked_sub(1).filter(|&n| n > 0).ok_or(
                    SpatialConductionResponseError::Invalid(
                        "at least two mixed-drive electrodes required",
                    ),
                )?
            }
        };
        let seeds = self.native.conductors();
        let mut applied = seeds.to_vec();
        let mut field_solves = 0usize;
        let mut relaxation_queries = 0usize;
        let (mut maximum_sigma, mut maximum_current, mut absolute_power) =
            (0.0_f64, 0.0_f64, 0.0_f64);
        for iteration in 1..=controls.maximum_iterations {
            check(cancelled)?;
            let needed = field_solves.checked_add(solves_per_iteration);
            if needed.is_none_or(|n| n > controls.maximum_field_solves) {
                return Err(SpatialConductionResponseError::FieldSolveBudget {
                    completed: field_solves,
                    required: solves_per_iteration,
                    limit: controls.maximum_field_solves,
                });
            }
            // All applied bundles come from their original admitted source
            // (construction or the preceding actual-field query).
            let native = SpatialConduction::assemble_admitted(
                Arc::clone(&self.native.geometry),
                applied,
                assembly,
                cancelled,
                |_| true,
            )?
            .bind_contacts(
                self.native.contacts.clone(),
                self.native.contact_liquid_mass_fraction_axes.clone(),
                assembly,
                cancelled,
            )?;
            let field = match drive {
                SpatialConductionResponseDrive::Boundary(boundary) => {
                    native.solve(boundary, field_controls, cancelled)?
                }
                SpatialConductionResponseDrive::MixedTerminals {
                    terminals,
                    current_terminals,
                    voltages_v,
                    currents_a,
                    response,
                    mixed,
                } => {
                    let matrix = native.conductance_matrix(
                        terminals,
                        response,
                        field_controls,
                        cancelled,
                    )?;
                    let mixed =
                        matrix.prepare_current_control(current_terminals, mixed, cancelled)?;
                    mixed.field_at_controls(
                        voltages_v,
                        currents_a,
                        field_controls.current_tolerance_a,
                        field_controls.power_tolerance_w,
                        cancelled,
                    )?
                }
            };
            field_solves += solves_per_iteration;
            let mut evaluated = Vec::with_capacity(seeds.len());
            let mut next = Vec::with_capacity(seeds.len());
            let mut converged = true;
            maximum_sigma = 0.0;
            maximum_current = 0.0;
            absolute_power = 0.0;
            let mut power_correction = 0.0;
            for (tetrahedron, (source, seed)) in self.sources.iter().zip(seeds).enumerate() {
                check(cancelled)?;
                let e = field.electric_field_v_m[tetrahedron];
                let j = field.current_density_a_m2[tetrahedron];
                let conductor = source.at_spatial_fields(seed, e, j)?;
                let sigma = finite(conductor.resistivity_ohm_m().recip())?;
                let used = native.conductors[tetrahedron].resistivity_ohm_m().recip();
                let conductivity_defect_s_m = finite(sigma - used)?;
                let defect = [
                    finite(j[0] - sigma * e[0])?,
                    finite(j[1] - sigma * e[1])?,
                    finite(j[2] - sigma * e[2])?,
                ];
                let current_defect = finite(defect[0].hypot(defect[1]).hypot(defect[2]))?;
                let volume = native.cell_volume_m3(tetrahedron).ok_or(
                    SpatialConductionResponseError::Invalid("original cell volume absent"),
                )?;
                let power_defect_w = finite(volume * dot(e, defect))?;
                let tolerance = finite(
                    controls.absolute_conductivity_tolerance_s_m
                        + controls.relative_conductivity_tolerance * sigma.max(used),
                )?;
                maximum_sigma = maximum_sigma.max(conductivity_defect_s_m.abs());
                maximum_current = maximum_current.max(current_defect);
                let adjusted = finite(power_defect_w.abs() - power_correction)?;
                let sum = finite(absolute_power + adjusted)?;
                power_correction = finite((sum - absolute_power) - adjusted)?;
                absolute_power = sum;
                converged &= conductivity_defect_s_m.abs() <= tolerance
                    && current_defect <= controls.current_density_tolerance_a_m2;
                next.push(conductor.clone());
                evaluated.push(SpatialConductionResponseSample {
                    tetrahedron,
                    conductor,
                    conductivity_defect_s_m,
                    current_density_defect_a_m2: defect,
                    power_defect_w,
                });
            }
            check(cancelled)?;
            if converged && absolute_power <= field_controls.power_tolerance_w {
                let report = SpatialConductionResponseReport {
                    field,
                    evaluated,
                    iterations: iteration,
                    field_solves,
                    source_evaluations: iteration
                        .checked_mul(seeds.len())
                        .and_then(|count| count.checked_add(relaxation_queries))
                        .ok_or(SpatialConductionResponseError::Invalid(
                            "source evaluation count overflows",
                        ))?,
                };
                return Ok(SpatialConductionResponseCandidate {
                    problem: Self {
                        native,
                        sources: self.sources.clone(),
                        relaxation: self.relaxation,
                    },
                    report,
                });
            }
            if iteration < controls.maximum_iterations && self.relaxation < 1.0 {
                for ((source, current), following) in
                    self.sources.iter().zip(native.conductors()).zip(&mut next)
                {
                    check(cancelled)?;
                    if !source.axes().is_empty() {
                        *following =
                            source.relax_operating_point(current, following, self.relaxation)?;
                        relaxation_queries = relaxation_queries.checked_add(1).ok_or(
                            SpatialConductionResponseError::Invalid(
                                "source evaluation count overflows",
                            ),
                        )?;
                    }
                }
            }
            applied = next;
        }
        Err(SpatialConductionResponseError::NotConverged {
            iterations: controls.maximum_iterations,
            maximum_conductivity_defect_s_m: maximum_sigma,
            maximum_current_density_defect_a_m2: maximum_current,
            absolute_power_defect_w: absolute_power,
        })
    }
}
