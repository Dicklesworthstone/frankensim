//! Dry-basis moisture transport and source-bound mechanical consumers.
//!
//! A homogeneous, isothermal slab stores kg water/kg dry solid in reference
//! finite volumes. Fick transport uses the selected material diffusivity;
//! backward Euler uses the existing checked-CG owner by default. Optional
//! implicit diffusivity resolves the original source at the nonlinear endpoint
//! with bounded Picard iteration and the native tridiagonal solver. Otherwise
//! diffusivity is lagged at the accepted start, then re-evaluated at every
//! proposed endpoint. No clipped water, inferred wood diffusivity,
//! instantaneous bulk equilibration or second material evaluator is used.
//!
//! Transport uses fixed reference dry-skeleton volumes and dilute bound moisture.
//! The solid consumer resolves quasistatic swelling in that reference chart.
//! The slab has no transverse transport; `cartesian` supplies three-dimensional
//! principal-axis diffusion and the same source-bound native solid mechanics.
//! Sorption heat, stress-assisted diffusion, free liquid and hysteresis are not
//! modeled. Ordinary surface conditions are
//! infinite reservoirs; `atmosphere` instead balances an actual finite vapor
//! inventory against this same diffusion owner. The plate consumer retains the
//! fixed geometry and rebuilds actual stiffness AND wet inertia.

/// Whole-duration temporal accuracy around the original conservative endpoints.
pub mod adaptive;
/// Finite vapor inventory, source sorption and the original diffusion endpoint.
pub mod atmosphere;
mod binding;
/// Three-dimensional principal-axis transport on actual Cartesian finite volumes.
pub mod cartesian;
mod implicit;
mod reservoir;
pub mod plate;
/// Directional swelling, compatible shape, stress and native vibration tangent.
pub mod solid;
/// Finite aqueous solvent inventory and activity-controlled bound-water exchange.
pub mod solution;
/// Source-backed finite surface exchange in the dry-basis flux convention.
pub mod surface;
pub use implicit::{MoistureNonlinearConfig, MoistureNonlinearReport};

use surface::ResolvedMoistureFace;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::{
    MaterialError,
    state_point::{
        MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::{
    Density, Dims, QuantitySpec,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
};
use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, SolveReport, checked_cg};
use fs_sparse::precond::Precond;

/// Dimensions of effective moisture diffusivity [m²/s].
pub const MOISTURE_DIFFUSIVITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);
/// Hard per-solve bound, additionally narrowed by the caller's budget.
pub const MAX_MOISTURE_CELLS: usize = 65_536;

/// Original source/solver/consumer refusal. Failed proposals publish no state.
#[derive(Debug)]
pub enum MoistureError {
    /// Native source-backed stiffness or solid material admission refused.
    SolidMaterial(fs_solid::linear3::source::SolidMaterialError),
    /// Native eigenstrain, equilibrium, geometry or prestress admission refused.
    Solid(fs_solid::linear3::TetElasticError),
    /// Malformed geometry, quantities, source applicability or limits.
    Invalid(&'static str),
    /// Shared material resolver refused.
    Material(MaterialStatePointError),
    /// Shared constitutive source refused.
    Constitutive(MaterialError),
    /// Native DKT admission failed.
    Plate(fs_plate::PlateError),
    /// Existing plate material lowering refused.
    Consumer(crate::acoustic_realize::AcousticRealizeError),
    /// Native linear solver did not converge.
    NotConverged(SolveReport),
    /// Native tridiagonal factorization or substitution refused.
    Tridiagonal(fs_la::tridiagonal::TridiagonalError),
    /// Endpoint-dependent diffusion failed its nonlinear residual or water budget.
    NonlinearNotConverged(MoistureNonlinearReport),
    /// Spatial endpoint-dependent transport exhausted its explicit work budget.
    SpatialNotConverged {
        iterations: usize,
        linear_iterations: usize,
        relative_residual: f64,
        water_balance_residual_kg: f64,
    },
    /// A valid cell endpoint exceeded the caller's numerical change limit.
    ContentIncrement {
        cell: usize,
        increment: f64,
        maximum: f64,
    },
    /// Sum of stored water did not match the two measured boundary transfers.
    MassBalance { residual_kg: f64, tolerance_kg: f64 },
    /// No endpoint is published after cancellation.
    Cancelled,
}
impl core::fmt::Display for MoistureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moisture transport: {self:?}")
    }
}
impl std::error::Error for MoistureError {}
impl From<MaterialStatePointError> for MoistureError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<MaterialError> for MoistureError {
    fn from(e: MaterialError) -> Self {
        Self::Constitutive(e)
    }
}
impl From<fs_solid::linear3::source::SolidMaterialError> for MoistureError {
    fn from(e: fs_solid::linear3::source::SolidMaterialError) -> Self {
        match e {
            fs_solid::linear3::source::SolidMaterialError::Solid(
                fs_solid::linear3::TetElasticError::Cancelled,
            ) => Self::Cancelled,
            e => Self::SolidMaterial(e),
        }
    }
}
impl From<fs_solid::linear3::TetElasticError> for MoistureError {
    fn from(e: fs_solid::linear3::TetElasticError) -> Self {
        match e {
            fs_solid::linear3::TetElasticError::Cancelled => Self::Cancelled,
            e => Self::Solid(e),
        }
    }
}
impl From<fs_plate::PlateError> for MoistureError {
    fn from(e: fs_plate::PlateError) -> Self {
        Self::Plate(e)
    }
}
impl From<crate::acoustic_realize::AcousticRealizeError> for MoistureError {
    fn from(e: crate::acoustic_realize::AcousticRealizeError) -> Self {
        Self::Consumer(e)
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll(cx: &Cx<'_>) -> Result<(), MoistureError> {
    cx.checkpoint().map_err(|_| MoistureError::Cancelled)
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

/// Surface moisture is distinct from ambient RH and from interior moisture.
#[derive(Clone, Debug, PartialEq)]
pub enum MoistureFace {
    /// Zero normal water flux.
    Sealed,
    /// Explicit reservoir value in kg water/kg dry solid.
    PrescribedDryBasis(f64),
    /// Pinned source sorption at prescribed ambient RH. A surface transfer law
    /// attached to `law` limits uptake/desorption; otherwise the surface is at
    /// equilibrium. This boundary is an infinite isothermal vapor reservoir.
    SourceSorption {
        /// Original equilibrium isotherm and optional surface-emission source.
        law: atmosphere::DryBasisSorption,
        /// Prescribed environmental vapor activity in [0,1].
        relative_humidity: f64,
    },
    /// Source-evaluated FPL sorption boundary at a declared T and RH.
    WoodEquilibrium {
        /// Immutable coefficient card consumed by the original FPL law.
        card: MaterialCard,
        /// Absolute environment temperature [K]; must match the slab.
        temperature_k: f64,
        /// Vapor partial-pressure ratio, never specimen mass fraction.
        relative_humidity: f64,
    },
}
impl MoistureFace {
    fn resolve(&self, temperature: f64) -> Result<ResolvedMoistureFace, MoistureError> {
        let value = match self {
            Self::Sealed => return Ok(ResolvedMoistureFace::default()),
            Self::PrescribedDryBasis(value) => *value,
            Self::SourceSorption {
                law,
                relative_humidity,
            } => {
                return law.resolve(temperature, *relative_humidity);
            }
            Self::WoodEquilibrium {
                card,
                temperature_k,
                relative_humidity,
            } => {
                if temperature_k.to_bits() != temperature.to_bits()
                    || !relative_humidity.is_finite()
                    || !(0.0..=1.0).contains(relative_humidity)
                {
                    return Err(MoistureError::Invalid(
                        "isothermal source boundary needs matching T and RH in [0,1]",
                    ));
                }
                fs_material::wood_moisture::equilibrium_moisture_content(
                    card,
                    *temperature_k,
                    *relative_humidity,
                )?
            }
        };
        if !value.is_finite() || value < 0.0 {
            return Err(MoistureError::Invalid(
                "surface moisture must be finite nonnegative dry-basis content",
            ));
        }
        Ok(ResolvedMoistureFace {
            equilibrium: Some(value),
            ..ResolvedMoistureFace::default()
        })
    }
}

/// Geometry and physical applicability of the fixed dry skeleton.
#[derive(Clone, Copy, Debug)]
pub struct MoistureSlabGeometry {
    /// Diffusion length [m]; plate consumption maps it to the in-plane x axis.
    pub length_m: f64,
    /// Cross-sectional area normal to transport [m²].
    pub cross_section_m2: f64,
    /// Maximum supported bound-water content, kg water/kg dry solid.
    /// This is an authored applicability limit, not an inferred saturation law.
    pub maximum_moisture: f64,
    /// Allocation and per-operator work limit.
    pub maximum_cells: usize,
}

/// Explicit linear work, conservation and step-admission limits.
#[derive(Clone, Copy, Debug)]
pub struct MoistureStepConfig {
    /// Native checked-CG controls for the lagged slab and Cartesian schemes.
    /// In the slab's implicit mode,
    /// tolerance is the true linear residual target, max_iterations bounds
    /// tridiagonal solves per Picard iteration, and max_corrections bounds
    /// additional defect solves using the same factorization.
    pub linear: CheckedCgConfig,
    /// Absolute allowed water balance residual [kg].
    pub water_balance_tolerance_kg: f64,
    /// Largest accepted cell-content change [kg water/kg dry solid].
    /// This bounds a proposal; it is not an estimated time-discretization error.
    pub maximum_content_increment: f64,
}

/// One accepted moisture field. Clone is a complete physical checkpoint.
#[derive(Clone, Debug)]
pub struct MoistureSlab {
    card: MaterialCard,
    reference: ResolvedMaterialStatePoint,
    transport: ResolvedMaterialStatePoint,
    diffusivity_property: String,
    geometry: MoistureSlabGeometry,
    dry_density_kg_m3: f64,
    dry_mass_per_cell_kg: f64,
    moisture: Vec<f64>,
    time_s: f64,
    identity: ContentHash,
    implicit_diffusivity: Option<MoistureNonlinearConfig>,
}

/// Actual source evaluations and the water transferred over one interval.
#[derive(Clone, Debug)]
pub struct MoistureStepReport {
    /// Start and accepted endpoint time [s].
    pub interval_s: [f64; 2],
    /// Signed water received through the left and right prescribed surfaces
    /// [kg]. The finite-atmosphere wrapper debits that same external inventory.
    pub boundary_water_kg: [f64; 2],
    /// Actual surface dry-basis content. Finite exchange permits a difference
    /// from reservoir equilibrium; None denotes an explicitly sealed face.
    pub surface_moisture: [Option<f64>; 2],
    /// Reservoir equilibrium content; unchanged from the prescribed value in
    /// the instantaneous surface-equilibrium limit.
    pub reservoir_moisture: [Option<f64>; 2],
    /// Pinned source sorption at the actual environmental RH, when used.
    pub surface_sorption: [Option<ResolvedMaterialStatePoint>; 2],
    /// Pinned surface-emission source at the same RH, when exchange is finite.
    pub surface_transfer: [Option<ResolvedMaterialStatePoint>; 2],
    /// Water stored before and after [kg].
    pub stored_water_kg: [f64; 2],
    /// Stored change minus external transfer [kg].
    pub water_balance_residual_kg: f64,
    /// Complete accepted-start material evaluations. The default scheme uses
    /// these coefficients; implicit mode uses the accepted endpoint sources.
    pub transport_start: Vec<ResolvedMaterialStatePoint>,
    /// All endpoint sources re-admitted before publication.
    pub transport_end: Vec<ResolvedMaterialStatePoint>,
    /// Original true-residual solver report.
    pub linear: SolveReport,
    /// Endpoint nonlinear residual/work, or None for accepted-start lagging.
    /// `linear` describes the last frozen-coefficient solve; this report tests
    /// the complete endpoint-dependent equations independently of that solve.
    pub nonlinear: Option<MoistureNonlinearReport>,
}

impl MoistureSlab {
    /// Bind a homogeneous source and a finite, possibly nonuniform initial field.
    ///
    /// The reference query must explicitly type temperature as absolute K and
    /// moisture_content as dry-mass moisture. Reference density is measured per
    /// reference volume at that moisture; rho_dry = rho_ref/(1+M_ref). The density
    /// claim is never re-queried to create/destroy dry skeleton or water.
    /// Diffusivity may depend on the actual local moisture and fixed environment.
    /// Selection pins density and the named diffusivity at the reference state.
    pub fn new(
        card: MaterialCard,
        reference_point: QueryPoint,
        diffusivity_property: &str,
        selection: MaterialPropertySelection,
        geometry: MoistureSlabGeometry,
        initial_moisture: Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        poll(cx)?;
        let n = initial_moisture.len();
        if n == 0
            || n > geometry.maximum_cells
            || n > MAX_MOISTURE_CELLS
            || !positive(geometry.length_m)
            || !positive(geometry.cross_section_m2)
            || !positive(geometry.maximum_moisture)
            || diffusivity_property == "density"
        {
            return Err(MoistureError::Invalid(
                "bounded nonempty cells, positive geometry/domain and a distinct diffusivity key required",
            ));
        }
        if initial_moisture.iter()
            .any(|m| !m.is_finite() || *m < 0.0 || *m > geometry.maximum_moisture)
        {
            return Err(MoistureError::Invalid(
                "initial/reference moisture outside declared bound-water domain",
            ));
        }
        let (reference, transport, dry_density_kg_m3) = binding::bind_source(
            &card, &reference_point, &[diffusivity_property], selection,
            geometry.maximum_moisture,
        )?;
        let dry_mass_per_cell_kg =
            dry_density_kg_m3 * geometry.cross_section_m2 * (geometry.length_m / n as f64);
        if !positive(dry_density_kg_m3)
            || !positive(dry_mass_per_cell_kg)
            || !positive(dry_mass_per_cell_kg * n as f64)
        {
            return Err(MoistureError::Invalid(
                "dry reference inventory is unrepresentable",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-couple.moisture-slab.v1");
        hash.update(reference.identity().as_bytes());
        for v in [
            geometry.length_m,
            geometry.cross_section_m2,
            geometry.maximum_moisture,
            dry_mass_per_cell_kg,
        ] {
            hash.update(&v.to_bits().to_le_bytes());
        }
        hash.update(&(n as u64).to_le_bytes());
        let slab = Self {
            card,
            reference,
            transport,
            diffusivity_property: diffusivity_property.into(),
            geometry,
            dry_density_kg_m3,
            dry_mass_per_cell_kg,
            moisture: initial_moisture,
            time_s: 0.0,
            identity: hash.finalize(),
            implicit_diffusivity: None,
        };
        if !positive(slab.wet_mass_kg()) || !slab.water_mass_kg().is_finite() {
            return Err(MoistureError::Invalid(
                "initial wet inventory is unrepresentable",
            ));
        }
        slab.transport_states(&slab.moisture, cx)?;
        Ok(slab)
    }

    /// Immutable material/geometry/reference-inventory identity.
    pub const fn definition_identity(&self) -> ContentHash {
        self.identity
    }
    /// Accepted local kg water/kg dry solid.
    pub fn moisture(&self) -> &[f64] {
        &self.moisture
    }
    /// Accepted physical time [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Conserved dry skeleton mass [kg].
    pub fn dry_mass_kg(&self) -> f64 {
        self.dry_mass_per_cell_kg * self.moisture.len() as f64
    }
    /// Actual water inventory [kg].
    pub fn water_mass_kg(&self) -> f64 {
        self.moisture
            .iter()
            .map(|m| m * self.dry_mass_per_cell_kg)
            .sum()
    }
    /// Wet physical mass = dry skeleton plus explicitly transported water [kg].
    pub fn wet_mass_kg(&self) -> f64 {
        self.dry_mass_kg() + self.water_mass_kg()
    }
    /// Immutable shared material source.
    pub const fn card(&self) -> &MaterialCard {
        &self.card
    }
    /// Reference material resolution, including the actual reference wet density.
    pub const fn reference(&self) -> &ResolvedMaterialStatePoint {
        &self.reference
    }

    fn field(&self) -> binding::MoistureField<'_> {
        binding::MoistureField {
            card: &self.card,
            reference: &self.reference,
            moisture: &self.moisture,
            dry_mass_per_cell_kg: self.dry_mass_per_cell_kg,
            time_s: self.time_s,
            identity: self.identity,
        }
    }

    /// Resolve moisture-dependent diffusivity at each backward-Euler endpoint.
    /// The choice propagates through existing finite-atmosphere, plate and solid
    /// advances because they all consume this same slab. It changes numerical
    /// integration, not the material/geometry identity or accepted physical state.
    pub fn with_implicit_diffusivity(
        mut self,
        config: MoistureNonlinearConfig,
    ) -> Result<Self, MoistureError> {
        config.validate()?;
        self.implicit_diffusivity = Some(config);
        Ok(self)
    }

    fn point_at(&self, moisture: f64) -> Result<QueryPoint, MoistureError> {
        self.transport
            .query()
            .with("moisture_content", moisture)
            .map_err(|_| MoistureError::Invalid("unrepresentable dry-basis moisture coordinate"))
    }
    fn transport_states(
        &self,
        moisture: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<ResolvedMaterialStatePoint>, MoistureError> {
        moisture
            .iter()
            .map(|&m| {
                poll(cx)?;
                Ok(self
                    .transport
                    .at_query_point(&self.card, &self.point_at(m)?)?)
            })
            .collect()
    }

    /// Propose one conservative interval without mutating the accepted slab.
    /// Face values are constant over this interval; cut intervals at RH events.
    /// Interior face diffusivity is the harmonic mean of adjacent cell values,
    /// at the accepted start by default or endpoint in implicit mode. The
    /// boundary lies half a cell from its cell center; optional
    /// finite surface exchange adds its resistance in series. Numerical
    /// coefficients are refreshed on EVERY call; no stale material operator is reused.
    pub fn advance(
        &self,
        dt_s: f64,
        faces: &[MoistureFace; 2],
        config: MoistureStepConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureStepReport), MoistureError> {
        poll(cx)?;
        let temperature = self
            .reference
            .query_point()
            .iter()
            .find(|(name, _)| name == "temperature")
            .expect("admitted temperature")
            .1;
        let boundary = [
            faces[0].resolve(temperature)?,
            faces[1].resolve(temperature)?,
        ];
        self.advance_resolved(dt_s, boundary, config, cx)
    }

    fn advance_resolved(
        &self,
        dt_s: f64,
        faces: [ResolvedMoistureFace; 2],
        config: MoistureStepConfig,
        cx: &Cx<'_>,
    ) -> Result<(Self, MoistureStepReport), MoistureError> {
        poll(cx)?;
        let end = self.time_s + dt_s;
        if !positive(dt_s)
            || !end.is_finite()
            || end <= self.time_s
            || !positive(config.water_balance_tolerance_kg)
            || !positive(config.maximum_content_increment)
        {
            return Err(MoistureError::Invalid(
                "representable forward interval and positive admission tolerances required",
            ));
        }
        let surface = [faces[0].equilibrium, faces[1].equilibrium];
        for m in surface.into_iter().flatten() {
            if m > self.geometry.maximum_moisture {
                return Err(MoistureError::Invalid(
                    "surface exceeds the bound-water applicability domain",
                ));
            }
            // Admission at the actual interface state, not just the interior.
            self.transport
                .at_query_point(&self.card, &self.point_at(m)?)?;
        }
        let start = self.transport_states(&self.moisture, cx)?;
        let (solution, endpoint, coupling, nonlinear) =
            if let Some(controls) = self.implicit_diffusivity {
                self.solve_implicit(dt_s, &faces, &start, config, controls, cx)?
            } else {
                let (operator, rhs, coupling) = self.transport_system(dt_s, &faces, &start, cx)?;
                let solution = checked_cg(&operator, &operator, &rhs, config.linear, |_| poll(cx))
                    .map_err(|e| match e {
                        CheckedCgError::Interrupted(error) => error,
                        CheckedCgError::InvalidInput(what) => MoistureError::Invalid(what),
                    })?;
                if !solution.report.converged_euclidean() {
                    return Err(MoistureError::NotConverged(solution.report));
                }
                let endpoint = self.transport_states(&solution.x, cx)?;
                (solution, endpoint, coupling, None)
            };
        let n = self.moisture.len();
        let lower = self
            .moisture
            .iter()
            .copied()
            .chain(surface.into_iter().flatten())
            .fold(f64::INFINITY, f64::min);
        let upper = self
            .moisture
            .iter()
            .copied()
            .chain(surface.into_iter().flatten())
            .fold(0.0, f64::max);
        let roundoff = 128.0 * f64::EPSILON * upper.max(1.0);
        for (cell, (&old, &new)) in self.moisture.iter().zip(&solution.x).enumerate() {
            if !new.is_finite()
                || new < 0.0
                || new > self.geometry.maximum_moisture
                || new < lower - roundoff
                || new > upper + roundoff
            {
                return Err(MoistureError::Invalid(
                    "candidate violates moisture maximum principle or domain",
                ));
            }
            let increment = (new - old).abs();
            if increment > config.maximum_content_increment {
                return Err(MoistureError::ContentIncrement {
                    cell,
                    increment,
                    maximum: config.maximum_content_increment,
                });
            }
        }
        let mut actual_surface = [None, None];
        for (side, i) in [0, n - 1].into_iter().enumerate() {
            poll(cx)?;
            if let Some(equilibrium) = surface[side] {
                let fraction = coupling[side].1;
                let value = if fraction == 1.0 {
                    equilibrium
                } else {
                    fraction.mul_add(equilibrium - solution.x[i], solution.x[i])
                };
                // The real surface, not just the reservoir and cell center,
                // must be inside the original transport source domain.
                self.transport
                    .at_query_point(&self.card, &self.point_at(value)?)?;
                actual_surface[side] = Some(value);
            }
        }
        let transfer = self.boundary_transfer(&solution.x, &faces, &coupling);
        let water_before = self.water_mass_kg();
        let water_after: f64 = solution
            .x
            .iter()
            .map(|m| self.dry_mass_per_cell_kg * m)
            .sum();
        if !positive(self.dry_mass_kg() + water_after) {
            return Err(MoistureError::Invalid(
                "proposed wet inventory is unrepresentable",
            ));
        }
        let residual = water_after - water_before - transfer[0] - transfer[1];
        if !residual.is_finite() || residual.abs() > config.water_balance_tolerance_kg {
            return Err(MoistureError::MassBalance {
                residual_kg: residual,
                tolerance_kg: config.water_balance_tolerance_kg,
            });
        }
        let [left, right] = faces;
        let report = MoistureStepReport {
            interval_s: [self.time_s, end],
            boundary_water_kg: transfer,
            surface_moisture: actual_surface,
            reservoir_moisture: surface,
            surface_sorption: [left.sorption, right.sorption],
            surface_transfer: [left.transfer, right.transfer],
            stored_water_kg: [water_before, water_after],
            water_balance_residual_kg: residual,
            transport_start: start,
            transport_end: endpoint,
            linear: solution.report,
            nonlinear,
        };
        let mut next = self.clone();
        next.time_s = end;
        next.moisture = solution.x;
        poll(cx)?;
        Ok((next, report))
    }

    /// One common finite-volume operator for both coefficient time levels.
    fn transport_system(
        &self,
        dt_s: f64,
        faces: &[ResolvedMoistureFace; 2],
        states: &[ResolvedMaterialStatePoint],
        cx: &Cx<'_>,
    ) -> Result<(DiffusionStep, Vec<f64>, [(f64, f64); 2]), MoistureError> {
        poll(cx)?;
        let surface = [faces[0].equilibrium, faces[1].equilibrium];
        let d: Vec<_> = states
            .iter()
            .map(|s| {
                s.property(&self.diffusivity_property)
                    .expect("pinned diffusivity")
                    .value_si()
            })
            .collect();
        let n = d.len();
        let dx = self.geometry.length_m / n as f64;
        let scale = dt_s / dx / dx;
        if !positive(scale) {
            return Err(MoistureError::Invalid(
                "diffusion interval/length scale is unrepresentable",
            ));
        }
        let mut operator = DiffusionStep {
            diagonal: vec![1.0; n],
            off: Vec::with_capacity(n - 1),
        };
        for pair in d.windows(2) {
            // Stable harmonic mean avoids D_i*D_j overflow.
            let lo = pair[0].min(pair[1]);
            let hi = pair[0].max(pair[1]);
            operator.off.push(scale * (2.0 * lo / (1.0 + lo / hi)));
        }
        for (i, &a) in operator.off.iter().enumerate() {
            operator.diagonal[i] += a;
            operator.diagonal[i + 1] += a;
        }
        let coupling = [
            faces[0].coupling(dt_s, dx, d[0], scale)?,
            faces[1].coupling(dt_s, dx, d[n - 1], scale)?,
        ];
        let boundary = [coupling[0].0, coupling[1].0];
        let mut rhs = self.moisture.clone();
        for (side, i) in [0, n - 1].into_iter().enumerate() {
            if let Some(m) = surface[side] {
                operator.diagonal[i] += boundary[side];
                rhs[i] += boundary[side] * m;
            }
        }
        if operator.diagonal.iter().any(|v| !positive(*v))
            || operator.off.iter().any(|v| !positive(*v))
            || rhs.iter().any(|v| !v.is_finite())
            || boundary.iter().any(|v| !v.is_finite())
            || (0..2).any(|side| {
                surface[side].is_some()
                    && faces[side].coefficient_m_s != Some(0.0)
                    && !positive(boundary[side])
            })
        {
            return Err(MoistureError::Invalid(
                "diffusion coefficients overflow or underflow",
            ));
        }
        Ok((operator, rhs, coupling))
    }

    fn boundary_transfer(
        &self,
        moisture: &[f64],
        faces: &[ResolvedMoistureFace; 2],
        coupling: &[(f64, f64); 2],
    ) -> [f64; 2] {
        let n = moisture.len();
        [
            self.dry_mass_per_cell_kg
                * coupling[0].0
                * (faces[0].equilibrium.unwrap_or(moisture[0]) - moisture[0]),
            self.dry_mass_per_cell_kg
                * coupling[1].0
                * (faces[1].equilibrium.unwrap_or(moisture[n - 1]) - moisture[n - 1]),
        ]
    }
}

struct DiffusionStep {
    diagonal: Vec<f64>,
    off: Vec<f64>,
}
impl LinearOp for DiffusionStep {
    fn n(&self) -> usize {
        self.diagonal.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for ((y, &a), &x) in y.iter_mut().zip(&self.diagonal).zip(x) {
            *y = a * x;
        }
        for (i, &a) in self.off.iter().enumerate() {
            y[i] -= a * x[i + 1];
            y[i + 1] -= a * x[i];
        }
    }
}
impl Precond for DiffusionStep {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((z, &r), &a) in z.iter_mut().zip(r).zip(&self.diagonal) {
            *z = r / a;
        }
    }
}
