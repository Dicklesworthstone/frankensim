//! Prescribed interference heating on the original stationary P1 enthalpy owner.
//!
//! The authored nodal quantity is a P1 volumetric thermal-power density, W/m3,
//! not a cell power or an assumed nodal control-volume power. Analytic interval
//! integrals enter the ORIGINAL nonlinear residual as dual nodal joules.
//! Reference mass, caloric charts, conductivity, contact and radiation retain
//! their original owners. No thermal evolution or constitutive law is replaced.
//!
//! Source-ripple admission bounds the positive sum of spatially integrated
//! local primitive deviations; cancelling sources cannot hide local forcing.
//! The independent maximum thermal step controls the first-order endpoint
//! transport approximation. Neither gate is a nonlinear time-error certificate.
//! A prescribed expected power drives a deterministic thermal response; it
//! does not establish the ensemble-mean temperature of a nonlinear system.

use std::collections::BTreeSet;
use std::fmt;

use fs_exec::Cx;
use fs_fft::harmonic_window::{PreparedHarmonicProductSpectrum, ProductSpectrumError};
use fs_math::dd::Dd;
use fs_solver::NewtonReport;

use super::{
    ASSEMBLY_TILE, DofMap, EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig,
    PointwiseBoundary, StepContext, energy_balance, poll,
};
use crate::{
    ConductionProblem, InterfaceFlux, MaterialPointwiseAmbientRadiation,
    MaterialPointwiseEnclosureRadiation, PointwiseAmbientRadiation, PointwiseEnclosureRadiation,
    PointwiseEnclosureReport, PointwiseRadiationReport, RobinFlux, ThermalInterfaces,
};

/// One explicit P1 vertex value. Unbound vertices have exactly zero source.
/// Each name and vertex must be unique within one prepared owner.
#[derive(Clone, Copy, Debug)]
pub struct NodalPowerSpectrum<'a> {
    pub name: &'a str,
    pub vertex: usize,
    pub spectrum: &'a PreparedHarmonicProductSpectrum,
    /// Caller-declared units of the supplied observable.
    pub input_unit: &'a str,
    /// Positive conversion to W/m3, applied BEFORE analytic time integration.
    pub watts_per_m3_per_input_unit: f64,
}

/// Additional source work/storage admission; native spatial and Newton budgets
/// remain those of the original storage owner and EnthalpyStepConfig.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInterferenceBudget {
    pub maximum_bindings: usize,
    /// Sum of retained DC and tone coefficients across all bindings.
    /// Each interval makes two bounded spectral passes, plus P1 mesh traversal.
    pub maximum_spectral_terms: usize,
    /// Bounds this adapter's source, quadrature and audit buffers; excludes
    /// the original nonlinear solver's independently admitted workspace.
    pub maximum_source_workspace_bytes: usize,
}

/// Independent temporal controls for one requested physical interval.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInterferenceControls {
    pub maximum_thermal_step_s: f64,
    /// Joules: positive spatial L1 bound on deviations of local source
    /// primitives from their interval chords. This is not a bound on T or h.
    pub maximum_source_ripple_j: f64,
}

/// Optional original radiation owner. Geometry, boundary and original material
/// chart identities are checked exactly as in the native enthalpy entry points.
#[derive(Clone, Copy)]
pub enum SpatialInterferenceRadiation<'a> {
    Ambient(&'a PointwiseAmbientRadiation<'a>),
    MaterialAmbient(&'a MaterialPointwiseAmbientRadiation<'a>),
    Enclosure(&'a PointwiseEnclosureRadiation<'a>),
    MaterialEnclosure(&'a MaterialPointwiseEnclosureRadiation<'a>),
}

#[derive(Debug)]
pub enum SpatialInterferenceError {
    Enthalpy(EnthalpyError),
    Spectrum {
        vertex: usize,
        source: ProductSpectrumError,
    },
    SourceRipple {
        bound_j: f64,
        tolerance_j: f64,
    },
}

impl fmt::Display for SpatialInterferenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Enthalpy(error) => error.fmt(f),
            Self::Spectrum { vertex, source } => {
                write!(f, "spatial heat spectrum at vertex {vertex}: {source}")
            }
            Self::SourceRipple {
                bound_j,
                tolerance_j,
            } => write!(
                f,
                "spatial input primitive ripple {bound_j} J exceeds {tolerance_j} J; shorten the interval"
            ),
        }
    }
}

impl std::error::Error for SpatialInterferenceError {}

impl From<EnthalpyError> for SpatialInterferenceError {
    fn from(value: EnthalpyError) -> Self {
        Self::Enthalpy(value)
    }
}

/// Signed physical accounts of THIS accepted interval. Source and storage are
/// retained in DD; boundary powers use the original endpoint quadratures.
/// Internal contact exchanges are not added to the external accounts.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInterferenceEnergy {
    pub source_j: Dd,
    pub stored_energy_change_j: Dd,
    pub neumann_out_j: Dd,
    pub robin_out_j: Dd,
    pub radiation_out_j: Dd,
    /// Actual constrained-node reaction, inward, including spectral source
    /// and any first-step fixed-temperature jump. Never inferred from mean W.
    pub dirichlet_in_j: Dd,
    pub residual_j: Dd,
}

/// An owned endpoint, with no borrowed source and no fictitious mean-power
/// report. The returned h field and DD end clock can start the next interval.
#[derive(Debug)]
pub struct SpatialInterferenceStep {
    pub start_clock_s: Dd,
    pub end_clock_s: Dd,
    pub duration_s: f64,
    pub source_ripple_bound_j: f64,
    pub specific_enthalpy_j_kg: Vec<f64>,
    pub temperature: Vec<f64>,
    pub liquid_mass_fraction: Vec<f64>,
    /// Assembled dual load, not a nodal density, in canonical vertex order.
    pub source_load_j: Vec<Dd>,
    /// Independent physical residual on free nodes; constrained rows are zero.
    /// Their complete physical balance is in dirichlet_in_j instead.
    pub nodal_residual_j: Vec<Dd>,
    pub energy: SpatialInterferenceEnergy,
    pub robin_fluxes: Vec<RobinFlux>,
    pub contact_fluxes: Vec<InterfaceFlux>,
    pub pointwise_radiation: Option<PointwiseRadiationReport>,
    pub pointwise_enclosure_radiation: Option<PointwiseEnclosureReport>,
    pub newton: NewtonReport,
}

/// Original storage, mesh and ordered source bindings. Labels never establish
/// geometric identity: queries must use the exact mesh retained by storage.
/// Sources are prescribed and immutable throughout an interval. Temperature
/// feedback requires rebuilding them through the original field/material owner.
pub struct PreparedSpatialInterference<'a, 'm, 'c> {
    storage: &'a EnthalpyBackwardEuler<'m, 'c>,
    bindings: Vec<NodalPowerSpectrum<'a>>,
    source_clock: &'a str,
    spectral_terms: usize,
}

impl<'a, 'm, 'c> PreparedSpatialInterference<'a, 'm, 'c> {
    pub fn new(
        cx: &Cx<'_>,
        storage: &'a EnthalpyBackwardEuler<'m, 'c>,
        bindings: &[NodalPowerSpectrum<'a>],
        source_clock: &'a str,
        budget: SpatialInterferenceBudget,
    ) -> Result<Self, SpatialInterferenceError> {
        poll(cx, 0)?;
        if source_clock.trim().is_empty() {
            return Err(
                EnthalpyError::InvalidInput("a common source clock must be declared").into(),
            );
        }
        admit_budget("spectral bindings", bindings.len(), budget.maximum_bindings)?;
        let n = storage.masses.len();
        // Density and dual-load vectors, DD audit scratch, scalar bounds and
        // duplicate bookkeeping. Native fields/assembly remain separately owned.
        let bytes = n
            .checked_mul(
                4 * std::mem::size_of::<Dd>()
                    + std::mem::size_of::<f64>()
                    + std::mem::size_of::<bool>(),
            )
            .and_then(|v| {
                bindings
                    .len()
                    .checked_mul(
                        std::mem::size_of::<NodalPowerSpectrum<'_>>()
                            + 8 * std::mem::size_of::<usize>(),
                    )
                    .and_then(|b| v.checked_add(b))
            })
            .filter(|&v| isize::try_from(v).is_ok())
            .ok_or(EnthalpyError::InvalidInput(
                "source workspace size overflow",
            ))?;
        admit_budget(
            "spectral source workspace bytes",
            bytes,
            budget.maximum_source_workspace_bytes,
        )?;
        let mut vertices = vec![false; n];
        let mut names = BTreeSet::new();
        let mut terms = 0usize;
        for (index, binding) in bindings.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            if binding.vertex >= n
                || binding.name.trim().is_empty()
                || binding.input_unit.trim().is_empty()
                || !binding.watts_per_m3_per_input_unit.is_finite()
                || binding.watts_per_m3_per_input_unit <= 0.0
            {
                return Err(EnthalpyError::InvalidInput(
                    "each source needs an original vertex, name, unit and positive W/m3 conversion",
                )
                .into());
            }
            if vertices[binding.vertex] || !names.insert(binding.name) {
                return Err(EnthalpyError::InvalidInput(
                    "duplicate spectral vertex or source name",
                )
                .into());
            }
            vertices[binding.vertex] = true;
            terms = terms
                .checked_add(binding.spectrum.tones().len())
                .and_then(|v| v.checked_add(1))
                .ok_or(EnthalpyError::InvalidInput("spectral work overflow"))?;
        }
        admit_budget("spectral terms", terms, budget.maximum_spectral_terms)?;
        terms
            .checked_mul(2)
            .ok_or(EnthalpyError::InvalidInput("spectral work overflow"))?;
        let mut bindings = bindings.to_vec();
        bindings.sort_by_key(|binding| binding.vertex);
        poll(cx, bindings.len())?;
        Ok(Self {
            storage,
            bindings,
            source_clock,
            spectral_terms: terms,
        })
    }

    pub fn bindings(&self) -> &[NodalPowerSpectrum<'a>] {
        &self.bindings
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub const fn spectral_terms(&self) -> usize {
        self.spectral_terms
    }

    /// Advance exactly this interval using the existing nonlinear residual,
    /// tangent, globalization and Cx cancellation. The caller supplies the
    /// original mesh/material/boundary/interface owners and immutable history.
    ///
    /// The static volumetric source must be zero: every prescribed volumetric
    /// contribution belongs in these explicit spectra, including DC. Native
    /// Neumann, Robin, fixed-temperature reactions and contact remain active.
    /// Fixed/declared contact resistance is supported; temperature-contact
    /// feedback owners require their own integrated-source composition.
    ///
    /// Both free-node and total physical joule gates use the unchanged
    /// config.energy_tolerance_j. Source ripple and maximum thermal step are
    /// separate admission checks. A refusal returns no endpoint or source debit.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_interval(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        start_clock_s: Dd,
        duration_s: f64,
        config: EnthalpyStepConfig,
        controls: SpatialInterferenceControls,
        radiation: Option<SpatialInterferenceRadiation<'_>>,
    ) -> Result<SpatialInterferenceStep, SpatialInterferenceError> {
        poll(cx, 0)?;
        if !controls.maximum_thermal_step_s.is_finite()
            || controls.maximum_thermal_step_s <= 0.0
            || !controls.maximum_source_ripple_j.is_finite()
            || controls.maximum_source_ripple_j < 0.0
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || duration_s > controls.maximum_thermal_step_s
        {
            return Err(EnthalpyError::InvalidInput(
                "positive interval within the thermal step cap and finite nonnegative ripple allowance required"
            ).into());
        }
        let start = checked_dd(Dd::from_pair(start_clock_s.hi, start_clock_s.lo))?;
        let end = checked_dd(start + Dd::from_f64(duration_s))?;
        if !positive(end - start) {
            return Err(EnthalpyError::InvalidInput(
                "interval does not advance the DD source clock",
            )
            .into());
        }
        self.storage
            .admit_step(cx, problem, old_specific_h, duration_s, config)?;
        problem
            .source
            .validate("spectral baseline source", self.storage.masses.len())
            .map_err(EnthalpyError::from)?;
        for vertex in 0..self.storage.masses.len() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if problem.source.at(vertex) != 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "static volumetric source must be zero; bind all volumetric input exactly once",
                )
                .into());
            }
        }
        let pointwise_radiation = self.admit_radiation(cx, problem, radiation)?;
        let n = self.storage.masses.len();
        let mut density_j_m3 = vec![Dd::ZERO; n];
        let mut ripple_j_m3 = vec![0.0; n];
        for (index, binding) in self.bindings.iter().enumerate() {
            poll(cx, index)?;
            let source_error = |source| SpatialInterferenceError::Spectrum {
                vertex: binding.vertex,
                source,
            };
            let raw_bound = binding
                .spectrum
                .input_ripple_bound_cancellable(duration_s, || cx.checkpoint().is_err())
                .map_err(source_error)?;
            ripple_j_m3[binding.vertex] =
                outward_product(raw_bound, binding.watts_per_m3_per_input_unit);
            density_j_m3[binding.vertex] = binding
                .spectrum
                .causal_integral_from_clock_dd_scaled_cancellable(
                    start,
                    duration_s,
                    0.0,
                    binding.watts_per_m3_per_input_unit,
                    || cx.checkpoint().is_err(),
                )
                .map_err(source_error)?
                .integral_dd();
        }
        let source_ripple_bound_j = spatial_ripple_bound(cx, self.storage.mesh, &ripple_j_m3)?;
        if source_ripple_bound_j > controls.maximum_source_ripple_j {
            return Err(SpatialInterferenceError::SourceRipple {
                bound_j: source_ripple_bound_j,
                tolerance_j: controls.maximum_source_ripple_j,
            });
        }
        let (source_load_j, source_j) = integrate_p1_source(cx, self.storage.mesh, &density_j_m3)?;
        let context = StepContext {
            storage: self.storage,
            cx,
            problem,
            interfaces,
            old: old_specific_h,
            dt: duration_s,
            pointwise_radiation,
        };
        let (h, newton) = context.solve_newton(config, Some(&source_load_j))?;
        finish_interval(
            &context,
            h,
            newton,
            source_load_j,
            source_j,
            start,
            end,
            source_ripple_bound_j,
            config.energy_tolerance_j,
        )
        .map_err(SpatialInterferenceError::from)
    }

    fn admit_radiation<'r>(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        radiation: Option<SpatialInterferenceRadiation<'r>>,
    ) -> Result<Option<PointwiseBoundary<'r>>, EnthalpyError> {
        let Some(radiation) = radiation else {
            return Ok(None);
        };
        let (matches, native) = match radiation {
            SpatialInterferenceRadiation::Ambient(owner) => (
                owner.matches(problem.mesh, problem.boundary),
                PointwiseBoundary::Ambient(owner),
            ),
            SpatialInterferenceRadiation::Enclosure(owner) => (
                owner.matches(problem.mesh, problem.boundary),
                PointwiseBoundary::Enclosure(owner),
            ),
            SpatialInterferenceRadiation::MaterialAmbient(owner) => (
                owner.matches(problem.mesh, problem.boundary)
                    && owner.matches_curves(cx, |vertex| self.storage.curve_for_vertex(vertex))?,
                PointwiseBoundary::MaterialAmbient(owner),
            ),
            SpatialInterferenceRadiation::MaterialEnclosure(owner) => (
                owner.matches(problem.mesh, problem.boundary)
                    && owner.matches_curves(cx, |vertex| self.storage.curve_for_vertex(vertex))?,
                PointwiseBoundary::MaterialEnclosure(owner),
            ),
        };
        if !matches {
            return Err(EnthalpyError::InvalidInput(
                "radiation must retain the exact spatial mesh, boundary and original caloric bindings",
            ));
        }
        Ok(Some(native))
    }
}

fn admit_budget(
    resource: &'static str,
    required: usize,
    limit: usize,
) -> Result<(), EnthalpyError> {
    if required > limit {
        Err(EnthalpyError::Budget {
            resource,
            required,
            limit,
        })
    } else {
        Ok(())
    }
}

fn checked_dd(value: Dd) -> Result<Dd, EnthalpyError> {
    if value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite() {
        Ok(value)
    } else {
        Err(EnthalpyError::NonFiniteArithmetic)
    }
}

fn is_zero(value: Dd) -> bool {
    value.hi == 0.0 && value.lo == 0.0
}
fn positive(value: Dd) -> bool {
    value.hi > 0.0 || (value.hi == 0.0 && value.lo > 0.0)
}
fn absolute(value: Dd) -> Dd {
    if positive(-value) { -value } else { value }
}

fn product(left: Dd, right: Dd) -> Result<Dd, EnthalpyError> {
    let value = checked_dd(left * right)?;
    if !is_zero(left) && !is_zero(right) && is_zero(value) {
        Err(EnthalpyError::NonFiniteArithmetic)
    } else {
        Ok(value)
    }
}

fn divided(value: Dd, divisor: f64) -> Result<Dd, EnthalpyError> {
    let result = checked_dd(value / Dd::from_f64(divisor))?;
    if !is_zero(value) && is_zero(result) {
        Err(EnthalpyError::NonFiniteArithmetic)
    } else {
        Ok(result)
    }
}

fn residual_dd(
    mass: f64,
    new_h: f64,
    old_h: f64,
    dt: f64,
    transport_w: f64,
    boundary_load_w: f64,
    source_j: Dd,
) -> Result<Dd, EnthalpyError> {
    let storage = product(
        Dd::from_f64(mass),
        checked_dd(Dd::from_f64(new_h) - Dd::from_f64(old_h))?,
    )?;
    let transport = product(
        Dd::from_f64(dt),
        checked_dd(Dd::from_f64(transport_w) - Dd::from_f64(boundary_load_w))?,
    )?;
    checked_dd(storage + transport - source_j)
}

// Original Newton receives the physical source subtraction in joules, before
// its final f64 residual interface. Signed cancellation is legitimate.
pub(super) fn physical_residual_j(
    mass: f64,
    new_h: f64,
    old_h: f64,
    dt: f64,
    transport_w: f64,
    boundary_load_w: f64,
    source_j: Dd,
) -> Result<f64, EnthalpyError> {
    let value = residual_dd(
        mass,
        new_h,
        old_h,
        dt,
        transport_w,
        boundary_load_w,
        source_j,
    )?;
    let rounded = value.to_f64();
    if !is_zero(value) && rounded == 0.0 {
        Err(EnthalpyError::NonFiniteArithmetic)
    } else {
        Ok(rounded)
    }
}

// Exact P1 source rule on the authoritative mesh volumes. The physical volume
// integral is computed independently of the assembled dual row sum.
fn integrate_p1_source(
    cx: &Cx<'_>,
    mesh: &crate::ConductionMesh,
    density: &[Dd],
) -> Result<(Vec<Dd>, Dd), EnthalpyError> {
    let mut load = vec![Dd::ZERO; mesh.vertex_count()];
    let mut total = Dd::ZERO;
    for (element, vertices) in mesh.complex().tets.iter().enumerate() {
        if element % ASSEMBLY_TILE == 0 {
            poll(cx, element)?;
        }
        let volume = Dd::from_f64(mesh.element_volume(element));
        let mut sum = Dd::ZERO;
        for &vertex in vertices {
            sum = checked_dd(sum + density[vertex as usize])?;
        }
        total = checked_dd(total + divided(product(volume, sum)?, 4.0)?)?;
        for &vertex in vertices {
            let vertex = vertex as usize;
            let row = divided(product(volume, checked_dd(sum + density[vertex])?)?, 20.0)?;
            load[vertex] = checked_dd(load[vertex] + row)?;
        }
    }
    poll(cx, mesh.element_count())?;
    Ok((load, total))
}

// Outward f64 arithmetic is only for nonnegative source-ripple bounds.
// Preserve exact darkness, and round positive underflow UP, never to darkness.
fn outward_product(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        0.0
    } else {
        (a * b).next_up()
    }
}

fn outward_sum(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        (a + b).next_up()
    }
}

fn spatial_ripple_bound(
    cx: &Cx<'_>,
    mesh: &crate::ConductionMesh,
    nodal: &[f64],
) -> Result<f64, EnthalpyError> {
    let mut bound = 0.0;
    for (element, vertices) in mesh.complex().tets.iter().enumerate() {
        if element % ASSEMBLY_TILE == 0 {
            poll(cx, element)?;
        }
        let weight = outward_product(mesh.element_volume(element), 0.25);
        for &vertex in vertices {
            bound = outward_sum(bound, outward_product(weight, nodal[vertex as usize]));
        }
    }
    poll(cx, mesh.element_count())?;
    Ok(bound)
}

#[allow(clippy::too_many_arguments)]
fn finish_interval(
    context: &StepContext<'_, '_, '_, '_>,
    h: Vec<f64>,
    newton: NewtonReport,
    source_load_j: Vec<Dd>,
    source_j: Dd,
    start_clock_s: Dd,
    end_clock_s: Dd,
    source_ripple_bound_j: f64,
    tolerance_j: f64,
) -> Result<SpatialInterferenceStep, EnthalpyError> {
    context.check_prescribed_values(&h)?;
    let temperature = context.temperatures(&h)?;
    let scales = context.phase_scales(&h)?;
    let system = context.assemble(&temperature, scales.as_deref())?;
    let radiation = context
        .pointwise_radiation
        .map(|owner| owner.load_and_report(context.cx, &temperature, &h))
        .transpose()?;
    let mut storage = Dd::ZERO;
    let mut dirichlet = Dd::ZERO;
    let mut nodal_residual_j = Vec::with_capacity(h.len());
    let mut liquid_mass_fraction = Vec::with_capacity(h.len());
    let tolerance = Dd::from_f64(tolerance_j);
    for (vertex, &value) in h.iter().enumerate() {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(context.cx, vertex)?;
        }
        storage = checked_dd(
            storage
                + product(
                    Dd::from_f64(context.storage.masses[vertex]),
                    checked_dd(Dd::from_f64(value) - Dd::from_f64(context.old[vertex]))?,
                )?,
        )?;
        liquid_mass_fraction.push(
            context
                .storage
                .curve_for_vertex(vertex)
                .values_at_specific_enthalpy(value)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?
                .liquid_mass_fraction.clamp(0.0, 1.0),
        );
        let (columns, entries) = system.operator.row(vertex);
        let transport: f64 = columns
            .iter()
            .zip(entries)
            .map(|(&column, &entry)| entry * temperature[column])
            .sum();
        let transport = match &radiation {
            Some((load, _)) => super::finite(transport + load[vertex])?,
            None => transport,
        };
        let residual = residual_dd(
            context.storage.masses[vertex],
            value,
            context.old[vertex],
            context.dt,
            transport,
            system.load[vertex],
            source_load_j[vertex],
        )?;
        if context.is_prescribed(vertex) {
            dirichlet = checked_dd(dirichlet + residual)?;
            nodal_residual_j.push(Dd::ZERO);
        } else {
            if positive(absolute(residual) - tolerance) {
                return Err(EnthalpyError::EndpointResidual {
                    residual_j: absolute(residual).to_f64(),
                    tolerance_j,
                });
            }
            nodal_residual_j.push(residual);
        }
    }
    let dofs = DofMap::new_allow_fully_prescribed(context.problem.boundary, h.len())?;
    let (boundary, robin_fluxes) = energy_balance(
        context.storage.mesh,
        context.problem.boundary,
        context.problem.source,
        &system,
        &dofs,
        &temperature,
    );
    let dt = Dd::from_f64(context.dt);
    let neumann = product(dt, Dd::from_f64(boundary.neumann_out_w))?;
    let robin = product(dt, Dd::from_f64(boundary.robin_out_w))?;
    let radiation_out = match &radiation {
        Some((_, report)) => product(dt, Dd::from_f64(report.outward_heat_w()))?,
        None => Dd::ZERO,
    };
    let residual =
        checked_dd(storage - checked_dd(source_j + dirichlet - neumann - robin - radiation_out)?)?;
    if positive(absolute(residual) - tolerance) {
        return Err(EnthalpyError::EnergyBalance {
            residual_j: residual.to_f64(),
            tolerance_j,
        });
    }
    let contact_fluxes = match context.interfaces {
        Some(interfaces) => interfaces.fluxes(&temperature)?,
        None => Vec::new(),
    };
    let (pointwise_radiation, pointwise_enclosure_radiation) =
        radiation.map_or((None, None), |(_, report)| report.into_reports());
    poll(context.cx, h.len())?;
    Ok(SpatialInterferenceStep {
        start_clock_s,
        end_clock_s,
        duration_s: context.dt,
        source_ripple_bound_j,
        specific_enthalpy_j_kg: h,
        temperature,
        liquid_mass_fraction,
        source_load_j,
        nodal_residual_j,
        energy: SpatialInterferenceEnergy {
            source_j,
            stored_energy_change_j: storage,
            neumann_out_j: neumann,
            robin_out_j: robin,
            radiation_out_j: radiation_out,
            dirichlet_in_j: dirichlet,
            residual_j: residual,
        },
        robin_fluxes,
        contact_fluxes,
        pointwise_radiation,
        pointwise_enclosure_radiation,
        newton,
    })
}
