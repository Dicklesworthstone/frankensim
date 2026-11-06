//! Conserved-inventory equilibrium of the original spatial enthalpy owner.
//!
//! Every connected component of the admitted heat operator keeps its own
//! reference enthalpy. Its equilibrium temperature follows from original bulk
//! and bonded-film storage. A sharp latent interval remains an interval; a
//! scalar energy constraint does not determine the spatial phase distribution.
//! This computes stationary end states, not relaxation times or a replacement
//! for the original transient transport solve.

use std::{collections::BTreeSet, fmt};

use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_solver::norm2;
use fs_sparse::Csr;

use super::{EnthalpyBackwardEuler, EnthalpyError, StepContext, finite, poll};
use crate::{
    ConductionError, ConductionProblem, ThermalBc, ThermalInterfaces,
    assemble::{ASSEMBLY_TILE, full_residual},
    material::TemperatureSpan,
};

/// Explicit work and physical acceptance limits for isolated equilibrium.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsolatedEquilibriumConfig {
    /// Maximum caloric evaluations per connected component, including its two
    /// support endpoints and exact latent-temperature probes.
    pub max_temperature_evaluations: usize,
    /// Absolute conserved reference-enthalpy defect per component, J.
    pub energy_tolerance_j: f64,
    /// Allowed chart-resolved temperature mismatch from a component's common
    /// temperature, K. Does not broaden any material's support.
    pub temperature_tolerance_k: f64,
    /// Norm of ALL original physical heat rows, W, at a materialized endpoint.
    pub heat_residual_tolerance_w: f64,
    /// Allowed negative entropy roundoff per isolated component, J/K.
    pub entropy_tolerance_j_k: f64,
}

/// No equilibrium or partial physical state is published on refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum IsolatedEquilibriumError {
    /// Original mesh, material, chart, film, arithmetic or cancellation refusal.
    Original(EnthalpyError),
    /// A connected inventory has no common supported temperature.
    NoCommonTemperature { component: usize },
    /// Conserved energy lies outside the common caloric support.
    EnergyOutsideDomain {
        component: usize,
        energy_change_interval_j: [f64; 2],
    },
    /// The bounded scalar solve could not close the original energy.
    WorkLimit {
        component: usize,
        evaluations: usize,
    },
    /// No further temperature between these floating-point endpoints exists.
    TemperatureResolution {
        component: usize,
        bracket_k: [f64; 2],
    },
    /// A latent phase distribution is not determined by total energy alone.
    AllocationRequired,
    /// A proposed enthalpy is outside its original isothermal interval.
    AllocationOutsideInterval {
        vertex: usize,
        specific_enthalpy_j_kg: f64,
    },
    /// A resolved nodal temperature fails the common-temperature tolerance.
    TemperatureMismatch { vertex: usize, difference_k: f64 },
    /// One component failed its own original energy balance.
    EnergyBalance {
        component: usize,
        residual_j: f64,
        tolerance_j: f64,
    },
    /// An isolated component lost entropy beyond the declared roundoff allowance.
    EntropyDecrease {
        component: usize,
        change_j_k: f64,
        tolerance_j_k: f64,
    },
    /// The original physical heat equations are not stationary.
    HeatResidual { residual_w: f64, tolerance_w: f64 },
}

impl From<EnthalpyError> for IsolatedEquilibriumError {
    fn from(error: EnthalpyError) -> Self {
        Self::Original(error)
    }
}
impl From<ConductionError> for IsolatedEquilibriumError {
    fn from(error: ConductionError) -> Self {
        Self::Original(error.into())
    }
}
impl fmt::Display for IsolatedEquilibriumError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Original(error) => error.fmt(f),
            other => write!(f, "isolated thermal equilibrium: {other:?}"),
        }
    }
}
impl std::error::Error for IsolatedEquilibriumError {}

/// The caloric constraint for one thermally connected, isolated inventory.
#[derive(Debug, Clone, PartialEq)]
pub struct IsolatedEquilibriumComponent {
    /// Common absolute temperature, K.
    pub temperature_k: f64,
    /// Number of original reference-mass vertices in this component.
    pub vertices: usize,
    /// Original energy change at the all-lower and all-upper enthalpy states,
    /// J. A plateau allocation must close energy between these extremes.
    pub energy_change_interval_j: [f64; 2],
    /// Energy to distribute above the all-lower state, J; the constraint is
    /// sum(m_i * (h_i - lower_i)) = this value within energy_tolerance_j.
    /// Sensible films are already included in the all-lower energy.
    pub latent_energy_above_lower_j: f64,
    /// Actual bounded caloric evaluations used by this component.
    pub temperature_evaluations: usize,
}

/// A conserved-energy equilibrium family bound to the exact original storage,
/// mesh, source, boundary and initial history.
///
/// Every h_i must lie in specific_enthalpy_intervals_j_kg(), and each connected
/// component retains one independent extensive energy constraint. On an
/// isothermal latent plateau there can be many valid phase distributions.
/// materialize() requires the caller's explicit selection and checks the
/// original physics; unique_state() is available when every interval is a point.
#[derive(Debug)]
pub struct IsolatedEnthalpyEquilibrium<'s, 'm, 'c> {
    storage: &'s EnthalpyBackwardEuler<'m, 'c>,
    problem: ConductionProblem<'s>,
    interfaces: Option<&'s ThermalInterfaces>,
    old_h: &'s [f64],
    config: IsolatedEquilibriumConfig,
    component_of_vertex: Vec<usize>,
    components: Vec<IsolatedEquilibriumComponent>,
    intervals: Vec<[f64; 2]>,
}

/// An explicitly resolved stationary enthalpy state, checked on the native
/// transport operator. This can be passed to the existing transient owner.
#[derive(Debug, Clone, PartialEq)]
pub struct IsolatedEquilibriumState {
    /// Original-chart specific enthalpies in mesh order, J/kg.
    pub specific_enthalpy_j_kg: Vec<f64>,
    /// Actual chart-resolved temperatures in mesh order, K.
    pub temperature: Vec<f64>,
    /// Actual original-chart liquid mass fractions.
    pub liquid_mass_fraction: Vec<f64>,
    /// Independently evaluated original-to-final storage change per component, J.
    pub component_energy_residual_j: Vec<f64>,
    /// Original-to-final bulk and film entropy change per component, J/K.
    pub component_entropy_change_j_k: Vec<f64>,
    /// Norm of every original physical heat row, W.
    pub heat_residual_norm_w: f64,
}

impl<'m, 'c> EnthalpyBackwardEuler<'m, 'c> {
    /// Find the stationary temperature of each isolated connected inventory.
    ///
    /// The same reference nodal masses, caloric charts, sensible films and
    /// material transport owners used by transient evolution are retained.
    /// Components come from the ACTUAL admitted native transport matrix:
    /// disconnected bodies conserve energy separately, while finite contacts
    /// join their inventories. No artificial reservoir or large timestep is used.
    ///
    /// All volume sources and all boundary fluxes must be exactly zero;
    /// prescribed temperatures and Robin reservoirs refuse. Balanced but
    /// nonzero external heating/cooling is not an isolated system. Radiation,
    /// moving geometry, pressure work and changing internal variables are not
    /// part of this endpoint. The conserved quantity is the original model's
    /// reference-mass enthalpy plus film storage, not reinterpreted internal energy.
    ///
    /// Constitutive monotonicity makes the component caloric energy monotone.
    /// Exact latent temperatures are searched before bounded safeguarded
    /// interpolation in sensible intervals. Original source/transport support
    /// is checked at the initial state and again at a materialized endpoint.
    /// Caloric families alone are not accepted transport states.
    ///
    /// # Errors
    /// Original admission/cancellation, open boundaries or sources, disjoint
    /// caloric support, unresolvable energy, or exhausted scalar work.
    #[allow(clippy::too_many_lines)]
    pub fn isolated_equilibrium<'s>(
        &'s self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'s>,
        interfaces: Option<&'s ThermalInterfaces>,
        old_h: &'s [f64],
        config: IsolatedEquilibriumConfig,
    ) -> Result<IsolatedEnthalpyEquilibrium<'s, 'm, 'c>, IsolatedEquilibriumError> {
        poll(cx, 0)?;
        validate_config(config)?;
        let n = self.masses.len();
        if !std::ptr::eq(self.mesh, problem.mesh) || old_h.len() != n {
            return Err(invalid(
                "storage and initial history must retain the exact transport mesh",
            ));
        }
        admit_isolated(cx, problem, n)?;
        if self
            .sourced_melting
            .as_ref()
            .is_some_and(|s| !s.matches_numerical_basis(problem))
            || self
                .sourced_phase
                .as_ref()
                .is_some_and(|s| !s.matches_numerical_basis(problem))
        {
            return Err(invalid(
                "source-bound conductivity requires its original numerical basis",
            ));
        }
        if let Some(source) = &self.sourced_phase {
            source.admit_enthalpies(cx, old_h)?;
        }
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_h,
            dt: 1.0,
            pointwise_radiation: None,
        };
        let old_temperature = context.temperatures(old_h)?;
        let conductivity = context.phase_scales(old_h)?;
        let system = context.assemble(&old_temperature, conductivity.as_ref())?;
        let (component_of_vertex, groups) = components(cx, &system.operator)?;
        let mut result = IsolatedEnthalpyEquilibrium {
            storage: self,
            problem,
            interfaces,
            old_h,
            config,
            component_of_vertex,
            components: Vec::with_capacity(groups.len()),
            intervals: vec![[0.0; 2]; n],
        };
        for (component, vertices) in groups.iter().enumerate() {
            poll(cx, component)?;
            let report = solve_component(
                cx,
                self,
                old_h,
                vertices,
                component,
                config,
                &mut result.intervals,
            )?;
            result.components.push(report);
        }
        poll(cx, n)?;
        Ok(result)
    }
}

impl IsolatedEnthalpyEquilibrium<'_, '_, '_> {
    /// Component index of each original mesh vertex.
    #[must_use]
    pub fn component_of_vertex(&self) -> &[usize] {
        &self.component_of_vertex
    }

    /// Independent common-temperature and energy constraints.
    #[must_use]
    pub fn components(&self) -> &[IsolatedEquilibriumComponent] {
        &self.components
    }

    /// Inclusive original caloric h intervals, J/kg, in mesh vertex order.
    #[must_use]
    pub fn specific_enthalpy_intervals_j_kg(&self) -> &[[f64; 2]] {
        &self.intervals
    }

    /// Whether a phase allocation remains to be supplied explicitly.
    #[must_use]
    pub fn requires_allocation(&self) -> bool {
        self.intervals.iter().any(|range| range[0] != range[1])
    }

    /// Materialize the unique sensible state without inventing a latent rule.
    ///
    /// # Errors
    /// AllocationRequired on a nontrivial latent family, or any materialize refusal.
    pub fn unique_state(
        &self,
        cx: &Cx<'_>,
    ) -> Result<IsolatedEquilibriumState, IsolatedEquilibriumError> {
        poll(cx, 0)?;
        if self.requires_allocation() {
            return Err(IsolatedEquilibriumError::AllocationRequired);
        }
        let h: Vec<_> = self.intervals.iter().map(|range| range[0]).collect();
        self.materialize(cx, &h)
    }

    /// Admit an explicit latent allocation on the original physical system.
    ///
    /// Both each component's energy and its entropy are checked independently;
    /// equal and opposite defects in disconnected bodies cannot cancel. The
    /// complete original heat operator is reassembled at the actual phase state,
    /// including conductivity changes and contacts. A calorically admissible
    /// phase allocation that leaves a source domain is refused here.
    ///
    /// No field, constitutive history or earlier accepted state is mutated.
    ///
    /// # Errors
    /// Wrong shape, out-of-interval enthalpy, unresolved common temperature,
    /// component energy/entropy failure, native source/heat refusal or cancellation.
    #[allow(clippy::too_many_lines)]
    pub fn materialize(
        &self,
        cx: &Cx<'_>,
        h: &[f64],
    ) -> Result<IsolatedEquilibriumState, IsolatedEquilibriumError> {
        poll(cx, 0)?;
        if h.len() != self.intervals.len() {
            return Err(invalid(
                "one equilibrium enthalpy per original vertex required",
            ));
        }
        for (vertex, (&value, range)) in h.iter().zip(&self.intervals).enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if !value.is_finite() || value < range[0] || value > range[1] {
                return Err(IsolatedEquilibriumError::AllocationOutsideInterval {
                    vertex,
                    specific_enthalpy_j_kg: value,
                });
            }
        }
        if let Some(source) = &self.storage.sourced_phase {
            source.admit_enthalpies(cx, h)?;
        }
        let context = StepContext {
            storage: self.storage,
            cx,
            problem: self.problem,
            interfaces: self.interfaces,
            old: self.old_h,
            dt: 1.0,
            pointwise_radiation: None,
        };
        let temperature = context.temperatures(h)?;
        let mut energy = vec![Dd::from_f64(0.0); self.components.len()];
        let mut entropy = vec![Dd::from_f64(0.0); self.components.len()];
        let mut fraction = Vec::with_capacity(h.len());
        for (vertex, ((&value, &old), &t)) in h.iter().zip(self.old_h).zip(&temperature).enumerate()
        {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            let component = self.component_of_vertex[vertex];
            let difference_k = finite(t - self.components[component].temperature_k)?;
            if difference_k.abs() > self.config.temperature_tolerance_k {
                return Err(IsolatedEquilibriumError::TemperatureMismatch {
                    vertex,
                    difference_k,
                });
            }
            energy[component] =
                energy[component] + self.storage.storage_change_dd(vertex, value, old, t)?;
            entropy[component] = entropy[component]
                + Dd::from_f64(
                    self.storage
                        .nodal_entropy_change_j_k(vertex, old, value, t)?,
                );
            fraction.push(
                self.storage
                    .curve_for_vertex(vertex)
                    .values_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?
                    .liquid_mass_fraction,
            );
        }
        let energy: Vec<_> = energy.into_iter().map(dd_value).collect::<Result<_, _>>()?;
        let entropy: Vec<_> = entropy
            .into_iter()
            .map(dd_value)
            .collect::<Result<_, _>>()?;
        for component in 0..self.components.len() {
            poll(cx, component)?;
            if energy[component].abs() > self.config.energy_tolerance_j {
                return Err(IsolatedEquilibriumError::EnergyBalance {
                    component,
                    residual_j: energy[component],
                    tolerance_j: self.config.energy_tolerance_j,
                });
            }
            if entropy[component] < -self.config.entropy_tolerance_j_k {
                return Err(IsolatedEquilibriumError::EntropyDecrease {
                    component,
                    change_j_k: entropy[component],
                    tolerance_j_k: self.config.entropy_tolerance_j_k,
                });
            }
        }
        let conductivity = context.phase_scales(h)?;
        let system = context.assemble(&temperature, conductivity.as_ref())?;
        let heat_residual_norm_w = finite(norm2(&full_residual(&system, &temperature)))?;
        if heat_residual_norm_w > self.config.heat_residual_tolerance_w {
            return Err(IsolatedEquilibriumError::HeatResidual {
                residual_w: heat_residual_norm_w,
                tolerance_w: self.config.heat_residual_tolerance_w,
            });
        }
        poll(cx, h.len())?;
        Ok(IsolatedEquilibriumState {
            specific_enthalpy_j_kg: h.to_vec(),
            temperature,
            liquid_mass_fraction: fraction,
            component_energy_residual_j: energy,
            component_entropy_change_j_k: entropy,
            heat_residual_norm_w,
        })
    }
}

fn validate_config(config: IsolatedEquilibriumConfig) -> Result<(), IsolatedEquilibriumError> {
    if config.max_temperature_evaluations < 2
        || [
            config.energy_tolerance_j,
            config.temperature_tolerance_k,
            config.heat_residual_tolerance_w,
            config.entropy_tolerance_j_k,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(invalid(
            "at least two caloric evaluations and finite nonnegative tolerances required",
        ));
    }
    Ok(())
}

fn admit_isolated(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    n: usize,
) -> Result<(), IsolatedEquilibriumError> {
    if problem.boundary.face_count() != problem.mesh.boundary().len()
        || !problem.boundary.dirichlet().is_empty()
        || problem.boundary.has_robin()
    {
        return Err(invalid(
            "isolated equilibrium requires the original unanchored adiabatic boundary",
        ));
    }
    problem.source.validate("isolated volume source", n)?;
    for vertex in 0..n {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(cx, vertex)?;
        }
        if problem.source.at(vertex) != 0.0 {
            return Err(invalid(
                "isolated equilibrium requires zero volume source at every vertex",
            ));
        }
    }
    for (slot, face) in problem.mesh.boundary().iter().enumerate() {
        if slot % ASSEMBLY_TILE == 0 {
            poll(cx, slot)?;
        }
        if let Some(condition) = problem.boundary.condition_for(slot) {
            condition.validate(n)?;
            match condition {
                ThermalBc::Neumann { outward_flux } => {
                    if face
                        .vertices
                        .iter()
                        .any(|&v| outward_flux.at(v as usize) != 0.0)
                    {
                        return Err(invalid(
                            "isolated equilibrium requires zero flux on every boundary face",
                        ));
                    }
                }
                _ => return Err(invalid("an external temperature reservoir is not isolated")),
            }
        }
    }
    Ok(())
}

fn components(
    cx: &Cx<'_>,
    operator: &Csr,
) -> Result<(Vec<usize>, Vec<Vec<usize>>), IsolatedEquilibriumError> {
    let n = operator.nrows();
    let mut parent: Vec<_> = (0..n).collect();
    for row in 0..n {
        if row % ASSEMBLY_TILE == 0 {
            poll(cx, row)?;
        }
        let (columns, values) = operator.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            finite(value)?;
            if column != row && value != 0.0 {
                let a = root(&mut parent, row);
                let b = root(&mut parent, column);
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    let mut labels = vec![usize::MAX; n];
    let mut component = Vec::with_capacity(n);
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for vertex in 0..n {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(cx, vertex)?;
        }
        let r = root(&mut parent, vertex);
        if labels[r] == usize::MAX {
            labels[r] = groups.len();
            groups.push(Vec::new());
        }
        let label = labels[r];
        component.push(label);
        groups[label].push(vertex);
    }
    Ok((component, groups))
}

fn root(parent: &mut [usize], mut index: usize) -> usize {
    while parent[index] != index {
        parent[index] = parent[parent[index]];
        index = parent[index];
    }
    index
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn solve_component(
    cx: &Cx<'_>,
    storage: &EnthalpyBackwardEuler<'_, '_>,
    old: &[f64],
    vertices: &[usize],
    component: usize,
    config: IsolatedEquilibriumConfig,
    intervals: &mut [[f64; 2]],
) -> Result<IsolatedEquilibriumComponent, IsolatedEquilibriumError> {
    let mut lo = 0.0_f64;
    let mut hi = f64::INFINITY;
    let mut latent = Vec::new();
    // Addresses are used only for membership, never iteration order: identical
    // charts are scanned once per component in original vertex order.
    let mut seen = BTreeSet::new();
    for (at, &vertex) in vertices.iter().enumerate() {
        if at % ASSEMBLY_TILE == 0 {
            poll(cx, at)?;
        }
        let curve = storage.curve_for_vertex(vertex);
        let knots = curve.knots();
        lo = lo.max(knots[0].temperature_k);
        hi = hi.min(knots[knots.len() - 1].temperature_k);
        if let TemperatureSpan::Sampled { low, high } = storage.film_temperature_span(vertex) {
            lo = lo.max(low);
            hi = hi.min(high);
        }
        if seen.insert(std::ptr::from_ref(curve) as usize) {
            for (index, pair) in knots.windows(2).enumerate() {
                if index % ASSEMBLY_TILE == 0 {
                    poll(cx, index)?;
                }
                if pair[0].temperature_k == pair[1].temperature_k {
                    latent.push(pair[0].temperature_k);
                }
            }
        }
    }
    if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || lo > hi {
        return Err(IsolatedEquilibriumError::NoCommonTemperature { component });
    }
    latent.retain(|&t| t > lo && t < hi);
    latent.sort_by(f64::total_cmp);
    latent.dedup();
    let mut evaluations = 0;
    let mut evaluate = |temperature| -> Result<[f64; 2], IsolatedEquilibriumError> {
        poll(cx, evaluations)?;
        if evaluations >= config.max_temperature_evaluations {
            return Err(IsolatedEquilibriumError::WorkLimit {
                component,
                evaluations,
            });
        }
        evaluations += 1;
        evaluate_energy(cx, storage, old, vertices, temperature, intervals)
    };
    let at_lo = evaluate(lo)?;
    if accepts(at_lo, config.energy_tolerance_j) {
        return finish_component(vertices, lo, at_lo, evaluations, intervals);
    }
    let at_hi = evaluate(hi)?;
    if accepts(at_hi, config.energy_tolerance_j) {
        return finish_component(vertices, hi, at_hi, evaluations, intervals);
    }
    if at_lo[0] > 0.0 || at_hi[1] < 0.0 {
        return Err(IsolatedEquilibriumError::EnergyOutsideDomain {
            component,
            energy_change_interval_j: [at_lo[0], at_hi[1]],
        });
    }
    let (mut lo_energy, mut hi_energy) = (at_lo[1], at_hi[0]);
    let (mut left, mut right) = (0, latent.len());
    while left < right {
        let middle = left + (right - left) / 2;
        let temperature = latent[middle];
        let energy = evaluate(temperature)?;
        if accepts(energy, config.energy_tolerance_j) {
            return finish_component(vertices, temperature, energy, evaluations, intervals);
        }
        if energy[1] < 0.0 {
            lo = temperature;
            lo_energy = energy[1];
            left = middle + 1;
        } else {
            hi = temperature;
            hi_energy = energy[0];
            right = middle;
        }
    }
    loop {
        // Normalize energies before forming the secant so opposite large
        // finite values cannot overflow their difference. The safeguard gives
        // bounded interval contraction even for sharply varying film capacity.
        let scale = lo_energy.abs().max(hi_energy.abs());
        let fraction = if scale > 0.0 {
            let a = -lo_energy / scale;
            let b = hi_energy / scale;
            (a / (a + b)).clamp(0.125, 0.875)
        } else {
            0.5
        };
        let temperature = fraction.mul_add(hi - lo, lo);
        if temperature <= lo || temperature >= hi {
            return Err(IsolatedEquilibriumError::TemperatureResolution {
                component,
                bracket_k: [lo, hi],
            });
        }
        let energy = evaluate(temperature)?;
        if accepts(energy, config.energy_tolerance_j) {
            return finish_component(vertices, temperature, energy, evaluations, intervals);
        }
        if energy[1] < 0.0 {
            lo = temperature;
            lo_energy = energy[1];
        } else {
            hi = temperature;
            hi_energy = energy[0];
        }
    }
}

fn accepts(energy: [f64; 2], tolerance: f64) -> bool {
    (energy[0] <= 0.0 && energy[1] >= 0.0)
        || (energy[0] == energy[1] && energy[0].abs() <= tolerance)
}

fn finish_component(
    vertices: &[usize],
    temperature_k: f64,
    mut energy: [f64; 2],
    evaluations: usize,
    intervals: &mut [[f64; 2]],
) -> Result<IsolatedEquilibriumComponent, IsolatedEquilibriumError> {
    // If the required energy is exactly an extreme, positivity of every
    // original mass forces every latent coordinate to that extreme.
    if energy[0] == 0.0 || energy[1] == 0.0 {
        let side = usize::from(energy[0] != 0.0);
        for &vertex in vertices {
            let h = intervals[vertex][side];
            intervals[vertex] = [h, h];
        }
        energy = [0.0; 2];
    }
    // No mass averaging or hidden phase allocation is performed.
    Ok(IsolatedEquilibriumComponent {
        temperature_k,
        vertices: vertices.len(),
        energy_change_interval_j: energy,
        latent_energy_above_lower_j: finite(-energy[0])?,
        temperature_evaluations: evaluations,
    })
}

#[allow(clippy::too_many_arguments)]
fn evaluate_energy(
    cx: &Cx<'_>,
    storage: &EnthalpyBackwardEuler<'_, '_>,
    old: &[f64],
    vertices: &[usize],
    temperature: f64,
    intervals: &mut [[f64; 2]],
) -> Result<[f64; 2], IsolatedEquilibriumError> {
    let mut sums = [Dd::from_f64(0.0); 2];
    for (at, &vertex) in vertices.iter().enumerate() {
        if at % ASSEMBLY_TILE == 0 {
            poll(cx, at)?;
        }
        let range = enthalpy_interval(storage, vertex, temperature)?;
        intervals[vertex] = range;
        if range[1] > range[0]
            && finite(storage.masses[vertex] * finite(range[1] - range[0])?)? <= 0.0
        {
            return Err(invalid("positive latent inventory is not representable"));
        }
        for side in 0..2 {
            let actual = storage
                .curve_for_vertex(vertex)
                .values_at_specific_enthalpy(range[side])
                .map_err(|source| EnthalpyError::Phase { vertex, source })?
                .temperature_k;
            sums[side] =
                sums[side] + storage.storage_change_dd(vertex, range[side], old[vertex], actual)?;
        }
    }
    Ok([dd_value(sums[0])?, dd_value(sums[1])?])
}

fn enthalpy_interval(
    storage: &EnthalpyBackwardEuler<'_, '_>,
    vertex: usize,
    temperature: f64,
) -> Result<[f64; 2], IsolatedEquilibriumError> {
    let knots = storage.curve_for_vertex(vertex).knots();
    let upper = knots.partition_point(|knot| knot.temperature_k < temperature);
    if upper == knots.len() || temperature < knots[0].temperature_k {
        return Err(invalid("common temperature left an original caloric chart"));
    }
    if knots[upper].temperature_k == temperature {
        let end = knots.partition_point(|knot| knot.temperature_k <= temperature) - 1;
        return Ok([
            knots[upper].specific_enthalpy_j_kg,
            knots[end].specific_enthalpy_j_kg,
        ]);
    }
    let a = knots[upper - 1];
    let b = knots[upper];
    let fraction = (temperature - a.temperature_k) / (b.temperature_k - a.temperature_k);
    let h = finite(fraction.mul_add(
        b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg,
        a.specific_enthalpy_j_kg,
    ))?;
    Ok([h, h])
}

fn dd_value(value: Dd) -> Result<f64, IsolatedEquilibriumError> {
    finite(value.hi)?;
    finite(value.lo)?;
    Ok(finite(value.to_f64())?)
}
fn invalid(reason: &'static str) -> IsolatedEquilibriumError {
    EnthalpyError::InvalidInput(reason).into()
}
