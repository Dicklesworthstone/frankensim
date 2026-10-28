//! Spatial Ohmic conduction and local material-temperature/phase feedback.
//!
//! The original FEEC field owner solves div(sigma(T,f_liquid) grad(V))=0 at each
//! proposed thermal endpoint. Its per-triangle E.J integrals feed the existing
//! backward-Euler enthalpy residual. No uniform-resistor approximation, second
//! heat integrator, or independent material-name lookup is introduced.

/// Spatial conductor terminals coupled to native circuit storage and finite heat.
pub mod circuit;
/// Retained dielectric charge and local conduction heat on this same material network.
pub mod dielectric;
// The complete sourced heat network enters this same electrical/thermal residual.
mod heat;

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, DriverAccount, ThermalDriver, cell_driver_candidate, finite,
    magnetic::same_material,
    material::{
        ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
    },
};
use fs_em::conduction::{
    ConductionControls, ConductionCurrentTrace, ConductionError, ConductionField, PlanarConduction,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    conductor::planar::PlanarConductorMaterial,
    conductor::{ConductorError, ResolvedConductorState},
    thermal::source::ThermalMaterialError,
};
use fs_solver::SolveReport;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// Held ideal electrical boundary for one electrothermal interval.
#[derive(Clone, Copy, Debug)]
pub enum ConductionDrive<'a> {
    /// Original prescribed exterior potentials [V]; other edges are insulated.
    Voltage(&'a [(usize, f64)]),
    /// Total current [A] entering the positive electrode and leaving the
    /// negative electrode. Each electrode is a disjoint union of complete
    /// exterior edges. Their potentials are solved from the current material
    /// field; the negative electrode is the zero-potential reference.
    Current {
        positive: &'a [usize],
        negative: &'a [usize],
        current_a: f64,
    },
    /// Distributed prescribed inward current density on exterior edges, with
    /// voltage traces supplying the reference and any return current. These
    /// current traces need not be equipotential. Native edge integration and
    /// the spatial field determine the source work and local heat together.
    CurrentDensity {
        voltage: &'a [(usize, f64)],
        currents: &'a [ConductionCurrentTrace],
    },
}

impl ConductionDrive<'_> {
    fn solve<'p, C: PlanarConductorMaterial>(
        self,
        problem: &'p PlanarConduction<C>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionField<'p, C>, Error> {
        match self {
            Self::Voltage(boundary) => problem
                .solve(boundary, controls, cancelled)
                .map_err(field_error),
            Self::CurrentDensity { voltage, currents } => problem
                .solve_with_current_traces(voltage, currents, controls, cancelled)
                .map_err(field_error),
            Self::Current {
                positive,
                negative,
                current_a,
            } => {
                if !current_a.is_finite() {
                    return Err(invalid("prescribed total electrode current must be finite"));
                }
                // The native two-terminal field extraction includes the exact
                // geometry and all regional sigma(T). Its unit solution also
                // admits both electrode traces, including at zero current.
                let resistance = problem
                    .resistance(positive, negative, controls, cancelled)
                    .map_err(field_error)?;
                let voltage = finite(current_a * resistance.resistance_ohm())?;
                if current_a != 0.0 && voltage == 0.0 {
                    return Err(invalid(
                        "prescribed current lost its representable electrode voltage",
                    ));
                }
                let field = resistance
                    .field_at_voltage(
                        voltage,
                        controls.current_tolerance_a,
                        controls.power_tolerance_w,
                        cancelled,
                    )
                    .map_err(field_error)?;
                // Check the actual recovered current independently of I*R.
                let mut actual = 0.0;
                for &vertex in positive {
                    if cancelled() {
                        return Err(CircuitEnthalpyError::Cancelled);
                    }
                    let current = field.nodal_boundary_current_a()[vertex].ok_or_else(|| {
                        invalid("current electrode lost its native terminal reaction")
                    })?;
                    actual = finite(actual + current)?;
                }
                if finite(actual - current_a)?.abs() > controls.current_tolerance_a {
                    return Err(invalid(
                        "recovered electrode current does not match its prescribed total",
                    ));
                }
                Ok(field)
            }
        }
    }
}

/// One explicit receiving thermal cell for an original electrical triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConductionRegion {
    pub triangle: usize,
    pub cell: CellId,
}

/// Complete original material/heat correspondence, with no electrical state of
/// its own. Several triangles may cover one thermal cell; every electrical
/// triangle has exactly one receiver and mapped cells must be covered in full.
pub struct ConductionThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub regions: &'a [ConductionRegion],
    /// Absolute agreement allowance [m3] for covered volume and mass/rho.
    pub volume_tolerance_m3: f64,
}

#[derive(Debug)]
pub enum ConductionThermalError {
    Material(MaterialNetworkError),
    Heat(super::contacts::heat::HeatNetworkError),
    Conductor(ConductorError),
    Field(ConductionError),
    Dielectric(fs_em::electrostatic::relaxation::DielectricRelaxationError),
    Invalid(&'static str),
}
impl core::fmt::Display for ConductionThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial conduction/enthalpy: {self:?}")
    }
}
impl std::error::Error for ConductionThermalError {}
type Error = CircuitEnthalpyError<ConductionThermalError>;
fn invalid(reason: &'static str) -> Error {
    CircuitEnthalpyError::Source(ConductionThermalError::Invalid(reason))
}
fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        error => CircuitEnthalpyError::Source(ConductionThermalError::Material(error)),
    }
}
fn field_error(error: ConductionError) -> Error {
    match error {
        ConductionError::Cancelled => CircuitEnthalpyError::Cancelled,
        error => CircuitEnthalpyError::Source(ConductionThermalError::Field(error)),
    }
}
fn check(cancelled: &RefCell<impl FnMut() -> bool>) -> Result<(), Error> {
    if (*cancelled.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}

/// Observations copied from the ORIGINAL accepted spatial solve. Electrical
/// conduction has no stored energy on this quasistatic lane. Source work is
/// independently recovered terminal power integrated by backward Euler.
#[derive(Clone, Debug)]
pub struct ConductionThermalStep {
    pub start_s: f64,
    pub end_s: f64,
    pub voltage_v: Vec<f64>,
    pub electric_field_v_m: Vec<[f64; 2]>,
    pub current_density_a_m2: Vec<[f64; 2]>,
    pub cell_power_w: Vec<f64>,
    /// Positive INTO the conductor at each prescribed boundary vertex.
    pub nodal_boundary_current_a: Vec<Option<f64>>,
    /// Authored current-density trace loads integrated by the native field
    /// owner [A]. At a voltage/current corner these are included in the combined
    /// boundary current above, never an additional external energy supply.
    pub nodal_applied_current_a: Vec<Option<f64>>,
    pub boundary_current_a: f64,
    pub maximum_free_current_residual_a: f64,
    pub power_balance_residual_w: f64,
    pub heat_j: f64,
    pub supplied_work_j: f64,
    pub linear: Option<SolveReport>,
}

struct FieldCandidate<C> {
    problem: PlanarConduction<C>,
    step: ConductionThermalStep,
}
impl<C: PlanarConductorMaterial> ThermalDriver for FieldCandidate<C> {
    type Model = PlanarConduction<C>;
    type Step = ConductionThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.problem, self.step)
    }
}

/// Jointly accepted heat history and its endpoint spatial material operator.
/// Consume them together to continue; all original inputs remain unchanged.
pub struct ConductionEnthalpyCandidate<'t, C = ResolvedConductorState, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    problem: PlanarConduction<C>,
    report: CircuitEnthalpyReport<S, ConductionThermalStep>,
}
impl<'t, C: PlanarConductorMaterial, S> ConductionEnthalpyCandidate<'t, C, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn problem(&self) -> &PlanarConduction<C> {
        &self.problem
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<S, ConductionThermalStep> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        PlanarConduction<C>,
        CircuitEnthalpyReport<S, ConductionThermalStep>,
    ) {
        (self.thermal, self.problem, self.report)
    }
}

impl ConductionThermalNetwork<'_> {
    /// Original prescribed-voltage entry point, sharing the same complete
    /// thermal/material solve with current-controlled electrodes.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        boundary_voltage_v: &[(usize, f64)],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionEnthalpyCandidate<'t, C>, Error> {
        self.step_driven_candidate(
            thermal,
            target_s,
            problem,
            ConductionDrive::Voltage(boundary_voltage_v),
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Solve actual spatial current redistribution, local heating and thermal
    /// transport at ONE common endpoint. Either the electrode voltages or the
    /// total two-terminal current or distributed current density is held for
    /// the interval. Total-current drive uses
    /// the native material/geometry resistance extraction at EVERY trial; its
    /// voltage and local heat therefore evolve with the same accepted material.
    /// The existing enthalpy owner supplies the nonlinear residual,
    /// event rules, finite-difference/Newton budgets and final energy acceptance.
    /// `maximum_circuit_trials` bounds these field solves, including the final
    /// fresh solve; `controls.circuit` is unused because there is no circuit DAE.
    ///
    /// Source choices, cards, coordinate quantities and fixed external values
    /// remain pinned. Each binding's temperature and explicitly declared liquid
    /// mass fraction evolve together. Stationary, confined solid/mushy/liquid
    /// charts are admitted only when their density matches the fixed electrical
    /// volume throughout the complete chart; a non-solid chart requires the
    /// liquid-fraction coordinate even on an isothermal latent-heat plateau.
    /// The supplied electrical law owns the effective phase response: no mixing
    /// rule or phase-interface impedance is inferred. Geometry motion, phase
    /// density changes, thermoelectric transfer, induction, contacts and dielectric
    /// storage require their corresponding coupled owners. This is first-order
    /// thermal/electrical time discretization, not unconditional stability or
    /// measured-material validation. Refine mesh, time and solve tolerances.
    /// Scalar, oriented orthotropic and mixed cells share this same residual;
    /// every principal electrical source and its fixed orientation is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn step_driven_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        drive: ConductionDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        self.admit_step(
            thermal,
            target_s,
            problem,
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |source| source,
            &cancel,
        )
    }

    // Admission precedes material preparation and is shared with complete heat.
    fn admit_step<C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'_>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        assembly: WeightedAssemblyLimits,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<f64, Error> {
        super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let dt = finite(target_s - thermal.time_s())?;
        if dt <= 0.0 || self.cells.len() != thermal.cells().len() {
            return Err(invalid(
                "positive interval and complete thermal material bindings required",
            ));
        }
        if problem.mesh().faces().len() > assembly.max_cells
            || problem.mesh().vertices().len() > assembly.max_dofs
            || problem.mesh().edges().len() > assembly.max_dofs
            || problem
                .mesh()
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > assembly.max_triplets)
        {
            return Err(invalid(
                "electrical geometry exceeds the declared assembly capacity",
            ));
        }
        Ok(dt)
    }

    // One spatial current, source refresh, Joule deposition and nonlinear owner.
    // Both concrete public entry points retain their original source samples.
    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, C: PlanarConductorMaterial, S>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        drive: ConductionDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'_>) -> Result<CircuitThermalCoefficients<S>, Error>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<ConductionEnthalpyCandidate<'t, C, S>, Error> {
        check(cancel)?;
        let dt = self.admit_step(
            thermal,
            target_s,
            problem,
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
        let routes = self.admit_routes(thermal, problem, cancel)?;
        let accepted = resolve(thermal)?;
        let initial = self.materials_at(problem, sample(&accepted.source), &routes, cancel)?;
        if initial.iter().zip(problem.conductors()).any(|(a, b)| {
            a.source_basis() != b.source_basis()
                || a.frame() != b.frame()
                || !same_material(a.material(), b.material())
        }) {
            return Err(invalid(
                "accepted conductor sources must match the actual thermal state",
            ));
        }
        let (thermal, problem, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check(cancel)?;
                let coefficients = resolve(state)?;
                let sources =
                    self.materials_at(problem, sample(&coefficients.source), &routes, cancel)?;
                let next =
                    PlanarConduction::new(problem.mesh().clone(), sources, assembly, &mut || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(field_error)?;
                let field = drive.solve(&next, field_controls, &mut || (*cancel.borrow_mut())())?;
                let heat_j = finite(dt * field.power_w())?;
                let supplied_work_j =
                    finite(dt * finite(field.power_w() - field.power_balance_residual_w())?)?;
                if (field.power_w() > 0.0 && heat_j == 0.0)
                    || (heat_j - supplied_work_j).abs() > controls.total_energy_tolerance_j
                {
                    return Err(invalid(
                        "spatial Joule heat and terminal source work do not close",
                    ));
                }
                let step = ConductionThermalStep {
                    start_s: thermal.time_s(),
                    end_s: target_s,
                    voltage_v: field.voltage_v().to_vec(),
                    electric_field_v_m: field.electric_field_v_m().to_vec(),
                    current_density_a_m2: field.current_density_a_m2().to_vec(),
                    cell_power_w: field.cell_power_w().to_vec(),
                    nodal_boundary_current_a: field.nodal_boundary_current_a().to_vec(),
                    nodal_applied_current_a: field.nodal_applied_current_a().to_vec(),
                    boundary_current_a: field.boundary_current_a(),
                    maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                    power_balance_residual_w: field.power_balance_residual_w(),
                    heat_j,
                    supplied_work_j,
                    linear: field.report().cloned(),
                };
                Ok(DrivenCoefficients {
                    circuit: FieldCandidate {
                        problem: next,
                        step,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |candidate, _dt, checkpoint| {
                let mut powers = BTreeMap::new();
                for (triangle, &power) in candidate.step.cell_power_w.iter().enumerate() {
                    checkpoint()?;
                    let cell = routes[triangle];
                    let previous = *powers.get(&cell).unwrap_or(&0.0);
                    let combined = finite(previous + power)?;
                    if power != 0.0 && combined == previous {
                        return Err(invalid(
                            "local conductor heat disappeared in receiving-cell sum",
                        ));
                    }
                    powers.insert(cell, combined);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel)?;
        Ok(ConductionEnthalpyCandidate {
            thermal,
            problem,
            report,
        })
    }

    fn admit_routes<C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'_>,
        problem: &PlanarConduction<C>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<CellId>, Error> {
        if !self.volume_tolerance_m3.is_finite()
            || self.volume_tolerance_m3 < 0.0
            || self.regions.len() != problem.mesh().faces().len()
        {
            return Err(invalid(
                "one heat receiver per electrical triangle and a finite volume tolerance required",
            ));
        }
        let mut routes = vec![None; self.regions.len()];
        let mut volumes = BTreeMap::new();
        for region in self.regions {
            check(cancel)?;
            let target = routes
                .get_mut(region.triangle)
                .ok_or_else(|| invalid("foreign electrical triangle"))?;
            if target.replace(region.cell).is_some() || !thermal.cells().contains_key(&region.cell)
            {
                return Err(invalid(
                    "duplicate electrical triangle or inactive thermal receiver",
                ));
            }
            let volume = problem
                .mesh()
                .face_measure(region.triangle)
                .ok_or_else(|| invalid("electrical triangle lacks physical volume"))?;
            let sum = volumes.entry(region.cell).or_insert(0.0);
            *sum = finite(*sum + volume)?;
        }
        for (id, volume) in volumes {
            check(cancel)?;
            let cell = &thermal.cells()[&id];
            if (volume - cell.geometry.volume_m3).abs() > self.volume_tolerance_m3 {
                return Err(invalid(
                    "electrical triangles must cover the whole receiving thermal volume",
                ));
            }
            let binding = self
                .cells
                .iter()
                .find(|s| s.cell == id)
                .ok_or_else(|| invalid("electrical receiver lacks a caloric binding"))?;
            for knot in binding.curve.knots() {
                check(cancel)?;
                if knot.liquid_mass_fraction != 0.0 && binding.point.liquid_fraction_axis.is_none()
                {
                    return Err(invalid(
                        "non-solid electrical charts require an explicit liquid mass fraction coordinate",
                    ));
                }
                if (finite(cell.mass_kg / knot.bulk_density_kg_m3)? - cell.geometry.volume_m3).abs()
                    > self.volume_tolerance_m3
                {
                    return Err(invalid(
                        "fixed electrical geometry requires matching constant density through every phase",
                    ));
                }
            }
        }
        routes
            .into_iter()
            .map(|r| r.ok_or_else(|| invalid("electrical triangle has no heat receiver")))
            .collect()
    }

    fn materials_at<C: PlanarConductorMaterial>(
        &self,
        problem: &PlanarConduction<C>,
        sample: &MaterialNetworkSample,
        routes: &[CellId],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<C>, Error> {
        let bindings: BTreeMap<_, _> = self.cells.iter().map(|c| (c.cell, c)).collect();
        let phases: BTreeMap<_, _> = sample.cells.iter().map(|c| (c.cell, c.phase)).collect();
        let mut sources = Vec::with_capacity(routes.len());
        for (&cell, original) in routes.iter().zip(problem.conductors()) {
            check(cancel)?;
            let binding = bindings
                .get(&cell)
                .ok_or_else(|| invalid("electrical receiver lacks a material binding"))?;
            if original.liquid_mass_fraction_axis().is_some()
                && original.liquid_mass_fraction_axis() != binding.point.liquid_fraction_axis
            {
                return Err(invalid(
                    "electrical phase coordinate must be driven by the same caloric liquid mass fraction",
                ));
            }
            let phase = *phases
                .get(&cell)
                .ok_or_else(|| invalid("electrical receiver lacks endpoint phase"))?;
            if !phase.is_fully_solid() && binding.point.liquid_fraction_axis.is_none() {
                return Err(invalid(
                    "non-solid electrical states require an explicit liquid mass fraction coordinate",
                ));
            }
            let point = binding
                .point
                .at(phase)
                .map_err(|e| material(MaterialNetworkError::Material(e)))?;
            let snapshot = original
                .material()
                .at_query_point(binding.point.material, &point)
                .map_err(|e| {
                    material(MaterialNetworkError::Material(
                        ThermalMaterialError::Material(e),
                    ))
                })?;
            let conductor = original
                .with_material(snapshot)
                .map_err(|e| CircuitEnthalpyError::Source(ConductionThermalError::Conductor(e)))?;
            for property in conductor.source_properties() {
                binding
                    .point
                    .admit_static(property)
                    .map_err(|e| material(MaterialNetworkError::Material(e)))?;
            }
            sources.push(conductor);
        }
        Ok(sources)
    }
}
