//! Conservative finite-volume osmosis with pressure and dilution feedback.
//!
//! The original `fs_material::osmosis` constitutive law owns water activity,
//! selectivity and exponentially fitted solute transport. This module evolves
//! actual solvent volumes and solute amounts, including initially absent
//! species. Mechanical pressure is a separate, explicit volume response.
//!
//! Each membrane is a fully implicit two-compartment exchange. At a proposed
//! water transfer, the solute equations have an analytic conservative solution;
//! only one scalar water balance needs iteration. Ordered membrane sweeps are
//! first-order Lie splitting. One full sweep versus two half sweeps controls
//! local temporal error; the two-half-step endpoint is retained. No uniqueness
//! theorem is claimed for the nonlinear permeating-solute residual, nor is the
//! local error estimate a certified global error bound.
//!
//! Compartments are well mixed and isothermal, with negligible solute volume.
//! A common declared species ordering, solvent and temperature is mandatory.
//! Ideal selective membranes can use quasi-steady concentration films, which
//! reconstruct actual membrane-face concentrations without hidden inventories.
//! Nonideal compartments use their shared free-energy pressure with fully
//! retained solutes and bare membranes; leaky and film transport need matching
//! chemical-potential laws. Separate owners provide finite-storage films in
//! [`transient`] and ideal fast-ion water exchange in [`donnan`]. None supplies
//! solvent inertia, a heat balance or mechanical inertia.
//! Successful advances allocate no memory. All workspace is prepared up front;
//! failed, cancelled or exhausted advances publish neither state nor time.

/// Nonideal Maxwell–Stefan diffusion at fixed solvent volume on this timeline.
pub mod diffusion;

/// Finite-storage concentration films coupled to membrane water transport.
pub mod transient;

/// Finite water exchange with conservative, fast-ion Donnan equilibrium.
pub mod donnan;

#[cfg(test)]
mod nonideal_tests;

#[cfg(test)]
use fs_material::osmosis::IdealDiluteSolution;
use fs_material::osmosis::polarization::PolarizedMembraneTransport;
use fs_material::osmosis::{AqueousSolution, MembraneTransport, OsmosisError};

/// Maximum independently conserved dissolved particle species.
pub const MAX_OSMOTIC_SPECIES: usize = 64;
const MAX_COMPARTMENTS: usize = 65_536;
const MAX_MEMBRANES: usize = 262_144;

/// Shared constitutive refusal or finite-volume integration failure.
#[derive(Debug)]
pub enum OsmoticNetworkError {
    /// Malformed geometry, layout, limits or incompatible thermodynamic state.
    Invalid(&'static str),
    /// The original solution or membrane law refused the state.
    Material(OsmosisError),
    /// Source-owned solution thermodynamics or original MS friction refused.
    SolutionDiffusion(fs_material::osmosis::diffusion::SolutionDiffusionError),
    /// A smaller physical diffusion step may satisfy withdrawal or passivity.
    DiffusionStep(diffusion::SolutionDiffusionStepRefusal),
    /// Actual closed-window species inventory failed its absolute molar budget.
    DiffusionBalance { species: usize, residual_mol: f64 },
    /// Closed fixed-volume mixing free energy increased beyond its whole-window budget.
    DiffusionFreeEnergy { change_j: f64, maximum_j: f64 },
    /// Prepared numerical workspace exceeds the caller's byte budget.
    Workspace {
        required_bytes: usize,
        budget_bytes: usize,
    },
    /// No admissible scalar water root was bracketed at this step size.
    Unbracketed,
    /// Scalar solve exhausted its iteration bound.
    NonlinearConvergence,
    /// Time integration exhausted its accepted-step or rejection budget.
    WorkBudget,
    /// Further halving cannot advance floating-point time.
    StepUnderflow,
    /// Cancellation left accepted state and clock unchanged.
    Cancelled,
}

impl core::fmt::Display for OsmoticNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-volume osmosis: {self:?}")
    }
}
impl std::error::Error for OsmoticNetworkError {}
impl From<OsmosisError> for OsmoticNetworkError {
    fn from(value: OsmosisError) -> Self {
        Self::Material(value)
    }
}

/// Mechanical pressure `p_ref + elastance * (V - V_ref)`.
///
/// A zero elastance is a prescribed pressure. An open vertical reservoir of
/// constant cross-section A has elastance rho*g/A; an elastic chamber can
/// provide its tangent pressure/volume response. This pressure excludes
/// osmotic pressure, which the membrane law accounts for exactly once.
/// Pressures may use any COMMON reference; this law does not admit cavitation.
#[derive(Clone, Copy, Debug)]
pub struct OsmoticPressureVolume {
    reference_volume_m3: f64,
    reference_pressure_pa: f64,
    elastance_pa_per_m3: f64,
}

impl OsmoticPressureVolume {
    /// Admit a finite reference pressure, positive volume and nonnegative stiffness.
    pub fn new(
        reference_volume_m3: f64,
        reference_pressure_pa: f64,
        elastance_pa_per_m3: f64,
    ) -> Result<Self, OsmoticNetworkError> {
        if !positive(reference_volume_m3)
            || !reference_pressure_pa.is_finite()
            || !nonnegative(elastance_pa_per_m3)
        {
            return Err(OsmoticNetworkError::Invalid(
                "invalid mechanical pressure/volume law",
            ));
        }
        Ok(Self {
            reference_volume_m3,
            reference_pressure_pa,
            elastance_pa_per_m3,
        })
    }

    /// Mechanical pressure at the supplied solvent volume [Pa].
    pub fn pressure_pa(self, volume_m3: f64) -> Result<f64, OsmoticNetworkError> {
        let pressure = self.reference_pressure_pa
            + self.elastance_pa_per_m3 * (volume_m3 - self.reference_volume_m3);
        if !positive(volume_m3) || !pressure.is_finite() {
            return Err(OsmoticNetworkError::Invalid(
                "nonfinite pressure/volume response",
            ));
        }
        Ok(pressure)
    }
}

/// One finite, well-mixed solution and its mechanical volume response.
#[derive(Clone, Debug)]
pub struct OsmoticCompartment {
    solution: AqueousSolution,
    pressure: OsmoticPressureVolume,
    minimum_volume_m3: f64,
}

impl OsmoticCompartment {
    /// The explicit positive volume floor is an applicability boundary, never
    /// a clipped state or invented residual water inventory.
    pub fn new(
        solution: impl Into<AqueousSolution>,
        pressure: OsmoticPressureVolume,
        minimum_volume_m3: f64,
    ) -> Result<Self, OsmoticNetworkError> {
        let solution = solution.into();
        if !positive(minimum_volume_m3) || minimum_volume_m3 >= solution.solvent_volume_m3() {
            return Err(OsmoticNetworkError::Invalid(
                "volume floor must be below initial volume",
            ));
        }
        pressure.pressure_pa(solution.solvent_volume_m3())?;
        Ok(Self {
            solution,
            pressure,
            minimum_volume_m3,
        })
    }
    /// Actual finite solvent and dissolved inventories.
    pub const fn solution(&self) -> &AqueousSolution {
        &self.solution
    }
    /// Current mechanical pressure [Pa].
    pub fn pressure_pa(&self) -> Result<f64, OsmoticNetworkError> {
        self.pressure.pressure_pa(self.solution.solvent_volume_m3())
    }
}

/// One oriented membrane with a fixed physical area and the original material law.
#[derive(Clone, Debug)]
pub struct OsmoticMembrane {
    left: usize,
    right: usize,
    area_m2: f64,
    material: MembraneBoundary,
}

#[derive(Clone, Debug)]
enum MembraneBoundary {
    WellMixed(MembraneTransport),
    Polarized(PolarizedMembraneTransport),
}

impl OsmoticMembrane {
    /// Positive transfer is left to right. Species follow the network ordering.
    pub fn new(
        left: usize,
        right: usize,
        area_m2: f64,
        material: MembraneTransport,
    ) -> Result<Self, OsmoticNetworkError> {
        Self::from_boundary(left, right, area_m2, MembraneBoundary::WellMixed(material))
    }

    /// Add species-specific concentration films on either membrane face.
    ///
    /// Films have negligible storage and prescribed mass-transfer coefficients.
    /// The original membrane owns selectivity and hydraulic resistance. Its
    /// osmotic force uses reconstructed face concentrations, while the same
    /// species flux debits and credits the finite bulk inventories. Film viscous
    /// losses and transient boundary-layer storage are outside this model.
    pub fn new_polarized(
        left: usize,
        right: usize,
        area_m2: f64,
        material: PolarizedMembraneTransport,
    ) -> Result<Self, OsmoticNetworkError> {
        Self::from_boundary(left, right, area_m2, MembraneBoundary::Polarized(material))
    }

    fn from_boundary(
        left: usize,
        right: usize,
        area_m2: f64,
        material: MembraneBoundary,
    ) -> Result<Self, OsmoticNetworkError> {
        if left == right || !positive(area_m2) {
            return Err(OsmoticNetworkError::Invalid(
                "membrane needs distinct cells and positive area",
            ));
        }
        Ok(Self {
            left,
            right,
            area_m2,
            material,
        })
    }

    fn membrane_law(&self) -> &MembraneTransport {
        match &self.material {
            MembraneBoundary::WellMixed(material) => material,
            MembraneBoundary::Polarized(material) => material.membrane(),
        }
    }

    fn heap_bytes(&self) -> usize {
        match &self.material {
            MembraneBoundary::WellMixed(material) => core::mem::size_of_val(material.solutes()),
            MembraneBoundary::Polarized(material) => material.heap_bytes(),
        }
    }
}

/// Explicit temporal accuracy and work limits; all quantities use SI units.
#[derive(Clone, Copy, Debug)]
pub struct OsmoticStepControls {
    /// Largest attempted substep [s].
    pub maximum_step_s: f64,
    /// Relative local step-doubling tolerance, in (0,1).
    pub relative_tolerance: f64,
    /// Absolute solvent-volume tolerance [m³].
    pub absolute_volume_tolerance_m3: f64,
    /// Absolute per-species amount tolerance [mol].
    pub absolute_amount_tolerance_mol: f64,
    /// Maximum accepted substeps across the requested interval.
    pub maximum_steps: usize,
    /// Maximum total rejected trials across the requested interval.
    pub maximum_rejections: usize,
    /// Maximum bracket-search and safeguarded scalar iterations per exchange.
    pub maximum_root_iterations: usize,
}

impl Default for OsmoticStepControls {
    fn default() -> Self {
        Self {
            maximum_step_s: 1.0,
            relative_tolerance: 1e-5,
            absolute_volume_tolerance_m3: 1e-18,
            absolute_amount_tolerance_mol: 1e-18,
            maximum_steps: 10_000,
            maximum_rejections: 100,
            maximum_root_iterations: 96,
        }
    }
}

/// Numerical work and actual closed-system inventory balances.
#[derive(Clone, Debug)]
pub struct OsmoticStepReport {
    /// Original and accepted time [s].
    pub interval_s: [f64; 2],
    /// Accepted two-half-step endpoints.
    pub accepted_steps: usize,
    /// Trials rejected for source, root or error-control reasons.
    pub rejected_steps: usize,
    /// Scalar water residuals or solution-diffusion face evaluations, including rejected trials.
    pub residual_evaluations: usize,
    /// Largest accepted step-doubling error divided by its requested scale.
    pub maximum_scaled_error: f64,
    /// Final minus initial sum of compartment solvent volumes [m³].
    pub water_balance_residual_m3: f64,
    /// Final minus initial amount for each species [mol]; only the first
    /// `species_count` entries are active.
    pub solute_balance_residual_mol: [f64; MAX_OSMOTIC_SPECIES],
    /// Number of active entries in the amount balance.
    pub species_count: usize,
    /// Reference-cancelled original-to-final mixing free-energy change [J].
    /// Some only for fixed-volume solution diffusion; membrane water work is
    /// not silently dropped from a mixed mechanical/chemical energy account.
    pub diffusion_free_energy_change_j: Option<f64>,
}

/// Prepared closed network. Public compartments are read-only checkpoints;
/// cloning the network retains complete physical and numerical state.
#[derive(Clone, Debug)]
pub struct OsmoticNetwork {
    species: Vec<String>,
    cells: Vec<OsmoticCompartment>,
    membranes: Vec<OsmoticMembrane>,
    diffusion: Vec<diffusion::SolutionDiffusionLink>,
    diffusion_controls: Option<diffusion::SolutionDiffusionControls>,
    work: Vec<OsmoticCompartment>,
    coarse: Vec<OsmoticCompartment>,
    fine: Vec<OsmoticCompartment>,
    workspace_bytes: usize,
    time_s: f64,
}

impl OsmoticNetwork {
    /// Prepare all storage, binding one explicit ordering of independent
    /// dissolved particles. Named ions must already be separate species; this
    /// function never infers dissociation or a van't Hoff multiplier.
    pub fn new(
        species: Vec<String>,
        cells: Vec<OsmoticCompartment>,
        membranes: Vec<OsmoticMembrane>,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        if cells.is_empty()
            || cells.len() > MAX_COMPARTMENTS
            || membranes.len() > MAX_MEMBRANES
            || species.len() > MAX_OSMOTIC_SPECIES
            || species
                .iter()
                .enumerate()
                .any(|(i, name)| name.trim().is_empty() || species[..i].contains(name))
        {
            return Err(OsmoticNetworkError::Invalid(
                "invalid compartment, membrane or species layout",
            ));
        }
        let first = &cells[0].solution;
        if cells.iter().any(|cell| {
            cell.solution.solute_moles().len() != species.len()
                || cell.solution.temperature_k().to_bits() != first.temperature_k().to_bits()
                || cell.solution.solvent_molar_volume_m3_mol().to_bits()
                    != first.solvent_molar_volume_m3_mol().to_bits()
        }) || membranes.iter().any(|m| {
            m.left >= cells.len()
                || m.right >= cells.len()
                || m.membrane_law().solutes().len() != species.len()
        }) {
            return Err(OsmoticNetworkError::Invalid(
                "all cells need one solvent, temperature and species ordering",
            ));
        }
        for membrane in &membranes {
            let nonideal = !cells[membrane.left].solution.is_ideal()
                || !cells[membrane.right].solution.is_ideal();
            if nonideal
                && (!matches!(membrane.material, MembraneBoundary::WellMixed(_))
                    || membrane
                        .membrane_law()
                        .solutes()
                        .iter()
                        .any(|(reflection, permeance)| {
                            reflection.value() != 1.0 || permeance.value() != 0.0
                        }))
            {
                return Err(OsmoticNetworkError::Invalid(
                    "nonideal solutions require fully retained solutes and a bare membrane; leaking and film transport need matching chemical-potential laws",
                ));
            }
        }
        // Four complete state buffers; input excess capacities and immutable
        // material/species data are counted conservatively. Bounded stack
        // scratch in the scalar solve is additional (under 4 KiB).
        let cell_bytes = core::mem::size_of::<OsmoticCompartment>();
        let required = cells
            .capacity()
            .checked_add(3 * cells.len())
            .and_then(|n| n.checked_mul(cell_bytes))
            .and_then(|n| {
                let heaps = cells.iter().try_fold(0usize, |total, cell| {
                    total.checked_add(cell.solution.heap_bytes())
                })?;
                n.checked_add(heaps.checked_mul(4)?)
            })
            .and_then(|n| {
                n.checked_add(membranes.capacity() * core::mem::size_of::<OsmoticMembrane>())
            })
            .and_then(|n| n.checked_add(membranes.iter().map(OsmoticMembrane::heap_bytes).sum()))
            .and_then(|n| n.checked_add(species.capacity() * core::mem::size_of::<String>()))
            .and_then(|n| n.checked_add(species.iter().map(String::capacity).sum()))
            .ok_or(OsmoticNetworkError::Invalid(
                "workspace byte count overflow",
            ))?;
        if required > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: required,
                budget_bytes: maximum_workspace_bytes,
            });
        }
        let (water, amounts) = totals(&cells, species.len());
        if !water.is_finite() || amounts[..species.len()].iter().any(|n| !n.is_finite()) {
            return Err(OsmoticNetworkError::Invalid(
                "network inventory total is not representable",
            ));
        }
        Ok(Self {
            species,
            work: cells.clone(),
            coarse: cells.clone(),
            fine: cells.clone(),
            cells,
            membranes,
            diffusion: Vec::new(),
            diffusion_controls: None,
            workspace_bytes: required,
            time_s: 0.0,
        })
    }
    /// Declared ordered species names.
    pub fn species(&self) -> &[String] {
        &self.species
    }
    /// Current accepted solution compartments.
    pub fn compartments(&self) -> &[OsmoticCompartment] {
        &self.cells
    }
    /// Current accepted time [s].
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Admitted numerical and immutable-input storage estimate [bytes].
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Advance exactly the requested interval with conservative, positive
    /// exchange and bounded adaptive substeps. Membrane stages are implicit;
    /// fixed-solvent diffusion stages use explicit paired fluxes. Pass
    /// `|| cx.checkpoint().is_err()` to use the shared execution scope.
    /// Cancellation is polled per scalar iteration, membrane and publication.
    pub fn advance(
        &mut self,
        duration_s: f64,
        controls: OsmoticStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OsmoticStepReport, OsmoticNetworkError> {
        admit_controls(duration_s, controls)?;
        let end_time = self.time_s + duration_s;
        if !end_time.is_finite() || end_time <= self.time_s {
            return Err(OsmoticNetworkError::StepUnderflow);
        }
        poll(&mut cancelled)?;
        copy_state(&self.cells, &mut self.work)?;
        let (water_before, amounts_before) = totals(&self.cells, self.species.len());
        let mut report = OsmoticStepReport {
            interval_s: [self.time_s, end_time],
            accepted_steps: 0,
            rejected_steps: 0,
            residual_evaluations: 0,
            maximum_scaled_error: 0.0,
            water_balance_residual_m3: 0.0,
            solute_balance_residual_mol: [0.0; MAX_OSMOTIC_SPECIES],
            species_count: self.species.len(),
            diffusion_free_energy_change_j: None,
        };
        let mut elapsed = 0.0;
        let mut step = duration_s.min(controls.maximum_step_s);
        while elapsed < duration_s {
            poll(&mut cancelled)?;
            if report.accepted_steps >= controls.maximum_steps {
                return Err(OsmoticNetworkError::WorkBudget);
            }
            let remaining = duration_s - elapsed;
            step = step.min(remaining);
            if elapsed + step <= elapsed || 0.5 * step == 0.0 {
                return Err(OsmoticNetworkError::StepUnderflow);
            }
            copy_state(&self.work, &mut self.coarse)?;
            copy_state(&self.work, &mut self.fine)?;
            let trial = (|| {
                sweep_all(
                    &mut self.coarse,
                    &self.membranes,
                    &mut self.diffusion,
                    self.diffusion_controls,
                    step,
                    controls,
                    &mut report.residual_evaluations,
                    &mut cancelled,
                )?;
                sweep_all(
                    &mut self.fine,
                    &self.membranes,
                    &mut self.diffusion,
                    self.diffusion_controls,
                    0.5 * step,
                    controls,
                    &mut report.residual_evaluations,
                    &mut cancelled,
                )?;
                sweep_all(
                    &mut self.fine,
                    &self.membranes,
                    &mut self.diffusion,
                    self.diffusion_controls,
                    0.5 * step,
                    controls,
                    &mut report.residual_evaluations,
                    &mut cancelled,
                )
            })();
            let error = match trial {
                Ok(()) => scaled_error(&self.work, &self.coarse, &self.fine, controls),
                Err(OsmoticNetworkError::Cancelled) => return Err(OsmoticNetworkError::Cancelled),
                Err(OsmoticNetworkError::DiffusionStep(_)) => f64::INFINITY,
                Err(error) if self.diffusion_controls.is_some() => return Err(error),
                Err(_) => f64::INFINITY,
            };
            if !error.is_finite() || error > 1.0 {
                if report.rejected_steps >= controls.maximum_rejections {
                    return Err(OsmoticNetworkError::WorkBudget);
                }
                report.rejected_steps += 1;
                step *= 0.5;
                continue;
            }
            core::mem::swap(&mut self.work, &mut self.fine);
            report.accepted_steps += 1;
            report.maximum_scaled_error = report.maximum_scaled_error.max(error);
            elapsed = if step == remaining {
                duration_s
            } else {
                elapsed + step
            };
            // Fixed binary growth preserves deterministic decision boundaries;
            // a first-order method's local difference scales as step².
            if error < 0.125 {
                step = (2.0 * step).min(controls.maximum_step_s);
            }
        }
        let (water_after, amounts_after) = totals(&self.work, self.species.len());
        if !water_after.is_finite()
            || amounts_after[..self.species.len()]
                .iter()
                .any(|n| !n.is_finite())
        {
            return Err(OsmoticNetworkError::Invalid(
                "network inventory total is not representable",
            ));
        }
        report.water_balance_residual_m3 = water_after - water_before;
        for i in 0..self.species.len() {
            report.solute_balance_residual_mol[i] = amounts_after[i] - amounts_before[i];
        }
        if let Some(physical) = self.diffusion_controls {
            // Sum actual endpoint changes, avoiding subtraction of two large
            // inventory totals. The same absolute budget covers the complete
            // window, independent of accepted/rejected stage counts.
            report.solute_balance_residual_mol.fill(0.0);
            let mut correction = [0.0; MAX_OSMOTIC_SPECIES];
            for (old, next) in self.cells.iter().zip(&self.work) {
                poll(&mut cancelled)?;
                for i in 0..self.species.len() {
                    compensated_add(
                        &mut report.solute_balance_residual_mol[i],
                        &mut correction[i],
                        next.solution.solute_moles()[i] - old.solution.solute_moles()[i],
                    );
                }
            }
            for i in 0..self.species.len() {
                let residual = report.solute_balance_residual_mol[i];
                if !residual.is_finite() || residual.abs() > physical.amount_balance_tolerance_mol {
                    return Err(OsmoticNetworkError::DiffusionBalance {
                        species: i,
                        residual_mol: residual,
                    });
                }
            }
            let change = diffusion::free_energy_change(&self.cells, &self.work, &mut cancelled)?;
            if change > physical.free_energy_tolerance_j {
                return Err(OsmoticNetworkError::DiffusionFreeEnergy {
                    change_j: change,
                    maximum_j: physical.free_energy_tolerance_j,
                });
            }
            report.diffusion_free_energy_change_j = Some(change);
        }
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.cells, &mut self.work);
        self.time_s = end_time;
        Ok(report)
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), OsmoticNetworkError> {
    if cancelled() {
        Err(OsmoticNetworkError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_controls(duration: f64, c: OsmoticStepControls) -> Result<(), OsmoticNetworkError> {
    if !positive(duration)
        || !positive(c.maximum_step_s)
        || !positive(c.relative_tolerance)
        || c.relative_tolerance >= 1.0
        || !positive(c.absolute_volume_tolerance_m3)
        || !positive(c.absolute_amount_tolerance_mol)
        || c.maximum_steps == 0
        || c.maximum_root_iterations == 0
        || c.maximum_root_iterations > 256
    {
        return Err(OsmoticNetworkError::Invalid(
            "invalid temporal accuracy or work controls",
        ));
    }
    Ok(())
}
fn copy_state(
    from: &[OsmoticCompartment],
    to: &mut [OsmoticCompartment],
) -> Result<(), OsmoticNetworkError> {
    for (source, dest) in from.iter().zip(to) {
        dest.solution.replace_inventories(
            source.solution.solvent_volume_m3(),
            source.solution.solute_moles(),
        )?;
    }
    Ok(())
}
fn scaled_error(
    old: &[OsmoticCompartment],
    coarse: &[OsmoticCompartment],
    fine: &[OsmoticCompartment],
    controls: OsmoticStepControls,
) -> f64 {
    let mut error = 0.0_f64;
    for ((old, coarse), fine) in old.iter().zip(coarse).zip(fine) {
        let before = old.solution.solvent_volume_m3();
        let after = fine.solution.solvent_volume_m3();
        error = error.max(
            (after - coarse.solution.solvent_volume_m3()).abs()
                / (controls.absolute_volume_tolerance_m3
                    + controls.relative_tolerance * before.abs().max(after.abs())),
        );
        for ((before, coarse), after) in old
            .solution
            .solute_moles()
            .iter()
            .zip(coarse.solution.solute_moles())
            .zip(fine.solution.solute_moles())
        {
            error = error.max(
                (after - coarse).abs()
                    / (controls.absolute_amount_tolerance_mol
                        + controls.relative_tolerance * before.abs().max(after.abs())),
            );
        }
    }
    error
}
fn totals(cells: &[OsmoticCompartment], species: usize) -> (f64, [f64; MAX_OSMOTIC_SPECIES]) {
    let mut volume = 0.0;
    let mut vc = 0.0;
    let mut amounts = [0.0; MAX_OSMOTIC_SPECIES];
    let mut corrections = [0.0; MAX_OSMOTIC_SPECIES];
    for cell in cells {
        compensated_add(&mut volume, &mut vc, cell.solution.solvent_volume_m3());
        for i in 0..species {
            compensated_add(
                &mut amounts[i],
                &mut corrections[i],
                cell.solution.solute_moles()[i],
            );
        }
    }
    (volume, amounts)
}
fn compensated_add(sum: &mut f64, correction: &mut f64, value: f64) {
    let corrected = value - *correction;
    let next = *sum + corrected;
    *correction = (next - *sum) - corrected;
    *sum = next;
}

#[allow(clippy::too_many_arguments)]
fn sweep_all(
    cells: &mut [OsmoticCompartment],
    membranes: &[OsmoticMembrane],
    links: &mut [diffusion::SolutionDiffusionLink],
    physical: Option<diffusion::SolutionDiffusionControls>,
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    match physical {
        Some(physical) => diffusion::sweep(cells, links, dt, physical, evaluations, cancelled),
        None => sweep(cells, membranes, dt, controls, evaluations, cancelled),
    }
}

fn sweep(
    cells: &mut [OsmoticCompartment],
    membranes: &[OsmoticMembrane],
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    for membrane in membranes {
        poll(cancelled)?;
        let endpoint = exchange(
            &cells[membrane.left],
            &cells[membrane.right],
            membrane,
            dt,
            controls,
            evaluations,
            cancelled,
        )?;
        let count = cells[membrane.left].solution.solute_moles().len();
        cells[membrane.left]
            .solution
            .replace_inventories(endpoint.volumes[0], &endpoint.left[..count])?;
        cells[membrane.right]
            .solution
            .replace_inventories(endpoint.volumes[1], &endpoint.right[..count])?;
    }
    Ok(())
}

struct PairEndpoint {
    volumes: [f64; 2],
    left: [f64; MAX_OSMOTIC_SPECIES],
    right: [f64; MAX_OSMOTIC_SPECIES],
    residual_m3: f64,
}

/// Eliminate implicit solute transfers analytically at a proposed water transfer.
fn pair_endpoint(
    left: &OsmoticCompartment,
    right: &OsmoticCompartment,
    membrane: &OsmoticMembrane,
    area_time: f64,
    water: f64,
) -> Result<PairEndpoint, OsmoticNetworkError> {
    let mut endpoint = PairEndpoint {
        volumes: [
            left.solution.solvent_volume_m3() - water,
            right.solution.solvent_volume_m3() + water,
        ],
        left: [0.0; MAX_OSMOTIC_SPECIES],
        right: [0.0; MAX_OSMOTIC_SPECIES],
        residual_m3: 0.0,
    };
    let [vl, vr] = endpoint.volumes;
    if !positive(vl) || !positive(vr) || vl < left.minimum_volume_m3 || vr < right.minimum_volume_m3
    {
        return Err(OsmoticNetworkError::Unbracketed);
    }
    let mut cl = [0.0; MAX_OSMOTIC_SPECIES];
    let mut cr = [0.0; MAX_OSMOTIC_SPECIES];
    let count = left.solution.solute_moles().len();
    let volume_flux_m_s = water / area_time;
    for i in 0..count {
        // Each film/membrane exponential fit is evaluated once for this trial.
        // The cached positive reconstruction is reused after the conservative
        // implicit species update. No nested water solve or film inventory.
        let polarization = match &membrane.material {
            MembraneBoundary::WellMixed(_) => None,
            MembraneBoundary::Polarized(material) => {
                Some(material.coefficients(i, volume_flux_m_s)?)
            }
        };
        let (kl, kr) = match &polarization {
            Some(coefficients) => coefficients.donor_coefficients(),
            None => membrane
                .membrane_law()
                .solute_stream_coefficients(i, volume_flux_m_s)?,
        };
        let a = area_time * kl / vl;
        let b = area_time * kr / vr;
        if !nonnegative(a) || !nonnegative(b) {
            return Err(OsmoticNetworkError::Invalid(
                "implicit solute coefficient overflow",
            ));
        }
        let nl = left.solution.solute_moles()[i];
        let nr = right.solution.solute_moles()[i];
        let scale = 1.0_f64.max(a).max(b);
        let a = a / scale;
        let b = b / scale;
        let inverse_scale = 1.0 / scale;
        let denominator = inverse_scale + a + b;
        let transfer = (a * nl - b * nr) / denominator;
        // Use a shared signed transfer for resolved changes. When nearly all
        // of a donor transfers, subtraction would erase its small but finite
        // remainder (and possibly its significant concentration). Recover that
        // remainder from positive implicit weights, then debit the SAME amount
        // on the other side. No endpoint is independently clipped.
        if transfer > 0.5 * nl {
            endpoint.left[i] = ((inverse_scale + b) / denominator) * nl + (b / denominator) * nr;
            endpoint.right[i] = nr + (nl - endpoint.left[i]);
        } else if -transfer > 0.5 * nr {
            endpoint.right[i] = ((inverse_scale + a) / denominator) * nr + (a / denominator) * nl;
            endpoint.left[i] = nl + (nr - endpoint.right[i]);
        } else {
            endpoint.left[i] = nl - transfer;
            endpoint.right[i] = nr + transfer;
        }
        let bulk_left = endpoint.left[i] / vl;
        let bulk_right = endpoint.right[i] / vr;
        (cl[i], cr[i]) = match &polarization {
            Some(coefficients) => {
                coefficients.membrane_face_concentrations(bulk_left, bulk_right)?
            }
            None => (bulk_left, bulk_right),
        };
    }
    let volume_flux = if left.solution.is_ideal() && right.solution.is_ideal() {
        membrane
            .membrane_law()
            .water_flux(
                left.solution.temperature_k(),
                left.pressure.pressure_pa(vl)?,
                right.pressure.pressure_pa(vr)?,
                &cl[..count],
                &cr[..count],
            )?
            .volume_flux_m_s
    } else {
        // Admission guarantees complete retention and no dilute Fick films.
        // The same water chemical potential supplies pressure and activity;
        // independent ionic/solute leakage is not inferred from this force.
        let osmotic_drop = left
            .solution
            .osmotic_pressure_at(vl, &endpoint.left[..count])?
            - right
                .solution
                .osmotic_pressure_at(vr, &endpoint.right[..count])?;
        membrane.membrane_law().hydraulic_permeance().value()
            * (left.pressure.pressure_pa(vl)? - right.pressure.pressure_pa(vr)? - osmotic_drop)
    };
    endpoint.residual_m3 = water - area_time * volume_flux;
    if !endpoint.residual_m3.is_finite() {
        return Err(OsmoticNetworkError::Invalid(
            "nonfinite implicit water residual",
        ));
    }
    Ok(endpoint)
}

fn exchange(
    left: &OsmoticCompartment,
    right: &OsmoticCompartment,
    membrane: &OsmoticMembrane,
    dt: f64,
    controls: OsmoticStepControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PairEndpoint, OsmoticNetworkError> {
    let area_time = dt * membrane.area_m2;
    if !positive(area_time) {
        return Err(OsmoticNetworkError::StepUnderflow);
    }
    // The solve tolerance is tighter than temporal error control. Relative
    // residual scaling follows transfer, not compartment size: tiny exchange
    // must not disappear merely because the reservoirs are large.
    let atol = controls.absolute_volume_tolerance_m3 * 0.01;
    let rtol = controls.relative_tolerance.min(1e-8) * 0.01;
    let mut evaluate = |water| {
        poll(cancelled)?;
        *evaluations += 1;
        pair_endpoint(left, right, membrane, area_time, water)
    };
    let at_zero = evaluate(0.0)?;
    if at_zero.residual_m3.abs() <= atol {
        return Ok(at_zero);
    }
    // Continue outwards only along the initial residual's flow direction.
    // Starting at a distant volume floor can select a remote root or overflow
    // the constitutive law even when a nearby physical root is admissible.
    let direction = if at_zero.residual_m3 < 0.0 { 1.0 } else { -1.0 };
    let available = if direction > 0.0 {
        left.solution.solvent_volume_m3() - left.minimum_volume_m3
    } else {
        right.solution.solvent_volume_m3() - right.minimum_volume_m3
    }
    .next_down();
    let mut distance = at_zero.residual_m3.abs().min(0.125 * available);
    let mut previous = 0.0;
    let mut previous_residual = at_zero.residual_m3;
    let mut refused_bound = None;
    let mut bracket = None;
    for _ in 0..controls.maximum_root_iterations {
        if distance <= previous || distance > available {
            return Err(OsmoticNetworkError::Unbracketed);
        }
        match evaluate(direction * distance) {
            Ok(endpoint) => {
                let residual = endpoint.residual_m3;
                if residual.abs() <= atol + rtol * distance {
                    return Ok(endpoint);
                }
                if direction * residual >= 0.0 {
                    bracket = Some(if direction > 0.0 {
                        (previous, distance, previous_residual, residual)
                    } else {
                        (-distance, -previous, residual, previous_residual)
                    });
                    break;
                }
                previous = distance;
                previous_residual = residual;
                distance = refused_bound.map_or_else(
                    || (2.0 * distance).min(available),
                    |upper| distance + 0.5 * (upper - distance),
                );
            }
            Err(OsmoticNetworkError::Cancelled) => return Err(OsmoticNetworkError::Cancelled),
            Err(_) => {
                refused_bound = Some(distance);
                distance = previous + 0.5 * (distance - previous);
            }
        }
    }
    let (mut lo, mut hi, mut flo, mut fhi) = bracket.ok_or(OsmoticNetworkError::Unbracketed)?;
    let mut bisect_next = false;
    for _ in 0..controls.maximum_root_iterations {
        let width = hi - lo;
        let midpoint = lo + 0.5 * width;
        // A narrow sign-changing bracket controls the transferred volume even
        // when a steep polarized osmotic force amplifies floating-point noise
        // in the residual. This is a state-error stopping rule, not permission
        // to publish an arbitrary last iterate after work exhaustion.
        if width <= 2.0 * (atol + rtol * midpoint.abs()) {
            return evaluate(midpoint);
        }
        // Scaled secant weights avoid overflow in fhi-flo. A midpoint follows
        // any interpolation that failed to halve the bracket, so interpolation
        // cannot stagnate at one endpoint. Keep the near-origin bracket itself:
        // permeating-solute models do not have a global unique-root guarantee.
        let scale = (-flo).max(fhi);
        let low_weight = -flo / scale;
        let high_weight = fhi / scale;
        let fraction = low_weight / (low_weight + high_weight);
        let interpolate = !bisect_next && fraction > 0.01 && fraction < 0.99;
        let water = if interpolate {
            lo + fraction * width
        } else {
            midpoint
        };
        let endpoint = evaluate(water)?;
        let residual = endpoint.residual_m3;
        if residual.abs() <= atol + rtol * water.abs() {
            return Ok(endpoint);
        }
        if water == lo || water == hi {
            return Err(OsmoticNetworkError::NonlinearConvergence);
        }
        if residual < 0.0 {
            lo = water;
            flo = residual;
        } else {
            hi = water;
            fhi = residual;
        }
        bisect_next = interpolate && hi - lo > 0.5 * width;
        // Both endpoint signs stay explicit; no unbracketed Newton jump is used.
        debug_assert!(flo <= 0.0 && fhi >= 0.0);
    }
    Err(OsmoticNetworkError::NonlinearConvergence)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod polarization_tests;
