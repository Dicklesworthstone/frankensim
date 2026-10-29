//! Conservative heat-driven condensed-phase decomposition.
//!
//! The same deposition cells, heat links and boundary operator carry a
//! composition-dependent virgin/char caloric mixture. The reaction is
//! `virgin -> y_char char + (1-y_char) volatile` by mass. Its finite-rate
//! backward-Euler extent and temperature close one local energy equation;
//! product formation enthalpy is carried once, never added as another heater.
//! The base owner has fixed geometry and instantaneous removal of volatiles at
//! their local production temperature. The `implicit` extension couples endpoint
//! conduction and reaction without the explicit conduction timestep limit.
//! The `slab` extension adds homogeneous noncharring recession and conservative
//! enthalpy rezoning at an explicitly supplied constant condensed density.
//! This reduced first-order model does not model pore
//! gas transport, oxidation, cracking, flame chemistry or a validated
//! PMMA/wood parameter set. Full cell exhaustion and exhausted work budgets
//! refuse instead of retaining an artificial seed mass.
//!
//! Governing-law source: FDS Technical Reference Guide, Solid Phase chapter,
//! "Solid Phase Reaction Rates" and the solid mass/energy conservation appendix.
//! Pinned chapter: <https://github.com/firemodels/fds/blob/a26081f3b24915552fab7ebf078185af817d9043/Manuals/FDS_Technical_Reference_Guide/Solid_Chapter.tex>.
//! Caller-owned material data, temperature domains and calibration remain
//! explicit. The endothermic single-step model is an estimate, not a measured
//! specimen or a general nonlinear/entropy stability certificate.

use super::{
    Boundary, CellId, DepositedCell, DepositionError, DepositionSimulation, EnergyAudit, Event,
    MaterialId, ThermalLaw, finite, nonnegative, positive,
};
use std::collections::BTreeMap;

/// Coupled backward-Euler heat transport and finite-rate decomposition.
#[path = "pyrolysis/implicit.rs"]
pub mod implicit;

/// Conservative shrinking and rezoning of homogeneous noncharring slabs.
#[path = "pyrolysis/slab.rs"]
pub mod slab;

/// SI universal molar gas constant, J/(mol K).
const GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// Original condensed caloric law with its source-owned inverse.
///
/// `h(T)` must be continuous and strictly increasing throughout the declared
/// material domain and invert the ORIGINAL `ThermalLaw::temperature_k`.
/// Virgin and char may have distinct, nonlinear heat capacities and formation
/// datums. Latent plateaus need a separate phase state and are not admitted here.
pub trait CondensedCaloric: ThermalLaw {
    /// Specific enthalpy in J/kg on the same original energy reference.
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, DepositionError>;
}

/// Product enthalpy carried out with each kilogram of released volatile.
///
/// The chart must be finite, continuous and nondecreasing on the material's
/// temperature domain and use the same formation reference as the solids.
/// A gas species/NASA owner can implement this directly; no gas Cp is invented.
/// All charts remain immutable across candidates and accepted advances.
pub trait VolatileEnthalpy {
    /// Carried specific enthalpy, J/kg, including its formation reference.
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, DepositionError>;
}
impl<F: Fn(f64) -> Result<f64, DepositionError>> VolatileEnthalpy for F {
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, DepositionError> {
        self(temperature_k)
    }
}

/// First-order, non-oxidative Arrhenius kinetics; its prefactor has units 1/s.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arrhenius {
    /// Nonnegative pre-exponential rate, 1/s; zero explicitly disables reaction.
    pub prefactor_s_inv: f64,
    /// Nonnegative activation energy, J/mol.
    pub activation_energy_j_mol: f64,
}
impl Arrhenius {
    /// Evaluate the actual first-order rate, preserving the shared math owner.
    pub fn rate_constant_s_inv(&self, temperature_k: f64) -> Result<f64, DepositionError> {
        if !nonnegative(self.prefactor_s_inv)
            || !nonnegative(self.activation_energy_j_mol)
            || !positive(temperature_k)
        {
            return Err(DepositionError::InvalidInput(
                "invalid first-order Arrhenius law",
            ));
        }
        let activation_temperature = self.activation_energy_j_mol / GAS_CONSTANT;
        finite(self.prefactor_s_inv * fs_math::det::exp(-activation_temperature / temperature_k))
    }

    /// The SAME backward-Euler inventory update used by the coupled solid.
    /// This isothermal query does not supply its thermostat's required heat.
    /// Small rate-duration products are applied as stable mass decrements.
    /// Consumption is then reconciled to the ACTUAL remaining inventory, so
    /// an unrepresentable decrement cannot emit fictitious volatile material.
    pub fn isothermal_step(
        &self,
        virgin_mass_kg: f64,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<KineticStep, DepositionError> {
        if !nonnegative(virgin_mass_kg) || !positive(duration_s) {
            return Err(DepositionError::InvalidInput(
                "invalid pyrolysis inventory or duration",
            ));
        }
        let rate = self.rate_constant_s_inv(temperature_k)?;
        let scaled = finite(duration_s * rate)?;
        let requested_extent = virgin_mass_kg * (scaled / (1.0 + scaled));
        let remaining = if scaled <= 1.0 {
            // 1+scaled can round to 1 even when this mass decrement is
            // representable; do not erase that physical reaction by division.
            virgin_mass_kg - requested_extent
        } else {
            // Avoid subtracting nearly the entire inventory for a stiff step.
            virgin_mass_kg / (1.0 + scaled)
        };
        if virgin_mass_kg > 0.0 && remaining == 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        let consumed = virgin_mass_kg - remaining;
        Ok(KineticStep {
            virgin_mass_kg: remaining,
            consumed_mass_kg: consumed,
            inventory_roundoff_kg: consumed - requested_extent,
        })
    }
}

/// Nonnegative finite-rate endpoint; no independent extent/heat clipping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticStep {
    /// Unreacted virgin material, kg.
    pub virgin_mass_kg: f64,
    /// Virgin material converted during the requested interval, kg.
    pub consumed_mass_kg: f64,
    /// Actual consumed minus stably computed backward-Euler extent, kg.
    /// This numerical discrepancy is neither emitted mass nor reaction heat.
    pub inventory_roundoff_kg: f64,
}

/// One immutable material closure, bound to the timeline's original virgin law.
#[derive(Clone, Copy)]
pub struct PyrolysisMaterial<'a> {
    /// Original timeline material identity.
    pub material_id: MaterialId,
    /// The exact original borrowed virgin caloric object.
    pub virgin: &'a dyn CondensedCaloric,
    /// Independent char calorics; required if char yield is positive.
    pub char_caloric: Option<&'a dyn CondensedCaloric>,
    /// Outgoing volatile enthalpy with a compatible formation reference.
    pub volatile: &'a dyn VolatileEnthalpy,
    /// Condensed residue kilograms per kilogram of virgin material consumed.
    pub char_mass_yield: f64,
    /// First-order finite-rate kinetics.
    pub kinetics: Arrhenius,
    /// Inclusive common temperature domain of every supplied chart, K.
    pub temperature_range_k: [f64; 2],
    /// Source/version and parameter applicability, supplied by the caller.
    pub source: &'a str,
}

/// Whole-advance time and algebraic budgets. They do not certify model accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PyrolysisControls {
    /// Requested maximum explicit conduction/reaction interval, s.
    pub maximum_dt_s: f64,
    /// Bounded scalar energy evaluations per cell and interval, 4 through 256.
    pub maximum_temperature_iterations: usize,
    /// Absolute local energy closure and kinetic energy-roundoff budget, J.
    pub cell_energy_tolerance_j: f64,
    /// Whole-advance absolute energy closure plus kinetic roundoff budget, J.
    pub total_energy_tolerance_j: f64,
    /// Combined mass closure plus accumulated absolute kinetic roundoff, kg.
    pub mass_tolerance_kg: f64,
}
impl Default for PyrolysisControls {
    fn default() -> Self {
        Self {
            maximum_dt_s: 0.02,
            maximum_temperature_iterations: 96,
            cell_energy_tolerance_j: 1e-8,
            total_energy_tolerance_j: 1e-6,
            mass_tolerance_kg: 1e-12,
        }
    }
}

/// Actual retained reactive inventory and common solid temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PyrolysisCell {
    /// Unconverted virgin material, kg.
    pub virgin_mass_kg: f64,
    /// Retained char, kg; it is not gas and contributes no ideal-gas pressure.
    pub char_mass_kg: f64,
    /// Common thermal-equilibrium temperature of these two condensed components.
    pub temperature_k: f64,
}

/// Time-integrated outgoing material from one original cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct VolatileRelease {
    /// Volatile mass emitted during this advance, kg.
    pub mass_kg: f64,
    /// Its actual carried enthalpy, J. This already contains reaction enthalpy.
    pub enthalpy_j: f64,
    /// Virgin mass consumed; char retention accounts for the difference.
    pub consumed_virgin_mass_kg: f64,
}

/// Coupled material/energy account from the published endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct PyrolysisAudit {
    /// Original thermal storage and external-source account. Its residual also
    /// includes outgoing volatile enthalpy, as required for this open solid.
    pub thermal: EnergyAudit,
    /// Source-cell-resolved volatile transfers, suitable for a gas coupling.
    pub releases: BTreeMap<CellId, VolatileRelease>,
    /// Total outgoing volatile mass, kg.
    pub volatile_mass_kg: f64,
    /// Total outgoing volatile enthalpy, J; do not add heat of pyrolysis again.
    pub volatile_enthalpy_j: f64,
    /// Actual final solid minus initial solid plus emitted gas mass, kg.
    pub mass_residual_kg: f64,
    /// Largest actual local stored-plus-emitted-minus-input energy defect, J.
    pub maximum_cell_energy_residual_j: f64,
    /// Sum of absolute local defects over all cells and accepted substeps, J.
    /// Implicit/remapping extensions also include their separately reported
    /// heat-transfer summation and remap discrepancies in this budget.
    pub absolute_energy_residual_j: f64,
    /// Sum of absolute consumed-mass roundoff against the stably computed BE
    /// extent. It is numerical discrepancy, not an extra physical mass stream.
    pub kinetic_mass_roundoff_kg: f64,
    /// Same discrepancy times local reaction enthalpy, summed in absolute J.
    /// No energy is injected or emitted for this account.
    pub kinetic_energy_roundoff_j: f64,
    /// Actual bounded scalar energy evaluations.
    pub temperature_evaluations: usize,
}

/// Fixed-geometry reactive solid on the original cells and transport operator.
#[derive(Clone)]
pub struct PyrolysisSimulation<'a> {
    thermal: DepositionSimulation<'a>,
    materials: BTreeMap<MaterialId, PyrolysisMaterial<'a>>,
    inventories: BTreeMap<CellId, PyrolysisCell>,
}

impl PyrolysisMaterial<'_> {
    fn enthalpies(&self, t: f64) -> Result<[f64; 3], DepositionError> {
        if !self.temperature_range_k[0].le(&t) || !t.le(&self.temperature_range_k[1]) {
            return Err(DepositionError::MaterialDomain);
        }
        let virgin = original_h(self.virgin, t)?;
        let char_h = if let Some(law) = self.char_caloric {
            original_h(law, t)?
        } else {
            virgin
        };
        let gas = finite(self.volatile.specific_enthalpy_j_kg(t)?)?;
        // Differences of the same datum cancel before yield weighting.
        let reaction = finite(
            self.char_mass_yield * (char_h - virgin)
                + (1.0 - self.char_mass_yield) * (gas - virgin),
        )?;
        if reaction < 0.0 {
            return Err(DepositionError::InvalidInput(
                "this pyrolysis lane requires endothermic reaction enthalpy throughout its domain",
            ));
        }
        Ok([virgin, char_h, gas])
    }

    fn admit(&self) -> Result<(), DepositionError> {
        let [lo, hi] = self.temperature_range_k;
        if !positive(lo)
            || !positive(hi)
            || lo >= hi
            || self.source.trim().is_empty()
            || !nonnegative(self.char_mass_yield)
            || self.char_mass_yield >= 1.0
            || (self.char_mass_yield > 0.0 && self.char_caloric.is_none())
            || !positive(self.virgin.maximum_temperature_slope())
            || self
                .char_caloric
                .is_some_and(|x| !positive(x.maximum_temperature_slope()))
        {
            return Err(DepositionError::InvalidInput(
                "invalid pyrolysis source, domain, residue yield or caloric slope",
            ));
        }
        self.kinetics.rate_constant_s_inv(lo)?;
        let a = self.enthalpies(lo)?;
        let b = self.enthalpies(lo + 0.5 * (hi - lo))?;
        let c = self.enthalpies(hi)?;
        if !(a[0] < b[0] && b[0] < c[0])
            || (self.char_mass_yield > 0.0 && !(a[1] < b[1] && b[1] < c[1]))
            || !(a[2] <= b[2] && b[2] <= c[2])
        {
            return Err(DepositionError::InvalidInput(
                "non-monotone pyrolysis caloric chart",
            ));
        }
        Ok(())
    }

    fn mixture_slope(&self, inventory: PyrolysisCell) -> Result<f64, DepositionError> {
        let mut capacity = inventory.virgin_mass_kg / self.virgin.maximum_temperature_slope();
        if inventory.char_mass_kg > 0.0 {
            capacity += inventory.char_mass_kg
                / self
                    .char_caloric
                    .ok_or(DepositionError::MaterialDomain)?
                    .maximum_temperature_slope();
        }
        let slope = (inventory.virgin_mass_kg + inventory.char_mass_kg) / capacity;
        if !positive(slope) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(slope)
    }
}

fn original_h(law: &dyn CondensedCaloric, t: f64) -> Result<f64, DepositionError> {
    let h = finite(law.specific_enthalpy_j_kg(t)?)?;
    let inverse = law.temperature_k(h)?;
    if !positive(inverse) || (inverse - t).abs() > 128.0 * f64::EPSILON * t.max(inverse) {
        return Err(DepositionError::InvalidInput(
            "pyrolysis h(T) differs from original caloric law",
        ));
    }
    Ok(h)
}

impl<'a> PyrolysisSimulation<'a> {
    /// Adopt an activated all-virgin thermal state without copying its physics.
    ///
    /// Future birth/contact events are intentionally refused: mass-changing
    /// solids need explicit geometry ownership. Temperatures, identities, cell
    /// energy references and original borrowed material laws must agree.
    pub fn new(
        thermal: DepositionSimulation<'a>,
        materials: &[PyrolysisMaterial<'a>],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if thermal.cells.is_empty()
            || !thermal.pending_events().is_empty()
            || materials.len() > thermal.budget.max_cells
        {
            return Err(DepositionError::InvalidInput(
                "pyrolysis requires activated cells, bounded materials and no pending geometry events",
            ));
        }
        let mut map = BTreeMap::new();
        for &material in materials {
            checkpoint()?;
            material.admit()?;
            if !thermal.uses_material_law(material.material_id, material.virgin)
                || map.insert(material.material_id, material).is_some()
            {
                return Err(DepositionError::InvalidInput(
                    "duplicate pyrolysis material or original virgin-law mismatch",
                ));
            }
        }
        let mut inventories = BTreeMap::new();
        for (&id, cell) in &thermal.cells {
            checkpoint()?;
            let model = map
                .get(&cell.material_id)
                .ok_or(DepositionError::InvalidInput(
                    "missing active pyrolysis material",
                ))?;
            let t = thermal.temperature_k(id)?;
            let h = model.enthalpies(t)?[0];
            if (h - cell.specific_h_j_kg).abs()
                > 128.0 * f64::EPSILON * h.abs().max(cell.specific_h_j_kg.abs()).max(1.0)
            {
                return Err(DepositionError::InvalidInput(
                    "original virgin enthalpy reference mismatch",
                ));
            }
            inventories.insert(
                id,
                PyrolysisCell {
                    virgin_mass_kg: cell.mass_kg,
                    char_mass_kg: 0.0,
                    temperature_k: t,
                },
            );
        }
        checkpoint()?;
        Ok(Self {
            thermal,
            materials: map,
            inventories,
        })
    }

    /// Original accepted clock, s.
    pub fn time_s(&self) -> f64 {
        self.thermal.time_s
    }
    /// Original cell identity, geometry, boundary and actual mixture energy.
    pub fn cells(&self) -> &BTreeMap<CellId, DepositedCell> {
        &self.thermal.cells
    }
    /// Actual condensed inventory; no volatile material remains in these cells.
    pub fn inventories(&self) -> &BTreeMap<CellId, PyrolysisCell> {
        &self.inventories
    }
    /// Sum of actual remaining condensed masses, kg.
    pub fn mass_kg(&self) -> f64 {
        self.thermal.mass_kg()
    }
    /// Sum of actual original cell `mass * mixture_h`, J.
    pub fn energy_j(&self) -> f64 {
        self.thermal.energy_j()
    }
    /// Composition-dependent common temperature, K, not a virgin-only Cp proxy.
    pub fn temperature_k(&self, id: CellId) -> Result<f64, DepositionError> {
        self.inventories
            .get(&id)
            .map(|c| c.temperature_k)
            .ok_or(DepositionError::InvalidInput("inactive pyrolysis cell"))
    }

    /// Atomically change one authored thermal boundary using original admission.
    /// The caller can supply convection or a separately resolved radiation flux.
    /// This endpoint change advances no time and does not apply its power twice.
    pub fn set_boundary(
        &mut self,
        id: CellId,
        boundary: Boundary,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        let mut trial = self.thermal.clone();
        trial.apply(Event::Boundary(id, boundary), &mut empty_thermal_audit())?;
        checkpoint()?;
        self.thermal = trial;
        Ok(())
    }

    /// Atomically replace supplied link conductances through the original owner.
    /// Material/porosity-dependent conductivity must be explicitly resolved by
    /// the caller; changing composition never silently invents a mixing rule.
    pub fn set_conductances(
        &mut self,
        conductances_w_k: &BTreeMap<u64, f64>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        let mut trial = self.thermal.clone();
        trial.set_step_conductances(conductances_w_k, &mut checkpoint)?;
        checkpoint()?;
        self.thermal = trial;
        Ok(())
    }

    /// Advance actual heat transport and thermally coupled finite-rate reaction.
    ///
    /// Conduction uses the original explicit stability bound with the CURRENT
    /// condensed mass and mixture heat-capacity bound. Each local reaction is
    /// implicit in its final temperature and extent. This is first-order time
    /// integration, with no unconditional stability claim. Near complete loss
    /// a fixed conductive mesh can become stiff; the original whole-advance
    /// substep budget refuses rather than hiding this with a seed mass.
    pub fn advance_to(
        &self,
        target_s: f64,
        controls: PyrolysisControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, PyrolysisAudit), DepositionError> {
        checkpoint()?;
        if !target_s.is_finite()
            || target_s < self.time_s()
            || !positive(controls.maximum_dt_s)
            || !(4..=256).contains(&controls.maximum_temperature_iterations)
            || !nonnegative(controls.cell_energy_tolerance_j)
            || !nonnegative(controls.total_energy_tolerance_j)
            || !nonnegative(controls.mass_tolerance_kg)
        {
            return Err(DepositionError::InvalidInput(
                "invalid pyrolysis time or numerical controls",
            ));
        }
        let mut next = self.clone();
        let mut audit = PyrolysisAudit {
            thermal: empty_thermal_audit(),
            releases: self
                .thermal
                .cells
                .keys()
                .map(|&id| (id, VolatileRelease::default()))
                .collect(),
            volatile_mass_kg: 0.0,
            volatile_enthalpy_j: 0.0,
            mass_residual_kg: 0.0,
            maximum_cell_energy_residual_j: 0.0,
            absolute_energy_residual_j: 0.0,
            kinetic_mass_roundoff_kg: 0.0,
            kinetic_energy_roundoff_j: 0.0,
            temperature_evaluations: 0,
        };
        let mut energy = Sum::default();
        let mut absolute_energy = Sum::default();
        let mut stored = Sum::default();
        let mut mass = Sum::default();
        let mut external = Sum::default();
        let mut kinetic_mass = Sum::default();
        let mut kinetic_energy = Sum::default();
        while next.time_s() < target_s {
            checkpoint()?;
            if audit.thermal.substeps >= self.thermal.budget.max_substeps {
                return Err(DepositionError::Budget(
                    "pyrolysis thermal substeps; refine geometry or shorten advance",
                ));
            }
            let remaining_s = target_s - next.time_s();
            let requested = controls.maximum_dt_s.min(remaining_s);
            let rates = next.thermal.heat_rates_with_calorics(
                requested,
                |id, cell| {
                    let inventory = next.inventories[&id];
                    let model = next.materials[&cell.material_id];
                    Ok((inventory.temperature_k, model.mixture_slope(inventory)?))
                },
                &mut checkpoint,
            )?;
            let dt = rates.explicit_dt_s;
            let endpoint = if dt == remaining_s {
                target_s
            } else {
                next.time_s() + dt
            };
            if !endpoint.is_finite() || endpoint <= next.time_s() {
                return Err(DepositionError::Arithmetic);
            }
            for (&id, cell) in &mut next.thermal.cells {
                checkpoint()?;
                let previous = next.inventories[&id];
                let input_j = finite(dt * rates.watts[&id])?;
                let step = solve_cell(
                    previous,
                    next.materials[&cell.material_id],
                    dt,
                    input_j,
                    controls,
                    &mut checkpoint,
                )?;
                let old_h = cell.specific_h_j_kg;
                let old_mass = cell.mass_kg;
                let next_mass = step.inventory.virgin_mass_kg + step.inventory.char_mass_kg;
                let actual_relative_storage = finite(next_mass * (step.specific_h_j_kg - old_h))?;
                let actual_mass_loss = old_mass - next_mass;
                // Cell-relative coordinates cancel a common formation datum
                // before arithmetic; reported gross energies retain that datum.
                let residual = finite(
                    actual_relative_storage
                        + step.release.mass_kg * (step.volatile_h_j_kg - old_h)
                        + (step.release.mass_kg - actual_mass_loss) * old_h
                        - input_j,
                )?;
                if residual.abs() > controls.cell_energy_tolerance_j {
                    return Err(DepositionError::EnergyBalance(residual));
                }
                if residual.abs() + step.kinetic_energy_roundoff_j
                    > controls.cell_energy_tolerance_j
                {
                    return Err(DepositionError::Budget(
                        "pyrolysis local reaction-energy roundoff",
                    ));
                }
                let local_mass = finite(step.release.mass_kg - actual_mass_loss)?;
                mass.add(local_mass)?;
                energy.add(residual)?;
                absolute_energy.add(residual.abs())?;
                kinetic_mass.add(step.kinetic_mass_roundoff_kg)?;
                kinetic_energy.add(step.kinetic_energy_roundoff_j)?;
                stored.add(finite(actual_relative_storage - actual_mass_loss * old_h)?)?;
                audit.maximum_cell_energy_residual_j =
                    audit.maximum_cell_energy_residual_j.max(residual.abs());
                audit.temperature_evaluations += step.evaluations;
                let release = audit
                    .releases
                    .get_mut(&id)
                    .expect("original cell set retained");
                release.mass_kg = finite(release.mass_kg + step.release.mass_kg)?;
                release.enthalpy_j = finite(release.enthalpy_j + step.release.enthalpy_j)?;
                release.consumed_virgin_mass_kg =
                    finite(release.consumed_virgin_mass_kg + step.release.consumed_virgin_mass_kg)?;
                cell.mass_kg = next_mass;
                cell.specific_h_j_kg = step.specific_h_j_kg;
                next.inventories.insert(id, step.inventory);
            }
            external.add(finite(dt * rates.external_w)?)?;
            next.thermal.time_s = endpoint;
            audit.thermal.substeps += 1;
        }
        audit.thermal.external_heat_j = external.value();
        audit.thermal.stored_energy_change_j = stored.value();
        audit.thermal.residual_j = energy.value();
        audit.absolute_energy_residual_j = absolute_energy.value();
        audit.mass_residual_kg = mass.value();
        audit.kinetic_mass_roundoff_kg = kinetic_mass.value();
        audit.kinetic_energy_roundoff_j = kinetic_energy.value();
        let mut emitted_mass = Sum::default();
        let mut emitted_enthalpy = Sum::default();
        for release in audit.releases.values() {
            emitted_mass.add(release.mass_kg)?;
            emitted_enthalpy.add(release.enthalpy_j)?;
        }
        audit.volatile_mass_kg = emitted_mass.value();
        audit.volatile_enthalpy_j = emitted_enthalpy.value();
        if audit.absolute_energy_residual_j > controls.total_energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(
                audit.absolute_energy_residual_j,
            ));
        }
        if audit.absolute_energy_residual_j + audit.kinetic_energy_roundoff_j
            > controls.total_energy_tolerance_j
        {
            return Err(DepositionError::Budget(
                "pyrolysis accumulated reaction-energy roundoff",
            ));
        }
        if audit.mass_residual_kg.abs() > controls.mass_tolerance_kg {
            return Err(DepositionError::InvalidInput(
                "pyrolysis solid/volatile mass balance",
            ));
        }
        if audit.mass_residual_kg.abs() + audit.kinetic_mass_roundoff_kg
            > controls.mass_tolerance_kg
        {
            return Err(DepositionError::Budget(
                "pyrolysis accumulated kinetic inventory roundoff",
            ));
        }
        checkpoint()?;
        Ok((next, audit))
    }
}

fn empty_thermal_audit() -> EnergyAudit {
    EnergyAudit {
        incoming_mass_kg: 0.0,
        incoming_energy_j: 0.0,
        external_heat_j: 0.0,
        stored_energy_change_j: 0.0,
        residual_j: 0.0,
        substeps: 0,
        events: 0,
    }
}

#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        let adjusted = finite(value - self.correction)?;
        let next = finite(self.sum + adjusted)?;
        self.correction = (next - self.sum) - adjusted;
        self.sum = next;
        Ok(())
    }
    fn value(&self) -> f64 {
        self.sum
    }
}

struct CellStep {
    inventory: PyrolysisCell,
    specific_h_j_kg: f64,
    volatile_h_j_kg: f64,
    release: VolatileRelease,
    evaluations: usize,
    kinetic_mass_roundoff_kg: f64,
    kinetic_energy_roundoff_j: f64,
}

fn solve_cell(
    initial: PyrolysisCell,
    model: PyrolysisMaterial<'_>,
    dt: f64,
    input_j: f64,
    controls: PyrolysisControls,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<CellStep, DepositionError> {
    let initial_h = model.enthalpies(initial.temperature_k)?;
    let evaluate = |t: f64| -> Result<(CellStep, f64), DepositionError> {
        let h = model.enthalpies(t)?;
        let kinetic = model
            .kinetics
            .isothermal_step(initial.virgin_mass_kg, t, dt)?;
        let consumed = kinetic.consumed_mass_kg;
        let char_mass = finite(initial.char_mass_kg + model.char_mass_yield * consumed)?;
        let out_mass = finite((1.0 - model.char_mass_yield) * consumed)?;
        let remaining_mass = finite(kinetic.virgin_mass_kg + char_mass)?;
        if !positive(remaining_mass) {
            return Err(DepositionError::InvalidInput(
                "pyrolysis cell exhausted; geometry transition required",
            ));
        }
        let h_mix = finite(h[0] + (char_mass / remaining_mass) * (h[1] - h[0]))?;
        let reaction_h =
            model.char_mass_yield * (h[1] - h[0]) + (1.0 - model.char_mass_yield) * (h[2] - h[0]);
        let residual = finite(
            initial.virgin_mass_kg * (h[0] - initial_h[0])
                + initial.char_mass_kg * (h[1] - initial_h[1])
                + consumed * reaction_h
                - input_j,
        )?;
        Ok((
            CellStep {
                inventory: PyrolysisCell {
                    virgin_mass_kg: kinetic.virgin_mass_kg,
                    char_mass_kg: char_mass,
                    temperature_k: t,
                },
                specific_h_j_kg: h_mix,
                volatile_h_j_kg: h[2],
                release: VolatileRelease {
                    mass_kg: out_mass,
                    enthalpy_j: finite(out_mass * h[2])?,
                    consumed_virgin_mass_kg: consumed,
                },
                evaluations: 0,
                kinetic_mass_roundoff_kg: kinetic.inventory_roundoff_kg.abs(),
                kinetic_energy_roundoff_j: finite(
                    kinetic.inventory_roundoff_kg.abs() * reaction_h,
                )?,
            },
            residual,
        ))
    };
    let [mut lo, mut hi] = model.temperature_range_k;
    checkpoint()?;
    let (left, mut f_lo) = evaluate(lo)?;
    checkpoint()?;
    let (right, mut f_hi) = evaluate(hi)?;
    if f_lo > 0.0 || f_hi < 0.0 {
        return Err(DepositionError::MaterialDomain);
    }
    // Local solve is tighter than final storage admission so representational
    // roundoff in actual mass*h has its own space within the caller's budget.
    let tolerance = 0.125 * controls.cell_energy_tolerance_j;
    for (mut step, residual) in [(left, f_lo), (right, f_hi)] {
        if residual.abs() <= tolerance {
            step.evaluations = 2;
            return Ok(step);
        }
    }
    for iteration in 2..controls.maximum_temperature_iterations {
        checkpoint()?;
        // Safeguarded secant accepts the middle 90%; periodic bisection makes
        // bracket progress independent of a sharply activated Arrhenius law.
        let fraction = -f_lo / (f_hi - f_lo);
        let t = if iteration % 4 != 0 && (0.05..=0.95).contains(&fraction) {
            lo + fraction * (hi - lo)
        } else {
            lo + 0.5 * (hi - lo)
        };
        if t <= lo || t >= hi {
            return Err(DepositionError::EnergyBalance(f_lo.abs().min(f_hi.abs())));
        }
        let (mut step, residual) = evaluate(t)?;
        if residual.abs() <= tolerance {
            step.evaluations = iteration + 1;
            return Ok(step);
        }
        if residual < 0.0 {
            lo = t;
            f_lo = residual;
        } else {
            hi = t;
            f_hi = residual;
        }
    }
    Err(DepositionError::Budget("pyrolysis temperature solve"))
}

#[cfg(test)]
#[path = "pyrolysis/tests.rs"]
mod tests;
