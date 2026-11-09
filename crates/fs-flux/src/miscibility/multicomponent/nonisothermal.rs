//! Closed local-temperature liquids with fully coupled Maxwell–Stefan mixing.
//!
//! Every species inventory and every cell temperature remain independent. The
//! original A+B/T liquid potential and pure Cp curves supply both the face
//! enthalpy carriage and the local caloric inverse. Thermal force uses the
//! source-conjugate composition effort: unequal temperatures alone cause no
//! diffusion unless an explicit heat-of-transport contrast is supplied. That
//! contrast also contributes reciprocal Dufour heat exactly once.
//!
//! One synchronous explicit step evaluates all faces at the accepted state,
//! applies paired species/energy transfers, and recovers local temperatures.
//! Only a complete candidate passing independent inventory, energy and entropy
//! checks can be published. Positive constitutive entropy production is not
//! substituted for the finite-step entropy check. Equal fixed partial molar
//! volumes, a fixed static pressure and positive common support are required.
//! The whole declared temperature bracket must satisfy a sufficient convexity
//! bound; nonconvex phase interfaces require their own energetic owner.

use super::{
    MulticomponentMiscibilityError, MulticomponentMiscibilityField,
    MulticomponentMiscibilityNetwork, MulticomponentMiscibilitySolve, Sum, finite, nonnegative,
    positive,
};
use fs_qty::Temperature;
use fs_thermochem::solution::caloric::BinaryCaloricErrorV1;
use fs_thermochem::solution::multicomponent::{
    MulticomponentSolutionErrorV1,
    caloric::{
        MulticomponentCaloricErrorV1, MulticomponentCaloricStateV1,
        MulticomponentLiquidCaloricModelV1,
        transport::{MulticomponentLiquidTransportPointV1, MulticomponentLiquidTransportSolveV1},
    },
};
use fs_thermochem::transport::maxwell_stefan::MaxwellStefanMolarSolveV1;
use std::sync::Arc;

mod precision;
pub(super) mod solve;

/// Coupled composition and local-temperature accuracy over a physical interval.
pub mod adaptive;

/// Conservative prescribed-volume stirring with original parcel enthalpy.
pub mod advection;
/// Stirring, full molecular diffusion and heat exchange over one shared clock.
pub mod transport;

/// Original material sources refreshed at every accepted local temperature.
#[cfg(any(feature = "continuum", feature = "experimental-miscibility-material"))]
pub mod material;

/// Independent kinetic data in the original network's canonical face order.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalFace {
    /// Conductivity at zero diffusive species flux, in W/(m K).
    pub thermal_conductivity_w_m_k: f64,
    /// Explicit canonical species heats of transport in J/mol. Only contrasts
    /// matter in the zero-net-molar frame. An all-zero vector excludes Soret
    /// and Dufour coupling while retaining original partial-enthalpy carriage.
    pub heat_of_transport_j_per_mol: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentNonisothermalError {
    Invalid(&'static str),
    Spatial(MulticomponentMiscibilityError),
    Caloric(MulticomponentCaloricErrorV1),
    Workspace {
        required: usize,
        maximum: usize,
    },
    /// Carried amount low words contain more original caloric energy than the
    /// declared budget permits the representable thermodynamic state to omit.
    UnresolvedCaloricInventory {
        cell: usize,
        bound_j: f64,
        tolerance_j: f64,
    },
    /// The original parcel enthalpy cannot be transported at the requested
    /// numerical energy resolution without amplifying rounding into heat.
    UnresolvedParcelEnthalpy {
        bound_j: f64,
        tolerance_j: f64,
    },
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
    Entropy {
        change_j_per_k: f64,
        tolerance_j_per_k: f64,
    },
    Cancelled,
}
impl core::fmt::Display for MulticomponentNonisothermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonisothermal multicomponent miscibility: {self:?}")
    }
}
impl std::error::Error for MulticomponentNonisothermalError {}
impl From<MulticomponentMiscibilityError> for MulticomponentNonisothermalError {
    fn from(error: MulticomponentMiscibilityError) -> Self {
        if error == MulticomponentMiscibilityError::Cancelled {
            Self::Cancelled
        } else {
            Self::Spatial(error)
        }
    }
}
impl From<MulticomponentCaloricErrorV1> for MulticomponentNonisothermalError {
    fn from(error: MulticomponentCaloricErrorV1) -> Self {
        match error {
            MulticomponentCaloricErrorV1::Cancelled
            | MulticomponentCaloricErrorV1::Caloric(BinaryCaloricErrorV1::Cancelled)
            | MulticomponentCaloricErrorV1::Thermodynamics(
                MulticomponentSolutionErrorV1::Cancelled,
            ) => Self::Cancelled,
            other => Self::Caloric(other),
        }
    }
}
type Result<T> = core::result::Result<T, MulticomponentNonisothermalError>;

/// Complete-interval physical budgets and bounded original caloric inversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentNonisothermalSolve {
    pub maxwell_stefan: MaxwellStefanMolarSolveV1,
    pub maximum_caloric_iterations: usize,
    pub temperature_tolerance_k: f64,
    pub amount_balance_tolerance_mol: f64,
    pub energy_tolerance_j: f64,
    pub entropy_tolerance_j_per_k: f64,
}
impl Default for MulticomponentNonisothermalSolve {
    fn default() -> Self {
        Self {
            maxwell_stefan: MaxwellStefanMolarSolveV1::default(),
            maximum_caloric_iterations: 96,
            temperature_tolerance_k: 1e-9,
            amount_balance_tolerance_mol: 1e-12,
            energy_tolerance_j: 1e-8,
            entropy_tolerance_j_per_k: 1e-10,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Owner {
    caloric: MulticomponentLiquidCaloricModelV1,
    faces: Box<[MulticomponentNonisothermalFace]>,
    temperature_bracket: [Temperature; 2],
    convexity_lower_bound: f64,
    maximum_inverse_temperature_magnitude_k: f64,
    workspace_bytes: usize,
}

/// Accepted local temperatures and independent original species inventories.
/// The isothermal state is deliberately private: its reference-temperature
/// free-energy gate is not the admissibility condition for this thermal owner.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalField {
    owner: Arc<Owner>,
    state: MulticomponentMiscibilityField,
    temperatures: Vec<Temperature>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalStep {
    pub duration_s: f64,
    /// Canonical face then canonical species, positive left to right.
    pub transferred_mol: Vec<Vec<f64>>,
    /// Original enthalpy transfer, including Fourier and reciprocal Dufour heat.
    pub transferred_energy_j: Vec<f64>,
    pub maxwell_stefan_evaluations: usize,
    pub maximum_gradient_residual_per_m: f64,
    pub maximum_molar_frame_residual_mol_m2_s: f64,
    pub caloric_iterations: u64,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
    /// Initial constitutive rate; the actual finite-step entropy is audited
    /// separately above. This is not a time-integrated dissipation bound.
    pub entropy_production_rate_w_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalCandidate {
    field: MulticomponentNonisothermalField,
    step: MulticomponentNonisothermalStep,
}
impl MulticomponentNonisothermalCandidate {
    pub const fn field(&self) -> &MulticomponentNonisothermalField {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentNonisothermalStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentNonisothermalField,
        MulticomponentNonisothermalStep,
    ) {
        (self.field, self.step)
    }
}

impl MulticomponentMiscibilityNetwork {
    /// Initialize local temperatures without losing accepted inventories, low
    /// words or elapsed time. This assignment is a thermal initialization, not
    /// a heat-transfer event. Subsequent evolution audits the original energy.
    /// Face coefficient rows follow `self.faces()` after canonical sorting.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_nonisothermal_field(
        &self,
        initial: &MulticomponentMiscibilityField,
        caloric: MulticomponentLiquidCaloricModelV1,
        temperatures: Vec<Temperature>,
        faces: Vec<MulticomponentNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentNonisothermalField> {
        self.prepare_nonisothermal_field(
            initial,
            caloric,
            temperatures,
            faces,
            temperature_bracket,
            maximum_workspace_bytes,
            true,
            cancelled,
        )
    }

    /// Shared original-source initialization. Only the explicit positive
    /// interface-energy owner may retain a nonconvex bulk branch; that owner
    /// never exposes this raw homogeneous thermal field.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare_nonisothermal_field(
        &self,
        initial: &MulticomponentMiscibilityField,
        caloric: MulticomponentLiquidCaloricModelV1,
        temperatures: Vec<Temperature>,
        faces: Vec<MulticomponentNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        require_convexity: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentNonisothermalField> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        if (!Arc::ptr_eq(&self.owner, &initial.network.owner) && self != initial.network())
            || self.model() != caloric.thermodynamics()
            || temperatures.len() != self.cells().len()
            || faces.len() != self.faces().len()
            || faces.iter().any(|f| {
                !nonnegative(f.thermal_conductivity_w_m_k)
                    || f.heat_of_transport_j_per_mol.len() != n
                    || f.heat_of_transport_j_per_mol.iter().any(|q| !q.is_finite())
            })
        {
            return Err(MulticomponentNonisothermalError::Invalid(
                "same original network/caloric law, one local temperature per cell and canonical finite thermal face data required",
            ));
        }
        // The original closure first validates both source endpoints and
        // evaluates both curvature bounds. Only its explicit convexity refusal
        // may be retained by the interfacial owner; no source, arithmetic,
        // stability or cancellation error is suppressed.
        let convexity_lower_bound = match caloric
            .admit_transport_temperature_interval(temperature_bracket, &mut cancelled)
        {
            Ok(bound) => bound,
            Err(MulticomponentCaloricErrorV1::TransportConvexity { lower_bound })
                if !require_convexity =>
            {
                lower_bound
            }
            Err(error) => return Err(error.into()),
        };
        // Covers simultaneous coarse/fine fields, original caloric states,
        // face rates/receipts, shared friction workspace and vector headers.
        // Original material-card storage remains outside the numerical budget.
        let extra_cells = self
            .cells()
            .len()
            .checked_mul(n)
            .and_then(|v| v.checked_mul(512));
        let extra_faces = self
            .faces()
            .len()
            .checked_mul(n)
            .and_then(|v| v.checked_mul(512));
        let workspace_bytes = extra_cells
            .and_then(|v| extra_faces.and_then(|w| v.checked_add(w)))
            .and_then(|v| {
                self.cells()
                    .len()
                    .checked_mul(2048)
                    .and_then(|w| v.checked_add(w))
            })
            // Donor stirring retains one original caloric state per cell.
            // Its inline MAX-species storage does not shrink with actual N;
            // account for that payload and its independent amount allocation.
            .and_then(|v| {
                n.checked_mul(core::mem::size_of::<f64>())
                    .and_then(|w| w.checked_add(core::mem::size_of::<MulticomponentCaloricStateV1>()))
                    .and_then(|w| self.cells().len().checked_mul(w))
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| {
                self.faces()
                    .len()
                    .checked_mul(512)
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| v.checked_add(self.workspace_bytes()))
            .and_then(|v| v.checked_add(32768 + core::mem::size_of::<Owner>()))
            .ok_or(MulticomponentNonisothermalError::Invalid(
                "thermal workspace size overflow",
            ))?;
        if workspace_bytes > maximum_workspace_bytes {
            return Err(MulticomponentNonisothermalError::Workspace {
                required: workspace_bytes,
                maximum: maximum_workspace_bytes,
            });
        }
        let field = MulticomponentNonisothermalField {
            owner: Arc::new(Owner {
                maximum_inverse_temperature_magnitude_k: caloric
                    .thermodynamics()
                    .chi_inverse_temperature_k()
                    .iter()
                    .fold(0.0_f64, |maximum, &b| maximum.max(b.abs())),
                caloric,
                faces: faces.into_boxed_slice(),
                temperature_bracket,
                convexity_lower_bound,
                workspace_bytes,
            }),
            state: initial.clone(),
            temperatures,
        };
        field.qualify(&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(field)
    }
}

impl MulticomponentNonisothermalField {
    /// Original geometry and molecular pair data. Its reference temperature is
    /// not a replacement for this field's actual local temperatures.
    pub fn network(&self) -> &MulticomponentMiscibilityNetwork {
        self.state.network()
    }
    pub fn caloric_model(&self) -> &MulticomponentLiquidCaloricModelV1 {
        &self.owner.caloric
    }
    pub fn temperatures(&self) -> &[Temperature] {
        &self.temperatures
    }
    pub fn component_moles(&self) -> &[f64] {
        self.state.component_moles()
    }
    pub fn mole_fractions(&self) -> &[f64] {
        self.state.mole_fractions()
    }
    pub fn cell_component_moles(&self, cell: usize) -> Option<&[f64]> {
        self.state.cell_component_moles(cell)
    }
    pub fn cell_mole_fractions(&self, cell: usize) -> Option<&[f64]> {
        self.state.cell_mole_fractions(cell)
    }
    pub fn elapsed_s(&self) -> f64 {
        self.state.elapsed_s()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.owner.workspace_bytes
    }
    pub fn temperature_bracket(&self) -> [Temperature; 2] {
        self.owner.temperature_bracket
    }
    pub fn transport_faces(&self) -> &[MulticomponentNonisothermalFace] {
        &self.owner.faces
    }
    pub fn convexity_lower_bound(&self) -> f64 {
        self.owner.convexity_lower_bound
    }

    /// One full coupled explicit physical interval. No clipping, dependent
    /// species reconstruction or hidden physical substeps are introduced.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentNonisothermalSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentNonisothermalCandidate> {
        solve::step(self, duration_s, controls, cancelled)
    }

    fn point(&self, cell: usize) -> MulticomponentLiquidTransportPointV1<'_> {
        MulticomponentLiquidTransportPointV1 {
            temperature: self.temperatures[cell],
            mole_fractions: self
                .cell_mole_fractions(cell)
                .expect("qualified thermal cell"),
        }
    }
    fn caloric_state(
        &self,
        cell: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MulticomponentCaloricStateV1> {
        Ok(self.caloric_model().homogeneous(
            self.temperatures[cell],
            self.cell_component_moles(cell)
                .expect("qualified thermal cell"),
            cancelled,
        )?)
    }
    fn qualify(&self, cancelled: &mut impl FnMut() -> bool) -> Result<()> {
        let [lo, hi] = self.temperature_bracket().map(|t| t.value());
        for (cell, &temperature) in self.temperatures.iter().enumerate() {
            poll(cancelled)?;
            if !temperature.value().is_finite()
                || temperature.value() < lo
                || temperature.value() > hi
            {
                return Err(MulticomponentNonisothermalError::Invalid(
                    "local temperature outside declared original-source bracket",
                ));
            }
            self.caloric_state(cell, cancelled)?;
        }
        Ok(())
    }
    fn balances(
        &self,
        next: &Self,
        controls: MulticomponentNonisothermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64, f64, f64)> {
        let balances = self.caloric_balances(next, controls, cancelled)?;
        admit_energy_entropy(balances.2, balances.3, controls)?;
        Ok(balances)
    }

    /// Original inventory and caloric increments, including the carried-amount
    /// energy-resolution guard. The caller applies its complete storage law:
    /// bulk enthalpy here, bulk plus interface energy in the phase-field owner.
    pub(super) fn caloric_balances(
        &self,
        next: &Self,
        controls: MulticomponentNonisothermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64, f64, f64)> {
        let (species, capacity) = self.state.inventory_balances(
            &next.state,
            controls.amount_balance_tolerance_mol,
            cancelled,
        )?;
        let mut energy = Sum::default();
        let mut entropy = Sum::default();
        for cell in 0..self.temperatures.len() {
            poll(cancelled)?;
            let old_state = self.caloric_state(cell, cancelled)?;
            let new_state = next.caloric_state(cell, cancelled)?;
            precision::admit_inventory_resolution(
                self,
                next,
                cell,
                &old_state,
                &new_state,
                controls.energy_tolerance_j / self.temperatures.len() as f64,
                cancelled,
            )?;
            let change = self.caloric_model().homogeneous_transport_change(
                &old_state,
                &new_state,
                &mut *cancelled,
            )?;
            energy.add(change.enthalpy_change_j)?;
            entropy.add(change.entropy_change_j_per_k)?;
        }
        Ok((species, capacity, energy.value()?, entropy.value()?))
    }
}

/// Apply the closed-window laws to the caller's complete stored energy.
pub(super) fn admit_energy_entropy(
    energy: f64,
    entropy: f64,
    controls: MulticomponentNonisothermalSolve,
) -> Result<()> {
    if energy.abs() > controls.energy_tolerance_j {
        return Err(MulticomponentNonisothermalError::Energy {
            residual_j: energy,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    if entropy < -controls.entropy_tolerance_j_per_k {
        return Err(MulticomponentNonisothermalError::Entropy {
            change_j_per_k: entropy,
            tolerance_j_per_k: controls.entropy_tolerance_j_per_k,
        });
    }
    Ok(())
}

pub(super) fn admit_controls(c: MulticomponentNonisothermalSolve) -> Result<()> {
    super::admit_controls(MulticomponentMiscibilitySolve {
        maxwell_stefan: c.maxwell_stefan,
        amount_balance_tolerance_mol: c.amount_balance_tolerance_mol,
        free_energy_tolerance_j: 0.0,
    })?;
    if !(1..=256).contains(&c.maximum_caloric_iterations)
        || !positive(c.temperature_tolerance_k)
        || !nonnegative(c.energy_tolerance_j)
        || !nonnegative(c.entropy_tolerance_j_per_k)
    {
        return Err(MulticomponentNonisothermalError::Invalid(
            "bounded caloric inverse and finite nonnegative energy/entropy allowances required",
        ));
    }
    Ok(())
}
pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentNonisothermalError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
