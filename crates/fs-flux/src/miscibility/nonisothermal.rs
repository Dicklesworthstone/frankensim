//! Closed, locally nonisothermal binary-liquid heat and species transport.
//!
//! Each finite volume retains its own temperature and both component amounts.
//! Equal fixed partial molar volumes make counterdiffusion volume preserving;
//! the same original homogeneous calorics supply energy, entropy and temperature
//! recovery. This lane requires kappa=0, positive common species support and a
//! convex homogeneous law throughout the declared temperature interval. It does
//! not assign a local temperature to unresolved square-gradient interfacial energy.
//!
//! A canonical face sweep splits the positive Onsager block into its chemical
//! rank-one part and Fourier part. Both use bounded implicit scalar solves.
//! Original source coefficients are frozen only during their pair solve and
//! refreshed at subsequent pairs/stages. Every debit has the same opposite
//! credit; whole-window species, energy and entropy checks precede publication.
//! The sweep is first order. The adaptive owner compares one sweep with two
//! half sweeps using both composition and temperature.

use super::{
    BinaryMiscibilityError, BinaryMiscibilityNetwork, BinaryMiscibilityState, Sum, finite,
    nonnegative, positive,
};
use fs_qty::Temperature;
use fs_thermochem::solution::caloric::{
    BinaryCaloricErrorV1, BinaryCaloricStateV1, BinaryLiquidCaloricModelV1,
    transport::{
        BinaryLiquidTransportCoefficientsV1, BinaryLiquidTransportPointV1,
        BinaryLiquidTransportSolveV1,
    },
};
use std::sync::Arc;

pub mod adaptive;
mod solve;

/// Explicit face data in the original network's canonical face order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryNonisothermalFace {
    /// Conductivity at zero diffusive species flux [W/(m K)].
    pub thermal_conductivity_w_m_k: f64,
    /// Ordered A-minus-B molar heat of transport [J/mol]. This authored Soret
    /// coefficient also enters reciprocal Dufour heat. Zero excludes that
    /// cross effect, but retains source-derived partial-enthalpy carriage.
    pub heat_of_transport_j_per_mol: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BinaryNonisothermalError {
    Invalid(&'static str),
    Spatial(BinaryMiscibilityError),
    Caloric(BinaryCaloricErrorV1),
    Workspace {
        required: usize,
        maximum: usize,
    },
    PairIterationLimit,
    PairOutsideTemperatureDomain,
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
    Entropy {
        change_j_per_k: f64,
        tolerance_j_per_k: f64,
    },
    Species {
        residual_mol: [f64; 2],
        tolerance_mol: f64,
    },
    Cancelled,
}
impl core::fmt::Display for BinaryNonisothermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonisothermal binary diffusion: {self:?}")
    }
}
impl std::error::Error for BinaryNonisothermalError {}
impl From<BinaryMiscibilityError> for BinaryNonisothermalError {
    fn from(value: BinaryMiscibilityError) -> Self {
        if value == BinaryMiscibilityError::Cancelled {
            Self::Cancelled
        } else {
            Self::Spatial(value)
        }
    }
}
impl From<BinaryCaloricErrorV1> for BinaryNonisothermalError {
    fn from(value: BinaryCaloricErrorV1) -> Self {
        if value == BinaryCaloricErrorV1::Cancelled {
            Self::Cancelled
        } else {
            Self::Caloric(value)
        }
    }
}
type Result<T> = core::result::Result<T, BinaryNonisothermalError>;

/// Numerical accuracy budgets, applied to the complete accepted interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryNonisothermalSolve {
    pub maximum_root_iterations: usize,
    pub absolute_transfer_tolerance_mol: f64,
    pub relative_transfer_tolerance: f64,
    pub temperature_tolerance_k: f64,
    pub amount_balance_tolerance_mol: f64,
    pub energy_tolerance_j: f64,
    pub entropy_tolerance_j_per_k: f64,
}
impl Default for BinaryNonisothermalSolve {
    fn default() -> Self {
        Self {
            maximum_root_iterations: 96,
            absolute_transfer_tolerance_mol: 1e-13,
            relative_transfer_tolerance: 1e-11,
            temperature_tolerance_k: 1e-9,
            amount_balance_tolerance_mol: 1e-12,
            energy_tolerance_j: 1e-8,
            entropy_tolerance_j_per_k: 1e-10,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Owner {
    network: BinaryMiscibilityNetwork,
    caloric: BinaryLiquidCaloricModelV1,
    faces: Arc<[BinaryNonisothermalFace]>,
    temperature_bracket: [Temperature; 2],
    workspace_bytes: usize,
}

/// Immutable accepted composition, local temperatures and original sources.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalField {
    owner: Arc<Owner>,
    state: BinaryMiscibilityState,
    temperatures: Vec<Temperature>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalStep {
    pub duration_s: f64,
    /// Signed A transfer from each canonical left cell to right cell [mol].
    /// B has the exact opposite transfer.
    pub transferred_a_mol: Vec<f64>,
    /// Total paired enthalpy/internal-energy transfer in the same orientation.
    pub transferred_energy_j: Vec<f64>,
    pub root_iterations: u64,
    pub caloric_iterations: u64,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalCandidate {
    field: BinaryNonisothermalField,
    step: BinaryNonisothermalStep,
}
impl BinaryNonisothermalCandidate {
    pub const fn field(&self) -> &BinaryNonisothermalField {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalStep {
        &self.step
    }
    pub fn into_parts(self) -> (BinaryNonisothermalField, BinaryNonisothermalStep) {
        (self.field, self.step)
    }
}

impl BinaryMiscibilityNetwork {
    /// Retain the accepted inventories and clock while assigning each cell's
    /// actual temperature. The caller explicitly supplies original liquid
    /// calorics and qualified conductivity/heat-of-transport face data. A bound
    /// material diffusivity continues to use its original claim and context.
    ///
    /// This is a thermal initialization, not an energy-conserving temperature
    /// event. Subsequent closed evolution conserves total energy. Original
    /// material-card storage/query receipts are outside the numerical workspace.
    pub fn bind_nonisothermal_field(
        &self,
        initial: &BinaryMiscibilityState,
        caloric: BinaryLiquidCaloricModelV1,
        temperatures: Vec<Temperature>,
        faces: Vec<BinaryNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryNonisothermalField> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        if self.model().gradient_coefficient_j_m() != 0.0
            || self.model().thermodynamics() != caloric.thermodynamics()
            || temperatures.len() != self.cells().len()
            || faces.len() != self.owner.faces.len()
            || faces.iter().any(|f| {
                !nonnegative(f.thermal_conductivity_w_m_k)
                    || !f.heat_of_transport_j_per_mol.is_finite()
            })
        {
            return Err(BinaryNonisothermalError::Invalid(
                "same original convex binary calorics, kappa=0, one temperature per cell and valid canonical transport faces required",
            ));
        }
        caloric.admit_transport_temperature_interval(temperature_bracket)?;
        let workspace_bytes = self
            .cells()
            .len()
            .checked_mul(2048)
            .and_then(|n| {
                self.owner
                    .faces
                    .len()
                    .checked_mul(512)
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| n.checked_add(self.workspace_bytes()))
            .and_then(|n| n.checked_add(core::mem::size_of::<Owner>()))
            .ok_or(BinaryNonisothermalError::Invalid(
                "nonisothermal workspace size overflow",
            ))?;
        if workspace_bytes > maximum_workspace_bytes {
            return Err(BinaryNonisothermalError::Workspace {
                required: workspace_bytes,
                maximum: maximum_workspace_bytes,
            });
        }
        let field = BinaryNonisothermalField {
            owner: Arc::new(Owner {
                network: self.clone(),
                caloric,
                faces: faces.into(),
                temperature_bracket,
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

impl BinaryNonisothermalField {
    /// Original geometry and homogeneous pair; local temperatures belong to
    /// this field and must not be replaced by the network's reference T.
    pub fn network(&self) -> &BinaryMiscibilityNetwork {
        &self.owner.network
    }
    pub fn caloric_model(&self) -> &BinaryLiquidCaloricModelV1 {
        &self.owner.caloric
    }
    pub fn component_moles(&self) -> &[[f64; 2]] {
        &self.state.amounts
    }
    pub fn mole_fractions(&self) -> &[[f64; 2]] {
        &self.state.fractions
    }
    pub fn temperatures(&self) -> &[Temperature] {
        &self.temperatures
    }
    pub fn elapsed_s(&self) -> f64 {
        self.state.elapsed_s
    }
    pub fn workspace_bytes(&self) -> usize {
        self.owner.workspace_bytes
    }
    pub fn transport_faces(&self) -> &[BinaryNonisothermalFace] {
        &self.owner.faces
    }
    pub fn temperature_bracket(&self) -> [Temperature; 2] {
        self.owner.temperature_bracket
    }

    /// One implicit canonical face sweep, with all physical state held until
    /// every local root and whole-field conservation/entropy check succeeds.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: BinaryNonisothermalSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryNonisothermalCandidate> {
        solve::step(self, duration_s, controls, cancelled)
    }

    fn caloric_state(&self, cell: usize) -> Result<BinaryCaloricStateV1> {
        Ok(self
            .caloric_model()
            .homogeneous(self.temperatures[cell], self.state.amounts[cell])?)
    }
    fn point(&self, cell: usize) -> BinaryLiquidTransportPointV1 {
        BinaryLiquidTransportPointV1 {
            temperature: self.temperatures[cell],
            mole_fractions: self.state.fractions[cell],
        }
    }
    fn face_diffusivity(&self, face: usize, cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
        let physical = self.network().owner.faces[face].physical;
        #[cfg(any(feature = "continuum", feature = "experimental-miscibility-material"))]
        if let Some(source) = &self.network().owner.diffusivity_source {
            let [i, j] = physical.cells;
            return Ok(source
                .resolve_nonisothermal_face(
                    [self.temperatures[i], self.temperatures[j]],
                    [self.state.fractions[i], self.state.fractions[j]],
                    cancelled,
                )
                .map_err(super::diffusivity::source_error)?
                .diffusivity_m2_s);
        }
        let _ = cancelled;
        Ok(physical.diffusivity_m2_s)
    }
    fn qualify(&self, cancelled: &mut impl FnMut() -> bool) -> Result<()> {
        let [lo, hi] = self.owner.temperature_bracket.map(|t| t.value());
        for (cell, &t) in self.temperatures.iter().enumerate() {
            poll(cancelled)?;
            if !t.value().is_finite()
                || t.value() < lo
                || t.value() > hi
                || self.state.fractions[cell].iter().any(|&x| !positive(x))
            {
                return Err(BinaryNonisothermalError::Invalid(
                    "source-supported temperatures and positive common species support required",
                ));
            }
            self.caloric_state(cell)?;
            #[cfg(any(feature = "continuum", feature = "experimental-miscibility-material"))]
            if let Some(source) = &self.network().owner.diffusivity_source {
                if self.network().owner.adjacent_offsets[cell]
                    == self.network().owner.adjacent_offsets[cell + 1]
                {
                    source
                        .resolve_face(t, [self.state.fractions[cell]; 2], &mut *cancelled)
                        .map_err(super::diffusivity::source_error)?;
                }
            }
        }
        for face in 0..self.owner.faces.len() {
            poll(cancelled)?;
            self.face_diffusivity(face, cancelled)?;
        }
        Ok(())
    }
    fn balances(
        &self,
        next: &Self,
        controls: BinaryNonisothermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; 2], f64, f64, f64)> {
        let mut amounts = [Sum::default(); 2];
        let mut energy = Sum::default();
        let mut entropy = Sum::default();
        let mut capacity = 0.0_f64;
        for cell in 0..self.temperatures.len() {
            poll(cancelled)?;
            let mut residual = Sum::default();
            for s in 0..2 {
                amounts[s].add(next.state.amounts[cell][s] - self.state.amounts[cell][s])?;
                amounts[s].add(
                    next.state.amount_remainders[cell][s] - self.state.amount_remainders[cell][s],
                )?;
                residual.add(next.state.amounts[cell][s])?;
                residual.add(next.state.amount_remainders[cell][s])?;
            }
            residual.add(-self.network().owner.capacities_mol[cell])?;
            capacity = capacity.max(residual.value()?.abs());
            let change = self.caloric_model().homogeneous_transport_change(
                &self.caloric_state(cell)?,
                &next.caloric_state(cell)?,
            )?;
            energy.add(change.enthalpy_change_j)?;
            entropy.add(change.entropy_change_j_per_k)?;
        }
        let species = [amounts[0].value()?, amounts[1].value()?];
        if capacity > controls.amount_balance_tolerance_mol
            || species
                .iter()
                .any(|r| r.abs() > controls.amount_balance_tolerance_mol)
        {
            return Err(BinaryNonisothermalError::Species {
                residual_mol: species,
                tolerance_mol: controls.amount_balance_tolerance_mol,
            });
        }
        let energy = energy.value()?;
        if energy.abs() > controls.energy_tolerance_j {
            return Err(BinaryNonisothermalError::Energy {
                residual_j: energy,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        let entropy = entropy.value()?;
        if entropy < -controls.entropy_tolerance_j_per_k {
            return Err(BinaryNonisothermalError::Entropy {
                change_j_per_k: entropy,
                tolerance_j_per_k: controls.entropy_tolerance_j_per_k,
            });
        }
        Ok((species, capacity, energy, entropy))
    }
}

fn admit_controls(c: BinaryNonisothermalSolve) -> Result<()> {
    if !(1..=256).contains(&c.maximum_root_iterations)
        || !positive(c.temperature_tolerance_k)
        || [
            c.absolute_transfer_tolerance_mol,
            c.relative_transfer_tolerance,
            c.amount_balance_tolerance_mol,
            c.energy_tolerance_j,
            c.entropy_tolerance_j_per_k,
        ]
        .into_iter()
        .any(|x| !nonnegative(x))
        || c.absolute_transfer_tolerance_mol + c.relative_transfer_tolerance == 0.0
    {
        return Err(BinaryNonisothermalError::Invalid(
            "bounded roots and finite nonnegative physical error budgets required",
        ));
    }
    Ok(())
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(BinaryNonisothermalError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
