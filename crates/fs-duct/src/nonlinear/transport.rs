//! Nonlinear compressible propagation with conservative molecular transport.
//!
//! The existing geometric image feeds one branch of the shared gas composition
//! owner. Its molecular/advection/molecular step evolves species, momentum and
//! total energy together. Viscous stress work and diffusing species enthalpy
//! travel in the energy flux; adding a separate viscous-heating source would
//! count the same energy twice. Geometry and thermochemistry are not resampled.
//!
//! This is longitudinal Navier--Stokes--Fourier/Maxwell--Stefan transport in a
//! constant-area duct. It does not resolve radial boundary layers or supply the
//! linear acoustic Zwikker--Kosten wall memory. It therefore must not be used
//! as evidence that realistic bore-wall attenuation or physical shock thickness
//! has converged. The full nonlinear DG entropy/physical-loss program remains
//! separate.

use super::{FrozenGas, NonlinearDuct, NonlinearDuctError, NonlinearDuctSample};
use crate::Duct;
use fs_gas::multicomponent::dissipation::TransportLaw;
use fs_gas::multicomponent::network::MixtureNetwork1d;
pub use fs_gas::multicomponent::network::dissipating::{
    DissipatingNetwork1d, DissipatingNetworkError, DissipatingNetworkReport, DissipatingStepControl,
};
use fs_gas::multicomponent::{Conserved, Error as MixtureError, MixtureEuler1d, MixtureGas};
use fs_qty::Time;

/// Source, geometry or complete coupled-step refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum DissipatingDuctError {
    /// The existing nonlinear geometric image refused.
    Duct(NonlinearDuctError),
    /// Original composition-dependent caloric or conservative owner refused.
    Mixture(MixtureError),
    /// Original molecular/advection composition refused without publication.
    Transport(DissipatingNetworkError),
    /// Source identity, initialization or observable is unsupported.
    Invalid(&'static str),
}
impl core::fmt::Display for DissipatingDuctError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear molecular duct: {self:?}")
    }
}
impl std::error::Error for DissipatingDuctError {}
impl From<NonlinearDuctError> for DissipatingDuctError {
    fn from(error: NonlinearDuctError) -> Self {
        Self::Duct(error)
    }
}
impl From<MixtureError> for DissipatingDuctError {
    fn from(error: MixtureError) -> Self {
        Self::Mixture(error)
    }
}
impl From<DissipatingNetworkError> for DissipatingDuctError {
    fn from(error: DissipatingNetworkError) -> Self {
        Self::Transport(error)
    }
}

/// One geometric bore advanced by the original conservative gas owners.
///
/// The admitted transport law is evaluated at the evolving face temperature,
/// pressure and composition. All candidate fields, source workspaces and
/// reports are prepared once by [`DissipatingNetwork1d`]. Successful steps do
/// not rebuild the duct, clone the flow or allocate a new numerical workspace.
pub struct DissipatingDuct<L: TransportLaw> {
    geometry: Duct,
    initial_gas: FrozenGas,
    solver: DissipatingNetwork1d<L>,
    centers_m: Vec<f64>,
    center_areas_m2: Vec<f64>,
}

impl NonlinearDuct {
    /// Attach molecular physics to a newly initialized constant-area bore.
    ///
    /// The law must retain exactly the same original NASA-9 species models,
    /// canonical order and standard-state conventions as this duct. Its
    /// elemental/charge declaration becomes the composition owner's original
    /// source. Initial mass fractions come from the retained frozen composition;
    /// no species is invented, dropped or equilibrated. Both owners use sensible
    /// energy relative to the same Tmin, so momentum and stored total energy
    /// transfer unchanged. Physical formation-energy totals agree to roundoff;
    /// subtracting their large rounded references would corrupt sensible heat.
    ///
    /// Attach before the first accepted step. This explicit initialization
    /// boundary prevents resetting an evolved clock or claiming that previously
    /// inviscid motion already included molecular loss. Tapers, reservoirs and
    /// unsupported molecular geometry retain the gas owner's named refusals.
    pub fn with_molecular_transport<L: TransportLaw>(
        self,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingDuct<L>, DissipatingDuctError> {
        poll(&mut cancelled)?;
        if self.time().value() != 0.0 {
            return Err(DissipatingDuctError::Invalid(
                "attach molecular transport before advancing the nonlinear duct",
            ));
        }
        let original = self.solver.gas();
        if original.source().components() != law.source().models() {
            return Err(DissipatingDuctError::Invalid(
                "molecular law and nonlinear duct have different original species models",
            ));
        }
        let gas = MixtureGas::new(law.source().clone(), original.domain(), &mut cancelled)?;
        let mole = original.source().mole_composition().fractions();
        let masses = original.source().components();
        let mean_mass = compensated(
            mole.iter()
                .zip(masses)
                .map(|(&x, species)| x * species.molar_mass().value()),
        );
        if !mean_mass.is_finite() || mean_mass <= 0.0 {
            return Err(DissipatingDuctError::Invalid(
                "initial mixture molar mass is not representable",
            ));
        }
        let mut fractions = reserve(mole.len())?;
        for (&x, species) in mole.iter().zip(masses) {
            fractions.push(x * (species.molar_mass().value() / mean_mass));
        }
        // Assign the final remainder to the largest constituent, so a trace
        // species never becomes the subtraction of two bulk densities.
        let largest = fractions
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .expect("admitted nonempty source")
            .0;
        let mut states = reserve(self.solver.states().len())?;
        for state in self.solver.states() {
            poll(&mut cancelled)?;
            let mut densities = reserve(fractions.len())?;
            for &fraction in &fractions {
                densities.push(state.density * fraction);
            }
            densities[largest] = state.density
                - compensated(
                    densities
                        .iter()
                        .enumerate()
                        .filter(|&(i, _)| i != largest)
                        .map(|(_, &density)| density),
                );
            if densities.iter().any(|x| !x.is_finite() || *x <= 0.0) {
                return Err(DissipatingDuctError::Invalid(
                    "initial partial density is not positive and representable",
                ));
            }
            states.push(Conserved {
                partial_densities: densities,
                momentum: state.momentum,
                energy: state.energy,
            });
        }
        let branch = MixtureEuler1d::new_on_mesh(
            gas,
            self.solver.mesh().clone(),
            states,
            self.solver.boundaries(),
            &mut cancelled,
        )?;
        let mut branches = reserve(1)?;
        branches.push(branch);
        let network = MixtureNetwork1d::new(branches, Vec::new(), Vec::new(), &mut cancelled)?;
        let solver = DissipatingNetwork1d::new(network, law, Vec::new(), &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(DissipatingDuct {
            geometry: self.geometry,
            initial_gas: original.clone(),
            solver,
            centers_m: self.centers_m,
            center_areas_m2: self.center_areas_m2,
        })
    }
}

impl<L: TransportLaw> DissipatingDuct<L> {
    /// Original immutable cylinder/profile geometry.
    pub const fn geometry(&self) -> &Duct {
        &self.geometry
    }
    /// Exact initial frozen source and composition, before molecular evolution.
    pub const fn initial_gas(&self) -> &FrozenGas {
        &self.initial_gas
    }
    /// Accepted source-bearing fields, totals and complete-step reports.
    pub const fn solver(&self) -> &DissipatingNetwork1d<L> {
        &self.solver
    }
    /// Transfer the accepted molecular/advection owner without reconstructing it.
    pub fn into_solver(self) -> DissipatingNetwork1d<L> {
        self.solver
    }
    /// Accepted clock after both molecular halves and the wave step complete.
    pub fn time(&self) -> Time {
        self.solver.time()
    }
    /// Advance one common CFL-limited physical interval.
    ///
    /// Molecular substeps cover the exact half-duration of the accepted wave
    /// step. A changed wave duration restarts the full composition. Refusal or
    /// cancellation preserves the previous gas field, time and accepted report.
    /// The report retains molecular and advective changes separately, together
    /// with complete species and physical-energy balances.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: DissipatingStepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DissipatingNetworkReport, DissipatingDuctError> {
        Ok(self.solver.step(maximum_duration, controls, cancelled)?)
    }
    /// Evaluate one accepted physical cell using its current composition.
    pub fn sample(&self, cell: usize) -> Result<NonlinearDuctSample, DissipatingDuctError> {
        let state = self.solver.branches()[0]
            .cell(cell)
            .ok_or(DissipatingDuctError::Invalid(
                "sample cell is outside the duct mesh",
            ))?;
        let primitive = self.solver.gas().primitive(state, || false)?;
        let area_m2 = self.center_areas_m2[cell];
        let volume_velocity_m3_s = area_m2 * primitive.velocity;
        let mass_flow_kg_s = primitive.density * volume_velocity_m3_s;
        if !volume_velocity_m3_s.is_finite() || !mass_flow_kg_s.is_finite() {
            return Err(DissipatingDuctError::Invalid(
                "duct volume or mass flow is unrepresentable",
            ));
        }
        Ok(NonlinearDuctSample {
            axial_position_m: self.centers_m[cell],
            area_m2,
            pressure_pa: primitive.pressure,
            density_kg_m3: primitive.density,
            temperature_k: primitive.temperature,
            velocity_m_s: primitive.velocity,
            sound_speed_m_s: primitive.sound_speed,
            volume_velocity_m3_s,
            mass_flow_kg_s,
        })
    }
    /// Snapshot accepted axial observables in original mesh order.
    pub fn samples(&self) -> Result<Vec<NonlinearDuctSample>, DissipatingDuctError> {
        let mut samples = reserve(self.centers_m.len())?;
        for cell in 0..self.centers_m.len() {
            samples.push(self.sample(cell)?);
        }
        Ok(samples)
    }
}

fn reserve<T>(capacity: usize) -> Result<Vec<T>, DissipatingDuctError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| MixtureError::Allocation)?;
    Ok(values)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipatingDuctError> {
    if cancelled() {
        Err(MixtureError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn compensated(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut correction) = (0.0, 0.0);
    for value in values {
        let adjusted = value - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    sum
}
