//! Fugacity-driven molecular diffusion on original homogeneous PR gas states.
//!
//! At a common face temperature, g_i = Delta ln(x_i p phi_i)/distance.
//! The mass-frame force is d_i = x_i g_i - Y_i sum_j(x_j g_j). Finite
//! fugacity differences need not satisfy Gibbs--Duhem: this projection removes
//! the mechanical mode explicitly. Adding lambda*M_i to g leaves d unchanged.
//! With sum_i(M_i N_i)=0, -R sum_i(N_i g_i) equals the original positive
//! Maxwell--Stefan pair friction. Every original row and this entropy pairing
//! are independently checked; independent Fick laws are not substituted.
//!
//! Both PR mixtures are evaluated at the common T and their actual densities,
//! retaining their original homogeneous branches without selecting a new PT
//! root. A nonzero authored temperature allowance is an ISOTHERMAL face
//! approximation, not Soret, Fourier or thermal-transpiration physics. Sources
//! admit the actual endpoints and the held-T evaluations separately. Spatial
//! and time refinement remain required. Mechanical stability is not a proof
//! of single-phase thermodynamic stability or absence of phase separation.
//!
//! Binary MS coefficients are explicit original material claims; no ideal pD,
//! viscosity inference or dense-fluid correction is manufactured. See Krishna
//! & van Baten, Chem. Eng. Sci. 153 (2016) 174--187,
//! doi:10.1016/j.ces.2016.07.025; Bothe, arXiv:1007.1775, equations 13--18.

mod force;
/// Conservative simultaneous diffusion on finite real-gas networks.
pub mod network;

use super::{GasEnergyError, IsochoricRealGasParcel, invalid, poll, positive};
use crate::state_point::ScalarAdmissibility;
use crate::thermochemical_gas::energy::parcel::transfer::mixed_model_from_portions;
use crate::thermochemical_gas::energy::parcel::transfer::species::{
    diffusion::{
        BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
        GasDiffusionGeometry, MAX_DIFFUSION_SPECIES, SpeciesDiffusionFlux,
        thermal::{GasThermodiffusionMaterial, GasThermodiffusionPairSample, PinnedScalar},
    },
    entropy::log_ratio,
};
use crate::thermochemical_gas::energy::updated_point;
use crate::thermochemical_gas::nonideal::energy::parcel::transfer::mixing;
use crate::thermochemical_gas::nonideal::{
    ResolvedThermochemicalRealGasState, ThermochemicalPengRobinsonModel,
};
use fs_math::dd::Dd;
use fs_thermochem::SpeciesId;
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanErrorV1, MaxwellStefanSolveV1, MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};

/// Explicit local isothermal approximation and original-equation tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealGasDiffusionSolve {
    /// Zero requires exactly equal actual endpoint temperatures.
    pub maximum_temperature_difference_k: f64,
    /// Bounds both actual and held-face-temperature pressure differences.
    pub maximum_pressure_difference_pa: f64,
    pub diffusion: MaxwellStefanSolveV1,
    /// Absolute allowance in reconstructed force power versus pair friction.
    /// The MS relative tolerance times their scale is added separately.
    pub entropy_pairing_tolerance_w_m3_k: f64,
}
impl RealGasDiffusionSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        if ![
            self.maximum_temperature_difference_k,
            self.maximum_pressure_difference_pa,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
            || ![
                self.diffusion.relative_tolerance,
                self.diffusion.absolute_gradient_tolerance_per_m,
                self.diffusion.absolute_mass_flux_tolerance_kg_m2_s,
                self.entropy_pairing_tolerance_w_m3_k,
            ]
            .iter()
            .all(|v| positive(*v))
        {
            return Err(invalid(
                "invalid real-gas diffusion approximation or residual controls",
            ));
        }
        Ok(())
    }
}

/// Prepared immutable original PR/NASA catalog and selected binary MS claims.
/// Every prepared species must be positive in both endpoints; no trace floors
/// replace one-sided logarithmic zeros or qualify a missing interaction pair.
/// Binary claims currently support temperature/pressure dependence only. The
/// original claim must qualify that composition-independent coefficient law
/// over the intended mixtures; PR fugacity does not supply or validate a dense
/// mixture diffusivity. Composition/humidity-dependent claims refuse here.
#[derive(Clone, Debug)]
pub struct RealGasDiffusionModel<'a> {
    reference: ThermochemicalPengRobinsonModel,
    diffusion: MaxwellStefanV1,
    pairs: Vec<(usize, usize, PinnedScalar<'a>)>,
}

/// Reusable dense friction scratch, without physical state or accepted history.
#[derive(Clone, Debug)]
pub struct RealGasDiffusionWorkspace {
    binary: Vec<f64>,
    diffusion: MaxwellStefanWorkspaceV1,
}
impl RealGasDiffusionWorkspace {
    pub fn workspace_bytes(&self) -> usize {
        self.binary.capacity() * core::mem::size_of::<f64>() + self.diffusion.workspace_bytes()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasDiffusionFlux {
    pub species: Vec<SpeciesDiffusionFlux>,
    /// Same three-point, original-claim receipts used by thermal diffusion.
    pub pairs: Vec<GasThermodiffusionPairSample>,
    pub face_temperature_k: f64,
    /// Mean of the actual physical endpoint pressures used for binary data.
    pub coefficient_pressure_pa: f64,
    /// Actual-density PR states at the common temperature used by fugacity.
    pub isothermal_traces: [ResolvedThermochemicalRealGasState; 2],
    /// Arithmetic mean of actual endpoint molar densities, never p/(RT).
    pub molar_concentration_mol_m3: f64,
    pub log_fugacity_gradients_per_m: Vec<f64>,
    /// Remove (M_i/M_pivot)*g_pivot with the largest-mass species as pivot.
    /// This eliminates only the mechanical gauge before entropy arithmetic.
    pub reduced_log_fugacity_gradients_per_m: Vec<f64>,
    pub projected_driving_forces_per_m: Vec<f64>,
    pub maximum_gradient_residual_per_m: f64,
    pub mass_frame_residual_kg_m2_s: f64,
    /// Original positive pair friction; not extra heat deposited into gas.
    pub drag_entropy_production_w_m3_k: f64,
    /// Independently reconstructed -R sum_i(N_i reduced_grad ln f_i).
    /// The separately admitted mass-frame residual carries no mechanical work.
    pub force_entropy_production_w_m3_k: f64,
}

fn diffusion_error(error: MaxwellStefanErrorV1) -> GasEnergyError {
    if error == MaxwellStefanErrorV1::Cancelled {
        GasEnergyError::Cancelled
    } else {
        source_error(error)
    }
}
fn source_error(error: impl core::fmt::Display) -> GasEnergyError {
    crate::thermochemical_gas::failure("real-gas molecular diffusion", error).into()
}
fn checked(value: Dd) -> Result<f64, GasEnergyError> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(invalid("real-gas diffusion reduction overflowed"))
    }
}

impl<'a> RealGasDiffusionModel<'a> {
    pub fn new(
        reference: &ThermochemicalPengRobinsonModel,
        pairs: &[BinaryDiffusionMaterial<'a>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasEnergyError> {
        poll(&mut cancelled)?;
        let species = reference.source().thermodynamics().components();
        let n = species.len();
        if !(1..=MAX_DIFFUSION_SPECIES).contains(&n) || pairs.len() != n * (n - 1) / 2 {
            return Err(invalid(
                "real-gas diffusion needs bounded complete original species/pair sources",
            ));
        }
        let diffusion = MaxwellStefanV1::new(
            species
                .iter()
                .map(|s| (s.species().clone(), s.molar_mass()))
                .collect(),
        )
        .map_err(diffusion_error)?;
        let mut pinned = Vec::with_capacity(pairs.len());
        for i in 0..n {
            for j in i + 1..n {
                poll(&mut cancelled)?;
                let mut selected = pairs.iter().filter(|pair| {
                    pair.species
                        == [
                            diffusion.species()[i].clone(),
                            diffusion.species()[j].clone(),
                        ]
                        || pair.species
                            == [
                                diffusion.species()[j].clone(),
                                diffusion.species()[i].clone(),
                            ]
                });
                let pair = selected
                    .next()
                    .ok_or(invalid("missing real-gas binary MS source"))?;
                if selected.next().is_some() {
                    return Err(invalid("duplicate real-gas binary MS source"));
                }
                pinned.push((
                    i,
                    j,
                    PinnedScalar::new(
                        GasThermodiffusionMaterial {
                            material: pair.material,
                            context: pair.context,
                            selection: pair.selection.clone(),
                        },
                        BINARY_DIFFUSIVITY_PROPERTY,
                        BINARY_DIFFUSIVITY_DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )?,
                ));
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            reference: reference.clone(),
            diffusion,
            pairs: pinned,
        })
    }

    pub fn species(&self) -> &[SpeciesId] {
        self.diffusion.species()
    }

    pub fn workspace(&self) -> Result<RealGasDiffusionWorkspace, GasEnergyError> {
        let n = self.species().len();
        Ok(RealGasDiffusionWorkspace {
            binary: vec![0.0; n * n],
            diffusion: self.diffusion.workspace().map_err(diffusion_error)?,
        })
    }

    // Reuse the original source joins: no equality of species labels alone can
    // admit changed NASA/transport/critical/interaction/pack authority.
    pub(super) fn admit_endpoint(
        &self,
        parcel: &IsochoricRealGasParcel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        poll(cancelled)?;
        let components = parcel.model().source().thermodynamics().components();
        if !components.iter().map(|s| s.species()).eq(self.species())
            || parcel
                .model()
                .mechanics()
                .mole_fractions()
                .iter()
                .any(|x| !positive(*x))
        {
            return Err(invalid(
                "real-gas diffusion requires the prepared positive species catalog on both sides",
            ));
        }
        let joined = mixed_model_from_portions(
            &[
                (self.reference.source(), 1.0),
                (parcel.model().source(), 1.0),
            ],
            2.0,
        )?;
        mixing::model_for_source(&[&self.reference, parcel.model()], joined)?;
        for (_, _, pair) in &self.pairs {
            poll(cancelled)?;
            pair.at(parcel.state().temperature_k(), parcel.state().pressure_pa())?;
        }
        Ok(())
    }

    pub(super) fn admit_pair(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        controls: RealGasDiffusionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        controls.validate()?;
        if (left.state().temperature_k() - right.state().temperature_k()).abs()
            > controls.maximum_temperature_difference_k
            || (left.state().pressure_pa() - right.state().pressure_pa()).abs()
                > controls.maximum_pressure_difference_pa
        {
            return Err(invalid(
                "real-gas endpoints exceed declared isothermal face approximation",
            ));
        }
        self.admit_endpoint(left, cancelled)?;
        self.admit_endpoint(right, cancelled)
    }

    /// Sample one open face without changing either physical inventory.
    /// Omit closed faces; positive area/distance and a model-sized workspace
    /// are explicit. The source coefficients and original MS matrix own all
    /// transport; PR fugacity supplies its thermodynamic driving force only.
    #[allow(clippy::too_many_arguments)]
    pub fn face_flux(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        geometry: GasDiffusionGeometry,
        controls: RealGasDiffusionSolve,
        work: &mut RealGasDiffusionWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasDiffusionFlux, GasEnergyError> {
        poll(&mut cancelled)?;
        let n = self.species().len();
        if !positive(geometry.area_m2)
            || !positive(geometry.distance_m)
            || work.binary.len() != n * n
        {
            return Err(invalid(
                "real-gas diffusion needs positive geometry and matching workspace",
            ));
        }
        self.admit_pair(left, right, controls, &mut cancelled)?;
        let temperatures = [left.state().temperature_k(), right.state().temperature_k()];
        let pressures = [left.state().pressure_pa(), right.state().pressure_pa()];
        let temperature = 0.5 * temperatures[0] + 0.5 * temperatures[1];
        let pressure = 0.5 * pressures[0] + 0.5 * pressures[1];
        let trace = |parcel: &IsochoricRealGasParcel| {
            let point = updated_point(parcel.point(), temperature, parcel.state().pressure_pa())?;
            parcel
                .model()
                .evaluate_density_temperature(&point, parcel.density_kg_m3())
                .map_err(GasEnergyError::from)
        };
        let traces = [trace(left)?, trace(right)?];
        if (traces[1].pressure_pa() - traces[0].pressure_pa()).abs()
            > controls.maximum_pressure_difference_pa
        {
            return Err(invalid(
                "held-temperature PR traces exceed the pressure approximation",
            ));
        }
        poll(&mut cancelled)?;
        let xl = left.model().mechanics().mole_fractions();
        let xr = right.model().mechanics().mole_fractions();
        let mut x: Vec<_> = xl.iter().zip(xr).map(|(a, b)| 0.5 * a + 0.5 * b).collect();
        fs_thermal::enthalpy::advection::species::close_generated_simplex(&mut x)
            .map_err(source_error)?;
        let concentration = 0.5
            * (left.density_kg_m3() / left.model().mechanics().molar_mass().value())
            + 0.5 * (right.density_kg_m3() / right.model().mechanics().molar_mass().value());
        if !positive(concentration) {
            return Err(invalid("real-gas molar concentration is unrepresentable"));
        }
        let log_pressure = log_ratio(traces[1].pressure_pa(), traces[0].pressure_pa());
        let mut gradients = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let increment = checked(
                Dd::from_f64(log_ratio(xr[i], xl[i]))
                    + Dd::from_f64(log_pressure)
                    + Dd::from_f64(traces[1].root().log_fugacity_coefficients()[i])
                    - Dd::from_f64(traces[0].root().log_fugacity_coefficients()[i]),
            )?;
            let value = increment / geometry.distance_m;
            if !value.is_finite() || (increment != 0.0 && value == 0.0) {
                return Err(invalid("real-gas fugacity force is unrepresentable"));
            }
            gradients.push(value);
        }
        work.binary.fill(0.0);
        let mut pairs = Vec::with_capacity(self.pairs.len());
        for (i, j, source) in &self.pairs {
            let sample = source.sample(temperatures, pressures, &mut cancelled)?;
            work.binary[i * n + j] = sample.value_si;
            work.binary[j * n + i] = sample.value_si;
            pairs.push(GasThermodiffusionPairSample {
                species: [self.species()[*i].clone(), self.species()[*j].clone()],
                sample,
            });
        }
        let evaluated = force::evaluate(
            &self.diffusion,
            &x,
            &gradients,
            concentration,
            &work.binary,
            controls.diffusion,
            controls.entropy_pairing_tolerance_w_m3_k,
            &mut work.diffusion,
            &mut cancelled,
        )
        .map_err(diffusion_error)?;
        let mut species = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let molar = evaluated.diffusion.molar_fluxes()[i];
            let mass = evaluated.diffusion.mass_fluxes()[i] * geometry.area_m2;
            if !mass.is_finite() || (evaluated.diffusion.mass_fluxes()[i] != 0.0 && mass == 0.0) {
                return Err(invalid(
                    "real-gas diffusion area-integrated mass is unrepresentable",
                ));
            }
            species.push(SpeciesDiffusionFlux {
                species: self.species()[i].clone(),
                molar_flux_mol_m2_s: molar,
                mass_rate_left_to_right_kg_s: mass,
            });
        }
        poll(&mut cancelled)?;
        Ok(RealGasDiffusionFlux {
            species,
            pairs,
            face_temperature_k: temperature,
            coefficient_pressure_pa: pressure,
            isothermal_traces: traces,
            molar_concentration_mol_m3: concentration,
            log_fugacity_gradients_per_m: gradients,
            reduced_log_fugacity_gradients_per_m: evaluated.reduced_gradients,
            projected_driving_forces_per_m: evaluated.forces,
            maximum_gradient_residual_per_m: evaluated.diffusion.maximum_gradient_residual_per_m(),
            mass_frame_residual_kg_m2_s: evaluated.diffusion.mass_frame_residual_kg_m2_s(),
            drag_entropy_production_w_m3_k: evaluated.drag_entropy,
            force_entropy_production_w_m3_k: evaluated.force_entropy,
        })
    }
}

#[cfg(test)]
mod tests;
