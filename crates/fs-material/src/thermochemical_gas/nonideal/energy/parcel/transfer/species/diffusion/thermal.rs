//! Source-bound nonisothermal PR diffusion with reciprocal heat transport.
//!
//! Actual Delta(mu_i/T) and a single chosen h*_i define the species effort.
//! h* combines the original NASA divided caloric with symmetric endpoint PR
//! partial residual enthalpy. The latter is a finite-face approximation, not
//! an exact Gibbs--Helmholtz integral over a finite temperature jump.
//! The native species owner carries donor partial enthalpy; additional heat
//! converts that carriage to the SAME h*, and adds Fourier plus reciprocal
//! Dufour heat. At equal T the species law reduces to the original PR face,
//! but centered residual enthalpy can differ from upwind energy carriage.
//!
//! Positive source-bound binary MS coefficients and zero-species-flux thermal
//! conductivity form a positive Schur complement on the mass-flux subspace.
//! Signed H_i are separately authored heats of transport, never inferred from
//! species mass, enthalpy, viscosity or identity. Their common offset is removed.
//! Every coefficient remains T/p-only and must qualify the intended mixtures.
//! PR thermodynamics does not create missing dense-mixture kinetic authority.
//!
//! Convention: Van-Brunt, Farrell & Monroe, arXiv:2109.05082v2, Eq.2.8 and
//! Sections 3,4.5. Local entropy pairing does not replace finite-window audits.

mod constitutive;
/// Simultaneous conservative real-gas species/heat transport.
pub mod network;

use super::{
    BinaryDiffusionMaterial, GasDiffusionGeometry, GasEnergyError, IsochoricRealGasParcel,
    RealGasDiffusionModel, RealGasDiffusionSolve, RealGasDiffusionWorkspace, SpeciesDiffusionFlux,
    SpeciesId, checked, diffusion_error, invalid, log_ratio, poll, positive, source_error,
};
use crate::state_point::{
    ScalarAdmissibility, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use crate::thermochemical_gas::{
    energy::parcel::transfer::species::diffusion::thermal::{
        GasThermodiffusionMaterial, GasThermodiffusionPairSample, GasThermodiffusionSample,
        HEAT_OF_TRANSPORT_DIMS, HEAT_OF_TRANSPORT_PROPERTY, PinnedScalar,
        SpeciesHeatOfTransportMaterial, SpeciesHeatOfTransportSample,
    },
    nonideal::ThermochemicalPengRobinsonModel,
};
use fs_math::dd::Dd;

/// The same numerical/endpoint bounds as real-gas diffusion. Here the declared
/// T difference bounds nonisothermal face resolution, not an isothermal model.
pub type RealGasThermodiffusionSolve = RealGasDiffusionSolve;
pub type RealGasThermodiffusionWorkspace = RealGasDiffusionWorkspace;

#[derive(Clone, Debug)]
pub struct RealGasThermodiffusionModel<'a> {
    molecular: RealGasDiffusionModel<'a>,
    heats: Vec<PinnedScalar<'a>>,
    conductivity: PinnedScalar<'a>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct RealGasThermodiffusionFlux {
    pub species: Vec<SpeciesDiffusionFlux>,
    pub pairs: Vec<GasThermodiffusionPairSample>,
    pub heats: Vec<SpeciesHeatOfTransportSample>,
    /// Conductivity at zero diffusive species flux, not bulk measured k with
    /// an unspecified diffusion constraint.
    pub conductivity: GasThermodiffusionSample,
    pub face_temperature_k: f64,
    pub coefficient_pressure_pa: f64,
    pub molar_concentration_mol_m3: f64,
    /// Actual endpoint logarithmic fugacity increments, including T dependence.
    pub log_fugacity_changes: Vec<f64>,
    /// Chosen symmetric PR partial residual face enthalpy [J/mol].
    pub face_residual_partial_enthalpy_j_mol: Vec<f64>,
    pub projected_driving_forces_per_m: Vec<f64>,
    pub maximum_gradient_residual_per_m: f64,
    pub mass_frame_residual_kg_m2_s: f64,
    pub fourier_heat_flux_w_m2: f64,
    pub reciprocal_heat_flux_w_m2: f64,
    pub enthalpy_correction_flux_w_m2: f64,
    /// Add once to original donor species enthalpy, before endpoint inversion.
    pub additional_heat_rate_left_to_right_w: f64,
    /// Positive original pair friction plus Fourier entropy.
    pub entropy_production_w_k: f64,
    /// Independent actual donor-energy minus chemical-potential source power.
    pub source_entropy_production_w_k: f64,
}

impl<'a> RealGasThermodiffusionModel<'a> {
    /// All three coefficient families must be explicit qualified material data.
    /// The conductivity is the zero-species-flux Schur complement; H_i are J/kg.
    pub fn new(
        reference: &ThermochemicalPengRobinsonModel,
        pairs: &[BinaryDiffusionMaterial<'a>],
        conductivity: GasThermodiffusionMaterial<'a>,
        heats: &[SpeciesHeatOfTransportMaterial<'a>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasEnergyError> {
        let molecular = RealGasDiffusionModel::new(reference, pairs, &mut cancelled)?;
        if heats.len() != molecular.species().len() {
            return Err(invalid(
                "one original heat-of-transport source per PR species required",
            ));
        }
        let mut selected = Vec::with_capacity(heats.len());
        for species in molecular.species() {
            poll(&mut cancelled)?;
            let mut matches = heats.iter().filter(|s| &s.species == species);
            let source = matches
                .next()
                .ok_or(invalid("missing PR heat-of-transport source"))?;
            if matches.next().is_some() {
                return Err(invalid("duplicate PR heat-of-transport source"));
            }
            selected.push(PinnedScalar::new(
                source.source.clone(),
                HEAT_OF_TRANSPORT_PROPERTY,
                HEAT_OF_TRANSPORT_DIMS,
                ScalarAdmissibility::Finite,
            )?);
        }
        let conductivity = PinnedScalar::new(
            conductivity,
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            molecular,
            heats: selected,
            conductivity,
        })
    }
    pub fn species(&self) -> &[SpeciesId] {
        self.molecular.species()
    }
    pub fn workspace(&self) -> Result<RealGasThermodiffusionWorkspace, GasEnergyError> {
        self.molecular.workspace()
    }
    fn admit_pair(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        controls: RealGasThermodiffusionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        controls.validate()?;
        if (left.state().temperature_k() - right.state().temperature_k()).abs()
            > controls.maximum_temperature_difference_k
            || (left.state().pressure_pa() - right.state().pressure_pa()).abs()
                > controls.maximum_pressure_difference_pa
        {
            return Err(invalid(
                "real-gas thermal endpoints exceed authored face support",
            ));
        }
        self.molecular.admit_endpoint(left, cancelled)?;
        self.molecular.admit_endpoint(right, cancelled)
    }
    /// Qualify every actual final endpoint and central coefficient state.
    pub fn admit_endpoints(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        controls: RealGasThermodiffusionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        self.admit_pair(left, right, controls, &mut cancelled)?;
        let t = [left.state().temperature_k(), right.state().temperature_k()];
        let p = [left.state().pressure_pa(), right.state().pressure_pa()];
        for (_, _, source) in &self.molecular.pairs {
            source.sample(t, p, &mut cancelled)?;
        }
        for source in &self.heats {
            source.sample(t, p, &mut cancelled)?;
        }
        self.conductivity.sample(t, p, &mut cancelled)?;
        poll(&mut cancelled)
    }
    /// Original density, composition, partial PR calorics and fugacity remain
    /// attached to each endpoint. No held-T ideal surrogate or new root is used.
    #[allow(clippy::too_many_arguments)]
    pub fn face_flux(
        &self,
        left: &IsochoricRealGasParcel,
        right: &IsochoricRealGasParcel,
        geometry: GasDiffusionGeometry,
        controls: RealGasThermodiffusionSolve,
        work: &mut RealGasThermodiffusionWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RealGasThermodiffusionFlux, GasEnergyError> {
        poll(&mut cancelled)?;
        let n = self.species().len();
        if !positive(geometry.area_m2)
            || !positive(geometry.distance_m)
            || work.binary.len() != n * n
        {
            return Err(invalid(
                "positive real-gas thermal geometry and matching workspace required",
            ));
        }
        self.admit_pair(left, right, controls, &mut cancelled)?;
        let temperatures = [left.state().temperature_k(), right.state().temperature_k()];
        let pressures = [left.state().pressure_pa(), right.state().pressure_pa()];
        let xl = left.model().mechanics().mole_fractions();
        let xr = right.model().mechanics().mole_fractions();
        let mut x: Vec<_> = xl.iter().zip(xr).map(|(a, b)| 0.5 * a + 0.5 * b).collect();
        fs_thermal::enthalpy::advection::species::close_generated_simplex(&mut x)
            .map_err(source_error)?;
        let concentration = 0.5
            * (left.density_kg_m3() / left.model().mechanics().molar_mass().value())
            + 0.5 * (right.density_kg_m3() / right.model().mechanics().molar_mass().value());
        let mut log_fugacity_changes = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            log_fugacity_changes.push(checked(
                Dd::from_f64(log_ratio(xr[i], xl[i]))
                    + Dd::from_f64(log_ratio(pressures[1], pressures[0]))
                    + Dd::from_f64(right.state().root().log_fugacity_coefficients()[i])
                    - Dd::from_f64(left.state().root().log_fugacity_coefficients()[i]),
            )?);
        }
        let partial_left = left
            .model()
            .partial_calorics_at_density(left.density_kg_m3(), temperatures[0])?;
        poll(&mut cancelled)?;
        let partial_right = right
            .model()
            .partial_calorics_at_density(right.density_kg_m3(), temperatures[1])?;
        let mut residual = [[0.0; 32]; 2];
        for (side, partial) in [&partial_left, &partial_right].into_iter().enumerate() {
            for (i, component) in partial.departure().components().iter().enumerate() {
                residual[side][i] = component.residual_enthalpy_j_mol();
            }
        }
        work.binary.fill(0.0);
        let mut pairs = Vec::with_capacity(self.molecular.pairs.len());
        for (i, j, source) in &self.molecular.pairs {
            let sample = source.sample(temperatures, pressures, &mut cancelled)?;
            work.binary[i * n + j] = sample.value_si;
            work.binary[j * n + i] = sample.value_si;
            pairs.push(GasThermodiffusionPairSample {
                species: [self.species()[*i].clone(), self.species()[*j].clone()],
                sample,
            });
        }
        let conductivity = self
            .conductivity
            .sample(temperatures, pressures, &mut cancelled)?;
        let mut heats = Vec::with_capacity(n);
        let mut h = [0.0; 32];
        for (i, source) in self.heats.iter().enumerate() {
            let sample = source.sample(temperatures, pressures, &mut cancelled)?;
            h[i] = sample.value_si;
            heats.push(SpeciesHeatOfTransportSample {
                species: self.species()[i].clone(),
                sample,
                centered_heat_j_kg: 0.0,
            });
        }
        let value = constitutive::evaluate(
            &self.molecular.diffusion,
            self.molecular
                .reference
                .source()
                .thermodynamics()
                .components(),
            &x,
            temperatures,
            &log_fugacity_changes,
            [&residual[0][..n], &residual[1][..n]],
            concentration,
            &work.binary,
            &h[..n],
            conductivity.value_si,
            geometry.distance_m,
            controls.diffusion,
            controls.entropy_pairing_tolerance_w_m3_k,
            &mut work.diffusion,
            &mut cancelled,
        )
        .map_err(diffusion_error)?;
        let mut species = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            heats[i].centered_heat_j_kg = value.centered_heats[i];
            let flux = value.kinetic.diffusion.mass_fluxes()[i];
            let rate = checked(Dd::from_f64(flux) * Dd::from_f64(geometry.area_m2))?;
            if flux != 0.0 && rate == 0.0 {
                return Err(invalid("real-gas thermal area discarded a trace mass flux"));
            }
            species.push(SpeciesDiffusionFlux {
                species: self.species()[i].clone(),
                molar_flux_mol_m2_s: value.kinetic.diffusion.molar_fluxes()[i],
                mass_rate_left_to_right_kg_s: rate,
            });
        }
        let heat = checked(Dd::from_f64(value.additional) * Dd::from_f64(geometry.area_m2))?;
        let entropy =
            checked(Dd::from_f64(value.entropy_from_dissipation) * Dd::from_f64(geometry.area_m2))?;
        let source_entropy =
            checked(Dd::from_f64(value.entropy_from_sources) * Dd::from_f64(geometry.area_m2))?;
        if (value.additional != 0.0 && heat == 0.0)
            || (value.entropy_from_dissipation > 0.0 && entropy == 0.0)
        {
            return Err(invalid("real-gas thermal area discarded heat or entropy"));
        }
        poll(&mut cancelled)?;
        Ok(RealGasThermodiffusionFlux {
            species,
            pairs,
            heats,
            conductivity,
            face_temperature_k: 0.5 * temperatures[0] + 0.5 * temperatures[1],
            coefficient_pressure_pa: 0.5 * pressures[0] + 0.5 * pressures[1],
            molar_concentration_mol_m3: concentration,
            log_fugacity_changes,
            face_residual_partial_enthalpy_j_mol: value.face_residual_enthalpies,
            projected_driving_forces_per_m: value.kinetic.forces,
            maximum_gradient_residual_per_m: value
                .kinetic
                .diffusion
                .maximum_gradient_residual_per_m(),
            mass_frame_residual_kg_m2_s: value.kinetic.diffusion.mass_frame_residual_kg_m2_s(),
            fourier_heat_flux_w_m2: value.fourier,
            reciprocal_heat_flux_w_m2: value.reciprocal,
            enthalpy_correction_flux_w_m2: value.correction,
            additional_heat_rate_left_to_right_w: heat,
            entropy_production_w_k: entropy,
            source_entropy_production_w_k: source_entropy,
        })
    }
}
