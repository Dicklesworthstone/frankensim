//! Molecular Rayleigh scattering of the same gas used by flow and acoustics.
//!
//! Optical data are explicit: a reference refractivity `n - 1`, molecular
//! number density and wavelength, a depolarization ratio, and a source identity.
//! No chemical name implies optical constants. Component cross sections add
//! by **number fraction**; averaging refractive indices before squaring would
//! introduce unphysical coherent cross terms between independent molecules.
//!
//! The reference cross section is the Lorentz–Lorenz/King expression
//! `24 pi³ [(n²-1)/(n²+2)]² F / (N² lambda⁴)`, with
//! `F = (6+3 delta)/(6-7 delta)`. See Bodhaine et al. (1999),
//! <https://journals.ametsoc.org/view/journals/atot/16/11/1520-0426_1999_016_1854_orodc_2_0_co_2.xml>.
//! The unpolarized phase uses the same depolarization convention as the total
//! cross section. Mixture phases are weighted by scattering cross section,
//! not by bulk mass fraction or an averaged depolarization ratio.
//! Phase convention: Lin et al. (2022), Eqs. 42–53,
//! <https://doi.org/10.3389/frsen.2022.880768>.
//!
//! Scope: independent, dilute molecules, elastic unpolarized scattering and
//! **constant molecular polarizability** over 380–780 nm. The exact lambda^-4
//! law is not a dispersion model, a resonance model, or a measured error bound.
//! Depolarization changes the scalar angular law; this module does not transport
//! Stokes polarization. It supplies no absorption or thermal emission.
//!
//! The same molecular references also determine the coherent phase index:
//! `L = (n²-1)/(n²+2) = rho K`, where
//! `K = (N_A/M_mix) sum_i x_i L_ref,i/N_ref,i`. Polarizabilities add before
//! evaluating the phase index; independent scattering cross sections add
//! separately. The King factor is not a correction to the phase index.
//! See Shevchenko et al., *Investigation of diffusive and optical properties of
//! vapour-air mixtures: The benefits of interferometry*, §3.3.1,
//! <https://doi.org/10.1016/j.ces.2020.116433>.

use fs_blake3::{ContentHash, DomainHasher};
use fs_qty::{MolarMass, semantic::Composition};

use crate::MaterialError;
use crate::gas::{ConductivityModel, GasSpec, GasState, R_USSA_1976};
use crate::gas_mixture::{MAX_GAS_MIXTURE_COMPONENTS, mix_gas_states};

/// Exact SI Avogadro constant [mol^-1].
pub const AVOGADRO_PER_MOL: f64 = 6.022_140_76e23;
/// Common reference for the prepared mixture, in vacuum nanometres.
pub const RAYLEIGH_REFERENCE_NM: f64 = 550.0;
/// Declared dilute-gas use boundary, shared by reference and local refractivity.
/// This numerical domain is not an experimental accuracy certificate.
pub const MAX_DILUTE_REFRACTIVITY: f64 = 0.01;

/// Explicit reference optical data for one species or declared pseudocomponent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularRayleigh {
    refractivity: f64,
    reference_number_density_m3: f64,
    reference_wavelength_nm: f64,
    depolarization_ratio: f64,
    source_identity: ContentHash,
    reference_cross_section_m2: f64,
}

impl MolecularRayleigh {
    /// Admit optical reference data, retaining `n - 1` directly so dilute-gas
    /// precision is not lost by adding one and subtracting it again.
    ///
    /// `delta` is the unpolarized 90-degree parallel/perpendicular intensity
    /// ratio. Its passive molecular range is `[0, 6/7)`. The refractivity
    /// window `[0, 0.01]` is an explicit dilute-reference use policy, not an
    /// accuracy certificate. Number density is molecules/m³, not mol/m³.
    pub fn try_new(
        refractivity: f64,
        reference_number_density_m3: f64,
        reference_wavelength_nm: f64,
        depolarization_ratio: f64,
        source_identity: ContentHash,
    ) -> Result<Self, GasScatteringError> {
        if !refractivity.is_finite()
            || !(0.0..=MAX_DILUTE_REFRACTIVITY).contains(&refractivity)
            || !reference_number_density_m3.is_finite()
            || reference_number_density_m3 <= 0.0
            || !depolarization_ratio.is_finite()
            || !(0.0..6.0 / 7.0).contains(&depolarization_ratio)
            || source_identity == ContentHash([0; 32])
        {
            return Err(GasScatteringError::OpticalReference);
        }
        wavelength(reference_wavelength_nm)?;
        let n2_minus_one = refractivity * (2.0 + refractivity);
        let lorentz_lorenz = n2_minus_one / (3.0 + n2_minus_one);
        let lambda = reference_wavelength_nm * 1.0e-9;
        let scaled = (lorentz_lorenz / reference_number_density_m3) / (lambda * lambda);
        let king = (6.0 + 3.0 * depolarization_ratio) / (6.0 - 7.0 * depolarization_ratio);
        let pi = core::f64::consts::PI;
        let reference_cross_section_m2 = (24.0 * pi * pi * pi * king) * scaled * scaled;
        finite_nonnegative(reference_cross_section_m2)?;
        if refractivity > 0.0 && reference_cross_section_m2 == 0.0 {
            return Err(GasScatteringError::Unrepresentable);
        }
        Ok(Self {
            refractivity,
            reference_number_density_m3,
            reference_wavelength_nm,
            depolarization_ratio,
            source_identity,
            reference_cross_section_m2,
        })
    }

    /// Total elastic scattering cross section per molecule [m²].
    pub fn cross_section_m2(self, wavelength_nm: f64) -> Result<f64, GasScatteringError> {
        inverse_fourth(
            self.reference_cross_section_m2,
            self.reference_wavelength_nm,
            wavelength_nm,
        )
    }

    /// Original caller-supplied source identity. Identity does not prove accuracy.
    #[must_use]
    pub const fn source_identity(self) -> ContentHash {
        self.source_identity
    }

    /// The same molecular anisotropy used by the King correction.
    #[must_use]
    pub const fn depolarization_ratio(self) -> f64 {
        self.depolarization_ratio
    }

    /// Molecular Lorentz–Lorenz volume `[(n²-1)/(n²+2)]/N` [m³/molecule].
    /// It is constant across the visible band in this explicit approximation.
    #[must_use]
    pub fn lorentz_lorenz_volume_m3(self) -> f64 {
        let n2_minus_one = self.refractivity * (2.0 + self.refractivity);
        (n2_minus_one / (3.0 + n2_minus_one)) / self.reference_number_density_m3
    }

    fn hash_into(self, hash: &mut DomainHasher) {
        hash.update(self.source_identity.as_bytes());
        for value in [
            self.refractivity,
            self.reference_number_density_m3,
            self.reference_wavelength_nm,
            self.depolarization_ratio,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
    }
}

/// Re-evaluable molecular optics attached to the original ideal-gas specs.
/// Component order is shared by the gas, composition and optical reference rows.
#[derive(Clone, Debug, PartialEq)]
pub struct MolecularRayleighMixture {
    specs: Vec<GasSpec>,
    composition: Composition,
    optical: Vec<MolecularRayleigh>,
    water_component: Option<usize>,
}

impl MolecularRayleighMixture {
    /// Retain the exact shared gas input; no species, composition or missing
    /// optical datum is inferred. Thermodynamic admission occurs at evaluation.
    pub fn try_new(
        specs: &[GasSpec],
        composition: Composition,
        optical: &[MolecularRayleigh],
        water_component: Option<usize>,
    ) -> Result<Self, GasScatteringError> {
        if specs.is_empty()
            || specs.len() > MAX_GAS_MIXTURE_COMPONENTS
            || specs.len() != optical.len()
            || specs.len() != composition.fractions().len()
            || water_component.is_some_and(|index| index >= specs.len())
        {
            return Err(GasScatteringError::Components);
        }
        Ok(Self {
            specs: specs.to_vec(),
            composition,
            optical: optical.to_vec(),
            water_component,
        })
    }

    /// Evaluate original component gases and their original Wilke/WMS mixture,
    /// then attach independent molecular scattering to that exact accepted state.
    /// This retains the original gas law's temperature, pressure and phase limits.
    pub fn evaluate(
        &self,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<RayleighGasState, GasScatteringError> {
        let components = self
            .specs
            .iter()
            .map(|spec| GasState::try_new(spec, temperature_k, pressure_pa))
            .collect::<Result<Vec<_>, _>>()
            .map_err(GasScatteringError::Gas)?;
        let gas = mix_gas_states(&components, &self.composition, self.water_component)
            .map_err(GasScatteringError::Gas)?;
        // Recover the same molar masses as the existing mixture owner. This
        // also keeps mass-to-mole conversion on its exact original arithmetic.
        let masses: Vec<_> = components
            .iter()
            .map(|state| MolarMass::new(R_USSA_1976 / state.specific_gas_constant))
            .collect();
        let mole = self
            .composition
            .to_mole_fractions(&masses)
            .map_err(|_| GasScatteringError::Components)?;
        let mut sigma = 0.0;
        let mut isotropic_sigma = 0.0;
        let mut refractive_volume = 0.0;
        let mut molar_mass = 0.0;
        let mut hash = DomainHasher::new("org.frankensim.material.rayleigh-gas.v2");
        hash.update(&(self.specs.len() as u64).to_le_bytes());
        hash.update(&(self.water_component.map_or(u64::MAX, |i| i as u64)).to_le_bytes());
        for (((&fraction, mass), &optical), spec) in mole
            .fractions()
            .iter()
            .zip(&masses)
            .zip(&self.optical)
            .zip(&self.specs)
        {
            let weighted = fraction * optical.cross_section_m2(RAYLEIGH_REFERENCE_NM)?;
            sigma += weighted;
            refractive_volume += fraction * optical.lorentz_lorenz_volume_m3();
            // Depolarized Rayleigh = w_iso * isotropic + (1-w_iso)*Rayleigh.
            // Weight the PHASE by scattering, never average delta directly.
            let delta = optical.depolarization_ratio;
            isotropic_sigma += weighted * (3.0 * delta / (2.0 + delta));
            molar_mass += fraction * mass.value();
            optical.hash_into(&mut hash);
            for value in [
                fraction,
                spec.molar_mass,
                spec.gamma,
                spec.sutherland_beta,
                spec.sutherland_s,
            ] {
                hash.update(&value.to_bits().to_le_bytes());
            }
            hash.update(&[match spec.conductivity {
                ConductivityModel::Ussa1976AirFit => 0,
                ConductivityModel::Eucken => 1,
            }]);
        }
        finite_nonnegative(sigma)?;
        if !molar_mass.is_finite() || molar_mass <= 0.0 {
            return Err(GasScatteringError::Unrepresentable);
        }
        let mass_scattering_reference_m2_kg = (sigma / molar_mass) * AVOGADRO_PER_MOL;
        finite_nonnegative(mass_scattering_reference_m2_kg)?;
        if sigma > 0.0 && mass_scattering_reference_m2_kg == 0.0 {
            return Err(GasScatteringError::Unrepresentable);
        }
        let mass_refractivity_m3_kg = (refractive_volume / molar_mass) * AVOGADRO_PER_MOL;
        finite_nonnegative(mass_refractivity_m3_kg)?;
        if refractive_volume > 0.0 && mass_refractivity_m3_kg == 0.0 {
            return Err(GasScatteringError::Unrepresentable);
        }
        let isotropic = if sigma == 0.0 {
            0.0
        } else {
            isotropic_sigma / sigma
        };
        let depolarization_ratio = 2.0 * isotropic / (3.0 - isotropic);
        if !depolarization_ratio.is_finite() || !(0.0..6.0 / 7.0).contains(&depolarization_ratio) {
            return Err(GasScatteringError::Unrepresentable);
        }
        for value in [
            gas.temperature,
            gas.pressure,
            gas.density,
            gas.sound_speed,
            gas.dynamic_viscosity,
            gas.thermal_conductivity,
            gas.gamma,
            gas.specific_gas_constant,
            gas.specific_heat_cp,
            gas.prandtl,
            gas.characteristic_impedance,
            gas.water_mole_fraction,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        let result = RayleighGasState {
            gas,
            mass_scattering_reference_m2_kg,
            depolarization_ratio,
            molecular_cross_section_reference_m2: sigma,
            mass_refractivity_m3_kg,
            molar_mass_kg_mol: molar_mass,
            identity: hash.finalize(),
        };
        // Admission proves representability over the entire visible band,
        // including the maximum coefficient at the short-wavelength endpoint.
        result.scattering_per_m(380.0)?;
        result.scattering_per_m(780.0)?;
        Ok(result)
    }

    /// Attach optics to an existing accepted gas only if the original gas owner
    /// reproduces its entire numerical state exactly. A differing transport or
    /// density refuses. `GasState` alone carries no species identity: the caller
    /// must use these same declared component inputs in the simulation. Numeric
    /// equality is not evidence that two thermodynamically identical species
    /// have the same optical properties.
    pub fn bind_gas_state(
        &self,
        accepted: &GasState,
    ) -> Result<RayleighGasState, GasScatteringError> {
        accepted
            .validate_coherence()
            .map_err(GasScatteringError::Gas)?;
        let result = self.evaluate(accepted.temperature, accepted.pressure)?;
        if result.gas != *accepted {
            return Err(GasScatteringError::StateMismatch);
        }
        Ok(result)
    }
}

/// Prepared optics plus the same gas state consumed by flow and acoustics.
/// Spectral evaluation is allocation-free and constant work after preparation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighGasState {
    gas: GasState,
    mass_scattering_reference_m2_kg: f64,
    molecular_cross_section_reference_m2: f64,
    mass_refractivity_m3_kg: f64,
    molar_mass_kg_mol: f64,
    depolarization_ratio: f64,
    identity: ContentHash,
}

impl RayleighGasState {
    /// Original accepted medium, with no optical override of its mass density.
    #[must_use]
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
    /// Source references, ordering, composition, gas-law inputs and accepted state.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }
    /// Mean molar mass using the gas owner's exact mole-fraction conversion.
    #[must_use]
    pub const fn molar_mass_kg_mol(&self) -> f64 {
        self.molar_mass_kg_mol
    }
    /// Number-weighted molecular cross section at 550 nm [m²/molecule].
    #[must_use]
    pub const fn molecular_cross_section_reference_m2(&self) -> f64 {
        self.molecular_cross_section_reference_m2
    }
    /// Prepared mass scattering coefficient at 550 nm [m²/kg].
    #[must_use]
    pub const fn mass_scattering_reference_m2_kg(&self) -> f64 {
        self.mass_scattering_reference_m2_kg
    }
    /// Scattering-weighted effective depolarization, not a fraction-weighted mean.
    #[must_use]
    pub const fn depolarization_ratio(&self) -> f64 {
        self.depolarization_ratio
    }

    /// Prepared `K` in `(n²-1)/(n²+2) = rho K`, in m³/kg.
    /// Composition and all original optical references are retained in the
    /// state identity; no measured wavelength dispersion is inferred.
    #[must_use]
    pub const fn mass_refractivity_m3_kg(&self) -> f64 {
        self.mass_refractivity_m3_kg
    }

    /// Coherent phase index of this exact accepted thermodynamic state.
    pub fn refractive_index(
        &self,
        wavelength_nm: f64,
    ) -> Result<GasRefractiveIndex, GasScatteringError> {
        self.refractive_index_at_density(self.gas.density, wavelength_nm)
    }

    /// Coherent phase index and exact density derivative for a local density
    /// [kg/m³] with the SAME frozen molecular composition and optical sources.
    /// This does not accept a new thermodynamic state or extrapolate a measured
    /// temperature law. A smooth density field may use the derivative to form
    /// `grad n = (dn/d rho) grad rho`; cell averages provide no such gradient.
    ///
    /// Constant molecular polarizability makes phase and group indices equal
    /// within this model. The original visible-band limit and dilute-use
    /// policy `0 <= n-1 <= 0.01` apply at every queried point. Vacuum is admitted.
    pub fn refractive_index_at_density(
        &self,
        density_kg_m3: f64,
        wavelength_nm: f64,
    ) -> Result<GasRefractiveIndex, GasScatteringError> {
        wavelength(wavelength_nm)?;
        if !density_kg_m3.is_finite() || density_kg_m3 < 0.0 {
            return Err(GasScatteringError::RefractivityDomain);
        }
        let lorentz_lorenz = density_kg_m3 * self.mass_refractivity_m3_kg;
        let max_n2_minus_one = MAX_DILUTE_REFRACTIVITY * (2.0 + MAX_DILUTE_REFRACTIVITY);
        let limit = max_n2_minus_one / (3.0 + max_n2_minus_one);
        // One outward ULP admits an exactly-on-limit source after converting
        // between molecular and mass density. It is only boundary roundoff.
        if !lorentz_lorenz.is_finite() || lorentz_lorenz > limit.next_up() {
            return Err(GasScatteringError::RefractivityDomain);
        }
        if density_kg_m3 > 0.0 && self.mass_refractivity_m3_kg > 0.0 && lorentz_lorenz == 0.0 {
            return Err(GasScatteringError::Unrepresentable);
        }
        let denominator = 1.0 - lorentz_lorenz;
        let n2_minus_one = 3.0 * lorentz_lorenz / denominator;
        let phase_index = fs_math::det::sqrt(1.0 + n2_minus_one);
        // Rationalization preserves n-1 even when n rounds to exactly one.
        let refractivity = n2_minus_one / (phase_index + 1.0);
        let derivative_density_m3_kg =
            1.5 * self.mass_refractivity_m3_kg / (phase_index * denominator * denominator);
        if !derivative_density_m3_kg.is_finite() {
            return Err(GasScatteringError::Unrepresentable);
        }
        Ok(GasRefractiveIndex {
            phase_index,
            refractivity,
            derivative_density_m3_kg,
        })
    }
    /// Exact inverse-fourth visible-band mass scattering coefficient [m²/kg].
    pub fn mass_scattering_m2_kg(&self, wavelength_nm: f64) -> Result<f64, GasScatteringError> {
        inverse_fourth(
            self.mass_scattering_reference_m2_kg,
            RAYLEIGH_REFERENCE_NM,
            wavelength_nm,
        )
    }
    /// Local attenuation by molecular scattering [m^-1], at the accepted density.
    pub fn scattering_per_m(&self, wavelength_nm: f64) -> Result<f64, GasScatteringError> {
        let value = self.gas.density * self.mass_scattering_m2_kg(wavelength_nm)?;
        finite_nonnegative(value)?;
        if self.mass_scattering_reference_m2_kg > 0.0 && value == 0.0 {
            return Err(GasScatteringError::Unrepresentable);
        }
        Ok(value)
    }
    /// Normalized unpolarized phase density [sr^-1], with both vectors taken
    /// along physical travel. This is the scalar molecular law, not Stokes data.
    pub fn phase_per_sr(&self, cosine: f64) -> Result<f64, GasScatteringError> {
        if !cosine.is_finite() || !(-1.0..=1.0).contains(&cosine) {
            return Err(GasScatteringError::Angle);
        }
        let delta = self.depolarization_ratio;
        Ok(3.0 * ((1.0 + delta) + (1.0 - delta) * cosine * cosine)
            / (8.0 * core::f64::consts::PI * (2.0 + delta)))
    }
}

/// Local coherent optics of the declared dilute, constant-polarizability model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasRefractiveIndex {
    /// Dimensionless phase index; also the group index in this nondispersive model.
    pub phase_index: f64,
    /// `n - 1`, retained without subtractive cancellation for very dilute gas.
    pub refractivity: f64,
    /// Exact `dn/d rho` at this density [m³/kg].
    pub derivative_density_m3_kg: f64,
}

/// Atomic refusal before a gas/scattering state is published.
#[derive(Clone, Debug, PartialEq)]
pub enum GasScatteringError {
    /// Invalid reference number density, refractivity, anisotropy or source identity.
    OpticalReference,
    /// Missing or misaligned gas/optical rows or unsupported composition basis.
    Components,
    /// The original thermodynamic/transport owner refused the gas state.
    Gas(MaterialError),
    /// The accepted gas did not equal the original model's complete evaluation.
    StateMismatch,
    /// Wavelength is not finite or is outside 380–780 nm.
    Wavelength,
    /// Scattering cosine is outside [-1, 1] or nonfinite.
    Angle,
    /// Invalid density or local refractivity outside the explicit dilute domain.
    RefractivityDomain,
    /// A physical coefficient is not representable in finite f64 arithmetic.
    Unrepresentable,
}

impl core::fmt::Display for GasScatteringError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "molecular Rayleigh scattering refused: {self:?}")
    }
}
impl std::error::Error for GasScatteringError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Gas(error) => Some(error),
            _ => None,
        }
    }
}

fn wavelength(value: f64) -> Result<(), GasScatteringError> {
    if value.is_finite() && (380.0..=780.0).contains(&value) {
        Ok(())
    } else {
        Err(GasScatteringError::Wavelength)
    }
}
fn finite_nonnegative(value: f64) -> Result<(), GasScatteringError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(GasScatteringError::Unrepresentable)
    }
}
fn inverse_fourth(
    reference: f64,
    reference_nm: f64,
    wavelength_nm: f64,
) -> Result<f64, GasScatteringError> {
    wavelength(wavelength_nm)?;
    let ratio = reference_nm / wavelength_nm;
    let square = ratio * ratio;
    let value = reference * square * square;
    finite_nonnegative(value)?;
    if reference > 0.0 && value == 0.0 {
        return Err(GasScatteringError::Unrepresentable);
    }
    Ok(value)
}

#[cfg(test)]
#[path = "gas_scattering/tests.rs"]
mod tests;
