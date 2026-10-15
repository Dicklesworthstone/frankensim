//! Named, near-ambient gas data over the existing ideal-gas/Wilke/WMS owners.
//!
//! These are explicit engineering models, NOT source-resolved material cards
//! or accuracy certificates. Heat capacities are frozen at 298.15 K. The
//! admitted window [273.15, 323.15] K and [10000, 200000] Pa deliberately excludes
//! cryogenic/high-pressure extrapolation. It is a model-use policy, not a
//! measured error bound. Water partial pressure is checked by the existing
//! liquid-water saturation owner; reactions and evolving composition are absent.
//!
//! Molar masses follow the committed `data/matdb/seed-v1` species data.
//! Molar cp values for N2, O2 and CO2 are NIST-JANAF's 298.15 K rows:
//! <https://janaf.nist.gov/tables/N-023.html>,
//! <https://janaf.nist.gov/tables/O-029.html>,
//! <https://janaf.nist.gov/tables/C-095.html>. Argon is calorically perfect
//! monatomic gas. Sutherland reference values are Table 5-2 in
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.46.html>.
//! Conductivity uses GasSpec's existing Eucken approximation, NOT that source's
//! separate conductivity fit. H2O reuses the existing water-vapor GasSpec.
//! No arbitrary mixture is thereby eligible for the air-only ISO 9613 law.

use crate::gas::{ConductivityModel, GasSpec, GasState, R_USSA_1976, saturation_pressure_water_pa};
use crate::MaterialError;
use fs_qty::semantic::Composition;
use std::str::FromStr;

/// Explicit species identities supported by the bounded convenience model.
/// Dry air is intentionally not a species: its legacy pseudocomponent model
/// remains available through `GasState::try_new_moist_air`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GasSpecies {
    /// Molecular nitrogen, N2.
    Nitrogen,
    /// Molecular oxygen, O2.
    Oxygen,
    /// Monatomic argon, Ar.
    Argon,
    /// Carbon dioxide, CO2.
    CarbonDioxide,
    /// Water vapor, H2O; not liquid water or ice.
    WaterVapor,
}

impl GasSpecies {
    /// Complete supported set; also bounds cold-path input allocation.
    pub const ALL: [Self; 5] = [Self::Nitrogen, Self::Oxygen, Self::Argon,
        Self::CarbonDioxide, Self::WaterVapor];

    /// Case-sensitive, canonical chemical symbol for configuration round trips.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::Nitrogen => "N2", Self::Oxygen => "O2", Self::Argon => "Ar",
            Self::CarbonDioxide => "CO2", Self::WaterVapor => "H2O",
        }
    }

    /// Declared approximate coefficients, not an evaluated or phase-admitted
    /// state. Use [`resolve_gas_mixture`] for bounded named-species admission.
    #[must_use]
    pub fn spec(self) -> GasSpec {
        let (molar_mass, gamma, mu0, sutherland_s) = match self {
            Self::Nitrogen => (28.01340e-3, 29.124 / (29.124 - R_USSA_1976), 1.663e-5, 107.0),
            Self::Oxygen => (31.99880e-3, 29.376 / (29.376 - R_USSA_1976), 1.919e-5, 139.0),
            Self::Argon => (39.948e-3, 5.0 / 3.0, 2.125e-5, 114.0),
            Self::CarbonDioxide => (44.00950e-3, 37.129 / (37.129 - R_USSA_1976), 1.370e-5, 222.0),
            Self::WaterVapor => return GasSpec::water_vapor_nist(),
        };
        let t0: f64 = 273.0;
        GasSpec {
            molar_mass, gamma, sutherland_s,
            // Convert the source's reference parameterization to GasSpec's
            // beta. Evaluation remains solely in GasState::try_new.
            sutherland_beta: mu0 * (t0 + sutherland_s) / (t0 * t0.sqrt()),
            conductivity: ConductivityModel::Eucken,
        }
    }
}

impl FromStr for GasSpecies {
    type Err = MaterialError;
    fn from_str(name: &str) -> Result<Self, Self::Err> {
        match name {
            "N2" | "nitrogen" => Ok(Self::Nitrogen),
            "O2" | "oxygen" => Ok(Self::Oxygen),
            "Ar" | "argon" => Ok(Self::Argon),
            "CO2" | "carbon-dioxide" => Ok(Self::CarbonDioxide),
            "H2O" | "water-vapor" => Ok(Self::WaterVapor),
            _ => Err(invalid("unsupported gas species; expected N2, O2, Ar, CO2 or H2O")),
        }
    }
}

fn invalid(what: &'static str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}

/// Resolve ordered, distinct named species at one temperature [K] and total
/// static pressure [Pa], returning the SAME medium used throughout physics.
/// Composition uses the shared mass/mole-fraction authority, never percentages,
/// balance-gas inference or parser-side normalization. Omit zero components.
///
/// All thermodynamic and transport evaluation delegates to
/// [`GasState::try_new_mixture`]. This wrapper selects declared data and checks
/// the named model's domain and water phase; it does not mix properties itself.
/// Legacy dry/humid-air evaluation and the general supplied-GasSpec API retain
/// their existing, separately documented domains and arithmetic.
///
/// # Errors
/// Empty/oversized/duplicate species, unsupported composition, invalid fractions,
/// temperatures/pressures outside the named model window, supersaturated water,
/// or any refusal of the shared material/composition owners.
pub fn resolve_gas_mixture(
    species: &[GasSpecies], composition: &Composition, temperature_k: f64, pressure_pa: f64,
) -> Result<GasState, MaterialError> {
    if species.is_empty() || species.len() > GasSpecies::ALL.len() {
        return Err(invalid("named gas mixture requires 1..=5 distinct supported species"));
    }
    for (i, component) in species.iter().enumerate() {
        if species[..i].contains(component) {
            return Err(invalid("duplicate gas species; supply each species exactly once"));
        }
    }
    if !(273.15..=323.15).contains(&temperature_k) {
        return Err(invalid("named gas model requires temperature in [273.15, 323.15] K"));
    }
    if !(10000.0..=200000.0).contains(&pressure_pa) {
        return Err(invalid("named gas model requires total static pressure in [10000, 200000] Pa"));
    }
    let specs: Vec<_> = species.iter().map(|species| species.spec()).collect();
    let water = species.iter().position(|species| *species == GasSpecies::WaterVapor);
    let state = GasState::try_new_mixture(&specs, composition, temperature_k, pressure_pa, water)?;
    if water.is_some() && state.water_mole_fraction * pressure_pa > saturation_pressure_water_pa(temperature_k)? {
        return Err(invalid("water partial pressure exceeds saturation; condensation is not modeled"));
    }
    Ok(state)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::{MolarMass, semantic::CompositionBasis};
    use super::GasSpecies::{Argon, CarbonDioxide, Nitrogen, Oxygen, WaterVapor};

    fn mole(fractions: Vec<f64>) -> Composition {
        Composition::new(CompositionBasis::MoleFraction, fractions).unwrap()
    }

    #[test]
    fn names_and_reference_viscosities_are_explicit() {
        for (species, mu0) in [(Nitrogen, 1.663e-5), (Oxygen, 1.919e-5),
            (Argon, 2.125e-5), (CarbonDioxide, 1.370e-5)] {
            assert_eq!(species.symbol().parse::<GasSpecies>().unwrap(), species);
            let state = GasState::try_new(&species.spec(), 273.0, 100000.0).unwrap();
            assert!((state.dynamic_viscosity / mu0 - 1.0).abs() < 1e-14);
        }
        assert_eq!("water-vapor".parse::<GasSpecies>().unwrap(), WaterVapor);
        assert_eq!(WaterVapor.spec(), GasSpec::water_vapor_nist());
        for name in ["air", "He", "N", "o2", "unknown"] {
            assert!(name.parse::<GasSpecies>().is_err(), "accepted {name}");
        }
    }

    #[test]
    fn named_resolution_is_the_existing_mixture_owner_bit_for_bit() {
        let species = [Nitrogen, Oxygen, Argon, CarbonDioxide, WaterVapor];
        let composition = mole(vec![0.75, 0.20, 0.03, 0.01, 0.01]);
        let specs: Vec<_> = species.iter().map(|s| s.spec()).collect();
        let expected = GasState::try_new_mixture(&specs, &composition, 300.0, 100000.0, Some(4)).unwrap();
        let actual = resolve_gas_mixture(&species, &composition, 300.0, 100000.0).unwrap();
        assert_eq!(actual, expected);
        actual.validate_coherence().unwrap();
    }

    #[test]
    fn mass_and_mole_inputs_describe_the_same_medium() {
        let species = [Nitrogen, Oxygen, CarbonDioxide];
        let mole = mole(vec![0.6, 0.3, 0.1]);
        let masses: Vec<_> = species.iter().map(|s| MolarMass::new(s.spec().molar_mass)).collect();
        let mass = mole.to_mass_fractions(&masses).unwrap();
        let a = resolve_gas_mixture(&species, &mole, 293.15, 100000.0).unwrap();
        let b = resolve_gas_mixture(&species, &mass, 293.15, 100000.0).unwrap();
        for (x, y) in [(a.density, b.density), (a.sound_speed, b.sound_speed),
            (a.specific_heat_cp, b.specific_heat_cp), (a.dynamic_viscosity, b.dynamic_viscosity),
            (a.thermal_conductivity, b.thermal_conductivity)] {
            assert!((x / y - 1.0).abs() < 1e-13);
        }
    }

    #[test]
    fn pressure_changes_density_not_ideal_gas_sound_speed() {
        let composition = mole(vec![0.7, 0.3]);
        for temperature in [273.15, 298.15, 323.15] {
            let a = resolve_gas_mixture(&[Nitrogen, Oxygen], &composition, temperature, 50000.0).unwrap();
            let b = resolve_gas_mixture(&[Nitrogen, Oxygen], &composition, temperature, 100000.0).unwrap();
            assert_eq!(2.0 * a.density, b.density);
            assert_eq!(2.0 * a.characteristic_impedance, b.characteristic_impedance);
            assert_eq!(a.sound_speed, b.sound_speed);
            assert_eq!(a.dynamic_viscosity, b.dynamic_viscosity);
            assert_eq!(a.thermal_conductivity, b.thermal_conductivity);
        }
        let n = resolve_gas_mixture(&[Nitrogen], &mole(vec![1.0]), 298.15, 100000.0).unwrap();
        let c = resolve_gas_mixture(&[CarbonDioxide], &mole(vec![1.0]), 298.15, 100000.0).unwrap();
        assert!(c.sound_speed < n.sound_speed && c.density > n.density);
        assert_ne!(c.prandtl, n.prandtl);
    }

    #[test]
    fn phase_and_domain_refusals_do_not_repair_inputs() {
        let one = mole(vec![1.0]);
        for (t, p) in [(273.14, 1e5), (323.16, 1e5), (f64::NAN, 1e5),
            (300.0, 0.0), (300.0, 9999.0), (300.0, 200001.0), (300.0, f64::INFINITY)] {
            assert!(resolve_gas_mixture(&[Nitrogen], &one, t, p).is_err());
        }
        assert!(resolve_gas_mixture(&[], &one, 300.0, 1e5).is_err());
        assert!(resolve_gas_mixture(&[Nitrogen, Nitrogen], &mole(vec![0.5, 0.5]), 300.0, 1e5).is_err());
        assert!(resolve_gas_mixture(&[Nitrogen, Oxygen], &mole(vec![1.0, 0.0]), 300.0, 1e5).is_err());
        assert!(resolve_gas_mixture(&[Nitrogen, Oxygen], &one, 300.0, 1e5).is_err());
        let wet = mole(vec![0.98, 0.02]);
        assert!(resolve_gas_mixture(&[Nitrogen, WaterVapor], &wet, 293.15, 100000.0).is_ok());
        assert!(resolve_gas_mixture(&[Nitrogen, WaterVapor], &wet, 293.15, 200000.0).is_err());
        assert!(resolve_gas_mixture(&[WaterVapor], &one, 293.15, 100000.0).is_err());
        let steam = resolve_gas_mixture(&[WaterVapor], &one, 323.15, 10000.0).unwrap();
        assert_eq!(steam.water_mole_fraction, 1.0);
    }
}
