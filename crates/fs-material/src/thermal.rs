//! Constitutive adapters for the shared conservative thermal timeline.
//!
//! The material remains the only owner of the enthalpy/phase relation. This
//! implementation makes the same immutable chart available to electrical,
//! thermal and mechanical consumers without a dependency on a FEM solver.

/// Shared material/phase conditions for conduction, electricity and radiation.
pub mod source;
/// Ordered two-surface thermal conditions shared by interface physics.
pub mod interface;
/// Prepared source-bound density, heat capacity and directional heat transport.
pub mod transport;
/// Original caloric charts restricted to their actual all-solid phase interval.
pub mod solid;
pub use solid::PureSolidThermalLaw;
/// Integrated temperature-dependent conduction through finite planar layers.
pub mod conduction;

use crate::phase::EquilibriumEnthalpyPhaseCurve;
use crate::phase::vaporization::IsobaricWaterVaporization;
use crate::phase::water::IsobaricWater;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

impl ThermalLaw for IsobaricWater {
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.state_at_specific_enthalpy(h, || false)
            .map(|answer| Some(answer.state().bulk_density_kg_m3()))
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.state_at_specific_enthalpy(specific_h_j_kg, || false)
            .map(|answer| answer.state().temperature_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        IsobaricWater::maximum_temperature_slope(self)
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        IsobaricWater::specific_entropy_change_j_kg_k(self, initial_h, final_h, || false)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

// The external constant-pressure load permits expansion/contraction. The
// thermal kernel's fixed geometry is reference heat-exchange geometry, not a
// rigid vessel or an automatically updated vapor volume.
impl ThermalLaw for IsobaricWaterVaporization {
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.state_at_specific_enthalpy(h)
            .map(|state| Some(state.bulk_density_kg_m3()))
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.state_at_specific_enthalpy(specific_h_j_kg)
            .map(|state| state.temperature_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 { 0.0 }

    fn specific_entropy_change_j_kg_k(
        &self, initial_h_j_kg: f64, final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        // Exact integral on the represented h interval. Original-source
        // absolute entropy endpoints can differ by source roundoff.
        self.state_at_specific_enthalpy(initial_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let final_state = self.state_at_specific_enthalpy(final_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(Some((final_h_j_kg - initial_h_j_kg) / final_state.temperature_k()))
    }
}

impl ThermalLaw for EquilibriumEnthalpyPhaseCurve {
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.values_at_specific_enthalpy(h)
            .map(|state| Some(state.bulk_density_kg_m3))
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity().as_bytes())
    }

    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.temperature_at_specific_enthalpy(specific_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn specific_entropy_change_j_kg_k(
        &self, initial_h_j_kg: f64, final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        EquilibriumEnthalpyPhaseCurve::specific_entropy_change_j_kg_k(
            self, initial_h_j_kg, final_h_j_kg)
            .map(Some).map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        // Original affine Cp is positive throughout each sensible interval,
        // so 1/min(Cp endpoints) bounds its nonlinear T(h) tangent. Authored
        // enthalpy chords and latent plateaus keep their original slopes.
        self.knots()
            .windows(2)
            .enumerate()
            .map(|(index, pair)| {
                self.sensible_heat_capacities_on_interval(index).map_or_else(
                    || (pair[1].temperature_k - pair[0].temperature_k)
                        / (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg),
                    |[a, b]| 1.0 / a.min(b),
                )
            })
            .fold(0.0, f64::max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phase::EnthalpyPhaseKnot;
    use fs_blake3::ContentHash;

    #[test]
    fn original_cp_bounds_whole_solid_tangent_and_temperature_inverse() {
        use crate::phase::{HeatCapacityKnot, SolidLiquidPhase};
        let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
            ContentHash([0x61; 32]), SolidLiquidPhase::Solid, 0.0,
            &[
                HeatCapacityKnot { temperature_k: 300.0, specific_heat_capacity_j_kg_k: 600.0, bulk_density_kg_m3: 1000.0 },
                HeatCapacityKnot { temperature_k: 600.0, specific_heat_capacity_j_kg_k: 1200.0, bulk_density_kg_m3: 900.0 },
            ], 1e6,
        ).unwrap();
        let solid = curve.pure_solid_thermal_law().unwrap();
        assert_eq!(ThermalLaw::maximum_temperature_slope(&curve), 1.0 / 600.0);
        assert_eq!(solid.maximum_temperature_slope(), 1.0 / 600.0);
        for temperature in [300.0, 350.25, 479.5, 600.0] {
            let h = solid.specific_enthalpy_at_temperature(temperature).unwrap();
            let expected = (temperature - 300.0) * (temperature + 300.0);
            assert!((h - expected).abs() < 1e-9);
            assert!((solid.temperature_k(h).unwrap() - temperature).abs() < 2e-12);
        }
    }

    #[test]
    fn phase_entropy_adapter_retains_sensible_and_latent_storage() {
        let knot = |h, t, fraction| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h, temperature_k: t,
            liquid_mass_fraction: fraction, bulk_density_kg_m3: 1000.0,
        };
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x76; 32]),
            vec![knot(0.0,300.0,0.0), knot(30_000.0,600.0,0.0),
                knot(55_000.0,600.0,1.0), knot(95_000.0,800.0,1.0)]).unwrap();
        let entropy = |a, b| ThermalLaw::specific_entropy_change_j_kg_k(&curve,a,b)
            .unwrap().unwrap();
        let exact = 100.0*2.0_f64.ln() + 25_000.0/600.0
            + 200.0*(800.0_f64/600.0).ln();
        assert!((entropy(0.0,95_000.0)-exact).abs() < 1e-12);
        assert_eq!(entropy(40_000.0,50_000.0),10_000.0/600.0);
        assert_eq!(entropy(95_000.0,0.0),-entropy(0.0,95_000.0));
        assert_eq!(entropy(40_000.0,40_000.0),0.0);
        assert!(ThermalLaw::specific_entropy_change_j_kg_k(&curve,-1.0,0.0).is_err());
    }
}
