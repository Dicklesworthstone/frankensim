//! Constitutive adapters for the shared conservative thermal timeline.
//!
//! The material remains the only owner of the enthalpy/phase relation. This
//! implementation makes the same immutable chart available to electrical,
//! thermal and mechanical consumers without a dependency on a FEM solver.

/// Shared material/phase conditions for conduction, electricity and radiation.
pub mod source;
/// Ordered two-surface thermal conditions shared by interface physics.
pub mod interface;

use crate::phase::EquilibriumEnthalpyPhaseCurve;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

impl ThermalLaw for EquilibriumEnthalpyPhaseCurve {
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
        // The chart is piecewise linear in h; plateaus contribute zero. Unlike
        // a local derivative query this remains valid across all chart knots.
        self.knots()
            .windows(2)
            .map(|pair| {
                (pair[1].temperature_k - pair[0].temperature_k)
                    / (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
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
