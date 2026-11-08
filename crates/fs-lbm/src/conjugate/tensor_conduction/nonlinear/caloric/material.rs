//! Original caloric sources for the native full-tensor storage coordinate.
//!
//! Capacity assemblies retain their fixed amount and mass normalization.
//! Pure-solid views retain their original restricted phase domain. Neither
//! adapter prepares a replacement caloric table or invents a conductivity.

use super::TensorCaloricLaw;
use fs_conduction::enthalpy::{DepositionError, ThermalLaw};
use fs_material::thermal::PureSolidThermalLaw;
use fs_material::thermal::transport::storage::HeatStorageThermalLaw;

impl TensorCaloricLaw for HeatStorageThermalLaw<'_> {
    fn specific_energy_bounds_j_kg(&self) -> [f64; 2] {
        self.specific_storage_span_j_kg()
    }

    fn temperature_slope_kg_k_j(&self, energy: f64) -> Result<f64, DepositionError> {
        Ok(self.temperature_and_slope_at_specific_energy(energy)?[1])
    }

    fn temperature_and_slope_at_specific_energy(
        &self,
        energy: f64,
    ) -> Result<[f64; 2], DepositionError> {
        // One original inverse per query; C_total includes every constituent
        // amount already. The reference mass converts its derivative to J/kg.
        let temperature = self.temperature_k(energy)?;
        let capacity = self
            .storage()
            .heat_capacity_j_k(temperature)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let slope = self.reference_mass_kg() / capacity;
        if !slope.is_finite() || slope <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok([temperature, slope])
    }

    fn validate_reference_cell(
        &self,
        reference_mass_kg: f64,
        reference_volume_m3: f64,
    ) -> Result<(), DepositionError> {
        let original_mass = self.reference_mass_kg();
        if !reference_mass_kg.is_finite()
            || reference_mass_kg <= 0.0
            || !reference_volume_m3.is_finite()
            || reference_volume_m3 <= 0.0
            || (reference_mass_kg - original_mass).abs()
                > 64.0 * f64::EPSILON * reference_mass_kg.max(original_mass)
        {
            return Err(DepositionError::InvalidInput(
                "tensor cell mass differs from its original heat-storage assembly",
            ));
        }
        // The source can contain mass and reference-volume capacities. Its
        // supplied amounts remain authoritative; no density is inferred.
        Ok(())
    }
}

impl TensorCaloricLaw for PureSolidThermalLaw<'_> {
    fn specific_energy_bounds_j_kg(&self) -> [f64; 2] {
        self.enthalpy_range_j_kg()
    }

    fn temperature_slope_kg_k_j(&self, energy: f64) -> Result<f64, DepositionError> {
        Ok(self.temperature_and_slope_at_specific_energy(energy)?[1])
    }

    fn temperature_and_slope_at_specific_energy(
        &self,
        energy: f64,
    ) -> Result<[f64; 2], DepositionError> {
        let temperature = self.temperature_k(energy)?;
        let curve = self.curve();
        let upper = self.enthalpy_range_j_kg()[1];
        let slope = if energy == upper {
            // The full chart takes the right derivative at an interior knot.
            // At this restricted solid endpoint that could be the excluded
            // latent branch. Use the exact left derivative of the last
            // admitted source interval, including retained affine Cp(T).
            let knots = curve.knots();
            let end = knots.partition_point(|k| k.specific_enthalpy_j_kg <= upper);
            let index = end.checked_sub(2).ok_or(DepositionError::MaterialDomain)?;
            let a = knots[index];
            let b = knots[index + 1];
            let difference = b.temperature_k - a.temperature_k;
            let slope = curve
                .sensible_heat_capacities_on_interval(index)
                .map_or_else(
                    || difference / (b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg),
                    |[_, hot_capacity]| hot_capacity.recip(),
                );
            if difference > 0.0 && slope == 0.0 {
                return Err(DepositionError::MaterialDomain);
            }
            slope
        } else {
            curve
                .temperature_derivative_at_specific_enthalpy(energy)
                .map_err(|_| DepositionError::MaterialDomain)?
        };
        if !slope.is_finite() || slope < 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok([temperature, slope])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_blake3::ContentHash;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
        UncertaintyModel,
    };
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    use fs_material::thermal::transport::storage::{
        HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage,
    };
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_qty::{Dims, QuantitySpec, Temperature};

    #[test]
    fn solid_endpoint_keeps_the_admitted_sensible_derivative() {
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x43; 32]),
            [
                (0.0, 300.0, 0.0),
                (100.0, 400.0, 0.0),
                (200.0, 400.0, 1.0),
                (300.0, 500.0, 1.0),
            ]
            .map(|(h, t, fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: fraction,
                bulk_density_kg_m3: 1.0,
            })
            .to_vec(),
        )
        .unwrap();
        let solid = curve.pure_solid_thermal_law().unwrap();
        assert_eq!(solid.specific_energy_bounds_j_kg(), [0.0, 100.0]);
        assert_eq!(
            curve
                .temperature_derivative_at_specific_enthalpy(100.0)
                .unwrap(),
            0.0
        );
        assert_eq!(
            solid
                .temperature_and_slope_at_specific_energy(100.0)
                .unwrap(),
            [400.0, 1.0]
        );
        assert!(solid.temperature_slope_kg_k_j(100.0_f64.next_up()).is_err());
        assert!(solid.temperature_slope_kg_k_j(-1.0).is_err());
    }

    fn original_capacity_assembly() -> PreparedHeatStorage {
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let dims = Dims([2, 0, -2, -1, 0, 0]);
        let capacity = PropertyKey::new("sensible-capacity", dims);
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(PropertyClaim {
                key: capacity.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, 2.0), (320.0, 6.0)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    300.0,
                    320.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured affine-capacity integration fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic body".into(),
                phase: "solid".into(),
                process: "fixed amount".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let lower = QueryPoint::new()
            .with_quantity("T", temperature, 300.0)
            .unwrap();
        let upper = QueryPoint::new()
            .with_quantity("T", temperature, 320.0)
            .unwrap();
        PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
            id: 1,
            card: &card,
            capacity: &capacity,
            amount: HeatStorageAmount::MassKg(2.0),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }])
        .unwrap()
    }

    #[test]
    fn assembly_uses_original_capacity_and_preserves_its_mass_binding() {
        let storage = original_capacity_assembly();
        let law = HeatStorageThermalLaw::try_new(&storage, 2.0, 300.0).unwrap();
        assert_eq!(law.specific_energy_bounds_j_kg(), [0.0, 80.0]);
        // Original c(T)=2+0.2*(T-300), hence e=2*x+0.1*x².
        for (temperature, energy) in [(300.0, 0.0), (310.0, 30.0), (320.0, 80.0)] {
            let [actual, slope] = law
                .temperature_and_slope_at_specific_energy(energy)
                .unwrap();
            assert!((actual - temperature).abs() < 2e-12);
            assert!((slope - 1.0 / (2.0 + 0.2 * (temperature - 300.0))).abs() < 2e-14);
        }
        assert!(law.validate_reference_cell(2.0, 0.1).is_ok());
        assert!(law.validate_reference_cell(1.0, 0.1).is_err());
        assert!(law.validate_reference_cell(4.0, 0.1).is_err());
        assert!(law.validate_reference_cell(2.0, 0.0).is_err());
        assert!(law.temperature_slope_kg_k_j(81.0).is_err());
    }
}
