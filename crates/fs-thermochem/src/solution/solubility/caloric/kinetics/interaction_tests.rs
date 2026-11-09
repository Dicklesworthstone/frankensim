use super::*;
use super::tests::{component, controls, model, near, state};
use crate::solution::caloric::MolarHeatCapacityCurveV1;
use crate::solution::solubility::IdealSolubilityV1;
use crate::solution::solubility::kinetics::PrecipitationPolicyV1;
use crate::{MolarThermalQuantityV1, SpeciesId};
use fs_qty::AmountConcentration;

// The original independent pair law is itself a consistent zero-excess source.
// Selecting it exercises the shared frozen-phase/heat split against a known
// thermodynamic limit without a second numerical implementation.
#[derive(Debug)]
struct OriginalPair {
    law: IdealSolubilityV1,
    domain: [Temperature; 2],
}
impl SolubilityThermalInteractionV1 for OriginalPair {
    fn temperature_domain(&self) -> [Temperature; 2] {
        self.domain
    }
    fn admit_binding(
        &self,
        caloric: &SolubilityCaloricModelV1,
        kinetics: &[SolubilityKineticsV1],
    ) -> Result<(), SolubilityCaloricErrorV1> {
        if kinetics.len() != 1
            || caloric.components().len() != 1
            || kinetics[0].law() != &self.law
            || caloric.components()[0].law() != &self.law
        {
            return Err(SolubilityCaloricErrorV1::Invalid("original test pair binding"));
        }
        Ok(())
    }
    fn conditional_saturation(
        &self,
        temperature: Temperature,
        _: Volume,
        _: &[SolubilityKineticStateV1],
        component: usize,
    ) -> Result<AmountConcentration, SolubilityCaloricErrorV1> {
        assert_eq!(component, 0);
        Ok(self.law.saturation_concentration(temperature)?)
    }
    fn enthalpy_increment_j(
        &self,
        _: Volume,
        _: Temperature,
        _: &[SolubilityKineticStateV1],
        _: Temperature,
        _: &[SolubilityKineticStateV1],
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        Ok(0.0)
    }
}

#[test]
fn interaction_split_uses_native_phase_trajectories_and_closes_original_heat() {
    let base = model(
        vec![component("A", 600.0, 3000.0, 30.0, 30.0)],
        &[1e-4],
        PrecipitationPolicyV1::ExistingSolidOnly,
    );
    let source: Arc<dyn SolubilityThermalInteractionV1> = Arc::new(OriginalPair {
        law: base.kinetics()[0].law().clone(),
        domain: [Temperature::new(280.0), Temperature::new(340.0)],
    });
    let bound = base.with_interactions(Arc::clone(&source)).unwrap();
    let initial = state(&bound, &[(1.0, 0.1)]);
    let mut ctl = controls(0.1);
    // One accepted comparison interval: reconstruct its TWO native half stages
    // and independent constant-Cp first-law temperature exactly.
    ctl.relative_tolerance = 1.0;
    let result = bound.advance(&initial, Time::new(0.1), 0.0, ctl, || false).unwrap();
    assert_eq!(result.accepted_steps(), 1);
    let capacity = 50.0 * 50.0 + 30.0;
    let first = base.kinetics()[0].advance(
        Temperature::new(300.0), Volume::new(1e-3), &initial.phases()[0], Time::new(0.05),
    ).unwrap();
    let first_t = 300.0 - 3000.0 * first.transferred_to_solution().value() / capacity;
    let second = base.kinetics()[0].advance(
        Temperature::new(first_t), Volume::new(1e-3), first.state(), Time::new(0.05),
    ).unwrap();
    let expected_t = 300.0 - 3000.0 * (second.state().dissolved().value() - 0.1) / capacity;
    near(result.final_state().temperature().value(), expected_t, 2e-10);
    near(result.final_state().phases()[0].dissolved().value(), second.state().dissolved().value(), 1e-13);
    near(bound.enthalpy_change_j(&initial, result.final_state()).unwrap(), 0.0, 1e-7);
    assert_eq!(result.final_state().phases()[0].total_solute(), initial.phases()[0].total_solute());
    assert!(Arc::ptr_eq(bound.interaction_source().unwrap(), &source));
    let snapshot = initial.clone();
    assert_eq!(bound.advance(&initial, Time::new(1.0), 0.0, controls(0.1), || true), Err(SolubilityCaloricErrorV1::Cancelled));
    assert_eq!(initial, snapshot);
}

#[test]
fn interaction_source_identity_and_temperature_support_remain_binding() {
    let base = model(
        vec![component("A", 600.0, 0.0, 30.0, 30.0)],
        &[1e-4],
        PrecipitationPolicyV1::ExistingSolidOnly,
    );
    let create = || -> Arc<dyn SolubilityThermalInteractionV1> { Arc::new(OriginalPair {
        law: base.kinetics()[0].law().clone(),
        domain: [Temperature::new(290.0), Temperature::new(310.0)],
    }) };
    let source = create();
    let bound = base.with_interactions(Arc::clone(&source)).unwrap();
    assert_eq!(bound.with_interactions(source).unwrap(), bound);
    assert!(bound.with_interactions(create()).is_err());
    let initial = state(&bound, &[(1.0, 0.1)]);
    assert!(base.advance(&initial, Time::new(1.0), 0.0, controls(0.1), || false).is_err());
    // The default broader root bracket may not extrapolate the narrower source.
    assert!(bound.advance(&initial, Time::new(1.0), 0.0, controls(0.1), || false).is_err());
    let mut ctl = controls(0.1);
    ctl.root.temperature_bracket = bound.temperature_domain();
    let result = bound.advance(&initial, Time::new(1.0), 0.0, ctl, || false).unwrap();
    let exact = base.kinetics()[0].advance(
        Temperature::new(300.0), Volume::new(1e-3), &initial.phases()[0], Time::new(1.0),
    ).unwrap();
    near(result.final_state().phases()[0].dissolved().value(), exact.state().dissolved().value(), 1e-14);
    assert_eq!(result.final_state().temperature(), initial.temperature());
}

#[test]
fn heat_capacity_bounds_include_interior_source_extrema_without_extrapolation() {
    let curve = MolarHeatCapacityCurveV1::new(
        SpeciesId::new("bounded-species").unwrap(), Temperature::new(300.0),
        [(280.0, 30.0), (300.0, 90.0), (310.0, 20.0), (340.0, 50.0)]
            .map(|(t, cp)| (Temperature::new(t), MolarThermalQuantityV1::new(cp))).to_vec(),
    ).unwrap();
    let bounds = |lo, hi| curve.heat_capacity_bounds([Temperature::new(lo), Temperature::new(hi)]).unwrap().map(|value| value.value());
    assert_eq!(bounds(290.0, 330.0), [20.0, 90.0]);
    assert_eq!(bounds(280.0, 290.0), [30.0, 60.0]);
    assert_eq!(bounds(305.0, 305.0), [55.0, 55.0]);
    assert!(curve.heat_capacity_bounds([Temperature::new(279.0), Temperature::new(300.0)]).is_err());
    assert!(curve.heat_capacity_bounds([Temperature::new(310.0), Temperature::new(300.0)]).is_err());
}
