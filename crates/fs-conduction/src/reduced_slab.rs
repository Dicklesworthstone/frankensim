//! First-mode one-dimensional slab cooling screen.
//!
//! This deliberately small reduced model evaluates
//! `T(t) = T_boundary + (T_initial - T_boundary) exp(-t/τ)` with
//! `τ = thickness²/(π² α)`. It admits only monotone cooling to a fixed,
//! uniform boundary temperature. It is not the transient finite-element rung,
//! does not infer convection coefficients or contact resistance, and does not
//! call a glass-transition crossing "solidification".

use core::f64::consts::PI;
use fs_matdb::{
    EvaluationDecision, MaterialCard, PropertyKey, PropertyUsageReceipt, PropertyValue,
    SelectionPolicy,
};

use crate::ConductionError;
use crate::transient::MaterialHeatCapacity;

/// Geometry and temperatures for a source-bound, constant-property screen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaterialFirstModeSlabInput {
    /// Full slab thickness [m]; conductivity is along this heat-flow direction.
    pub thickness_m: f64,
    /// Initial temperature [K].
    pub initial_temperature_k: f64,
    /// Fixed boundary temperature [K].
    pub boundary_temperature_k: f64,
    /// Strictly interior cooling threshold [K].
    pub threshold_temperature_k: f64,
}

/// Source-bound slab screening refusal, retaining the original material error.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialFirstModeSlabError {
    /// Same-card, source-coordinate, quantity, or temperature-domain refusal.
    Material(ConductionError),
    /// A curve or nonconstant source cannot supply a constant diffusivity.
    NonconstantConductivity,
    /// Invalid geometry, temperatures, elapsed time or arithmetic.
    Screen(FirstModeSlabError),
}

/// The existing exponential screen with diffusivity derived from one material.
///
/// The full cooling interval uses the original, constant scalar rho, Cp and k
/// sources. This is a first-mode screening approximation, not an exact transient
/// slab solution or an extension of any source's temperature support.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialFirstModeSlabCooling {
    source: MaterialHeatCapacity,
    conductivity_receipts: Vec<PropertyUsageReceipt>,
    input: FirstModeSlabInput,
    step: FirstModeSlabStep,
}

impl MaterialFirstModeSlabCooling {
    /// Bind exact conductivity metadata to the original heat-capacity card and
    /// its unchanged external conditions. Selection pins one claim over the
    /// complete prepared interval; even a flat conductivity curve refuses.
    ///
    /// Conductivity is the component along the one-dimensional slab heat-flow
    /// direction. Choosing a directional key does not make that material isotropic.
    pub fn from_card(
        source: &MaterialHeatCapacity,
        card: &MaterialCard,
        conductivity: &PropertyKey,
        selection: SelectionPolicy,
        controls: MaterialFirstModeSlabInput,
    ) -> Result<Self, MaterialFirstModeSlabError> {
        let model = source
            .isotropic_conductivity_with_property(card, conductivity, selection)
            .map_err(MaterialFirstModeSlabError::Material)?;
        let conductivity_receipts: Vec<_> = model.receipts().into_iter().cloned().collect();
        if conductivity_receipts.is_empty()
            || conductivity_receipts.iter().any(|receipt| {
                receipt.decision != EvaluationDecision::ConstantWithinValidity
                    || !card
                        .claims()
                        .claim(receipt.selected)
                        .is_some_and(|claim| matches!(claim.value, PropertyValue::Scalar { .. }))
            })
        {
            return Err(MaterialFirstModeSlabError::NonconstantConductivity);
        }
        for temperature in [
            controls.initial_temperature_k,
            controls.boundary_temperature_k,
            controls.threshold_temperature_k,
        ] {
            source
                .temperature_span()
                .check(temperature)
                .map_err(MaterialFirstModeSlabError::Material)?;
        }
        let k = model
            .tensor_at(controls.initial_temperature_k)
            .map_err(MaterialFirstModeSlabError::Material)?[0][0];
        let input = FirstModeSlabInput {
            thickness_m: controls.thickness_m,
            thermal_diffusivity_m2_s: k / source.value_j_per_m3_k(),
            initial_temperature_k: controls.initial_temperature_k,
            boundary_temperature_k: controls.boundary_temperature_k,
            threshold_temperature_k: controls.threshold_temperature_k,
        };
        let step =
            step_first_mode_slab_cooling(input).map_err(MaterialFirstModeSlabError::Screen)?;
        if step.time_constant_s <= 0.0 || step.time_to_threshold_s <= 0.0 {
            return Err(MaterialFirstModeSlabError::Screen(
                FirstModeSlabError::NonFiniteResult,
            ));
        }
        Ok(Self {
            source: source.clone(),
            conductivity_receipts,
            input,
            step,
        })
    }

    /// Original density/Cp binding, including its reference and endpoint receipts.
    #[must_use]
    pub fn source(&self) -> &MaterialHeatCapacity {
        &self.source
    }

    /// Original conductivity selections at the prepared interval's endpoints.
    #[must_use]
    pub fn conductivity_receipts(&self) -> &[PropertyUsageReceipt] {
        &self.conductivity_receipts
    }

    /// Admitted SI inputs, including source-derived diffusivity.
    #[must_use]
    pub const fn input(&self) -> FirstModeSlabInput {
        self.input
    }

    /// Time constant and threshold crossing from the existing screen.
    #[must_use]
    pub const fn step(&self) -> FirstModeSlabStep {
        self.step
    }

    /// Evaluate using only this screen's admitted source, geometry and temperatures.
    pub fn temperature_at_s(&self, elapsed_s: f64) -> Result<f64, MaterialFirstModeSlabError> {
        self.step
            .temperature_at_s(self.input, elapsed_s)
            .map_err(MaterialFirstModeSlabError::Screen)
    }
}

/// SI inputs to the declared first-mode cooling screen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FirstModeSlabInput {
    /// Slab thickness in metres.
    pub thickness_m: f64,
    /// Thermal diffusivity in square metres per second.
    pub thermal_diffusivity_m2_s: f64,
    /// Initial uniform temperature in kelvin.
    pub initial_temperature_k: f64,
    /// Fixed boundary temperature in kelvin.
    pub boundary_temperature_k: f64,
    /// Screening threshold in kelvin.
    pub threshold_temperature_k: f64,
}

/// Admitted first-mode cooling result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FirstModeSlabStep {
    /// First-mode time constant in seconds.
    pub time_constant_s: f64,
    /// Time at which the exponential screen crosses the declared threshold.
    pub time_to_threshold_s: f64,
}

impl FirstModeSlabStep {
    /// Evaluate the declared exponential screen at non-negative time.
    pub fn temperature_at_s(
        self,
        input: FirstModeSlabInput,
        elapsed_s: f64,
    ) -> Result<f64, FirstModeSlabError> {
        if !elapsed_s.is_finite() {
            return Err(FirstModeSlabError::NonFinite { field: "elapsed_s" });
        }
        if elapsed_s < 0.0 {
            return Err(FirstModeSlabError::NegativeTime);
        }
        Ok(input.boundary_temperature_k
            + (input.initial_temperature_k - input.boundary_temperature_k)
                * (-elapsed_s / self.time_constant_s).exp())
    }
}

/// Typed refusal at the reduced thermal boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FirstModeSlabError {
    /// A named input was NaN or infinite.
    NonFinite {
        /// Stable input field name.
        field: &'static str,
    },
    /// A named geometric or material quantity was zero or negative.
    NonPositive {
        /// Stable input field name.
        field: &'static str,
    },
    /// Temperatures do not describe monotone cooling through the threshold.
    ThresholdOutsideCoolingInterval,
    /// A requested evaluation time was negative.
    NegativeTime,
    /// The admitted calculation became non-finite.
    NonFiniteResult,
}

/// Evaluate the declared fixed-boundary first-mode cooling screen.
pub fn step_first_mode_slab_cooling(
    input: FirstModeSlabInput,
) -> Result<FirstModeSlabStep, FirstModeSlabError> {
    for (field, value) in [
        ("thickness_m", input.thickness_m),
        ("thermal_diffusivity_m2_s", input.thermal_diffusivity_m2_s),
        ("initial_temperature_k", input.initial_temperature_k),
        ("boundary_temperature_k", input.boundary_temperature_k),
        ("threshold_temperature_k", input.threshold_temperature_k),
    ] {
        if !value.is_finite() {
            return Err(FirstModeSlabError::NonFinite { field });
        }
    }
    for (field, value) in [
        ("thickness_m", input.thickness_m),
        ("thermal_diffusivity_m2_s", input.thermal_diffusivity_m2_s),
        ("initial_temperature_k", input.initial_temperature_k),
        ("boundary_temperature_k", input.boundary_temperature_k),
        ("threshold_temperature_k", input.threshold_temperature_k),
    ] {
        if value <= 0.0 {
            return Err(FirstModeSlabError::NonPositive { field });
        }
    }
    if !(input.initial_temperature_k > input.threshold_temperature_k
        && input.threshold_temperature_k > input.boundary_temperature_k)
    {
        return Err(FirstModeSlabError::ThresholdOutsideCoolingInterval);
    }

    let time_constant_s = input.thickness_m.powi(2) / (PI.powi(2) * input.thermal_diffusivity_m2_s);
    let temperature_ratio = (input.initial_temperature_k - input.boundary_temperature_k)
        / (input.threshold_temperature_k - input.boundary_temperature_k);
    let time_to_threshold_s = time_constant_s * temperature_ratio.ln();
    if !time_constant_s.is_finite() || !time_to_threshold_s.is_finite() {
        return Err(FirstModeSlabError::NonFiniteResult);
    }
    Ok(FirstModeSlabStep {
        time_constant_s,
        time_to_threshold_s,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, Provenance, QueryPoint,
        UncertaintyModel,
    };
    use fs_material::state_point::{
        DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
        SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    };
    use fs_qty::{Density, Temperature};

    const INPUT: FirstModeSlabInput = FirstModeSlabInput {
        thickness_m: 0.2e-3,
        thermal_diffusivity_m2_s: 0.082e-6,
        initial_temperature_k: 498.15,
        boundary_temperature_k: 298.15,
        threshold_temperature_k: 378.15,
    };

    #[test]
    fn crosses_the_threshold_at_the_reported_time() {
        let step = step_first_mode_slab_cooling(INPUT).expect("admitted cooling screen");
        let temperature = step
            .temperature_at_s(INPUT, step.time_to_threshold_s)
            .expect("non-negative crossing time");
        assert!((temperature - INPUT.threshold_temperature_k).abs() < 1.0e-10);
        assert!(step.time_to_threshold_s < step.time_constant_s);
    }

    #[test]
    fn thickness_has_the_expected_quadratic_scaling() {
        let base = step_first_mode_slab_cooling(INPUT).expect("base screen");
        let thick = step_first_mode_slab_cooling(FirstModeSlabInput {
            thickness_m: INPUT.thickness_m * 2.0,
            ..INPUT
        })
        .expect("thick screen");
        assert!((thick.time_constant_s / base.time_constant_s - 4.0).abs() < 1.0e-12);
        assert!((thick.time_to_threshold_s / base.time_to_threshold_s - 4.0).abs() < 1.0e-12);
    }

    #[test]
    fn refuses_unordered_thresholds_and_invalid_time() {
        assert_eq!(
            step_first_mode_slab_cooling(FirstModeSlabInput {
                threshold_temperature_k: INPUT.initial_temperature_k,
                ..INPUT
            }),
            Err(FirstModeSlabError::ThresholdOutsideCoolingInterval)
        );
        let step = step_first_mode_slab_cooling(INPUT).expect("base screen");
        assert_eq!(
            step.temperature_at_s(INPUT, -0.1),
            Err(FirstModeSlabError::NegativeTime)
        );
    }

    fn source_card(k: f64, flat_curve: bool, k_span: [f64; 2]) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, value, dims) in [
            (DENSITY_PROPERTY, 800.0, Density::DIMS),
            (
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                1500.0,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            (THERMAL_CONDUCTIVITY_PROPERTY, k, THERMAL_CONDUCTIVITY_DIMS),
        ] {
            let conductivity = name == THERMAL_CONDUCTIVITY_PROPERTY;
            let curve = conductivity && flat_curve;
            let [low, high] = if conductivity { k_span } else { [280.0, 360.0] };
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: if curve {
                        PropertyValue::Curve {
                            abscissa: "T".into(),
                            abscissa_dims: Temperature::DIMS,
                            knots: vec![(low, value), (high, value)],
                            dims,
                        }
                    } else {
                        PropertyValue::Scalar { value, dims }
                    },
                    validity: ValidityDomain::unconstrained().with("T", low, high),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: if curve {
                        InterpolationPolicy::LinearInside
                    } else {
                        InterpolationPolicy::ConstantWithinValidity
                    },
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: "manufactured thermal screening coefficients; not measured data"
                            .into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic slab".into(),
                phase: "solid".into(),
                process: "numerical regression".into(),
                revision: 0,
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    fn source_capacity(card: &MaterialCard) -> MaterialHeatCapacity {
        MaterialHeatCapacity::from_card(
            card,
            &QueryPoint::new().with("T", 300.0).unwrap(),
            280.0,
            360.0,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }

    const SOURCE_INPUT: MaterialFirstModeSlabInput = MaterialFirstModeSlabInput {
        thickness_m: 0.01,
        initial_temperature_k: 350.0,
        boundary_temperature_k: 290.0,
        threshold_temperature_k: 320.0,
    };

    fn conductivity_key() -> PropertyKey {
        PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS)
    }

    #[test]
    fn g1_source_bound_cooling_uses_all_three_material_coefficients_and_retains_receipts() {
        for k in [0.2, 0.4] {
            let card = source_card(k, false, [280.0, 360.0]);
            let source = source_capacity(&card);
            let screen = MaterialFirstModeSlabCooling::from_card(
                &source,
                &card,
                &conductivity_key(),
                SelectionPolicy::SingleClaimOnly,
                SOURCE_INPUT,
            )
            .unwrap();
            let expected_tau = 0.01_f64.powi(2) * 800.0 * 1500.0 / (PI.powi(2) * k);
            assert!((screen.step().time_constant_s / expected_tau - 1.0).abs() < 1e-14);
            let crossing = screen
                .temperature_at_s(expected_tau * 2.0_f64.ln())
                .unwrap();
            assert!((crossing - 320.0).abs() < 1e-12);
            assert_eq!(screen.temperature_at_s(0.0).unwrap(), 350.0);
            assert_eq!(screen.source(), &source);
            let receipts = screen.conductivity_receipts();
            assert!(!receipts.is_empty());
            for receipt in receipts {
                assert_eq!(receipt.selected, receipts[0].selected);
                card.claims().verify_receipt(receipt).unwrap();
            }
            assert_eq!(
                screen.temperature_at_s(-1.0),
                Err(MaterialFirstModeSlabError::Screen(
                    FirstModeSlabError::NegativeTime
                ))
            );
        }
    }

    #[test]
    fn g0_source_bound_cooling_rejects_changed_card_curves_and_unsupported_temperatures() {
        let card = source_card(0.2, false, [280.0, 360.0]);
        let source = source_capacity(&card);
        let changed = source_card(0.4, false, [280.0, 360.0]);
        assert!(
            MaterialFirstModeSlabCooling::from_card(
                &source,
                &changed,
                &conductivity_key(),
                SelectionPolicy::SingleClaimOnly,
                SOURCE_INPUT
            )
            .is_err()
        );
        for changed_input in [
            MaterialFirstModeSlabInput {
                initial_temperature_k: 370.0,
                ..SOURCE_INPUT
            },
            MaterialFirstModeSlabInput {
                boundary_temperature_k: 270.0,
                ..SOURCE_INPUT
            },
            MaterialFirstModeSlabInput {
                threshold_temperature_k: 361.0,
                ..SOURCE_INPUT
            },
        ] {
            assert!(
                MaterialFirstModeSlabCooling::from_card(
                    &source,
                    &card,
                    &conductivity_key(),
                    SelectionPolicy::SingleClaimOnly,
                    changed_input
                )
                .is_err()
            );
        }
        let curve = source_card(0.2, true, [280.0, 360.0]);
        assert_eq!(
            MaterialFirstModeSlabCooling::from_card(
                &source_capacity(&curve),
                &curve,
                &conductivity_key(),
                SelectionPolicy::SingleClaimOnly,
                SOURCE_INPUT
            ),
            Err(MaterialFirstModeSlabError::NonconstantConductivity)
        );
        let narrow = source_card(0.2, false, [300.0, 340.0]);
        assert!(
            MaterialFirstModeSlabCooling::from_card(
                &source_capacity(&narrow),
                &narrow,
                &conductivity_key(),
                SelectionPolicy::SingleClaimOnly,
                SOURCE_INPUT
            )
            .is_err()
        );
        assert!(
            MaterialFirstModeSlabCooling::from_card(
                &source,
                &card,
                &PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
                SelectionPolicy::SingleClaimOnly,
                SOURCE_INPUT
            )
            .is_err()
        );
    }
}
