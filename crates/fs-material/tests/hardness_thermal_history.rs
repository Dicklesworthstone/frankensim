//! Synthetic original-source/native-history regressions, not material calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::hardness::thermal_history::{
    ThermoplasticHardnessSource, ThermoplasticHardnessStepError,
};
use fs_material::phase::{
    HeatCapacitySourceQuery, PhaseStateError, ResolvedHeatCapacityCurve, SolidLiquidPhase,
    UniformEnthalpyStep, UniformEnthalpyStepInput,
};
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY,
};
use fs_material::plastic::thermal::{CaloricPoint, FeedbackConfig, FeedbackError, SourceQuery};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, QtyAny, QuantitySpec, Temperature};

const CP: f64 = 500.0;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    hardness: [PropertyKey; 2],
}

impl Fixture {
    fn new(maximum_alpha: f64, maximum_hardness_temperature: f64, revision: u32) -> Self {
        let provenance = Provenance {
            source: "synthetic thermoplastic indentation specimen".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        };
        let mut claims = ClaimSet::new();
        for (name, dims, value, slope) in [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9, -2e6),
            (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2, 0.0),
            (YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6, -1e4),
            (
                ISOTROPIC_HARDENING_MODULUS_PROPERTY,
                Pressure::DIMS,
                1e8,
                -1e6,
            ),
            (
                KINEMATIC_HARDENING_MODULUS_PROPERTY,
                Pressure::DIMS,
                1.6e8,
                0.0,
            ),
            ("Cp", CP_DIMS, CP, 0.0),
            (DENSITY_PROPERTY, Density::DIMS, 1000.0, 0.0),
            (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, 6e6, 0.0),
            (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, 300.0, 0.0),
            ("fast_C", Pressure::DIMS, 4e8, -8e6),
            ("slow_C", Pressure::DIMS, 1e8, 2e6),
            ("fast_gamma", Dims::NONE, 120.0, 0.0),
            ("slow_gamma", Dims::NONE, 15.0, 0.0),
            (
                OVERSTRESS_VISCOSITY_PROPERTY,
                OVERSTRESS_VISCOSITY_DIMS,
                2e6,
                0.0,
            ),
            (OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 1e6, 0.0),
            (OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, 3.0, 0.0),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: if slope == 0.0 {
                        PropertyValue::Scalar { value, dims }
                    } else {
                        PropertyValue::Curve {
                            abscissa: "T".into(),
                            abscissa_dims: Temperature::DIMS,
                            knots: [280.0, 300.0, 320.0]
                                .map(|t| (t, value + slope * (t - 300.0)))
                                .to_vec(),
                            dims,
                        }
                    },
                    validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: if slope == 0.0 {
                        InterpolationPolicy::ConstantWithinValidity
                    } else {
                        InterpolationPolicy::LinearInside
                    },
                    observations: vec![],
                    provenance: provenance.clone(),
                })
                .unwrap();
        }
        let observation = claims
            .register_observation(ObservationDataset {
                specimen: "synthetic solid at declared temperature and accumulated strain".into(),
                method: "synthetic force/projected-area indentation".into(),
                artifact: fs_blake3::hash_domain("fixture", b"thermal history hardness"),
                caveats: "numerical regression only".into(),
                provenance: provenance.clone(),
            })
            .unwrap();
        let hardness = [10.0, 20.0].map(|load| {
            PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
                .with_hardness_test(
                    HardnessTestContext::new(
                        "synthetic spherical indenter, radius 1 mm",
                        vec![
                            HardnessLoadStep::new(
                                QtyAny::new(load, Dims([1, 1, -2, 0, 0, 0])),
                                QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
                            )
                            .unwrap(),
                        ],
                        "synthetic projected-area test revision 1",
                        observation,
                    )
                    .unwrap(),
                )
                .unwrap()
        });
        for (index, key) in hardness.iter().enumerate() {
            let scale = 1.0 + index as f64;
            claims
                .insert_claim(PropertyClaim {
                    key: key.clone(),
                    value: PropertyValue::Surface2D {
                        axes: [ALPHA.into(), "T".into()],
                        axis_dims: [Dims::NONE, Temperature::DIMS],
                        samples: vec![
                            ([0.0, 280.0], scale * 9e8),
                            ([maximum_alpha, 280.0], scale * 1.2e9),
                            ([maximum_alpha, maximum_hardness_temperature], scale * 6e8),
                            ([0.0, maximum_hardness_temperature], scale * 4e8),
                        ],
                        triangles: vec![[0, 1, 2], [0, 2, 3]],
                        dims: Pressure::DIMS,
                    },
                    validity: ValidityDomain::unconstrained()
                        .with("T", 280.0, maximum_hardness_temperature)
                        .with_quantity(
                            ALPHA,
                            QuantitySpec::dimensional(Dims::NONE),
                            0.0,
                            maximum_alpha,
                        ),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![observation],
                    provenance: provenance.clone(),
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermoplastic solid".into(),
                phase: "solid".into(),
                process: "synthetic thermal hardness fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let card = match revision {
            0 => card,
            1 => MaterialCard::supersede(&card, card.claims().clone(), vec![]).unwrap(),
            _ => panic!("fixture supports its original card and one successor"),
        };
        let lower = QueryPoint::new().with("T", 280.0).unwrap();
        let upper = QueryPoint::new().with("T", 320.0).unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("Cp", CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * CP,
            1e-6,
        )
        .unwrap();
        Self {
            card,
            caloric,
            hardness,
        }
    }

    fn query(&self) -> SourceQuery<'_> {
        SourceQuery {
            card: &self.card,
            caloric: &self.caloric,
            context: QueryPoint::new().with("T", 300.0).unwrap(),
            temperature_axis: "T".into(),
            selection: SELECT,
            hardening_basis: HardeningPropertyBasis::MixedPlasticStrainModuli,
            thermal_expansion: None,
            strain_limit: 0.04,
        }
    }
    fn point(&self) -> CaloricPoint<'_> {
        CaloricPoint::try_new(self.query(), 0.0, 1e-6, 1e-3).unwrap()
    }
}

fn config() -> FeedbackConfig {
    FeedbackConfig {
        max_iterations: 64,
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        relaxation: 1.0,
    }
}

fn transport(
    input: UniformEnthalpyStepInput<'_>,
    external_heat_j: f64,
) -> Result<UniformEnthalpyStep<f64>, PhaseStateError> {
    Ok(UniformEnthalpyStep {
        state: input.curve.advance_specific_energy(
            input.initial,
            (input.internal_heat_j + external_heat_j) / input.mass_kg,
        )?,
        external_heat_j,
        energy_residual_tolerance_j: 1e-8,
        report: external_heat_j,
    })
}

#[test]
fn g1_actual_caloric_temperature_and_native_accumulation_share_the_original_hardness_test() {
    let fixture = Fixture::new(0.04, 320.0, 0);
    let mut point = fixture.point();
    let binding =
        ThermoplasticHardnessSource::try_new(&point, fixture.hardness[0].clone(), SELECT).unwrap();
    let other =
        ThermoplasticHardnessSource::try_new(&point, fixture.hardness[1].clone(), SELECT).unwrap();
    let selected = binding
        .hardness_source()
        .property()
        .answer()
        .receipt
        .selected;
    let mut heat_j = 0.0;
    for (axial, external_heat) in [(0.012, 6.0), (0.0118, -3.0), (-0.009, -3.0)] {
        let strain = [axial, -0.5 * axial, -0.5 * axial, 0.001, 0.0, 0.0];
        let previous = binding.at_accepted(&point).unwrap();
        let native = point
            .advance(
                strain,
                0.1,
                config(),
                |input| transport(input, external_heat),
                || false,
            )
            .unwrap();
        let step = binding
            .advance(
                &point,
                strain,
                0.1,
                config(),
                |input| transport(input, external_heat),
                || false,
            )
            .unwrap();
        assert_eq!(step.feedback.mechanical, native.mechanical);
        assert_eq!(step.feedback.next.history(), native.next.history());
        assert_eq!(
            step.feedback.next.caloric_state(),
            native.next.caloric_state()
        );
        assert_eq!(step.feedback.next.resolved(), native.next.resolved());
        assert_eq!(step.feedback.transport.report, native.transport.report);
        assert_eq!(
            [
                step.feedback.mechanical_work_j,
                step.feedback.physical_heat_j,
                step.feedback.reversible_heat_j,
                step.feedback.coefficient_storage_change_j,
                step.feedback.thermoelastic_storage_change_j,
                step.feedback.combined_residual_j
            ],
            [
                native.mechanical_work_j,
                native.physical_heat_j,
                native.reversible_heat_j,
                native.coefficient_storage_change_j,
                native.thermoelastic_storage_change_j,
                native.combined_residual_j
            ]
        );
        assert_eq!(
            step.hardness.resolved().state_coordinate("T"),
            Some(native.next.caloric_state().temperature_k())
        );
        assert_eq!(
            step.hardness.resolved().state_coordinate(ALPHA),
            Some(native.next.history().alpha)
        );
        assert_eq!(step.hardness.property().answer().receipt.selected, selected);
        assert_eq!(step.hardness.key(), &fixture.hardness[0]);
        assert_eq!(
            step.hardness_change_pa,
            step.hardness.hardness_pa() - previous.hardness_pa()
        );
        let expected = fixture
            .card
            .claims()
            .query_pinned_typed(
                &fixture.hardness[0],
                &step.hardness.resolved().query(),
                selected,
            )
            .unwrap();
        assert_eq!(step.hardness.property().answer(), &expected);
        let second_test = other.at_accepted(&native.next).unwrap();
        assert_eq!(second_test.hardness_pa(), 2.0 * step.hardness.hardness_pa());
        heat_j += external_heat + step.feedback.physical_heat_j + step.feedback.reversible_heat_j;
        point = step.feedback.next;
    }
    assert!(point.history().alpha > 0.01);
    // The external heat sums to zero. The isotropic hardening modulus falls by
    // 1 % per K, so the stored hardening carries entropy: its reversible heat,
    // about -T |dH/dT| alpha^2 / 2 per volume (-0.074 J at alpha = 0.022),
    // outweighs the plastic heat (about Y alpha, 0.066 J). The specimen ends
    // a few hundredths of a kelvin below 300 K. Its actual caloric temperature
    // is the first-law endpoint of the reported heats on the constant-Cp chart
    // (h = Cp (T - 300 K)).
    let expected = 300.0 + heat_j / (point.mass_kg() * CP);
    let temperature = point.caloric_state().temperature_k();
    assert!(
        (temperature - expected).abs() <= 1e-9,
        "caloric temperature {temperature} K, first-law endpoint {expected} K"
    );
}

#[test]
fn g3_viscous_thermal_chaboche_keeps_all_native_backstresses_and_recovery_heat() {
    let fixture = Fixture::new(0.04, 320.0, 0);
    let mut point = fixture
        .point()
        .with_voce_hardening(SELECT)
        .unwrap()
        .with_power_law_overstress(SELECT)
        .unwrap()
        .with_chaboche_branches(
            SELECT,
            [
                BackstressProperties::new("fast_C", "fast_gamma"),
                BackstressProperties::new("slow_C", "slow_gamma"),
            ],
        )
        .unwrap();
    let binding =
        ThermoplasticHardnessSource::try_new(&point, fixture.hardness[0].clone(), SELECT).unwrap();
    for (axial, heat, duration) in [(0.012, 0.0, 2e-4), (0.012, 5.0, 4e-3), (-0.009, -2.0, 2e-3)] {
        let strain = [axial, -0.5 * axial, -0.5 * axial, 0.001, 0.0, 0.0];
        let native = point
            .advance(
                strain,
                duration,
                config(),
                |input| transport(input, heat),
                || false,
            )
            .unwrap();
        let step = binding
            .advance(
                &point,
                strain,
                duration,
                config(),
                |input| transport(input, heat),
                || false,
            )
            .unwrap();
        assert_eq!(step.feedback.next.history(), native.next.history());
        assert_eq!(step.feedback.mechanical, native.mechanical);
        assert_eq!(step.feedback.physical_heat_j, native.physical_heat_j);
        assert_eq!(
            step.feedback.coefficient_storage_change_j,
            native.coefficient_storage_change_j
        );
        assert_eq!(
            step.feedback.next.caloric_state(),
            native.next.caloric_state()
        );
        assert_ne!(
            native.next.history().backstress_pa[0],
            native.next.history().backstress_pa[1]
        );
        assert_eq!(
            step.hardness.resolved().state_coordinate(ALPHA),
            Some(native.next.history().plastic.alpha)
        );
        assert_eq!(
            step.hardness.resolved().state_coordinate("T"),
            Some(native.next.caloric_state().temperature_k())
        );
        point = step.feedback.next;
    }
    assert!(point.history().plastic.alpha > 0.0);
}

#[test]
fn g0_missing_history_or_temperature_coverage_rejects_the_whole_candidate() {
    for (maximum_alpha, maximum_temperature, heat) in [(1e-5, 320.0, 0.0), (0.04, 305.0, 5.0)] {
        let fixture = Fixture::new(maximum_alpha, maximum_temperature, 0);
        let point = fixture.point();
        let binding =
            ThermoplasticHardnessSource::try_new(&point, fixture.hardness[0].clone(), SELECT)
                .unwrap();
        let old_history = point.history().clone();
        let old_thermal = *point.caloric_state();
        let strain = [0.012, -0.006, -0.006, 0.0, 0.0, 0.0];
        assert!(
            point
                .advance(
                    strain,
                    0.1,
                    config(),
                    |input| transport(input, heat),
                    || false
                )
                .is_ok()
        );
        assert!(matches!(
            binding.advance(
                &point,
                strain,
                0.1,
                config(),
                |input| transport(input, heat),
                || false
            ),
            Err(ThermoplasticHardnessStepError::Hardness(_))
        ));
        assert_eq!(point.history(), &old_history);
        assert_eq!(point.caloric_state(), &old_thermal);
        assert_eq!(point.time_s(), 0.0);
        assert_eq!(
            binding.at_accepted(&point).unwrap(),
            *binding.hardness_source()
        );
    }
}

#[test]
fn g0_original_sources_and_native_cancellation_cannot_be_replaced_by_equal_state_coordinates() {
    let fixture = Fixture::new(0.04, 320.0, 0);
    let point = fixture.point();
    let binding =
        ThermoplasticHardnessSource::try_new(&point, fixture.hardness[0].clone(), SELECT).unwrap();
    let replacement = Fixture::new(0.04, 320.0, 1);
    assert!(binding.at_accepted(&replacement.point()).is_err());
    let mut changed_query = fixture.query();
    changed_query.hardening_basis = HardeningPropertyBasis::PlasticStrainModulus;
    let changed = CaloricPoint::try_new(changed_query, 0.0, 1e-6, 1e-3).unwrap();
    assert_eq!(point.history(), changed.history());
    assert_eq!(point.caloric_state(), changed.caloric_state());
    assert!(binding.at_accepted(&changed).is_err());
    assert!(
        ThermoplasticHardnessSource::try_new(&point, PropertyKey::new("H", Pressure::DIMS), SELECT)
            .is_err()
    );
    let mut calls = 0;
    let result = binding.advance(
        &point,
        [0.0; 6],
        0.1,
        config(),
        |input| {
            calls += 1;
            transport(input, 0.0)
        },
        || true,
    );
    assert!(matches!(
        result,
        Err(ThermoplasticHardnessStepError::Feedback(
            FeedbackError::Cancelled
        ))
    ));
    assert_eq!(calls, 0);
    assert_eq!(point.time_s(), 0.0);
}
