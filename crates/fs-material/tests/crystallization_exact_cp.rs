//! G1/G3: exact original Cp, phase density, live properties and elastic release.
//! Manufactured sources establish numerical behavior, not measured validity.

use std::sync::LazyLock;

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::phase::crystallization::mechanics::{
    CrystallizationLock, CrystallizationStrainLimits,
};
use fs_material::phase::crystallization::properties::{
    CrystallizationMaterialPoint, CrystallizationPropertySource,
};
use fs_material::phase::crystallization::source::CrystallizationSourcePair;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermal::transport::expansion::{
    FreeExpansionSourceQuery, FreeExpansionStrainSources, PreparedFreeExpansionDensity,
};
use fs_material::thermal::transport::isobaric::crystallization::{
    IsobaricCrystallizationPair, IsobaricCrystallizationQuery,
};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_thermal::crystallization::pressure::IsobaricCrystallizationControls;
use fs_thermal::crystallization::pressure::parcel::IsobaricCrystallizationParcel;
use fs_thermal::crystallization::{
    CrystallizationError, CrystallizationPopulation, CrystallizationSolve,
};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const ENERGY: Dims = Dims([2, 0, -2, 0, 0, 0]);
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

// Pair preparation borrows one query lifetime alongside the original Cp laws.
static REFERENCE: LazyLock<QueryPoint> = LazyLock::new(|| point(400.0, 100_000.0));
static LATENT_KEY: LazyLock<PropertyKey> = LazyLock::new(|| PropertyKey::new("latent", ENERGY));

fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), p)
        .unwrap()
}

fn key(prefix: &str, phase: usize, dims: Dims) -> PropertyKey {
    PropertyKey::new(format!("{prefix}_{phase}"), dims)
}

fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut insert = |key: PropertyKey, value: PropertyValue, validity: ValidityDomain| {
        let interpolation = match &value {
            PropertyValue::Scalar { .. } => InterpolationPolicy::ConstantWithinValidity,
            _ => InterpolationPolicy::LinearInside,
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured exact-Cp phase/property/shrinkage reference".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let support = |low, high| {
        ValidityDomain::unconstrained()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), low, high)
            .with_quantity(
                "p",
                QuantitySpec::dimensional(Pressure::DIMS),
                100_000.0,
                100_000.0,
            )
    };
    for phase in 0..2 {
        insert(
            key("cp", phase, CP),
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, [900.0, 600.0][phase]),
                    (500.0, [1500.0, 1000.0][phase]),
                ],
                dims: CP,
            },
            support(300.0, 500.0),
        );
        insert(
            key("rho", phase, DENSITY),
            PropertyValue::Scalar {
                value: [1000.0, 1050.0][phase],
                dims: DENSITY,
            },
            support(400.0, 400.0),
        );
        insert(
            key("length", phase, Dims::NONE),
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, [-0.01, -0.005][phase]),
                    (500.0, [0.01, 0.005][phase]),
                ],
                dims: Dims::NONE,
            },
            support(300.0, 500.0),
        );
    }
    insert(
        PropertyKey::new("latent", ENERGY),
        PropertyValue::Scalar {
            value: 80_000.0,
            dims: ENERGY,
        },
        support(400.0, 400.0),
    );
    for (name, dims, axis, knots) in [
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            "X",
            vec![(0.0, 1.0), (1.0, 5.0)],
        ),
        (
            "volume_observation",
            Dims::NONE,
            "V",
            vec![(0.0, 0.0), (1.0, 1.0)],
        ),
    ] {
        insert(
            PropertyKey::new(name, dims),
            PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: Dims::NONE,
                knots,
                dims,
            },
            support(300.0, 500.0).with_quantity(
                axis,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                1.0,
            ),
        );
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic one-composition fixture".into(),
            phase: "parent and crystalline product".into(),
            process: "fixed pressure".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn caloric(card: &MaterialCard, phase: usize) -> IsobaricHeatCapacityLaw {
    IsobaricHeatCapacityLaw::try_new(
        card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &key("cp", phase, CP),
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &point(300.0, 100_000.0),
            upper: &point(500.0, 100_000.0),
            selection: SELECTION,
        },
        400.0,
    )
    .unwrap()
}

fn density(card: &MaterialCard, phase: usize) -> PreparedFreeExpansionDensity {
    PreparedFreeExpansionDensity::try_new(
        card,
        FreeExpansionSourceQuery {
            density: &key("rho", phase, DENSITY),
            density_selection: SELECTION,
            reference: &point(400.0, 100_000.0),
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &key("length", phase, Dims::NONE),
                selection: SELECTION,
            },
            temperature_axis: "T",
            lower: &point(300.0, 100_000.0),
            upper: &point(500.0, 100_000.0),
        },
    )
    .unwrap()
}

fn pair<'a>(
    card: &MaterialCard,
    cp: &'a [IsobaricHeatCapacityLaw; 2],
    rho: &'a [PreparedFreeExpansionDensity; 2],
) -> IsobaricCrystallizationPair<'a> {
    IsobaricCrystallizationPair::try_new_with_free_expansion(
        card,
        IsobaricCrystallizationQuery {
            parent: &cp[0],
            crystal: &cp[1],
            enthalpy_gap: &LATENT_KEY,
            reference: &REFERENCE,
            selection: SELECTION,
            temperature_range_k: [300.0, 500.0],
        },
        &rho[0],
        &rho[1],
    )
    .unwrap()
}

fn properties<'a>(
    card: &'a MaterialCard,
    query: &'a QueryPoint,
) -> CrystallizationMaterialPoint<'a> {
    CrystallizationMaterialPoint {
        material: card,
        context: query,
        temperature_axis: "T",
        crystalline_mass_fraction_axis: Some("X"),
        crystalline_volume_fraction_axis: Some("V"),
    }
}

fn requirements() -> [ScalarPropertyRequirement; 2] {
    [
        ScalarPropertyRequirement::try_new(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        ScalarPropertyRequirement::try_new(
            "volume_observation",
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
    ]
}

fn context(p: f64) -> QueryPoint {
    point(350.0, p)
        .with_quantity("X", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
        .with_quantity("V", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
}

fn rate(_: f64) -> Result<f64, CrystallizationError> {
    Ok(0.3)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn specific_volume(t: f64, x: f64) -> f64 {
    (1.0 - x) * (1.0 + 1e-4 * (t - 400.0)).powi(3) / 1000.0
        + x * (1.0 + 5e-5 * (t - 400.0)).powi(3) / 1050.0
}
fn limits() -> CrystallizationStrainLimits {
    CrystallizationStrainLimits {
        maximum_absolute_total_mandel_strain: 0.03,
        maximum_absolute_elastic_mandel_strain: 0.03,
        maximum_absolute_free_isotropic_strain: 0.02,
    }
}
fn controls() -> IsobaricCrystallizationControls {
    IsobaricCrystallizationControls {
        caloric: CrystallizationSolve::default(),
        energy_tolerance_j: 1e-6,
        volume_tolerance_m3: 1e-15,
    }
}

#[test]
fn g1_exact_cp_recalescence_drives_original_mass_and_volume_property_coordinates() {
    let card = card();
    let cp = [caloric(&card, 0), caloric(&card, 1)];
    let rho = [density(&card, 0), density(&card, 1)];
    let pair = pair(&card, &cp, &rho);
    let law = pair
        .law(&rate, ContentHash([31; 32]), 1.0, 0.8, 550.0)
        .unwrap();
    let initial = law.state_at_temperature(350.0, 0.0).unwrap();
    let query = context(100_000.0);
    let source = CrystallizationPropertySource::new(
        &pair,
        law,
        initial,
        properties(&card, &query),
        &requirements(),
        MaterialPropertySelection::SingleClaimOnly,
        CrystallizationSolve::default(),
    )
    .unwrap();
    let grown = law
        .advance_specific_enthalpy(
            initial,
            initial.specific_enthalpy_j_kg(),
            1.0,
            CrystallizationSolve::default(),
            || false,
        )
        .unwrap()
        .state;
    let response = source.evaluate(grown).unwrap();
    let x = 0.8 * (1.0 - (-0.3_f64).exp());
    // h(T,X)=(1.5-X/2)*T²-240000, exactly from Cp_parent=3T and Cp_crystal=2T.
    let t = (1.5 * 350.0_f64.powi(2) / (1.5 - 0.5 * x)).sqrt();
    close(grown.temperature_k(), t, 1e-8);
    close(grown.crystalline_mass_fraction(), x, 1e-14);
    close(
        response.volume().specific_volume_m3_kg,
        specific_volume(t, x),
        1e-14,
    );
    let phi = (x * (1.0 + 5e-5 * (t - 400.0)).powi(3) / 1050.0) / specific_volume(t, x);
    close(response.volume().crystalline_volume_fraction, phi, 1e-11);
    close(
        response
            .material()
            .property(THERMAL_CONDUCTIVITY_PROPERTY)
            .unwrap()
            .value_si(),
        1.0 + 4.0 * x,
        1e-12,
    );
    close(
        response
            .material()
            .property("volume_observation")
            .unwrap()
            .value_si(),
        phi,
        1e-11,
    );
    assert_eq!(response.volume().caloric_pair_identity, pair.identity());
    assert_eq!(response.material().state_coordinate("p"), Some(100_000.0));
    for property in response.material().properties() {
        assert_eq!(
            property.answer().receipt.selected,
            source
                .original_source()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    let population = CrystallizationPopulation::new(
        law,
        &[(1.0, initial), (2.0, grown)],
        2,
        CrystallizationSolve::default(),
        || false,
    )
    .unwrap();
    let population_source = CrystallizationPropertySource::for_population(
        &pair,
        &population,
        properties(&card, &query),
        &requirements(),
        MaterialPropertySelection::SingleClaimOnly,
        CrystallizationSolve::default(),
        || false,
    )
    .unwrap();
    let answer = population_source
        .evaluate_population(&population, || false)
        .unwrap();
    close(
        answer.material().state_coordinate("X").unwrap(),
        population.crystalline_mass_fraction(),
        0.0,
    );
    close(
        answer.volume().specific_volume_m3_kg,
        specific_volume(
            population.temperature_k(),
            population.crystalline_mass_fraction(),
        ),
        1e-14,
    );
}

#[test]
fn g1_exact_density_drives_free_and_constrained_shrinkage_without_replacing_calorics() {
    let card = card();
    let cp = [caloric(&card, 0), caloric(&card, 1)];
    let rho = [density(&card, 0), density(&card, 1)];
    let pair = pair(&card, &cp, &rho);
    let lock = CrystallizationLock::new(&pair, 350.0, 0.0, [0.0; 6]).unwrap();
    let elastic =
        ElasticAcousticState::from_bulk_and_shear(1000.0, 1e10, 4e9, ContentHash([32; 32]), 0.03)
            .unwrap();
    let frozen = lock.freeze(&pair, 350.0, 0.4, &elastic, limits()).unwrap();
    let shrinkage = (specific_volume(350.0, 0.4) / specific_volume(350.0, 0.0)).cbrt() - 1.0;
    close(frozen.free_isotropic_strain(), shrinkage, 2e-15);
    assert!(shrinkage < 0.0);
    for stress in frozen.stress(frozen.free_strain(), &()) {
        close(stress, 0.0, 0.0);
    }
    for stress in &frozen.stress(&[0.0; 6], &())[..3] {
        close(*stress, -3e10 * shrinkage, 1e-4);
    }
    assert!(lock.freeze(&pair, 350.0, -0.1, &elastic, limits()).is_err());
    close(
        pair.latent_enthalpy_j_kg(350.0).unwrap(),
        0.5 * 350.0_f64.powi(2),
        1e-9,
    );
}

#[test]
fn g1_exact_cp_unlocking_closes_elastic_heat_and_original_pressure_work() {
    let card = card();
    let cp = [caloric(&card, 0), caloric(&card, 1)];
    let rho = [density(&card, 0), density(&card, 1)];
    let pair = pair(&card, &cp, &rho);
    let law = pair
        .law(&rate, ContentHash([33; 32]), 1.0, 0.8, 550.0)
        .unwrap();
    let initial = law.state_at_temperature(350.0, 0.6).unwrap();
    let parcel = pair
        .isobaric_parcel(&law, initial, 2.0, 0.002, controls(), || false)
        .unwrap();
    let elastic =
        ElasticAcousticState::from_bulk_and_shear(1000.0, 1e10, 4e9, ContentHash([34; 32]), 0.03)
            .unwrap();
    let lock = CrystallizationLock::new(
        &pair,
        initial.temperature_k(),
        initial.crystalline_mass_fraction(),
        [0.0; 6],
    )
    .unwrap();
    let release = || {
        lock.freeze(
            &pair,
            initial.temperature_k(),
            initial.crystalline_mass_fraction(),
            &elastic,
            limits(),
        )
        .unwrap()
        .unlock_at_fixed_strain([0.01, -0.01, 0.0, 0.0, 0.0, 0.0], 0.001, 0.75)
        .unwrap()
    };
    let (warm, report) = release()
        .thermalize_isobaric(&pair, &parcel, controls(), || false)
        .unwrap();
    close(report.removed_elastic_energy_j, 1600.0, 1e-9);
    close(report.thermalized_energy_j, 1200.0, 1e-9);
    close(report.exported_energy_j, 400.0, 1e-9);
    close(
        report.thermal.internal_energy_change_j
            + report.thermal.boundary_work_j
            + report.exported_energy_j,
        report.removed_elastic_energy_j,
        1e-6,
    );
    close(
        warm.volume_m3(),
        2.0 * specific_volume(
            warm.state().temperature_k(),
            initial.crystalline_mass_fraction(),
        ),
        1e-14,
    );
    assert!(report.thermal.boundary_work_j > 0.0);
    assert_eq!(warm.time_s(), parcel.time_s());
    assert_eq!(warm.state().exposure(), initial.exposure());
    let relock = report
        .state
        .relock_at(
            &pair,
            warm.state().temperature_k(),
            warm.state().crystalline_mass_fraction(),
            [0.0; 6],
        )
        .unwrap();
    assert_eq!(relock.caloric_pair_identity(), pair.identity());
    let wrong_pressure = IsobaricCrystallizationParcel::new(
        &law,
        initial,
        2.0,
        0.002,
        200_000.0,
        controls(),
        || false,
    )
    .unwrap();
    assert!(
        release()
            .thermalize_isobaric(&pair, &wrong_pressure, controls(), || false)
            .is_err()
    );
    assert!(
        release()
            .thermalize_isobaric(&pair, &parcel, controls(), || true)
            .is_err()
    );
    assert_eq!(parcel.state(), initial);
}

#[test]
fn g0_missing_density_and_mismatched_exact_source_context_refuse() {
    let card = card();
    let cp = [caloric(&card, 0), caloric(&card, 1)];
    let rho = [density(&card, 0), density(&card, 1)];
    let pair = pair(&card, &cp, &rho);
    let law = pair
        .law(&rate, ContentHash([35; 32]), 1.0, 0.8, 550.0)
        .unwrap();
    let initial = law.state_at_temperature(350.0, 0.0).unwrap();
    let missing_pressure = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 350.0)
        .unwrap()
        .with_quantity("X", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
        .with_quantity("V", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    for query in [context(200_000.0), missing_pressure] {
        assert!(
            CrystallizationPropertySource::new(
                &pair,
                law,
                initial,
                properties(&card, &query),
                &requirements(),
                MaterialPropertySelection::SingleClaimOnly,
                CrystallizationSolve::default()
            )
            .is_err()
        );
    }
    let plain = IsobaricCrystallizationPair::try_new(
        &card,
        IsobaricCrystallizationQuery {
            parent: &cp[0],
            crystal: &cp[1],
            enthalpy_gap: &LATENT_KEY,
            reference: &REFERENCE,
            selection: SELECTION,
            temperature_range_k: [300.0, 500.0],
        },
    )
    .unwrap();
    assert!(
        CrystallizationSourcePair::from(&plain)
            .volume_state(350.0, 0.0)
            .is_err()
    );
    assert!(CrystallizationLock::new(&plain, 350.0, 0.0, [0.0; 6]).is_err());
    let query = context(100_000.0);
    assert!(
        CrystallizationPropertySource::new(
            &plain,
            law,
            initial,
            properties(&card, &query),
            &requirements(),
            MaterialPropertySelection::SingleClaimOnly,
            CrystallizationSolve::default()
        )
        .is_err()
    );
    assert!(
        CrystallizationSourcePair::from(&pair)
            .volume_state(299.0, 0.0)
            .is_err()
    );
}

#[test]
fn g1_exact_cp_contact_solves_live_bulk_and_interface_crystallinity_together() {
    use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
    use fs_material::thermal::interface::crystallization::{
        CrystallizationBulkConduction, CrystallizationInterfaceContact,
        CrystallizationInterfacePhaseAxes,
    };
    use fs_material::thermal::interface::{
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        ThermalInterfacePoint,
    };
    use fs_thermal::crystallization::pressure::exchange::{
        IsobaricExchangeControls, IsobaricExchangeSideControls,
    };

    let card = card();
    let cp = [caloric(&card, 0), caloric(&card, 1)];
    let rho = [density(&card, 0), density(&card, 1)];
    let pair = pair(&card, &cp, &rho);
    let law = pair
        .law(&rate, ContentHash([36; 32]), 1.0, 0.8, 550.0)
        .unwrap();
    let bodies = [400.0, 350.0].map(|temperature| {
        pair.isobaric_parcel(
            &law,
            law.state_at_temperature(temperature, 0.0).unwrap(),
            1.0,
            0.001,
            controls(),
            || false,
        )
        .unwrap()
    });
    let bulk_context = context(100_000.0);
    let bulk = CrystallizationBulkConduction::new(
        [properties(&card, &bulk_context); 2],
        [&pair; 2],
        bodies.each_ref(),
        [0.001; 2],
        [
            MaterialPropertySelection::SingleClaimOnly,
            MaterialPropertySelection::SingleClaimOnly,
        ],
        [CrystallizationSolve::default(); 2],
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "X_b".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 1e-4), (0.8, 2.6e-4)],
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ["T_a", "T_b"]
                .into_iter()
                .fold(ValidityDomain::unconstrained(), |domain, axis| {
                    domain.with_quantity(
                        axis,
                        QuantitySpec::dimensional(Temperature::DIMS),
                        300.0,
                        500.0,
                    )
                })
                .with_quantity("X_b", QuantitySpec::dimensional(Dims::NONE), 0.0, 0.8)
                .with_quantity("V_b", QuantitySpec::dimensional(Dims::NONE), 0.0, 0.8),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured live crystalline contact resistance".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |label: &str| SurfaceSpec {
        material: card.id().clone(),
        texture_frame: label.into(),
    };
    let interface = InterfaceSystemCard::assemble(
        surface("left"),
        surface("right"),
        SystemContext {
            medium: "fixed contact layer".into(),
            third_body: None,
            environment: "fixed area and layer distances".into(),
            history: "uniform bulk crystalline fractions at the interface".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let interface_context = QueryPoint::new()
        .with_quantity("T_a", QuantitySpec::dimensional(Temperature::DIMS), 400.0)
        .unwrap()
        .with_quantity("T_b", QuantitySpec::dimensional(Temperature::DIMS), 350.0)
        .unwrap()
        .with_quantity("X_b", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap()
        .with_quantity("V_b", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    let contact = CrystallizationInterfaceContact::new_with_crystallinity(
        ThermalInterfacePoint {
            interface: &interface,
            materials: [&card; 2],
            context: &interface_context,
            temperature_axes: ["T_a", "T_b"],
            liquid_fraction_axes: [None; 2],
            mean_temperature_axis: None,
        },
        [&pair; 2],
        bodies.each_ref(),
        0.01,
        None,
        None,
        CrystallizationInterfacePhaseAxes {
            crystalline_mass_fraction: [None, Some("X_b")],
            crystalline_volume_fraction: [None, Some("V_b")],
        },
        MaterialPropertySelection::SingleClaimOnly,
        Some(bulk),
    )
    .unwrap();
    let side = IsobaricExchangeSideControls {
        caloric: controls(),
        ..Default::default()
    };
    let duration = 0.5;
    let (left, right, report) = contact
        .exchange_heat(
            &bodies[0],
            &bodies[1],
            duration,
            side,
            side,
            IsobaricExchangeControls {
                heat_bracket_j: [0.0, 2000.0],
                maximum_evaluations: 64,
                residual_tolerance_j: 1e-6,
                energy_tolerance_j: 1e-5,
            },
            || false,
        )
        .unwrap();

    // Independent analytic calorics and kinetics close the same two-body
    // backward-Euler heat equation; both conductivity and resistance change.
    let x = 0.8 * (1.0 - (-0.3_f64 * duration).exp());
    let conductivity = 1.0 + 4.0 * x;
    let resistance = 1e-4 * (1.0 + 2.0 * x);
    let conductance = 0.01 / (resistance + 0.002 / conductivity);
    let temperatures = |q: f64| {
        [
            ((1.5 * 400.0_f64.powi(2) - q) / (1.5 - 0.5 * x)).sqrt(),
            ((1.5 * 350.0_f64.powi(2) + q) / (1.5 - 0.5 * x)).sqrt(),
        ]
    };
    let mut bracket = [0.0, 2000.0];
    for _ in 0..64 {
        let heat = 0.5 * (bracket[0] + bracket[1]);
        let t = temperatures(heat);
        if heat > conductance * duration * (t[0] - t[1]) {
            bracket[1] = heat;
        } else {
            bracket[0] = heat;
        }
    }
    let expected_heat = 0.5 * (bracket[0] + bracket[1]);
    let expected_t = temperatures(expected_heat);
    close(report.thermal.heat_left_to_right_j, expected_heat, 2e-6);
    for (body, temperature) in [&left, &right].into_iter().zip(expected_t) {
        close(body.state().temperature_k(), temperature, 1e-8);
        close(body.state().crystalline_mass_fraction(), x, 1e-14);
    }
    close(report.interface.resistance_m2_k_w, resistance, 1e-16);
    close(report.interface.conductance_w_k, conductance, 1e-12);
    close(
        report.interface.source.state_coordinate("X_b").unwrap(),
        x,
        1e-14,
    );
    let volume_fraction =
        x * (1.0 + 5e-5 * (expected_t[1] - 400.0)).powi(3) / 1050.0
            / specific_volume(expected_t[1], x);
    close(
        report.interface.source.state_coordinate("V_b").unwrap(),
        volume_fraction,
        1e-11,
    );
    for k in report.interface.bulk.as_ref().unwrap().conductivity_w_m_k {
        close(k, conductivity, 1e-12);
    }
    assert_eq!(
        report.interface.source.properties()[0].answer().receipt.selected,
        contact.original_source().properties()[0].answer().receipt.selected
    );
    close(report.thermal.combined_enthalpy_change_j, 0.0, 1e-5);
    close(
        report.thermal.combined_internal_energy_change_j
            + report.thermal.combined_boundary_work_j,
        0.0,
        1e-5,
    );
}
