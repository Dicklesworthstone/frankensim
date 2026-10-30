//! Source curves -> film properties -> actual shared finite-cell heat storage.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    UncertaintyModel,
};
use fs_material::fluid_transport::{ThermalFluidTransport, ThermalFluidTransportState};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::state_point::*;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

// Synthetic source curves, deliberately not presented as measured fluid data.
fn liquid_card(density_slope: f64, viscosity_slope: f64, upper_k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (DENSITY_PROPERTY, Density::DIMS, 1000.0, -density_slope),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            4000.0,
            2.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.5,
            0.002,
        ),
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            0.001,
            -0.001 * viscosity_slope,
        ),
    ] {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [280.0, 320.0]
                        .map(|t| (t, value + slope * (t - 300.0)))
                        .to_vec(),
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, upper_k)
                    .with("p", 100000.0, 100000.0)
                    .with("solute_mass_fraction", 0.2, 0.2),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored film-convection fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "source-only substitutions".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn liquid_law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    let point = QueryPoint::new()
        .with_quantity("T", temperature_quantity(), 300.0)
        .unwrap()
        .with("p", 100000.0)
        .unwrap()
        .with("solute_mass_fraction", 0.2)
        .unwrap();
    LiquidTemperatureLaw::try_new(
        card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn film(cell: CellId) -> SourceConvectiveCell {
    SourceConvectiveCell {
        cell,
        reference_temperature: Temperature::new(300.0),
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: SourceCellConvectionLaw::VerticalPlate { gravity_m_s2: 9.81 },
    }
}
// Independent Churchill-Chu algebra, with film-evaluated material curves but
// the FULL wall-to-reservoir temperature difference driving buoyancy and heat.
fn expected_h(density_slope: f64, viscosity_slope: f64, wall: f64, property: f64) -> f64 {
    let delta = property - 300.0;
    let rho = 1000.0 - density_slope * delta;
    let cp = 4000.0 + 2.0 * delta;
    let k = 0.5 + 0.002 * delta;
    let mu = 0.001 * (1.0 - viscosity_slope * delta);
    let pr = mu * cp / k;
    let beta = density_slope / rho;
    let ra =
        (9.81 * beta * (wall - 300.0) * 0.1_f64.powi(3) / ((mu / rho) * (k / (rho * cp)))).abs();
    (0.825 + 0.387 * ra.powf(1.0 / 6.0) / (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0))
        .powi(2)
        * k
        / 0.1
}
fn endpoint(density_slope: f64, viscosity_slope: f64, old: f64, dt: f64) -> f64 {
    let (mut lower, mut upper) = (old.min(300.0), old.max(300.0));
    for _ in 0..80 {
        let t = 0.5 * lower + 0.5 * upper;
        let power =
            0.008 * expected_h(density_slope, viscosity_slope, t, 0.5 * t + 150.0) * (300.0 - t);
        if MASS * CP * (t - old) - dt * power > 0.0 {
            upper = t;
        } else {
            lower = t;
        }
    }
    0.5 * lower + 0.5 * upper
}

#[test]
fn g1_film_transport_and_signed_expansion_keep_original_material_and_bath() {
    let chart = curve();
    let heat = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    for density_slope in [0.5, -0.5] {
        let card = liquid_card(density_slope, 0.01, 320.0);
        let law = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let mut queries = Vec::new();
        let result = evaluate_source_cell_convection(
            &heat,
            &IDS.map(film),
            2,
            |binding, t| {
                queries.push((binding.cell, t.value()));
                law.at_temperature(t)
            },
            || false,
        )
        .unwrap();
        assert_eq!(queries, vec![(CellId(10), 295.0), (CellId(30), 305.0)]);
        for (source, sample) in result.materials().iter().zip(result.convection().samples()) {
            assert_eq!(source.cell, sample.cell);
            assert_eq!(
                source.material.material().card_identity(),
                card.content_hash()
            );
            assert_eq!(source.material.material().query().axes()["p"], 100000.0);
            assert_eq!(
                source.material.material().query().axes()["solute_mass_fraction"],
                0.2
            );
            assert_eq!(
                source.material.material().query().axes()["T"],
                source.property_temperature_k
            );
            assert_eq!(source.thermal.transport(), sample.transport);
            assert_eq!(
                source
                    .thermal
                    .volumetric_expansion_per_k()
                    .unwrap()
                    .signum(),
                density_slope.signum()
            );
            for (actual, reference) in source
                .material
                .material()
                .properties()
                .iter()
                .zip(law.reference_state().material().properties())
            {
                assert_eq!(
                    actual.answer().receipt.selected,
                    reference.answer().receipt.selected
                );
            }
            let wall = heat.temperature_k(sample.cell).unwrap();
            let h = expected_h(density_slope, 0.01, wall, source.property_temperature_k);
            near(
                sample
                    .evaluation
                    .as_ref()
                    .unwrap()
                    .coefficient
                    .value
                    .value(),
                h,
                1e-8,
            );
            near(sample.heat_into_cell_w, 0.008 * h * (300.0 - wall), 1e-8);
            assert_eq!(sample.reference_temperature_k, 300.0);
        }
        let reservoir = IDS.map(|id| SourceConvectiveCell {
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            ..film(id)
        });
        let frozen = evaluate_source_cell_convection(
            &heat,
            &reservoir,
            2,
            |_, t| law.at_temperature(t),
            || false,
        )
        .unwrap();
        assert!(
            (frozen.convection().sources_w()[&IDS[0]] - result.convection().sources_w()[&IDS[0]])
                .abs()
                > 1e-4
        );
    }
}

#[test]
fn g3_source_viscosity_curve_changes_the_actual_accepted_thermal_endpoint() {
    let chart = curve();
    let heat = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let saved = heat.cells().clone();
    let dt = 0.2;
    let mut endpoints = Vec::new();
    for viscosity_slope in [0.01, 0.03] {
        let card = liquid_card(0.5, viscosity_slope, 320.0);
        let law = liquid_law(&card)
            .with_density_expansion(Temperature::new(280.0), Temperature::new(320.0))
            .unwrap();
        let step = BackwardEulerStep::new(&heat, dt, || Ok(())).unwrap();
        let increments = step
            .cell_ids()
            .iter()
            .map(|&id| {
                let old = heat.temperature_k(id).unwrap();
                MASS * CP * (endpoint(0.5, viscosity_slope, old, dt) - old)
            })
            .collect::<Vec<_>>();
        let trial = step
            .evaluate(
                &increments,
                |state| {
                    let source = evaluate_source_cell_convection(
                        state,
                        &IDS.map(film),
                        2,
                        |_, t| law.at_temperature(t),
                        || false,
                    )?;
                    Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                        sources_w: source.convection().sources_w().clone(),
                        conductances_w_k: BTreeMap::new(),
                        payload: source,
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, source, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        for sample in source.convection().samples() {
            let old = heat.temperature_k(sample.cell).unwrap();
            near(
                sample.wall_temperature_k,
                endpoint(0.5, viscosity_slope, old, dt),
                1e-10,
            );
            near(
                MASS * CP * (sample.wall_temperature_k - old),
                dt * sample.heat_into_cell_w,
                1e-8,
            );
        }
        near(
            next.energy_j() - heat.energy_j(),
            audit.external_heat_j,
            1e-8,
        );
        near(
            audit.external_heat_j,
            dt * source.convection().total_heat_into_cells_w(),
            1e-10,
        );
        endpoints.push(next.temperature_k(IDS[0]).unwrap());
    }
    assert!(endpoints[1] < endpoints[0] && endpoints[1] > 300.0);
    assert_eq!(heat.cells(), &saved);
}

#[test]
fn g0_missing_expansion_support_mismatch_duplicate_and_cancel_leave_history_intact() {
    let chart = curve();
    let heat = initial(&chart, [310.0, 290.0], 0.0, 0.002);
    let saved = heat.cells().clone();
    let card = liquid_card(0.5, 0.01, 320.0);
    let law = liquid_law(&card);
    let b = film(IDS[0]);
    assert!(matches!(
        evaluate_source_cell_convection(&heat, &[b], 1, |_, t| law.at_temperature(t), || false),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Invalid(
                "natural convection requires sourced isobaric volumetric expansion"
            )
        ))
    ));
    let forced = SourceConvectiveCell {
        law: SourceCellConvectionLaw::Forced {
            correlation: CorrelationId::FlatPlateLaminarAverage,
            speed_m_s: 1.0,
            settings: ForcedConvectionSettings::default(),
        },
        ..b
    };
    evaluate_source_cell_convection(&heat, &[forced], 1, |_, t| law.at_temperature(t), || false)
        .unwrap();
    assert!(matches!(
        evaluate_source_cell_convection(
            &heat,
            &[forced],
            1,
            |_, _| law.at_temperature(Temperature::new(300.0)),
            || false
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Invalid(_)
        ))
    ));
    let mut queries = 0;
    assert!(
        evaluate_source_cell_convection(
            &heat,
            &[forced, forced],
            2,
            |_, t| {
                queries += 1;
                law.at_temperature(t)
            },
            || false
        )
        .is_err()
    );
    assert_eq!(queries, 0);
    let limited_card = liquid_card(0.5, 0.01, 302.0);
    let limited = liquid_law(&limited_card);
    assert!(matches!(
        evaluate_source_cell_convection(
            &heat,
            &[forced],
            1,
            |_, t| limited.at_temperature(t),
            || false
        ),
        Err(SourceCellConvectionError::Source {
            cell: CellId(30),
            property_temperature_k: 305.0,
            ..
        })
    ));
    let mut polls = 0;
    evaluate_source_cell_convection(
        &heat,
        &[forced],
        1,
        |_, t| law.at_temperature(t),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut at = 0;
    assert!(matches!(
        evaluate_source_cell_convection(
            &heat,
            &[forced],
            1,
            |_, t| law.at_temperature(t),
            || {
                at += 1;
                at == polls
            }
        ),
        Err(SourceCellConvectionError::Convection(
            ConvectionSourceError::Cancelled
        ))
    ));
    // Distinct views of one owner cannot silently carry different transport.
    struct Incoherent(ThermalFluidTransportState);
    impl FluidTransport for Incoherent {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            FluidTransportState::try_new(1.0, 1000.0, 0.02, 2e-5)
        }
    }
    impl ThermalFluidTransport for Incoherent {
        fn thermal_transport_state(
            &self,
        ) -> Result<ThermalFluidTransportState, fs_material::MaterialError> {
            Ok(self.0)
        }
    }
    assert!(
        evaluate_source_cell_convection(
            &heat,
            &[forced],
            1,
            |_, t| {
                Ok::<_, ()>(Incoherent(
                    ThermalFluidTransportState::try_new(
                        FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap(),
                        t.value(),
                        None,
                    )
                    .unwrap(),
                ))
            },
            || false
        )
        .is_err()
    );
    assert_eq!(heat.cells(), &saved);
}
