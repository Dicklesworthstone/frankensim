use super::*;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn ambient(density: f64) -> FireAmbientV1 {
    FireAmbientV1::new(
        Temperature::new(300.0),
        Density::new(density),
        SpecificHeatCapacityV1::new(1000.0),
        Acceleration::new(9.0),
    )
    .unwrap()
}

fn reference_fire() -> HeskestadFireV1 {
    HeskestadFireV1::new(
        Power::new(900_000.0),
        Power::new(112_500.0),
        Length::new(1.0),
        ambient(1.0),
    )
    .unwrap()
}

#[test]
fn heskestad_source_equations_and_density_similarity_g1_g3() {
    // NIST TN1889v1 Eqs3.15--3.17,3.20--3.21:
    // https://nvlpubs.nist.gov/nistpubs/TechnicalNotes/NIST.TN.1889v1.pdf
    // Test-owned dimensional case: Q*=1, S=Qc/(rho Cp T sqrt(g))=1/8,
    // L=2.68, z0=.38. At z-z0=8, the exact powers are S^(1/3)=1/2
    // and (z-z0)^(5/3)=32. These are equation checks, not experimental data.
    let fire = reference_fire();
    close(fire.dimensionless_heat_release(), 1.0, 2e-14);
    close(fire.mean_flame_height().value(), 2.68, 2e-14);
    close(fire.virtual_origin().unwrap().value(), 0.38, 2e-14);
    let entrainment = fire.entrainment(Length::new(8.38)).unwrap();
    assert_eq!(entrainment.value().region, HeskestadRegionV1::AboveFlame);
    // .196*1*3*(1/2)*32 + .196*2.9*(112500/(1000*300))
    close(entrainment.value().mass_rate.value(), 9.62115, 2e-13);
    let centerline = fire.centerline(Length::new(8.38)).unwrap();
    close(
        centerline.value().temperature_rise.value(),
        21.328125,
        2e-13,
    );
    close(centerline.value().temperature.value(), 321.328125, 2e-13);
    close(centerline.value().velocity.value(), 2.55, 2e-13);
    assert_eq!(
        entrainment.card().correlation,
        FireCorrelationIdV1::HeskestadEntrainment
    );
    assert_eq!(
        centerline.card().correlation,
        FireCorrelationIdV1::HeskestadCenterline
    );

    // Double actual ambient density and both HRRs. Q* and S remain fixed,
    // so only the entrained mass doubles. A hidden standard-air density
    // or a total/convective-HRR substitution breaks this dimensional law.
    let denser = HeskestadFireV1::new(
        Power::new(1_800_000.0),
        Power::new(225_000.0),
        Length::new(1.0),
        ambient(2.0),
    )
    .unwrap();
    close(denser.mean_flame_height().value(), 2.68, 2e-14);
    close(
        denser
            .entrainment(Length::new(8.38))
            .unwrap()
            .value()
            .mass_rate
            .value(),
        19.2423,
        2e-13,
    );
    let doubled = denser.centerline(Length::new(8.38)).unwrap().value();
    close(doubled.temperature_rise.value(), 21.328125, 2e-13);
    close(doubled.velocity.value(), 2.55, 2e-13);
    let farther = fire.centerline(Length::new(64.38)).unwrap().value();
    close(farther.temperature_rise.value(), 21.328125 / 32.0, 2e-13);
    close(farther.velocity.value(), 2.55 / 2.0, 2e-13);
}

#[test]
fn flaming_entrainment_joins_plume_and_centerline_refuses_inside_flame_g1() {
    let fire = reference_fire();
    let height = fire.mean_flame_height().value();
    let boundary = fire.entrainment(Length::new(height)).unwrap().value();
    let middle = fire.entrainment(Length::new(0.5 * height)).unwrap().value();
    let ground = fire.entrainment(Length::new(0.0)).unwrap().value();
    assert_eq!(boundary.region, HeskestadRegionV1::AboveFlame);
    assert_eq!(middle.region, HeskestadRegionV1::Flaming);
    assert_eq!(ground.region, HeskestadRegionV1::Flaming);
    close(
        middle.mass_rate.value(),
        0.5 * boundary.mass_rate.value(),
        2e-14,
    );
    assert_eq!(ground.mass_rate.value().to_bits(), 0.0_f64.to_bits());
    assert!(fire.centerline(Length::new(height)).is_ok());
    assert!(matches!(
        fire.centerline(Length::new(0.5 * height)),
        Err(FireCorrelationErrorV1::Regime {
            correlation: FireCorrelationIdV1::HeskestadCenterline,
            ..
        })
    ));
    // A source with zero convective fraction has no buoyant plume; its total
    // HRR still fixes the same fitted flame geometry.
    let radiation_only = HeskestadFireV1::new(
        Power::new(900_000.0),
        Power::new(0.0),
        Length::new(1.0),
        ambient(1.0),
    )
    .unwrap();
    assert_eq!(radiation_only.mean_flame_height(), fire.mean_flame_height());
    assert_eq!(
        radiation_only
            .entrainment(Length::new(8.38))
            .unwrap()
            .value()
            .mass_rate
            .value(),
        0.0
    );
    let no_convection = radiation_only
        .centerline(Length::new(8.38))
        .unwrap()
        .value();
    assert_eq!(no_convection.temperature_rise.value(), 0.0);
    assert_eq!(no_convection.velocity.value(), 0.0);
}

#[test]
fn dormant_source_is_exact_and_invalid_or_unrepresented_states_refuse_g1() {
    let dormant = HeskestadFireV1::new(
        Power::new(0.0),
        Power::new(0.0),
        Length::new(1.0),
        ambient(1.0),
    )
    .unwrap();
    assert_eq!(dormant.mean_flame_height().value(), 0.0);
    assert_eq!(dormant.dimensionless_heat_release(), 0.0);
    assert_eq!(dormant.virtual_origin(), None);
    for height in [0.0, 1.0, 10.0] {
        let mass = dormant.entrainment(Length::new(height)).unwrap().value();
        assert_eq!(mass.region, HeskestadRegionV1::Dormant);
        assert_eq!(mass.mass_rate.value().to_bits(), 0.0_f64.to_bits());
        let center = dormant.centerline(Length::new(height)).unwrap().value();
        assert_eq!(center.temperature.value(), 300.0);
        assert_eq!(center.temperature_rise.value().to_bits(), 0.0_f64.to_bits());
        assert_eq!(center.velocity.value().to_bits(), 0.0_f64.to_bits());
    }
    for (total, convective, diameter) in [
        (-1.0, 0.0, 1.0),
        (0.0, 1.0, 1.0),
        (900_000.0, -1.0, 1.0),
        (900_000.0, 900_001.0, 1.0),
        (900_000.0, 100_000.0, 0.0),
        (f64::INFINITY, 1.0, 1.0),
    ] {
        assert!(matches!(
            HeskestadFireV1::new(
                Power::new(total),
                Power::new(convective),
                Length::new(diameter),
                ambient(1.0)
            ),
            Err(FireCorrelationErrorV1::InvalidInput { .. })
        ));
    }
    assert!(matches!(
        HeskestadFireV1::new(
            Power::new(1.0),
            Power::new(0.7),
            Length::new(1.0),
            ambient(1.0)
        ),
        Err(FireCorrelationErrorV1::Regime {
            correlation: FireCorrelationIdV1::HeskestadFlame,
            ..
        })
    ));
    let payload = 0x7ff8_0000_0000_0042;
    assert_eq!(
        FireAmbientV1::new(
            Temperature::new(f64::from_bits(payload)),
            Density::new(1.0),
            SpecificHeatCapacityV1::new(1000.0),
            Acceleration::new(9.0)
        )
        .unwrap_err(),
        FireCorrelationErrorV1::InvalidInput {
            field: "ambient temperature K",
            bits: payload
        }
    );
    let fire = reference_fire();
    assert!(fire.entrainment(Length::new(-1.0)).is_err());
    assert!(fire.centerline(Length::new(f64::NAN)).is_err());
    assert!(matches!(
        fire.entrainment(Length::new(f64::MAX)),
        Err(FireCorrelationErrorV1::Arithmetic(_))
    ));
    assert!(matches!(
        fire.centerline(Length::new(f64::MAX)),
        Err(FireCorrelationErrorV1::Arithmetic(_))
    ));
}

#[test]
fn mqh_reproduces_published_nist_sp1169_table_2_1_g1() {
    // Independent published worked example, not an experimental validation:
    // https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.1169.pdf
    // Table2.1, Eqs2.1--2.4: room6x4x3m, vent.76x2.06m, Q200kW,
    // Tamb33C, k=.000463 kW/(m K), rho1607kg/m3, c=1kJ/(kg K).
    // All retained times precede thermal penetration through the .10m wall.
    // h=sqrt(k*rho*c/t) is calculated explicitly here in W/(m2 K).
    let opening = FireOpeningV1::new(Area::new(0.76 * 2.06), Length::new(2.06)).unwrap();
    let area = Area::new(106.4344);
    for (time, published_celsius) in [
        (60.0, 111.45),
        (120.0, 121.05),
        (180.0, 127.21),
        (600.0, 148.14),
        (1200.0, 162.24),
        (1800.0, 171.28),
        (2400.0, 178.07),
        (3000.0, 183.57),
    ] {
        let wall_h = (0.463_f64 * 1607.0 * 1000.0 / time).sqrt();
        let result = mqh_hot_layer_temperature(
            Power::new(200_000.0),
            Temperature::new(306.15),
            opening,
            area,
            HeatTransferCoefficient::new(wall_h),
        )
        .unwrap();
        let celsius = result.value().temperature.value() - 273.15;
        // 0.006K accommodates only the source table's 0.01K printing precision.
        // It is not a physical uncertainty or empirical acceptance envelope.
        assert!(
            (celsius - published_celsius).abs() < 0.006,
            "t={time}: actual={celsius}, published={published_celsius}"
        );
        assert_eq!(result.card().correlation, FireCorrelationIdV1::MqhHotLayer);
    }
    let zero = mqh_hot_layer_temperature(
        Power::new(0.0),
        Temperature::new(306.15),
        opening,
        area,
        HeatTransferCoefficient::new(20.0),
    )
    .unwrap()
    .value();
    assert_eq!(zero.temperature.value(), 306.15);
    assert_eq!(zero.temperature_rise.value().to_bits(), 0.0_f64.to_bits());
    assert!(
        mqh_hot_layer_temperature(
            Power::new(200_000.0),
            Temperature::new(306.15),
            opening,
            area,
            HeatTransferCoefficient::new(0.0)
        )
        .is_err()
    );
}

#[test]
fn flashover_and_ventilation_capacity_use_explicit_si_and_closed_limit_g1() {
    // NIST SP921 sections5.13--5.14 and Defining Flashover (1999).
    // Ao=2m2, Ho=4m -> Ao sqrt(Ho)=4m^(5/2), AT=100m2:
    // Thomas=(7.8*100+378*4)kW; Babrauskas=750*4kW;
    // air=.5*4kg/s; heat=.8*3MJ/kg*2kg/s.
    let opening = FireOpeningV1::new(Area::new(2.0), Length::new(4.0)).unwrap();
    let thomas = thomas_flashover_heat_release(Area::new(100.0), opening).unwrap();
    let babrauskas = babrauskas_flashover_heat_release(opening).unwrap();
    assert_eq!(thomas.value().value(), 2_292_000.0);
    assert_eq!(babrauskas.value().value(), 3_000_000.0);
    let capacity =
        kawagoe_ventilation_capacity(opening, SpecificEnergyV1::new(3_000_000.0), 0.8).unwrap();
    assert_eq!(capacity.value().air_mass_rate.value(), 2.0);
    assert_eq!(capacity.value().heat_release_capacity.value(), 4_800_000.0);
    let no_combustion =
        kawagoe_ventilation_capacity(opening, SpecificEnergyV1::new(3_000_000.0), 0.0)
            .unwrap()
            .value();
    assert_eq!(no_combustion.air_mass_rate.value(), 2.0);
    assert_eq!(no_combustion.heat_release_capacity.value(), 0.0);
    for efficiency in [-0.01, 1.01, f64::NAN] {
        assert!(
            kawagoe_ventilation_capacity(opening, SpecificEnergyV1::new(3_000_000.0), efficiency)
                .is_err()
        );
    }
    let closed = FireOpeningV1::new(Area::new(0.0), Length::new(4.0)).unwrap();
    let zero = kawagoe_ventilation_capacity(closed, SpecificEnergyV1::new(3_000_000.0), 0.8)
        .unwrap()
        .value();
    assert_eq!(zero.air_mass_rate.value().to_bits(), 0.0_f64.to_bits());
    assert_eq!(
        zero.heat_release_capacity.value().to_bits(),
        0.0_f64.to_bits()
    );
    assert!(thomas_flashover_heat_release(Area::new(100.0), closed).is_err());
    assert!(babrauskas_flashover_heat_release(closed).is_err());
    assert!(
        mqh_hot_layer_temperature(
            Power::new(200_000.0),
            Temperature::new(300.0),
            closed,
            Area::new(100.0),
            HeatTransferCoefficient::new(20.0)
        )
        .is_err()
    );
    assert!(FireOpeningV1::new(Area::new(1.0), Length::new(0.0)).is_err());
}

#[test]
fn all_fire_model_cards_preserve_nominal_status_and_unknown_discrepancy() {
    for id in FireCorrelationIdV1::ALL {
        let card = id.card();
        assert_eq!(card.correlation, id);
        assert_eq!(card.version, 1);
        assert_eq!(card.status, FirePredictionStatusV1::Nominal);
        assert_eq!(card.discrepancy, FireDiscrepancyV1::Unquantified);
        assert_eq!(card.quantitative_validity_domain, None);
        assert!(!card.source.citation.is_empty());
        assert!(!card.source.identifier.is_empty());
        assert!(card.source_url.starts_with("https://"));
        assert!(!card.domain_assumptions.is_empty());
        assert!(!card.known_failures.is_empty());
    }
}
