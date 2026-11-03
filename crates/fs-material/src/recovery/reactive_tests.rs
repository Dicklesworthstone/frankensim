use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::{Dims, QuantitySpec};

fn card(rate: f64) -> LocalTransferCard {
    let description =
        format!("authored transient local rate {rate} /s; fixed pH 9..11; no calibration");
    let artifact = hash_domain("reactive-cascade-test", description.as_bytes());
    LocalTransferCard::from_card(ConstitutiveModelCard {
        law: LawId(TRANSIENT_TRANSFER_LAW.into()),
        law_version: 1,
        parameters: [(
            TRANSFER_RATE_PARAMETER.into(),
            LawParameter {
                value: rate,
                dims: Frequency::DIMS,
            },
        )]
        .into(),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "pH",
            QuantitySpec::dimensional(Dims::NONE),
            9.0,
            11.0,
        ),
        sources: vec![artifact],
        provenance: Provenance {
            source: description,
            license: "authored test".into(),
            artifact: Some(artifact),
        },
    })
    .unwrap()
}
fn inventory_pair(u: f64, d: f64) -> ConstituentInventory {
    ConstituentInventory {
        undissolved: Mass::new(u),
        dissolved: Mass::new(d),
    }
}
fn train(rates: &[f64], volume: f64, inventories: &[(f64, f64)]) -> ReactiveTankCascade {
    ReactiveTankCascade::new(
        "Cu",
        Volume::new(volume),
        rates.iter().map(|r| card(*r)).collect(),
        inventories
            .iter()
            .map(|(u, d)| inventory_pair(*u, *d))
            .collect(),
    )
    .unwrap()
}
fn points(n: usize) -> Vec<QueryPoint> {
    vec![
        QueryPoint::new()
            .with_quantity("pH", QuantitySpec::dimensional(Dims::NONE), 10.0)
            .unwrap();
        n
    ]
}
fn advance(train: &mut ReactiveTankCascade, t: f64, q: f64, u: f64, d: f64) -> ReactiveInterval {
    let point = points(train.tanks().len());
    train
        .advance(
            Time::new(t),
            VolumetricFlowRate::new(q),
            Density::new(u),
            Density::new(d),
            &point,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn balance(interval: &ReactiveInterval) {
    close(interval.undissolved_residual.value(), 0.0);
    close(interval.dissolved_residual.value(), 0.0);
    close(
        interval.final_inventory.undissolved.value()
            + interval.final_inventory.dissolved.value()
            + interval.outlet.undissolved.value()
            + interval.outlet.dissolved.value(),
        interval.initial.undissolved.value()
            + interval.initial.dissolved.value()
            + interval.inlet.undissolved.value()
            + interval.inlet.dissolved.value(),
    );
}

#[test]
fn one_tank_matches_independently_integrated_linear_balances() {
    let mut model = train(&[0.25], 2.0, &[(2.0, 1.0)]);
    let result = advance(&mut model, 4.0, 0.5, 5.0, 1.0);
    // a=Q/V=.25, a+k=.5; Uss=5 kg, total steady inventory=12 kg.
    let u = 5.0 - 3.0 * (-2.0_f64).exp();
    let total = 12.0 - 9.0 * (-1.0_f64).exp();
    let integral_u = 20.0 - 6.0 * (1.0 - (-2.0_f64).exp());
    let integral_total = 48.0 - 36.0 * (1.0 - (-1.0_f64).exp());
    close(model.tanks()[0].undissolved.value(), u);
    close(model.tanks()[0].dissolved.value(), total - u);
    close(result.outlet.undissolved.value(), 0.25 * integral_u);
    close(
        result.outlet.dissolved.value(),
        0.25 * (integral_total - integral_u),
    );
    close(result.converted.value(), 0.25 * integral_u);
    assert_eq!(result.inlet, inventory_pair(10.0, 2.0));
    balance(&result);
    assert_eq!(model.constituent(), "Cu");
    assert_eq!(model.laws()[0].source(), card(0.25).source());
}

#[test]
fn stopped_flow_converts_inventory_without_exporting_or_destroying_constituent() {
    let mut model = train(&[0.1, 0.2], 4.0, &[(10.0, 1.0), (20.0, 2.0)]);
    let result = advance(&mut model, 3.0, 0.0, 999.0, 999.0);
    let remaining = 10.0 * (-0.3_f64).exp() + 20.0 * (-0.6_f64).exp();
    close(result.final_inventory.undissolved.value(), remaining);
    close(result.final_inventory.dissolved.value(), 33.0 - remaining);
    close(result.converted.value(), 30.0 - remaining);
    assert_eq!(result.outlet, inventory_pair(0.0, 0.0));
    assert_eq!(result.inlet, inventory_pair(0.0, 0.0));
    balance(&result);
    let restart = advance(&mut model, 0.2, 1.0, 0.0, 0.0);
    assert!(restart.outlet.dissolved.value() > 0.0);
    assert!(restart.outlet.undissolved.value() > 0.0);
    balance(&restart);
}

#[test]
fn zero_rates_reproduce_each_phase_of_the_real_hydraulic_model() {
    use super::super::residence::StirredTankCascade;
    let mut model = train(&[0.0; 3], 6.0, &[(2.0, 0.2), (0.0, 0.3), (1.0, 0.4)]);
    let source = card(0.0).source().provenance.clone();
    let mut u = StirredTankCascade::new(
        Volume::new(6.0),
        vec![Density::new(1.0), Density::new(0.0), Density::new(0.5)],
        source.clone(),
    )
    .unwrap();
    let mut d = StirredTankCascade::new(
        Volume::new(6.0),
        vec![Density::new(0.1), Density::new(0.15), Density::new(0.2)],
        source,
    )
    .unwrap();
    let result = advance(&mut model, 1.7, 0.4, 4.0, 5.0);
    let ur = u
        .advance(
            Time::new(1.7),
            VolumetricFlowRate::new(0.4),
            Density::new(4.0),
        )
        .unwrap();
    let dr = d
        .advance(
            Time::new(1.7),
            VolumetricFlowRate::new(0.4),
            Density::new(5.0),
        )
        .unwrap();
    close(result.outlet.undissolved.value(), ur.outlet_mass.value());
    close(result.outlet.dissolved.value(), dr.outlet_mass.value());
    for (i, tank) in model.tanks().iter().enumerate() {
        close(
            tank.undissolved.value(),
            2.0 * u.tank_concentrations()[i].value(),
        );
        close(
            tank.dissolved.value(),
            2.0 * d.tank_concentrations()[i].value(),
        );
    }
    assert_eq!(result.converted.value(), 0.0);
    balance(&result);
}

#[test]
fn unequal_local_rates_retain_the_independent_sequential_steady_solution() {
    // a=.5/s and feed U=2 kg/m3, D=.4 kg/m3. Total per tank=4.8 kg.
    let us = [8.0 / 3.0, 4.0 / 3.0, 4.0 / 9.0];
    let mut model = train(&[0.25, 0.5, 1.0], 6.0, &us.map(|u| (u, 4.8 - u)));
    let before = model.tanks().to_vec();
    let result = advance(&mut model, 20.0, 1.0, 2.0, 0.4);
    for (actual, expected) in model.tanks().iter().zip(before) {
        close(actual.undissolved.value(), expected.undissolved.value());
        close(actual.dissolved.value(), expected.dissolved.value());
    }
    close(result.outlet.undissolved.value(), 40.0 / 9.0);
    close(result.outlet.dissolved.value(), 48.0 - 40.0 / 9.0);
    close(result.converted.value(), 40.0 - 40.0 / 9.0);
    balance(&result);
}

#[test]
fn feed_pulse_flow_change_and_split_intervals_preserve_phase_balances() {
    let mut one = train(&[0.1, 0.3, 0.0, 0.2], 8.0, &[(0.0, 0.0); 4]);
    let mut split = one.clone();
    for (duration, flow, u, d) in [
        (5.0, 0.5, 3.0, 0.1),
        (3.0, 0.0, 0.0, 0.0),
        (20.0, 0.8, 0.0, 0.0),
    ] {
        let whole = advance(&mut one, duration, flow, u, d);
        let a = advance(&mut split, duration / 2.0, flow, u, d);
        let b = advance(&mut split, duration / 2.0, flow, u, d);
        close(
            whole.outlet.undissolved.value(),
            a.outlet.undissolved.value() + b.outlet.undissolved.value(),
        );
        close(
            whole.outlet.dissolved.value(),
            a.outlet.dissolved.value() + b.outlet.dissolved.value(),
        );
        close(
            whole.converted.value(),
            a.converted.value() + b.converted.value(),
        );
        for (left, right) in one.tanks().iter().zip(split.tanks()) {
            close(left.undissolved.value(), right.undissolved.value());
            close(left.dissolved.value(), right.dissolved.value());
        }
        balance(&whole);
        balance(&a);
        balance(&b);
    }
}

#[test]
fn source_validity_numeric_failures_and_work_limits_leave_state_unchanged() {
    let mut model = train(&[0.25], 2.0, &[(2.0, 1.0)]);
    let before = model.clone();
    for (dt, q, u, d) in [
        (f64::NAN, 1.0, 1.0, 0.0),
        (1.0, -1.0, 1.0, 0.0),
        (1.0, 1.0, -1.0, 0.0),
        (1.0, 1.0, 1.0, f64::INFINITY),
        (10000.0, 1.0, 1.0, 0.0),
        (1.0, f64::MAX, 1.0, 0.0),
    ] {
        assert!(
            model
                .advance(
                    Time::new(dt),
                    VolumetricFlowRate::new(q),
                    Density::new(u),
                    Density::new(d),
                    &points(1)
                )
                .is_err()
        );
        assert_eq!(model, before);
    }
    assert!(
        model
            .advance(
                Time::new(1.0),
                VolumetricFlowRate::new(1.0),
                Density::new(1.0),
                Density::new(0.0),
                &[QueryPoint::new()]
            )
            .is_err()
    );
    assert_eq!(model, before);
    let zero = advance(&mut model, 0.0, 0.0, 0.0, 0.0);
    assert_eq!(zero.substeps, 0);
    assert_eq!(model, before);
}

#[test]
fn source_scaling_retains_tiny_exposure_and_underflow_refuses_transactionally() {
    let mut tiny = train(&[1.0], 1.0, &[(1e200, 0.0)]);
    let result = advance(&mut tiny, 1e-200, 0.0, 0.0, 0.0);
    close(result.converted.value(), 1.0);
    close(result.final_inventory.dissolved.value(), 1.0);
    assert_eq!(result.outlet, inventory_pair(0.0, 0.0));

    let mut exhausted = train(&[1.0], 1.0, &[(1.0, 0.0)]);
    let original = exhausted.clone();
    assert!(matches!(
        exhausted.advance(
            Time::new(1000.0),
            VolumetricFlowRate::new(0.0),
            Density::new(0.0),
            Density::new(0.0),
            &points(1)
        ),
        Err(ReactiveError::Arithmetic(_))
    ));
    assert_eq!(exhausted, original);
}

#[test]
fn stationary_and_batch_sources_cannot_authorize_transient_transfer() {
    for law in [
        super::super::SINGLE_RATE_LAW,
        super::super::continuous::CONTINUOUS_SINGLE_RATE_LAW,
    ] {
        let mut source = card(0.1).source().clone();
        source.law = LawId(law.into());
        assert!(LocalTransferCard::from_card(source).is_err());
    }
    let mut source = card(0.1).source().clone();
    source
        .parameters
        .get_mut(TRANSFER_RATE_PARAMETER)
        .unwrap()
        .dims = Dims::NONE;
    assert!(LocalTransferCard::from_card(source).is_err());
    let mut source = card(0.1).source().clone();
    source.validity = source.validity.with_quantity(
        "elapsed-time",
        fs_qty::QuantitySpec::dimensional(Time::DIMS),
        0.0,
        1.0,
    );
    assert!(LocalTransferCard::from_card(source).is_err());
    assert!(
        ReactiveTankCascade::new(
            "",
            Volume::new(1.0),
            vec![card(0.1)],
            vec![inventory_pair(0.0, 0.0)]
        )
        .is_err()
    );
    assert!(
        ReactiveTankCascade::new(
            "Cu",
            Volume::new(1.0),
            vec![card(0.1)],
            vec![inventory_pair(-1.0, 0.0)]
        )
        .is_err()
    );
}
