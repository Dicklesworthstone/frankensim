//! Independent material calorimetry through the original gas/network owners.
use super::*;
use crate::multicomponent::network::dissipation::tests::owned;
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallCaloricConstraint, exchange_cell, exchange_cell_with_heat_storage,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use std::sync::Arc;

// Synthetic source C(T)=a+bT J/K, using one explicit kilogram. This is not a
// measured material or a claimed phase-stability interval.
fn source(a: f64, b: f64, low: f64, high: f64) -> MaterialWallCaloric<'static> {
    let key = PropertyKey::new("constant_volume_specific_heat", Dims([2, 0, -2, -1, 0, 0]));
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(low, a + b * low), (high, a + b * high)],
                dims: key.dims(),
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic wall caloric fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "stationary solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |temperature| {
        QueryPoint::new()
            .with_quantity(
                "T",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                temperature,
            )
            .unwrap()
    };
    let storage = PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
        id: 0,
        card: &card,
        capacity: &key,
        amount: HeatStorageAmount::MassKg(1.0),
        temperature_axis: "T",
        lower: &point(low),
        upper: &point(high),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }])
    .unwrap();
    MaterialWallCaloric::new(
        Arc::new(storage),
        WallCaloricConstraint::IsochoricInternalEnergy,
    )
}

#[test]
fn source_capacity_integral_matches_quadratic_reference_and_preserves_drag_heat() {
    let g = gas(1, false, 1e12);
    let q = g.conserved(&[1.0], 1.0, 2.0).unwrap();
    let (a, b, volume, dt, conductance) = (1.0, 0.2, 0.2, 0.2, 0.3);
    // Initial gas T=2 lies OUTSIDE this wall chart. The valid root is inside;
    // checking unsupported arbitrary trial endpoints must not reject it.
    let law = source(a, b, 4.8, 5.2);
    let contact = law
        .contact(5.0, conductance, WallDrag::Linear { rate_per_s: 2.0 })
        .unwrap();
    let update = exchange_cell_with_heat_storage(
        &g,
        q.as_cell(),
        volume,
        Time::new(dt),
        contact,
        &law,
        controls(),
        || false,
    )
    .unwrap();
    let gas_t = 2.0 + (1.0 - (-4.0_f64 * dt).exp()) / 5.0;
    let gas_capacity = 0.5;
    let transfer = dt * conductance / (1.0 + dt * conductance / gas_capacity);
    // Cwall=a+b*T: .5*b*Tw²+(a+transfer)*Tw
    // = a*5 + .5*b*25 + transfer*Tg_after_drag.
    let constant = a * 5.0 + 0.5 * b * 25.0 + transfer * gas_t;
    let linear = a + transfer;
    let expected_wall = 2.0 * constant / (linear + (linear * linear + 2.0 * b * constant).sqrt());
    let expected_heat = transfer * (expected_wall - gas_t);
    let actual = crate::multicomponent::Cell {
        partial_densities: &q.partial_densities,
        momentum: update.momentum,
        energy: update.energy,
    };
    close(wall_temperature(update.contact), expected_wall, 2e-11);
    close(
        g.primitive(actual, || false).unwrap().temperature,
        gas_t + expected_heat / gas_capacity,
        2e-11,
    );
    close(update.report.heat_to_gas_j, expected_heat, 2e-12);
    let wall_change = a * (expected_wall - 5.0) + 0.5 * b * (expected_wall * expected_wall - 25.0);
    close(update.report.wall_energy_change_j, wall_change, 2e-11);
    close(
        update.report.drag_heating_j,
        0.5 * volume * (1.0 - (-4.0_f64 * dt).exp()),
        2e-12,
    );
    close(update.report.energy_residual_j, 0.0, 2e-11);
    assert_eq!(update.report.external_heat_to_gas_j, 0.0);
    assert_eq!(
        law.constraint(),
        WallCaloricConstraint::IsochoricInternalEnergy
    );
    // Constant-source parity uses the same native gas drag and scalar solve.
    let constant = source(1.0, 0.0, 1.0, 10.0);
    let contact = constant.contact(5.0, conductance, contact.drag).unwrap();
    let native = exchange_cell(
        &g,
        q.as_cell(),
        volume,
        Time::new(dt),
        contact,
        controls(),
        || false,
    )
    .unwrap();
    let bound = exchange_cell_with_heat_storage(
        &g,
        q.as_cell(),
        volume,
        Time::new(dt),
        contact,
        &constant,
        controls(),
        || false,
    )
    .unwrap();
    close(native.energy, bound.energy, 2e-11);
    close(
        wall_temperature(native.contact),
        wall_temperature(bound.contact),
        2e-11,
    );
}

#[test]
fn material_substitution_changes_actual_network_and_keeps_one_wall_state() {
    let g = gas(1, false, 0.0);
    let mut outcomes = Vec::new();
    let mut identities = Vec::new();
    for slope in [0.0, 0.5] {
        let law = source(1.0, slope, 1.0, 10.0);
        identities.push(law.storage().unwrap().identity());
        let mut n = network(&g, &[0.0], true);
        let initial = n.totals().absolute_energy(&g).unwrap();
        let contact = law.contact(5.0, 0.3, WallDrag::None).unwrap();
        let mut wall = uniform_operator(&n, contact)
            .with_material_heat_storage(
                vec![vec![Some(law.clone()); 2]],
                vec![Some(law.clone())],
                || false,
            )
            .unwrap();
        assert!(core::ptr::eq(
            wall.branch_heat_storage(0).unwrap()[0]
                .as_ref()
                .unwrap()
                .storage()
                .unwrap(),
            law.storage().unwrap()
        ));
        assert_eq!(wall.plenum_heat_storage().len(), 1);
        let mut reported = 0.0;
        for _ in 0..3 {
            let report = wall
                .advance_for(&mut n, Time::new(0.2), controls(), || false)
                .unwrap();
            reported += report.wall_energy_change_j;
            assert_eq!(report.external_heat_to_gas_j, 0.0);
        }
        let mut actual_wall_change = 0.0;
        for contact in wall
            .branch_contacts(0)
            .unwrap()
            .iter()
            .chain(wall.plenum_contacts())
        {
            let t = wall_temperature(*contact);
            actual_wall_change += (t - 5.0) + 0.5 * slope * (t * t - 25.0);
            if let WallThermal::Lumped {
                heat_capacity_j_per_k,
                ..
            } = contact.thermal
            {
                close(heat_capacity_j_per_k, 1.0 + slope * t, 2e-12);
            }
        }
        close(actual_wall_change, reported, 3e-11);
        close(
            n.totals().absolute_energy(&g).unwrap() - initial + actual_wall_change,
            0.0,
            4e-11,
        );
        outcomes.push(
            g.primitive(n.branches()[0].cell(0).unwrap(), || false)
                .unwrap()
                .temperature,
        );
    }
    assert_ne!(identities[0], identities[1]);
    assert!(outcomes[1] > outcomes[0] + 1e-3);
}

#[test]
fn reacting_wall_domain_exit_and_late_cancellation_preserve_all_accepted_state() {
    use crate::chemistry::IsochoricKineticsSolve;
    use crate::multicomponent::network::reactive::{ReactingNetwork1d, ReactiveNetworkError};
    use crate::multicomponent::reactive::ReactiveStepControl;
    use crate::multicomponent::tests::{bank, source as gas_source};
    use crate::multicomponent::{CaloricDomain, Conserved};
    use fs_thermochem::kinetics::{
        ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
    };
    use fs_thermochem::{
        MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    };
    let make = |upper| {
        let thermodynamics = bank(vec![
            gas_source("A", R, 3.5, 0.0, 0.0),
            gas_source("B", R, 3.5, 0.0, -10.0),
        ]);
        let reaction = ElementaryReactionV1::new(
            ReactionId::new("wall-heat-isomers").unwrap(),
            vec![(SpeciesId::new("A").unwrap(), 1)],
            vec![(SpeciesId::new("B").unwrap(), 1)],
            ArrheniusRateV1::new(
                Temperature::new(1.0),
                ReactionProgressRateQuantityV1::new(100.0),
                0.0,
                MolarEnergyQuantityV1::new(0.0),
            )
            .unwrap(),
        )
        .unwrap();
        let mechanism = IdealGasKineticsV1::new(
            thermodynamics.clone(),
            "wall-heat-isomers",
            1,
            vec![reaction],
            1e-12,
        )
        .unwrap();
        let gas = MixtureGas::new(
            thermodynamics,
            CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
            || false,
        )
        .unwrap();
        let branch = MixtureEuler1d::new(
            gas.clone(),
            Length::new(0.5),
            0.4,
            vec![gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap(); 2],
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let n = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
        let law = source(1.0, 0.0, 1.9, upper);
        let wall = uniform_operator(&n, law.contact(2.0, 100.0, WallDrag::None).unwrap())
            .with_material_heat_storage(vec![vec![Some(law.clone()); 2]], vec![], || false)
            .unwrap();
        let chemistry = IsochoricKineticsSolve {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            energy_tolerance_j: 1e-10,
            initial_step_s: 1e-5,
            ..IsochoricKineticsSolve::default()
        };
        (
            ReactingNetwork1d::new(n, mechanism, chemistry, || false).unwrap(),
            wall,
        )
    };
    let states = |flow: &ReactingNetwork1d| -> Vec<Conserved> {
        (0..2)
            .map(|i| owned(flow.branches()[0].cell(i).unwrap()))
            .collect()
    };
    let step = ReactiveStepControl {
        maximum_coupling_retries: 0,
        ..ReactiveStepControl::default()
    };
    let (mut valid, mut valid_wall) = make(10.0);
    let mut polls = 0;
    valid
        .step_with_wall(&mut valid_wall, Time::new(1e-3), step, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(valid.branches()[0].cell(0).unwrap().partial_densities[1] > 0.0);
    assert!(wall_temperature(valid_wall.branch_contacts(0).unwrap()[0]) > 2.00001);
    for domain_exit in [true, false] {
        let (mut flow, mut wall) = make(if domain_exit { 2.00001 } else { 10.0 });
        let before = states(&flow);
        let contacts = wall.branch_contacts(0).unwrap().to_vec();
        let totals = flow.totals();
        let mut visits = 0;
        let error = flow
            .step_with_wall(&mut wall, Time::new(1e-3), step, controls(), || {
                visits += 1;
                !domain_exit && visits >= polls
            })
            .unwrap_err();
        fn wall_refusal(e: &ReactiveNetworkError) -> bool {
            match e {
                ReactiveNetworkError::Wall(_) => true,
                ReactiveNetworkError::RetryLimit { last, .. } => wall_refusal(last),
                _ => false,
            }
        }
        if domain_exit {
            assert!(wall_refusal(&error), "{error:?}");
        } else {
            assert_eq!(error, ReactiveNetworkError::Cancelled);
        }
        assert_eq!(states(&flow), before);
        assert_eq!(flow.totals(), totals);
        assert_eq!(flow.time(), Time::new(0.0));
        assert_eq!(flow.last_report(), None);
        assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
        assert_eq!(wall.last_report(), None);
    }
}
