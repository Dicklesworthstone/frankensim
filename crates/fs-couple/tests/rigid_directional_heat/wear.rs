//! Actual path wear, original thermal hardness and atomic mechanical/heat state.
use super::*;
use fs_contact::interface_binding::wear::SourceArchardWear;
use fs_couple::friction_heat::directional::trajectory::apply_source_directional_trajectory_heat;
use fs_couple::friction_heat::directional::wear::*;
use fs_material::state_point::resolve_interface_state_point;
use fs_qty::Pressure;
use fs_tribo::WearState;
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;

const HARDNESS: [f64; 2] = [1e9, 2e9];

fn wear_close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}
fn hardness_key() -> PropertyKey {
    PropertyKey::new("indentation-hardness", Pressure::DIMS)
}

fn material_with_hardness(side: usize, limit: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = material(side);
    let mut claims = original.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: hardness_key(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, 1.2 * HARDNESS[side]), (400.0, 0.8 * HARDNESS[side])],
                dims: Pressure::DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                200.0,
                limit,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new("density", Density::DIMS);
    let chart = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &rho,
            temperature_axis: "T",
            lower: &material_query(200.0),
            upper: &material_query(400.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap();
    (card, chart)
}

fn wear_interface(cards: &[MaterialCard; 2], k: [f64; 2]) -> InterfaceSystemCard {
    let original = interface(cards, [scalar(0.2), scalar(0.4)], 400.0);
    let mut claims = original.claims().clone();
    for side in 0..2 {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(["wear-a", "wear-b"][side], Dims::NONE),
                value: curve("path", LENGTH, [(0.0, k[side]), (10.0, 2.0 * k[side])]),
                validity: original.claims_for("kinetic-x")[0].1.validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

fn wear_state(
    law: &SourceThermalEllipticCoulomb,
    cards: &[MaterialCard; 2],
    card: &InterfaceSystemCard,
) -> SourceArchardWear {
    let keys = ["wear-a", "wear-b"].map(|name| PropertyKey::new(name, Dims::NONE));
    let requirements = keys.each_ref().map(|k| {
        ScalarPropertyRequirement::try_with_key(k, ScalarAdmissibility::NonNegative).unwrap()
    });
    let coefficients = resolve_interface_state_point(
        card,
        &law.principal_sources()[0].initial_source().query(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let h = ScalarPropertyRequirement::try_with_key(
        &hardness_key(),
        ScalarAdmissibility::StrictlyPositive,
    )
    .unwrap();
    let hardness = cards.each_ref().map(|card| {
        resolve_material_state_point(
            card,
            &material_query(300.0),
            std::slice::from_ref(&h),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    });
    SourceArchardWear::new(
        law,
        coefficients,
        keys,
        hardness,
        [hardness_key(), hardness_key()],
        ["T".into(), "T".into()],
        [WearState::default(); 2],
    )
    .unwrap()
}
fn heat_controls() -> DirectionalHeatSolve {
    let mut solve = controls();
    solve.motion.distance_bracket_m = [0.0, 2.0];
    solve
}
fn trajectory_controls() -> EllipticRelaxationControls {
    EllipticRelaxationControls {
        time_tolerance_s: 1e-10,
        travel_tolerance_m: 1e-10,
        maximum_evaluations: 40_000,
        maximum_root_iterations: 96,
    }
}

#[test]
fn curved_path_wears_both_original_materials_at_actual_temperature_without_changing_heat() {
    let (a, ca) = material_with_hardness(0, 400.0);
    let (b, cb) = material_with_hardness(1, 400.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let k = [1e-6, 3e-6];
    let card = wear_interface(&cards, k);
    let law = law(&card, &q, [0.6, 0.8]);
    let mut state = law.initial_state();
    let mut worn = wear_state(&law, &cards, &card);
    let original_wear = worn.clone();
    let mut heat = thermal(&curves, [290.0, 310.0]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = mechanics(true);
    let spare = bodies[2];
    let mut plain_bodies = bodies;
    let mut plain_state = state.clone();
    let mut plain_heat = heat.clone();
    let plain = apply_source_directional_trajectory_heat(
        &mut plain_bodies,
        &props,
        &mut plain_state,
        &law,
        contact,
        &mut plain_heat,
        &route,
        binding(&cards, &curves, &card, &q),
        2.0,
        heat_controls(),
        trajectory_controls(),
        &[0.0; 2],
        305.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    let result = apply_source_directional_trajectory_heat_with_wear(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        &mut worn,
        2.0,
        heat_controls(),
        trajectory_controls(),
        &[0.0; 2],
        305.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    let report = result.phase;
    assert!(
        report.actual_distance_m - fs_math::det::hypot(report.net_slip_m[0], report.net_slip_m[1])
            > 1e-4
    );
    assert_eq!(bodies, plain_bodies);
    assert_eq!(state, plain_state);
    assert_eq!(bodies[2], spare);
    assert_eq!(heat.cells(), plain_heat.cells());
    assert_eq!(heat.time_s(), plain_heat.time_s());
    assert_eq!(report.heat.dissipated_j, plain.heat.dissipated_j);
    assert_eq!(report.actual_distance_m, plain.actual_distance_m);
    let path = 0.25 + report.actual_distance_m;
    for i in 0..2 {
        let temperature = report.phases[i].temperature_k();
        let expected_h = HARDNESS[i] * (1.0 - 0.002 * (temperature - 300.0));
        let expected_k = k[i] * (1.0 + 0.1 * path);
        let expected = expected_k * 10.0 * report.actual_distance_m / expected_h;
        wear_close(result.wear[i].increment_m3, expected);
        wear_close(worn.volumes_m3()[i], expected);
        close(
            worn.hardness_sources()[i].state_coordinate("T").unwrap(),
            temperature,
        );
        assert_eq!(
            worn.hardness_sources()[i].properties()[0]
                .answer()
                .receipt
                .selected,
            original_wear.hardness_sources()[i].properties()[0]
                .answer()
                .receipt
                .selected
        );
        assert_eq!(
            worn.coefficient_source().properties()[i]
                .answer()
                .receipt
                .selected,
            original_wear.coefficient_source().properties()[i]
                .answer()
                .receipt
                .selected
        );
    }
    close(report.combined_energy_residual_j, 0.0);
}

#[test]
fn reversal_accumulates_unsigned_wear_then_arrest_adds_none_and_zero_k_stays_zero() {
    let (a, ca) = material_with_hardness(0, 400.0);
    let (b, cb) = material_with_hardness(1, 400.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = wear_interface(&cards, [1e-6, 0.0]);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let mut worn = wear_state(&law, &cards, &card);
    let mut heat = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (original, props, contact) = mechanics(false);
    let mut bodies = original;
    let mut total = 0.0;
    for phase in 0..3 {
        if phase == 1 {
            // A separately prescribed drive restarts the original pair in reverse.
            for i in 0..2 {
                bodies[i] = RigidBodyState::new(
                    original[i].pose(),
                    original[i].linear_momentum_world().scale(-1.0),
                    original[i].angular_momentum_body().scale(-1.0),
                )
                .unwrap();
            }
        }
        let before = worn.volumes_m3();
        // Warm-start from the actual accepted thermal state. At equilibrium,
        // a stale temperature guess can ask for unrepresentable tiny cell heat.
        let temperature_guess =
            0.5 * (heat.temperature_k(CELLS[0]).unwrap() + heat.temperature_k(CELLS[1]).unwrap());
        let result = apply_source_directional_trajectory_heat_with_wear(
            &mut bodies,
            &props,
            &mut state,
            &law,
            contact,
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            &mut worn,
            0.8,
            heat_controls(),
            trajectory_controls(),
            &[0.0; 2],
            temperature_guess,
            |_, p| transport(&cards, p),
            &CancelGate::new(),
        )
        .unwrap_or_else(|error| panic!("wear phase {phase}: {error:?}"));
        let ds = result.phase.actual_distance_m;
        if phase < 2 {
            close(ds, 30.0 / 101.0);
            assert!(result.phase.net_slip_m[0] * if phase == 0 { 1.0 } else { -1.0 } > 0.0);
            assert!(result.wear[0].increment_m3 > 0.0);
            assert!(worn.volumes_m3()[0] > before[0]);
        } else {
            // Recoil can leave a roundoff contact velocity after finite arrest.
            // Keep its native quadratic path; do not invent a wear deadband.
            assert!(ds <= 64.0 * f64::EPSILON * f64::EPSILON);
            assert_eq!(worn.volumes_m3(), before);
        }
        total += result.wear[0].increment_m3;
        wear_close(worn.volumes_m3()[0], total);
        assert_eq!(result.wear[1].increment_m3, 0.0);
        assert_eq!(worn.volumes_m3()[1], 0.0);
        close(result.phase.combined_energy_residual_j, 0.0);
    }
    close(heat.time_s(), 2.4);
}

#[test]
fn late_original_hardness_domain_failure_rolls_back_bodies_history_heat_and_both_wear_volumes() {
    let (a, ca) = material_with_hardness(0, 400.0);
    let (b, cb) = material_with_hardness(1, 300.1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = wear_interface(&cards, [1e-6, 3e-6]);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let mut worn = wear_state(&law, &cards, &card);
    let mut heat = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = mechanics(false);
    let old_bodies = bodies;
    let old_state = state.clone();
    let old_heat = heat.clone();
    let old_wear = worn.clone();
    let calls = std::cell::Cell::new(0);
    let result = apply_source_directional_trajectory_heat_with_wear(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        &mut worn,
        0.8,
        heat_controls(),
        trajectory_controls(),
        &[0.0; 2],
        305.0,
        |_, p| {
            calls.set(calls.get() + 1);
            transport(&cards, p)
        },
        &CancelGate::new(),
    );
    assert!(matches!(result, Err(DirectionalWearHeatError::Wear(_))));
    assert!(calls.get() > 0);
    assert_eq!(bodies, old_bodies);
    assert_eq!(state, old_state);
    assert_eq!(worn, old_wear);
    assert_eq!(heat.time_s(), old_heat.time_s());
    assert_eq!(heat.cells(), old_heat.cells());
    assert_eq!(heat.links(), old_heat.links());
    assert_eq!(heat.energy_j(), old_heat.energy_j());
}

#[test]
fn source_wear_refuses_swapped_roles_changed_owner_and_evolution_of_fixed_conditions() {
    let (a, ca) = material_with_hardness(0, 400.0);
    let (b, cb) = material_with_hardness(1, 400.0);
    let cards = [a, b];
    let q = query()
        .with_quantity("ambient-pressure", static_pressure(), 1e5)
        .unwrap();
    let card = wear_interface(&cards, [1e-6, 3e-6]);
    let friction = law(&card, &q, [1.0, 0.0]);
    let state = wear_state(&friction, &cards, &card);
    let keys = ["wear-a", "wear-b"].map(|name| PropertyKey::new(name, Dims::NONE));
    let mut swapped = state.hardness_sources().clone();
    swapped.swap(0, 1);
    assert!(
        SourceArchardWear::new(
            &friction,
            state.coefficient_source().clone(),
            keys,
            swapped,
            [hardness_key(), hardness_key()],
            ["T".into(), "T".into()],
            [WearState::default(); 2]
        )
        .is_err()
    );
    let changed = wear_interface(&cards, [2e-6, 3e-6]);
    let other = law(&changed, &q, [1.0, 0.0]);
    assert!(state.admit_friction(&other).is_err());
    let phases = [ca, cb]
        .each_ref()
        .map(|curve| curve.curve().state_at_specific_enthalpy(100.0).unwrap());
    assert!(
        state
            .advance_candidate(&friction, &q, [&cards[0], &cards[1]], phases, 10.0, 0.1)
            .is_ok()
    );
    let changed_conditions = q
        .with_quantity("ambient-pressure", static_pressure(), 2e5)
        .unwrap();
    assert!(
        state
            .advance_candidate(
                &friction,
                &changed_conditions,
                [&cards[0], &cards[1]],
                phases,
                10.0,
                0.1
            )
            .is_err()
    );
}
