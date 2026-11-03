//! G0/G1/G3/G4: finite gas, neutral acid/base mass action and shared pH.
use super::*;
use crate::stoichiometry::thermodynamics::{
    acid_base::gas_liquid::VolatileAcidBaseSystem, gas_liquid::IdealGasLiquidPartition,
};
use fs_qty::Volume;

const AG: &str = "authored-HA(g)";
const BH: &str = "authored-BH[+](aq)";
const B: &str = "authored-B(aq)";
const BG: &str = "authored-B(g)";
const RT: f64 = crate::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE;

fn append(fixture: &mut Fixture, definition: Definition, phase: &str) {
    let at = point(phase, definition.charge != 0, fixture.reference);
    fixture
        .packs
        .push(pack(definition.name, definition.gibbs, &at, None));
    fixture.points.push(at);
    fixture.definitions.push(definition);
}

fn add_base(fixture: &mut Fixture, ka: f64) -> AcidFamily {
    append(
        fixture,
        Definition {
            name: BH,
            atoms: [4, 0, 1, 0],
            charge: 1,
            gibbs: 0.0,
        },
        "source-phase-aqueous",
    );
    append(
        fixture,
        Definition {
            name: B,
            atoms: [3, 0, 1, 0],
            charge: 0,
            gibbs: -RT * ka.ln(),
        },
        "source-phase-aqueous",
    );
    AcidFamily::new(
        "base",
        vec![fixture.step("authored-base-dissociation", BH, B).unwrap()],
    )
    .unwrap()
}

fn partition(
    fixture: &mut Fixture,
    gas: &'static str,
    neutral: &'static str,
    log_k: f64,
) -> IdealGasLiquidPartition {
    let definition = &fixture.definitions[fixture.index(neutral)];
    let atoms = definition.atoms;
    let gibbs = definition.gibbs + RT * log_k;
    append(
        fixture,
        Definition {
            name: gas,
            atoms,
            charge: 0,
            gibbs,
        },
        "source-phase-gas",
    );
    let equation = fixture.reaction(gas, &[(gas, -1), (neutral, 1)]);
    IdealGasLiquidPartition::from_sources(
        equation.clone(),
        id(gas),
        id(neutral),
        &fixture.sources(&equation),
    )
    .unwrap()
}

fn volume(ratio: f64) -> Volume {
    // For K=1 and n_ref=1 mol, n_g/n_neutral = V p0/(R T).
    Volume::new(ratio * RT / 100_000.0)
}

#[test]
fn g1_finite_volatile_acid_and_base_match_independent_manufactured_balances() {
    let mut fixture = Fixture::new(1e-14, 1e-5, [1e-6, 1e-10]);
    let acid_partition = partition(&mut fixture, AG, HA, 0.0);
    let acid =
        VolatileAcidBaseSystem::new(fixture.system(vec![fixture.mono()]), vec![acid_partition])
            .unwrap();
    // At h=1e-4, A-=h-Kw/h, HA=h*A-/Ka, gas=2*HA.
    let anion = 1e-4 - 1e-10;
    let neutral = 10.0 * anion;
    let total_amount = anion + 3.0 * neutral;
    let result = acid
        .speciate(
            &[total("monoprotic", total_amount)],
            &[],
            Amount::new(1.0),
            volume(2.0),
            options(),
            || false,
        )
        .unwrap();
    near(result.aqueous.ph.value(), 4.0, 0.0, 2e-12);
    let phase = &result.volatile_families[0];
    near(phase.gas_amount.value(), 2.0 * neutral, 5e-12, 0.0);
    near(phase.aqueous_amount.value(), anion + neutral, 5e-12, 0.0);
    near(
        phase.partial_pressure.value(),
        neutral * 100_000.0,
        5e-12,
        0.0,
    );
    assert!(phase.component_balance_residual.value().abs() < 1e-17);
    assert!(phase.log_equilibrium_residual.unwrap().value().abs() < 1e-12);

    let base_family = add_base(&mut fixture, 1e-9);
    let base_partition = partition(&mut fixture, BG, B, 0.0);
    let base = VolatileAcidBaseSystem::new(fixture.system(vec![base_family]), vec![base_partition])
        .unwrap();
    // At h=1e-10, BH+=Kw/h-h, B=Ka*BH+/h, gas=2*B.
    let result = base
        .speciate(
            &[total("base", total_amount)],
            &[],
            Amount::new(1.0),
            volume(2.0),
            options(),
            || false,
        )
        .unwrap();
    near(result.aqueous.ph.value(), 10.0, 0.0, 2e-12);
    near(
        result.volatile_families[0].gas_amount.value(),
        2.0 * neutral,
        5e-12,
        0.0,
    );
    near(
        result.aqueous.families[0].forms[0].amount.value(),
        anion,
        5e-12,
        0.0,
    );
    assert_eq!(result.model().partitions().len(), 1);
    assert_eq!(result.model().reaction_network().reactions().len(), 3);
}

#[test]
fn g0_declared_base_enhances_acid_uptake_and_depletes_the_same_finite_gas() {
    let mut fixture = Fixture::new(1e-14, 1e-5, [1e-6, 1e-10]);
    let phase = partition(&mut fixture, AG, HA, 0.0);
    let model =
        VolatileAcidBaseSystem::new(fixture.system(vec![fixture.mono()]), vec![phase]).unwrap();
    let run = |spectators: &[SpectatorIon]| {
        model
            .speciate(
                &[total("monoprotic", 0.02)],
                spectators,
                Amount::new(1.0),
                volume(2.0),
                options(),
                || false,
            )
            .unwrap()
    };
    let initial = run(&[]);
    let with_base = run(&[ion("explicit-Na+", 1, 0.015)]);
    assert!(with_base.aqueous.ph > initial.aqueous.ph);
    assert!(
        with_base.volatile_families[0].aqueous_amount > initial.volatile_families[0].aqueous_amount
    );
    assert!(with_base.volatile_families[0].gas_amount < initial.volatile_families[0].gas_amount);
    let p = &with_base.volatile_families[0];
    near(
        p.gas_amount.value() + p.aqueous_amount.value(),
        0.02,
        5e-15,
        0.0,
    );
    assert!(with_base.aqueous.charge_balance_residual.value().abs() < 1e-13);
    assert!(with_base.aqueous.max_log_mass_action_residual.value() < 1e-12);
    assert!(p.log_equilibrium_residual.unwrap().value().abs() < 1e-12);
}

#[test]
fn g1_g3_multiple_gases_share_ph_nonvolatile_buffer_and_extensive_scaling() {
    let mut fixture = Fixture::new(1e-14, 1e-5, [1e-6, 1e-10]);
    let base_family = add_base(&mut fixture, 1e-9);
    let acid_partition = partition(&mut fixture, AG, HA, 0.0);
    let base_partition = partition(&mut fixture, BG, B, 0.0);
    let aqueous = fixture.system(vec![fixture.mono(), fixture.poly(), base_family]);
    let model = VolatileAcidBaseSystem::new(
        aqueous.clone(),
        vec![acid_partition.clone(), base_partition.clone()],
    )
    .unwrap();
    let reversed =
        VolatileAcidBaseSystem::new(aqueous, vec![base_partition, acid_partition]).unwrap();
    let h = 1e-8;
    // Independent explicit weights at the manufactured shared pH 8.
    let acid_neutral = 0.05 / (3.0 + 1e-5 / h);
    let acid_anion = acid_neutral * 1e-5 / h;
    let base_neutral = 0.03 / (3.0 + h / 1e-9);
    let base_cation = base_neutral * h / 1e-9;
    let poly_neutral = 0.02 / (1.0 + 1e-6 / h + 1e-16 / (h * h));
    let spectator =
        1e-14 / h - h + acid_anion + poly_neutral * (1e-6 / h + 2e-16 / (h * h)) - base_cation;
    let run = |model: &VolatileAcidBaseSystem, scale: f64| {
        model
            .speciate(
                &[
                    total("monoprotic", 0.05 * scale),
                    total("base", 0.03 * scale),
                    total("diprotic", 0.02 * scale),
                ],
                &[ion("explicit-Na+", 1, spectator * scale)],
                Amount::new(scale),
                Volume::new(volume(2.0).value() * scale),
                options(),
                || false,
            )
            .unwrap()
    };
    let answer = run(&model, 1.0);
    let canonical = run(&reversed, 1.0);
    let scaled = run(&model, 37.0);
    near(answer.aqueous.ph.value(), 8.0, 0.0, 2e-12);
    assert_eq!(answer.aqueous.ph, canonical.aqueous.ph);
    near(
        scaled.aqueous.ph.value(),
        answer.aqueous.ph.value(),
        0.0,
        2e-12,
    );
    assert_eq!(answer.aqueous.families[1].id, "diprotic");
    assert_eq!(answer.aqueous.families[1].total.value(), 0.02);
    for (index, neutral) in [base_neutral, acid_neutral].into_iter().enumerate() {
        let a = &answer.volatile_families[index];
        near(a.gas_amount.value(), 2.0 * neutral, 1e-11, 0.0);
        assert_eq!(a.gas_amount, canonical.volatile_families[index].gas_amount);
        near(
            scaled.volatile_families[index].gas_amount.value(),
            a.gas_amount.value() * 37.0,
            1e-11,
            0.0,
        );
        near(
            scaled.volatile_families[index].partial_pressure.value(),
            a.partial_pressure.value(),
            1e-11,
            0.0,
        );
        assert!(a.log_equilibrium_residual.unwrap().value().abs() < 1e-12);
    }
}

#[test]
fn g0_trace_gas_and_aqueous_inventories_are_not_subtraction_remainders() {
    for log_k in [-500.0, 500.0] {
        let mut fixture = Fixture::new(1e-14, 1e-5, [1e-6, 1e-10]);
        let phase = partition(&mut fixture, AG, HA, log_k);
        let model =
            VolatileAcidBaseSystem::new(fixture.system(vec![fixture.mono()]), vec![phase]).unwrap();
        let result = model
            .speciate(
                &[total("monoprotic", 0.01)],
                &[],
                Amount::new(1.0),
                volume(2.0),
                options(),
                || false,
            )
            .unwrap();
        let phase = &result.volatile_families[0];
        assert!(phase.gas_amount.value() > 0.0 && phase.aqueous_amount.value() > 0.0);
        assert!(phase.gas_amount.value().min(phase.aqueous_amount.value()) < 1e-200);
        near(
            phase.gas_amount.value() + phase.aqueous_amount.value(),
            0.01,
            1e-14,
            0.0,
        );
        assert!(phase.log_equilibrium_residual.unwrap().value().abs() < 2e-12);
    }
}

#[test]
fn g0_g4_empty_pools_and_invalid_or_cancelled_calls_publish_no_partial_result() {
    let mut fixture = Fixture::new(1e-14, 1e-5, [1e-6, 1e-10]);
    let phase = partition(&mut fixture, AG, HA, 0.0);
    let aqueous = fixture.system(vec![fixture.mono()]);
    assert!(VolatileAcidBaseSystem::new(aqueous.clone(), vec![]).is_err());
    assert!(
        VolatileAcidBaseSystem::new(aqueous.clone(), vec![phase.clone(), phase.clone()]).is_err()
    );
    let model = VolatileAcidBaseSystem::new(aqueous, vec![phase]).unwrap();
    let empty = model
        .speciate(
            &[total("monoprotic", 0.0)],
            &[],
            Amount::new(1.0),
            volume(2.0),
            options(),
            || false,
        )
        .unwrap();
    near(empty.aqueous.ph.value(), 7.0, 0.0, 2e-12);
    assert_eq!(empty.volatile_families[0].gas_amount.value(), 0.0);
    assert_eq!(empty.volatile_families[0].aqueous_amount.value(), 0.0);
    assert!(
        empty.volatile_families[0]
            .log_equilibrium_residual
            .is_none()
    );
    assert!(matches!(
        model.speciate(
            &[total("monoprotic", 0.01)],
            &[],
            Amount::new(1.0),
            volume(2.0),
            options(),
            || true,
        ),
        Err(AcidBaseError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        model.speciate(
            &[total("monoprotic", 0.01)],
            &[],
            Amount::new(1.0),
            volume(2.0),
            options(),
            || {
                polls += 1;
                polls > 12
            },
        ),
        Err(AcidBaseError::Cancelled)
    ));
    let mut short = options();
    short.max_iterations = 1;
    assert!(matches!(
        model.speciate(
            &[total("monoprotic", 0.01)],
            &[],
            Amount::new(1.0),
            volume(2.0),
            short,
            || false,
        ),
        Err(AcidBaseError::Budget)
    ));
    for invalid_volume in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .speciate(
                    &[total("monoprotic", 0.01)],
                    &[],
                    Amount::new(1.0),
                    Volume::new(invalid_volume),
                    options(),
                    || false,
                )
                .is_err()
        );
    }
    assert!(
        model
            .speciate(
                &[total("monoprotic", 0.01)],
                &[ion(AG, 1, 0.01)],
                Amount::new(1.0),
                volume(2.0),
                options(),
                || false,
            )
            .is_err()
    );
    assert!(
        model
            .speciate(&[], &[], Amount::new(1.0), volume(2.0), options(), || false,)
            .is_err()
    );
}

#[test]
fn g0_mismatched_shared_neutral_source_and_unbound_family_are_refused() {
    let mut fixture = Fixture::new(1e-14, 1e-5, [1e-6, 1e-10]);
    let aqueous = fixture.system(vec![fixture.mono()]);
    let index = fixture.index(HA);
    fixture.packs[index] = pack(HA, 123.0, &fixture.points[index], None);
    let phase = partition(&mut fixture, AG, HA, 0.0);
    assert!(VolatileAcidBaseSystem::new(aqueous, vec![phase.clone()]).is_err());
    assert!(
        VolatileAcidBaseSystem::new(fixture.system(vec![fixture.poly()]), vec![phase]).is_err()
    );
}
