//! Independent analytical ideal-model limits, not measured recovery fixtures.
use super::*;
use crate::stoichiometry::thermodynamics::equilibrium::{
    EquilibriumExtentStatus as Status, EquilibriumInventory, EquilibriumSolveOptions,
    IdealActivityModel as Law, SingleReactionEquilibrium,
};
use fs_qty::Amount;

fn options() -> EquilibriumSolveOptions {
    EquilibriumSolveOptions {
        extent_tolerance: Amount::new(1e-12),
        log_residual_tolerance: Dimensionless::new(1e-11),
        max_iterations: 256,
    }
}

struct Fixture {
    reaction: ConservedReaction,
    points: [QueryPoint; 2],
    packs: [NormalizedPack; 2],
    inventory: Vec<EquilibriumInventory>,
}
impl Fixture {
    fn new(scale: i128, aqueous: [bool; 2], amounts: [f64; 2], gibbs: f64) -> Self {
        let reaction = exchange(scale, 0);
        let points = aqueous.map(|aq| {
            point(
                if aq {
                    "source-phase-aqueous"
                } else {
                    "source-phase-solid"
                },
                false,
            )
        });
        let packs = [
            pack(0.0, &points[0], None, PROPERTY),
            pack(gibbs, &points[1], None, PROPERTY),
        ];
        let inventory = reaction
            .terms()
            .iter()
            .enumerate()
            .map(|(i, term)| EquilibriumInventory {
                species: term.species.clone(),
                amount: Amount::new(amounts[i]),
                activity: if aqueous[i] {
                    Law::IdealDilute {
                        amount_at_unit_activity: Amount::new(1.0),
                    }
                } else {
                    Law::PureCondensed
                },
            })
            .collect();
        Self {
            reaction,
            points,
            packs,
            inventory,
        }
    }

    fn solve(&self) -> Result<SingleReactionEquilibrium> {
        self.reaction.equilibrate_ideal(
            &sources(&self.reaction, &self.packs, &self.points),
            &self.inventory,
            options(),
        )
    }
}

fn approximate(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn ideal_solution_matches_analytical_equilibrium_in_both_directions() {
    // Independent K=10 isomer balance: xi=(K*n_A0-n_B0)/(1+K).
    // Literal standard energy is -R*T*ln(10), not fitted from the solver.
    let gibbs = -5708.009_502_535_259;
    let forward = Fixture::new(1, [true, true], [4.0, 1.0], gibbs)
        .solve()
        .unwrap();
    let reversed = Fixture::new(-1, [true, true], [4.0, 1.0], gibbs)
        .solve()
        .unwrap();
    let tripled = Fixture::new(3, [true, true], [4.0, 1.0], gibbs)
        .solve()
        .unwrap();
    for result in [&forward, &reversed, &tripled] {
        assert_eq!(result.status, Status::Interior);
        approximate(result.final_inventory[0].1.value(), 5.0 / 11.0);
        approximate(result.final_inventory[1].1.value(), 50.0 / 11.0);
        assert!(result.log_quotient_over_equilibrium.unwrap().value().abs() < 1e-11);
        assert!(result.extent_bracket[1].value() - result.extent_bracket[0].value() <= 1e-12);
        assert_eq!(result.standard.terms().len(), 2);
    }
    approximate(forward.extent.value(), 39.0 / 11.0);
    approximate(reversed.extent.value(), -39.0 / 11.0);
    approximate(tripled.extent.value(), 13.0 / 11.0);
}

#[test]
fn mineral_dissolution_reaches_saturation_or_inventory_limit() {
    // Authored S(s) <=> S(aq), K=1 and n_ref=1: saturation is exactly 1 mol.
    let saturated = Fixture::new(1, [false, true], [2.0, 0.0], 0.0)
        .solve()
        .unwrap();
    assert_eq!(saturated.status, Status::Interior);
    assert_eq!(saturated.extent.value(), 1.0);
    assert_eq!(saturated.final_inventory[0].1.value(), 1.0);
    assert_eq!(saturated.final_inventory[1].1.value(), 1.0);

    let depleted = Fixture::new(1, [false, true], [0.25, 0.0], 0.0)
        .solve()
        .unwrap();
    assert_eq!(depleted.status, Status::ForwardInventoryLimit);
    assert_eq!(depleted.extent.value(), 0.25);
    assert_eq!(depleted.final_inventory[0].1.value(), 0.0);
    approximate(
        depleted.log_quotient_over_equilibrium.unwrap().value(),
        -1.386_294_361_119_890_6,
    );
}

#[test]
fn absent_mineral_precipitates_only_when_supersaturated() {
    let supersaturated = Fixture::new(1, [false, true], [0.0, 2.0], 0.0)
        .solve()
        .unwrap();
    assert_eq!(supersaturated.status, Status::Interior);
    assert_eq!(supersaturated.extent.value(), -1.0);
    assert_eq!(supersaturated.final_inventory[0].1.value(), 1.0);
    assert_eq!(supersaturated.final_inventory[1].1.value(), 1.0);

    let undersaturated = Fixture::new(1, [false, true], [0.0, 0.25], 0.0)
        .solve()
        .unwrap();
    assert_eq!(undersaturated.status, Status::ForwardInventoryLimit);
    assert_eq!(undersaturated.extent.value(), 0.0);
    assert_eq!(undersaturated.final_inventory[0].1.value(), 0.0);
    assert_eq!(undersaturated.final_inventory[1].1.value(), 0.25);
    assert!(
        undersaturated
            .log_quotient_over_equilibrium
            .unwrap()
            .value()
            < 0.0
    );
}

#[test]
fn pure_phase_equality_is_nonunique_and_inequalities_choose_the_correct_end() {
    let indifferent = Fixture::new(1, [false, false], [2.0, 3.0], 0.0)
        .solve()
        .unwrap();
    assert_eq!(indifferent.status, Status::NonUnique);
    assert_eq!(indifferent.extent.value(), 0.0);
    assert_eq!(indifferent.extent_bracket.map(|v| v.value()), [-3.0, 2.0]);
    let forward = Fixture::new(1, [false, false], [2.0, 3.0], -100.0)
        .solve()
        .unwrap();
    assert_eq!(forward.status, Status::ForwardInventoryLimit);
    assert_eq!(forward.extent.value(), 2.0);
    let reverse = Fixture::new(1, [false, false], [2.0, 3.0], 100.0)
        .solve()
        .unwrap();
    assert_eq!(reverse.status, Status::ReverseInventoryLimit);
    assert_eq!(reverse.extent.value(), -3.0);
}

#[test]
fn zero_stocks_are_blocked_without_log_zero_or_invented_material() {
    let result = Fixture::new(1, [true, true], [0.0, 0.0], -100.0)
        .solve()
        .unwrap();
    assert_eq!(result.status, Status::InventoryBlocked);
    assert_eq!(result.extent.value(), 0.0);
    assert_eq!(result.log_quotient_over_equilibrium, None);
    assert!(result.final_inventory.iter().all(|(_, n)| n.value() == 0.0));
}

#[test]
fn fractional_limit_retains_fma_remainder_without_overspending() {
    let result = Fixture::new(3, [false, true], [1.0, 0.0], -10000.0)
        .solve()
        .unwrap();
    assert_eq!(result.status, Status::ForwardInventoryLimit);
    assert!(result.extent.value().mul_add(3.0, -1.0) <= 0.0);
    assert_eq!(
        result.final_inventory[0].1.value(),
        (-3.0_f64).mul_add(result.extent.value(), 1.0)
    );
    assert!(result.final_inventory[0].1.value() > 0.0);
    approximate(
        result.final_inventory.iter().map(|(_, n)| n.value()).sum(),
        1.0,
    );
}

#[test]
fn finite_logs_survive_unrepresentable_activity_ratios_and_wide_brackets() {
    let mut fixture = Fixture::new(1, [true, true], [1e300, 1e300], 0.0);
    for entry in &mut fixture.inventory {
        entry.activity = Law::IdealDilute {
            amount_at_unit_activity: Amount::new(1e-300),
        };
    }
    let result = fixture.solve().unwrap();
    assert_eq!(result.status, Status::Interior);
    assert_eq!(result.extent.value(), 0.0);

    let mut fixture = Fixture::new(1, [true, true], [1e-250, 0.0], 0.0);
    for entry in &mut fixture.inventory {
        entry.activity = Law::IdealDilute {
            amount_at_unit_activity: Amount::new(1e250),
        };
    }
    let result = fixture
        .reaction
        .equilibrate_ideal(
            &sources(&fixture.reaction, &fixture.packs, &fixture.points),
            &fixture.inventory,
            EquilibriumSolveOptions {
                extent_tolerance: Amount::new(1e-265),
                ..options()
            },
        )
        .unwrap();
    assert_eq!(result.status, Status::Interior);
    assert_eq!(result.extent.value(), 5e-251);
}

#[test]
fn complete_inputs_are_canonical_and_original_source_receipts_remain_valid() {
    let mut fixture = Fixture::new(1, [true, true], [4.0, 1.0], -5708.009_502_535_259);
    let original = fixture.solve().unwrap();
    fixture.inventory.reverse();
    let mut inputs = sources(&fixture.reaction, &fixture.packs, &fixture.points);
    inputs.reverse();
    let permuted = fixture
        .reaction
        .equilibrate_ideal(&inputs, &fixture.inventory, options())
        .unwrap();
    assert_eq!(original.extent, permuted.extent);
    assert_eq!(original.final_inventory, permuted.final_inventory);
    assert_eq!(original.initial, permuted.initial);
    for (index, term) in original.standard.terms().iter().enumerate() {
        fixture.packs[index]
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
}

#[test]
fn incomplete_duplicate_nonfinite_and_mismatched_phase_models_refuse() {
    let fixture = Fixture::new(1, [true, true], [1.0, 0.0], 0.0);
    let src = sources(&fixture.reaction, &fixture.packs, &fixture.points);
    assert!(
        fixture
            .reaction
            .equilibrate_ideal(&src, &fixture.inventory[..1], options())
            .is_err()
    );
    let duplicate = [fixture.inventory[0].clone(), fixture.inventory[0].clone()];
    assert!(
        fixture
            .reaction
            .equilibrate_ideal(&src, &duplicate, options())
            .is_err()
    );
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        let mut inventory = fixture.inventory.clone();
        inventory[0].amount = Amount::new(bad);
        assert!(
            fixture
                .reaction
                .equilibrate_ideal(&src, &inventory, options())
                .is_err()
        );
    }
    let mut inventory = fixture.inventory.clone();
    inventory[0].activity = Law::PureCondensed;
    assert!(
        fixture
            .reaction
            .equilibrate_ideal(&src, &inventory, options())
            .is_err()
    );
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        inventory[0].activity = Law::IdealDilute {
            amount_at_unit_activity: Amount::new(scale),
        };
        assert!(
            fixture
                .reaction
                .equilibrate_ideal(&src, &inventory, options())
                .is_err()
        );
    }
    let mut solid = Fixture::new(1, [false, true], [1.0, 0.0], 0.0);
    solid.inventory[0].activity = Law::IdealDilute {
        amount_at_unit_activity: Amount::new(1.0),
    };
    assert!(solid.solve().is_err());
}

#[test]
fn bounded_work_and_unresolved_interior_equilibrium_refuse_without_changing_inputs() {
    let fixture = Fixture::new(1, [true, true], [4.0, 1.0], -5708.009_502_535_259);
    let before = fixture.inventory.clone();
    let src = sources(&fixture.reaction, &fixture.packs, &fixture.points);
    assert!(
        fixture
            .reaction
            .equilibrate_ideal(
                &src,
                &fixture.inventory,
                EquilibriumSolveOptions {
                    max_iterations: 1,
                    ..options()
                }
            )
            .is_err()
    );
    assert_eq!(before, fixture.inventory);
    for bad in [0, 4097] {
        assert!(
            fixture
                .reaction
                .equilibrate_ideal(
                    &src,
                    &fixture.inventory,
                    EquilibriumSolveOptions {
                        max_iterations: bad,
                        ..options()
                    }
                )
                .is_err()
        );
    }
    // Equilibrium demands a residual reactant below the representable extent
    // spacing near xi=1. It must not be relabelled complete conversion.
    let impossible = Fixture::new(1, [true, true], [1.0, 0.0], -250_000.0);
    assert!(matches!(
        impossible.solve(),
        Err(ReactionBalanceError::Arithmetic(_))
    ));
}
