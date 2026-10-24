//! G1/G3 material-driven heat allocation through the real coupled thermal solve.
use super::*;
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::Density;

fn source(cp: f64, knots: Vec<(f64, f64)>) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = common::material(1.0, cp, 0.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: knots.clone(),
                dims: THERMAL_CONDUCTIVITY_DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &common::point(200.0),
            upper: &common::point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, caloric)
}

fn material_binding<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    interface: &'a InterfaceSystemCard,
    interface_point: &'a QueryPoint,
    bulk_point: &'a QueryPoint,
    depths: [f64; 2],
) -> SlidingInterface<'a> {
    let mut data = binding(
        cards,
        curves,
        interface,
        interface_point,
        HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
    );
    data.partition = SlidingHeatPartition::SemiInfiniteEffusivity(SemiInfiniteHeatPartition {
        surfaces: [0, 1].map(|i| ThermalMaterialPoint {
            material: &cards[i],
            context: bulk_point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        }),
        selection: [
            MaterialPropertySelection::SingleClaimOnly,
            MaterialPropertySelection::SingleClaimOnly,
        ],
        thermal_depth_m: depths,
        maximum_penetration_fraction: 0.5,
    });
    data.outside_destination = None;
    data
}

fn constant_interface(cards: &[MaterialCard; 2]) -> InterfaceSystemCard {
    interface(
        cards,
        "Ta",
        Temperature::DIMS,
        [(200.0, 0.2), (2000.0, 0.2)],
        2000.0,
    )
}
fn drive(distance: f64) -> SlidingDrive {
    SlidingDrive {
        normal_load_n: 100.0,
        relative_increment_m: [distance, 0.0],
    }
}

#[test]
fn conductivity_and_caloric_swaps_change_actual_heat_stores_without_retuning() {
    let point = context();
    let bulk = common::point(777.0); // stale authoring temperature must be replaced
    for cp_a in [4.0_f64, 16.0] {
        for k_a in [1e-6_f64, 9e-6] {
            let (a, ca) = common::material_with_conductivity(1.0, cp_a, 0.0, k_a);
            let (b, cb) = common::material_with_conductivity(1.0, 4.0, 0.0, 1e-6);
            let cards = [a, b];
            let curves = [ca, cb];
            let thermal = common::initial(&curves, 0.0, None);
            let card = constant_interface(&cards);
            let contact = MaterialSlidingContact::new(
                &thermal,
                material_binding(&cards, &curves, &card, &point, &bulk, [0.01; 2]),
                100.0,
                1.0,
                InputAuthority::SyntheticFixture,
            )
            .unwrap();
            let candidate = contact
                .step_candidate(
                    &thermal,
                    0.2,
                    drive(0.2),
                    controls(),
                    &[0.0; 2],
                    other,
                    || false,
                )
                .unwrap();
            let report = candidate.report();
            // Independent half-space impedance ratio, not implementation output.
            let fraction = (k_a * cp_a).sqrt() / ((k_a * cp_a).sqrt() + (1e-6_f64 * 4.0).sqrt());
            close(report.drive_work_j, 4.0);
            close(report.partition.surface_a_heat_j(), 4.0 * fraction);
            close(report.partition.surface_b_heat_j(), 4.0 * (1.0 - fraction));
            close(
                candidate.thermal().temperature_k(CELLS[0]).unwrap(),
                300.0 + 4.0 * fraction / cp_a,
            );
            close(
                candidate.thermal().temperature_k(CELLS[1]).unwrap(),
                300.0 + 1.0 - fraction,
            );
            close(report.total_energy_residual_j, 0.0);
            let source = report.material.effusivity.as_ref().unwrap();
            for side in 0..2 {
                close(
                    source.conductivity[side].query().axes()["T"],
                    candidate.thermal().temperature_k(CELLS[side]).unwrap(),
                );
                assert_eq!(source.conductivity[side].query().axes()["lot"], 7.0);
            }
            close(source.specific_heat_capacity_j_kg_k[0], cp_a);
            close(source.contact_age_s, 0.2);
            assert_eq!(contact.work().dissipated_work_j(), 0.0);
        }
    }
}

#[test]
fn endpoint_conductivity_feedback_uses_the_same_temperature_as_heat_storage() {
    let (a, ca) = source(
        4.0,
        vec![
            (200.0, 1e-6),
            (300.0, 1e-6),
            (302.0, 21e-6),
            (2000.0, 21e-6),
        ],
    );
    let (b, cb) = common::material_with_conductivity(1.0, 4.0, 0.0, 1e-6);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let bulk = common::point(700.0);
    let thermal = common::initial(&curves, 0.0, None);
    let card = constant_interface(&cards);
    let contact = MaterialSlidingContact::new(
        &thermal,
        material_binding(&cards, &curves, &card, &point, &bulk, [0.01; 2]),
        100.0,
        1.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let candidate = contact
        .step_candidate(
            &thermal,
            0.2,
            drive(0.2),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    // T_a=300+beta, k_a/k_b=1+10 beta. Solve the independent scalar balance.
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    for _ in 0..80 {
        let beta = 0.5 * (low + high);
        let e = (1.0 + 10.0 * beta).sqrt();
        if beta < e / (e + 1.0) {
            low = beta;
        } else {
            high = beta;
        }
    }
    let fraction = 0.5 * (low + high);
    assert!(fraction > 0.7);
    close(candidate.report().material.partition.surface_a(), fraction);
    close(
        candidate.thermal().temperature_k(CELLS[0]).unwrap(),
        300.0 + fraction,
    );
    close(
        candidate.thermal().temperature_k(CELLS[1]).unwrap(),
        301.0 - fraction,
    );
    close(candidate.report().total_energy_residual_j, 0.0);
    let SlidingHeatPartition::SemiInfiniteEffusivity(pinned) =
        &candidate.contact().binding().partition
    else {
        panic!("selected effusivity model")
    };
    assert!(matches!(
        pinned.selection[0],
        MaterialPropertySelection::PinnedByProperty(_)
    ));
}

#[test]
fn diffusion_envelope_accumulates_across_steps_and_cloned_restarts() {
    let (a, ca) = common::material_with_conductivity(1.0, 4.0, 0.0, 9e-6);
    let (b, cb) = common::material_with_conductivity(1.0, 4.0, 0.0, 1e-6);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let bulk = common::point(300.0);
    let thermal = common::initial(&curves, 0.0, None);
    let card = constant_interface(&cards);
    // The A-side bound is reached at 0.0625 s. Each individual 0.04 s step
    // fits, but the original contact's accumulated 0.08 s exposure does not.
    let data = material_binding(&cards, &curves, &card, &point, &bulk, [0.0015, 0.01]);
    let contact = MaterialSlidingContact::new(
        &thermal,
        data.clone(),
        100.0,
        1.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let first = contact
        .step_candidate(
            &thermal,
            0.04,
            drive(0.04),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let restarted = first.contact().clone();
    let old_cells = first.thermal().cells().clone();
    let old_work = restarted.work().dissipated_work_j();
    let rejected = restarted.step_candidate(
        first.thermal(),
        0.08,
        drive(0.04),
        controls(),
        &[0.0; 2],
        other,
        || false,
    );
    assert!(
        matches!(rejected,Err(SlidingError::Invalid(message)) if message.contains("diffusion-depth"))
    );
    assert_eq!(first.thermal().cells(), &old_cells);
    assert_eq!(restarted.work().dissipated_work_j(), old_work);
    assert_eq!(restarted.time_s(), 0.04);
    let mut foreign = data;
    let SlidingHeatPartition::SemiInfiniteEffusivity(model) = &mut foreign.partition else {
        panic!("selected effusivity model")
    };
    model.surfaces.swap(0, 1);
    assert!(matches!(
        MaterialSlidingContact::new(
            &thermal,
            foreign,
            100.0,
            1.0,
            InputAuthority::SyntheticFixture
        ),
        Err(SlidingError::Invalid(_))
    ));
}

#[test]
fn a_conductivity_source_exit_refuses_the_entire_thermal_and_sliding_candidate() {
    let (a, ca) = source(4.0, vec![(200.0, 1e-6), (300.01, 1e-6)]);
    let (b, cb) = common::material_with_conductivity(1.0, 4.0, 0.0, 1e-6);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let bulk = common::point(300.0);
    let thermal = common::initial(&curves, 0.0, None);
    let card = constant_interface(&cards);
    let contact = MaterialSlidingContact::new(
        &thermal,
        material_binding(&cards, &curves, &card, &point, &bulk, [0.01; 2]),
        100.0,
        1.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let original = thermal.cells().clone();
    assert!(matches!(
        contact.step_candidate(
            &thermal,
            0.2,
            drive(0.2),
            controls(),
            &[0.0; 2],
            other,
            || false
        ),
        Err(SlidingError::Material(_))
    ));
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(contact.state().slip_m(), [0.0; 2]);
    assert_eq!(contact.work().dissipated_work_j(), 0.0);
}
