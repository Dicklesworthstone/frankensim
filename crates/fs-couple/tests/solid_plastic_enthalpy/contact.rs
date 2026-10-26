//! Closed pressure/thermal/constitutive feedback on an actual native tetrahedron.
#[path = "contact/dynamics.rs"]
mod dynamics;
use super::*;
use fs_conduction::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RDIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RNAME,
};
use fs_couple::solid_enthalpy::plastic::contact::{
    ContactError, ContactSolve, ThermalSupportContact,
};
use fs_matdb::{ClaimId, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use std::convert::Infallible;

fn controls() -> ContactSolve {
    ContactSolve {
        maximum_cells: 2,
        maximum_mechanical_trials: 400,
        maximum_newton_iterations: 25,
        difference_step_j: 1e-7,
        cell_energy_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-10,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn other(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}
fn point() -> QueryPoint {
    // Deliberately stale pressure: the adapter must replace it with solved stress.
    QueryPoint::new()
        .with("Ta", 300.0)
        .unwrap()
        .with("Tb", 300.0)
        .unwrap()
        .with("pressure", 7e6)
        .unwrap()
}
fn interface(
    a: &Fixture,
    b: &Fixture,
    scale: f64,
    maximum_pressure: f64,
) -> (InterfaceSystemCard, ClaimId) {
    let mut claim = scalar(RNAME, RDIMS, 0.002 * scale);
    claim.value = PropertyValue::Curve {
        abscissa: "pressure".into(),
        abscissa_dims: Pressure::DIMS,
        knots: vec![(0.0, 0.002 * scale), (maximum_pressure, 0.0004 * scale)],
        dims: RDIMS,
    };
    claim.validity = ValidityDomain::unconstrained()
        .with("pressure", 0.0, maximum_pressure)
        .with("Ta", 280.0, 320.0)
        .with("Tb", 280.0, 320.0);
    claim.interpolation = InterpolationPolicy::LinearInside;
    let mut claims = ClaimSet::new();
    let pin = claims.insert_claim(claim).unwrap();
    (
        InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: a.card.id().clone(),
                texture_frame: "specimen finish".into(),
            },
            SurfaceSpec {
                material: b.card.id().clone(),
                texture_frame: "support finish".into(),
            },
            SystemContext {
                medium: "synthetic dry contact".into(),
                third_body: None,
                environment: "numerical fixture".into(),
                history: "fixed reference footprint".into(),
            },
            claims,
            vec![],
        )
        .unwrap(),
        pin,
    )
}
fn binding<'a>(
    a: &'a Fixture,
    b: &'a Fixture,
    card: &'a InterfaceSystemCard,
    claim: ClaimId,
    context: &'a QueryPoint,
) -> ThermalSupportContact<'a> {
    ThermalSupportContact {
        link: 0,
        face_nodes: [2, 0, 1],
        point: ThermalInterfacePoint {
            interface: card,
            materials: [&a.card, &b.card],
            context,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        caloric: [a.caloric.curve(), b.caloric.curve()],
        claim,
        pressure_axis: "pressure",
    }
}
fn prescribed(compression: f64, shear: f64) -> Vec<(usize, f64)> {
    (0..12)
        .map(|d| {
            (
                d,
                -compression * NODES[d / 3][d % 3]
                    + if d % 3 == 0 {
                        2.0 * shear * NODES[d / 3][1]
                    } else {
                        0.0
                    },
            )
        })
        .collect()
}

#[test]
fn solved_pressure_and_interface_data_recover_the_two_finite_body_oracle() {
    context(|cx| {
        let a = Fixture::new(500.0);
        let b = Fixture::new(1000.0);
        let thermal = initial(&[&a, &b], &[(0, 1, 0.0)], &[0.2, 0.0]);
        let solid = bind(&thermal, &[&a], &NODES, &TETS, cx);
        let original = solid.solution().clone();
        let old_cells = thermal.cells().clone();
        let context = point();
        let mut temperatures = Vec::new();
        for (compression, scale) in [(0.001, 1.0), (0.002, 1.0), (0.001, 3.0)] {
            let (card, pin) = interface(&a, &b, scale, 8e6);
            let contacts = [binding(&a, &b, &card, pin, &context)];
            let supports = prescribed(compression, 0.0);
            let next = solid
                .advance_with_thermal_contacts(
                    &thermal,
                    &contacts,
                    1.0,
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    other,
                    cx,
                )
                .unwrap();
            let report = next.report();
            let sample = &report.contacts[0];
            // sigma=-3K*compression*I; finite Cp stores use the exact BE solution.
            let pressure = 2e9 * compression;
            let r = scale * (0.002 - 0.0016 * pressure / 8e6);
            let g = 0.5 * L * L / r;
            let ca = RHO * V * 500.0;
            let cb = RHO * V * 1000.0;
            let ta = 300.0 + 0.2 / (ca + g * cb / (cb + g));
            let tb = 300.0 + g / (cb + g) * (ta - 300.0);
            close(sample.compression_pa, pressure, 1e-6);
            assert_eq!(sample.outward_normal, [0.0, 0.0, -1.0]);
            close(sample.resistance.value_m2_k_per_w(), r, 1e-15);
            assert_eq!(sample.resistance.card_identity(), card.content_hash());
            assert_eq!(sample.resistance.receipt().selected, pin);
            close(next.thermal().temperature_k(CellId(100)).unwrap(), ta, 1e-7);
            close(next.thermal().temperature_k(CellId(101)).unwrap(), tb, 1e-7);
            close(sample.heat_a_to_b_j, cb * (tb - 300.0), 1e-9);
            close(report.balance.external_heat_j, 0.2, 1e-10);
            close(
                report.support_work_j,
                3.0 * pressure * compression * V,
                1e-9,
            );
            close(report.balance.combined_residual_j, 0.0, 1e-8);
            close(report.balance.physical_heat_j, 0.0, 1e-12);
            close(report.balance.applied_heat_j, 0.0, 1e-12);
            assert!(report.balance.numerical_loss_j > 0.0);
            assert_eq!(next.solid().solution().element_states[0].alpha, 0.0);
            temperatures.push(tb);
        }
        assert!(temperatures[1] > temperatures[0]); // More compression, better contact.
        assert!(temperatures[2] < temperatures[0]); // Only interface source data changed.
        assert_eq!(solid.solution(), &original);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn material_expansion_changes_contact_cooling_and_closes_plastic_thermal_support_work() {
    context(|cx| {
        let mut pressures = Vec::new();
        for alpha in [1e-5, 4e-5] {
            let a = Fixture::with_expansion(500.0, Some(alpha));
            let b = Fixture::new(1000.0);
            let thermal = initial(&[&a, &b], &[(0, 1, 0.0)], &[0.2, 0.0]);
            let solid = bind(&thermal, &[&a], &NODES, &TETS, cx);
            let (card, pin) = interface(&a, &b, 1.0, 8e6);
            let point = point();
            let contacts = [binding(&a, &b, &card, pin, &point)];
            let supports = prescribed(0.001, 0.005);
            let next = solid
                .advance_with_thermal_contacts(
                    &thermal,
                    &contacts,
                    1.0,
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    other,
                    cx,
                )
                .unwrap();
            let r = next.report();
            let balance = &r.balance;
            let sample = &r.contacts[0];
            let ta = next.thermal().temperature_k(CellId(100)).unwrap();
            let tb = next.thermal().temperature_k(CellId(101)).unwrap();
            close(
                sample.compression_pa,
                2e9 * (0.001 + alpha * (ta - 300.0)),
                1e-5,
            );
            close(sample.phases[0].temperature_k(), ta, 1e-12);
            close(sample.phases[1].temperature_k(), tb, 1e-12);
            assert!(ta > tb && tb > 300.0);
            assert!(balance.physical_heat_j > 0.0);
            assert!(balance.reversible_heat_j > 0.0);
            assert!(balance.numerical_loss_j > 0.0);
            assert!(next.solid().solution().element_states[0].alpha > 0.0);
            close(
                balance.applied_heat_j,
                balance.physical_heat_j + balance.reversible_heat_j,
                1e-10,
            );
            close(
                next.thermal().energy_j() - thermal.energy_j()
                    + balance.stored_mechanical_change_j
                    + balance.thermoelastic_storage_change_j
                    + balance.numerical_loss_j
                    - r.support_work_j,
                0.2,
                1e-8,
            );
            close(sample.heat_a_to_b_j, RHO * V * 1000.0 * (tb - 300.0), 1e-9);
            pressures.push(sample.compression_pa);
            // Hold the accepted displacement: no old loading work can be applied again.
            let held = next
                .solid()
                .advance_with_thermal_contacts(
                    next.thermal(),
                    &contacts,
                    2.0,
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    other,
                    cx,
                )
                .unwrap();
            close(held.report().support_work_j, 0.0, 1e-14);
            close(held.report().balance.combined_residual_j, 0.0, 1e-8);
            assert_eq!(held.report().contacts[0].resistance.receipt().selected, pin);
            assert_eq!(solid.solution().element_states[0].alpha, 0.0);
            assert_eq!(next.solid().time_s(), 1.0);
        }
        assert!(pressures[1] > pressures[0]);
    });
}

#[test]
fn late_contact_source_refusal_tension_and_cancellation_preserve_both_owners() {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 0x435450,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let a = Fixture::new(500.0);
        let b = Fixture::new(1000.0);
        let thermal = initial(&[&a, &b], &[(0, 1, 0.0)], &[0.2, 0.0]);
        let solid = bind(&thermal, &[&a], &NODES, &TETS, &cx);
        let old_cells = thermal.cells().clone();
        let old_solid = solid.solution().clone();
        let point = point();
        for (maximum_pressure, compression) in [(0.5e6, 0.001), (8e6, -0.001)] {
            let (card, pin) = interface(&a, &b, 1.0, maximum_pressure);
            let contacts = [binding(&a, &b, &card, pin, &point)];
            let result = solid.advance_with_thermal_contacts(
                &thermal,
                &contacts,
                1.0,
                &prescribed(compression, 0.005),
                &[],
                IncrementSettings::default(),
                controls(),
                &[0.0; 2],
                other,
                &cx,
            );
            if compression > 0.0 {
                assert!(matches!(result, Err(ContactError::Conduction(_))));
            } else {
                assert!(matches!(result, Err(ContactError::Invalid(_))));
            }
        }
        let (card, pin) = interface(&a, &b, 1.0, 8e6);
        let copied = b.caloric.curve().clone();
        let mut wrong = binding(&a, &b, &card, pin, &point);
        wrong.caloric[1] = &copied;
        assert!(matches!(
            solid.advance_with_thermal_contacts(
                &thermal,
                &[wrong],
                1.0,
                &prescribed(0.001, 0.0),
                &[],
                IncrementSettings::default(),
                controls(),
                &[0.0; 2],
                other,
                &cx
            ),
            Err(ContactError::Invalid(_))
        ));
        let calls = std::cell::Cell::new(0);
        let contacts = [binding(&a, &b, &card, pin, &point)];
        let cancelled = solid.advance_with_thermal_contacts(
            &thermal,
            &contacts,
            1.0,
            &prescribed(0.001, 0.005),
            &[],
            IncrementSettings::default(),
            controls(),
            &[0.0; 2],
            |endpoint| {
                calls.set(calls.get() + 1);
                gate.request();
                other(endpoint)
            },
            &cx,
        );
        assert_eq!(calls.get(), 1); // Cancellation after a real mechanical proposal.
        assert!(matches!(cancelled, Err(ContactError::Cancelled)));
        assert_eq!(solid.solution(), &old_solid);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
