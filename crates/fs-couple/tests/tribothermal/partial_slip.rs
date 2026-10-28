//! G1/G3 analytical references use explicit synthetic stiffness and material
//! sources; these are implementation checks, not measured contact calibration.
use super::common;
use fs_couple::tribothermal::partial_slip::*;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::HeatPartition;
use fs_tribo::partial_slip::*;
use std::collections::BTreeMap;
use std::convert::Infallible;

const DESTINATIONS: [CellId; 2] = [CellId(30), CellId(10)];

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-7, "{a} != {b}");
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}
fn parameters(mu: f64) -> PartialSlipParameters {
    PartialSlipParameters {
        static_mu: mu,
        kinetic_mu: mu,
        tangential_stiffness_n_per_m: 1000.0,
        torsional_stiffness_nm_per_rad: 10.0,
        torsional_capacity_factor: 1.0,
        partial_slip_onset_fraction: 0.2,
        partial_slip_hardening_fraction: 0.25,
    }
}
fn model(mu: f64, slopes: [StiffnessDerivative; 2]) -> AffineContactStiffness {
    AffineContactStiffness::new(
        PartialSlipLaw::new(
            "declared affine contact",
            "synthetic free-energy fixture",
            parameters(mu),
        )
        .unwrap(),
        [300.0; 2],
        slopes,
        [[200.0, 600.0]; 2],
    )
    .unwrap()
}
fn slopes() -> [StiffnessDerivative; 2] {
    [
        StiffnessDerivative {
            tangential_n_per_m_k: -1.0,
            torsional_nm_per_rad_k: -0.01,
        },
        StiffnessDerivative {
            tangential_n_per_m_k: -0.5,
            torsional_nm_per_rad_k: -0.02,
        },
    ]
}
fn zero_slopes() -> [StiffnessDerivative; 2] {
    [StiffnessDerivative {
        tangential_n_per_m_k: 0.0,
        torsional_nm_per_rad_k: 0.0,
    }; 2]
}
fn initial<'a>(
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    temperatures: [f64; 2],
    conductance: f64,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for i in 0..2 {
        let knots = curves[i].curve().knots();
        let cp = (knots[1].specific_enthalpy_j_kg - knots[0].specific_enthalpy_j_kg)
            / (knots[1].temperature_k - knots[0].temperature_k);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: DESTINATIONS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: (temperatures[i] - 200.0) * cp,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: DESTINATIONS[0],
            b: DESTINATIONS[1],
            area_m2: 1.0,
            conductance_w_k: conductance,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn contact<'a>(
    thermal: &DepositionSimulation<'_>,
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    model: AffineContactStiffness,
    strain: [f64; 3],
) -> CaloricPartialSlip<'a> {
    let temperatures = DESTINATIONS.map(|id| thermal.temperature_k(id).unwrap());
    let checkpoint = PartialSlipCheckpoint::new(
        NormalPatchView::new(
            "patch",
            "declared held footprint",
            "synthetic patch",
            NormalPatchAuthority::SyntheticFixture,
            10.0,
            0.1,
            0.1,
            0.005,
        )
        .unwrap(),
        PartialSlipInterface::new(
            "ordered A/B",
            "loaded history",
            "synthetic pair",
            NormalPatchAuthority::SyntheticFixture,
        )
        .unwrap(),
        model.law_at::<Infallible>(temperatures).unwrap(),
        PartialSlipState::from_region_checkpoint(
            [[strain[0], strain[1]]; 2],
            [strain[2]; 2],
            [false; 2],
            Some(frame()),
            1,
        )
        .unwrap(),
    )
    .unwrap();
    CaloricPartialSlip::new(
        thermal,
        model,
        checkpoint,
        [curves[0].curve(), curves[1].curve()],
        DESTINATIONS,
        HeatPartition::new(0.25, 0.75, 0.0).unwrap(),
    )
    .unwrap()
}
fn controls() -> CaloricContactSolve {
    CaloricContactSolve {
        maximum_cells: 2,
        maximum_contact_trials: 600,
        maximum_newton_iterations: 32,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 2e-8,
        thermal_energy_tolerance_j: 4e-8,
        total_energy_tolerance_j: 1e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn motion(du: f64) -> PartialSlipKinematics {
    PartialSlipKinematics {
        creepage: [du, 0.0],
        rolling_speed_mps: 1.0,
        torsional_spin_rad_per_s: 0.0,
        contact_spin_rad_per_s: 0.0,
        dt_s: 1.0,
    }
}
fn work(name: &str) -> GeneralizedWorkOwnership {
    GeneralizedWorkOwnership::new("patch", name, "longitudinal", "lateral", "twist").unwrap()
}
fn other(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}
fn curves() -> [ResolvedHeatCapacityCurve; 2] {
    [
        common::material(1.0, 10.0, 0.0).1,
        common::material(1.0, 20.0, 0.0).1,
    ]
}

#[test]
fn retained_elastic_loading_cools_each_finite_body_by_its_own_entropy_change() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = contact(&thermal, &curves, model(1000.0, slopes()), [0.1, 0.0, 0.0]);
    let step = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("elastic"),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let r = step.report();
    let entropy_change = [
        0.5 * (0.2f64.powi(2) - 0.1f64.powi(2)),
        0.25 * (0.2f64.powi(2) - 0.1f64.powi(2)),
    ];
    for (i, capacity) in [10.0, 20.0].into_iter().enumerate() {
        // Independent caloric/entropy equation: C(T-T0)=-T*Delta(s).
        let expected = 300.0 * capacity / (capacity + entropy_change[i]);
        near(
            step.thermal().temperature_k(DESTINATIONS[i]).unwrap(),
            expected,
        );
        near(r.reversible_heat_j[i], capacity * (expected - 300.0));
        near(
            r.thermodynamics.entropy_j_k[i] - old.thermodynamics().entropy_j_k[i],
            entropy_change[i],
        );
        near(r.physical_heat_j[i], 0.0);
    }
    let t = step.contact().temperature_k();
    let k = 1000.0 - (t[0] - 300.0) - 0.5 * (t[1] - 300.0);
    near(r.contact.tangent_force_n[0], -k * 0.2);
    near(r.thermodynamics.free_energy_j, 0.5 * k * 0.2f64.powi(2));
    // U uses k-T_a*k_a'-T_b*k_b'=1450 N/m, not current k.
    near(
        r.thermodynamics.internal_energy_j,
        0.5 * 1450.0 * 0.2f64.powi(2),
    );
    near(
        r.contact_internal_energy_change_j,
        0.5 * 1450.0 * (0.2f64.powi(2) - 0.1f64.powi(2)),
    );
    near(r.total_energy_residual_j, 0.0);
    assert!(r.contact.generalized_work.numerical_dissipation_j > 0.0);
    assert_eq!(old.checkpoint().state().accepted_steps(), 1);
    assert_eq!(step.contact().checkpoint().state().accepted_steps(), 2);
}

#[test]
fn heating_a_loaded_stationary_contact_does_not_turn_coefficient_work_into_heat() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = contact(&thermal, &curves, model(1000.0, slopes()), [0.1, 0.0, 0.2]);
    let step = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.0),
            &work("heat"),
            controls(),
            &[200.0, 100.0],
            |_| {
                Ok::<_, Infallible>(StepCoefficients {
                    sources_w: BTreeMap::from([(DESTINATIONS[0], 100.0), (DESTINATIONS[1], 200.0)]),
                    conductances_w_k: BTreeMap::new(),
                    payload: "independent heaters",
                })
            },
            || false,
        )
        .unwrap();
    let r = step.report();
    near(step.contact().temperature_k()[0], 310.0);
    near(step.contact().temperature_k()[1], 310.0);
    assert!(r.contact.generalized_work.coefficient_work_j < 0.0);
    for i in 0..2 {
        near(r.reversible_heat_j[i], 0.0);
        near(r.deposited_heat_j[i], 0.0);
    }
    near(r.contact_internal_energy_change_j, 0.0);
    near(r.contact.generalized_work.work_into_interface_j, 0.0);
    near(r.contact.generalized_work.numerical_dissipation_j, 0.0);
    near(r.external_heat_j, 300.0);
    near(r.thermal.stored_energy_change_j, 300.0);
    near(r.total_energy_residual_j, 0.0);
    assert_eq!(r.source, "independent heaters");
    assert_eq!(
        old.checkpoint().state().region_elastic_displacement_m(),
        step.contact()
            .checkpoint()
            .state()
            .region_elastic_displacement_m()
    );
}

#[test]
fn only_physical_two_region_slip_loss_heats_the_bodies() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = contact(&thermal, &curves, model(0.4, zero_slopes()), [0.0; 3]);
    let step = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("slide"),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let r = step.report();
    // Independent Jenkins sliders: Ci=(0.6,3.4) N, ki=(750,250) N/m.
    // Each finishes at Ci/ki and slips by imposed travel minus this strain.
    let heat = 0.6 * (0.1 - 0.6 / 750.0) + 3.4 * (0.1 - 3.4 / 250.0);
    let numerical = 0.5 * (0.6f64.powi(2) / 750.0 + 3.4f64.powi(2) / 250.0);
    near(r.contact.dissipation.heat_j, heat);
    near(
        r.contact.generalized_work.numerical_dissipation_j,
        numerical,
    );
    near(r.thermal.stored_energy_change_j, heat);
    near(r.physical_heat_j[0], 0.25 * heat);
    near(r.physical_heat_j[1], 0.75 * heat);
    near(
        step.contact().temperature_k()[0],
        300.0 + 0.25 * heat / 10.0,
    );
    near(
        step.contact().temperature_k()[1],
        300.0 + 0.75 * heat / 20.0,
    );
    near(r.total_energy_residual_j, 0.0);
    assert!(numerical > 0.0);
    assert_eq!(
        step.contact().checkpoint().state().sliding_regions(),
        [true; 2]
    );
}

#[test]
fn unequal_surface_temperatures_conduct_through_the_original_network() {
    let curves = curves();
    let thermal = initial(&curves, [400.0, 300.0], 1.0);
    let old = contact(&thermal, &curves, model(1000.0, slopes()), [0.1, 0.05, 0.2]);
    let mut k = motion(0.0);
    k.contact_spin_rad_per_s = 0.3;
    let step = old
        .step_candidate(
            &thermal,
            1.0,
            frame().rotated(0.4).unwrap(),
            k,
            &work("conduct"),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    // BE conduction: d_new=d_old/(1+G*dt*(1/Ca+1/Cb)).
    let transfer = 100.0 / (1.0 + 1.0 / 10.0 + 1.0 / 20.0);
    near(step.contact().temperature_k()[0], 400.0 - transfer / 10.0);
    near(step.contact().temperature_k()[1], 300.0 + transfer / 20.0);
    for q in step.report().deposited_heat_j {
        near(q, 0.0);
    }
    near(step.report().contact_internal_energy_change_j, 0.0);
    near(step.report().thermal.stored_energy_change_j, 0.0);
    near(step.report().total_energy_residual_j, 0.0);
}

#[test]
fn refusal_and_cancellation_preserve_loaded_history_and_replay() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = contact(&thermal, &curves, model(1000.0, slopes()), [0.1, 0.0, 0.0]);
    let accepted = old.checkpoint().clone();
    let mut budget = controls();
    budget.maximum_contact_trials = 1;
    assert!(
        old.step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("budget"),
            budget,
            &[0.0; 2],
            other,
            || false
        )
        .is_err()
    );
    assert!(matches!(
        old.step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("cancel"),
            controls(),
            &[0.0; 2],
            other,
            || true
        ),
        Err(CaloricContactError::Cancelled)
    ));
    let calls = std::cell::Cell::new(0);
    let rejected = old.step_candidate(
        &thermal,
        1.0,
        frame(),
        motion(0.1),
        &work("late"),
        controls(),
        &[0.0; 2],
        |state| {
            calls.set(calls.get() + 1);
            if calls.get() > 1 {
                Err("source refused later thermal trial")
            } else {
                Ok(other(state).unwrap())
            }
        },
        || false,
    );
    assert!(matches!(
        rejected,
        Err(CaloricContactError::Source(
            "source refused later thermal trial"
        ))
    ));
    assert_eq!(old.checkpoint(), &accepted);
    near(thermal.time_s(), 0.0);
    let a = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("replay"),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let b = old
        .clone()
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.1),
            &work("replay"),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    assert_eq!(a.contact().checkpoint(), b.contact().checkpoint());
    assert_eq!(a.contact().temperature_k(), b.contact().temperature_k());
    let bad = old.step_candidate(
        &thermal,
        1.0,
        frame(),
        motion(0.0),
        &work("domain"),
        controls(),
        &[0.0; 2],
        |_| {
            Ok::<_, Infallible>(StepCoefficients {
                sources_w: BTreeMap::from([(DESTINATIONS[0], 10000.0)]),
                conductances_w_k: BTreeMap::new(),
                payload: (),
            })
        },
        || false,
    );
    assert!(bad.is_err());
    assert_eq!(old.checkpoint(), &accepted);
}

#[test]
fn affine_domain_and_original_caloric_identity_are_enforced() {
    let mut bad = slopes();
    bad[0].tangential_n_per_m_k = -100.0;
    assert!(
        AffineContactStiffness::new(
            PartialSlipLaw::new("affine", "declared", parameters(100.0)).unwrap(),
            [300.0; 2],
            bad,
            [[200.0, 600.0]; 2]
        )
        .is_err()
    );
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let old = contact(&thermal, &curves, model(100.0, slopes()), [0.1, 0.0, 0.0]);
    let foreign = common::material(1.0, 10.0, 0.0).1;
    assert!(
        CaloricPartialSlip::new(
            &thermal,
            model(100.0, slopes()),
            old.checkpoint().clone(),
            [foreign.curve(), curves[1].curve()],
            DESTINATIONS,
            HeatPartition::new(0.5, 0.5, 0.0).unwrap()
        )
        .is_err()
    );
    assert!(old.model().law_at::<Infallible>([601.0, 300.0]).is_err());
    assert!(old.model().law_at::<Infallible>([f64::NAN, 300.0]).is_err());
}
