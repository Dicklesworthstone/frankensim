//! Original three-body inertia and three finite caloric cells; no effective bodies.
use fs_couple::tribothermal::partial_slip::resolved::ResolvedCaloricContactSolve;
use fs_couple::tribothermal::partial_slip::rigid::network::*;
use fs_couple::tribothermal::partial_slip::*;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::*;
use std::cell::Cell;
use std::collections::BTreeMap;
use std::convert::Infallible;

const AREA: f64 = core::f64::consts::PI * 0.25;
const CAP: [f64; 3] = [10.0, 20.0, 30.0];
const IDS: [CellId; 3] = [CellId(70), CellId(10), CellId(30)];
const LINKS: [u64; 2] = [9, 4];
const PATCH: [&str; 2] = ["lower", "upper"];
const ARMS: [Vec3; 2] = [Vec3::new(0.0, 0.0, 0.5), Vec3::new(0.0, 0.0, -0.5)];
const DERIVATIVE: [[f64; 2]; 2] = [[-100.0, -50.0], [-50.0, -100.0]];
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a:.16e} != {b:.16e}");
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}
fn curves() -> [EquilibriumEnthalpyPhaseCurve; 3] {
    std::array::from_fn(|i| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            fs_blake3::hash_bytes(&[1 + i as u8]),
            SolidLiquidPhase::Solid,
            [200.0, 400.0]
                .map(|t| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: CAP[i] * (t - 200.0),
                    temperature_k: t,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1.0,
                })
                .to_vec(),
        )
        .unwrap()
    })
}
fn thermal(curves: &[EquilibriumEnthalpyPhaseCurve; 3]) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for i in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 100.0 * CAP[i],
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64 - 1.5],
                    end_m: [0.0, 0.0, i as f64 - 0.5],
                    volume_m3: AREA,
                    surface_area_m2: 2.0 * AREA + core::f64::consts::PI,
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
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINKS[i],
                a: IDS[i],
                b: IDS[i + 1],
                area_m2: AREA,
                conductance_w_k: AREA,
            }),
        });
    }
    let laws: Vec<_> = curves
        .iter()
        .enumerate()
        .map(|(i, c)| (MaterialId(i as u64), c as &dyn ThermalLaw))
        .collect();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 32,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.02,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn model(index: usize, mu: f64, variable: bool) -> AffineContactStiffness {
    AffineContactStiffness::new(
        PartialSlipLaw::new(
            "affine/network",
            "synthetic/two-contact-free-energy",
            PartialSlipParameters {
                static_mu: mu,
                kinetic_mu: mu,
                tangential_stiffness_n_per_m: 1000.0,
                torsional_stiffness_nm_per_rad: 10.0,
                torsional_capacity_factor: 0.25,
                partial_slip_onset_fraction: 0.5,
                partial_slip_hardening_fraction: 0.4,
            },
        )
        .unwrap(),
        [300.0; 2],
        std::array::from_fn(|i| StiffnessDerivative {
            tangential_n_per_m_k: if variable { DERIVATIVE[index][i] } else { 0.0 },
            torsional_nm_per_rad_k: if variable {
                0.01 * DERIVATIVE[index][i]
            } else {
                0.0
            },
        }),
        [[295.0, 305.0]; 2],
    )
    .unwrap()
}
fn contacts<'a>(
    thermal: &DepositionSimulation<'_>,
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 3],
    mu: f64,
    variable: bool,
    strain: [f64; 2],
) -> [CaloricPartialSlip<'a>; 2] {
    std::array::from_fn(|i| {
        let model = model(i, mu, variable);
        let checkpoint = PartialSlipCheckpoint::new(
            NormalPatchView::new(
                PATCH[i],
                "held ellipse",
                "synthetic/normal",
                NormalPatchAuthority::SyntheticFixture,
                10.0,
                0.5,
                0.5,
                0.125,
            )
            .unwrap(),
            PartialSlipInterface::new(
                PATCH[i],
                "retained native island",
                "synthetic/interface",
                NormalPatchAuthority::SyntheticFixture,
            )
            .unwrap(),
            model.reference_law().clone(),
            PartialSlipState::from_region_checkpoint(
                [[strain[i], 0.0]; 2],
                [0.0; 2],
                [false; 2],
                Some(frame()),
                1,
            )
            .unwrap(),
        )
        .unwrap();
        CaloricPartialSlip::new_with_resolved_interface(
            thermal,
            model,
            checkpoint,
            [&curves[i], &curves[i + 1]],
            LINKS[i],
            [1.0; 2],
        )
        .unwrap()
    })
}
fn mechanics() -> ([RigidBodyState; 3], [MassProperties; 3]) {
    let properties = [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(); 3];
    let bodies = std::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                Vec3::new(0.0, 0.0, i as f64 - 1.0),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::new(2.0 * [0.6, 0.2, -0.3][i], 0.0, 0.0),
            Vec3::ZERO,
        )
        .unwrap()
    });
    (bodies, properties)
}
fn work() -> [GeneralizedWorkOwnership; 2] {
    PATCH.map(|patch| {
        GeneralizedWorkOwnership::new(patch, "native island", "x", "y", "spin").unwrap()
    })
}
fn inputs<'s, 'a>(
    contacts: &'s [CaloricPartialSlip<'a>; 2],
    work: &'s [GeneralizedWorkOwnership; 2],
) -> [RigidCaloricNetworkContact<'s, 'a>; 2] {
    std::array::from_fn(|i| RigidCaloricNetworkContact {
        contact: &contacts[i],
        bodies: [i, i + 1],
        arms_body_m: ARMS,
        frame: frame(),
        creepage_reference_speed_mps: 1.0,
        ownership: &work[i],
        impulse_guess_n_s: [0.0; 3],
        interface_temperature_guess_k: 300.0,
    })
}
fn controls() -> RigidCaloricNetworkSolve {
    RigidCaloricNetworkSolve {
        contact: ResolvedCaloricContactSolve {
            contact: CaloricContactSolve {
                maximum_cells: 3,
                maximum_contact_trials: 1200,
                maximum_newton_iterations: 32,
                difference_step_j: 1e-4,
                cell_energy_tolerance_j: 2e-8,
                thermal_energy_tolerance_j: 4e-8,
                total_energy_tolerance_j: 2e-7,
                newton: NewtonKrylovConfig {
                    absolute_tolerance: 2e-10,
                    relative_tolerance: 1e-12,
                    forcing_minimum: 1e-4,
                    linear_restart: 12,
                    ..Default::default()
                },
            },
            temperature_scale_j_per_k: 1.0,
            interface_energy_tolerance_j: 2e-8,
        },
        maximum_bodies: 3,
        maximum_contacts: 2,
        impulse_coordinate_speed_mps: 1.0,
        torsional_coordinate_length_m: 0.1,
        geometry_tolerance_m: 1e-12,
        linear_momentum_tolerance_n_s: 2e-8,
        angular_momentum_tolerance_nm_s: 2e-8,
        mechanical_energy_tolerance_j: 2e-7,
    }
}
fn source(links: &[u64], external_w: [f64; 3]) -> RigidCaloricNetworkSource<()> {
    RigidCaloricNetworkSource {
        interfaces: links
            .iter()
            .map(|id| CaloricNetworkInterfaceSource {
                link_id: *id,
                conductivity_w_m_k: [1.0; 2],
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: IDS.into_iter().zip(external_w).collect(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}

// Independent tridiagonal three-capacity balance including BOTH contacts'
// surface entropy in the middle cell. This never calls a production solver.
fn caloric_reference(dt: f64, ds: [f64; 3], supplied_j: [f64; 3]) -> [f64; 3] {
    let g = dt * AREA;
    let mut b = [
        CAP[0] + g + ds[0],
        CAP[1] + 2.0 * g + ds[1],
        CAP[2] + g + ds[2],
    ];
    let mut r = std::array::from_fn::<_, 3, _>(|i| CAP[i] * 300.0 + supplied_j[i]);
    b[1] -= g * g / b[0];
    r[1] += g * r[0] / b[0];
    b[2] -= g * g / b[1];
    r[2] += g * r[1] / b[1];
    let t2 = r[2] / b[2];
    let t1 = (r[1] + g * t2) / b[1];
    [(r[0] + g * t1) / b[0], t1, t2]
}
fn native_reference(bodies: &[RigidBodyState; 3], impulses: [f64; 2]) -> ([f64; 3], [f64; 3]) {
    let [a, b] = impulses;
    (
        [
            bodies[0].linear_momentum_world().x + a,
            bodies[1].linear_momentum_world().x - a + b,
            bodies[2].linear_momentum_world().x - b,
        ],
        [0.5 * a, 0.5 * (a + b), 0.5 * b],
    )
}

#[test]
fn g1_shared_native_inertia_and_affine_entropy_match_independent_three_cell_reference() {
    let curves = curves();
    let thermal = thermal(&curves);
    let old_strain = [0.001, -0.0007];
    let contacts = contacts(&thermal, &curves, 10.0, true, old_strain);
    let (bodies, properties) = mechanics();
    let work = work();
    let inputs = inputs(&contacts, &work);
    let dt = 0.02;
    let external_w = [100.0, 250.0, -50.0];
    let calls = Cell::new(0);
    let candidate = step_rigid_caloric_network(
        &bodies,
        &properties,
        &inputs,
        &thermal,
        dt,
        controls(),
        &[0.0; 3],
        |_, links, temperatures| -> Result<_, Infallible> {
            calls.set(calls.get() + 1);
            assert_eq!(links, [4, 9]);
            assert_eq!(temperatures.len(), 2);
            Ok(source(links, external_w))
        },
        || false,
    )
    .unwrap();
    let mut t = [300.0; 3];
    let mut impulse = [0.0; 2];
    let mut z = [0.0; 2];
    for _ in 0..100 {
        let k: [f64; 2] = std::array::from_fn(|i| {
            1000.0 + DERIVATIVE[i][0] * (t[i] - 300.0) + DERIVATIVE[i][1] * (t[i + 1] - 300.0)
        });
        // W = [[1.5,-.25],[-.25,1.5]] comes from all ORIGINAL masses,
        // transverse inertias and half-metre contact arms, including the shared body.
        let a = 1.0 + 0.75 * dt * dt * k[0];
        let b = -0.125 * dt * dt * k[0];
        let c = -0.125 * dt * dt * k[1];
        let d = 1.0 + 0.75 * dt * dt * k[1];
        let rhs = [
            -dt * k[0] * (old_strain[0] + 0.4 * dt),
            -dt * k[1] * (old_strain[1] + 0.5 * dt),
        ];
        impulse = [
            (d * rhs[0] - b * rhs[1]) / (a * d - b * c),
            (a * rhs[1] - c * rhs[0]) / (a * d - b * c),
        ];
        z = std::array::from_fn(|i| -impulse[i] / (dt * k[i]));
        let dq = std::array::from_fn::<_, 2, _>(|i| z[i] * z[i] - old_strain[i] * old_strain[i]);
        let ds = [50.0 * dq[0], 25.0 * (dq[0] + dq[1]), 50.0 * dq[1]];
        t = caloric_reference(dt, ds, external_w.map(|p| dt * p));
    }
    let r = candidate.report();
    assert_eq!(calls.get(), r.trials);
    near(r.external_heat_j, 6.0, 1e-12);
    let (momentum, angular) = native_reference(&bodies, impulse);
    for i in 0..3 {
        near(
            candidate.thermal().temperature_k(IDS[i]).unwrap(),
            t[i],
            2e-7,
        );
        near(
            candidate.bodies()[i].linear_momentum_world().x,
            momentum[i],
            2e-9,
        );
        near(
            candidate.bodies()[i].angular_momentum_body().y,
            angular[i],
            2e-9,
        );
        assert_eq!(candidate.bodies()[i].pose(), bodies[i].pose());
    }
    for (report, accepted) in r.interfaces.iter().zip(candidate.contacts()) {
        let i = LINKS
            .iter()
            .position(|id| *id == report.heat.link_id)
            .unwrap();
        assert_eq!(report.contact.state, PartialSlipStateKind::Sticking);
        for region in 0..2 {
            near(
                report.contact.next_state.region_elastic_displacement_m()[region][0],
                z[i],
                2e-9,
            );
        }
        near(report.applied_impulse_on_a_n_s.x, impulse[i], 2e-9);
        assert!(report.reversible_heat_j.iter().all(|q| *q < 0.0));
        assert_eq!(report.contact.dissipation.heat_j, 0.0);
        assert!(report.contact.generalized_work.coefficient_work_j.abs() > 1e-6);
        assert_eq!(accepted.time_s(), dt);
    }
    near(r.mechanical_energy_residual_j, 0.0, 2e-7);
    near(r.total_energy_residual_j, 0.0, 2e-7);
    assert!(r.world_angular_momentum_change_nm_s.norm_squared() < 1e-24);
    // Unchanged inputs demonstrate that accepting requires the whole candidate.
    assert_eq!(thermal.time_s(), 0.0);
    assert!(contacts.iter().all(|c| c.time_s() == 0.0));
}

#[test]
fn g1_two_gross_sliding_contacts_partition_only_physical_heat_into_shared_cell() {
    let curves = curves();
    let thermal = thermal(&curves);
    let contacts = contacts(&thermal, &curves, 0.04, true, [0.0; 2]);
    let (bodies, properties) = mechanics();
    let work = work();
    let inputs = inputs(&contacts, &work);
    let dt = 0.02;
    let candidate = step_rigid_caloric_network(
        &bodies,
        &properties,
        &inputs,
        &thermal,
        dt,
        controls(),
        &[0.0; 3],
        |_, links, _| -> Result<_, Infallible> { Ok(source(links, [0.0; 3])) },
        || false,
    )
    .unwrap();
    let impulse = [-0.4 * dt; 2];
    let (momentum, angular) = native_reference(&bodies, impulse);
    let slip = [0.4 + 0.625 * impulse[0], 0.5 + 0.625 * impulse[1]];
    let mut t = [300.0; 3];
    let mut heat = [0.0; 2];
    let mut q = [0.0; 2];
    for _ in 0..80 {
        for i in 0..2 {
            let k =
                1000.0 + DERIVATIVE[i][0] * (t[i] - 300.0) + DERIVATIVE[i][1] * (t[i + 1] - 300.0);
            // Independent two Jenkins regions: stiffness weights .6/.4,
            // Coulomb shares .3/.7, and positive plastic slip in both regions.
            let z = [0.3 * 0.4 / (0.6 * k), 0.7 * 0.4 / (0.4 * k)];
            q[i] = 0.6 * z[0] * z[0] + 0.4 * z[1] * z[1];
            heat[i] = 0.4 * (dt * slip[i] - 0.3 * z[0] - 0.7 * z[1]);
        }
        t = caloric_reference(
            dt,
            [50.0 * q[0], 25.0 * (q[0] + q[1]), 50.0 * q[1]],
            [0.5 * heat[0], 0.5 * (heat[0] + heat[1]), 0.5 * heat[1]],
        );
    }
    for i in 0..3 {
        near(
            candidate.thermal().temperature_k(IDS[i]).unwrap(),
            t[i],
            2e-7,
        );
        near(
            candidate.bodies()[i].linear_momentum_world().x,
            momentum[i],
            2e-9,
        );
        near(
            candidate.bodies()[i].angular_momentum_body().y,
            angular[i],
            2e-9,
        );
    }
    for report in &candidate.report().interfaces {
        let i = LINKS
            .iter()
            .position(|id| *id == report.heat.link_id)
            .unwrap();
        assert_eq!(report.contact.state, PartialSlipStateKind::GrossSlide);
        near(report.contact.dissipation.heat_j, heat[i], 1e-9);
        for side in 0..2 {
            near(report.physical_heat_j[side], 0.5 * heat[i], 1e-9);
        }
        assert!(report.contact.generalized_work.numerical_dissipation_j > 0.0);
        assert!(report.thermodynamics.internal_energy_j > report.thermodynamics.free_energy_j);
    }
    near(candidate.report().total_energy_residual_j, 0.0, 2e-7);
}

#[test]
fn g1_retained_joint_time_steps_converge_to_the_original_coupled_oscillator() {
    let curves = curves();
    let (initial, properties) = mechanics();
    let work = work();
    let time = 0.04;
    let symmetric_frequency = (1250.0_f64).sqrt();
    let antisymmetric_frequency = (1750.0_f64).sqrt();
    let reference = [
        0.45 * (symmetric_frequency * time).cos() - 0.05 * (antisymmetric_frequency * time).cos(),
        0.45 * (symmetric_frequency * time).cos() + 0.05 * (antisymmetric_frequency * time).cos(),
    ];
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut thermal = thermal(&curves);
        let mut contacts = contacts(&thermal, &curves, 1000.0, false, [0.0; 2]);
        let mut bodies = initial.to_vec();
        for step in 1..=steps {
            let input = inputs(&contacts, &work);
            let candidate = step_rigid_caloric_network(
                &bodies,
                &properties,
                &input,
                &thermal,
                time * step as f64 / steps as f64,
                controls(),
                &[0.0; 3],
                |_, links, _| -> Result<_, Infallible> { Ok(source(links, [0.0; 3])) },
                || false,
            )
            .unwrap();
            assert!(
                candidate
                    .report()
                    .interfaces
                    .iter()
                    .all(|r| r.contact.dissipation.heat_j == 0.0)
            );
            let (b, c, t, _) = candidate.into_parts();
            bodies = b;
            thermal = t;
            // Returned contacts are in link order [4,9]; physical lower/upper
            // order is reconstructed by original IDs, never by solve position.
            let mut c = c.into_iter();
            let upper = c.next().unwrap();
            let lower = c.next().unwrap();
            contacts = [lower, upper];
        }
        let actual: [f64; 2] = std::array::from_fn(|i| {
            bodies[i]
                .point_kinematics(properties[i], ARMS[0])
                .unwrap()
                .point_velocity_world
                .x
                - bodies[i + 1]
                    .point_kinematics(properties[i + 1], ARMS[1])
                    .unwrap()
                    .point_velocity_world
                    .x
        });
        errors.push((actual[0] - reference[0]).hypot(actual[1] - reference[1]));
        for id in IDS {
            near(thermal.temperature_k(id).unwrap(), 300.0, 1e-10);
        }
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (1.7..2.3).contains(&ratio),
            "first-order coupled velocity errors: {errors:?}"
        );
    }
    assert!(
        errors[2] < 0.009,
        "finest native velocity error: {errors:?}"
    );
}

#[test]
fn g0_permutation_late_refusals_and_retry_preserve_all_original_owners() {
    let curves = curves();
    let thermal = thermal(&curves);
    let contacts = contacts(&thermal, &curves, 0.04, true, [0.0; 2]);
    let (bodies, properties) = mechanics();
    let work = work();
    let input = inputs(&contacts, &work);
    let resolve = |_: &DepositionSimulation<'_>,
                   links: &[u64],
                   _: &[f64]|
     -> Result<_, &'static str> { Ok(source(links, [0.0; 3])) };
    let run = || {
        step_rigid_caloric_network(
            &bodies,
            &properties,
            &input,
            &thermal,
            0.02,
            controls(),
            &[0.0; 3],
            resolve,
            || false,
        )
        .unwrap()
    };
    let accepted = run();
    let calls = Cell::new(0);
    let failed = step_rigid_caloric_network(
        &bodies,
        &properties,
        &input,
        &thermal,
        0.02,
        controls(),
        &[0.0; 3],
        |_, links, _| {
            calls.set(calls.get() + 1);
            if calls.get() == 4 {
                Err("late shared source refusal")
            } else {
                Ok(source(links, [0.0; 3]))
            }
        },
        || false,
    );
    assert!(failed.is_err());
    assert_eq!(calls.get(), 4);
    let calls = Cell::new(0);
    let cancelled = step_rigid_caloric_network(
        &bodies,
        &properties,
        &input,
        &thermal,
        0.02,
        controls(),
        &[0.0; 3],
        |_, links, _| -> Result<_, Infallible> {
            calls.set(calls.get() + 1);
            Ok(source(links, [0.0; 3]))
        },
        || calls.get() >= 4,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls.get(), 4);
    let mut limited = controls();
    limited.contact.contact.maximum_contact_trials = 1;
    assert!(
        step_rigid_caloric_network(
            &bodies,
            &properties,
            &input,
            &thermal,
            0.02,
            limited,
            &[0.0; 3],
            resolve,
            || false
        )
        .is_err()
    );
    let mut reordered = inputs(&contacts, &work);
    reordered.swap(0, 1);
    let permutation = step_rigid_caloric_network(
        &bodies,
        &properties,
        &reordered,
        &thermal,
        0.02,
        controls(),
        &[0.0; 3],
        resolve,
        || false,
    )
    .unwrap();
    let retry = run();
    for next in [&permutation, &retry] {
        assert_eq!(next.bodies(), accepted.bodies());
        for id in IDS {
            assert_eq!(
                next.thermal().cells()[&id].specific_h_j_kg,
                accepted.thermal().cells()[&id].specific_h_j_kg
            );
        }
        for (a, b) in next.contacts().iter().zip(accepted.contacts()) {
            assert_eq!(a.checkpoint(), b.checkpoint());
            assert_eq!(a.thermodynamics(), b.thermodynamics());
        }
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert!(contacts.iter().all(|c| c.time_s() == 0.0));
}
