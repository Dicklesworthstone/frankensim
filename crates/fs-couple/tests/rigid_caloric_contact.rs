//! Independent native inertia and affine free-energy references.
//! All constitutive coefficients below are explicit analytical fixtures.
use fs_couple::tribothermal::material::{KINETIC_FRICTION_PROPERTY, STATIC_FRICTION_PROPERTY};
use fs_couple::tribothermal::partial_slip::resolved::*;
use fs_couple::tribothermal::partial_slip::rigid::*;
use fs_couple::tribothermal::partial_slip::*;
use fs_couple::tribothermal::{MaterialSlidingContact, SlidingHeatPartition, SlidingInterface};
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::*;
use fs_tribo::{HeatPartition, InputAuthority};
use std::collections::BTreeMap;
use std::convert::Infallible;

const DT: f64 = 0.02;
const AREA: f64 = core::f64::consts::PI * 0.5 * 0.5;
const CELL: [CellId; 2] = [CellId(30), CellId(10)];
const CAP: [f64; 2] = [10.0, 20.0];
const ARMS: [Vec3; 2] = [Vec3::new(0.0, 0.0, 0.5), Vec3::new(0.0, 0.0, -0.5)];
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a:.16e} != {b:.16e}");
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}
fn curves() -> [EquilibriumEnthalpyPhaseCurve; 2] {
    std::array::from_fn(|i| material_caloric(i).1)
}
fn thermal<'a>(curves: &'a [EquilibriumEnthalpyPhaseCurve; 2]) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 100.0 * CAP[i],
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64 - 1.0],
                    end_m: [0.0, 0.0, i as f64],
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
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CELL[0],
            b: CELL[1],
            area_m2: AREA,
            conductance_w_k: 1.5 * AREA,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &curves[0] as &dyn ThermalLaw),
            (MaterialId(1), &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn model(mu: f64, variable: bool) -> AffineContactStiffness {
    AffineContactStiffness::new(
        PartialSlipLaw::new(
            "affine/rigid",
            "synthetic/free-energy",
            PartialSlipParameters {
                static_mu: mu,
                kinetic_mu: mu,
                tangential_stiffness_n_per_m: 1000.0,
                torsional_stiffness_nm_per_rad: 10.0,
                torsional_capacity_factor: 0.1 / 0.125_f64.sqrt(),
                partial_slip_onset_fraction: 0.5,
                partial_slip_hardening_fraction: 0.4,
            },
        )
        .unwrap(),
        [300.0; 2],
        std::array::from_fn(|i| StiffnessDerivative {
            tangential_n_per_m_k: if variable { [-100.0, -50.0][i] } else { 0.0 },
            torsional_nm_per_rad_k: if variable { [-1.0, -0.5][i] } else { 0.0 },
        }),
        [[295.0, 305.0]; 2],
    )
    .unwrap()
}
fn contact<'a>(
    thermal: &DepositionSimulation<'_>,
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    model: AffineContactStiffness,
    strain: [f64; 3],
) -> CaloricPartialSlip<'a> {
    let checkpoint = PartialSlipCheckpoint::new(
        NormalPatchView::new(
            "patch",
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
            "A->B",
            "retained rigid",
            "synthetic/interface",
            NormalPatchAuthority::SyntheticFixture,
        )
        .unwrap(),
        model.law_at::<Infallible>([300.0; 2]).unwrap(),
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
    CaloricPartialSlip::new_with_resolved_interface(
        thermal,
        model,
        checkpoint,
        [&curves[0], &curves[1]],
        7,
        [3.0, 1.0],
    )
    .unwrap()
}
fn mechanics() -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let properties = [
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.5, 1.0)).unwrap(),
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(2.0, 2.5, 1.0)).unwrap(),
    ];
    let bodies = std::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                Vec3::new(0.0, 0.0, if i == 0 { -0.5 } else { 0.5 }),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            if i == 0 {
                Vec3::new(1.0, 0.6, 0.0)
            } else {
                Vec3::ZERO
            },
            Vec3::new(0.0, 0.0, if i == 0 { 0.2 } else { -0.2 }),
        )
        .unwrap()
    });
    (bodies, properties)
}
fn controls() -> RigidCaloricContactSolve {
    let newton = NewtonKrylovConfig {
        absolute_tolerance: 1e-9,
        relative_tolerance: 1e-12,
        forcing_minimum: 1e-4,
        linear_restart: 3,
        ..Default::default()
    };
    RigidCaloricContactSolve {
        contact: ResolvedCaloricContactSolve {
            contact: CaloricContactSolve {
                maximum_cells: 2,
                maximum_contact_trials: 800,
                maximum_newton_iterations: 32,
                difference_step_j: 1e-4,
                cell_energy_tolerance_j: 2e-8,
                thermal_energy_tolerance_j: 4e-8,
                total_energy_tolerance_j: 1e-7,
                newton,
            },
            temperature_scale_j_per_k: 1.0,
            interface_energy_tolerance_j: 2e-8,
        },
        maximum_motion_trials: 300,
        maximum_newton_iterations: 24,
        difference_step_n_s: 1e-4,
        torsional_coordinate_length_m: 0.1,
        geometry_tolerance_m: 1e-12,
        linear_momentum_tolerance_n_s: 2e-8,
        angular_momentum_tolerance_nm_s: 2e-8,
        energy_tolerance_j: 2e-7,
        newton,
    }
}
fn source(
    _: &DepositionSimulation<'_>,
    _: f64,
) -> Result<ResolvedCaloricContactSource<()>, Infallible> {
    Ok(ResolvedCaloricContactSource {
        conductivity_w_m_k: [3.0, 1.0],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    })
}
fn work(name: &str) -> GeneralizedWorkOwnership {
    GeneralizedWorkOwnership::new("patch", name, "x", "y", "spin").unwrap()
}

#[test]
fn g1_native_shear_spin_and_variable_stiffness_match_independent_entropy_reference() {
    let curves = curves();
    let thermal = thermal(&curves);
    let strain = [0.03, 0.01, 0.05];
    let contact = contact(&thermal, &curves, model(1000.0, true), strain);
    let (bodies, properties) = mechanics();
    let candidate = contact
        .step_rigid_resolved_candidate(
            &bodies,
            properties,
            ARMS,
            &thermal,
            DT,
            frame(),
            1.0,
            &work("elastic"),
            controls(),
            [0.0; 3],
            &[0.0; 2],
            300.0,
            source,
            || false,
        )
        .unwrap();
    // Independent elimination of the native diagonal tangent/spin mobility,
    // followed by a two-capacity caloric solve with -T_i Delta s_i.
    let mobility = [
        5.0 / 6.0 + 0.25 * (1.0 / 1.5 + 1.0 / 2.5),
        5.0 / 6.0 + 0.25 * (1.0 + 0.5),
        2.0,
    ];
    let speed = [0.5, 0.3, 0.4];
    let mut t = [300.0; 2];
    let mut next = [0.0; 3];
    for _ in 0..100 {
        let k = 1000.0 - 100.0 * (t[0] - 300.0) - 50.0 * (t[1] - 300.0);
        let stiffness = [k, k, 0.01 * k];
        next = std::array::from_fn(|i| {
            (strain[i] + DT * speed[i]) / (1.0 + 0.5 * DT * DT * mobility[i] * stiffness[i])
        });
        let dq =
            next[0] * next[0] + next[1] * next[1] - strain[0] * strain[0] - strain[1] * strain[1];
        let dr = next[2] * next[2] - strain[2] * strain[2];
        let ds = [50.0 * dq + 0.5 * dr, 25.0 * dq + 0.25 * dr];
        let g = DT * 1.5 * AREA;
        let a = CAP[0] + g + ds[0];
        let b = CAP[1] + g + ds[1];
        let det = a * b - g * g;
        t = [
            (300.0 * CAP[0] * b + g * 300.0 * CAP[1]) / det,
            (300.0 * CAP[1] * a + g * 300.0 * CAP[0]) / det,
        ];
    }
    let report = candidate.report();
    let native = &report.contact.caloric.contact;
    assert_eq!(native.state, PartialSlipStateKind::Sticking);
    for i in 0..2 {
        near(candidate.contact().temperature_k()[i], t[i], 2e-7);
        assert!(t[i] < 300.0);
    }
    for region in 0..2 {
        let actual = native.next_state.region_elastic_displacement_m()[region];
        for i in 0..2 {
            near(actual[i], next[i], 2e-9);
        }
        near(
            native.next_state.region_elastic_twist_rad()[region],
            next[2],
            2e-9,
        );
    }
    let k = 1000.0 - 100.0 * (t[0] - 300.0) - 50.0 * (t[1] - 300.0);
    near(native.tangent_force_n[0], -k * next[0], 2e-6);
    near(native.tangent_force_n[1], -k * next[1], 2e-6);
    near(native.torsional_moment_nm, -0.01 * k * next[2], 2e-7);
    assert!(native.generalized_work.coefficient_work_j.abs() > 1e-3);
    assert!(
        report
            .contact
            .caloric
            .reversible_heat_j
            .iter()
            .all(|q| *q < 0.0)
    );
    assert_eq!(native.dissipation.heat_j, 0.0);
    near(report.mechanical_energy_residual_j, 0.0, 2e-7);
    near(report.total_energy_residual_j, 0.0, 2e-7);
    assert!(report.world_angular_momentum_change_nm_s.norm_squared() < 1e-20);
    for i in 0..2 {
        assert_eq!(candidate.bodies()[i].pose(), bodies[i].pose());
    }
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn g1_gross_slide_heats_only_physical_loss_with_native_off_center_impulses() {
    let curves = curves();
    let thermal = thermal(&curves);
    let contact = contact(&thermal, &curves, model(0.04, false), [0.0; 3]);
    let (bodies, properties) = mechanics();
    let candidate = contact
        .step_rigid_resolved_candidate(
            &bodies,
            properties,
            ARMS,
            &thermal,
            DT,
            frame(),
            1.0,
            &work("sliding"),
            controls(),
            [0.0; 3],
            &[0.0; 2],
            300.0,
            source,
            || false,
        )
        .unwrap();
    let r = candidate.report();
    let step = &r.contact.caloric.contact;
    assert_eq!(step.state, PartialSlipStateKind::GrossSlide);
    assert!(step.dissipation.heat_j > 0.0);
    near(
        r.contact.caloric.thermal.stored_energy_change_j,
        step.dissipation.heat_j,
        1e-7,
    );
    assert!(-r.kinetic_energy_change_j > step.dissipation.heat_j);
    near(r.total_energy_residual_j, 0.0, 2e-7);
    // Original body moments include both contact-arm torque and spin couple.
    for i in 0..2 {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let expected = ARMS[i]
            .cross(r.applied_impulse_on_a_n_s)
            .add(r.applied_spin_impulse_on_a_nm_s)
            .scale(sign);
        let actual = candidate.bodies()[i]
            .angular_momentum_body()
            .sub(bodies[i].angular_momentum_body());
        near(actual.sub(expected).norm_squared(), 0.0, 1e-20);
    }
}

#[test]
fn g0_refusal_keeps_all_original_histories_and_retry_matches() {
    let curves = curves();
    let thermal = thermal(&curves);
    let contact = contact(&thermal, &curves, model(1000.0, true), [0.03, 0.01, 0.05]);
    let (bodies, properties) = mechanics();
    let original = contact.checkpoint().clone();
    for mode in 0..4 {
        let mut config = controls();
        if mode == 0 {
            config.maximum_motion_trials = 1;
        }
        let mut arms = ARMS;
        if mode == 1 {
            arms[0].z += 0.01;
        }
        let calls = std::cell::Cell::new(0usize);
        let failed = contact.step_rigid_resolved_candidate(
            &bodies,
            properties,
            arms,
            &thermal,
            DT,
            frame(),
            1.0,
            &work("retry"),
            config,
            [0.0; 3],
            &[0.0; 2],
            300.0,
            |s, t| {
                calls.set(calls.get() + 1);
                if mode == 2 && calls.get() > 8 {
                    Err("material refused")
                } else {
                    source(s, t).map_err(|never| match never {})
                }
            },
            || mode == 3 && calls.get() > 8,
        );
        assert!(failed.is_err());
        assert_eq!(contact.checkpoint(), &original);
        assert_eq!(contact.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
    }
    let a = contact
        .step_rigid_resolved_candidate(
            &bodies,
            properties,
            ARMS,
            &thermal,
            DT,
            frame(),
            1.0,
            &work("retry"),
            controls(),
            [0.0; 3],
            &[0.0; 2],
            300.0,
            source,
            || false,
        )
        .unwrap();
    let b = contact
        .step_rigid_resolved_candidate(
            &bodies,
            properties,
            ARMS,
            &thermal,
            DT,
            frame(),
            1.0,
            &work("retry"),
            controls(),
            [0.0; 3],
            &[0.0; 2],
            300.0,
            source,
            || false,
        )
        .unwrap();
    assert_eq!(a.bodies(), b.bodies());
    assert_eq!(a.contact().checkpoint(), b.contact().checkpoint());
    assert_eq!(a.contact().temperature_k(), b.contact().temperature_k());
}

fn material_caloric(index: usize) -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
            .unwrap()
            .with("lot", 7.0)
            .unwrap()
    };
    let heat = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
    let density = PropertyKey::new("density", Density::DIMS);
    let mut claims = ClaimSet::new();
    for (key, value) in [(heat.clone(), CAP[index]), (density.clone(), 1.0)] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (600.0, value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        200.0,
                        600.0,
                    )
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic native caloric fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: format!("authored contact cell {index}"),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &heat,
            density: &density,
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(600.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, caloric.curve().clone())
}

#[test]
fn g1_pinned_mu_and_affine_storage_match_native_two_temperature_reference() {
    let (a, ca) = material_caloric(0);
    let (b, cb) = material_caloric(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = thermal(&curves);
    const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
    let point = QueryPoint::new()
        .with_quantity("Ta", QuantitySpec::dimensional(Temperature::DIMS), 290.0)
        .unwrap()
        .with_quantity("Tb", QuantitySpec::dimensional(Temperature::DIMS), 310.0)
        .unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 7.0)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.01)
        .unwrap()
        .with("lot", 7.0)
        .unwrap();
    let mut claims = ClaimSet::new();
    for name in [STATIC_FRICTION_PROPERTY, KINETIC_FRICTION_PROPERTY] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Curve {
                    abscissa: "Ta".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(295.0, 0.65), (305.0, 0.15)],
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "Ta",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        295.0,
                        305.0,
                    )
                    .with_quantity(
                        "Tb",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        295.0,
                        305.0,
                    )
                    .with_quantity(
                        "pressure",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        0.0,
                        1000.0,
                    )
                    .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0, 10.0)
                    .with("lot", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic temperature weakening".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "A finish".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "B finish".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic thermal fixture".into(),
            history: "retained ordered surfaces".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let data = SlidingInterface {
        point: ThermalInterfacePoint {
            interface: &interface,
            materials: [&cards[0], &cards[1]],
            context: &point,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        caloric: [&curves[0], &curves[1]],
        cells: CELL,
        nominal_area_m2: AREA,
        pressure_axis: "pressure",
        speed_axis: "speed",
        accumulated_distance_axis: None,
        tangent_frame_id: "native material frame",
        selection: MaterialPropertySelection::SingleClaimOnly,
        // Resolved routing must derive its actual 3:1 conductance partition;
        // this legacy metadata deliberately supplies the opposite fraction.
        partition: SlidingHeatPartition::Prescribed(HeatPartition::new(0.25, 0.75, 0.0).unwrap()),
        outside_destination: None,
    };
    let sampler = MaterialSlidingContact::new(
        &thermal,
        data.clone(),
        10.0,
        5.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let bound = sampler.interface().interface();
    let stiffness = model(0.4, true);
    let checkpoint = PartialSlipCheckpoint::new(
        NormalPatchView::new(
            "patch",
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
            bound.ordered_system_id(),
            bound.history_id(),
            bound.provenance().source_id(),
            NormalPatchAuthority::SyntheticFixture,
        )
        .unwrap(),
        stiffness.reference_law().clone(),
        PartialSlipState::from_region_checkpoint(
            [[0.001, 0.0]; 2],
            [0.0; 2],
            [false; 2],
            Some(frame()),
            7,
        )
        .unwrap(),
    )
    .unwrap();
    let contact = CaloricPartialSlip::new_with_resolved_interface(
        &thermal,
        stiffness,
        checkpoint,
        [&curves[0], &curves[1]],
        7,
        [3.0, 1.0],
    )
    .unwrap()
    .with_material_friction(&thermal, data, 5.0, InputAuthority::SyntheticFixture)
    .unwrap();
    let (mut bodies, properties) = mechanics();
    for i in 0..2 {
        bodies[i] = RigidBodyState::new(
            bodies[i].pose(),
            if i == 0 {
                Vec3::new(10.0, 0.0, 0.0)
            } else {
                Vec3::ZERO
            },
            Vec3::ZERO,
        )
        .unwrap();
    }
    let external = |_: &DepositionSimulation<'_>,
                    _: f64|
     -> Result<ResolvedCaloricContactSource<()>, Infallible> {
        Ok(ResolvedCaloricContactSource {
            conductivity_w_m_k: [3.0, 1.0],
            transport: StepCoefficients {
                sources_w: BTreeMap::from([(CELL[0], 1000.0)]),
                conductances_w_k: BTreeMap::new(),
                payload: (),
            },
        })
    };
    let next = contact
        .step_rigid_resolved_candidate(
            &bodies,
            properties,
            ARMS,
            &thermal,
            DT,
            frame(),
            1.0,
            &work("pinned material"),
            controls(),
            [0.0; 3],
            &[0.0; 2],
            300.0,
            external,
            || false,
        )
        .unwrap();
    // Independent scalar native mobility, two fully sliding spring caps and
    // two caloric balances. No constitutive or thermal trial calls enter this
    // reference. The source uses Ta, not the massless interface temperature.
    let mobility = 5.0 / 6.0 + 0.25 * (1.0 / 1.5 + 1.0 / 2.5);
    let cap_square_weight = 0.3 * 0.3 / 0.6 + 0.7 * 0.7 / 0.4;
    let mut t = [300.0; 2];
    for _ in 0..100 {
        let force = 10.0 * (0.4 - 0.05 * (t[0] - 300.0));
        let k = 1000.0 - 100.0 * (t[0] - 300.0) - 50.0 * (t[1] - 300.0);
        let travel = DT * (5.0 - 0.5 * DT * mobility * force);
        let heat = force * (0.001 + travel) - force * force * cap_square_weight / k;
        let dq = force * force * cap_square_weight / (k * k) - 0.001 * 0.001;
        let ds = [50.0 * dq, 25.0 * dq];
        let g = DT * 1.5 * AREA;
        let aa = CAP[0] + g + ds[0];
        let bb = CAP[1] + g + ds[1];
        let determinant = aa * bb - g * g;
        let rhs = [
            300.0 * CAP[0] + 0.75 * heat + 20.0,
            300.0 * CAP[1] + 0.25 * heat,
        ];
        t = [
            (rhs[0] * bb + g * rhs[1]) / determinant,
            (rhs[1] * aa + g * rhs[0]) / determinant,
        ];
    }
    let mu = 0.4 - 0.05 * (t[0] - 300.0);
    let force = 10.0 * mu;
    let k = 1000.0 - 100.0 * (t[0] - 300.0) - 50.0 * (t[1] - 300.0);
    let speed = 5.0 - 0.5 * DT * mobility * force;
    let heat = force * (0.001 + DT * speed) - force * force * cap_square_weight / k;
    let r = next.report();
    let step = &r.contact.caloric.contact;
    assert_eq!(step.state, PartialSlipStateKind::GrossSlide);
    for i in 0..2 {
        near(next.contact().temperature_k()[i], t[i], 2e-7);
    }
    near(step.tangent_force_n[0], -force, 2e-7);
    near(step.dissipation.heat_j, heat, 2e-8);
    near(r.midpoint_relative_velocity_mps[0], speed, 2e-8);
    near(
        next.contact().checkpoint().law().parameters().static_mu,
        mu,
        2e-8,
    );
    assert!(mu < 0.31);
    assert!(k < 810.0);
    assert!(step.generalized_work.coefficient_work_j.abs() > 1e-5);
    near(r.total_energy_residual_j, 0.0, 2e-7);
    near(r.contact.caloric.external_heat_j, 20.0, 2e-8);
    let sample = r.contact.caloric.material.as_ref().unwrap();
    near(sample.speed_mps, speed, 2e-8);
    near(sample.resolved.state_coordinate("Ta").unwrap(), t[0], 2e-7);
    near(sample.resolved.state_coordinate("Tb").unwrap(), t[1], 2e-7);
    near(sample.partition.surface_a(), 0.75, 1e-12);
    near(sample.partition.surface_b(), 0.25, 1e-12);
    assert!((r.contact.interface_temperature_k - t[0]).abs() > 0.1);
    let MaterialPropertySelection::PinnedByProperty(pins) =
        &contact.material_friction_binding().unwrap().selection
    else {
        panic!("material claims must remain pinned")
    };
    for property in sample.resolved.properties() {
        assert_eq!(
            pins.iter()
                .find(|(name, _)| name == property.requirement().name())
                .unwrap()
                .1,
            property.answer().receipt.selected
        );
    }
    assert_eq!(contact.checkpoint().state().accepted_steps(), 7);
    assert_eq!(next.contact().checkpoint().state().accepted_steps(), 8);
    // Actual spin remains outside this source's stated applicability even
    // though the general native owner can integrate shear and spin together.
    bodies[0] = RigidBodyState::new(
        bodies[0].pose(),
        bodies[0].linear_momentum_world(),
        Vec3::new(0.0, 0.0, 0.1),
    )
    .unwrap();
    assert!(
        contact
            .step_rigid_resolved_candidate(
                &bodies,
                properties,
                ARMS,
                &thermal,
                DT,
                frame(),
                1.0,
                &work("unowned source spin"),
                controls(),
                [0.0; 3],
                &[0.0; 2],
                300.0,
                external,
                || false
            )
            .is_err()
    );
}
