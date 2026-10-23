//! Convergent material fibers plus the actual string/body/friction consumer.
use super::*;
use fs_couple::fiber_ribbon::{
    FiberDirection, FiberRibbonGeometry, FiberRibbonResolution, TensionedFiberRibbon,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, axial::YOUNG_VISCOSITY_PROPERTY,
};

fn fiber_material(rho: f64, young: f64, eta: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, young),
        (ETA, DynViscosity::DIMS, eta),
        (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, eta),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", 0.0, 1e8),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "manufactured fiber for analytic/wiring checks, not horsehair measurements"
                            .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-fiber".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn fibers(
    cells: usize,
    modes: usize,
    rho: f64,
    young: f64,
    eta: f64,
    tension: f64,
    rate: u32,
    gate: &CancelGate,
) -> TensionedFiberRibbon {
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 40e-6,
            material: fiber_material(rho, young, eta),
        }],
    )
    .unwrap();
    TensionedFiberRibbon::try_new(
        FiberRibbonGeometry {
            length_m: 0.6,
            width_m: 0.01,
            fiber_count: 180,
            tension_n: tension,
        },
        FiberRibbonResolution {
            cells,
            axial_modes: modes,
            transverse_modes: modes,
        },
        section
            .clone()
            .with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY)
            .unwrap(),
        section.with_kelvin_voigt_bending(ETA).unwrap(),
        rate,
        ModalAcousticTimeBudget::audible_reference(),
        gate,
    )
    .unwrap()
}

fn static_compliance(
    ribbon: &TensionedFiberRibbon,
    station: f64,
    direction: FiberDirection,
) -> f64 {
    ribbon
        .displacement_map(station, direction)
        .unwrap()
        .iter()
        .zip(ribbon.cell_model().modes())
        .map(|(shape, mode)| shape.powi(2) / mode.angular_frequency_rad_s.powi(2))
        .sum()
}

#[test]
fn axial_and_tensioned_normal_continua_converge_to_independent_static_green_functions() {
    let gate = CancelGate::new_clock_free();
    let mut previous = [f64::INFINITY; 2];
    for modes in [7, 15, 31] {
        let ribbon = fibers(3, modes, 1300.0, 3e9, 0.0, 40.0, 1_000_000, &gate);
        let geometry = ribbon.geometry();
        let section = ribbon.axial_material().elastic();
        let share = geometry.fiber_count as f64 / ribbon.resolution().cells as f64;
        let ea = share * section.axial_stiffness_n();
        let ei = share * section.bending_stiffness_n_m2();
        let t = geometry.tension_n / ribbon.resolution().cells as f64;
        let l = geometry.length_m;
        let k = (t / ei).sqrt();
        // Midspan fixed-end axial bar and simply-supported beam under tension.
        let exact = [
            l / (4.0 * ea),
            l / (4.0 * t) - (0.5 * k * l).tanh() / (2.0 * t * k),
        ];
        for (i, direction) in [FiberDirection::Axial, FiberDirection::Normal]
            .into_iter()
            .enumerate()
        {
            let actual = static_compliance(&ribbon, 0.5, direction);
            let relative = (actual / exact[i] - 1.0).abs();
            assert!(
                relative < 0.56 * previous[i],
                "{modes} modes: {direction:?} {relative}"
            );
            assert!(actual < exact[i]); // positive omitted high-frequency compliance
            previous[i] = relative;
        }
    }
    assert!(previous.iter().all(|error| *error < 0.014));
}

#[test]
fn width_refinement_preserves_mass_and_uniform_load_compliance() {
    let gate = CancelGate::new_clock_free();
    let reference = fibers(1, 3, 1300.0, 3e9, 1.0, 40.0, 48000, &gate);
    let reference_c = static_compliance(&reference, 0.37, FiberDirection::Axial);
    for cells in [3, 9, 27] {
        let ribbon = fibers(cells, 3, 1300.0, 3e9, 1.0, 40.0, 48000, &gate);
        assert_eq!(ribbon.mass_kg(), reference.mass_kg());
        assert_eq!(ribbon.cell_model().modes(), reference.cell_model().modes());
        let normalized_c = static_compliance(&ribbon, 0.37, FiberDirection::Axial) / cells as f64;
        assert!((normalized_c / reference_c - 1.0).abs() < 2e-14);
        assert!(
            (ribbon.cell_line_density_kg_m() * cells as f64 / reference.cell_line_density_kg_m()
                - 1.0)
                .abs()
                < 2e-14
        );
    }
}

#[test]
fn material_and_tension_change_the_correct_independent_wave_families() {
    let gate = CancelGate::new_clock_free();
    let a = fibers(3, 3, 1300.0, 3e9, 1.0, 40.0, 48000, &gate);
    let dense = fibers(3, 3, 5200.0, 3e9, 1.0, 40.0, 48000, &gate);
    let taut = fibers(3, 3, 1300.0, 3e9, 1.0, 160.0, 48000, &gate);
    for (original, changed) in a
        .cell_model()
        .modes()
        .iter()
        .zip(dense.cell_model().modes())
    {
        assert!(
            (changed.angular_frequency_rad_s / original.angular_frequency_rad_s - 0.5).abs()
                < 1e-14
        );
    }
    assert_eq!(
        &a.cell_model().modes()[..3],
        &taut.cell_model().modes()[..3]
    );
    for (original, changed) in a.cell_model().modes()[3..]
        .iter()
        .zip(&taut.cell_model().modes()[3..])
    {
        assert!(changed.angular_frequency_rad_s > 1.99 * original.angular_frequency_rad_s);
    }
    let axial = a.displacement_map(0.4, FiberDirection::Axial).unwrap();
    let normal = a.displacement_map(0.4, FiberDirection::Normal).unwrap();
    assert!(axial[..3].iter().any(|x| *x != 0.0));
    assert!(axial[3..].iter().all(|x| *x == 0.0));
    assert!(normal[..3].iter().all(|x| *x == 0.0));
    assert!(normal[3..].iter().any(|x| *x != 0.0));
}

#[test]
fn physical_fibers_receive_opposite_traction_and_change_string_body_pressure() {
    with_cx(|cx, gate| {
        let make = |rho| {
            base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(
                    fibers(3, 2, rho, 3e9, 50.0, 40.0, 48000, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
                .unwrap()
        };
        let mut a = make(1300.0);
        let mut b = make(5200.0);
        let source = a.source().fiber_ribbon().unwrap();
        let first = source.component_start();
        assert!(a.source().bank().material_cavity().is_some());
        for i in 0..3 {
            let surface = a.system().port(i).unwrap().right.as_ref().unwrap();
            assert_eq!(surface.component, first + i);
            assert_eq!(
                surface.shapes,
                source
                    .model()
                    .displacement_map(0.4, FiberDirection::Axial)
                    .unwrap()
            );
            assert_eq!(
                source.normal_stop_receivers()[i].component,
                first + i - a.cards().len()
            );
        }
        let zero = vec![0.0; a.system().mode_count()];
        let mut drives = vec![drive(0.0, 0.0); 3];
        a.source()
            .ribbon()
            .unwrap()
            .project_drive(drive(0.01, 0.01), &mut drives)
            .unwrap();
        let (mut pressure_difference, mut force_difference) = (0.0, 0.0);
        for _ in 0..256 {
            let fa = a.step(&zero, &drives).unwrap().clone();
            let fb = b.step(&zero, &drives).unwrap().clone();
            pressure_difference += (fa.observer_pressure_pa - fb.observer_pressure_pa).abs();
            force_difference += (fa.ports[0].traction_n - fb.ports[0].traction_n).abs();
            assert!(fa.energy_residual_j.abs() <= fa.energy_tolerance_j);
            assert!(fb.energy_residual_j.abs() <= fb.energy_tolerance_j);
        }
        assert!(pressure_difference > 1e-12 && force_difference > 1e-12);
        assert!(a.system().components()[first..].iter().any(|cell| {
            cell.states()[..2]
                .iter()
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0)
        }));
        // A refused contact trial cannot alter string, fiber or body state.
        let before = states(&a);
        let frame = a.system().last_frame().unwrap().clone();
        drives[2].normal_force_n = f64::NAN;
        assert!(a.step(&zero, &drives).is_err());
        assert_eq!(states(&a), before);
        assert_eq!(a.system().last_frame(), Some(&frame));
    });
}

#[test]
fn fiber_tension_participates_in_the_simultaneous_normal_load_and_friction_solve() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::contact_load::BowContactLoad;
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::{
        MaterialStringStop, StringStopSide,
    };
    use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{
        TransverseDirection, TransverseDirections,
    };
    use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
    use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
    scoped_normal(|cx, gate, tension| {
        let original = base(cx, gate, 7800.0);
        let mut members = original.source().strings().to_vec();
        members[0].source = members[0]
            .source
            .clone()
            .with_transverse_polarizations(
                TransverseDirections::try_new([1.0, 0.0], [1.0, 0.0]).unwrap(),
            )
            .unwrap();
        let raw = ReciprocalStringBank::with_body_connections(
            cx,
            members,
            *original.source().gas(),
            original.system().components()[2..].to_vec(),
            original.source().body_connections().to_vec(),
            law(),
            options(),
            gate,
        )
        .unwrap();
        let bank = raw
            .with_tensioned_fiber_ribbon(
                fibers(3, 2, 1300.0, 3e9, 50.0, tension, 48000, gate),
                0.4,
                ribbon_limits(),
                gate,
            )
            .unwrap();
        let normal = bank
            .source()
            .fiber_ribbon()
            .unwrap()
            .normal_stop_receivers();
        let stops = bank
            .source()
            .bows()
            .iter()
            .zip(normal)
            .map(|(bow, body)| MaterialStringStop {
                location: bow.location,
                side: StringStopSide::Along(TransverseDirection::try_new([0.0, 1.0]).unwrap()),
                body: Some(body),
                law: fs_dcontact::Obstacle::new(
                    vec![-1.0],
                    1,
                    1,
                    vec![-1e-6],
                    vec![1.0],
                    5000.0,
                    1.0,
                    "manufactured finite compliance; not measured rosin".into(),
                )
                .unwrap(),
                limits: ModalContactConfig {
                    max_iterations: 96,
                    maximum_force_n: 10.0,
                    maximum_penetration_m: 0.001,
                    force_absolute_tolerance_n: 1e-11,
                    force_relative_tolerance: 1e-9,
                },
            })
            .collect();
        bank.with_string_stops(
            stops,
            FrictionObstacleConfig {
                max_contacts: 3,
                max_setup_terms: 1_000_000,
            },
            gate,
        )
        .unwrap()
        .with_bow_contact_loads(
            (0..3)
                .map(|bow| BowContactLoad { bow, stop: bow })
                .collect(),
            gate,
        )
        .unwrap()
    });
}

fn scoped_normal(make: impl Fn(&Cx<'_>, &CancelGate, f64) -> MultiBowStringBank) {
    with_cx(|cx, gate| {
        let mut loose = make(cx, gate, 10.0);
        let mut taut = make(cx, gate, 100.0);
        let first = loose.source().fiber_ribbon().unwrap().component_start();
        let zero = vec![0.0; loose.system().mode_count()];
        let drives = vec![drive(0.01, 0.0); 3];
        let mut normal_difference = 0.0;
        for _ in 0..128 {
            let a = loose.step(&zero, &drives).unwrap().clone();
            let b = taut.step(&zero, &drives).unwrap().clone();
            for (port, normal) in a
                .ports
                .iter()
                .zip(loose.system().last_obstacle_frames().unwrap())
            {
                assert_eq!(port.normal_force_n, normal.normal_force_n);
            }
            normal_difference += (a.ports[0].normal_force_n - b.ports[0].normal_force_n).abs();
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
        }
        assert!(normal_difference > 1e-10);
        assert!(loose.system().components()[first..].iter().any(|cell| {
            cell.states()[2..]
                .iter()
                .any(|s| s.displacement_m_sqrt_kg != 0.0)
        }));
    });
}
