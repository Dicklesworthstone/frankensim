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

fn strain_config()
-> fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig {
    fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig {
        max_iterations: 128,
        maximum_generalized_force: 1e6,
        force_absolute_tolerance: 1e-10,
        force_relative_tolerance: 1e-9,
    }
}

fn strain_limits() -> fs_couple::fiber_ribbon::strain::FiberStrainLimits {
    fs_couple::fiber_ribbon::strain::FiberStrainLimits {
        maximum_slope: 0.3,
        maximum_strain: 0.1,
    }
}

#[test]
fn fiber_strain_matches_independent_continuum_energy_and_gradient() {
    // G1: all nonzero geometric harmonics below are retained, so numerical
    // integration of EA/2 * (u_x + y_x²/2)² is an independent continuum oracle.
    use fs_phs::Storage;
    let gate = CancelGate::new_clock_free();
    let ribbon = fibers(3, 4, 1300.0, 3e9, 0.0, 40.0, 192000, &gate)
        .with_axial_transverse_strain(strain_limits(), strain_config())
        .unwrap();
    let strain = ribbon.strain().unwrap();
    let q = [5e-9, -2e-9, 0.0, 0.0, 1e-6, -2e-6, 0.0, 0.0];
    let l = ribbon.geometry().length_m;
    let scale = (2.0 / (ribbon.cell_line_density_kg_m() * l)).sqrt();
    let ea = ribbon.axial_material().elastic().axial_stiffness_n()
        * ribbon.geometry().fiber_count as f64
        / ribbon.resolution().cells as f64;
    let mut integral = 0.0;
    for i in 0..2048 {
        let x = (i as f64 + 0.5) / 2048.0;
        let mut axial = 0.0;
        let mut normal = 0.0;
        for n in 1..=4 {
            let derivative = scale * n as f64 * std::f64::consts::PI / l
                * (n as f64 * std::f64::consts::PI * x).cos();
            axial += derivative * q[n - 1];
            normal += derivative * q[4 + n - 1];
        }
        integral += 0.5 * ea * (axial + 0.5 * normal * normal).powi(2) * l / 2048.0;
    }
    assert!((strain.hamiltonian(&q) / integral - 1.0).abs() < 2e-13);
    let mut gradient = [0.0; 8];
    strain.gradient(&q, &mut gradient);
    for i in 0..8 {
        let mut a = q;
        let mut b = q;
        let h = 2e-11;
        a[i] -= h;
        b[i] += h;
        let numerical = (strain.hamiltonian(&b) - strain.hamiltonian(&a)) / (2.0 * h);
        assert!(
            (gradient[i] - numerical).abs() < 1e-7 * gradient[i].abs() + 1e-12,
            "coordinate {i}: {} != {numerical}",
            gradient[i]
        );
    }
}

#[test]
fn axial_enrichment_converges_to_uncondensed_continuum_strain_energy() {
    use fs_phs::Storage;
    let gate = CancelGate::new_clock_free();
    let source = fibers(3, 4, 1300.0, 3e9, 0.0, 40.0, 192000, &gate);
    let l = source.geometry().length_m;
    let amplitudes = [2e-4, -1e-4, 5e-5, 2.5e-5];
    let ea = source.axial_material().elastic().axial_stiffness_n() * 60.0;
    let mut integral = 0.0;
    for i in 0..2048 {
        let x = (i as f64 + 0.5) / 2048.0;
        let slope: f64 = amplitudes
            .iter()
            .enumerate()
            .map(|(n, amplitude)| {
                let k = (n + 1) as f64 * std::f64::consts::PI;
                amplitude * k / l * (k * x).cos()
            })
            .sum();
        integral += 0.125 * ea * slope.powi(4) * l / 2048.0;
    }
    let mut previous_error = f64::INFINITY;
    for axial_modes in [2, 4, 8] {
        let ribbon = TensionedFiberRibbon::try_new(
            source.geometry(),
            FiberRibbonResolution {
                cells: 3,
                axial_modes,
                transverse_modes: 4,
            },
            source.axial_material().clone(),
            source.bending_material().clone(),
            192000,
            ModalAcousticTimeBudget::audible_reference(),
            &gate,
        )
        .unwrap()
        .with_axial_transverse_strain(strain_limits(), strain_config())
        .unwrap();
        let normalization = (0.5 * ribbon.cell_line_density_kg_m() * l).sqrt();
        let mut q = vec![0.0; axial_modes + amplitudes.len()];
        for (slot, amplitude) in q[axial_modes..].iter_mut().zip(amplitudes) {
            *slot = normalization * amplitude;
        }
        let error = (ribbon.strain().unwrap().hamiltonian(&q) / integral - 1.0).abs();
        assert!(error < previous_error, "{axial_modes} axial modes: {error}");
        previous_error = error;
    }
    assert!(previous_error < 2e-13);
}

#[test]
fn fiber_strain_width_refinement_preserves_physical_energy_and_tension() {
    use fs_phs::Storage;
    let gate = CancelGate::new_clock_free();
    let mut previous: Option<(f64, f64, f64)> = None;
    for cells in [1, 3, 9] {
        let ribbon = fibers(cells, 2, 1300.0, 3e9, 0.0, 40.0, 48000, &gate)
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .unwrap();
        // Identical physical displacements for each quadrature strip.
        let mass_scale =
            (0.5 * ribbon.cell_line_density_kg_m() * ribbon.geometry().length_m).sqrt();
        let q = [mass_scale * 1e-8, 0.0, mass_scale * 2e-4, 0.0];
        let strain = ribbon.strain().unwrap();
        let observation = strain.observe(&q).unwrap();
        let current = (
            cells as f64 * strain.hamiltonian(&q),
            cells as f64 * observation.minimum_tension_n,
            cells as f64 * observation.maximum_tension_n,
        );
        if let Some(old) = previous {
            assert!((current.0 / old.0 - 1.0).abs() < 2e-14);
            assert!((current.1 / old.1 - 1.0).abs() < 2e-14);
            assert!((current.2 / old.2 - 1.0).abs() < 2e-14);
        }
        previous = Some(current);
    }
}

#[test]
fn normal_fiber_motion_drives_axial_waves_in_the_original_contact_bank() {
    // G3: the original linear rung cannot transfer normal motion into axial
    // waves. The conservative strain rung must, while retaining the string's
    // own nonlinear energy, all body/cavity states and the same work audit.
    with_cx(|cx, gate| {
        let make = |coupled| {
            let ribbon = fibers(3, 2, 1300.0, 3e9, 50.0, 40.0, 48000, gate);
            let ribbon = if coupled {
                ribbon
                    .with_axial_transverse_strain(strain_limits(), strain_config())
                    .unwrap()
            } else {
                ribbon
            };
            base(cx, gate, 7800.0)
                .with_geometric_stretching(0.2, strain_config(), gate)
                .unwrap()
                .with_tensioned_fiber_ribbon(ribbon, 0.4, ribbon_limits(), gate)
                .unwrap()
        };
        let mut linear = make(false);
        let mut coupled = make(true);
        assert!(coupled.source().bank().geometric_stretching().is_some());
        assert!(coupled.source().bank().material_cavity().is_some());
        let source = coupled.source().fiber_ribbon().unwrap();
        let first = source.component_start();
        let normal = source
            .model()
            .displacement_map(0.5, FiberDirection::Normal)
            .unwrap();
        let mut external = vec![0.0; coupled.system().mode_count()];
        let start: usize = coupled.system().components()[..first]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        for cell in 0..3 {
            for (mode, shape) in normal.iter().enumerate() {
                external[start + 4 * cell + mode] = 0.02 * shape;
            }
        }
        let drives = vec![drive(0.0, 0.0); 3];
        for sample in 0..128 {
            if sample == 64 {
                external.fill(0.0);
            }
            for bank in [&mut linear, &mut coupled] {
                let frame = bank.step(&external, &drives).unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
        assert!(linear.system().components()[first..].iter().all(|c| {
            c.states()[..2]
                .iter()
                .all(|s| s.velocity_m_sqrt_kg_per_s == 0.0)
        }));
        let axial_energy: f64 = coupled.system().components()[first..]
            .iter()
            .map(|c| {
                c.states()[..2]
                    .iter()
                    .map(|s| 0.5 * s.velocity_m_sqrt_kg_per_s.powi(2))
                    .sum::<f64>()
            })
            .sum();
        assert!(
            axial_energy > 1e-18,
            "normal motion did not reach the original axial wave state"
        );
        let q: Vec<_> = coupled.system().components()[first]
            .states()
            .iter()
            .map(|s| s.displacement_m_sqrt_kg)
            .collect();
        let observation = coupled
            .source()
            .fiber_ribbon()
            .unwrap()
            .model()
            .strain()
            .unwrap()
            .observe(&q)
            .unwrap();
        assert!(observation.mean_dynamic_strain > 0.0);
        // The same changed hair state must reach actual traction and the
        // existing string/body pressure observer when contact is engaged.
        let drives = vec![drive(0.01, 0.01); 3];
        let mut pressure_difference = 0.0;
        let mut traction_difference = 0.0;
        for _ in 0..64 {
            let a = linear.step(&external, &drives).unwrap().clone();
            let b = coupled.step(&external, &drives).unwrap().clone();
            pressure_difference += (a.observer_pressure_pa - b.observer_pressure_pa).abs();
            traction_difference += (a.ports[0].traction_n - b.ports[0].traction_n).abs();
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
        }
        assert!(pressure_difference > 1e-18 && traction_difference > 1e-18);
        let before = states(&coupled);
        let accepted = coupled.system().last_frame().unwrap().clone();
        external[start] = f64::NAN;
        assert!(coupled.step(&external, &drives).is_err());
        assert_eq!(states(&coupled), before);
        assert_eq!(coupled.system().last_frame(), Some(&accepted));
    });
}

#[test]
fn fiber_strain_refuses_invalid_motion_and_mismatched_joint_solve() {
    let gate = CancelGate::new_clock_free();
    let ribbon = fibers(3, 2, 1300.0, 3e9, 0.0, 40.0, 48000, &gate)
        .with_axial_transverse_strain(strain_limits(), strain_config())
        .unwrap();
    let law = ribbon.strain().unwrap();
    assert!(law.observe(&[0.0; 3]).is_err());
    assert!(law.observe(&[f64::NAN, 0.0, 0.0, 0.0]).is_err());
    assert!(law.observe(&[0.0, 0.0, 1.0, 0.0]).is_err());
    let normalization = (0.5 * ribbon.cell_line_density_kg_m() * ribbon.geometry().length_m).sqrt();
    let ea = ribbon.axial_material().elastic().axial_stiffness_n() * 60.0;
    let q =
        2.0 * (40.0 / 3.0) / ea * ribbon.geometry().length_m / std::f64::consts::PI * normalization;
    assert!(
        law.observe(&[q, 0.0, 0.0, 0.0]).is_err(),
        "a fiber with local compressive tension must refuse"
    );
    assert!(
        ribbon
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .is_err()
    );
    with_cx(|cx, gate| {
        let mut different = strain_config();
        different.max_iterations -= 1;
        let ribbon = fibers(3, 2, 1300.0, 3e9, 0.0, 40.0, 48000, gate)
            .with_axial_transverse_strain(strain_limits(), different)
            .unwrap();
        assert!(
            base(cx, gate, 7800.0)
                .with_geometric_stretching(0.2, strain_config(), gate)
                .unwrap()
                .with_tensioned_fiber_ribbon(ribbon, 0.4, ribbon_limits(), gate)
                .is_err()
        );
    });
}
