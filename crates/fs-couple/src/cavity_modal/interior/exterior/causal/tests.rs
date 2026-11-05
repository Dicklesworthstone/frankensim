use super::*;
use fs_bem::helmholtz::aperture::ApertureRow;
use fs_bem::reciprocal_radiation::Budget;
use fs_feec::acoustics::helmholtz::Medium as FeMedium;
use fs_feec::fixtures::kuhn_cube;
use fs_material::gas::GasSpec;
use fs_vfit::passivity::CertificateClass;

#[test]
fn original_box_bake_prepares_and_advances_its_own_full_field_and_volume_source() {
    use fs_bem::helmholtz::aperture::Aperture;
    use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousError;
    use fs_feec::acoustics::helmholtz::heterogeneous::time::exact::ExactTimeLimits;

    // One original 3 cm gas cavity, with the entire x=0 face open to the same
    // native closed exterior. The 2x2x2 P1 mesh has 27 pressure coordinates;
    // none is removed. Its 48 exterior triangles resolve the declared 2.54 kHz
    // band at more than six panels per wavelength along the longest facet.
    // These choices are fixed physical/numerical declarations, not a search
    // over fixtures or a weakening of the constructor's production admission.
    const SIDE: f64 = 0.03;
    const VOLUME: f64 = SIDE * SIDE * SIDE;
    let (complex, points) = kuhn_cube(2);
    let points: Vec<_> = points
        .iter()
        .map(|p| [SIDE * p[0], SIDE * p[1], SIDE * p[2]])
        .collect();
    let gas = gas();
    let fe = HelmholtzFe::new(&complex, &points, 1, medium(gas)).unwrap();
    let surface = boundary_surface(&fe);
    let faces: Vec<_> = fe
        .boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&v| points[v as usize][0] == 0.0))
        .collect();
    let mode = BoundaryMode::uniform(&faces);
    // An independent original FE trace integral for checking boundary work
    // from the accepted full-field pressure potential, below.
    let original_load = fe.boundary_linear_load(&mode.traces).unwrap();
    let mut opts = options();
    opts.angular_step_rad_s = 125.0;
    opts.relative_tolerance = 1e-3; // The declared native operator tolerance.
    let native_limits = load_limits(&opts);
    let dt = 2.0 * core::f64::consts::PI / native_limits.band_rad_s[1];
    let gate = CancelGate::new();
    let owner = BoundCausalRadiation::new(
        &fe,
        &[mode],
        &[0],
        &surface,
        &gas,
        3,
        spatial(),
        &route(),
        opts,
        dt,
        native_limits,
        &gate,
    )
    .unwrap();
    assert_eq!(owner.field().ndof(), fe.ndof());
    assert_eq!(owner.field().ndof(), 27);
    assert_eq!(owner.exterior_ports(), &[0]);
    assert_eq!(owner.fit_rungs().len(), 3);
    assert!(
        owner
            .causality_reports()
            .iter()
            .all(|report| report.tested > 0 && report.worst_ratio <= 1.0)
    );
    assert_eq!(owner.source().rows.len(), 128);
    assert_eq!(owner.validation_source().rows.len(), 127);

    // A second original BEM call constructs its own uniform physical opening,
    // rather than reading a fitted model or the adapter's aperture mapping.
    let independent_aperture = Aperture {
        shape: surface
            .triangles()
            .unwrap()
            .iter()
            .map(|triangle| {
                if triangle.iter().all(|point| point[0] == 0.0) {
                    1.0
                } else {
                    0.0
                }
            })
            .collect(),
        reference_point: [0.0; 3],
        normal: [-1.0, 0.0, 0.0],
    };
    let row = &owner.source().rows[63];
    let independent = aperture_impedance(
        &surface,
        &[independent_aperture],
        &[row.omega],
        Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        },
        &route(),
        &gate,
    )
    .unwrap();
    assert!(relative_gap(&row.impedance, &independent.rows[0].impedance).unwrap() < 1e-11);

    // This source borrows the exact operator created by the public constructor.
    // Its normalized uniform volume injection has physical amplitude m^3/s.
    let source = owner
        .field()
        .volume_velocity_source(
            3,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    assert!((source.integrated_shape() - 1.0).abs() < 2e-13);
    let time_limits = CausalTimeLimits {
        interior: ExactTimeLimits::default(),
        maximum_exterior_ports: 1,
        interface_relative_tolerance: 1e-8,
    };
    let mut time = owner
        .prepare(&[], &[&source], &[], time_limits, &gate)
        .unwrap();
    assert_eq!(time.state().pressure_pa().len(), fe.ndof());
    assert!(time.driven_ports().is_empty());
    let CausalExteriorLoad::Rational(load) = time.exterior() else {
        panic!("source-bound zero-extracted-delay construction changed its native owner");
    };
    assert_eq!(load.source_model(), owner.load().source_model());
    assert_eq!(load.reports().1.class, CertificateClass::GridOnly);

    let mut heard = false;
    for flow in [1e-7, 1e-7, 0.0, 0.0] {
        let before = time.state().pressure_potential_pa_s().to_vec();
        let step = time.step(&[], &[flow], &mut || false).unwrap();
        let integrated_force = original_load
            .iter()
            .zip(&before)
            .zip(time.state().pressure_potential_pa_s())
            .map(|((&weight, &a), &b)| weight * (b - a))
            .sum::<f64>();
        let mean_force = integrated_force / dt;
        let force_scale = mean_force
            .abs()
            .max(step.exterior_mean_force_n[0].abs())
            .max(f64::MIN_POSITIVE);
        assert!((mean_force - step.exterior_mean_force_n[0]).abs() < 2e-8 * force_scale);
        let work = step.outward_velocity_m_s[0] * integrated_force;
        let scale = step
            .work
            .delta_fluid_energy_j
            .abs()
            .max(step.work.source_work_j.abs())
            .max(step.work.exterior_exchange_work_j.abs())
            .max(f64::MIN_POSITIVE);
        assert!((work - step.work.exterior_exchange_work_j).abs() < 3e-8 * scale);
        assert!(step.work.balance_residual_j().abs() < 3e-8 * scale);
        assert!(step.interface_relative_residual <= time_limits.interface_relative_tolerance);
        heard |= time.state().pressure_pa().iter().any(|p| p.abs() > 1e-9);
    }
    assert!(heard);
    assert_eq!(time.accepted_steps(), 4);

    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        owner.prepare(&[], &[&source], &[], time_limits, &cancelled),
        Err(CausalRadiationError::Geometry(RadiationError::Cancelled))
    ));
}

fn options() -> RadiationRealizationOptions {
    RadiationRealizationOptions {
        samples: 128,
        angular_step_rad_s: 16.0,
        fits: [2, 4, 6].map(|order| FitOptions {
            fit_e: false,
            iterations: 20,
            ..FitOptions::new(order)
        }),
        causality: CausalityOptions {
            subtractions: 4,
            tail_bound_ratio: 2.0,
            relative_data_tolerance: 1e-12,
        },
        relative_tolerance: 1e-5,
        maximum_sample_bytes: 4 << 20,
        maximum_relocation_scratch_bytes: 64 << 20,
        maximum_relocation_work: 5_000_000_000,
        maximum_causality_work: 5_000_000_000,
    }
}
fn load_limits(options: &RadiationRealizationOptions) -> MimoLoadLimits {
    MimoLoadLimits {
        band_rad_s: options.band_rad_s(),
        sampled_passivity_points: 129,
        max_states: 64,
        max_memory_bytes: 4 << 20,
    }
}
fn spatial() -> RadiationBudget {
    RadiationBudget {
        maximum_pressure_dofs: 128,
        maximum_panels: 128,
        maximum_trace_entries: 1024,
        maximum_projection_terms: 32768,
    }
}
fn route() -> ApertureRoute {
    ApertureRoute::Galerkin {
        budget: Budget {
            maximum_unknowns: 256,
            maximum_cache_bytes: 32 << 20,
            maximum_work: 200_000_000,
            angular_refinement: 0,
        },
    }
}
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn medium(gas: GasState) -> FeMedium {
    FeMedium {
        density_kg_m3: gas.density,
        sound_speed_m_s: gas.sound_speed,
    }
}
fn boundary_surface(fe: &HelmholtzFe<'_>) -> SpherePanels {
    let triangles = fe
        .boundary_faces()
        .into_iter()
        .map(|face| {
            let mut t = face.map(|v| fe.positions()[v as usize]);
            let (_, normal, _) = fe.boundary_geometry(face).unwrap();
            let a: [f64; 3] = std::array::from_fn(|i| t[1][i] - t[0][i]);
            let b: [f64; 3] = std::array::from_fn(|i| t[2][i] - t[0][i]);
            let cross = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            if cross.iter().zip(normal).map(|(x, y)| x * y).sum::<f64>() < 0.0 {
                t.swap(1, 2);
            }
            t
        })
        .collect();
    SpherePanels::from_triangles(triangles).unwrap()
}

#[test]
fn signed_zero_net_flow_modes_retain_native_bem_force_coordinates() {
    let (complex, positions) = kuhn_cube(1);
    let positions: Vec<_> = positions
        .iter()
        .map(|p| [0.06 * p[0], 0.04 * p[1], 0.04 * p[2]])
        .collect();
    let gas = gas();
    let fe = HelmholtzFe::new(&complex, &positions, 2, medium(gas)).unwrap();
    let surface = boundary_surface(&fe);
    let faces: Vec<_> = fe
        .boundary_faces()
        .into_iter()
        .filter(|f| f.iter().all(|&v| positions[v as usize][0] == 0.0))
        .collect();
    let modes = [
        BoundaryMode::uniform(&faces),
        BoundaryMode {
            traces: faces
                .iter()
                .map(|&face| {
                    (
                        face,
                        face.map(|v| 2.0 * positions[v as usize][1] / 0.04 - 1.0),
                    )
                })
                .collect(),
        },
    ];
    let gate = CancelGate::new();
    let (_, apertures) =
        bind_openings(&fe, &modes, &[1, 0], &surface, &gas, spatial(), &gate).unwrap();
    // Original P2 pressure field remains complete while the velocity shapes
    // remain the declared linear traces; signed mode is first in source order.
    assert!(apertures[0].shape.iter().any(|&x| x < 0.0));
    assert!(apertures[0].shape.iter().any(|&x| x > 0.0));
    let source = aperture_impedance(
        &surface,
        &apertures,
        &[core::f64::consts::TAU * 480.0],
        Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        },
        &route(),
        &gate,
    )
    .unwrap();
    assert!(source.volume_gains[0].abs() < 1e-16);
    assert!((source.volume_gains[1] - 0.04 * 0.04).abs() < 1e-16);
    let converted = conjugated_rows(&source).unwrap();
    assert_eq!(converted[0].len(), 4);
    for (&a, &b) in converted[0].iter().zip(&source.rows[0].impedance) {
        assert_eq!(a, b.conj());
    }
    // Compare independent coherent native work readings in the two phasor
    // conventions, without dividing by the signed mode's zero volume gain.
    let velocity = [C64::new(0.03, 0.01), C64::new(-0.02, 0.04)];
    let work = |matrix: &[C64], v: &[C64]| {
        (0..2)
            .map(|a| {
                (0..2)
                    .map(|b| (v[a].conj() * matrix[2 * a + b] * v[b]).re)
                    .sum::<f64>()
            })
            .sum::<f64>()
            * 0.5
    };
    let a = work(&converted[0], &velocity);
    let b = work(&source.rows[0].impedance, &velocity.map(|x| x.conj()));
    assert!((a - b).abs() < 1e-13 * a.abs().max(b.abs()));
}

#[test]
fn opening_binding_refuses_substituted_gas_and_malformed_original_trace() {
    let (complex, positions) = kuhn_cube(1);
    let gas = gas();
    let fe = HelmholtzFe::new(&complex, &positions, 1, medium(gas)).unwrap();
    let surface = boundary_surface(&fe);
    let face = fe.boundary_faces()[0];
    let modes = [BoundaryMode::uniform(&[face])];
    let gate = CancelGate::new();
    let mut other = gas;
    other.density *= 1.01;
    assert!(bind_openings(&fe, &modes, &[0], &surface, &other, spatial(), &gate).is_err());
    let duplicate = [BoundaryMode::uniform(&[face, face])];
    assert!(bind_openings(&fe, &duplicate, &[0], &surface, &gas, spatial(), &gate).is_err());
    let empty = [BoundaryMode { traces: vec![] }];
    assert!(bind_openings(&fe, &empty, &[0], &surface, &gas, spatial(), &gate).is_err());
    let bad = [BoundaryMode {
        traces: vec![(face, [1.0, f64::NAN, 1.0])],
    }];
    assert!(bind_openings(&fe, &bad, &[0], &surface, &gas, spatial(), &gate).is_err());
}

#[test]
fn sampled_band_clock_and_native_work_limits_refuse_before_baking() {
    let opts = options();
    let limits = load_limits(&opts);
    let dt = 2.0 * core::f64::consts::PI / limits.band_rad_s[1];
    let (omega, validation) = admit(&opts, 2, dt, limits).unwrap();
    assert_eq!(omega[0], limits.band_rad_s[0]);
    assert_eq!(*omega.last().unwrap(), limits.band_rad_s[1]);
    assert_eq!(validation.len() + 1, omega.len());
    for (pair, &w) in omega.windows(2).zip(&validation) {
        assert!(pair[0] < w && w < pair[1]);
        assert!(!omega.contains(&w));
    }
    let mut wider = limits;
    wider.band_rad_s[1] *= 1.001;
    assert!(admit(&opts, 2, dt, wider).is_err());
    assert!(admit(&opts, 2, dt * 0.49, limits).is_err());
    let mut excessive = opts.clone();
    excessive.fits[2].fit_e = true;
    assert!(admit(&excessive, 2, dt, limits).is_err());
    excessive = opts.clone();
    excessive.maximum_causality_work = 1;
    assert!(admit(&excessive, 2, dt, limits).is_err());
    excessive = opts.clone();
    excessive.maximum_relocation_scratch_bytes = 1;
    assert!(admit(&excessive, 2, dt, limits).is_err());
    excessive = opts.clone();
    excessive.maximum_sample_bytes = 1;
    assert!(admit(&excessive, 2, dt, limits).is_err());
}

/// Independent reciprocal positive-real matrix: six physical relaxation
/// impedances k*s/(s+a), each acting through a real spatial work vector.
/// It contains all self and mutual terms and has no improper part.
fn relaxation_matrix(omega: f64) -> Vec<C64> {
    let s = C64::new(0.0, omega);
    let mut z = vec![C64::from_re(0.5), C64::ZERO, C64::ZERO, C64::from_re(0.8)];
    for (k, pole) in [30.0, 100.0, 300.0, 900.0, 2700.0, 8100.0]
        .into_iter()
        .enumerate()
    {
        let direction = [1.0, if k % 2 == 0 { 0.6 } else { -0.35 }];
        let relaxation = s * (s + C64::from_re(pole)).recip();
        for a in 0..2 {
            for b in 0..2 {
                z[a * 2 + b] = z[a * 2 + b]
                    + relaxation.scale((0.4 + 0.1 * k as f64) * direction[a] * direction[b]);
            }
        }
    }
    z
}

#[test]
fn original_native_order_ladder_and_independent_validation_feed_exact_zoh() {
    let opts = options();
    let limits = load_limits(&opts);
    let dt = 2.0 * core::f64::consts::PI / limits.band_rad_s[1];
    let (omega, validation_omega) = admit(&opts, 2, dt, limits).unwrap();
    let data: Vec<_> = omega.iter().map(|&w| relaxation_matrix(w)).collect();
    let validation: Vec<_> = validation_omega
        .iter()
        .map(|&w| relaxation_matrix(w))
        .collect();
    let gate = CancelGate::new();
    let (model, rungs) = fit_ladder(
        &omega,
        &data,
        &validation_omega,
        &validation,
        2,
        &opts,
        &gate,
    )
    .unwrap();
    assert_eq!(model.ports, 2);
    assert!(model.e.iter().all(|&x| x == 0.0));
    assert_eq!(rungs.len(), 3);
    assert!(rungs[0].previous_validation_change.is_none());
    assert!(
        rungs[1..]
            .iter()
            .all(|r| r.previous_validation_change.is_some())
    );
    assert!(rungs[2].fit.max_relative_error <= opts.relative_tolerance);
    assert!(rungs[2].validation_max_relative_error <= opts.relative_tolerance);
    let mut load = PreparedMimoLoad::new(&model, dt, limits).unwrap();
    assert!(load.reports().0.passive);
    assert_eq!(load.reports().1.class, CertificateClass::GridOnly);
    let mut force = [0.0; 2];
    for velocity in [[0.01, -0.02], [0.0, 0.0], [-0.01, 0.03]] {
        load.step_mean_into(&velocity, &mut force).unwrap();
        assert!(force.iter().all(|x| x.is_finite()));
    }
    // A good training fit cannot hide a different physical operator on heldout
    // frequencies: those samples never become extra least-squares fit rows.
    let wrong: Vec<Vec<_>> = validation
        .iter()
        .map(|row| row.iter().map(|z| z.scale(2.0)).collect())
        .collect();
    let result = fit_ladder(&omega, &data, &validation_omega, &wrong, 2, &opts, &gate);
    assert!(matches!(result, Err(CausalRadiationError::Accuracy { .. })));
}

#[test]
fn outgoing_bem_convention_is_conjugated_once_before_native_causality() {
    // Analytic pulsating-sphere radiation in N s/m. This tests only convention
    // and causality; a sphere is not misrepresented as a planar FE opening.
    let count = 1025;
    let step = core::f64::consts::TAU * 5000.0 / 1024.0;
    let omega: Vec<_> = (0..count).map(|j| (j as f64 + 0.5) * step).collect();
    let rows: Vec<_> = omega
        .iter()
        .map(|&w| {
            let s = C64::new(0.0, w * 0.05 / 343.0);
            let z = (s * (C64::ONE + s).recip()).scale(413.0);
            ApertureRow {
                omega: w,
                impedance: vec![z.conj()],
                route: "analytic-sphere",
                reciprocity_defect: 0.0,
                passivity_margin: z.re / z.abs(),
            }
        })
        .collect();
    let source = ApertureImpedance {
        declaration: 0,
        apertures: 1,
        volume_gains: vec![0.0],
        rows,
    };
    let data = conjugated_rows(&source).unwrap();
    let opts = options().causality;
    assert!(check_mimo_causality(&omega, &data, 1, &[0.0], &opts).is_ok());
    let wrong: Vec<Vec<_>> = data
        .iter()
        .map(|r| r.iter().map(|z| z.conj()).collect())
        .collect();
    assert!(matches!(
        check_mimo_causality(&omega, &wrong, 1, &[0.0], &opts),
        Err(CausalityError::NonCausal { .. })
    ));
}
