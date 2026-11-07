//! G1/G3: nonlinear composite slabs, wall heat observers and conservative storage.

use super::*;
use fs_lbm::conjugate::FaceContactResistance;

// Independent primitives for two positive, linear conductivity laws. These
// formulas are not obtained from the production table integration or root solve.
fn potential_a(t: f64) -> f64 {
    let x = t - 300.0;
    x + 0.005 * x * x
}

fn temperature_a(potential: f64) -> f64 {
    300.0 + 2.0 * potential / (1.0 + (1.0 + 0.02 * potential).sqrt())
}

fn potential_b(t: f64) -> f64 {
    let x = t - 300.0;
    4.0 * x - 0.005 * x * x
}

fn temperature_b(potential: f64) -> f64 {
    300.0 + 2.0 * potential / (4.0 + (16.0 - 0.02 * potential).sqrt())
}

fn materials() -> [SolidMaterial; 2] {
    [
        SolidMaterial::new("increasing conductivity", 1.0).with_conductivity_table(&[
            (250.0, 0.5),
            (330.0, 1.3),
            (350.0, 1.5),
            (370.0, 1.7),
            (450.0, 2.5),
        ]),
        SolidMaterial::new("decreasing conductivity", 4.0).with_conductivity_table(&[
            (250.0, 4.5),
            (330.0, 3.7),
            (350.0, 3.5),
            (370.0, 3.3),
            (450.0, 2.5),
        ]),
    ]
}

#[test]
fn nonlinear_material_interface_matches_composite_slab_in_each_direction() {
    let flux = 100.0;
    let (width_a, width_b) = (0.5, 0.8);
    let gate = CancelGate::new();
    for resistance in [0.0, 0.2] {
        // Manufacture a known flux and contact-temperature jump, then derive
        // the boundary and cell temperatures from independent quadratics.
        let surface_b = 340.0;
        let surface_a = surface_b + flux * resistance;
        let hot_wall = temperature_a(potential_a(surface_a) + flux * width_a);
        let cold_wall = temperature_b(potential_b(surface_b) - flux * width_b);
        let expected_a = temperature_a(potential_a(surface_a) + 0.5 * flux * width_a);
        let expected_b = temperature_b(potential_b(surface_b) - 0.5 * flux * width_b);
        for axis in 0..3 {
            for reversed in [false, true] {
                let mut widths = [vec![0.2], vec![0.3], vec![0.4]];
                widths[axis] = if reversed {
                    vec![width_b, width_a]
                } else {
                    vec![width_a, width_b]
                };
                let split = widths[axis][0];
                let domain = VoxelDomain::graded_from_fn(widths, |p| {
                    let a = (p[axis] < split) != reversed;
                    Voxel::Solid(if a { 0 } else { 1 })
                })
                .unwrap();
                let mut faces = [ThermalFace::Adiabatic; 6];
                faces[2 * axis] =
                    ThermalFace::Temperature(if reversed { cold_wall } else { hot_wall });
                faces[2 * axis + 1] =
                    ThermalFace::Temperature(if reversed { hot_wall } else { cold_wall });
                let mut setup = ThermalSetup::new(faces);
                setup.contacts.push(ContactResistance {
                    materials: (1, 0),
                    resistance_m2_k_w: resistance,
                });
                let solution = solve_energy(
                    &domain,
                    &unit_fluid(),
                    &materials(),
                    &FlowField::quiescent(&domain),
                    &setup,
                    &EnergyConfig::default(),
                    &gate,
                )
                .unwrap();
                let a = usize::from(reversed);
                let b = 1 - a;
                assert!(
                    (solution.temperature[a] - expected_a).abs() < 1e-7,
                    "A: axis={axis}, reversed={reversed}, R={resistance}, {solution:?}"
                );
                assert!(
                    (solution.temperature[b] - expected_b).abs() < 1e-7,
                    "B: axis={axis}, reversed={reversed}, R={resistance}, {solution:?}"
                );
                let through_a = (potential_a(hot_wall) - potential_a(solution.temperature[a]))
                    / (0.5 * width_a);
                let through_b = (potential_b(solution.temperature[b]) - potential_b(cold_wall))
                    / (0.5 * width_b);
                assert!((through_a - flux).abs() < 1e-6);
                assert!((through_b - flux).abs() < 1e-6);
                assert!(solution.report.balance.relative_residual < 1e-9);
                // The report must be the original law at these exact cells,
                // not the final linear solve's previous conductivity epoch.
                let ka = 1.0 + 0.01 * (solution.temperature[a] - 300.0);
                let kb = 4.0 - 0.01 * (solution.temperature[b] - 300.0);
                for actual in solution.conductivity[a] {
                    assert!((actual - ka).abs() < 2e-15 * ka);
                }
                for actual in solution.conductivity[b] {
                    assert!((actual - kb).abs() < 2e-15 * kb);
                }
                let boundary_heat = domain.face_area(a, axis) * (through_b - through_a);
                assert!((solution.report.balance.boundary_outflow_w - boundary_heat).abs() < 1e-10);
                assert!(solution.report.relative_residual <= EnergyConfig::default().tolerance);
                assert!(solution.report.iterations > 0);
            }
        }
    }
}

#[test]
fn nonlinear_endpoint_refuses_cancelling_local_heat_errors() {
    // Both reservoirs and the initial iterate are 350 K, but the cells have
    // opposite 100 W sources. That uniform field closes the global heat ledger
    // EXACTLY while violating both local equations by 100 W. Its row-scaled
    // algebraic residual is below 0.2 because that norm contains the 350 K
    // reservoir offset. No Krylov work is allowed to repair the field.
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let material = SolidMaterial::new("nonlinear slab", 1.0)
        .with_conductivity_table(&[(300.0, 1.0), (400.0, 2.0)])
        .with_heat_capacity(1.0);
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(350.0);
    faces[1] = ThermalFace::Temperature(350.0);
    let mut setup = ThermalSetup::new(faces);
    setup.power_w = vec![100.0, -100.0];
    let config = EnergyConfig { tolerance: 0.2, max_iterations: 0,
        ..EnergyConfig::default() };
    let error = solve_energy(
        &domain, &unit_fluid(), std::slice::from_ref(&material),
        &FlowField::quiescent(&domain), &setup, &config, &CancelGate::new(),
    ).unwrap_err();
    assert!(matches!(error, ChtError::SolverNotConverged { .. }), "{error:?}");

    // With funded Krylov work the same loose request can be repaired. The
    // acceptance is physical: independently integrate both boundary halves
    // and the common full-cell segment at the returned temperatures.
    let funded = EnergyConfig { max_iterations: 50, ..config };
    let accepted = solve_energy(
        &domain, &unit_fluid(), std::slice::from_ref(&material),
        &FlowField::quiescent(&domain), &setup, &funded, &CancelGate::new(),
    ).unwrap();
    let [a, b] = [accepted.temperature[0], accepted.temperature[1]];
    assert!(a > 350.0 && b < 350.0);
    let boundary = [2.0 * (potential_a(a) - potential_a(350.0)),
        2.0 * (potential_a(b) - potential_a(350.0))];
    let internal = potential_a(a) - potential_a(b);
    for (i, transfer) in [internal, -internal].into_iter().enumerate() {
        let residual = boundary[i] + transfer - setup.power_w[i];
        let scale = boundary[i].abs() + transfer.abs() + setup.power_w[i].abs();
        assert!(residual.abs() <= config.tolerance * scale);
    }
    assert!((accepted.report.balance.boundary_outflow_w - boundary.iter().sum::<f64>()).abs() < 1e-10);
    assert!(accepted.report.iterations > 0);

    // The exact same defect must not become a physical time step: frozen old
    // storage is zero at this candidate, as are every boundary/contact flux.
    let error = march_energy(
        &domain, &unit_fluid(), &[material], &FlowField::quiescent(&domain),
        &setup, &[350.0; 2], |_| 1.0,
        &TransientConfig { time_step_s: 0.1, steps: 1, energy: config },
        None, &CancelGate::new(),
    ).unwrap_err();
    assert!(matches!(error, ChtError::SolverNotConverged { .. }), "{error:?}");
}

#[test]
fn nonlinear_solid_fluid_heat_observer_matches_actual_contact_flux_and_sources() {
    let domain = VoxelDomain::graded_from_fn([vec![0.5, 0.8], vec![0.2], vec![0.3]], |p| {
        if p[0] < 0.5 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let area = 0.2 * 0.3;
    let fluid = FluidProperties {
        conductivity_w_m_k: 1.25,
        ..unit_fluid()
    };
    for sign in [-1.0, 1.0] {
        let flux = sign * 100.0;
        let surface_solid = 350.0 + sign * 10.0;
        let surface_fluid = 350.0 - sign * 10.0;
        let solid_cell = temperature_a(potential_a(surface_solid) + flux * 0.25);
        let fluid_cell = surface_fluid - flux * 0.4 / 2.0;
        let source_w = 3.0;
        let outward_w = area * flux + source_w;
        let outer_temperature = fluid_cell - outward_w / area * 0.4 / 2.0;
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[1] = ThermalFace::Temperature(outer_temperature);
        let mut setup = ThermalSetup::new(faces);
        setup.fixed_temperature = vec![(0, solid_cell)];
        setup.power_w = vec![0.0, source_w];
        // The fluid's admitted molecular and eddy values sum to 2 W/(m K).
        setup.eddy_conductivity_w_m_k = vec![0.0, 0.75];
        setup.face_contacts.push(FaceContactResistance {
            cells: [1, 0],
            resistance_m2_k_w: 0.2,
        });
        let solution = solve_energy(
            &domain,
            &fluid,
            &materials()[..1],
            &FlowField::quiescent(&domain),
            &setup,
            &EnergyConfig::default(),
            &CancelGate::new(),
        )
        .unwrap();
        assert!((solution.temperature[1] - fluid_cell).abs() < 1e-7);
        assert!((solution.solid_to_fluid_heat_w(&domain) - area * flux).abs() < 1e-7);
        let balance = solution.report.balance;
        assert!((balance.fixed_cell_injection_w - area * flux).abs() < 1e-7);
        assert!((balance.source_w - source_w).abs() < 1e-12);
        assert!((balance.boundary_outflow_w - outward_w).abs() < 1e-7);
        assert!(balance.residual_w.abs() < 1e-7);
    }
}

#[test]
fn nonlinear_material_contact_transient_matches_independent_endpoint_and_storage() {
    let domain = VoxelDomain::graded_from_fn([vec![0.5, 0.8], vec![0.2], vec![0.3]], |p| {
        Voxel::Solid(u16::from(p[0] >= 0.5))
    })
    .unwrap();
    let [a, b] = materials();
    let solids = [a.with_heat_capacity(100.0), b.with_heat_capacity(250.0)];
    let area = 0.2 * 0.3;
    let capacity = [100.0 * 0.5 * area, 250.0 * 0.8 * area];
    let initial = [390.0, 310.0];
    let dt = 1.0;
    let resistance = 0.2;
    // One scalar endpoint oracle eliminates both stored energies first, then
    // uses analytical conductivity primitives to recover the two surface
    // temperatures. It does not call the production interface primitive.
    let cells_at = |flux: f64| {
        [
            initial[0] - dt * area * flux / capacity[0],
            initial[1] + dt * area * flux / capacity[1],
        ]
    };
    let mismatch = |flux: f64| {
        let [ta, tb] = cells_at(flux);
        let sa = temperature_a(potential_a(ta) - 0.25 * flux);
        let sb = temperature_b(potential_b(tb) + 0.4 * flux);
        sa - sb - resistance * flux
    };
    let (mut lo, mut hi) = (0.0, (initial[0] - initial[1]) / resistance);
    assert!(mismatch(lo) > 0.0 && mismatch(hi) < 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if mismatch(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let flux = 0.5 * (lo + hi);
    let expected = cells_at(flux);
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: resistance,
    });
    let result = march_energy(
        &domain,
        &unit_fluid(),
        &solids,
        &FlowField::quiescent(&domain),
        &setup,
        &initial,
        |_| 1.0,
        &TransientConfig {
            time_step_s: dt,
            steps: 1,
            energy: EnergyConfig::default(),
        },
        None,
        &CancelGate::new(),
    )
    .unwrap();
    for (actual, expected) in result.temperature.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-7, "{result:?}");
    }
    let independently_stored: f64 = (0..2)
        .map(|c| capacity[c] * (result.temperature[c] - initial[c]))
        .sum();
    assert!(independently_stored.abs() < 1e-7);
    let record = result.records[0];
    assert!(record.stored_energy_change_j.abs() < 1e-7);
    assert_eq!(record.source_j, 0.0);
    assert_eq!(record.boundary_outflow_j, 0.0);
    assert!(record.closure_j.abs() < 1e-7);
}


// A sharp but strictly positive phase-region conductivity. The independent
// integral below is measured from the 300 K reservoir; no production table
// integration or interface solver is used by these endpoint oracles.
fn sharp_material() -> SolidMaterial {
    SolidMaterial::new("sharp conductivity transition", 1.0)
        .with_conductivity_table(&[(300.0, 1.0), (350.0, 1.0), (351.0, 100.0), (400.0, 100.0)])
        .with_heat_capacity(2.0)
}

fn sharp_potential(t: f64) -> f64 {
    let ramp = (t - 350.0).clamp(0.0, 1.0);
    t.min(350.0) - 300.0 + ramp + 49.5 * ramp * ramp
        + 100.0 * (t - 351.0).max(0.0)
}

#[test]
fn sharp_conductivity_heating_closes_steady_and_original_time_storage() {
    // Minmod previously cycled inside its nested reconstruction solve before
    // the conductivity owner's line search could inspect a single proposal.
    for scheme in [ConvectionScheme::PowerLaw, ConvectionScheme::EnthalpyMinmod] {
        let config = EnergyConfig { scheme, ..EnergyConfig::default() };
        let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let material = sharp_material();
        let flow = FlowField::quiescent(&domain);
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::Temperature(300.0);
        let mut setup = ThermalSetup::new(faces);
        let expected = 350.5;
        let outward = 2.0 * sharp_potential(expected); // half-cell Fourier integral
        setup.power_w = vec![outward];
        let solution = solve_energy(
            &domain, &unit_fluid(), std::slice::from_ref(&material), &flow,
            &setup, &config, &CancelGate::new(),
        ).unwrap();
        // Undamped Picard instead cycles between 362.875 and 303.069305609472 K.
        assert!((solution.temperature[0] - expected).abs() < 1e-7, "{solution:?}");
        let actual_outward = 2.0 * sharp_potential(solution.temperature[0]);
        assert!((actual_outward - outward).abs() < 1e-7);
        assert!((solution.report.balance.boundary_outflow_w - actual_outward).abs() < 1e-8);
        assert!(solution.report.balance.relative_residual < 1e-9);
        assert!(solution.report.relative_residual <= config.tolerance);
        for k in solution.conductivity[0] {
            assert!((k - 50.5).abs() < 1e-7);
        }

        // One backward-Euler endpoint has the same target with additional stored
        // energy. Every line-search trial must retain this ORIGINAL old state.
        let initial = [300.0];
        let dt = 0.25;
        let stored = 2.0 * (expected - initial[0]);
        setup.power_w[0] = outward + stored / dt;
        let transient = march_energy(
            &domain, &unit_fluid(), &[material], &flow, &setup, &initial, |_| 1.0,
            &TransientConfig { time_step_s: dt, steps: 1, energy: config },
            None, &CancelGate::new(),
        ).unwrap();
        assert!((transient.temperature[0] - expected).abs() < 1e-7, "{transient:?}");
        let record = transient.records[0];
        assert!((record.stored_energy_change_j - stored).abs() < 1e-7);
        assert!((record.boundary_outflow_j - dt * outward).abs() < 1e-7);
        assert!((record.source_j - (stored + dt * outward)).abs() < 1e-10);
        assert!(record.closure_j.abs() < 1e-7);
    }
}

#[test]
fn sharp_conductivity_graded_composite_preserves_contact_jump_in_each_axis() {
    let (width_a, width_b) = (0.5, 0.8);
    let (cell_a, surface_a, resistance, kb) = (350.5, 350.2, 0.2, 2.0);
    let flux = (sharp_potential(cell_a) - sharp_potential(surface_a)) / (0.5 * width_a);
    let surface_b = surface_a - resistance * flux;
    let cell_b = surface_b - flux * 0.5 * width_b / kb;
    let bath = surface_b - flux * width_b / kb;
    let solids = [sharp_material(), SolidMaterial::new("constant backing", kb)];
    for axis in 0..3 {
        for reversed in [false, true] {
            let mut widths = [vec![0.2], vec![0.3], vec![0.4]];
            widths[axis] = if reversed { vec![width_b, width_a] } else { vec![width_a, width_b] };
            let split = widths[axis][0];
            let domain = VoxelDomain::graded_from_fn(widths, |p| {
                Voxel::Solid(u16::from((p[axis] < split) == reversed))
            }).unwrap();
            let (a, b) = (usize::from(reversed), usize::from(!reversed));
            let area = domain.face_area(a, axis);
            let mut faces = [ThermalFace::Adiabatic; 6];
            faces[2 * axis + usize::from(!reversed)] = ThermalFace::Temperature(bath);
            let mut setup = ThermalSetup::new(faces);
            setup.power_w = vec![0.0; 2];
            setup.power_w[a] = flux * area;
            setup.contacts.push(ContactResistance {
                materials: (0, 1), resistance_m2_k_w: resistance,
            });
            let solution = solve_energy(
                &domain, &unit_fluid(), &solids, &FlowField::quiescent(&domain),
                &setup, &EnergyConfig::default(), &CancelGate::new(),
            ).unwrap();
            assert!((solution.temperature[a] - cell_a).abs() < 1e-7, "{solution:?}");
            assert!((solution.temperature[b] - cell_b).abs() < 1e-7, "{solution:?}");
            let actual_flux = kb * (solution.temperature[b] - bath) / (0.5 * width_b);
            assert!((actual_flux - flux).abs() < 1e-7);
            // Invert the independent quadratic in the crossed transition
            // interval to recover A's physical surface, then check the jump.
            let potential = sharp_potential(solution.temperature[a])
                - actual_flux * 0.5 * width_a - 50.0;
            let sa = 350.0 + 2.0 * potential / (1.0 + (1.0 + 198.0 * potential).sqrt());
            let sb = solution.temperature[b] + actual_flux * 0.5 * width_b / kb;
            assert!((sa - sb - resistance * actual_flux).abs() < 1e-7);
            assert!((solution.report.balance.boundary_outflow_w - area * flux).abs() < 1e-7);
            assert!(solution.report.balance.relative_residual < 1e-9);
        }
    }
}

#[test]
fn heterogeneous_graded_face_peclet_matches_the_actual_series_resistance() {
    let fluid = FluidProperties {
        density_kg_m3: 2.0,
        specific_heat_j_kg_k: 3.0,
        ..unit_fluid()
    };
    let field = [[0.5, 1.0, 2.0], [5.0, 10.0, 20.0]];
    for axis in 0..3 {
        let mut widths = [vec![0.3], vec![0.4], vec![0.5]];
        widths[axis] = vec![0.2, 0.8];
        let domain = VoxelDomain::graded(widths).unwrap();
        for speed in [-2.0_f64, 2.0] {
            let mut flow_faces = [FlowFace::Wall; 6];
            flow_faces[2 * axis] = FlowFace::Fixed;
            flow_faces[2 * axis + 1] = FlowFace::Free;
            let flow = FlowField::from_face_velocity(&domain, flow_faces, |component, _| {
                if component == axis { speed } else { 0.0 }
            }).unwrap();
            let inlet = 2 * axis + usize::from(speed < 0.0);
            let mut faces = [ThermalFace::Adiabatic; 6];
            faces[inlet] = ThermalFace::Inflow { temperature: 300.0 };
            faces[inlet ^ 1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
            let solution = EnergyConfig::default().solve_with_conductivity(
                &domain, &fluid, &[], &flow, &ThermalSetup::new(faces), &field, &CancelGate::new(),
            ).unwrap();
            // Area cancels between rho*Cp*Q and A/(h0/(2*k0)+h1/(2*k1)).
            let expected = fluid.density_kg_m3 * fluid.specific_heat_j_kg_k * speed.abs()
                * (0.1 / field[0][axis] + 0.4 / field[1][axis]);
            assert!((solution.report.max_cell_peclet - expected).abs() < 1e-12 * expected,
                "axis={axis}, speed={speed}: {:?}", solution.report);
            for t in solution.temperature {
                assert!((t - 300.0).abs() < 1e-8);
            }
        }
    }
}
