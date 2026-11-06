use fs_lbm::ThermalLbm;
use fs_lbm::d3q19::ThermalLbm3;
use fs_lbm::thermal::sources::{ThermalSourceError, ThermalSourceLimits};

fn limits() -> ThermalSourceLimits {
    ThermalSourceLimits { max_cells: 4096, temperature_range: [-2.0, 2.0], balance_tolerance: 1e-10 }
}
fn heat(cells: usize, plane: usize, q: f64) -> Vec<f64> {
    let mut h = vec![0.0; cells]; h[plane..cells-plane].fill(q); h
}
fn mass<const Q: usize>(f: &[[f64; Q]]) -> f64 { f.iter().flatten().sum() }

#[test]
fn g1_switching_held_forces_delivers_exactly_one_impulse_per_step() {
    let density = 1.7;
    let mut two = ThermalLbm::slab(3, 32, 0.9, 0.8, 0.0);
    let mut three = ThermalLbm3::slab(3, 2, 32, 0.9, 0.9, 0.0);
    let mut velocity2 = [0.01, -0.02];
    let mut velocity3 = [0.01, -0.02, 0.015];
    two.grid.f.fill(fs_lbm::equilibrium(density, velocity2[0], velocity2[1]));
    three.grid.f.fill(fs_lbm::equilibrium3(density, velocity3));
    // Central cells remain outside the walls' finite lattice domain of
    // dependence. Uniform force has no pressure or viscous momentum exchange:
    // rho * delta u = F * delta t, with delta t=1 in lattice units.
    for amplitude in [1e-4, 1e-4, 0.0, -2e-4, 1.5e-4] {
        let force2 = [amplitude, -0.5 * amplitude];
        let force3 = [amplitude, -0.5 * amplitude, 0.25 * amplitude];
        let force2 = vec![force2; two.gpop.len()];
        let force3 = vec![force3; three.gpop.len()];
        two.step_with_sources(&vec![0.0; force2.len()], &force2, limits()).unwrap();
        three.step_with_external_forces(&force3).unwrap();
        let actual2 = two.grid.moments(two.grid.idx(1, 16));
        let actual3 = three.grid.moments3(three.grid.idx(1, 1, 16));
        assert!((actual2.rho - density).abs() < 2e-14);
        assert!((actual3.0 - density).abs() < 2e-14);
        for a in 0..2 {
            velocity2[a] += force2[0][a] / density;
            assert!((actual2.u[a] - velocity2[a]).abs() < 2e-14);
        }
        for a in 0..3 {
            velocity3[a] += force3[0][a] / density;
            assert!((actual3.1[a] - velocity3[a]).abs() < 2e-14);
        }
    }
    // Returning to the ordinary zero-buoyancy step must also remove the last
    // force's coordinate frame without a physical kick.
    two.step();
    three.step();
    for (actual, expected) in two.grid.moments(two.grid.idx(1, 16)).u.iter().zip(velocity2) {
        assert!((actual - expected).abs() < 2e-14);
    }
    for (actual, expected) in three.grid.moments3(three.grid.idx(1, 1, 16)).1.iter().zip(velocity3) {
        assert!((actual - expected).abs() < 2e-14);
    }
}

#[test]
fn zero_heat_preserves_the_existing_forced_step_in_both_dimensions() {
    let mut a = ThermalLbm::slab(3, 8, 0.9, 0.8, 1e-5);
    let mut b = ThermalLbm::slab(3, 8, 0.9, 0.8, 1e-5);
    b.set_coriolis(fs_lbm::core2::frame::Coriolis2::default()).unwrap();
    let f = vec![[1e-6, -2e-6]; a.gpop.len()];
    for _ in 0..12 {
        a.step_with_external_forces(&f).unwrap();
        let balance = b.step_with_sources(&vec![0.0; f.len()], &f, limits()).unwrap();
        assert_eq!(balance.source, 0.0);
        assert_eq!(a.grid.f, b.grid.f); assert_eq!(a.gpop, b.gpop); assert_eq!(a.grid.fext, b.grid.fext);
    }
    let mut a = ThermalLbm3::slab(3, 2, 8, 0.9, 0.9, 1e-5);
    let mut b = ThermalLbm3::slab(3, 2, 8, 0.9, 0.9, 1e-5);
    b.set_coriolis(fs_lbm::d3q19::Coriolis3::default()).unwrap();
    let f = vec![[1e-6, -2e-6, 3e-6]; a.gpop.len()];
    for _ in 0..12 {
        a.step_with_external_forces(&f).unwrap();
        b.step_with_sources(&vec![0.0; f.len()], &f, limits()).unwrap();
        assert_eq!(a.grid.f, b.grid.f); assert_eq!(a.gpop, b.gpop); assert_eq!(a.grid.fext, b.grid.fext);
    }
}

#[test]
fn trt_source_steps_retain_native_flow_and_thermal_history() {
    use fs_lbm::core2::CollisionModel2;

    let mut two = ThermalLbm::slab(3, 6, 0.9, 0.8, 1e-5);
    let mut direct2 = ThermalLbm::slab(3, 6, 0.9, 0.8, 1e-5);
    let mut three = ThermalLbm3::slab(3, 1, 6, 0.9, 0.8, 1e-5);
    let mut direct3 = ThermalLbm3::slab(3, 1, 6, 0.9, 0.8, 1e-5);
    for slab in [&mut two, &mut direct2] {
        slab.grid.collision = CollisionModel2::Trt { magic: 3.0 / 16.0 };
        slab.gpop[10][0] += 0.01;
    }
    for slab in [&mut three, &mut direct3] {
        slab.grid.set_trt_magic(3.0 / 16.0).unwrap();
        slab.gpop[10][0] += 0.01;
    }
    for step in 0..12 {
        let load = if step < 5 { 2e-5 } else { -1e-5 };
        let force2 = [[load, 1e-5]; 24];
        let force3 = [[load, -0.5 * load, 1e-5]; 24];
        direct2.step_with_external_forces(&force2).unwrap();
        direct3.step_with_external_forces(&force3).unwrap();
        two.step_with_sources(&[0.0; 24], &force2, limits()).unwrap();
        three.step_with_sources(&[0.0; 24], &force3, limits()).unwrap();
        assert_eq!(two.grid.f, direct2.grid.f);
        assert_eq!(two.gpop, direct2.gpop);
        assert_eq!(three.grid.f, direct3.grid.f);
        assert_eq!(three.gpop, direct3.gpop);
        assert_eq!(three.grid.trt_magic(), Some(3.0 / 16.0));
    }
}

#[test]
fn trt_material_heat_and_force_use_the_original_physical_scaling() {
    use fs_lbm::material::{MaterialThermalSlab3, plan_scaling_from_material};
    use fs_lbm::thermal::sources::HeatSourceLimits;
    use fs_material::fluid_transport::FluidTransportState;

    let material = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.001).unwrap();
    let input = MaterialThermalSlab3 {
        nx: 2, ny: 1, nz: 6,
        bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    };
    let mut actual = scale.thermal_slab3(input).unwrap();
    let mut direct = scale.thermal_slab3(input).unwrap();
    actual.lattice.grid.set_trt_magic(3.0 / 16.0).unwrap();
    direct.lattice.grid.set_trt_magic(3.0 / 16.0).unwrap();
    let native_force = [scale.force_density_lattice(20.0).unwrap(), 0.0, 0.0];
    for _ in 0..12 {
        actual.step_with_power_density(
            &[0.0; 16], &[[20.0, 0.0, 0.0]; 16],
            HeatSourceLimits {
                max_cells: 16, temperature_range_k: [290.0, 320.0],
                balance_tolerance_j: 1e-10,
            },
        ).unwrap();
        direct.lattice.step_with_external_forces(&[native_force; 16]).unwrap();
        assert_eq!(actual.lattice.grid.f, direct.lattice.grid.f);
        assert_eq!(actual.lattice.gpop, direct.lattice.gpop);
        assert_eq!(actual.lattice.grid.trt_magic(), Some(3.0 / 16.0));
    }
    assert!(actual.lattice.grid.moments3(actual.lattice.grid.idx(0, 0, 3)).1[0] > 0.0);
}

#[test]
fn heating_drives_buoyancy_in_the_same_step_without_erasing_external_loads() {
    let mut a = ThermalLbm::slab(3, 8, 0.9, 0.8, 0.01);
    let n = a.gpop.len(); let i = a.grid.idx(1, 4); let before = a.temperature(1, 4);
    let f = vec![[1e-5, -2e-5]; n]; let h = heat(n, 3, 0.002);
    let m = mass(&a.grid.f);
    let b = a.step_with_sources(&h, &f, limits()).unwrap();
    assert!((a.grid.fext[i][1] - (0.01 * (before + 0.001 - 0.5) - 2e-5)).abs() < 1e-14);
    assert_eq!(a.grid.fext[i][0], 1e-5);
    assert!((b.source - 24.0*0.002).abs() < 1e-14);
    assert!(b.defect.abs() < 1e-12); assert!((mass(&a.grid.f)-m).abs() < 1e-10);
    let mut a = ThermalLbm3::slab(3, 2, 8, 0.9, 0.9, 0.01);
    let n = a.gpop.len(); let i = a.grid.idx(1, 1, 4); let before = a.temperature(1, 1, 4);
    let f = vec![[1e-5, 2e-5, -3e-5]; n]; let h = heat(n, 6, 0.002);
    let b = a.step_with_sources(&h, &f, limits()).unwrap();
    assert!((a.grid.fext[i][2] - (0.01 * (before + 0.001 - 0.5) - 3e-5)).abs() < 1e-14);
    assert_eq!(&a.grid.fext[i][..2], &[1e-5, 2e-5]);
    assert!(b.defect.abs() < 1e-12);
}

#[test]
fn signed_sources_and_wall_exchange_close_over_a_heating_cooling_schedule() {
    let mut a = ThermalLbm3::slab(2, 3, 8, 0.9, 0.9, 1e-5);
    let n = a.gpop.len(); let f = vec![[0.0; 3]; n]; let m = mass(&a.grid.f);
    let mut stored = 0.0; let mut supplied = 0.0; let mut walls = 0.0;
    for step in 0..80 {
        let q = if step < 40 { 0.0002 } else { -0.0001 };
        let b = a.step_with_sources(&heat(n, 6, q), &f, limits()).unwrap();
        stored += b.stored_change; supplied += b.source; walls += b.bottom_wall + b.top_wall;
        assert!(b.defect.abs() < 1e-12);
    }
    assert!((stored-supplied-walls).abs() < 1e-10);
    assert!((mass(&a.grid.f)-m).abs() < 1e-10);
    let b = a.step_with_sources(&vec![0.0; n], &f, limits()).unwrap();
    assert_eq!(b.source, 0.0); // Switching off removes the input, not stored heat.
}

fn heated_profile_error(height: usize, three: bool) -> f64 {
    let q = 0.2 / (height*height) as f64;
    let (temperatures, wall) = if three {
        let mut s = ThermalLbm3::slab(1, 1, height, 0.9, 0.9, 0.0);
        s.t_bottom = 0.0; s.t_top = 0.0; s.gpop.fill([0.0; 7]);
        let n = s.gpop.len(); let h = heat(n, 1, q); let f = vec![[0.0; 3]; n];
        let mut wall = 0.0;
        for _ in 0..20*height*height {
            let b = s.step_with_sources(&h, &f, limits()).unwrap(); wall = b.bottom_wall + b.top_wall;
        }
        ((1..=height).map(|z| s.temperature(0, 0, z)).collect::<Vec<_>>(), wall)
    } else {
        let mut s = ThermalLbm::slab(1, height, 0.9, 0.8, 0.0);
        s.t_bottom = 0.0; s.t_top = 0.0; s.gpop.fill([0.0; 5]);
        let n = s.gpop.len(); let h = heat(n, 1, q); let f = vec![[0.0; 2]; n];
        let mut wall = 0.0;
        for _ in 0..20*height*height {
            let b = s.step_with_sources(&h, &f, limits()).unwrap(); wall = b.bottom_wall + b.top_wall;
        }
        ((1..=height).map(|y| s.temperature(0, y)).collect::<Vec<_>>(), wall)
    };
    assert!((wall + q*height as f64).abs() < 1e-8);
    temperatures.iter().enumerate().map(|(i,t)| {
        let z = (i as f64+0.5)/height as f64;
        (t-z*(1.0-z)).abs() // Independent -alpha*T'' = q, zero wall values.
    }).fold(0.0, f64::max)
}
#[test]
fn uniform_heating_converges_to_the_analytic_parabola_in_2d_and_3d() {
    for three in [false, true] {
        let coarse = heated_profile_error(8, three); let fine = heated_profile_error(16, three);
        assert!(coarse < 0.003);
        assert!((3.8..4.2).contains(&(coarse/fine)), "coarse={coarse}, fine={fine}");
    }
}

#[test]
fn refusals_preserve_both_fields_even_after_a_tentative_transport_step() {
    let mut a = ThermalLbm::slab(2, 8, 0.9, 0.8, 0.0);
    let n = a.gpop.len(); let f = vec![[0.0; 2]; n];
    let before = (a.grid.f.clone(), a.gpop.clone(), a.grid.fext.clone());
    let mut h = heat(n, 2, 0.1);
    // First kick stays below 1; the second exceeds the admitted upper bound.
    let mut tight = limits(); tight.temperature_range = [0.0, 1.0];
    assert!(matches!(a.step_with_sources(&h, &f, tight), Err(ThermalSourceError::Temperature { .. })));
    assert_eq!((&a.grid.f, &a.gpop, &a.grid.fext), (&before.0, &before.1, &before.2));
    h[0] = 0.1;
    assert!(a.step_with_sources(&h, &f, limits()).is_err());
    h[0] = 0.0;
    let mut small = limits(); small.max_cells = n-1;
    assert!(a.step_with_sources(&h, &f, small).is_err());
    h[3] = f64::NAN;
    assert!(a.step_with_sources(&h, &f, limits()).is_err());
    assert_eq!((&a.grid.f, &a.gpop, &a.grid.fext), (&before.0, &before.1, &before.2));
    let mut a = ThermalLbm3::slab(2, 2, 8, 0.9, 0.9, 0.0);
    let n = a.gpop.len(); let h = heat(n, 4, 0.001); let mut f = vec![[0.0; 3]; n];
    let before = (a.grid.f.clone(), a.gpop.clone(), a.grid.fext.clone());
    f[4][2] = f64::INFINITY;
    assert!(a.step_with_sources(&h, &f, limits()).is_err());
    assert!(a.step_with_sources(&h[..n-1], &f, limits()).is_err());
    assert_eq!((&a.grid.f, &a.gpop, &a.grid.fext), (&before.0, &before.1, &before.2));
}


#[test]
fn g1_rotating_thermal_sources_preserve_inertial_phase_without_coriolis_heat() {
    use fs_lbm::core2::CollisionModel2;
    use fs_lbm::core2::frame::Coriolis2;
    use fs_lbm::d3q19::Coriolis3;
    use fs_lbm::thermal::BuoyancyQuadrature;
    use fs_lbm::thermal::boundary::ThermalWall;

    let density = 1.7;
    let initial2 = [0.01, -0.02];
    let initial3 = [0.02, -0.015, 0.01];
    let spin2: f64 = 0.075;
    let spin3: [f64; 3] = [0.03, -0.025, 0.04];
    let norm3 = spin3.iter().map(|w| w * w).sum::<f64>().sqrt();
    let axis3 = spin3.map(|w| w / norm3);
    let axial: f64 = axis3.iter().zip(initial3).map(|(a, u)| a * u).sum();
    let cross = [
        axis3[1] * initial3[2] - axis3[2] * initial3[1],
        axis3[2] * initial3[0] - axis3[0] * initial3[2],
        axis3[0] * initial3[1] - axis3[1] * initial3[0],
    ];
    for quadrature in [
        BuoyancyQuadrature::OpeningTemperature,
        BuoyancyQuadrature::EndpointAverage,
    ] {
        for collision in 0..3 {
            let mut two = ThermalLbm::slab(2, 96, 0.9, 0.8, 0.0);
            let mut three = ThermalLbm3::slab(2, 1, 96, 0.9, 0.9, 0.0);
            if collision == 1 {
                two.grid.collision = CollisionModel2::Trt { magic: 3.0 / 16.0 };
                three.grid.set_trt_magic(3.0 / 16.0).unwrap();
            } else if collision == 2 {
                three.grid.set_central_moment_bulk_collision(1.1, 1.3).unwrap();
            }
            two.grid.f.fill(fs_lbm::equilibrium(density, initial2[0], initial2[1]));
            three.grid.f.fill(fs_lbm::equilibrium3(density, initial3));
            two.gpop.fill([1.0 / 6.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0]);
            three.gpop.fill([0.125, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625]);
            two.set_thermal_walls([ThermalWall::Adiabatic; 2]).unwrap();
            three.set_thermal_walls([ThermalWall::Adiabatic; 2]).unwrap();
            two.set_buoyancy_quadrature(quadrature);
            three.set_buoyancy_quadrature(quadrature);
            two.set_coriolis(Coriolis2 { angular_velocity: spin2 }).unwrap();
            three.set_coriolis(Coriolis3 { angular_velocity: spin3 }).unwrap();
            let count = two.gpop.len();
            let heat = vec![0.0; count];
            let force2 = vec![[0.0; 2]; count];
            let force3 = vec![[0.0; 3]; count];
            // The observed cells remain outside the walls' lattice domain of
            // dependence for all 32 ticks. This tests the homogeneous source,
            // not energy conservation of the spatial/boundary discretization.
            for tick in 1..=32 {
                let balance2 = two.step_with_sources(&heat, &force2, limits()).unwrap();
                let balance3 = three.step_with_sources(&heat, &force3, limits()).unwrap();
                for balance in [balance2, balance3] {
                    assert_eq!(balance.source, 0.0);
                    assert_eq!(balance.bottom_wall, 0.0);
                    assert_eq!(balance.top_wall, 0.0);
                    assert!(balance.stored_change.abs() < 1e-11);
                }
                let u2 = two.grid.moments(two.grid.idx(1, 48)).u;
                let u3 = three.grid.moments3(three.grid.idx(1, 0, 48)).1;
                let angle2 = -2.0 * f64::from(tick) * spin2.atan();
                let expected2 = [
                    angle2.cos() * initial2[0] - angle2.sin() * initial2[1],
                    angle2.sin() * initial2[0] + angle2.cos() * initial2[1],
                ];
                let angle3 = -2.0 * f64::from(tick) * norm3.atan();
                for d in 0..2 {
                    assert!((u2[d] - expected2[d]).abs() < 3e-12);
                }
                for d in 0..3 {
                    let expected = angle3.cos() * initial3[d]
                        + angle3.sin() * cross[d]
                        + (1.0 - angle3.cos()) * axial * axis3[d];
                    assert!((u3[d] - expected).abs() < 3e-12);
                }
                let actual_axial: f64 = axis3.iter().zip(u3).map(|(a, u)| a * u).sum();
                assert!((actual_axial - axial).abs() < 3e-12);
            }
            let before2 = two.grid.moments(two.grid.idx(1, 48)).u;
            let before3 = three.grid.moments3(three.grid.idx(1, 0, 48)).1;
            let snapshot2 = (two.grid.f.clone(), two.grid.fext.clone(), two.gpop.clone());
            let snapshot3 = (three.grid.f.clone(), three.grid.fext.clone(), three.gpop.clone());
            for invalid in [f64::NAN, f64::MAX] {
                assert!(two.set_coriolis(Coriolis2 { angular_velocity: invalid }).is_err());
                assert!(three.set_coriolis(Coriolis3 {
                    angular_velocity: [invalid, 0.0, 0.0],
                }).is_err());
            }
            assert_eq!(two.coriolis().angular_velocity, spin2);
            assert_eq!(three.coriolis().angular_velocity, spin3);
            two.set_coriolis(Coriolis2::default()).unwrap();
            three.set_coriolis(Coriolis3::default()).unwrap();
            assert_eq!((two.grid.f.clone(), two.grid.fext.clone(), two.gpop.clone()), snapshot2);
            assert_eq!((three.grid.f.clone(), three.grid.fext.clone(), three.gpop.clone()), snapshot3);
            two.step_with_sources(&heat, &force2, limits()).unwrap();
            three.step_with_sources(&heat, &force3, limits()).unwrap();
            for (actual, expected) in two.grid.moments(two.grid.idx(1, 48)).u.iter().zip(before2) {
                assert!((actual - expected).abs() < 3e-12);
            }
            for (actual, expected) in three.grid.moments3(three.grid.idx(1, 0, 48)).1.iter().zip(before3) {
                assert!((actual - expected).abs() < 3e-12);
            }
        }
    }
}
