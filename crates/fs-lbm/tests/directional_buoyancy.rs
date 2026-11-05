//! Native inclined-layer/vertical-channel buoyancy uses the same thermal and
//! momentum populations as the original horizontal Rayleigh–Bénard slab.
use fs_lbm::ThermalLbm;
use fs_lbm::d3q19::ThermalLbm3;
use fs_lbm::thermal::ThermalForcingError;
use fs_lbm::thermal::sources::ThermalSourceLimits;

fn limits() -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: 4096,
        temperature_range: [-2.0, 2.0],
        balance_tolerance: 1e-10,
    }
}

#[test]
fn direction_admission_is_scale_independent_and_atomic() {
    let mut two = ThermalLbm::slab(2, 4, 0.9, 0.8, 1e-5);
    let mut three = ThermalLbm3::slab(2, 2, 4, 0.9, 0.9, 1e-5);
    assert_eq!(two.buoyancy_direction(), [0.0, 1.0]);
    assert_eq!(three.buoyancy_direction(), [0.0, 0.0, 1.0]);
    two.set_buoyancy_direction([f64::MAX, 0.0]).unwrap();
    three
        .set_buoyancy_direction([0.0, -f64::from_bits(1), 0.0])
        .unwrap();
    assert_eq!(two.buoyancy_direction(), [1.0, 0.0]);
    assert_eq!(three.buoyancy_direction(), [0.0, -1.0, 0.0]);
    let before2 = (two.grid.f.clone(), two.gpop.clone(), two.grid.fext.clone());
    let before3 = (
        three.grid.f.clone(),
        three.gpop.clone(),
        three.grid.fext.clone(),
    );
    for bad in [[0.0, 0.0], [f64::NAN, 1.0], [1.0, f64::INFINITY]] {
        assert_eq!(
            two.set_buoyancy_direction(bad),
            Err(ThermalForcingError::InvalidDirection)
        );
    }
    for bad in [
        [0.0; 3],
        [1.0, f64::NAN, 0.0],
        [0.0, 0.0, f64::NEG_INFINITY],
    ] {
        assert_eq!(
            three.set_buoyancy_direction(bad),
            Err(ThermalForcingError::InvalidDirection)
        );
    }
    assert_eq!(two.buoyancy_direction(), [1.0, 0.0]);
    assert_eq!(three.buoyancy_direction(), [0.0, -1.0, 0.0]);
    // A finite external force can overflow on a tangential axis after adding
    // buoyancy. Admission must check that axis before mutating either field.
    two.gbeta = f64::MAX;
    three.gbeta = f64::MAX;
    two.t_ref = -0.125;
    three.t_ref = -0.125;
    assert_eq!(
        two.step_with_external_forces(&vec![[f64::MAX, 0.0]; two.gpop.len()]),
        Err(ThermalForcingError::NonFinite {
            cell: 2,
            component: 0
        }),
    );
    assert_eq!(
        three.step_with_external_forces(&vec![[0.0, -f64::MAX, 0.0]; three.gpop.len()]),
        Err(ThermalForcingError::NonFinite {
            cell: 4,
            component: 1
        }),
    );
    assert_eq!((two.grid.f, two.gpop, two.grid.fext), before2);
    assert_eq!((three.grid.f, three.gpop, three.grid.fext), before3);
}

#[test]
fn rotating_native_buoyancy_delivers_one_impulse_through_checked_sources() {
    let mut two = ThermalLbm::slab(3, 32, 0.9, 0.8, 2e-4);
    let mut three = ThermalLbm3::slab(3, 2, 32, 0.9, 0.9, 2e-4);
    // Uniform temperature and a central observation outside the walls' domain
    // of dependence: rho du/dt = g beta (T-T_ref) d, independently of direction.
    two.t_bottom = 1.0;
    two.t_top = 1.0;
    three.t_bottom = 1.0;
    three.t_top = 1.0;
    two.gpop
        .fill([1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]);
    three
        .gpop
        .fill([0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125]);
    let center2 = two.grid.idx(1, 16);
    let center3 = three.grid.idx(1, 1, 16);
    let zero2 = vec![[0.0; 2]; two.gpop.len()];
    let zero3 = vec![[0.0; 3]; three.gpop.len()];
    let heat2 = vec![0.0; two.gpop.len()];
    let heat3 = vec![0.0; three.gpop.len()];
    let mut velocity2 = [0.0; 2];
    let mut velocity3 = [0.0; 3];
    for (d2, d3) in [
        ([1.0, 0.0], [1.0, 0.0, 0.0]),
        ([3.0, 4.0], [2.0, -3.0, 6.0]),
        ([0.0, -1.0], [0.0, 0.0, -1.0]),
        ([-1.0, 0.0], [-1.0, 0.0, 0.0]),
    ] {
        let before2 = two.grid.moments(center2).u;
        let before3 = three.grid.moments3(center3).1;
        two.set_buoyancy_direction(d2).unwrap();
        three.set_buoyancy_direction(d3).unwrap();
        assert_eq!(two.grid.moments(center2).u, before2);
        assert_eq!(three.grid.moments3(center3).1, before3);
        let direction2 = two.buoyancy_direction();
        let direction3 = three.buoyancy_direction();
        two.step_with_sources(&heat2, &zero2, limits()).unwrap();
        three.step_with_sources(&heat3, &zero3, limits()).unwrap();
        assert_eq!(two.buoyancy_direction(), direction2);
        assert_eq!(three.buoyancy_direction(), direction3);
        for a in 0..2 {
            velocity2[a] += 1e-4 * direction2[a];
            assert!((two.grid.moments(center2).u[a] - velocity2[a]).abs() < 2e-14);
        }
        for a in 0..3 {
            velocity3[a] += 1e-4 * direction3[a];
            assert!((three.grid.moments3(center3).1[a] - velocity3[a]).abs() < 2e-14);
        }
    }
}

fn vertical_channel_error(height: usize) -> [f64; 2] {
    let nu = (0.9 - 0.5) / 3.0;
    let h = height as f64;
    let amplitude = 0.03;
    let gbeta = amplitude * nu / (h * h);
    let mut two = ThermalLbm::slab(2, height, 0.9, 0.9, gbeta);
    let mut three = ThermalLbm3::slab(2, 2, height, 0.9, 0.9, gbeta);
    two.set_buoyancy_direction([1.0, 0.0]).unwrap();
    // Rotate into the other tangential axis and reverse gravity. The thermal
    // forcing and continuum viscous profile must transform together.
    three.set_buoyancy_direction([0.0, -1.0, 0.0]).unwrap();
    for _ in 0..(3.0 * h * h / nu) as usize {
        two.step();
        three.step();
    }
    let mut squared = [0.0; 2];
    let mut norm = 0.0;
    for row in 1..=height {
        let s = (row as f64 - 0.5) / h;
        // T=1-s, T_ref=1/2; nu u'' + g beta (1/2-s)=0,
        // u(0)=u(1)=0. No fitted Nusselt correlation enters this solution.
        let exact = amplitude * (s * s * s / 6.0 - s * s / 4.0 + s / 12.0);
        let u2 = two.grid.moments(two.grid.idx(0, row)).u;
        let u3 = three.grid.moments3(three.grid.idx(0, 0, row)).1;
        squared[0] += (u2[0] - exact).powi(2);
        squared[1] += (u3[1] + exact).powi(2);
        norm += exact * exact;
        assert!(u2[1].abs() < 1e-11);
        assert!(u3[0].abs() < 1e-11 && u3[2].abs() < 1e-11);
        assert!((two.temperature(0, row) - (1.0 - s)).abs() < 1e-10);
        assert!((three.temperature(0, 0, row) - (1.0 - s)).abs() < 1e-10);
    }
    assert!((two.nusselt() - 1.0).abs() < 1e-8);
    assert!((three.nusselt() - 1.0).abs() < 1e-8);
    squared.map(|value| fs_math::det::sqrt(value / norm))
}

#[test]
fn heated_vertical_channel_converges_to_analytic_counterflow_in_2d_and_3d() {
    let coarse = vertical_channel_error(12);
    let fine = vertical_channel_error(24);
    eprintln!("vertical-channel relative L2 velocity error: H=12 {coarse:?}, H=24 {fine:?}");
    for dimension in 0..2 {
        assert!(
            fine[dimension] < 0.4 * coarse[dimension],
            "coarse={coarse:?}, fine={fine:?}"
        );
        assert!(fine[dimension] < 0.02, "coarse={coarse:?}, fine={fine:?}");
    }
}
