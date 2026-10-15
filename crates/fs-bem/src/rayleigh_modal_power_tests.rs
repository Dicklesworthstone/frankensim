use super::*;
use core::f64::consts::TAU;

fn air() -> Medium { Medium { density: 1.2, sound_speed: 343.0 } }

fn work(z: &[C64], velocity: &[C64]) -> f64 {
    let n = velocity.len();
    let mut result = C64::ZERO;
    for i in 0..n { for j in 0..n {
        result = result + velocity[i].conj() * z[i*n+j] * velocity[j];
    } }
    0.5 * result.re
}

#[test]
fn cancelling_motion_has_positive_power_without_a_diagonal_loss_floor() {
    // Three NONOVERLAPPING equivalent discs (diameter 11.3 mm, spacing
    // 20 mm). Monopole and dipole moments cancel. The old self-cell rule
    // makes this perfectly resolved quadrupole supply acoustic energy.
    let area = 1e-4;
    let distance = 0.02;
    let cells: Vec<_> = [-1.0, 0.0, 1.0].into_iter().zip([1.0, -2.0, 1.0])
        .map(|(x, shape)| AreaSample { position_m: [x*distance, 0.0],
            area_m2: area, mode_shape: vec![shape] }).collect();
    let surface = PreparedSurface::new(&cells).unwrap();
    for hz in [100.0, 400.0, 1000.0] {
        let omega = TAU*hz;
        let old = surface.impedance(omega, air()).unwrap();
        let new = surface.power_consistent_impedance(omega, air()).unwrap();
        assert!(old[0].re < 0.0, "fixture must expose the old anti-damping at {hz} Hz");
        assert!(new[0].re > 0.0, "physical quadrupole must radiate at {hz} Hz");
        assert_eq!(new[0].im, old[0].im, "no change to singular reactance");
        // Independent positive outgoing angular power. With centres -d,0,d,
        // the far-field amplitude is -4*A*sin(k*d*direction_x/2)^2.
        // Integrate its SQUARE directly; do not form a cancelling sinc sum.
        let intervals = 256;
        let mut integral = 0.0;
        for i in 0..=intervals {
            let x = -1.0 + 2.0*i as f64/intervals as f64;
            let weight = if i==0 || i==intervals {1.0} else if i%2==0 {2.0} else {4.0};
            integral += weight * (0.5*omega/air().sound_speed*distance*x).sin().powi(4);
        }
        integral *= 2.0/(3.0*intervals as f64);
        let expected = 2.0*air().density*omega*omega*area*area/(PI*air().sound_speed)*integral;
        assert!((0.5*new[0].re/expected - 1.0).abs() < 2e-6,
            "{hz} Hz: matrix {} vs outgoing {expected}", 0.5*new[0].re);
    }
}

fn surface(fields: usize) -> Vec<AreaSample> {
    (0..9).map(|i| {
        let x = (i%3) as f64 - 1.0;
        let y = (i/3) as f64 - 1.0;
        AreaSample { position_m: [0.04*x, 0.035*y], area_m2: 1e-4*(1.0+0.13*i as f64),
            mode_shape: (0..fields).map(|m| if m==0 {1.0} else {
                ((m+1) as f64*(0.7*x+0.31*y+0.15)).sin()
            }).collect() }
    }).collect()
}

fn far_power(cells: &[AreaSample], velocity: &[C64], omega: f64, medium: Medium) -> f64 {
    // Independent angular pressure integration, not the impedance kernel.
    // Simpson in cos(theta), periodic midpoint trapezoid in azimuth.
    let (polar, azimuth) = (128, 96);
    let local: Vec<_> = cells.iter().map(|cell| {
        velocity.iter().zip(&cell.mode_shape).fold(C64::ZERO,
            |sum, (&v, &shape)| sum + v.scale(cell.area_m2*shape))
    }).collect();
    let mut integral = 0.0;
    for p in 0..=polar {
        let mu = p as f64/polar as f64;
        let transverse = (1.0-mu*mu).sqrt();
        let weight = if p==0 || p==polar {1.0} else if p%2==0 {2.0} else {4.0};
        for a in 0..azimuth {
            let angle = TAU*(a as f64+0.5)/azimuth as f64;
            let direction = [transverse*angle.cos(), transverse*angle.sin()];
            let mut volume = C64::ZERO;
            for (cell, &v) in cells.iter().zip(&local) {
                let phase = -omega/medium.sound_speed*(direction[0]*cell.position_m[0]
                    + direction[1]*cell.position_m[1]);
                volume = volume + v*C64::new(phase.cos(), phase.sin());
            }
            integral += weight*(volume.re*volume.re+volume.im*volume.im);
        }
    }
    integral * TAU/(3.0*polar as f64*azimuth as f64)
        * medium.density*omega*omega/(2.0*medium.sound_speed*TAU*TAU)
}

#[test]
fn complex_modal_work_matches_outgoing_power_with_all_128_fields() {
    for fields in [3, 128] {
        let cells = surface(fields);
        let prepared = PreparedSurface::new(&cells).unwrap();
        let velocity: Vec<_> = (0..fields).map(|i|
            C64::new((0.31*i as f64+0.4).cos(), (0.73*i as f64-0.2).sin())
                .scale(1.0/fields as f64)).collect();
        for hz in [80.0, 700.0, 2000.0] {
            let omega = TAU*hz;
            let z = prepared.power_consistent_impedance(omega, air()).unwrap();
            assert_eq!(z.len(), fields*fields);
            for i in 0..fields {for j in 0..fields {assert_eq!(z[i*fields+j],z[j*fields+i]);}}
            let expected = far_power(&cells, &velocity, omega, air());
            let got = work(&z, &velocity);
            assert!(got > 0.0 && (got/expected-1.0).abs() < 2e-6,
                "{fields} fields at {hz} Hz: {got} vs {expected}");
        }
    }
}

#[test]
fn signed_field_projection_preserves_the_new_power_and_old_reactance() {
    let cells = surface(3);
    let omega = TAU*420.0;
    let p = PreparedSurface::new(&cells).unwrap();
    let new = p.power_consistent_impedance(omega, air()).unwrap();
    let old = p.impedance(omega, air()).unwrap();
    assert!(new.iter().zip(&old).all(|(a,b)|a.im.to_bits()==b.im.to_bits()));
    let columns = [[0.4, -0.7, 0.2], [-0.5, 0.1, 0.9]];
    let projected: Vec<_> = cells.into_iter().map(|mut cell| {
        cell.mode_shape = columns.iter().map(|q|
            q.iter().zip(&cell.mode_shape).map(|(x,y)|x*y).sum()).collect();
        cell
    }).collect();
    let direct = PreparedSurface::new(&projected).unwrap()
        .power_consistent_impedance(omega, air()).unwrap();
    for i in 0..2 { for j in 0..2 {
        let mut expected = C64::ZERO;
        for a in 0..3 { for b in 0..3 {
            expected = expected + new[a*3+b].scale(columns[i][a]*columns[j][b]);
        } }
        assert!((direct[i*2+j]-expected).abs() < 2e-13*new[0].abs());
    } }
}

#[test]
fn self_resistance_is_the_regular_limit_not_a_finite_disk_correction() {
    let cell = AreaSample { position_m: [0.2, -0.3], area_m2: 2e-4,
        mode_shape: vec![-2.0, 0.0, 0.7] };
    let p = PreparedSurface::new(&[cell.clone()]).unwrap();
    for hz in [1e-3, 50.0, 1000.0] {
        let omega = TAU*hz;
        let z = p.power_consistent_impedance(omega, air()).unwrap();
        let factor = air().density*omega*omega/(TAU*air().sound_speed)*cell.area_m2*cell.area_m2;
        for i in 0..3 { for j in 0..3 {
            let expected = factor*cell.mode_shape[i]*cell.mode_shape[j];
            assert!((z[i*3+j].re-expected).abs() < 2e-14*factor);
        } }
        assert!(z[2].re < 0.0, "signed cross entries must NOT be clipped");
        assert_eq!(z[4], C64::ZERO, "zero physical field receives no loss floor");
    }
}

#[test]
fn new_image_preserves_scaling_and_refuses_unresolved_inputs() {
    let cells = surface(3);
    let prepared = PreparedSurface::new(&cells).unwrap();
    let omega = TAU*500.0;
    let base = prepared.power_consistent_impedance(omega, air()).unwrap();
    let mut scaled = cells;
    for cell in &mut scaled { cell.position_m=cell.position_m.map(|v|3.0*v+0.1);cell.area_m2*=9.0; }
    let large = PreparedSurface::new(&scaled).unwrap().power_consistent_impedance(omega/3.0,air()).unwrap();
    let dense = prepared.power_consistent_impedance(omega,Medium {density:2.4,..air()}).unwrap();
    for i in 0..base.len() {
        assert!((large[i]-base[i].scale(9.0)).abs() < 2e-12*base[0].abs());
        assert!((dense[i]-base[i].scale(2.0)).abs() < 2e-13*base[0].abs());
    }
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, TAU*100_000.0] {
        assert!(prepared.power_consistent_impedance(omega,air()).is_err());
    }
    assert!(prepared.power_consistent_impedance(1.0,Medium {density:0.0,..air()}).is_err());
}
