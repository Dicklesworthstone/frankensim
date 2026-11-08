//! Independent continuous solutions for complete-spectrum physical damping.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::{
    CompleteSpectralSolid, CompleteSpectralSolidControls, CompleteSpectralSolidError,
    DampedSpectralSolid,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::{Coo, Csr};

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(&gate, arena, StreamKey { seed: 62, kernel_id: 35, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}
fn controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 128, maximum_preparation_bytes: 8 * 1024 * 1024,
        maximum_preparation_multiply_terms: 100_000_000, maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-12, maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12, maximum_work_error_j: 1e-10, maximum_steps: 1000,
    }
}
fn matrix(values: &[f64], n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n { for j in 0..n { coo.push(i, j, values[i * n + j]); } }
    coo.assemble()
}
fn pencil(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    TetElasticAssembly {
        stiffness: matrix(k, n), mass: matrix(m, n), free_dofs: (0..n).collect(),
        total_mass_kg: 1.0, minimum_scaled_jacobian: 1.0, element_volumes_m3: vec![1.0 / 6.0],
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:e} != {expected:e} +/- {tolerance:e}");
}
fn bits(s: &DampedSpectralSolid) -> Vec<u64> {
    s.displacement_m().iter().chain(s.velocity_mps()).copied()
        .chain([s.time_s(), s.stored_energy_j(), s.kinetic_energy_j(),
            s.cumulative_work_j(), s.cumulative_dissipated_j()]).map(f64::to_bits).collect()
}

// Closed-form particular solution plus independent homogeneous pole formulas.
fn exact(w: f64, c: f64, q: f64, v: f64, force: f64, slope: f64, t: f64) -> (f64, f64) {
    if w == 0.0 {
        let a = (force - slope / c) / c;
        let b = slope / c;
        let transient = v - a;
        return (q + a * t + 0.5 * b * t * t - transient * (-c * t).exp_m1() / c,
            a + b * t + transient * (-c * t).exp());
    }
    let b = slope / (w * w);
    let a = (force - c * b) / (w * w);
    let (x, y, decay) = (q - a, v - b, 0.5 * c);
    let (position, velocity) = if c == 2.0 * w {
        let e = (-decay * t).exp();
        (e * (x + (y + decay * x) * t), e * (y - decay * (y + decay * x) * t))
    } else if c < 2.0 * w {
        let frequency = (w * w - decay * decay).sqrt();
        let (sin, cos) = (frequency * t).sin_cos();
        let e = (-decay * t).exp();
        (e * (x * cos + (y + decay * x) * sin / frequency),
            e * (y * cos - (decay * y + w * w * x) * sin / frequency))
    } else {
        let fast = -decay - (decay * decay - w * w).sqrt();
        let slow = w * w / fast;
        let amplitude = (y - fast * x) / (slow - fast);
        let first = amplitude * (slow * t).exp();
        let second = (x - amplitude) * (fast * t).exp();
        (first + second, slow * first + fast * second)
    };
    (a + b * t + position, b + velocity)
}

#[test]
fn all_128_coupled_coordinates_match_affine_motion_work_and_heat() {
    with_cx(false, |cx| {
        let n = 128;
        let shape = [[1.0, 0.5], [0.25, 1.0]];
        let inverse = [[8.0 / 7.0, -4.0 / 7.0], [-2.0 / 7.0, 8.0 / 7.0]];
        let omega: Vec<_> = (0..n).map(|i| if i == 0 { 0.0 } else { 2.0 + 0.75 * (i % 7) as f64 }).collect();
        let rates: Vec<_> = omega.iter().enumerate().map(|(i, w)| {
            if i == 0 { 0.8 } else { [0.1, 2.0, 4.0][i % 3] * w }
        }).collect();
        let lambda: Vec<_> = omega.iter().map(|w| w * w).collect();
        let (mut k, mut m, mut c, mut basis) = (vec![0.0; n*n], vec![0.0; n*n], vec![0.0; n*n], vec![0.0; n*n]);
        for first in (0..n).step_by(2) { for a in 0..2 { for b in 0..2 {
            basis[(first+a)*n+first+b] = shape[a][b];
            for j in 0..2 {
                let product = inverse[j][a] * inverse[j][b];
                m[(first+a)*n+first+b] += product;
                k[(first+a)*n+first+b] += lambda[first+j] * product;
                c[(first+a)*n+first+b] += rates[first+j] * product;
            }
        } } }
        let transform = |values: &[f64], map: [[f64; 2]; 2]| -> Vec<f64> {
            (0..n).map(|i| map[i%2][0]*values[i/2*2] + map[i%2][1]*values[i/2*2+1]).collect()
        };
        let dual = [[inverse[0][0], inverse[1][0]], [inverse[0][1], inverse[1][1]]];
        let q0: Vec<_> = (0..n).map(|i| 0.001 * (1+i%5) as f64).collect();
        let v0: Vec<_> = (0..n).map(|i| 0.002 * (i%7) as f64 - 0.004).collect();
        let force: Vec<_> = (0..n).map(|i| 0.01 + 0.001 * (i%3) as f64).collect();
        let slope: Vec<_> = (0..n).map(|i| 0.003 * (i%5) as f64 - 0.006).collect();
        let source = |t: f64| transform(&(0..n).map(|i| force[i]+slope[i]*t).collect::<Vec<_>>(), dual);
        let mut solid = DampedSpectralSolid::prepare_with_basis(&pencil(&k,&m,n), &matrix(&c,n),
            &lambda, &basis, 0.037, &transform(&q0,shape), &transform(&v0,shape), controls(), cx).unwrap();
        for interval in 0..8 {
            let start = interval as f64 * 0.037;
            let end = start + 0.037;
            let frame = solid.advance_affine(&source(start), &source(end), cx).unwrap();
            let mut q = vec![0.0; n];
            let mut v = vec![0.0; n];
            let (mut work, mut heat) = (0.0, 0.0);
            for i in 0..n {
                (q[i],v[i]) = exact(omega[i],rates[i],q0[i],v0[i],force[i],slope[i],end);
                // Resolved Simpson integration of an independent analytic trajectory.
                for node in 0..=128 {
                    let t = start + (end-start)*node as f64/128.0;
                    let speed = exact(omega[i],rates[i],q0[i],v0[i],force[i],slope[i],t).1;
                    let weight = if node == 0 || node == 128 { 1.0 } else if node%2 == 0 { 2.0 } else { 4.0 };
                    work += weight * (force[i]+slope[i]*t) * speed * (end-start)/384.0;
                    heat += weight * rates[i] * speed*speed * (end-start)/384.0;
                }
            }
            for (actual, expected) in solid.displacement_m().iter().zip(transform(&q,shape)) { near(*actual,expected,2e-14); }
            for (actual, expected) in solid.velocity_mps().iter().zip(transform(&v,shape)) { near(*actual,expected,5e-14); }
            near(frame.energy.supplied_work_j,work,2e-13);
            near(frame.energy.dissipated_j,heat,2e-13);
            assert!(frame.energy.dissipated_j > 0.0);
            assert!(frame.energy.cumulative_balance_error_j.abs() < 1e-12);
        }
        assert_eq!(solid.damping_rates_per_s().len(),n);
    });
}

#[test]
fn rigid_drag_reaches_terminal_velocity_and_high_q_heat_survives() {
    with_cx(false, |cx| {
        let mut rigid = DampedSpectralSolid::prepare_with_basis(&pencil(&[0.0],&[4.0],1),
            &matrix(&[12.0],1), &[0.0], &[0.5], 1.7, &[0.7], &[0.2], controls(), cx).unwrap();
        let frame = rigid.advance(&[8.0],cx).unwrap();
        let (q,v) = exact(0.0,3.0,0.7,0.2,2.0,0.0,1.7);
        near(rigid.displacement_m()[0],q,2e-14);
        near(rigid.velocity_mps()[0],v,2e-14);
        near(frame.energy.supplied_work_j,8.0*(q-0.7),2e-13);
        let terminal: f64 = 2.0/3.0;
        let transient = 0.2-terminal;
        let heat = 12.0*(terminal*terminal*1.7 - 2.0*terminal*transient*(-5.1_f64).exp_m1()/3.0
            - transient*transient*(-10.2_f64).exp_m1()/6.0);
        near(frame.energy.dissipated_j,heat,2e-13);
        let mut high_q = DampedSpectralSolid::prepare_with_basis(&pencil(&[4.0],&[1.0],1),
            &matrix(&[1e-20],1), &[4.0], &[1.0], 0.3, &[0.7], &[-0.2], controls(), cx).unwrap();
        let loss = high_q.advance(&[0.0],cx).unwrap().energy.dissipated_j;
        assert!(loss > 1e-22 && loss < 1e-20);
        assert_eq!(high_q.cumulative_dissipated_j(),loss);
    });
}

#[test]
fn nonclassical_or_negative_viscosity_refuses_without_a_scale_floor() {
    with_cx(false, |cx| {
        let assembly = pencil(&[1.0,0.0,0.0,4.0],&[1.0,0.0,0.0,1.0],2);
        let prepare = |c: &[f64]| DampedSpectralSolid::prepare_with_basis(&assembly,&matrix(c,2),
            &[1.0,4.0], &[1.0,0.0,0.0,1.0], 0.1, &[0.0;2], &[0.0;2], controls(), cx);
        for c in [[1.0,0.1,0.1,1e12], [0.0,1e-100,1e-100,1.0], [-1e-30,0.0,0.0,1.0],
            [1.0,0.1,0.0,1.0], [1.0,f64::NAN,0.0,1.0]] { assert!(prepare(&c).is_err()); }
        let accepted = prepare(&[1.0,1e-14,1e-14,2.0]).unwrap();
        near(accepted.relative_damping_projection_error(),1e-14,1e-29);
    });
}

#[test]
fn impulses_keep_past_heat_and_refused_intervals_are_atomic() {
    with_cx(false, |cx| {
        let mut solid = DampedSpectralSolid::prepare_with_basis(&pencil(&[16.0],&[4.0],1),
            &matrix(&[2.0],1), &[4.0], &[0.5], 0.1, &[0.1], &[0.2], controls(), cx).unwrap();
        solid.advance(&[0.0],cx).unwrap();
        let (heat,old_velocity,time) = (solid.cumulative_dissipated_j(),solid.velocity_mps()[0],solid.time_s());
        let impulse = solid.apply_momentum_impulse(&[0.8],1e-13,cx).unwrap();
        near(solid.velocity_mps()[0],old_velocity+0.2,1e-15);
        near(impulse.energy.supplied_work_j,0.8*(old_velocity+0.1),1e-15);
        assert_eq!(impulse.energy.cumulative_dissipated_j,heat);
        assert_eq!(solid.time_s(),time);
        let before = bits(&solid);
        assert!(solid.advance(&[f64::MAX],cx).is_err());
        assert!(solid.advance(&[],cx).is_err());
        with_cx(true, |cancelled| assert!(matches!(solid.advance(&[0.0],cancelled),Err(CompleteSpectralSolidError::Cancelled))));
        assert_eq!(bits(&solid),before);
        solid.advance(&[0.0],cx).unwrap();
        assert!(solid.cumulative_dissipated_j() > heat);
    });
}

#[test]
fn zero_viscosity_preserves_conservative_motion_and_receipts_bitwise() {
    with_cx(false, |cx| {
        let assembly = pencil(&[4.0],&[1.0],1);
        let mut plain = CompleteSpectralSolid::prepare_with_basis(&assembly,&[4.0],&[1.0],0.1,&[0.7],&[-0.2],controls(),cx).unwrap();
        let mut damped = DampedSpectralSolid::prepare_with_basis(&assembly,&matrix(&[0.0],1),&[4.0],&[1.0],0.1,&[0.7],&[-0.2],controls(),cx).unwrap();
        for i in 0..20 {
            let a = [0.1*i as f64];
            let b = [-0.05*i as f64];
            let expected = plain.advance_affine(&a,&b,cx).unwrap();
            let actual = damped.advance_affine(&a,&b,cx).unwrap();
            assert_eq!(actual.energy,expected.energy);
            assert_eq!(actual.mean_velocity_mps,expected.mean_velocity_mps);
            assert_eq!(plain.displacement_m()[0].to_bits(),damped.displacement_m()[0].to_bits());
            assert_eq!(plain.velocity_mps()[0].to_bits(),damped.velocity_mps()[0].to_bits());
        }
    });
}
