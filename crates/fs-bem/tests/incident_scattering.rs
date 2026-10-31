//! G0/G1/G2/G3/G5: coherent finite-source scattering on the original BEM owner.
use fs_bem::helmholtz::bc::{
    self, Admittance, IncidentField, IncidentSource, IncidentTraces, Monopole, PlaneWave,
    PreparedBoundary, sample_incident_sources,
};
use fs_bem::helmholtz::{self, Formulation, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_math::{c64::C64, det};
use std::f64::consts::PI;

const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 340.0,
};

fn field(surface: &SpherePanels, k: f64, sources: &[IncidentSource]) -> IncidentTraces {
    IncidentTraces::from_sources(surface, k, MEDIUM, sources, 1e-5).unwrap()
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-20),
        "{a:?} != {b:?}"
    );
}

#[test]
fn monopole_retains_reactive_velocity_and_conserved_outward_power() {
    let k = 0.7;
    let source = IncidentSource::Monopole(Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::new(0.003, -0.001),
    });
    let rho_c = MEDIUM.density * MEDIUM.sound_speed;
    let expected_power = rho_c * k * k * (0.003_f64.powi(2) + 0.001_f64.powi(2)) / (8.0 * PI);
    for r in [0.01, 0.2, 3.0, 100.0] {
        let point = [0.0, 0.0, r];
        let value = sample_incident_sources(k, MEDIUM, &[source], point).unwrap();
        close(
            value.velocity_m_s[2] / value.pressure_pa,
            C64::new(1.0, 1.0 / (k * r)).scale(1.0 / rho_c),
            2e-14,
        );
        let power =
            0.5 * (value.pressure_pa * value.velocity_m_s[2].conj()).re * (4.0 * PI * r * r);
        assert!((power / expected_power - 1.0).abs() < 2e-13);
        // Independent central difference of pressure verifies the exp(-iwt)
        // Euler gradient and the sign of the spherical-wave near field.
        let h = r.min(1.0) * 1e-5;
        let plus = sample_incident_sources(k, MEDIUM, &[source], [0.0, 0.0, r + h]).unwrap();
        let minus = sample_incident_sources(k, MEDIUM, &[source], [0.0, 0.0, r - h]).unwrap();
        let velocity = (plus.pressure_pa - minus.pressure_pa) / C64::new(0.0, 2.0 * h * k * rho_c);
        close(velocity, value.velocity_m_s[2], 2e-8);
    }
}

#[test]
fn prepared_robin_batch_preserves_plane_wave_solution_and_replay() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let beta = vec![C64::new(0.4, -0.2); surface.areas().len()];
    let law = Admittance::Local {
        beta: &beta,
        declared_active: false,
    };
    let k = 0.8;
    let waves = [
        PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(1.0, 0.2),
        },
        PlaneWave {
            direction: [1.0, 0.0, 0.0],
            amplitude_pa: C64::new(-0.4, 0.7),
        },
    ];
    let fields: Vec<_> = waves
        .iter()
        .map(|&w| field(&surface, k, &[IncidentSource::PlaneWave(w)]))
        .collect();
    let prepared =
        PreparedBoundary::new(&surface, k, MEDIUM, law, Formulation::BurtonMiller).unwrap();
    let batch = prepared.scatter_many(&fields).unwrap();
    for ((wave, input), batch) in waves.iter().zip(&fields).zip(batch) {
        let legacy =
            bc::solve_scattering_bc(&surface, k, MEDIUM, law, *wave, Formulation::BurtonMiller)
                .unwrap();
        let again = prepared.scatter(input).unwrap();
        assert_eq!(batch.scattered.pressure, again.scattered.pressure);
        assert_eq!(batch.scattered.velocity, again.scattered.velocity);
        assert_eq!(
            batch.absorbed_power_w.to_bits(),
            again.absorbed_power_w.to_bits()
        );
        for (&a, &b) in batch
            .scattered
            .pressure
            .iter()
            .zip(&legacy.scattered.pressure)
        {
            close(a, b, 3e-13);
        }
        assert!((batch.absorbed_power_w / legacy.absorbed_power_w - 1.0).abs() < 3e-13);
    }
}

#[test]
fn prepared_rigid_and_zero_beta_preserve_original_radiation_bits() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let k = 0.8;
    let drive: Vec<_> = surface
        .normals()
        .iter()
        .map(|n| C64::new(n[0], 0.3 * n[2]))
        .collect();
    let original =
        helmholtz::solve_radiation(&surface, k, MEDIUM, &drive, Formulation::BurtonMiller).unwrap();
    let zero = vec![C64::ZERO; drive.len()];
    for law in [
        Admittance::Rigid,
        Admittance::Local {
            beta: &zero,
            declared_active: false,
        },
    ] {
        let prepared =
            PreparedBoundary::new(&surface, k, MEDIUM, law, Formulation::BurtonMiller).unwrap();
        let solution = prepared.solve(&drive).unwrap();
        assert_eq!(solution.radiation.pressure, original.pressure);
        assert_eq!(solution.radiation.velocity, original.velocity);
        assert_eq!(
            solution.radiation.radiated_power.to_bits(),
            original.radiated_power.to_bits()
        );
        assert_eq!(
            solution.radiation.condition_lower_bound.to_bits(),
            original.condition_lower_bound.to_bits()
        );
    }
}

#[test]
fn coherent_sources_obey_superposition_and_total_wall_absorption() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let k = 0.8;
    let beta = vec![C64::new(0.5, 0.1); surface.areas().len()];
    let law = Admittance::Local {
        beta: &beta,
        declared_active: false,
    };
    let sources = [
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, -1.0],
            amplitude_pa: C64::ONE,
        }),
        IncidentSource::Monopole(Monopole {
            position_m: [0.0, 0.0, 3.0],
            volume_velocity_m3_s: C64::new(0.1, -0.02),
        }),
    ];
    let inputs = [
        field(&surface, k, &sources[..1]),
        field(&surface, k, &sources[1..]),
        field(&surface, k, &sources),
    ];
    let prepared =
        PreparedBoundary::new(&surface, k, MEDIUM, law, Formulation::BurtonMiller).unwrap();
    let solutions = prepared.scatter_many(&inputs).unwrap();
    for i in 0..surface.areas().len() {
        close(
            solutions[0].scattered.pressure[i] + solutions[1].scattered.pressure[i],
            solutions[2].scattered.pressure[i],
            2e-13,
        );
    }
    let combined = &solutions[2];
    let transferred = IncidentTraces::from_traces(
        &surface,
        k,
        MEDIUM,
        inputs[2].pressure_pa(),
        inputs[2].normal_velocity_m_s(),
    )
    .unwrap();
    assert_eq!(
        prepared.scatter(&transferred).unwrap().scattered.pressure,
        combined.scattered.pressure
    );
    let observed = prepared
        .scatter_field(IncidentField::Sources(&sources), 1e-5)
        .unwrap();
    assert_eq!(
        observed.solution().scattered.pressure,
        combined.scattered.pressure
    );
    assert!(observed.total_far_field([0.0, 0.0, 1.0]).is_err());
    let total_p = combined.total_pressure_pa();
    let total_v = combined.total_normal_velocity_m_s();
    let rho_c = MEDIUM.density * MEDIUM.sound_speed;
    for i in 0..total_p.len() {
        close(
            total_v[i],
            -(beta[i] * total_p[i]).scale(1.0 / rho_c),
            2e-13,
        );
    }
    let absorption: f64 = total_p
        .iter()
        .zip(surface.areas())
        .map(|(&p, &area)| 0.5 * beta[0].re * p.abs().powi(2) * area / rho_c)
        .sum();
    assert!(absorption > 0.0);
    assert!((combined.absorbed_power_w / absorption - 1.0).abs() < 1e-14);
    assert!(
        (combined.absorbed_power_w - solutions[0].absorbed_power_w - solutions[1].absorbed_power_w)
            .abs()
            > 0.05 * combined.absorbed_power_w,
        "coherent interference must survive the power calculation"
    );

    let opposite = [
        sources[0],
        IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, -1.0],
            amplitude_pa: -C64::ONE,
        }),
    ];
    let zero = prepared.scatter(&field(&surface, k, &opposite)).unwrap();
    assert!(zero.scattered.pressure.iter().all(|&p| p == C64::ZERO));
    assert_eq!(zero.absorbed_power_w, 0.0);
}

// Independent spherical-Bessel series and upward Neumann recurrence.
fn spherical(x: f64, count: usize) -> (Vec<f64>, Vec<f64>) {
    let mut j = Vec::new();
    let mut leading = 1.0;
    for n in 0..count {
        if n > 0 {
            leading *= x / (2 * n + 1) as f64;
        }
        let (mut term, mut sum) = (leading, leading);
        for m in 1..60 {
            term *= -x * x / (2 * m * (2 * n + 2 * m + 1)) as f64;
            sum += term;
        }
        j.push(sum);
    }
    let mut y = vec![-det::cos(x) / x, -det::cos(x) / (x * x) - det::sin(x) / x];
    for n in 1..count - 1 {
        y.push((2 * n + 1) as f64 / x * y[n] - y[n - 1]);
    }
    (j, y)
}

// NIST DLMF https://dlmf.nist.gov/10.60 (10.60.1–2):
// e^(ik|x-s|)/|x-s| = ik Σ(2l+1) j_l(kr) h_l(kR) P_l(mu)
// for r<R. Enforce p'_total + ik*beta*p_total = 0 at the unit sphere, then
// take h_l(kr) ~ (-i)^(l+1)e^(ikr)/(kr). No BEM kernel enters this oracle.
fn monopole_sphere_far_field(k: f64, source_r: f64, q: C64, beta: C64, mu: f64) -> C64 {
    let (j, y) = spherical(k, 26);
    let (sj, sy) = spherical(k * source_r, 26);
    let (mut previous, mut p) = (1.0, 1.0);
    let mut phase = C64::new(0.0, -1.0);
    let mut sum = C64::ZERO;
    for n in 0..24 {
        if n == 1 {
            p = mu;
        } else if n > 1 {
            let next = ((2 * n - 1) as f64 * mu * p - (n - 1) as f64 * previous) / n as f64;
            previous = p;
            p = next;
        }
        let jd = if n == 0 {
            -j[1]
        } else {
            j[n - 1] - (n + 1) as f64 * j[n] / k
        };
        let yd = if n == 0 {
            -y[1]
        } else {
            y[n - 1] - (n + 1) as f64 * y[n] / k
        };
        let ib = C64::new(0.0, 1.0) * beta;
        let ratio =
            -(C64::from_re(jd) + ib.scale(j[n])) / (C64::new(jd, yd) + ib * C64::new(j[n], y[n]));
        sum = sum + C64::new(sj[n], sy[n]) * ratio * phase.scale((2 * n + 1) as f64 * p);
        phase = phase * C64::new(0.0, -1.0);
    }
    sum * q.scale(MEDIUM.density * MEDIUM.sound_speed * k / (4.0 * PI))
}

#[test]
fn finite_source_robin_sphere_converges_to_independent_partial_waves() {
    let k = 0.8;
    let q = C64::new(0.03, -0.01);
    let beta = C64::new(0.4, 0.15);
    let source = IncidentSource::Monopole(Monopole {
        position_m: [0.0, 0.0, 3.0],
        volume_velocity_m3_s: q,
    });
    let directions = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]];
    let expected = directions.map(|d| monopole_sphere_far_field(k, 3.0, q, beta, d[2]));
    let mut errors = Vec::new();
    for level in 0..=2 {
        let surface = SpherePanels::icosphere(1.0, level).unwrap();
        let betas = vec![beta; surface.areas().len()];
        let law = Admittance::Local {
            beta: &betas,
            declared_active: false,
        };
        let input = field(&surface, k, &[source]);
        // Below the first interior resonance, the original triangle-CBIE lane
        // supplies a separate geometric integration axis from incident sampling.
        let prepared =
            PreparedBoundary::new(&surface, k, MEDIUM, law, Formulation::PlainCbie).unwrap();
        let solution = prepared.scatter(&input).unwrap();
        let observed =
            helmholtz::far_field_integrated(&surface, &solution.scattered, MEDIUM, &directions);
        let error = (observed
            .iter()
            .zip(expected)
            .map(|(&a, b)| (a - b).abs().powi(2))
            .sum::<f64>()
            / expected.iter().map(|z| z.abs().powi(2)).sum::<f64>())
        .sqrt();
        eprintln!(
            "monopole Robin sphere level={level} panels={} relative_far_field_error={error:.9e}",
            surface.areas().len()
        );
        assert!(solution.absorbed_power_w > 0.0);
        errors.push(error);
    }
    assert!(
        errors[1] < 0.65 * errors[0] && errors[2] < 0.65 * errors[1],
        "{errors:?}"
    );
    assert!(
        errors[2] < 0.03,
        "finite-source scattering spatial ladder: {errors:?}"
    );
}

#[test]
fn incident_admission_rejects_interior_sources_singularities_and_cross_wiring() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let k = 0.8;
    for position in [[0.0; 3], surface.triangles().unwrap()[0][0]] {
        let source = IncidentSource::Monopole(Monopole {
            position_m: position,
            volume_velocity_m3_s: C64::ONE,
        });
        assert!(IncidentTraces::from_sources(&surface, k, MEDIUM, &[source], 1e-4).is_err());
        assert!(sample_incident_sources(k, MEDIUM, &[source], position).is_err());
    }
    let input = field(
        &surface,
        k,
        &[IncidentSource::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::ONE,
        })],
    );
    let other = SpherePanels::icosphere(1.1, 1).unwrap();
    let prepared = PreparedBoundary::new(
        &other,
        k,
        MEDIUM,
        Admittance::Rigid,
        Formulation::BurtonMiller,
    )
    .unwrap();
    assert!(prepared.scatter(&input).is_err());
    let prepared = PreparedBoundary::new(
        &surface,
        0.9,
        MEDIUM,
        Admittance::Rigid,
        Formulation::BurtonMiller,
    )
    .unwrap();
    assert!(prepared.scatter(&input).is_err());
    let other_medium = Medium {
        density: 1.5,
        ..MEDIUM
    };
    let prepared = PreparedBoundary::new(
        &surface,
        k,
        other_medium,
        Admittance::Rigid,
        Formulation::BurtonMiller,
    )
    .unwrap();
    assert!(prepared.scatter(&input).is_err());
    assert!(IncidentTraces::from_traces(&surface, k, MEDIUM, &[], &[]).is_err());
}

#[test]
fn prepared_source_fields_preserve_total_receivers_and_bounded_far_field() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let k = 0.8;
    let source = Monopole {
        position_m: [0.0, 0.0, 3.0],
        volume_velocity_m3_s: C64::new(0.001, -0.0002),
    };
    let finite = [source];
    let mixed = [IncidentSource::Monopole(source)];
    let prepared = PreparedBoundary::new(
        &surface,
        k,
        MEDIUM,
        Admittance::Rigid,
        Formulation::PlainCbie,
    )
    .unwrap();
    let reused = prepared
        .scatter_field(IncidentField::Sources(&mixed), 1e-5)
        .unwrap();
    let original = bc::incident::solve_scattering(
        &surface,
        k,
        MEDIUM,
        Admittance::Rigid,
        IncidentField::Monopoles(&finite),
        1e-5,
        Formulation::PlainCbie,
    )
    .unwrap();
    assert_eq!(
        reused.solution().scattered.pressure,
        original.solution().scattered.pressure
    );
    assert_eq!(
        reused.total_far_field([1.0, 0.0, 0.0]).unwrap(),
        original.total_far_field([1.0, 0.0, 0.0]).unwrap()
    );
    let points = [[2.0, 0.3, 0.1], [-1.5, 0.1, 0.7]];
    let a = reused
        .total_pressure_at(&points, 1e-5, fs_bem::near_field::Options::default())
        .unwrap();
    let b = original
        .total_pressure_at(&points, 1e-5, fs_bem::near_field::Options::default())
        .unwrap();
    assert_eq!(a.pressure, b.pressure);
}
