use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.15e} != {b:.15e} (tol {tolerance})"
    );
}
fn draw(state: &mut u64) -> f64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut x = *state;
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^= x >> 31;
    (f64::from((x >> 32) as u32) + 0.5) / 4_294_967_296.0
}
fn walk(model: SmithGgx, incident: Vec3, state: &mut u64) -> (u32, Vec3, f64) {
    let mut direction = incident;
    let mut depth = 0.0;
    let mut log_density = 0.0;
    for count in 0..256 {
        let step = model
            .sample(direction, depth, [draw(state), draw(state), draw(state)])
            .unwrap();
        log_density += step.log_density;
        match step.event {
            SmithEvent::Escape => {
                assert!(direction.z > 0.0);
                assert!(log_density.is_finite());
                return (count, direction, log_density);
            }
            SmithEvent::Reflection {
                depth: next_depth,
                direction: next,
                facet,
            } => {
                assert!(facet.z > 0.0 && direction.dot(facet) < 0.0);
                // An upward heightfield normal always increases the ray's z.
                assert!(next.z >= direction.z - 8.0 * f64::EPSILON);
                close(next.dot(next), 1.0, 8e-15);
                direction = next;
                depth = next_depth;
            }
        }
    }
    panic!("a geometric path exceeded the test's explicit 256-facet limit")
}

#[test]
fn heights_obey_smith_survival_atoms_and_signed_rates() {
    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
    // At alpha=1: downward mu=-0.5 has rate 1.5; upward mu=0.5
    // has rate 0.5. These are direct specializations of Smith Eq. (30).
    let radial = det::sqrt(0.75);
    let downward = model
        .sample(Vec3::new(radial, 0.0, -0.5), 2.0, [0.5, 0.3, 0.7])
        .unwrap();
    let SmithEvent::Reflection { depth, .. } = downward.event else {
        panic!("downward escape")
    };
    close(depth, 2.0 + det::ln(2.0) / 1.5, 2e-15);
    let upward = model
        .sample(Vec3::new(radial, 0.0, 0.5), 2.0, [0.5, 0.3, 0.7])
        .unwrap();
    let SmithEvent::Reflection { depth, .. } = upward.event else {
        panic!("early escape")
    };
    close(depth, 2.0 - det::ln(2.0) / 0.5, 2e-15);
    let escape = model
        .sample(Vec3::new(radial, 0.0, 0.5), 2.0, [0.1, 0.3, 0.7])
        .unwrap();
    assert_eq!(escape.event, SmithEvent::Escape);
    close(escape.log_density, -1.0, 2e-15);
    let normal = model
        .sample(Vec3::new(0.0, 0.0, 1.0), 1e300, [0.5, 0.3, 0.7])
        .unwrap();
    assert_eq!(
        normal,
        SmithStep {
            event: SmithEvent::Escape,
            log_density: 0.0
        }
    );
    let horizon = model
        .sample(Vec3::new(1.0, 0.0, 0.0), 2.0, [0.1, 0.3, 0.7])
        .unwrap();
    let SmithEvent::Reflection { depth, .. } = horizon.event else {
        panic!("horizontal escape")
    };
    assert_eq!(depth, 2.0);

    let mut escaped = 0;
    let count = 16_384;
    for i in 0..count {
        let u = (i as f64 + 0.5) / count as f64;
        if matches!(
            model
                .sample(Vec3::new(radial, 0.0, 0.5), 2.0, [u, 0.3, 0.7])
                .unwrap()
                .event,
            SmithEvent::Escape
        ) {
            escaped += 1;
        }
    }
    close(
        escaped as f64 / count as f64,
        det::exp(-1.0),
        1.0 / count as f64,
    );
}

#[test]
fn full_vndf_matches_angular_quadrature_on_both_sides_of_horizon() {
    let model = SmithGgx::try_new(0.4, 0.8).unwrap();
    for mu in [-0.8_f64, 0.0, 0.7] {
        let incoming = Vec3::new(det::sqrt(1.0 - mu * mu), 0.0, mu);
        let wo = incoming.scale(-1.0);
        let area = model.projected_area(wo);
        let mut quadrature = [0.0; 4];
        // Independent integration of D(half)/(4 projected_area) over outgoing
        // solid angle, including directions below the macro horizon.
        const NZ: usize = 320;
        const NP: usize = 640;
        let measure = 4.0 * PI / (NZ * NP) as f64;
        for iz in 0..NZ {
            let z = -1.0 + 2.0 * (iz as f64 + 0.5) / NZ as f64;
            let radial = det::sqrt(1.0 - z * z);
            for ip in 0..NP {
                let phi = 2.0 * PI * (ip as f64 + 0.5) / NP as f64;
                let outgoing = Vec3::new(radial * det::cos(phi), radial * det::sin(phi), z);
                let half = Vec3::new(wo.x + outgoing.x, wo.y + outgoing.y, wo.z + outgoing.z);
                if half.z <= 0.0 {
                    continue;
                }
                let facet = unit(half).unwrap();
                let density = model.ndf(facet) / (4.0 * area);
                let weight = density * measure;
                for (sum, observable) in quadrature.iter_mut().zip([
                    1.0,
                    outgoing.x,
                    outgoing.z,
                    outgoing.z * outgoing.z,
                ]) {
                    *sum += weight * observable;
                }
            }
        }
        close(quadrature[0], 1.0, 0.0035);
        let mut moments = [0.0; 3];
        const N: usize = 262_144;
        let mut state = 23_u64;
        for _ in 0..N {
            let facet = model
                .visible_normal(wo, draw(&mut state), draw(&mut state))
                .unwrap();
            let c = wo.dot(facet);
            let outgoing = Vec3::new(
                2.0 * c * facet.x - wo.x,
                2.0 * c * facet.y - wo.y,
                2.0 * c * facet.z - wo.z,
            );
            assert!(c > 0.0 && facet.z > 0.0);
            for (sum, value) in
                moments
                    .iter_mut()
                    .zip([outgoing.x, outgoing.z, outgoing.z * outgoing.z])
            {
                *sum += value / N as f64;
            }
        }
        for (sampled, integrated) in moments.into_iter().zip(quadrature[1..].iter()) {
            close(sampled, *integrated, 0.004);
        }
    }
}

#[test]
fn joint_path_density_integrates_to_correlated_single_scattering() {
    let model = SmithGgx::try_new(0.55, 0.8).unwrap();
    let incoming = unit(Vec3::new(0.4, -0.1, -0.9)).unwrap();
    let wo = incoming.scale(-1.0);
    let area = model.projected_area(wo);
    let incident_rate = area / wo.z;
    let mut integral = 0.0;
    let mut expected = 0.0;
    const N: usize = 16_384;
    const END: f64 = 32.0;
    for i in 0..N {
        let depth = END * (i as f64 + 0.5) / N as f64;
        let first = model
            .sample(incoming, 0.0, [det::exp(-incident_rate * depth), 0.2, 0.7])
            .unwrap();
        let SmithEvent::Reflection {
            depth: actual,
            facet,
            direction,
        } = first.event
        else {
            panic!("entry escaped")
        };
        close(actual, depth, 1e-13);
        assert!(direction.z > 0.0);
        let exit_rate = model.projected_area(direction.scale(-1.0)) / direction.z;
        let second = model
            .sample(
                direction,
                actual,
                [det::exp(-exit_rate * actual - 1.0), 0.3, 0.6],
            )
            .unwrap();
        assert_eq!(second.event, SmithEvent::Escape);
        integral += det::exp(first.log_density + second.log_density) * END / N as f64;
        expected = model.ndf(facet) / (4.0 * wo.z * (incident_rate + exit_rate));
    }
    close(integral, expected, 2e-7);
}

// At alpha=1 the reflected cosine nu is uniform on (mu,1). Integrating
// the first two exponential height transitions analytically leaves this
// independent 2D angular quadrature for exactly TWO facet reflections.
fn double_scattering_probability() -> f64 {
    const N: usize = 512;
    let mut total = 0.0;
    for i in 0..N {
        let mu = -1.0 + 2.0 * (i as f64 + 0.5) / N as f64;
        let lower = mu.max(0.0);
        let a = (1.0 - mu) / (2.0 * mu.abs());
        for j in 0..N {
            let nu = lower + (1.0 - lower) * (j as f64 + 0.5) / N as f64;
            let b = (1.0 - nu) / (2.0 * nu);
            let height_integral = if mu < 0.0 {
                a / ((1.0 + b) * (a + b))
            } else {
                a / ((1.0 + a) * (1.0 + b))
            };
            total += height_integral * (1.0 - lower) / ((1.0 - mu) * (N * N) as f64);
        }
    }
    total
}

#[test]
fn perfect_reflector_escapes_and_absorbing_orders_match_independent_reference() {
    let model = SmithGgx::try_new(1.0, 1.0).unwrap();
    let mut state = 97_u64;
    let mut first = 0;
    let mut second = 0;
    let mut reflected = 0.0;
    let mut mean_facets = 0.0;
    let mut max_facets = 0;
    const N: usize = 262_144;
    const RHO: f64 = 0.1;
    for _ in 0..N {
        let (count, _, _) = walk(model, Vec3::new(0.0, 0.0, -1.0), &mut state);
        first += usize::from(count == 1);
        second += usize::from(count == 2);
        reflected += det::powi(RHO, count as i32) / N as f64;
        mean_facets += count as f64 / N as f64;
        max_facets = max_facets.max(count);
    }
    // Every unit-Fresnel path escaped with unit power. The first-order
    // component must equal the original correlated Smith white furnace.
    let p1 = 1.0 - det::ln(2.0);
    let p2 = double_scattering_probability();
    close(first as f64 / N as f64, p1, 0.0035);
    close(second as f64 / N as f64, p2, 0.0035);
    // All remaining paths have >=3 reflections. With passive constant facet
    // reflectance rho, their total contribution is bounded by rho^3 P(N>=3).
    let lower = RHO * p1 + RHO * RHO * p2;
    let upper = lower + RHO * RHO * RHO * (1.0 - p1 - p2);
    assert!(
        reflected >= lower - 0.0004 && reflected <= upper + 0.0004,
        "absorbing albedo {reflected}, independent [{lower},{upper}]"
    );
    println!(
        "Smith alpha=1: P1={:.9}, P2={:.9} (reference {p2:.9}), rho=.1 albedo={reflected:.9}, reference=[{lower:.9},{upper:.9}], mean_facets={mean_facets:.6}, max={max_facets}",
        first as f64 / N as f64,
        second as f64 / N as f64
    );
}

#[test]
fn anisotropic_grazing_and_nearly_smooth_paths_retain_unit_fresnel_power() {
    for (ax, ay, mu) in [
        (0.0001, 0.0001, -1.0_f64),
        (0.2, 0.8, -0.3),
        (1.0, 0.1, -1e-8),
    ] {
        let model = SmithGgx::try_new(ax, ay).unwrap();
        let mut state = 41;
        let incident = Vec3::new(det::sqrt(1.0 - mu * mu), 0.0, mu);
        for _ in 0..8192 {
            let (events, escaped, _) = walk(model, incident, &mut state);
            assert!(events >= 1 && escaped.z > 0.0);
        }
    }
}

#[test]
fn invalid_inputs_refuse_without_a_dark_sample() {
    for a in [0.0, -0.1, f64::NAN, f64::INFINITY, 1.01] {
        assert_eq!(SmithGgx::try_new(a, 0.5), Err(SmithError::InvalidRoughness));
    }
    let model = SmithGgx::try_new(0.5, 0.5).unwrap();
    for draw in [
        [0.0, 0.5, 0.5],
        [1.0, 0.5, 0.5],
        [0.5, 1.0, 0.5],
        [0.5, 0.5, f64::NAN],
    ] {
        assert_eq!(
            model.sample(Vec3::new(0.0, 0.0, -1.0), 0.0, draw),
            Err(SmithError::InvalidSample)
        );
    }
    assert_eq!(
        model.sample(Vec3::new(0.0, 0.0, 0.0), 0.0, [0.5; 3]),
        Err(SmithError::InvalidDirection)
    );
    assert_eq!(
        model.sample(Vec3::new(0.0, 0.0, -1.0), -1.0, [0.5; 3]),
        Err(SmithError::InvalidDepth)
    );
}
