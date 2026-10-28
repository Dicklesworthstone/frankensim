//! Mixed surroundings must combine irradiation before reducing to one temperature.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ExposedFace, RadiationConfig, SolidMaterial, Voxel, VoxelDomain, escape_factors,
    radiated_power, radiative_sinks,
};

fn fixture() -> (VoxelDomain, [SolidMaterial; 1]) {
    let mut domain = VoxelDomain::graded([
        vec![0.008, 0.010, 0.014],
        vec![0.009, 0.012, 0.020],
        vec![0.015, 0.018, 0.010],
    ])
    .unwrap();
    domain.set(1, 1, 1, Voxel::Solid(0));
    (domain, [SolidMaterial::new("radiating block", 10.0)])
}

fn mixed_and_isolated(
    surroundings: [Option<f64>; 6],
) -> (VoxelDomain, Vec<ExposedFace>, Vec<ExposedFace>) {
    let (domain, solids) = fixture();
    let gate = CancelGate::new();
    let mut config = RadiationConfig::new(vec![0.8], surroundings);
    config.surface_exchange = false;
    config.rays_per_face = 512;
    let mixed = escape_factors(&domain, &solids, &config, &gate).unwrap();
    let mut isolated = Vec::new();
    for (side, temperature) in surroundings.into_iter().enumerate() {
        if temperature.is_some() {
            config.surroundings_k = [None; 6];
            config.surroundings_k[side] = temperature;
            isolated.extend(escape_factors(&domain, &solids, &config, &gate).unwrap());
        }
    }
    (domain, mixed, isolated)
}

fn close(actual: f64, expected: f64) {
    let tolerance = 128.0 * f64::EPSILON * actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn mixed_surroundings_power_equals_sum_of_isolated_openings() {
    for surroundings in [
        [220.0, 290.0, 370.0, 540.0, 760.0, 1100.0].map(Some),
        [Some(230.0), None, Some(820.0), None, Some(410.0), None],
    ] {
        // Identical deterministic rays reach the same boundary in each run.
        // Opening one boundary at a time therefore supplies independent
        // isothermal reservoir contributions, without sampling uncertainty.
        let (domain, mixed, isolated) = mixed_and_isolated(surroundings);
        assert!(!mixed.is_empty());
        let mut sees_multiple_temperatures = false;
        for face in &mixed {
            let contributions: Vec<_> = isolated
                .iter()
                .filter(|other| other.cell == face.cell && other.face == face.face)
                .collect();
            sees_multiple_temperatures |= contributions.len() > 1;
            let escape = contributions.iter().map(|f| f.escape).sum::<f64>();
            close(face.escape, escape);
            let fourth_moment = contributions
                .iter()
                .map(|f| f.escape * f.surroundings_k.powi(4))
                .sum::<f64>()
                / escape;
            close(face.surroundings_k.powi(4), fourth_moment);
        }
        assert!(sees_multiple_temperatures);
        for surface_k in [280.0, 650.0, 1400.0] {
            let temperature = vec![surface_k; domain.cell_count()];
            close(
                radiated_power(&domain, &mixed, &temperature),
                radiated_power(&domain, &isolated, &temperature),
            );
        }
    }
}

#[test]
fn mixed_surroundings_newton_tangent_matches_isolated_openings() {
    let (domain, mixed, isolated) =
        mixed_and_isolated([220.0, 290.0, 370.0, 540.0, 760.0, 1100.0].map(Some));
    for surface_k in [280.0, 650.0, 1400.0] {
        let temperature = vec![surface_k; domain.cell_count()];
        let sinks = radiative_sinks(&domain, &mixed, &temperature);
        let reference = radiative_sinks(&domain, &isolated, &temperature);
        let slope = sinks.iter().map(|sink| sink.conductance_w_k).sum::<f64>();
        close(
            slope,
            reference.iter().map(|sink| sink.conductance_w_k).sum(),
        );
        for change_k in [-3.0, 0.0, 4.0] {
            let linearized = |rows: &[fs_lbm::conjugate::CellSink]| {
                rows.iter()
                    .map(|sink| {
                        sink.conductance_w_k * (surface_k + change_k - sink.temperature_k)
                    })
                    .sum::<f64>()
            };
            close(linearized(&sinks), linearized(&reference));
            if change_k == 0.0 {
                close(
                    linearized(&sinks),
                    radiated_power(&domain, &isolated, &temperature),
                );
            }
        }
        let delta_k = 0.01;
        let plus = vec![surface_k + delta_k; domain.cell_count()];
        let minus = vec![surface_k - delta_k; domain.cell_count()];
        let finite_difference = (radiated_power(&domain, &mixed, &plus)
            - radiated_power(&domain, &mixed, &minus))
            / (2.0 * delta_k);
        assert!((finite_difference / slope - 1.0).abs() < 2e-8);
    }
}

#[test]
fn fourth_power_mean_preserves_uniform_and_unseen_temperature_limits() {
    let (domain, solids) = fixture();
    let gate = CancelGate::new();
    let mut config = RadiationConfig::new(vec![0.8], [None; 6]);
    config.surface_exchange = false;
    config.rays_per_face = 512;
    for ambient in [0.25, 303.25, 1e100, f64::MAX / 16.0] {
        config.surroundings_k = [Some(ambient); 6];
        let faces = escape_factors(&domain, &solids, &config, &gate).unwrap();
        assert_eq!(faces.len(), 6);
        for face in faces {
            assert_eq!(face.escape, 1.0);
            assert_eq!(face.surroundings_k.to_bits(), ambient.to_bits());
        }
    }
    // Rays from the block's negative-x face cannot reach positive-x.
    // A very hot unseen boundary must not destroy its cold irradiation
    // through underflow while normalizing the fourth-power mean.
    config.surroundings_k = [Some(0.25); 6];
    config.surroundings_k[1] = Some(1e300);
    let faces = escape_factors(&domain, &solids, &config, &gate).unwrap();
    let cold = faces.iter().find(|face| face.face == 0).unwrap();
    assert_eq!(cold.surroundings_k.to_bits(), 0.25_f64.to_bits());
    assert!(faces.iter().all(|face| face.surroundings_k.is_finite()));
    config.surroundings_k = [None; 6];
    assert!(escape_factors(&domain, &solids, &config, &gate).unwrap().is_empty());
}
