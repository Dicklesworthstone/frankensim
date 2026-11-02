//! Exact positive P1 reference-density row moments for native enthalpy storage.

/// Integral rho0*N_a over one admitted affine tetrahedron. Normalize before
/// summing densities, so representable mass does not require a representable
/// unscaled density sum. The caller owns mesh indexing, budgets and cancellation.
pub(super) fn p1_nodal_masses(volume: f64, density: [f64; 4]) -> Option<[f64; 4]> {
    if !volume.is_finite()
        || volume <= 0.0
        || density.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
    {
        return None;
    }
    let masses = if density
        .iter()
        .all(|rho| rho.to_bits() == density[0].to_bits())
    {
        // Match the original cellwise density operation tree exactly.
        [density[0] * volume / 4.0; 4]
    } else {
        let scale = density.into_iter().fold(0.0_f64, f64::max);
        let normalized = density.map(|rho| rho / scale);
        let sum: f64 = normalized.iter().sum();
        let mass = volume * (scale * (0.25 * sum));
        normalized.map(|rho| (mass / 5.0) * (1.0 + rho / sum))
    };
    masses
        .iter()
        .all(|mass| mass.is_finite() && *mass > 0.0)
        .then_some(masses)
}

/// Exact signed action of `integral(N_a N_b dV)` on a nodal density
/// direction. The same symmetric action contracts a storage multiplier back
/// to nodal density. Zero directions are valid; density admission belongs to
/// the positive reference profile, not its derivative.
pub(super) fn p1_nodal_mass_direction(volume: f64, direction: [f64; 4]) -> Option<[f64; 4]> {
    if !volume.is_finite() || volume <= 0.0 || direction.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let scale = direction.into_iter().map(f64::abs).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Some([0.0; 4]);
    }
    let normalized = direction.map(|value| value / scale);
    let sum: f64 = normalized.iter().sum();
    let scaled_volume = volume * scale;
    let action = normalized.map(|value| {
        let coefficient = (sum + value) / 20.0;
        if scaled_volume.is_finite() {
            // Preserve a subnormal direction when geometric scaling makes
            // its final contribution representable before applying the weight.
            scaled_volume * coefficient
        } else {
            volume * (scale * coefficient)
        }
    });
    action
        .iter()
        .all(|value| value.is_finite())
        .then_some(action)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2e-14 * expected.abs().max(1.0),
            "{actual:.17e} != {expected:.17e}"
        );
    }

    #[test]
    fn graded_row_masses_match_independent_degree_two_quadrature() {
        let large = 0.585_410_196_624_968_5;
        let small = 0.138_196_601_125_010_5;
        for (volume, density) in [
            (1.0 / 6.0, [3.0, 15.0, 3.0, 3.0]),
            (0.23, [0.01, 1.0, 8.0, 4.0]),
            (1.2, [14.0, 2.0, 7.0, 6.0]),
        ] {
            let actual = p1_nodal_masses(volume, density).unwrap();
            let mut integral = [0.0; 4];
            for point in 0..4 {
                let shape: [f64; 4] =
                    std::array::from_fn(|a| if a == point { large } else { small });
                let rho: f64 = shape.iter().zip(density).map(|(n, rho)| n * rho).sum();
                for a in 0..4 {
                    integral[a] += (volume / 4.0) * rho * shape[a];
                }
            }
            for a in 0..4 {
                close(actual[a], integral[a]);
            }
            close(
                actual.iter().sum(),
                volume * density.iter().sum::<f64>() / 4.0,
            );
        }
    }

    #[test]
    fn graded_storage_retains_total_enthalpy_under_exact_mesh_partitioning() {
        // On the unit right tet, rho=3+12x and h=40+20x-8y+12z.
        // Independent Cartesian monomials: integral 1=1/6, x=1/24,
        // x²=1/60, xy=1/120. Their product integrates to 47.4 J.
        let density = [3.0, 15.0, 3.0, 3.0];
        let enthalpy = [40.0, 60.0, 32.0, 52.0];
        let mass = p1_nodal_masses(1.0 / 6.0, density).unwrap();
        for (actual, expected) in mass.into_iter().zip([0.225, 0.325, 0.225, 0.225]) {
            close(actual, expected);
        }
        let energy: f64 = mass.iter().zip(enthalpy).map(|(m, h)| m * h).sum();
        close(energy, 47.4);
        let averaged: f64 = p1_nodal_masses(1.0 / 6.0, [6.0; 4])
            .unwrap()
            .into_iter()
            .zip(enthalpy)
            .map(|(m, h)| m * h)
            .sum();
        assert!((energy - averaged).abs() > 1.0);

        // Join the centroid to each original face: four volume-V/4 tets.
        // The original P1 density and h interpolate without changing fields.
        let mut refined_mass = 0.0;
        let mut refined_energy = 0.0;
        for opposite in 0..4 {
            let mut rho = density;
            let mut h = enthalpy;
            rho[opposite] = density.iter().sum::<f64>() / 4.0;
            h[opposite] = enthalpy.iter().sum::<f64>() / 4.0;
            let local = p1_nodal_masses(1.0 / 24.0, rho).unwrap();
            refined_mass += local.iter().sum::<f64>();
            refined_energy += local.iter().zip(h).map(|(m, h)| m * h).sum::<f64>();
        }
        close(refined_mass, mass.iter().sum());
        close(refined_energy, energy);
    }

    #[test]
    fn uniform_storage_bits_and_representable_extreme_scales_are_preserved() {
        let original: f64 = 7.3 * 0.17 / 4.0;
        assert!(
            p1_nodal_masses(0.17, [7.3; 4])
                .unwrap()
                .iter()
                .all(|mass| mass.to_bits() == original.to_bits())
        );
        let normalized = [0.9, 0.8, 0.7, 0.6];
        let large = normalized.map(|rho| rho * f64::MAX);
        assert!(!large.iter().sum::<f64>().is_finite());
        let actual = p1_nodal_masses(1.0 / f64::MAX, large).unwrap();
        let reference = p1_nodal_masses(1.0, normalized).unwrap();
        for a in 0..4 {
            close(actual[a], reference[a]);
        }
        for density in [
            [1.0, 0.0, 1.0, 1.0],
            [1.0, -1.0, 1.0, 1.0],
            [1.0, f64::NAN, 1.0, 1.0],
            [1.0, f64::INFINITY, 1.0, 1.0],
        ] {
            assert!(p1_nodal_masses(1.0, density).is_none());
        }
        assert!(p1_nodal_masses(f64::MAX, [2.0; 4]).is_none());
        assert!(p1_nodal_masses(f64::MIN_POSITIVE, [f64::MIN_POSITIVE; 4]).is_none());
    }

    #[test]
    fn signed_density_action_and_transpose_match_independent_quadrature() {
        let volume = 0.23;
        let direction = [0.7, -1.3, 0.2, -0.4];
        let multiplier = [-0.6, 0.8, 1.1, -0.2];
        let large = 0.585_410_196_624_968_5;
        let small = 0.138_196_601_125_010_5;
        let mut expected = [0.0; 4];
        for point in 0..4 {
            let shape: [f64; 4] = std::array::from_fn(|a| if a == point { large } else { small });
            let value: f64 = shape.iter().zip(direction).map(|(n, d)| n * d).sum();
            for a in 0..4 {
                expected[a] += (volume / 4.0) * shape[a] * value;
            }
        }
        let forward = p1_nodal_mass_direction(volume, direction).unwrap();
        let transpose = p1_nodal_mass_direction(volume, multiplier).unwrap();
        for a in 0..4 {
            close(forward[a], expected[a]);
        }
        close(
            forward.iter().zip(multiplier).map(|(a, b)| a * b).sum(),
            transpose.iter().zip(direction).map(|(a, b)| a * b).sum(),
        );
        let density = [2.0, 4.0, 1.0, 3.0];
        let perturbed =
            p1_nodal_masses(volume, std::array::from_fn(|a| density[a] + direction[a])).unwrap();
        let base = p1_nodal_masses(volume, density).unwrap();
        for a in 0..4 {
            close(perturbed[a] - base[a], forward[a]);
        }
        assert_eq!(p1_nodal_mass_direction(volume, [0.0; 4]), Some([0.0; 4]));
        assert!(p1_nodal_mass_direction(volume, [f64::NAN, 0.0, 0.0, 0.0]).is_none());
        assert!(p1_nodal_mass_direction(f64::MAX, [f64::MAX; 4]).is_none());
        let extreme =
            p1_nodal_mass_direction(1.0 / f64::MAX, [f64::MAX, -f64::MAX, f64::MAX, -f64::MAX])
                .unwrap();
        for (actual, expected) in extreme.into_iter().zip([0.05, -0.05, 0.05, -0.05]) {
            close(actual, expected);
        }
        for sign in [-1.0, 1.0] {
            let subnormal = sign * f64::from_bits(1);
            let expected = (1e300 * subnormal) / 4.0;
            assert_ne!(expected, 0.0);
            assert_eq!(
                p1_nodal_mass_direction(1e300, [subnormal; 4]),
                Some([expected; 4])
            );
        }
    }
}
