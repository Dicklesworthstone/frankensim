use super::*;

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g0_births_preserve_particle_number_and_mean_energy_on_original_cells() {
    let grid = AlphaEnergyGrid::new(
        vec![Energy::new(1.0), Energy::new(2.0), Energy::new(5.0)],
        || false,
    )
    .unwrap();
    assert_eq!(grid.widths_j(), &[1.5, 2.0, 3.0]);
    let cloned = grid.clone();
    assert_eq!(cloned.energies().as_ptr(), grid.energies().as_ptr());
    assert_eq!(cloned.widths_j().as_ptr(), grid.widths_j().as_ptr());

    let mut amounts = [0.0; MAX_ALPHA_BINS];
    grid.inject(&mut amounts, Amount::new(7e-6), Energy::new(1.0))
        .unwrap();
    assert_eq!(amounts[0], 7e-6);
    assert_eq!(amounts[1], 0.0);
    grid.inject(&mut amounts, Amount::new(4e-6), Energy::new(3.0))
        .unwrap();
    close(amounts.iter().sum(), 11e-6, 3e-15, 0.0);
    close(
        amounts
            .iter()
            .zip(grid.energies())
            .map(|(n, e)| n * e.value())
            .sum(),
        7e-6 + 3.0 * 4e-6,
        3e-15,
        0.0,
    );
    assert!(amounts[1] > 0.0 && amounts[2] > 0.0);
    assert!(amounts[grid.len()..].iter().all(|&n| n == 0.0));

    // A birth adjacent to a center must retain its small neighbor population.
    let mut tiny = [0.0; MAX_ALPHA_BINS];
    grid.inject(
        &mut tiny,
        Amount::new(1.0),
        Energy::new(f64::from_bits(1.0_f64.to_bits() + 1)),
    )
    .unwrap();
    assert!(tiny[1] > 0.0 && tiny[1] < 1e-15);
    close(tiny.iter().sum(), 1.0, 2e-15, 0.0);

    let before = amounts;
    for energy in [0.5, 5.1] {
        assert_eq!(
            grid.inject(&mut amounts, Amount::new(1.0), Energy::new(energy)),
            Err(DtFusionError::OutsideDomain)
        );
        assert_eq!(amounts, before);
    }
    assert!(matches!(
        grid.inject(&mut amounts, Amount::new(-1.0), Energy::new(3.0)),
        Err(DtFusionError::Invalid(_))
    ));
    assert_eq!(amounts, before);
}

#[test]
fn g1_grid_builders_retain_endpoints_and_explicit_coordinate_spacing() {
    let speed =
        AlphaEnergyGrid::uniform_speed(Energy::new(1.0), Energy::new(25.0), 5, || false).unwrap();
    assert_eq!(
        speed.energies(),
        &[1.0, 4.0, 9.0, 16.0, 25.0].map(Energy::new)
    );
    let log =
        AlphaEnergyGrid::logarithmic(Energy::new(1.0), Energy::new(16.0), 5, || false).unwrap();
    assert_eq!(log.energies()[0], Energy::new(1.0));
    assert_eq!(log.energies()[4], Energy::new(16.0));
    for (actual, expected) in log.energies().iter().zip([1.0, 2.0, 4.0, 8.0, 16.0]) {
        close(actual.value(), expected, 3e-15, 0.0);
    }
    assert_eq!(speed.len(), 5);
    assert!(!speed.is_empty());
}

#[test]
fn g4_grid_domains_budgets_and_cancellation_refuse_without_repair() {
    for centers in [
        vec![],
        vec![Energy::new(1.0)],
        vec![Energy::new(1.0); MAX_ALPHA_BINS + 1],
        vec![Energy::new(1.0), Energy::new(1.0)],
        vec![Energy::new(2.0), Energy::new(1.0)],
        vec![Energy::new(0.0), Energy::new(1.0)],
        vec![Energy::new(1.0), Energy::new(f64::INFINITY)],
    ] {
        assert!(matches!(
            AlphaEnergyGrid::new(centers, || false),
            Err(DtFusionError::Invalid(_))
        ));
    }
    assert_eq!(
        AlphaEnergyGrid::new(vec![Energy::new(1.0), Energy::new(2.0)], || true),
        Err(DtFusionError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        AlphaEnergyGrid::uniform_speed(Energy::new(1.0), Energy::new(25.0), 256, || {
            polls += 1;
            polls > 40
        }),
        Err(DtFusionError::Cancelled)
    );
    let grid = AlphaEnergyGrid::new(vec![Energy::new(1.0), Energy::new(2.0)], || false).unwrap();
    let mut bad = [0.0; MAX_ALPHA_BINS];
    bad[4] = 1.0;
    let before = bad;
    assert!(matches!(
        grid.inject(&mut bad, Amount::new(1.0), Energy::new(1.5)),
        Err(DtFusionError::Invalid(_))
    ));
    assert_eq!(bad, before);
}

#[test]
fn g1_scaled_cells_preserve_birth_number_physical_energy_and_reversible_geometry() {
    let grid = AlphaEnergyGrid::new(
        [1.1e-16, 2.3e-16, 5.7e-16]
            .into_iter()
            .map(Energy::new)
            .collect(),
        || false,
    )
    .unwrap();
    let unchanged = grid.clone();
    let noop = grid.scaled(1.0, || false).unwrap();
    assert_eq!(noop.energies().as_ptr(), grid.energies().as_ptr());
    assert_eq!(noop.widths_j().as_ptr(), grid.widths_j().as_ptr());

    // Multiplying coordinates by powers of two and reversing loses no bins,
    // grid resolution or represented energy. No population remap is involved.
    let compressed = grid.scaled(4.0, || false).unwrap();
    let restored = compressed.scaled(0.25, || false).unwrap();
    assert_eq!(restored, grid);
    for (original, actual) in grid.energies().iter().zip(compressed.energies()) {
        assert_eq!(actual.value(), 4.0 * original.value());
    }

    for scale in [0.25, 1.7, 4.0] {
        let physical = grid.scaled(scale, || false).unwrap();
        grid.admit_scale(scale).unwrap();
        let mut populations = [0.0; MAX_ALPHA_BINS];
        // An exact physical endpoint must be accepted even if division by the
        // scale would round outside the original coordinate's support.
        grid.inject_scaled(
            &mut populations,
            Amount::new(7e-6),
            physical.energies()[0],
            scale,
        )
        .unwrap();
        assert_eq!(populations[0], 7e-6);
        assert!(populations[1..].iter().all(|n| *n == 0.0));
        let birth = 0.3 * physical.energies()[1].value() + 0.7 * physical.energies()[2].value();
        grid.inject_scaled(
            &mut populations,
            Amount::new(4e-6),
            Energy::new(birth),
            scale,
        )
        .unwrap();
        close(populations.iter().sum(), 11e-6, 3e-15, 0.0);
        close(
            populations
                .iter()
                .zip(physical.energies())
                .map(|(n, e)| n * e.value())
                .sum(),
            7e-6 * physical.energies()[0].value() + 4e-6 * birth,
            3e-15,
            0.0,
        );
        assert!(populations[1] > 0.0 && populations[2] > 0.0);

        let before = populations;
        let outside = f64::from_bits(physical.energies()[2].value().to_bits() + 1);
        assert_eq!(
            grid.inject_scaled(
                &mut populations,
                Amount::new(1e-6),
                Energy::new(outside),
                scale,
            ),
            Err(DtFusionError::OutsideDomain)
        );
        assert_eq!(populations, before);
    }
    assert_eq!(grid, unchanged);
}

#[test]
fn g4_scaled_support_refuses_unrepresentable_geometry_and_cancellation() {
    let grid = AlphaEnergyGrid::new(vec![Energy::new(1.0), Energy::new(2.0)], || false).unwrap();
    let original = grid.clone();
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            grid.scaled(scale, || false),
            Err(DtFusionError::Invalid(_))
        ));
    }
    // The last center remains finite, but its extrapolated upper face does not.
    assert_eq!(
        grid.scaled(0.5 * f64::MAX, || false),
        Err(DtFusionError::OutsideDomain)
    );
    let tiny =
        AlphaEnergyGrid::new(vec![Energy::new(1e-300), Energy::new(2e-300)], || false).unwrap();
    assert_eq!(
        tiny.scaled(1e-30, || false),
        Err(DtFusionError::OutsideDomain)
    );
    assert_eq!(grid.scaled(1.0, || true), Err(DtFusionError::Cancelled));
    let mut polls = 0;
    assert_eq!(
        grid.scaled(2.0, || {
            polls += 1;
            polls > 3
        }),
        Err(DtFusionError::Cancelled)
    );
    assert_eq!(grid, original);
}
