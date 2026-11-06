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
