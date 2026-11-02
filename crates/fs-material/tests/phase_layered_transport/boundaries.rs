//! Two physical face baths retain separate heat and entropy at one phase store.
use super::*;

#[test]
fn g1_two_face_baths_drive_original_phase_conductivity_and_latent_energy() {
    let card = material(1.0, 4.0);
    let curve = caloric(&card, 400.0, 100.0);
    let law = prepare(&card, &curve);
    for shape in shapes() {
        let body = body(&[&law], &[0.5], &[1], &[], shape);
        let bound = body.bind(ids(), [Some(0.2), Some(0.05)]).unwrap();
        let before = initial(&bound, &[20.0], [Some(350.0), Some(500.0)]);
        let id = bound.cell_ids()[0];
        let original = before.cells().clone();
        let a = start(shape);
        let b = a + 0.5;
        let x = centroid(shape, a, b);
        let mass = volume(shape, a, b);
        let bulk = [metric(shape, a, x), metric(shape, x, b)];
        let film = [0.2 / area(shape, a), 0.05 / area(shape, b)];
        let dt = 0.01;
        // One original latent store, with its source k(h) refreshed at the
        // actual endpoint. Independent dimensional resistance integrals.
        let powers = |heat: f64| {
            let k = 1.0 + 3.0 * (20.0 + heat / mass) / 100.0;
            [
                -50.0 / (bulk[0] / k + film[0]),
                100.0 / (bulk[1] / k + film[1]),
            ]
        };
        let expected = root(0.0, 80.0 * mass, |heat| {
            let q = powers(heat);
            heat - dt * (q[0] + q[1])
        });
        let q = powers(expected);
        let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let (next, report) = chain
            .advance_to(&before, dt, controls(dt), &mut workspace, || Ok(()))
            .unwrap();
        close(
            next.cells()[&id].specific_h_j_kg,
            20.0 + expected / mass,
            2e-8,
        );
        close(next.temperature_k(id).unwrap(), 400.0, 1e-12);
        close(report.energy.stored_energy_change_j, expected, 1e-8);
        close(report.energy.external_heat_j, expected, 1e-8);
        close(report.energy.residual_j, 0.0, 1e-8);
        close(
            report.entropy.thermal_storage_change_j_k,
            expected / 400.0,
            1e-10,
        );
        close(
            report.entropy.reservoir_entropy_input_j_k,
            dt * (q[0] / 350.0 + q[1] / 500.0),
            1e-10,
        );
        close(
            report.entropy.conduction_generation_j_k,
            dt * (q[0] * (1.0 / 400.0 - 1.0 / 350.0) + q[1] * (1.0 / 400.0 - 1.0 / 500.0)),
            1e-10,
        );
        close(report.entropy.temporal_defect_j_k, 0.0, 1e-10);
        assert!(report.entropy.conduction_generation_j_k > 0.0);
        assert_eq!(next.cells()[&id].boundary, before.cells()[&id].boundary);
        assert_eq!(before.cells(), &original);
    }
}

#[test]
fn g1_equal_opposite_bath_heat_preserves_latent_state_and_produces_physical_entropy() {
    let card = material(2.0, 2.0);
    let curve = caloric(&card, 400.0, 100.0);
    let law = prepare(&card, &curve);
    let shape = LayerStackGeometry::Planar { area_m2: 2.0 };
    let body = body(&[&law], &[0.5], &[1], &[], shape);
    let bound = body.bind(ids(), [Some(0.1), Some(0.1)]).unwrap();
    let before = initial(&bound, &[35.0], [Some(350.0), Some(450.0)]);
    let original = before.cells().clone();
    let id = bound.cell_ids()[0];
    let dt = 0.02;
    let inward_hot = 50.0 / (0.25 / (2.0 * 2.0) + 0.1 / 2.0);
    let expected_bath_entropy = dt * inward_hot * (1.0 / 450.0 - 1.0 / 350.0);
    let observer = bound.conduction_observer(&before, || Ok(())).unwrap();
    let (explicit, _, energy, entropy) = before
        .advance_observed_with_entropy(
            &observer,
            dt,
            StepConfig {
                max_dt_s: dt,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap();
    close(explicit.cells()[&id].specific_h_j_kg, 35.0, 1e-12);
    close(energy.external_heat_j, 0.0, 1e-12);
    close(entropy.thermal_storage_change_j_k, 0.0, 1e-12);
    close(
        entropy.reservoir_entropy_input_j_k,
        expected_bath_entropy,
        1e-12,
    );
    close(
        entropy.conduction_generation_j_k,
        -expected_bath_entropy,
        1e-12,
    );
    assert!(entropy.conduction_generation_j_k > 0.0);

    let chain = bound.implicit_chain(&before, || Ok(())).unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let policy = AdaptiveImplicitConductionControls {
        step: controls(dt),
        initial_dt_s: dt,
        maximum_trials: 16,
        enthalpy_absolute_tolerance_j_kg: 1e-7,
        temperature_absolute_tolerance_k: 1e-7,
        phase_fraction_absolute_tolerance: Some(1e-9),
        relative_tolerance: 1e-7,
    };
    let (adaptive, report) = chain
        .advance_adaptive_to(&before, dt, policy, &mut workspace, || Ok(()))
        .unwrap();
    close(adaptive.cells()[&id].specific_h_j_kg, 35.0, 1e-12);
    close(report.solve.energy.external_heat_j, 0.0, 1e-12);
    close(
        report.solve.entropy.reservoir_entropy_input_j_k,
        expected_bath_entropy,
        1e-12,
    );
    close(
        report.solve.entropy.conduction_generation_j_k,
        -expected_bath_entropy,
        1e-12,
    );
    assert_eq!(adaptive.cells()[&id].boundary, before.cells()[&id].boundary);
    assert_eq!(before.cells(), &original);

    assert!(
        bound
            .birth_events(0.0, &[35.0], [Some(350.0), None], || Ok(()))
            .is_err()
    );
    assert!(matches!(
        bound.conduction_observer(&before, || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    assert!(
        chain
            .advance_to(&before, dt, controls(dt), &mut workspace, || Err(
                DepositionError::Cancelled
            ))
            .is_err()
    );
    assert_eq!(before.cells(), &original);
}
