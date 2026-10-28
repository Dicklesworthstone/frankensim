use super::*;

fn balance(e: &FrictionHeatEntropy) {
    close(
        e.thermal_storage_change_j_k,
        e.reservoir_entropy_input_j_k
            + e.network_conduction_generation_j_k
            + e.interface_conduction_generation_j_k
            + e.friction_generation_j_k
            + e.unclassified_source_at_cell_j_k
            + e.temporal_defect_j_k
            + e.interface_residual_entropy_j_k,
        2e-12,
    );
}

#[test]
fn g1_friction_entropy_uses_solved_interface_and_exact_body_storage() {
    let old = network(1, [1.0, 3.0], [2.0, 9.0], [2.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [2.0, 9.0]).unwrap();
    let dt = 0.2;
    let mut lu = factor(&old, dt, [4.0, 6.0], 0.0);
    let (energy, t) = endpoint(&old, &mut lu, dt, 10.0);
    let law = FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: 1.0,
    };
    let accepted = law.initial_state();
    let resolve = |_: &DepositionSimulation<'_>, temperature| {
        let friction =
            tangential_heat_candidate(&law, &accepted, input(dt), temperature, 1e-12).unwrap();
        Ok::<_, DepositionError>(source(friction.dissipated_j, [2.0, 9.0], ()))
    };
    let step = FrictionHeatStep::new(&old, dt, &route, || Ok(())).unwrap();
    let (plain, (), pr) = step
        .evaluate(&energy, t, &resolve, || Ok(()))
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    let (next, (), report) = step
        .evaluate_with_entropy(&energy, t, &resolve, || Ok(()))
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    assert!(pr.entropy.is_none());
    assert_eq!(next.cells(), plain.cells());
    assert_eq!(report.body_heat_j, pr.body_heat_j);
    assert_eq!(report.thermal, pr.thermal);
    let e = report.entropy.unwrap();
    let ta = next.temperature_k(CellId(1)).unwrap();
    let tb = next.temperature_k(CellId(2)).unwrap();
    close(
        e.thermal_storage_change_j_k,
        2.0 * (ta / 300.0).ln() + 6.0 * (tb / 300.0).ln(),
        1e-12,
    );
    close(e.friction_generation_j_k, 2.0 / t, 1e-15);
    close(
        e.interface_conduction_generation_j_k,
        report.body_heat_j[0] * (1.0 / ta - 1.0 / t) + report.body_heat_j[1] * (1.0 / tb - 1.0 / t),
        1e-15,
    );
    assert!(e.interface_conduction_generation_j_k > 0.0 && e.temporal_defect_j_k > 0.0);
    assert_eq!(e.network_conduction_generation_j_k, 0.0);
    assert_eq!(e.reservoir_entropy_input_j_k, 0.0);
    assert_eq!(e.unclassified_source_at_cell_j_k, 0.0);
    assert!(!e.has_unclassified_sources);
    balance(&e);
}

#[test]
fn g1_zero_slip_transfers_heat_without_inventing_friction_entropy() {
    let old = network(1, [1.0, 3.0], [2.0, 9.0], [2.0; 2], [0.0; 2], [10.0, 0.0]);
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [2.0, 9.0]).unwrap();
    let dt = 0.2;
    let mut lu = factor(&old, dt, [4.0, 6.0], 0.0);
    let (energy, t) = endpoint(&old, &mut lu, dt, 0.0);
    let step = FrictionHeatStep::new(&old, dt, &route, || Ok(())).unwrap();
    let (next, (), r) = step
        .evaluate_with_entropy(
            &energy,
            t,
            |_, _| Ok::<_, DepositionError>(source(0.0, [2.0, 9.0], ())),
            || Ok(()),
        )
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    let e = r.entropy.unwrap();
    assert!(r.body_heat_j[0] < 0.0 && r.body_heat_j[1] > 0.0);
    assert_eq!(e.friction_generation_j_k, 0.0);
    assert!(e.interface_conduction_generation_j_k > 0.0 && e.thermal_storage_change_j_k > 0.0);
    assert_eq!(e.network_conduction_generation_j_k, 0.0);
    assert!(!e.has_unclassified_sources);
    close(next.energy_j(), old.energy_j(), 1e-12);
    balance(&e);
}

#[test]
fn g1_other_source_entropy_stays_separate_from_known_friction_heat() {
    let old = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [1.0; 2]).unwrap();
    let step = FrictionHeatStep::new(&old, 1.0, &route, || Ok(())).unwrap();
    // Interface 303 K supplies [4,2] J. A distinct source removes 3 J from A.
    let (_, (), r) = step
        .evaluate_with_entropy(
            &[1.0, 2.0],
            303.0,
            |_, _| {
                let mut input = source(6.0, [1.0; 2], ());
                input.transport.sources_w.insert(CellId(1), -3.0);
                Ok::<_, DepositionError>(input)
            },
            || Ok(()),
        )
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    let e = r.entropy.unwrap();
    assert_eq!(r.body_heat_j, [4.0, 2.0]);
    close(r.external_heat_j, -3.0, 1e-15);
    close(e.unclassified_source_at_cell_j_k, -3.0 / 301.0, 1e-15);
    close(e.friction_generation_j_k, 6.0 / 303.0, 1e-15);
    assert!(e.has_unclassified_sources);
    assert_eq!(e.reservoir_entropy_input_j_k, 0.0);
    balance(&e);
}

#[test]
fn g1_latent_storage_uses_body_temperature_and_interface_generation_once() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if (0.0..=100.0).contains(&h) {
                Ok(300.0)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.0
        }
        fn specific_entropy_change_j_kg_k(
            &self,
            a: f64,
            b: f64,
        ) -> Result<Option<f64>, DepositionError> {
            self.temperature_k(a)?;
            self.temperature_k(b)?;
            Ok(Some((b - a) / 300.0))
        }
    }
    static LATENT: Plateau = Plateau;
    let old = network_with_law(
        1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [10.0; 2], &LATENT,
    );
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [1.0; 2]).unwrap();
    let step = FrictionHeatStep::new(&old, 1.0, &route, || Ok(())).unwrap();
    let (next, (), r) = step
        .evaluate_with_entropy(
            &[6.0; 2],
            303.0,
            |_, _| Ok::<_, DepositionError>(source(12.0, [1.0; 2], ())),
            || Ok(()),
        )
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    let e = r.entropy.unwrap();
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(next.temperature_k(CellId(2)).unwrap(), 300.0);
    close(e.thermal_storage_change_j_k, 12.0 / 300.0, 1e-15);
    close(e.friction_generation_j_k, 12.0 / 303.0, 1e-15);
    close(
        e.interface_conduction_generation_j_k,
        12.0 * (1.0 / 300.0 - 1.0 / 303.0),
        1e-15,
    );
    close(e.temporal_defect_j_k, 0.0, 1e-15);
    assert!(!e.has_unclassified_sources);
    balance(&e);
}

#[test]
fn g0_missing_entropy_and_cancellation_do_not_publish_partial_history() {
    struct Missing;
    impl ThermalLaw for Missing {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            SENSIBLE.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
    }
    static MISSING: Missing = Missing;
    let old = network_with_law(
        1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2], &MISSING,
    );
    let original = old.cells().clone();
    let route = FrictionHeatRoute::from_link(&old, CONTACT, [1.0; 2]).unwrap();
    let step = FrictionHeatStep::new(&old, 1.0, &route, || Ok(())).unwrap();
    let resolve =
        |_: &DepositionSimulation<'_>, _| Ok::<_, DepositionError>(source(0.0, [1.0; 2], ()));
    assert!(step.evaluate(&[0.0; 2], 300.0, &resolve, || Ok(())).is_ok());
    assert!(matches!(
        step.evaluate_with_entropy(&[0.0; 2], 300.0, &resolve, || Ok(())),
        Err(DepositionError::InvalidInput(_))
    ));
    assert!(matches!(
        step.evaluate_with_entropy(&[0.0; 2], 300.0, &resolve, || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(old.cells(), &original);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn g3_friction_numerical_entropy_reduces_under_time_refinement() {
    let defects: Vec<f64> = [10, 20, 40]
        .into_iter()
        .map(|count| {
            let mut thermal = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
            let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0; 2]).unwrap();
            let dt = 1.0 / f64::from(count);
            let (mut defect, mut storage) = (0.0, 0.0);
            for _ in 0..count {
                let t = thermal.temperature_k(CellId(1)).unwrap() + dt + 0.5;
                let step =
                    FrictionHeatStep::new(&thermal, thermal.time_s() + dt, &route, || Ok(()))
                        .unwrap();
                let (next, (), r) = step
                    .evaluate_with_entropy(
                        &[dt; 2],
                        t,
                        |_, _| Ok::<_, DepositionError>(source(2.0 * dt, [1.0; 2], ())),
                        || Ok(()),
                    )
                    .unwrap()
                    .accept(tolerances(), || Ok(()))
                    .unwrap();
                let e = r.entropy.unwrap();
                assert!(
                    e.friction_generation_j_k > 0.0 && e.interface_conduction_generation_j_k > 0.0
                );
                defect += e.temporal_defect_j_k;
                storage += e.thermal_storage_change_j_k;
                thermal = next;
            }
            close(storage, 2.0 * (301.0_f64 / 300.0).ln(), 2e-12);
            defect
        })
        .collect();
    for pair in defects.windows(2) {
        assert!((0.45..0.55).contains(&(pair[1] / pair[0])), "{defects:?}");
    }
}
