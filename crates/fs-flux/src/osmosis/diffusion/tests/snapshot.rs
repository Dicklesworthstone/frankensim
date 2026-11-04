//! Physical readouts are solved at the accepted state without consuming time.
use super::electroneutral::loop_network;
use super::*;

#[test]
fn accepted_snapshot_replays_current_native_fluxes_and_exposes_loop_potentials() {
    let mut system = loop_network(ElectroneutralDiffusionControls::default());
    system.advance(2.0, numerical(1e-5), || false).unwrap();
    let before = bits(&system);
    let work_before = system.electroneutral_diffusion.as_ref().unwrap().report();
    let snapshot = system.ionic_diffusion_snapshot(1 << 20, || false).unwrap();
    assert_eq!(bits(&system), before);
    assert_eq!(
        system.electroneutral_diffusion.as_ref().unwrap().report(),
        work_before
    );
    assert_eq!(snapshot.time_s, system.time_s());
    assert_eq!(snapshot.temperature_k, T);
    assert_eq!(snapshot.species_count, 2);
    assert_eq!(snapshot.connections().len(), 3);
    assert_eq!(snapshot.potential_gauges(), &[0, 0, 0]);
    assert_eq!(snapshot.dimensionless_potentials()[0], 0.0);
    assert_eq!(snapshot.work.potential_solves, 1);
    assert_eq!(snapshot.work.face_evaluations, 6);
    assert!(snapshot.solute_fluxes_mol_m2_s(3).is_none());
    let mut currents = [0.0; 3];
    for (index, connection) in snapshot.connections().iter().enumerate() {
        let link = &system.solution_diffusion_links()[connection.link_index];
        assert_eq!(connection.endpoints, link.endpoints());
        let [left, right] = link.endpoints();
        let difference =
            snapshot.dimensionless_potentials()[right] - snapshot.dimensionless_potentials()[left];
        assert_eq!(
            connection
                .transport
                .dimensionless_potential_right_minus_left,
            difference
        );
        let law = link.material();
        let mut actual = [0.0; 2];
        let direct = law
            .flux_at_potential_into(
                system.compartments()[left].solution(),
                system.compartments()[right].solution(),
                connection.distance_m,
                physical().constitutive,
                link.zero_current_diffusion().unwrap(),
                difference,
                &mut law.workspace().unwrap(),
                &mut actual,
                || false,
            )
            .unwrap();
        let observed = snapshot.solute_fluxes_mol_m2_s(index).unwrap();
        for i in 0..2 {
            assert!((actual[i] - observed[i]).abs() < 1e-18);
        }
        assert!(
            (direct.charge_flux_mol_m2_s - connection.transport.charge_flux_mol_m2_s).abs() < 1e-18
        );
        assert!(
            connection
                .transport
                .electrochemical_entropy_production_w_m3_k
                .unwrap()
                >= -1e-9
        );
        let charge = connection.area_m2 * (observed[0] - observed[1]);
        currents[left] += charge;
        currents[right] -= charge;
    }
    for (actual, observed) in currents.iter().zip(snapshot.node_charge_rates_mol_s()) {
        assert!((actual - observed).abs() < 1e-22);
        assert!(observed.abs() <= 1e-16);
    }
    assert!(snapshot.work.maximum_edge_charge_rate_mol_s > 1e-10);
}

#[test]
fn snapshot_work_budget_and_cancellation_do_not_change_accepted_evolution() {
    let mut system = loop_network(ElectroneutralDiffusionControls::default());
    let mut reference = system.clone();
    let before = bits(&system);
    assert!(matches!(
        system.ionic_diffusion_snapshot(0, || false),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    assert_eq!(bits(&system), before);
    let mut polls = 0;
    system
        .ionic_diffusion_snapshot(1 << 20, || {
            polls += 1;
            false
        })
        .unwrap();
    let work = system.electroneutral_diffusion.as_ref().unwrap().report();
    let mut seen = 0;
    assert!(matches!(
        system.ionic_diffusion_snapshot(1 << 20, || {
            seen += 1;
            seen == polls - 1
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&system), before);
    assert_eq!(
        system.electroneutral_diffusion.as_ref().unwrap().report(),
        work
    );
    let a = system.advance(1.0, numerical(1e-5), || false).unwrap();
    let b = reference.advance(1.0, numerical(1e-5), || false).unwrap();
    assert_eq!(bits(&system), bits(&reference));
    assert_eq!(a.electroneutral_diffusion, b.electroneutral_diffusion);

    let mut controls = ElectroneutralDiffusionControls::default();
    controls.maximum_face_evaluations = 5;
    let mut bounded = loop_network(controls);
    let saved = bits(&bounded);
    assert!(matches!(
        bounded.ionic_diffusion_snapshot(1 << 20, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&bounded), saved);
    assert_eq!(
        bounded
            .electroneutral_diffusion
            .as_ref()
            .unwrap()
            .report()
            .face_evaluations,
        0
    );
}
