//! G1/G3: original gas plumbing drives the mesh-derived mode without resets.
use super::*;
use fs_couple::cavity_phs::modal::exchange::{ModalGasNetwork, ModalGasNetworkError};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::SnapshotNozzleLink;
use std::collections::BTreeMap;

const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

fn reservoir(x: f64, pressure: f64) -> ClosedGasParcel {
    let gas = gas(x);
    let point = gas.point().clone().with_quantity("pressure",
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
        pressure).unwrap();
    ClosedGasParcel::new(gas.model().clone(), point, VOLUME).unwrap()
}
fn links() -> [SnapshotNozzleLink; 2] {
    [(5, [20, 10]), (9, [10, 30])].map(|(id, parcels)| SnapshotNozzleLink {
        id, parcels, nozzle: GasNozzle { area_m2: 2e-6, discharge_coefficient: 0.8 },
        nozzle_solve: controls().gas,
    })
}
fn network<'a>(reservoirs: &'a BTreeMap<u64, ClosedGasParcel>, links: &'a [SnapshotNozzleLink],
    endpoints: &'a BTreeMap<u64, GasEnergySolve>) -> ModalGasNetwork<'a> {
    ModalGasNetwork { chamber_id: 10, reservoirs, links, endpoints,
        exchange: GasSnapshotSolve { maximum_withdrawn_fraction: 0.1,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 },
        total_energy_tolerance_j: 1e-7 }
}
fn molar_mass(x: f64) -> f64 { x * 0.028 + (1.0 - x) * 0.004 }

// Independent ideal isentropic nozzle, constant molar heat capacities. This
// does not query the native throat, energy inversion or transport reports.
fn moved_moles(x: f64, p0: f64, pb: f64) -> f64 {
    let cp = (2.5 + x) * R;
    let gamma = (2.5 + x) / (1.5 + x);
    let pe = pb.max(p0 * (2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0)));
    let te = 400.0 * (pe / p0).powf((gamma - 1.0) / gamma);
    let velocity = (2.0 * cp * (400.0 - te) / molar_mass(x)).sqrt();
    0.8 * 2e-6 * pe / (R * te) * velocity / RATE as f64
}

#[test]
fn simultaneous_intake_exhaust_changes_composition_pressure_and_real_wall_motion() {
    with_cx(|cx| {
        let x0 = 0.9;
        let mut speeds = Vec::new();
        for xin in [0.1, 0.7] {
            let reservoirs = BTreeMap::from([(20, reservoir(xin, 180_000.0)),
                (30, reservoir(0.5, 60_000.0))]);
            let endpoints = BTreeMap::from([(10, controls().gas), (20, controls().gas), (30, controls().gas)]);
            let links = links();
            let run = GasValveMotion::new(gas(x0), motion(1000.0, 0.0, 0.0, AREA, RATE, cx), cx).unwrap();
            let next = run.network_candidate(network(&reservoirs, &links, &endpoints),
                0.0, EXTERIOR, controls(), cx).unwrap();
            let incoming = moved_moles(xin, 180_000.0, EXTERIOR);
            let outgoing = moved_moles(x0, EXTERIOR, 60_000.0);
            let n0 = EXTERIOR * VOLUME / (R * 400.0);
            let na = (n0 - outgoing) * x0 + incoming * xin;
            let nb = (n0 - outgoing) * (1.0 - x0) + incoming * (1.0 - xin);
            let capacity = (2.5 * na + 1.5 * nb) * R;
            let stored = n0 * (1.5 + x0) * R * 400.0
                + incoming * (2.5 + xin) * R * 400.0
                - outgoing * (2.5 + x0) * R * 400.0;
            let tstar = stored / capacity;
            let nr = (na + nb) * R;
            let gamma = 1.0 + nr / capacity;
            let stiffness = 2000.0 / 3.0;
            let omega = (stiffness / 0.6_f64).sqrt();
            let compliance = 2.0 * AREA.powi(2) * (omega / (2.0 * RATE as f64)).sin().powi(2) / stiffness;
            let a = 0.5 * gamma;
            let b = VOLUME + a * compliance * EXTERIOR;
            let drive = compliance * (nr * tstar - EXTERIOR * VOLUME);
            let dv = 2.0 * drive / (b + (b * b + 4.0 * a * drive).sqrt());
            let held = nr * tstar / (VOLUME + a * dv);
            let temperature = tstar - held * dv / capacity;
            let velocity = AREA * (held - EXTERIOR) * (omega / RATE as f64).sin() / (0.6 * omega);
            near(next.step().flow.links[0].transferred_mass_kg, incoming * molar_mass(xin), 1e-14);
            near(next.step().flow.links[1].transferred_mass_kg, outgoing * molar_mass(x0), 1e-14);
            near(next.motion().gas().mass_kg(), na * 0.028 + nb * 0.004, 1e-14);
            near(next.motion().gas().state().state().specific_gas_constant,
                nr / (na * 0.028 + nb * 0.004), 1e-8);
            near(next.motion().gas().state().state().temperature, temperature, 2e-7);
            near(next.motion().gas().volume_m3(), VOLUME + dv, 1e-14);
            near(next.motion().motion().nodal_motion(cx).unwrap().velocity_m_s[7][0], velocity, 1e-9);
            near(next.step().total_energy_residual_j, 0.0, 1e-7);
            near(next.step().total_mass_residual_kg, 0.0, 1e-14);
            assert_eq!(next.step().evolution.gas.supplied_heat_j, 0.0);
            assert_eq!(run.motion().accepted_steps(), 0);
            assert_eq!(next.motion().motion().accepted_steps(), 1);
            assert_eq!(next.motion().motion().valve().identity(), run.motion().valve().identity());
            speeds.push(next.motion().gas().state().state().sound_speed);
            // Renaming reporting IDs and permuting their input order cannot
            // change either inventories or subsequent native structural state.
            let mut renamed = links;
            renamed[0].id = 91; renamed[1].id = 3; renamed.reverse();
            let same = run.network_candidate(network(&reservoirs, &renamed, &endpoints),
                0.0, EXTERIOR, controls(), cx).unwrap();
            assert_eq!(same.motion().gas(), next.motion().gas());
            assert_eq!(same.motion().motion().model().states(), next.motion().motion().model().states());
            assert_eq!(same.reservoirs(), next.reservoirs());
        }
        assert!((speeds[0] - speeds[1]).abs() > 0.01);
    });
}

#[test]
fn closed_flow_is_native_motion_and_continuation_retains_the_anchored_volume() {
    with_cx(|cx| {
        let mut run = GasValveMotion::new(gas(0.9), motion(1000.0, 0.01, 0.05, AREA, RATE, cx), cx).unwrap();
        let identity = run.motion().valve().identity();
        let mut reservoirs = BTreeMap::from([(20, reservoir(0.1, 180_000.0)), (30, reservoir(0.5, 60_000.0))]);
        let endpoints = BTreeMap::from([(10, controls().gas), (20, controls().gas), (30, controls().gas)]);
        let mut closed = links();
        for link in &mut closed { link.nozzle.area_m2 = 0.0; }
        let mut expected = run.clone();
        expected.step(0.2, EXTERIOR, controls(), cx).unwrap();
        let next = run.network_candidate(network(&reservoirs, &closed, &endpoints),
            0.2, EXTERIOR, controls(), cx).unwrap();
        assert_eq!(next.motion().gas(), expected.gas());
        assert_eq!(next.motion().motion().model().states(), expected.motion().model().states());
        assert_eq!(next.reservoirs(), &reservoirs);
        let open = links();
        for tick in 1..=5 {
            let next = run.network_candidate(network(&reservoirs, &open, &endpoints),
                0.0, EXTERIOR, controls(), cx).unwrap();
            (run, reservoirs, _, _) = next.into_parts();
            assert_eq!(run.motion().accepted_steps(), tick);
            assert_eq!(run.motion().valve().identity(), identity);
            let actual = run.motion().nodal_motion(cx).unwrap();
            near(run.gas().volume_m3() - VOLUME, AREA * actual.displacement_m[7][0], 1e-14);
        }
    });
}

#[test]
fn refused_motion_after_real_flow_preserves_reservoirs_and_retry_history() {
    with_cx(|cx| {
        let run = GasValveMotion::new(gas(0.9), motion(1000.0, 0.0, 0.0, AREA, RATE, cx), cx).unwrap();
        let reservoirs = BTreeMap::from([(20, reservoir(0.1, 180_000.0)), (30, reservoir(0.5, 60_000.0))]);
        let before = reservoirs.clone();
        let initial_gas = run.gas().clone();
        let states = run.motion().model().states().to_vec();
        let endpoints = BTreeMap::from([(10, controls().gas), (20, controls().gas), (30, controls().gas)]);
        let links = links();
        let wrong = ModalGasSolve { volume_interval_m3: [0.0011, 0.0012], ..controls() };
        assert!(matches!(run.network_candidate(network(&reservoirs, &links, &endpoints),
            0.0, EXTERIOR, wrong, cx), Err(ModalGasNetworkError::Evolution(ModalGasError::Invalid(_)))));
        assert_eq!(reservoirs, before);
        assert_eq!(run.gas(), &initial_gas);
        assert_eq!(run.motion().model().states(), states);
        let retry = || run.network_candidate(network(&reservoirs, &links, &endpoints),
            0.0, EXTERIOR, controls(), cx).unwrap();
        let first = retry(); let second = retry();
        assert_eq!(first.motion().gas(), second.motion().gas());
        assert_eq!(first.reservoirs(), second.reservoirs());
        assert_eq!(first.motion().motion().model().states(), second.motion().motion().model().states());
    });
}
