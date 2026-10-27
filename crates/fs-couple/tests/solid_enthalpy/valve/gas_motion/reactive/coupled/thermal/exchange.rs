//! G1/G3/G4: finite plumbing feeds the same thermal/chemical/elastic endpoint.
use super::*;
use fs_couple::cavity_phs::modal::exchange::{ModalGasNetwork, ModalGasNetworkError};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::SnapshotNozzleLink;
use std::collections::BTreeMap;

fn reservoir(reaction: &GasReaction, a: f64, t: f64, p: f64) -> ClosedGasParcel {
    let model = reaction.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction,
        vec![a, 0.0, 1.0 - a]).unwrap()).unwrap();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new().with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap().with_quantity("pressure", quantity(QuantityKind::Pressure), p).unwrap();
    ClosedGasParcel::new(model, point, VOLUME).unwrap()
}
fn links() -> [SnapshotNozzleLink; 2] {
    [(5, [20, 10]), (9, [10, 30])].map(|(id, parcels)| SnapshotNozzleLink {
        id, parcels, nozzle: GasNozzle { area_m2: 2e-6, discharge_coefficient: 0.8 },
        nozzle_solve: chemistry().gas,
    })
}
fn network<'a>(reservoirs: &'a BTreeMap<u64, ClosedGasParcel>, links: &'a [SnapshotNozzleLink],
    endpoints: &'a BTreeMap<u64, GasEnergySolve>) -> ModalGasNetwork<'a> {
    ModalGasNetwork { chamber_id: 10, reservoirs, links, endpoints,
        exchange: GasSnapshotSolve { maximum_withdrawn_fraction: 0.1,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }, total_energy_tolerance_j: 1e-7 }
}

fn moved_moles(t0: f64, p0: f64, back: f64) -> f64 {
    let pe = back.max(p0 * (2.0_f64 / 2.4).powf(3.5));
    let te = t0 * (pe / p0).powf(2.0 / 7.0);
    let speed = (7.0 * R * (t0 - te) / 0.028).sqrt();
    0.8 * 2e-6 * pe / (R * te) * speed / RATE as f64
}

// Independent analytic constant-cp nozzle and BE isomer rate, six-tet modal
// response, and caloric/conduction elimination. No native endpoint is sampled.
fn reference(xin: f64) -> (f64, f64, f64, f64, f64) {
    let ni = moved_moles(600.0, 180_000.0, EXTERIOR);
    let no = moved_moles(TEMPERATURE, EXTERIOR, 60_000.0);
    let n = amount() + ni - no;
    let a = 0.5 * (amount() - no) + xin * ni;
    let ustar = amount() * 2.5 * R * TEMPERATURE + ni * 3.5 * R * 600.0 - no * 3.5 * R * TEMPERATURE;
    let stiffness = 2000.0 / 3.0;
    let omega = (stiffness / 0.6_f64).sqrt();
    let compliance = 2.0 * COUPLED_AREA.powi(2) * (omega / (2.0 * RATE as f64)).sin().powi(2) / stiffness;
    let gas = |q: f64| {
        let sample = |p: f64| {
            let dv = compliance * (p - OUTSIDE);
            let t = p * (VOLUME + dv) / (n * R);
            let beta = 3000.0 * R * t / (RATE as f64 * EXTERIOR);
            let xi = beta * a / (1.0 + beta * (1.0 + (-RELEASE / (R * t)).exp()));
            let residual = n * 2.5 * R * t - RELEASE * xi - ustar + p * dv - GAS_HEAT - q;
            (residual, t, xi)
        };
        let p = bracketed_root(0.5 * EXTERIOR, 1.5 * EXTERIOR, |p| sample(p).0);
        let (_, t, xi) = sample(p);
        (p, t, xi)
    };
    let q = bracketed_root(-2.0, 2.0, |q| {
        let (_, tg, _) = gas(q);
        let tb = BODY_T0 + (BODY_HEAT - q) / (BODY_MASS * 50.0);
        let kb = 2.0 * (1.0 + (tb - 200.0) / 600.0);
        q - 0.2 / (0.001 / 0.03 + 0.02 / kb) * (tb - tg) / RATE as f64
    });
    let (p, t, xi) = gas(q);
    (q, p, t, xi, BODY_T0 + (BODY_HEAT - q) / (BODY_MASS * 50.0))
}

#[test]
fn intake_composition_changes_reaction_heat_and_mesh_motion_with_one_inventory_debit() {
    with_cx(|cx| {
        let mut extents = Vec::new();
        for xin in [0.1, 0.9] {
            let (reaction, gas) = reacting_gas(3000.0, [0.5, 0.0, 0.5]);
            let reaction = implicit(reaction);
            let run = GasValveMotion::new(gas,
                motion(1000.0, 0.0, 0.0, COUPLED_AREA, RATE, cx), cx).unwrap();
            let reservoirs = BTreeMap::from([(20, reservoir(&reaction, xin, 600.0, 180_000.0)),
                (30, reservoir(&reaction, 0.5, 400.0, 60_000.0))]);
            let before = reservoirs.clone();
            let endpoints = BTreeMap::from([(10, chemistry().gas), (20, chemistry().gas), (30, chemistry().gas)]);
            let links = links();
            let (card, key, caloric) = thermal_material(2.0, 50.0);
            let point = thermal_point(BODY_T0);
            let input = body(caloric.curve(), 15_000.0);
            let next = run.reactive_thermal_network_layers_candidate(network(&reservoirs, &links, &endpoints),
                &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(), layers(&card, &key, &point), cx).unwrap();
            let (q, p, tg, xi, tb) = reference(xin);
            let report = next.step();
            close(report.evolution.exchanged_heat_j, q, 2e-9);
            close(next.motion().gas().state().state().pressure, p, 2e-5);
            close(next.motion().gas().state().state().temperature, tg, 4e-7);
            close(next.thermal().state.temperature_k(), tb, 1e-7);
            close(report.evolution.coupled.reaction.reaction.extent_mol, xi, 3e-13);
            close(report.flow.links[0].transferred_mass_kg, 0.028 * moved_moles(600.0, 180_000.0, EXTERIOR), 1e-14);
            close(reservoirs[&20].mass_kg() - next.reservoirs()[&20].mass_kg(), report.flow.links[0].transferred_mass_kg, 1e-14);
            close(report.total_mass_residual_kg, 0.0, 1e-14);
            close(report.total_energy_residual_j, 0.0, 1e-7);
            close(report.supplied_heat_j, GAS_HEAT + BODY_HEAT, 1e-15);
            let stiffness = 2000.0 / 3.0;
            let omega = (stiffness / 0.6_f64).sqrt();
            let velocity = COUPLED_AREA * (p - OUTSIDE) * (omega / RATE as f64).sin() / (0.6 * omega);
            close(next.motion().motion().nodal_motion(cx).unwrap().velocity_m_s[7][0], velocity, 1e-9);
            assert_eq!(reservoirs, before);
            assert_eq!(run.motion().accepted_steps(), 0);
            assert_eq!(next.motion().motion().accepted_steps(), 1);
            assert!(report.evolution.exchange_trials > 1);
            extents.push(xi);
        }
        assert!((extents[0] - extents[1]).abs() > 1e-7);
    });
}

#[test]
fn late_source_failure_cancels_every_gas_body_and_modal_candidate_then_retry_is_identical() {
    with_cx(|cx| {
        let (reaction, gas) = reacting_gas(3000.0, [0.5, 0.0, 0.5]);
        let reaction = implicit(reaction);
        let run = GasValveMotion::new(gas,
            motion(1000.0, 0.0, 0.0, COUPLED_AREA, RATE, cx), cx).unwrap();
        let initial = run.clone();
        let reservoirs = BTreeMap::from([(20, reservoir(&reaction, 0.9, 600.0, 180_000.0)),
            (30, reservoir(&reaction, 0.5, 400.0, 60_000.0))]);
        let original = reservoirs.clone();
        let endpoints = BTreeMap::from([(10, chemistry().gas), (20, chemistry().gas), (30, chemistry().gas)]);
        let links = links();
        let (card, key, caloric) = thermal_material(2.0, 50.0);
        let point = thermal_point(BODY_T0);
        let input = body(caloric.curve(), 15_000.0);
        let mut count = 0;
        let failed = run.reactive_thermal_network_candidate(network(&reservoirs, &links, &endpoints),
            &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(), |g, b| {
                count += 1;
                if count == 2 { return Err("late source refusal"); }
                let source = layers(&card, &key, &point).resolve(g, b).unwrap();
                Ok(ConductanceSample { conductance_w_k: source.conductance_w_k, source: () })
            }, cx);
        assert!(matches!(failed, Err(ModalGasNetworkError::Evolution(
            ReactiveModalThermalError::Interface("late source refusal")))));
        assert_eq!(count, 2);
        assert_eq!(reservoirs, original);
        assert_eq!(run.gas(), initial.gas());
        assert_eq!(run.motion().model().states(), initial.motion().model().states());
        let retry = || run.reactive_thermal_network_layers_candidate(network(&reservoirs, &links, &endpoints),
            &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(), layers(&card, &key, &point), cx).unwrap();
        let a = retry(); let b = retry();
        assert_eq!(a.motion().gas(), b.motion().gas());
        assert_eq!(a.reservoirs(), b.reservoirs());
        assert_eq!(a.motion().motion().model().states(), b.motion().motion().model().states());
        assert_eq!(a.thermal(), b.thermal());
    });
}

#[test]
fn cancellation_after_flow_and_actual_material_probe_publishes_no_endpoints() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 4, kernel_id: 18, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let (reaction, gas) = reacting_gas(3000.0, [0.5, 0.0, 0.5]);
        let reaction = implicit(reaction);
        let run = GasValveMotion::new(gas,
            motion(1000.0, 0.0, 0.0, COUPLED_AREA, RATE, &cx), &cx).unwrap();
        let reservoirs = BTreeMap::from([(20, reservoir(&reaction, 0.9, 600.0, 180_000.0))]);
        let original = reservoirs.clone();
        let endpoints = BTreeMap::from([(10, chemistry().gas), (20, chemistry().gas)]);
        let inlet = [links()[0]];
        let (card, key, caloric) = thermal_material(2.0, 50.0);
        let point = thermal_point(BODY_T0);
        let input = body(caloric.curve(), 15_000.0);
        let mut probed = false;
        let result = run.reactive_thermal_network_candidate(network(&reservoirs, &inlet, &endpoints),
            &reaction, GAS_HEAT, OUTSIDE, input, chemistry(), thermal_controls(), |g, b| {
                let sample = layers(&card, &key, &point).resolve(g, b)?;
                probed = true;
                gate.request();
                Ok::<_, LayerConductionError>(sample)
            }, &cx);
        assert!(probed);
        assert!(matches!(result, Err(ModalGasNetworkError::Cancelled)));
        assert_eq!(reservoirs, original);
        assert_eq!(run.motion().accepted_steps(), 0);
    });
}
