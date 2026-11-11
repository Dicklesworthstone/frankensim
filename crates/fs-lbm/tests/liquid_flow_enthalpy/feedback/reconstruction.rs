//! Live momentum changes use the original limited enthalpy traces and storage.

use super::*;
use fs_conduction::transient::deposition::advection::reconstruction::{
    AxialCellStencil, LimitedEnthalpyReconstruction,
};
use std::sync::Arc;

#[cfg(feature = "implicit-enthalpy-coupling")]
#[path = "reconstruction/implicit.rs"]
mod implicit;

fn plan(thermal: &DepositionSimulation<'_>) -> Arc<LimitedEnthalpyReconstruction> {
    Arc::new(
        LimitedEnthalpyReconstruction::try_new(
            thermal,
            &[
                AxialCellStencil {
                    previous: CellId(0),
                    cell: CellId(1),
                    next: CellId(2),
                },
                AxialCellStencil {
                    previous: CellId(1),
                    cell: CellId(2),
                    next: CellId(3),
                },
            ],
            2,
            || Ok(()),
        )
        .unwrap(),
    )
}

fn profile<'a>(
    domain: &VoxelDomain,
    rho: f64,
    caloric: &'a ResolvedHeatCapacityCurve,
    h: [f64; 4],
) -> DepositionSimulation<'a> {
    let original = initial(domain, rho, caloric);
    let events = original
        .cells()
        .values()
        .map(|cell| {
            let mut cell = cell.clone();
            cell.specific_h_j_kg = h[cell.id.0 as usize];
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MAT, caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 64,
            max_events: 4,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

// Independent uniform-cell minmod face values, including the original external
// donor traces. Source viscosity and density still determine the actual flux.
fn traces(h: [f64; 4], forward: bool, inlet_h: f64) -> [f64; 5] {
    let mut face = [0.0; 5];
    face[0] = if forward { inlet_h } else { h[0] };
    face[4] = if forward { h[3] } else { inlet_h };
    for (i, trace) in face.iter_mut().enumerate().take(4).skip(1) {
        let c = if forward { i - 1 } else { i };
        let slope = if c == 0 || c == 3 {
            0.0
        } else {
            let left = h[c] - h[c - 1];
            let right = h[c + 1] - h[c];
            if left * right > 0.0 {
                left.signum() * left.abs().min(right.abs())
            } else {
                0.0
            }
        };
        *trace = h[c] + if forward { 0.5 * slope } else { -0.5 * slope };
    }
    face
}

#[test]
fn live_reconstruction_tracks_source_viscosity_and_reversed_native_flow() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = feedback_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = profile(&domain, 1.0, &caloric, [20.0, 40.0, 80.0, 140.0]);
    let plan = plan(&old);
    let gate = CancelGate::new();
    for forward in [true, false] {
        let mut flow = porous_config();
        flow.internal_fans[0].blows_positive = forward;
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &flow, Heater::default(), &gate, 100)
            .unwrap()
            .with_reconstruction(Arc::clone(&plan))
            .unwrap();
        let (next, accepted, energy) = old
            .advance_observed(&observer, 0.5, config(), || Ok(()))
            .unwrap();
        let mut h = [20.0, 40.0, 80.0, 140.0];
        let mut stream_j = 0.0;
        let mut mass_kg = 0.0;
        let mut maximum_courant: f64 = 0.0;
        let mut last_q = 0.0;
        let mut changed_q = false;
        for &dt in &accepted.inner().durations {
            let q = if forward { 4.0 } else { -4.0 }
                / h.iter().map(|&h| viscosity(1.0, h)).sum::<f64>();
            changed_q |= last_q != 0.0 && (q - last_q).abs() > 1e-5;
            let face = traces(h, forward, 0.0);
            stream_j += dt * q * (face[0] - face[4]);
            mass_kg += dt * q.abs();
            maximum_courant = maximum_courant.max(dt * q.abs());
            for c in 0..4 {
                h[c] += dt * (if c == 1 { 100.0 } else { 0.0 } + q * (face[c] - face[c + 1]));
            }
            last_q = q;
        }
        assert!(changed_q);
        assert!(std::ptr::eq(
            accepted.reconstruction().unwrap(),
            plan.as_ref()
        ));
        assert_eq!(accepted.momentum_solves(), accepted.inner().durations.len());
        near(
            accepted.applied().unwrap().flow.field.flux_x(2, 0, 0),
            last_q,
        );
        for sample in accepted.endpoint() {
            near(next.cells()[&sample.cell].specific_h_j_kg, h[sample.voxel]);
            near(sample.temperature_k, 300.0 + h[sample.voxel] / 2.0);
            near(
                sample.transport.dynamic_viscosity_pa_s(),
                viscosity(1.0, h[sample.voxel]),
            );
            let receipt = &sample
                .material
                .material()
                .property(DYNAMIC_VISCOSITY_PROPERTY)
                .unwrap()
                .answer()
                .receipt;
            assert_eq!(
                receipt.selected,
                liquid
                    .reference_state()
                    .material()
                    .property(DYNAMIC_VISCOSITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
            assert!(receipt.query_point.contains(&("lot".into(), 7.0)));
        }
        near(accepted.audit().external_enthalpy_j, stream_j);
        near(accepted.audit().inlet_mass_kg, mass_kg);
        near(accepted.audit().outlet_mass_kg, mass_kg);
        near(accepted.audit().maximum_courant, maximum_courant);
        near(
            next.energy_j() - old.energy_j(),
            accepted.inner().supplied_j + stream_j,
        );
        near(
            energy.external_heat_j,
            accepted.inner().supplied_j + stream_j,
        );
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(observer.momentum_solves(), 0);
        assert!(observer.applied().is_none());
        assert_eq!(observer.inner(), &Heater::default());
    }
}

#[derive(Clone, Default)]
struct ReversingHeat {
    intervals: Vec<(f64, f64)>,
}

impl ReversingHeat {
    fn watts(time_s: f64) -> f64 {
        if time_s < 0.25 { 20.0 } else { -40.0 }
    }
}

impl DepositionObserver for ReversingHeat {
    type Error = DepositionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.intervals.clear();
        Ok(())
    }
    fn next_transition_time_s(
        &mut self,
        state: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        Ok((state.time_s() < 0.25).then_some(0.25))
    }
    fn sources(
        &mut self,
        state: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok((0..4)
            .map(|c| (CellId(c), Self::watts(state.time_s())))
            .collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        start: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.intervals.push((start, dt));
        Ok(())
    }
}

#[test]
fn reconstructed_buoyant_flow_retains_geometric_neighbors_at_rest_and_through_reversal() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let mut flow = porous_config();
    flow.internal_fans.clear();
    let card = buoyant_card(0.001);
    let liquid = LiquidTemperatureLaw::try_new(
        &card,
        &point(340.0, 100000.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = liquid_caloric(&card);
    let rho = liquid.reference_state().density_kg_m3();
    let old = profile(&domain, rho, &caloric, [80.0; 4]);
    let plan = plan(&old);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let mut binding = input(&ids, &liquid, &caloric);
    binding.inlet_specific_h_j_kg = [Some(80.0); 6];
    let mut observer = binding
        .feedback_boussinesq(
            &old,
            &domain,
            &flow,
            ReversingHeat::default(),
            &gate,
            100,
            buoyancy(0.1),
        )
        .unwrap()
        .with_reconstruction(Arc::clone(&plan))
        .unwrap();
    let mut current = old.clone();
    let mut h = [80.0; 4];
    let mut saw_forward = false;
    let mut saw_reverse = false;
    for target in [0.125, 0.25, 0.75] {
        let (next, accepted, energy) = current
            .advance_observed(&observer, target, config(), || Ok(()))
            .unwrap();
        let mut external_j = 0.0;
        let mut last_q = 0.0;
        for &(start, dt) in &accepted.inner().intervals {
            // Native series Darcy balance after subtracting hydrostatic rho0.
            let q = h
                .iter()
                .map(|h| (1.0 - 0.001 * h / 2.0 - rho) * -9.81)
                .sum::<f64>()
                / h.iter().map(|&h| viscosity(1.0, h)).sum::<f64>();
            saw_forward |= q > 1e-5;
            saw_reverse |= q < -1e-5;
            let face = traces(h, q >= 0.0, 80.0);
            let watts = ReversingHeat::watts(start);
            external_j += dt * (4.0 * watts + rho * q * (face[0] - face[4]));
            for c in 0..4 {
                h[c] += dt * (watts / rho + q * (face[c] - face[c + 1]));
            }
            last_q = q;
        }
        if target == 0.125 {
            near(last_q, 0.0);
            assert_eq!(accepted.momentum_solves(), 1);
        }
        near(
            accepted.applied().unwrap().flow.field.flux_x(2, 0, 0),
            last_q,
        );
        for (c, expected) in h.iter().enumerate() {
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, *expected);
        }
        near(next.energy_j() - current.energy_j(), external_j);
        near(energy.external_heat_j, external_j);
        near(next.mass_kg(), old.mass_kg());
        assert!(std::ptr::eq(
            accepted.reconstruction().unwrap(),
            plan.as_ref()
        ));
        current = next;
        observer = accepted;
    }
    assert!(saw_forward && saw_reverse);
    assert!(observer.applied().unwrap().flow.field.flux_x(2, 0, 0) < 0.0);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn live_empty_reconstruction_preserves_donor_arithmetic_and_late_rejection_is_atomic() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let flow = porous_config();
    let card = feedback_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let donor = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow, Heater::default(), &gate, 100)
        .unwrap();
    let empty = LimitedEnthalpyReconstruction::try_new(&old, &[], 0, || Ok(())).unwrap();
    let wrapped = donor.clone().with_reconstruction(empty).unwrap();
    let (a, a_observer, a_audit) = old
        .advance_observed(&donor, 0.5, config(), || Ok(()))
        .unwrap();
    let (b, b_observer, b_audit) = old
        .advance_observed(&wrapped, 0.5, config(), || Ok(()))
        .unwrap();
    assert_eq!(a.cells(), b.cells());
    assert_eq!(a_audit, b_audit);
    assert_eq!(a_observer.audit(), b_observer.audit());
    assert_eq!(a_observer.inner(), b_observer.inner());
    assert!(b_observer.reconstruction().is_none());
    assert!(a_observer.with_reconstruction(plan(&old)).is_err());

    let history = Heater {
        reject_step: Some(2),
        ..Heater::default()
    };
    let refused = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow, history.clone(), &gate, 100)
        .unwrap()
        .with_reconstruction(plan(&old))
        .unwrap();
    assert!(
        old.advance_observed(&refused, 0.5, config(), || Ok(()))
            .is_err()
    );
    assert_eq!(refused.inner(), &history);
    assert_eq!(refused.audit().substeps, 0);
    assert_eq!(refused.momentum_solves(), 0);
    assert!(refused.applied().is_none());
    assert_eq!(old.time_s(), 0.0);
    near(old.energy_j(), 0.0);
}

#[test]
fn stationary_wide_channel_retains_only_the_selected_axial_neighbor_connections() {
    let domain = VoxelDomain::new(4, 2, 1, 1.0).unwrap();
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..8).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    for budget in [2, 3] {
        let mut binding = input(&ids, &liquid, &caloric);
        binding.max_connections = budget;
        let observer = binding
            .feedback(&old, &domain, &flow, (), &gate, 4)
            .unwrap()
            .with_reconstruction(plan(&old))
            .unwrap();
        let result = old.advance_observed(&observer, 0.1, config(), || Ok(()));
        if budget == 2 {
            assert!(matches!(
                result,
                Err(MaterialFlowError::Thermal(DepositionError::Budget(
                    "flow/enthalpy connections"
                )))
            ));
        } else {
            let (next, accepted, audit) = result.unwrap();
            assert_eq!(next.cells(), old.cells());
            assert_eq!(accepted.momentum_solves(), 1);
            assert_eq!(accepted.audit().transported_mass_kg, 0.0);
            assert_eq!(audit.external_heat_j, 0.0);
        }
        assert!(observer.applied().is_none());
    }
}
