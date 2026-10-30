use super::*;
use crate::conjugate::FlowFace;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> { Ok(300.0 + h) }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
    fn specific_entropy_change_j_kg_k(&self, from: f64, to: f64)
        -> Result<Option<f64>, DepositionError> {
        Ok(Some(((300.0 + to) / (300.0 + from)).ln()))
    }
}

fn input<'a>(ids: &'a [Option<CellId>], law: &'a Caloric) -> FlowEnthalpyInput<'a> {
    FlowEnthalpyInput {
        cell_ids: ids, material_id: MaterialId(4), caloric: law,
        inlet_specific_h_j_kg: [Some(10.0), Some(11.0), Some(12.0), Some(13.0), Some(14.0), Some(15.0)],
        max_cells: ids.len(), max_connections: 1000,
    }
}

#[test]
fn graded_three_axis_fluxes_keep_direction_area_density_and_unique_face_identity() {
    let domain = VoxelDomain::graded([
        vec![0.5, 1.5, 2.0], vec![0.25, 0.75], vec![0.5, 1.0],
    ]).unwrap();
    let ids: Vec<_> = (0..domain.cell_count()).map(|c| Some(CellId(c as u64))).collect();
    let velocity = [-0.25, 0.125, -0.5];
    let field = FlowField::from_face_velocity(&domain, [FlowFace::Fixed; 6],
        |axis, _| velocity[axis]).unwrap();
    let (flows, boundaries) = mass_graph(&domain, &field, 8.0, input(&ids, &Caloric), &mut || Ok(())).unwrap();
    let mut identities = BTreeSet::new();
    let mut balance = vec![0.0; domain.cell_count()];
    // Each interior face appears exactly once; positive axis is independent of
    // actual flow direction. All three unequal face areas enter in SI units.
    assert_eq!(flows.len(), 20);
    assert_eq!(boundaries.len(), 32);
    for flow in &flows {
        assert!(identities.insert(flow.id));
        let a = flow.a.0 as usize;
        let b = flow.b.0 as usize;
        let ca = domain.coords(a);
        let cb = domain.coords(b);
        let axis = (0..3).find(|&axis| ca[axis] != cb[axis]).unwrap();
        assert_eq!(flow.kg_s, 8.0 * domain.face_area(a, axis) * velocity[axis]);
        assert_eq!(cb[axis], ca[axis] + 1);
        balance[a] -= flow.kg_s;
        balance[b] += flow.kg_s;
    }
    for boundary in &boundaries {
        match *boundary {
            MassBoundary::Inlet { id, cell, kg_s, specific_h_j_kg, .. } => {
                assert!(identities.insert(id));
                let face = (specific_h_j_kg - 10.0) as usize;
                assert!([1, 2, 5].contains(&face));
                assert_eq!(kg_s, 8.0 * domain.face_area(cell.0 as usize, face / 2) * velocity[face / 2].abs());
                balance[cell.0 as usize] += kg_s;
            }
            MassBoundary::Outlet { id, cell, kg_s } => {
                assert!(identities.insert(id));
                balance[cell.0 as usize] -= kg_s;
            }
        }
    }
    assert!(balance.iter().all(|&rate| rate == 0.0));
}

#[test]
fn closed_circulation_has_no_external_enthalpy_and_no_duplicate_faces() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let field = FlowField::from_face_velocity(&domain, [FlowFace::Wall; 6], |axis, centre| {
        if axis == 0 { if centre[1] < 1.0 { 0.25 } else { -0.25 } }
        else if axis == 1 { if centre[0] < 1.0 { -0.25 } else { 0.25 } }
        else { 0.0 }
    }).unwrap();
    let (flows, boundaries) = mass_graph(&domain, &field, 4.0, input(&ids, &Caloric), &mut || Ok(())).unwrap();
    assert_eq!(flows.len(), 4);
    assert!(boundaries.is_empty());
    assert!(flows.iter().all(|f| f.kg_s.abs() == 1.0));
}

#[test]
fn divergence_missing_backflow_and_connection_budget_refuse_before_history_creation() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let ids = [Some(CellId(0)), Some(CellId(1))];
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Fixed;
    faces[1] = FlowFace::Free;
    let divergent = FlowField::from_face_velocity(&domain, faces,
        |axis, centre| if axis == 0 { 0.125 * (1.0 + centre[0]) } else { 0.0 }).unwrap();
    assert!(matches!(mass_graph(&domain, &divergent, 1.0, input(&ids, &Caloric), &mut || Ok(())),
        Err(DepositionError::InvalidInput(_))));
    let reverse = FlowField::from_face_velocity(&domain, faces,
        |axis, _| if axis == 0 { -0.125 } else { 0.0 }).unwrap();
    let mut missing = input(&ids, &Caloric);
    missing.inlet_specific_h_j_kg[1] = None;
    assert!(matches!(mass_graph(&domain, &reverse, 1.0, missing, &mut || Ok(())),
        Err(DepositionError::InvalidInput(_))));
    let mut bounded = input(&ids, &Caloric);
    bounded.max_connections = 2; // One interior plus two boundary connections.
    assert!(matches!(mass_graph(&domain, &reverse, 1.0, bounded, &mut || Ok(())),
        Err(DepositionError::Budget(_))));
    assert!(matches!(mass_graph(&domain, &reverse, 1.0, input(&ids, &Caloric),
        &mut || Err(DepositionError::Cancelled)), Err(DepositionError::Cancelled)));
}

#[test]
fn advected_material_preserves_dynamic_surface_reservoir_heat_and_entropy() {
    use fs_conduction::enthalpy::conduction::{
        ConductivityLaw, ConductionObserver, ConductiveReservoir,
    };
    use fs_conduction::enthalpy::{
        Boundary, CellGeometry, DepositedCell, DepositionBudget, Reservoir, StepConfig, TimedEvent,
    };

    struct Conductivity;
    impl ConductivityLaw for Conductivity {
        fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
            Ok(1.0 + h / 100.0)
        }
    }
    // Observe the actual conduction owner and accepted intervals; this wrapper
    // supplies no independent flux, timestep or thermal update.
    #[derive(Clone)]
    struct Trace<'a> {
        conduction: ConductionObserver<'a>,
        pending_g: Option<f64>,
        steps: Vec<(f64, f64)>,
        cancel_after: Option<usize>,
    }
    impl DepositionObserver for Trace<'_> {
        type Error = DepositionError;
        fn begin(&mut self) -> Result<(), Self::Error> {
            self.conduction.begin()?;
            self.steps.clear();
            self.pending_g = None;
            Ok(())
        }
        fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
            check: &mut dyn FnMut() -> Result<(), DepositionError>)
            -> Result<BTreeMap<u64, f64>, Self::Error> {
            self.conduction.conductances(thermal, check)
        }
        fn reservoir_conductances(&mut self, thermal: &DepositionSimulation<'_>,
            check: &mut dyn FnMut() -> Result<(), DepositionError>)
            -> Result<BTreeMap<CellId, f64>, Self::Error> {
            let values = self.conduction.reservoir_conductances(thermal, check)?;
            self.pending_g = values.get(&CellId(0)).copied();
            Ok(values)
        }
        fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>)
            -> Result<(), Self::Error> {
            self.conduction.event(event, thermal)
        }
        fn step(&mut self, start: f64, dt: f64, temperatures: &BTreeMap<CellId, f64>,
            thermal: &DepositionSimulation<'_>,
            check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
            self.steps.push((dt, self.pending_g.take().expect("actual reservoir hook")));
            if self.cancel_after == Some(self.steps.len()) {
                return Err(DepositionError::Cancelled);
            }
            self.conduction.step(start, dt, temperatures, thermal, check)
        }
    }
    fn near(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-12 * b.abs().max(1.0), "{a} != {b}");
    }

    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Fixed;
    faces[1] = FlowFace::Free;
    let field = FlowField::from_face_velocity(&domain, faces,
        |axis, _| if axis == 0 { 0.25 } else { 0.0 }).unwrap();
    let caloric = Caloric;
    let conductivity = Conductivity;
    let source = FluidTransportState::try_new(2.0, 1.0, 1.0, 0.1).unwrap();
    let controls = StepConfig { max_dt_s: 0.125, energy_tolerance_j: 1e-10 };
    let events = (0..2).map(|i| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(i), segment_id: 0, body_id: 0, material_id: MaterialId(4),
            mass_kg: 2.0, specific_h_j_kg: 0.0,
            geometry: CellGeometry {
                start_m: [i as f64, 0.5, 0.5], end_m: [i as f64 + 1.0, 0.5, 0.5],
                volume_m3: 1.0, surface_area_m2: 6.0,
            },
            boundary: Boundary {
                ambient_k: 300.0, htc_w_m2_k: 0.0, source_w: 0.0,
                reservoir: (i == 0).then_some(Reservoir {
                    area_m2: 1.0, conductance_w_k: 0.0, temperature_k: 320.0,
                }),
            },
        }),
    }).collect();
    let old = DepositionSimulation::new(events, &[(MaterialId(4), &caloric)],
        DepositionBudget { max_cells: 2, max_links: 8, max_events: 2, max_substeps: 1000 })
        .unwrap().advance_to(0.0, controls, || Ok(())).unwrap().0;
    let conduction = ConductionObserver::new(&old, &[(MaterialId(4), &conductivity)], &[],
        &[ConductiveReservoir { cell_id: CellId(0), half_distance_m: 0.5,
            surface_resistance_m2_k_w: 0.25 }], || Ok(())).unwrap();
    let ids = [Some(CellId(0)), Some(CellId(1))];
    let mut binding = input(&ids, &caloric);
    binding.inlet_specific_h_j_kg = [Some(0.0); 6];
    let (flows, boundaries) = mass_graph(&domain, &field, 2.0, binding, &mut || Ok(())).unwrap();
    let observer = MaterialFlowAdvection {
        material: &source, transport: source, domain: &domain, field: &field, input: binding,
        advected: Advected::new(Trace { conduction, pending_g: None, steps: Vec::new(),
            cancel_after: None }, flows, 8).unwrap().with_boundaries(boundaries, 8).unwrap(),
    };
    let (next, accepted, energy, entropy) = old
        .advance_observed_with_entropy(&observer, 1.0, controls, || Ok(())).unwrap();
    let steps = &accepted.advection().inner().steps;
    assert!(steps.len() > 2);
    assert!(steps.last().unwrap().1 > steps[0].1);
    let mut h = [0.0; 2];
    let mut bath_heat = 0.0;
    let mut exported = 0.0;
    for &(dt, g) in steps {
        // Independent donor-cell energy recursion. Only the accepted durations
        // come from the owner; the half-cell/film law is checked independently.
        let expected_g = 1.0 / (0.5 / (1.0 + h[0] / 100.0) + 0.25);
        near(g, expected_g);
        let q = expected_g * (320.0 - (300.0 + h[0]));
        bath_heat += dt * q;
        exported += dt * 0.5 * h[1];
        let previous = h;
        h[0] += dt * (q - 0.5 * previous[0]) / 2.0;
        h[1] += dt * 0.5 * (previous[0] - previous[1]) / 2.0;
    }
    assert!(bath_heat > 0.0 && exported > 0.0 && h[1] > 0.0);
    for (i, expected) in h.into_iter().enumerate() {
        near(next.cells()[&CellId(i as u64)].specific_h_j_kg, expected);
        near(next.temperature_k(CellId(i as u64)).unwrap(), 300.0 + expected);
    }
    near(accepted.advection().audit().external_enthalpy_j, -exported);
    near(energy.external_heat_j, bath_heat - exported);
    near(energy.stored_energy_change_j, 2.0 * (h[0] + h[1]));
    near(entropy.reservoir_entropy_input_j_k, bath_heat / 320.0);
    near(next.mass_kg(), old.mass_kg());
    assert_eq!(next.cells()[&CellId(0)].boundary.reservoir.unwrap().conductance_w_k, 0.0);

    let mut cancelled = observer.clone();
    cancelled.advected = Advected::new(Trace {
        cancel_after: Some(2), ..observer.advection().inner().clone()
    }, observer.advection().flows().to_vec(), 8).unwrap()
        .with_boundaries(observer.advection().boundaries().to_vec(), 8).unwrap();
    assert!(matches!(old.advance_observed_with_entropy(&cancelled, 1.0, controls, || Ok(())),
        Err(DepositionError::Cancelled)));
    assert!(cancelled.advection().inner().steps.is_empty());
    assert_eq!(cancelled.advection().audit().substeps, 0);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.energy_j(), 0.0);
}
