//! G1/G3/G4: entropy on the original swept-contact thermal candidate.
//! Independent temperature-coordinate linear solves; no measured-material claim.
use fs_couple::friction_heat::{FrictionHeatEntropy, FrictionHeatTolerances, moving::*};
use fs_la::prepared_lu::LuWorkspace;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::cell::Cell;
use std::collections::BTreeMap;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > -300.0 { Ok(300.0 + h) }
        else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
    fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64)
        -> Result<Option<f64>, DepositionError> {
        let t = self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b - a) / t).ln_1p()))
    }
}
static SENSIBLE: Sensible = Sensible;
const MASS: [f64; 3] = [1.0, 1.0, 4.0];
const AREA: [f64; 3] = [6.0, 6.0, 10.0];
const AMBIENT: [f64; 3] = [290.0, 320.0, 295.0];
fn close(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol, "{a:.16e} != {b:.16e}; tolerance {tol:e}");
}
fn tolerance() -> FrictionHeatTolerances {
    FrictionHeatTolerances { cell_j: 1e-10, thermal_j: 1e-10,
        interface_j: 1e-10, total_j: 1e-10 }
}
fn position(x: f64) -> MovingContactPosition {
    MovingContactPosition { body_offsets_m: [0.0; 2], footprint_offset_m: x }
}
fn fixture(law: &'static dyn ThermalLaw, initial_h: [f64; 3], htc: [f64; 3])
    -> (DepositionSimulation<'static>, MovingFrictionRoute) {
    let mut events = Vec::new();
    for i in 0..3 {
        let x = if i == 2 { 1.0 } else { i as f64 + 0.5 };
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(i as u64), segment_id: i as u64,
            body_id: if i == 2 { 1 } else { 0 }, material_id: MaterialId(0),
            mass_kg: MASS[i], specific_h_j_kg: initial_h[i],
            geometry: CellGeometry { start_m: [x, 0.0, 0.0],
                end_m: [x, if i == 2 { 1.0 } else { -1.0 }, 0.0],
                volume_m3: if i == 2 { 2.0 } else { 1.0 }, surface_area_m2: AREA[i] },
            boundary: Boundary { ambient_k: AMBIENT[i], htc_w_m2_k: htc[i],
                reservoir: None, source_w: 0.0 },
        }) });
    }
    let state = DepositionSimulation::new(events, &[(MaterialId(0), law)],
        DepositionBudget { max_cells: 3, max_links: 0, max_events: 3, max_substeps: 1000 })
        .unwrap().advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-10 },
            || Ok(())).unwrap().0;
    let route = MovingFrictionRoute::new(&state,
        StripFrame { origin_m: [0.0; 3], tangent: [1.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0], width_m: 1.0 },
        [vec![StripSurfaceCell { cell: CellId(0), span_m: [0.0, 1.0] },
              StripSurfaceCell { cell: CellId(1), span_m: [1.0, 2.0] }],
         vec![StripSurfaceCell { cell: CellId(2), span_m: [0.0, 2.0] }]],
        [0.0, 1.0], MovingFrictionBudget { maximum_surface_cells: 3,
            maximum_overlap_tests: 4, maximum_patches: 4 }, || Ok(())).unwrap();
    (state, route)
}
fn input(work: &[f64], sources: [f64; 3]) -> MovingFrictionInput<()> {
    MovingFrictionInput { dissipated_j: work.to_vec(),
        conductivity_w_m_k: (0..3).map(|i| (CellId(i), 1.0)).collect(),
        transport: StepCoefficients {
            sources_w: sources.into_iter().enumerate().filter(|(_, q)| *q != 0.0)
                .map(|(i, q)| (CellId(i as u64), q)).collect(),
            conductances_w_k: BTreeMap::new(), payload: (),
        } }
}
fn step<'s>(state: &'s DepositionSimulation<'static>, route: &MovingFrictionRoute, dt: f64)
    -> MovingFrictionHeatStep<'s, 'static> {
    MovingFrictionHeatStep::new(state, state.time_s() + dt, route,
        MovingContactPath { start: position(0.0), end: position(1.0) }, || Ok(())).unwrap()
}
// Eliminate the two massless interfaces, then solve the THREE original stores.
// For the full sweep each overlap is 1/2 m2 and each half-conductance is 1 W/K.
fn endpoint(old: &DepositionSimulation<'_>, dt: f64, work: [f64; 2], sources: [f64; 3])
    -> (Vec<f64>, [f64; 2]) {
    let coverage = [0.5, 0.5, 1.0];
    let mut matrix = vec![0.0; 9];
    let mut rhs = vec![0.0; 3];
    for i in 0..3 {
        let c = &old.cells()[&CellId(i as u64)];
        let g = (AREA[i] - coverage[i]) * c.boundary.htc_w_m2_k;
        matrix[3*i+i] = MASS[i]/dt + g;
        rhs[i] = MASS[i]/dt * old.temperature_k(c.id).unwrap() + g*AMBIENT[i] + sources[i];
    }
    for i in 0..2 {
        matrix[3*i+i] += 0.5; matrix[8] += 0.5;
        matrix[3*i+2] -= 0.5; matrix[6+i] -= 0.5;
        rhs[i] += work[i]/(2.0*dt); rhs[2] += work[i]/(2.0*dt);
    }
    let mut lu = LuWorkspace::new(3).unwrap();
    lu.factor(&matrix).unwrap();
    let mut t = vec![0.0; 3];
    lu.solve_factored_into(&rhs, &mut t).unwrap();
    let energy = (0..3).map(|i| MASS[i]*(t[i]-old.temperature_k(CellId(i as u64)).unwrap())).collect();
    let interfaces = std::array::from_fn(|i| 0.5*(work[i]/dt+t[i]+t[2]));
    (energy, interfaces)
}
fn balance(e: FrictionHeatEntropy) {
    close(e.thermal_storage_change_j_k, e.reservoir_entropy_input_j_k
        + e.network_conduction_generation_j_k + e.interface_conduction_generation_j_k
        + e.friction_generation_j_k + e.unclassified_source_at_cell_j_k
        + e.interface_residual_entropy_j_k + e.temporal_defect_j_k, 2e-12);
}

#[test]
fn moving_mask_removes_the_actual_reservoir_entropy_and_resolves_source_once() {
    let (old, route) = fixture(&SENSIBLE, [0.0, 10.0, -5.0], [0.2, 0.3, 0.4]);
    let dt = 0.25;
    let step = step(&old, &route, dt);
    assert_eq!(step.patches().len(), 2);
    for p in step.patches() { close(p.mean_area_m2, 0.5, 0.0); }
    let work = [4.0, 7.0];
    let (energy, interfaces) = endpoint(&old, dt, work, [0.0; 3]);
    let (plain, (), plain_report) = step.evaluate(&energy, &interfaces,
        |_, _, _| Ok::<_, DepositionError>(input(&work, [0.0; 3])), || Ok(()))
        .unwrap().accept(tolerance(), || Ok(())).unwrap();
    let calls = Cell::new(0);
    let (next, (), report) = step.evaluate_with_entropy(&energy, &interfaces,
        |_, _, _| { calls.set(calls.get()+1);
            Ok::<_, DepositionError>(input(&work, [0.0; 3])) }, || Ok(()))
        .unwrap().accept(tolerance(), || Ok(())).unwrap();
    assert_eq!(calls.get(), 1);
    assert_eq!(next.cells(), plain.cells());
    assert_eq!(report.thermal, plain_report.thermal);
    assert!(plain_report.entropy.is_none());
    let e = report.entropy.unwrap();
    assert!(!e.has_unclassified_sources, "covered-area cancellation is NOT an unknown source");
    assert_eq!(e.unclassified_source_at_cell_j_k, 0.0);
    let mut storage = 0.0;
    let mut boundary = 0.0;
    let mut generation = 0.0;
    let mut temporal = 0.0;
    for i in 0..3 {
        let t = next.temperature_k(CellId(i as u64)).unwrap();
        let t0 = old.temperature_k(CellId(i as u64)).unwrap();
        let ds = MASS[i]*(t/t0).ln();
        storage += ds;
        temporal += ds-MASS[i]*(t-t0)/t;
        let q = dt * ([5.5, 5.5, 9.0][i]) * [0.2, 0.3, 0.4][i] * (AMBIENT[i]-t);
        boundary += q/AMBIENT[i];
        generation += q*(1.0/t-1.0/AMBIENT[i]);
    }
    close(e.thermal_storage_change_j_k, storage, 2e-12);
    close(e.reservoir_entropy_input_j_k, boundary, 2e-14);
    close(e.network_conduction_generation_j_k, generation, 2e-14);
    close(e.temporal_defect_j_k, temporal, 2e-12);
    close(e.friction_generation_j_k, work[0]/interfaces[0]+work[1]/interfaces[1], 1e-15);
    assert!(e.interface_conduction_generation_j_k > 0.0);
    assert_eq!(report.position, position(1.0));
    balance(e);
}

#[test]
fn zero_slip_redistributes_heat_without_friction_or_duplicate_shared_storage() {
    let (old, route) = fixture(&SENSIBLE, [20.0, -10.0, 0.0], [0.0; 3]);
    let step = step(&old, &route, 0.5);
    let (energy, interfaces) = endpoint(&old, 0.5, [0.0; 2], [0.0; 3]);
    let (next, (), r) = step.evaluate_with_entropy(&energy, &interfaces,
        |_, _, _| Ok::<_, DepositionError>(input(&[0.0; 2], [0.0; 3])), || Ok(()))
        .unwrap().accept(tolerance(), || Ok(())).unwrap();
    let e = r.entropy.unwrap();
    assert_eq!(e.friction_generation_j_k, 0.0);
    assert_eq!(e.reservoir_entropy_input_j_k, 0.0);
    assert!(e.interface_conduction_generation_j_k > 0.0);
    assert!(e.thermal_storage_change_j_k > 0.0);
    close(next.energy_j(), old.energy_j(), 1e-10);
    let expected: f64 = (0..3).map(|i| MASS[i]
        * (next.temperature_k(CellId(i as u64)).unwrap()
            / old.temperature_k(CellId(i as u64)).unwrap()).ln()).sum();
    close(e.thermal_storage_change_j_k, expected, 2e-12);
    balance(e);
}

#[test]
fn prescribed_source_is_not_confused_with_masking_or_friction() {
    let (old, route) = fixture(&SENSIBLE, [0.0; 3], [0.2, 0.3, 0.4]);
    let step = step(&old, &route, 0.25);
    let extra = [0.1, -0.2, 0.3];
    let work = [4.0, 7.0];
    let (energy, interfaces) = endpoint(&old, 0.25, work, extra);
    let (next, (), r) = step.evaluate_with_entropy(&energy, &interfaces,
        |_, _, _| Ok::<_, DepositionError>(input(&work, extra)), || Ok(()))
        .unwrap().accept(tolerance(), || Ok(())).unwrap();
    let e = r.entropy.unwrap();
    let expected: f64 = (0..3).map(|i| 0.25*extra[i]
        / next.temperature_k(CellId(i as u64)).unwrap()).sum();
    assert!(e.has_unclassified_sources);
    close(e.unclassified_source_at_cell_j_k, expected, 1e-15);
    balance(e);
}

#[test]
fn moving_latent_heat_stores_entropy_once_at_the_original_plateau() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if (0.0..=100.0).contains(&h) { Ok(300.0) }
            else { Err(DepositionError::MaterialDomain) }
        }
        fn maximum_temperature_slope(&self) -> f64 { 0.0 }
        fn specific_entropy_change_j_kg_k(&self, a: f64, b: f64)
            -> Result<Option<f64>, DepositionError> {
            self.temperature_k(a)?; self.temperature_k(b)?; Ok(Some((b-a)/300.0))
        }
    }
    static PLATEAU: Plateau = Plateau;
    let (old, route) = fixture(&PLATEAU, [10.0; 3], [0.0; 3]);
    let step = step(&old, &route, 1.0);
    let (next, (), r) = step.evaluate_with_entropy(&[2.0, 4.0, 6.0], &[302.0, 304.0],
        |_, _, _| Ok::<_, DepositionError>(input(&[4.0, 8.0], [0.0; 3])), || Ok(()))
        .unwrap().accept(tolerance(), || Ok(())).unwrap();
    let e = r.entropy.unwrap();
    for i in 0..3 { assert_eq!(next.temperature_k(CellId(i)).unwrap(), 300.0); }
    close(e.thermal_storage_change_j_k, 12.0/300.0, 1e-15);
    close(e.friction_generation_j_k, 4.0/302.0+8.0/304.0, 1e-15);
    close(e.temporal_defect_j_k, 0.0, 1e-15);
    assert!(!e.has_unclassified_sources);
    balance(e);
}

#[test]
fn missing_caloric_entropy_and_late_cancellation_preserve_the_original_history() {
    struct Missing;
    impl ThermalLaw for Missing {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> { SENSIBLE.temperature_k(h) }
        fn maximum_temperature_slope(&self) -> f64 { 1.0 }
    }
    static MISSING: Missing = Missing;
    let (old, route) = fixture(&MISSING, [0.0; 3], [0.0; 3]);
    let saved = old.cells().clone();
    let step = step(&old, &route, 0.25);
    let (energy, interfaces) = endpoint(&old, 0.25, [4.0, 7.0], [0.0; 3]);
    assert!(step.evaluate(&energy, &interfaces,
        |_, _, _| Ok::<_, DepositionError>(input(&[4.0, 7.0], [0.0; 3])), || Ok(())).is_ok());
    assert!(matches!(step.evaluate_with_entropy(&energy, &interfaces,
        |_, _, _| Ok::<_, DepositionError>(input(&[4.0, 7.0], [0.0; 3])), || Ok(())),
        Err(DepositionError::InvalidInput(_))));
    assert_eq!(old.cells(), &saved);
    let (other, route) = fixture(&SENSIBLE, [0.0; 3], [0.0; 3]);
    let candidate = MovingFrictionHeatStep::new(&other, 0.25, &route,
        MovingContactPath { start: position(0.0), end: position(1.0) }, || Ok(())).unwrap()
        .evaluate_with_entropy(&energy, &interfaces,
            |_, _, _| Ok::<_, DepositionError>(input(&[4.0, 7.0], [0.0; 3])), || Ok(())).unwrap();
    assert!(matches!(candidate.accept(tolerance(), || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)));
    assert_eq!(other.cells(), &saved);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(other.time_s(), 0.0);
}
