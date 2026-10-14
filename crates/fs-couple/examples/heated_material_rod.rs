//! Source data -> shared enthalpy transport -> shared material rod equilibrium.
//! The example is axially guided, quasistatic and one-way thermomechanical.
//! Mechanical energy is observed separately, NOT silently added to the heat audit.
#[path = "heated_material_rod/material.rs"]
mod material;

use fs_conduction::transient::deposition::{Boundary, CellGeometry, CellId, ContactPatchUpdate,
    DepositedCell, DepositionBudget, DepositionError, DepositionSimulation, EnergyAudit,
    Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent};
use fs_exec::CancelGate;
use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection,
    ResolvedMaterialStatePoint, THERMAL_CONDUCTIVITY_PROPERTY};
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::{SpatialMaterialRod, supports::RodSupport};
use fs_solid::rod::material_field::thermal::source::{RodLayerSource, RodMaterialSource, RodSegmentSource};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
const N: usize = 4;
const LENGTH: f64 = 0.04;
const RADIUS: f64 = 0.002;
const HTC: f64 = 10.0;
fn area() -> f64 { std::f64::consts::PI * RADIUS * RADIUS }
fn step_config() -> StepConfig { StepConfig { max_dt_s: 0.01, energy_tolerance_j: 1e-8 } }
fn poll(gate: &CancelGate) -> std::result::Result<(), DepositionError> {
    if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) }
}
fn k(state: &ResolvedMaterialStatePoint) -> f64 {
    state.property(THERMAL_CONDUCTIVITY_PROPERTY).expect("conductivity was required").value_si()
}
struct Experiment<'a> {
    material: &'a material::Material,
    source: RodMaterialSource<'a>,
    fixed_end: bool,
    power_w: f64,
    duration_s: f64,
}
#[derive(Clone)]
struct State<'a> {
    heat: DepositionSimulation<'a>,
    rod: SpatialMaterialRod,
    // Retain the actual source snapshots used to set each interval's links.
    applied_conductivity: Vec<ResolvedMaterialStatePoint>,
    external_heat_j: f64,
}
impl<'a> Experiment<'a> {
    fn new(material: &'a material::Material, fixed_end: bool, power_w: f64,
        duration_s: f64, gate: &CancelGate) -> Result<Self> {
        poll(gate)?;
        if !power_w.is_finite() || power_w < 0.0 || !duration_s.is_finite()
            || duration_s <= 0.0 || duration_s > 1000.0 {
            return Err("example needs nonnegative power and duration in (0,1000] seconds".into());
        }
        let mut segments = Vec::new();
        for _ in 0..N {
            segments.push(RodSegmentSource { inner_radius_m: 0.0, shear_correction: 0.8,
                layers: vec![RodLayerSource { outer_radius_m: RADIUS, card: &material.card,
                    reference_point: material::point(300.0)?,
                    mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
                    expansion_selection: MaterialPropertySelection::SingleClaimOnly }] });
        }
        let source = RodMaterialSource::try_new(LENGTH, segments, "shear_modulus", 0.03, gate)?;
        let reference = &source.reference_rod().sections()[0].elastic().layers()[0].material;
        let density_claim = reference.property(DENSITY_PROPERTY).unwrap().answer().receipt.selected;
        if material.caloric.curve().material_card_identity() != reference.card_identity()
            || material.caloric.receipts().iter().filter(|r| r.property == DENSITY_PROPERTY)
                .any(|r| r.selected != density_claim || r.axis_quantities != *reference.axis_quantities()
                    || r.query_point.iter().any(|(axis, _)| axis != "T")) {
            return Err("caloric and mechanical density must select the same card, claim and fixed coordinates".into());
        }
        Ok(Self { material, source, fixed_end, power_w, duration_s })
    }
    fn boundary(&self, cell: usize, powered: bool) -> Boundary {
        Boundary { ambient_k: 300.0, htc_w_m2_k: HTC, reservoir: None,
            source_w: if cell == 0 && powered { self.power_w } else { 0.0 } }
    }
    fn resolve_mechanics(&self, temperature: &[f64], gate: &CancelGate) -> Result<SpatialMaterialRod> {
        let rod = self.source.resolve_at(temperature, gate)?;
        // Explicit straight axial guides, NOT a claim about a buckling-free beam.
        let supports: Vec<_> = (1..=N).map(|node| RodSupport { node,
            position_m: [if self.fixed_end && node == N { Some(LENGTH) } else { None }, Some(0.0), Some(0.0)],
            frame: Some(self.source.reference_rod().frames()[node]),
        }).collect();
        let mut rod = rod.with_supports(&supports, gate)?;
        rod.solve_static(&TipLoad::default(), 1, 1e-8, gate)?;
        Ok(rod)
    }
    fn initial(&self, gate: &CancelGate) -> Result<State<'a>> {
        let dx = LENGTH/N as f64;
        let h = self.material.caloric.state_at_temperature(300.0)?.specific_enthalpy_j_kg();
        let mut events = Vec::new();
        let mut applied_conductivity = Vec::new();
        for cell in 0..N {
            poll(gate)?;
            applied_conductivity.push(self.material.conductivity_at(300.0)?);
            // Use source-resolved reference mass, not current rho(T)*fixed volume.
            let mass = dx*self.source.reference_rod().sections()[cell].elastic().line_density_kg_m();
            events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
                id: CellId(cell as u64), segment_id: cell as u64, body_id: 1, material_id: MaterialId(1),
                mass_kg: mass, specific_h_j_kg: h,
                geometry: CellGeometry { start_m: [cell as f64*dx, 0.0, 0.0],
                    end_m: [(cell+1) as f64*dx, 0.0, 0.0], volume_m3: area()*dx,
                    surface_area_m2: 2.0*std::f64::consts::PI*RADIUS*dx + 2.0*area() },
                boundary: self.boundary(cell, true),
            }) });
        }
        for cell in 1..N {
            events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
                id: cell as u64, a: CellId((cell-1) as u64), b: CellId(cell as u64),
                area_m2: area(), conductance_w_k: self.conductance(&applied_conductivity, cell)?,
            }) });
        }
        events.push(TimedEvent { time_s: self.duration_s/2.0,
            event: Event::Boundary(CellId(0), self.boundary(0, false)) });
        let law: &'a dyn ThermalLaw = self.material.caloric.curve();
        let heat = DepositionSimulation::new(events, &[(MaterialId(1), law)], DepositionBudget {
            max_cells: N, max_links: N-1, max_events: 2*N, max_substeps: 200_000,
        })?.advance_to(0.0, step_config(), || poll(gate))?.0;
        let rod = self.resolve_mechanics(&[300.0; N], gate)?;
        Ok(State { heat, rod, applied_conductivity, external_heat_j: 0.0 })
    }
    fn conductance(&self, states: &[ResolvedMaterialStatePoint], cell: usize) -> Result<f64> {
        let dx = LENGTH/N as f64;
        // Existing network receives the series resistance of two half cells.
        let g = area()/(0.5*dx/k(&states[cell-1]) + 0.5*dx/k(&states[cell]));
        if !g.is_finite() || g <= 0.0 { return Err("axial thermal conductance is unrepresentable".into()); }
        Ok(g)
    }
    fn advance(&self, state: &State<'a>, target_s: f64, gate: &CancelGate) -> Result<(State<'a>, EnergyAudit)> {
        poll(gate)?;
        let mut applied_conductivity = Vec::with_capacity(N);
        for cell in 0..N {
            let t = state.heat.temperature_k(CellId(cell as u64))?;
            let property = self.material.conductivity_at(t)?;
            // Necessary external-film check only; axial spatial accuracy still
            // needs refinement and k remains lagged within this coupling interval.
            if HTC*RADIUS/(2.0*k(&property)) > fs_conduction::lumped::LUMPED_BIOT_CEILING {
                return Err("external-film Biot exceeds the reduced-cell limit".into());
            }
            applied_conductivity.push(property);
        }
        let mut updates = Vec::with_capacity(N-1);
        for cell in 1..N {
            let original = state.heat.contact_patch(cell as u64).ok_or("missing adjacent thermal link")?;
            let mut next = original;
            next.link.conductance_w_k = self.conductance(&applied_conductivity, cell)?;
            updates.push(ContactPatchUpdate { id: cell as u64, expected: Some(original), replacement: Some(next) });
        }
        let candidate = state.heat.replace_contact_patches(&updates, || poll(gate))?;
        let (heat, audit) = candidate.advance_to(target_s, step_config(), || poll(gate))?;
        let mut temperatures = Vec::with_capacity(N);
        for cell in 0..N {
            let t = heat.temperature_k(CellId(cell as u64))?;
            // Endpoint source support is checked even on the last observation.
            self.material.conductivity_at(t)?;
            temperatures.push(t);
        }
        let rod = self.resolve_mechanics(&temperatures, gate)?;
        let external_heat_j = state.external_heat_j + audit.external_heat_j;
        if !external_heat_j.is_finite() { return Err("thermal account overflow".into()); }
        poll(gate)?;
        // No accepted heat/time/support state changes before ALL owners succeed.
        Ok((State { heat, rod, applied_conductivity, external_heat_j }, audit))
    }
}
fn run(args: &[String]) -> Result<()> {
    if args.len() != 6 {
        return Err("usage: heated_material_rod --demo CP_J_KG_K|--pack MATERIAL_PACK free|fixed POWER_W DURATION_S INTERVALS".into());
    }
    let material = match args[0].as_str() {
        "--demo" => material::demo(args[1].parse()?, 1.0, 500.0)?,
        "--pack" => material::Material::load(&args[1])?,
        _ => return Err("choose --demo (synthetic data) or --pack (canonical material card pack)".into()),
    };
    let fixed = match args[2].as_str() { "free" => false, "fixed" => true, _ => return Err("support must be free or fixed".into()) };
    let intervals: usize = args[5].parse()?;
    if intervals == 0 || intervals > 10_000 { return Err("interval count must lie in 1..=10000".into()); }
    let gate = CancelGate::new();
    let experiment = Experiment::new(&material, fixed, args[3].parse()?, args[4].parse()?, &gate)?;
    let mut state = experiment.initial(&gate)?;
    let initial_heat = state.heat.energy_j();
    eprintln!("source={}, length_m={LENGTH}, radius_m={RADIUS}, segments={N}, h_W_m2_K={HTC}, reference_K=300, caloric_span_K=[300,500], axial_guides=true, fixed_end={fixed}, reference_mass_kg={}; one-way quasistatic mechanics, no melting or total thermomechanical energy claim",
        material.card.content_hash().to_hex(), experiment.source.reference_mass_kg());
    println!("time_s,T0_K,T1_K,T2_K,T3_K,tip_extension_m,axial_force_N,elastic_energy_J,stored_heat_change_J,external_heat_J,heat_residual_J,applied_k0_W_m_K");
    for step in 1..=intervals {
        let target = experiment.duration_s * (step as f64 / intervals as f64);
        state = experiment.advance(&state, target, &gate)?.0;
        let temperatures = (0..N).map(|c| state.heat.temperature_k(CellId(c as u64)))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let heat_change = state.heat.energy_j()-initial_heat;
        println!("{target},{},{},{},{},{},{},{},{},{},{},{}", temperatures[0], temperatures[1], temperatures[2], temperatures[3],
            state.rod.positions_m()[N][0]-LENGTH, state.rod.segment_resultants(N-1)?.0[0],
            state.rod.energy_j()?, heat_change, state.external_heat_j, heat_change-state.external_heat_j,
            k(&state.applied_conductivity[0]));
    }
    Ok(())
}
fn main() -> Result<()> { run(&std::env::args().skip(1).collect::<Vec<_>>()) }

#[cfg(test)]
mod tests {
    use super::*;
    fn at_half(material: &material::Material, fixed: bool) -> State<'_> {
        let gate = CancelGate::new();
        let experiment = Experiment::new(material, fixed, 1.0, 4.0, &gate).unwrap();
        let mut state = experiment.initial(&gate).unwrap();
        for i in 1..=10 { state = experiment.advance(&state, 0.2*i as f64, &gate).unwrap().0; }
        state
    }
    #[test]
    fn heat_capacity_data_changes_actual_heating_and_equilibrium_not_just_a_getter() {
        let a = material::demo(500.0, 1.0, 500.0).unwrap();
        let b = material::demo(1000.0, 1.0, 500.0).unwrap();
        let sa = at_half(&a, false); let sb = at_half(&b, false);
        assert!(sa.heat.temperature_k(CellId(0)).unwrap() > sb.heat.temperature_k(CellId(0)).unwrap()+1.0);
        assert!(sa.rod.positions_m()[N][0] > sb.rod.positions_m()[N][0]);
        for s in [&sa, &sb] { assert!((s.heat.energy_j()-s.external_heat_j).abs() < 1e-8); }
    }
    #[test]
    fn conductivity_data_changes_spatial_heat_and_restrained_reaction() {
        let a = material::demo(500.0, 1.0, 500.0).unwrap();
        let b = material::demo(500.0, 4.0, 500.0).unwrap();
        let sa = at_half(&a, true); let sb = at_half(&b, true);
        assert!(sa.heat.temperature_k(CellId(0)).unwrap() > sb.heat.temperature_k(CellId(0)).unwrap());
        assert!(sb.heat.temperature_k(CellId(1)).unwrap() > sa.heat.temperature_k(CellId(1)).unwrap());
        assert_eq!(sa.rod.positions_m()[N][0], LENGTH);
        assert!(sa.rod.segment_resultants(N-1).unwrap().0[0] < -1e-5);
        assert!((sa.rod.segment_resultants(N-1).unwrap().0[0]-sb.rod.segment_resultants(N-1).unwrap().0[0]).abs() > 1e-8);
    }
    #[test]
    fn power_off_preserves_heat_history_and_cooling_changes_the_next_equilibrium() {
        let material = material::demo(500.0, 1.0, 500.0).unwrap();
        let gate = CancelGate::new();
        let experiment = Experiment::new(&material, false, 1.0, 4.0, &gate).unwrap();
        let middle = at_half(&material, false);
        let final_state = experiment.advance(&middle, 4.0, &gate).unwrap().0;
        assert!(final_state.heat.energy_j() < middle.heat.energy_j());
        assert!(final_state.rod.positions_m()[N][0] < middle.rod.positions_m()[N][0]);
        assert!(final_state.heat.temperature_k(CellId(0)).unwrap() > 300.0);
        assert!((final_state.heat.energy_j()-final_state.external_heat_j).abs() < 1e-8);
    }
    #[test]
    fn late_source_refusal_and_cancellation_leave_heat_time_and_pose_unchanged() {
        let material = material::demo(500.0, 1.0, 301.0).unwrap();
        let gate = CancelGate::new();
        let experiment = Experiment::new(&material, false, 1.0, 4.0, &gate).unwrap();
        let old = experiment.initial(&gate).unwrap();
        let cells = old.heat.cells().clone(); let positions = old.rod.positions_m().to_vec();
        assert!(experiment.advance(&old, 0.5, &gate).is_err());
        assert_eq!(old.heat.time_s(), 0.0); assert_eq!(old.heat.cells(), &cells);
        assert_eq!(old.rod.positions_m(), positions);
        gate.request(); assert!(experiment.advance(&old, 0.01, &gate).is_err());
        assert_eq!(old.heat.cells(), &cells);
    }
}
