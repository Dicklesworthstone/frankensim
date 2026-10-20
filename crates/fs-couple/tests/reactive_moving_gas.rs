//! Actual source-backed flow, chemistry and momentum; synthetic isomers only.
//! Reuse the existing synthetic card loader rather than inventing material data.
#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
use std::collections::BTreeMap;
use fs_couple::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use fs_couple::cavity_phs::piston::exchange::graph::{self, MovingGasDrive, MovingGasError, MovingGasSolve};
use fs_couple::cavity_phs::piston::exchange::graph::reactive::{self, ChamberReaction};
use fs_couple::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::{GasEnergySolve, parcel::ClosedGasParcel}};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{GasNozzle, network::snapshot::{self, SnapshotNozzleLink}};
use fs_ports::{ModuleWindowAudit, PortedModule};
use fs_thermochem::{Composition, CompositionBasis, SpeciesId, ElementId, ReactionId,
    ElementalMatrix, StoichiometricMatrix, ChargeVector, verify_conservation,
    Nasa9RegionV1, Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::equilibrium::{DeclaredStoichiometricMassActionV1, IdealGasReactionEquilibriumV1,
    ThermodynamicReverseRateClosureV1, ReactionProgressRateScaleV1, ReactionProgressRateQuantityV1,
    MassActionKineticOrderConventionV1};

fn law(release: f64, rate: f64) -> GasReaction {
    let base = source::load("A 0.25 0.028 3.5 0 2e-5 0.03\nB 0.25 0.028 4.0 0 3e-5 0.045\nI 0.5 0.004 2.5 0 1.4e-5 0.021\n",
        &format!("synthetic A/B isomers, declared product formation {} J/mol", -release)).unwrap();
    let mut components = base.thermodynamics().components().to_vec();
    let b = &components[1];
    let mut card = b.regions()[0].card().clone();
    card.parameters.get_mut("a7").unwrap().value = -release / R;
    components[1] = Nasa9StandardStateModelV1::new(b.species().clone(), b.molar_mass(),
        b.convention().clone(), vec![Nasa9RegionV1::from_card(card).unwrap()]).unwrap();
    let ids = vec![SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let elements = ElementalMatrix::new(vec![ElementId::new("X").unwrap()], ids.clone(), vec![vec![1, 1]]).unwrap();
    let reaction = ReactionId::new("synthetic-isomer-conversion").unwrap();
    let stoich = StoichiometricMatrix::new(ids.clone(), vec![reaction.clone()], vec![vec![-1], vec![1]]).unwrap();
    let charges = ChargeVector::new(ids, vec![0, 0]).unwrap();
    let conserved = verify_conservation(&elements, &stoich, &charges).unwrap();
    let equilibrium = IdealGasReactionEquilibriumV1::new(stoich, conserved, reaction, components[..2].to_vec()).unwrap();
    let kinetic = DeclaredStoichiometricMassActionV1::new(ThermodynamicReverseRateClosureV1::new(equilibrium,
        ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(rate)).unwrap()),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients);
    let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(components,
        base.thermodynamics().mole_composition().clone()).unwrap(), base.transport_materials().to_vec(), None).unwrap();
    GasReaction::new(model, kinetic).unwrap()
}
fn gas(law: &GasReaction, x: [f64; 3], p: f64, t: f64) -> ClosedGasParcel {
    let model = law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction, x.to_vec()).unwrap()).unwrap();
    ClosedGasParcel::new(model, source::point(t, p).unwrap(), 0.001).unwrap()
}
fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn chemical() -> GasReactionSolve {
    GasReactionSolve { gas: energy(), maximum_consumed_fraction: 0.2,
        species_tolerance_mol: 1e-13, entropy_tolerance_j_k: 1e-11 }
}
fn solve() -> MovingGasSolve {
    MovingGasSolve { exchange: GasSnapshotSolve { maximum_withdrawn_fraction: 0.1,
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }, total_energy_tolerance_j: 1e-7 }
}
fn drive(dt: f64) -> MovingGasDrive {
    MovingGasDrive { input: PistonInput { duration_s: dt, heat_j: 0.0, exterior_pressure_pa: 100000.0, external_force_n: 0.0 },
        solve: PistonSolve { volume_bracket_m3: [0.0007, 0.0015], gas: energy(),
            momentum_tolerance_n_s: 2e-7 * dt, energy_tolerance_j: 1e-7, maximum_iterations: 128 } }
}
fn links() -> Vec<SnapshotNozzleLink> {
    [(10, [0, 1], 1e-5), (20, [0, 2], 1e-5), (30, [1, 3], 1e-6), (40, [2, 3], 1e-6)]
        .map(|(id, parcels, area_m2)| SnapshotNozzleLink { id, parcels,
            nozzle: GasNozzle { area_m2, discharge_coefficient: 0.8 }, nozzle_solve: energy() }).to_vec()
}
fn owners(law: &GasReaction) -> (BTreeMap<u64, GasPiston>, BTreeMap<u64, ClosedGasParcel>) {
    let chamber = gas(law, [0.0, 0.0, 1.0], 100000.0, 600.0);
    (BTreeMap::from([(1, GasPiston::new(chamber.clone(), 1.0, 0.01, 0.0).unwrap()),
        (2, GasPiston::new(chamber, 1.0, 0.01, 0.0).unwrap())]),
        BTreeMap::from([(0, gas(law, [0.75, 0.0, 0.25], 300000.0, 800.0)),
            (3, gas(law, [0.0, 0.0, 1.0], 60000.0, 600.0))]))
}
fn endpoints() -> BTreeMap<u64, GasEnergySolve> { (0..4).map(|id| (id, energy())).collect() }
fn drives(dt: f64) -> BTreeMap<u64, MovingGasDrive> { BTreeMap::from([(1, drive(dt)), (2, drive(dt))]) }
fn bindings(law: &GasReaction) -> Vec<ChamberReaction<'_>> {
    [1, 2].map(|chamber| ChamberReaction { chamber, reaction: law, controls: chemical() }).to_vec()
}
fn amount(gas: &ClosedGasParcel, name: &str) -> f64 {
    let model = gas.model().thermodynamics();
    let n = gas.mass_kg() / gas.state().thermodynamics().properties().molar_mass().value();
    model.components().iter().zip(model.mole_composition().fractions())
        .find(|(s, _)| s.species().as_str() == name).map_or(0.0, |(_, x)| n * x)
}
fn delta(old: &ClosedGasParcel, new: &ClosedGasParcel) -> f64 {
    new.mass_kg() * (new.state().specific_internal_energy_j_kg() - old.state().specific_internal_energy_j_kg())
        + (new.mass_kg() - old.mass_kg()) * old.state().specific_internal_energy_j_kg()
}

#[test]
fn empty_chemistry_is_the_identical_existing_graph_path() {
    let law = law(8000.0, 1000.0);
    let (p, r) = owners(&law);
    let reference = graph::step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5), solve(), || false).unwrap();
    let actual = reactive::step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5), &[], solve(), || false).unwrap();
    assert_eq!(actual.pistons(), reference.pistons());
    assert_eq!(actual.reservoirs(), reference.reservoirs());
    assert_eq!(actual.step(), reference.step());
    assert!(actual.reactions().is_empty());
    for id in [1, 2] {
        assert_eq!(actual.pistons()[&id].last_window_transfers(), reference.pistons()[&id].last_window_transfers());
    }
}

#[test]
fn incoming_species_react_once_before_actual_wall_motion_without_a_heat_credit() {
    let law = law(8000.0, 1000.0);
    let (p, r) = owners(&law);
    let mut initial = r.clone();
    for (&id, piston) in &p { initial.insert(id, piston.gas().clone()); }
    let flow = snapshot::step_candidate(&initial, &links(), 5e-5, &endpoints(), solve().exchange, || false).unwrap();
    let actual = reactive::step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5), &bindings(&law), solve(), || false).unwrap();
    assert_eq!(&actual.step().flow, flow.step());
    assert_eq!(actual.step().supplied_heat_j, 0.0);
    assert!(actual.step().total_energy_residual_j.abs() < 1e-7);
    for id in [1, 2] {
        assert_eq!(amount(p[&id].gas(), "A"), 0.0);
        let chemical = law.parcel_candidate(&flow.parcels()[&id], 5e-5, chemical(), || false).unwrap();
        assert_eq!(&actual.reactions()[&id], chemical.step(), "motion probes must not react again");
        assert!(actual.reactions()[&id].extent_mol > 0.0);
        let mut expected = GasPiston::new(chemical.into_parts().0, 1.0, 0.01, 0.0).unwrap();
        let report = expected.step(drive(5e-5).input, drive(5e-5).solve, || false).unwrap();
        assert_eq!(&actual.pistons()[&id], &expected);
        assert_eq!(&actual.step().pistons[&id], &report);
        assert!(amount(actual.pistons()[&id].gas(), "B") > 0.0);
        assert!(matches!(actual.pistons()[&id].audit_last_window().unwrap(), ModuleWindowAudit::Unavailable { .. }));
    }
    assert_eq!(amount(&actual.reservoirs()[&3], "B"), 0.0, "new products cannot retroactively leave this snapshot");
    let (next_p, next_r, _, _) = actual.into_parts();
    let next = reactive::step_candidate(&next_p, &next_r, &links(), 5e-5, &endpoints(), &drives(5e-5), &bindings(&law), solve(), || false).unwrap();
    assert!(amount(&next.reservoirs()[&3], "B") > 0.0);
}

#[test]
fn chemistry_respects_symmetry_binding_order_and_original_link_order_independence() {
    let law = law(8000.0, 1000.0);
    let (p, r) = owners(&law);
    let run = |links: &[SnapshotNozzleLink], bindings: &[ChamberReaction<'_>]| {
        reactive::step_candidate(&p, &r, links, 5e-5, &endpoints(), &drives(5e-5), bindings, solve(), || false).unwrap()
    };
    let expected = run(&links(), &bindings(&law));
    let mut edges = links(); edges.reverse();
    let mut chemical = bindings(&law); chemical.reverse();
    assert_eq!(run(&edges, &chemical), expected);
    assert_eq!(expected.pistons()[&1], expected.pistons()[&2]);
    assert_eq!(expected.reactions()[&1], expected.reactions()[&2]);
    for (i, link) in edges.iter_mut().enumerate() { link.id = i as u64 + 200; }
    let renamed = run(&edges, &chemical);
    assert_eq!(renamed.pistons(), expected.pistons());
    assert_eq!(renamed.reservoirs(), expected.reservoirs());
    assert_eq!(renamed.reactions(), expected.reactions());
}

#[test]
fn closed_reacting_walls_and_finite_fixed_reactors_use_the_original_chemical_owners() {
    let law = law(8000.0, 1000.0);
    let piston = GasPiston::new(gas(&law, [0.75, 0.0, 0.25], 100000.0, 600.0), 1.0, 0.01, 0.0).unwrap();
    let fixed = gas(&law, [0.5, 0.0, 0.5], 100000.0, 600.0);
    let reference = piston.reactive_candidate(&law, drive(5e-5).input, chemical(), drive(5e-5).solve, || false).unwrap();
    let fixed_reference = law.parcel_candidate(&fixed, 5e-5, chemical(), || false).unwrap();
    let b: Vec<_> = [0, 1].map(|chamber| ChamberReaction { chamber, reaction: &law, controls: chemical() }).to_vec();
    let actual = reactive::step_candidate(&BTreeMap::from([(1, piston)]), &BTreeMap::from([(0, fixed)]), &[], 5e-5,
        &BTreeMap::from([(0, energy()), (1, energy())]), &BTreeMap::from([(1, drive(5e-5))]), &b, solve(), || false).unwrap();
    assert_eq!(&actual.pistons()[&1], reference.piston());
    assert_eq!(&actual.reactions()[&1], &reference.step().reaction);
    assert_eq!(&actual.reservoirs()[&0], fixed_reference.gas());
    assert_eq!(&actual.reactions()[&0], fixed_reference.step());
    let ModuleWindowAudit::Available(window) = actual.pistons()[&1].audit_last_window().unwrap() else {
        panic!("a closed reacting wall retains its complete single-interval port account");
    };
    assert!(window.violations().is_empty());
}

#[test]
fn repeated_flow_chemistry_and_motion_close_storage_and_chemical_species_budgets() {
    let law = law(8000.0, 1000.0);
    let (mut p, mut r) = owners(&law);
    let (initial_p, initial_r) = (p.clone(), r.clone());
    let mut progress = 0.0;
    let mut outside = 0.0;
    for i in 0..32 {
        // A physical change of rate policy: turn reactions off without resetting
        // material inventory, wall momentum, flow or the chemical product state.
        let b = if i < 16 { bindings(&law) } else { vec![] };
        let next = reactive::step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5), &b, solve(), || false).unwrap();
        progress += next.reactions().values().map(|s| s.extent_mol).sum::<f64>();
        outside += next.step().exterior_pressure_work_j;
        let (np, nr, _, _) = next.into_parts(); p = np; r = nr;
        let du: f64 = p.iter().map(|(id, w)| delta(initial_p[id].gas(), w.gas())).sum::<f64>()
            + r.iter().map(|(id, g)| delta(&initial_r[id], g)).sum::<f64>();
        let kinetic: f64 = p.values().map(|w| 0.5*w.mass_kg()*w.velocity_m_s().powi(2)).sum();
        assert!((du + kinetic + outside).abs() < 2e-7);
        let products: f64 = p.values().map(|w| amount(w.gas(), "B")).sum::<f64>()
            + r.values().map(|g| amount(g, "B")).sum::<f64>();
        assert!((products-progress).abs() < 1e-11);
    }
    assert!(progress > 0.0);
    assert!(amount(&r[&3], "B") > 0.0);
    assert_ne!(p[&1].gas().state().state().sound_speed, initial_p[&1].gas().state().state().sound_speed);
}

#[test]
fn invalid_bindings_late_reaction_refusal_and_final_cancellation_leave_all_inputs_intact() {
    let law = law(8000.0, 1000.0);
    let (p, r) = owners(&law);
    let original = (p.clone(), r.clone());
    let mut b = bindings(&law);
    let run = |b: &[ChamberReaction<'_>], cancel: &mut dyn FnMut() -> bool| {
        reactive::step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5), b, solve(), cancel)
    };
    b[1].chamber = 1;
    assert!(matches!(run(&b, &mut || false), Err(MovingGasError::Invalid(_))));
    b[1].chamber = 99;
    assert!(matches!(run(&b, &mut || false), Err(MovingGasError::Invalid(_))));
    b = bindings(&law);
    b[1].controls.maximum_consumed_fraction = 1e-12;
    assert!(matches!(run(&b, &mut || false), Err(MovingGasError::Reaction { id: 2, .. })));
    let mut polls = 0;
    run(&bindings(&law), &mut || { polls += 1; false }).unwrap();
    for stop in [1, polls] {
        let mut calls = 0;
        assert_eq!(run(&bindings(&law), &mut || { calls += 1; calls == stop }), Err(MovingGasError::Cancelled));
    }
    assert_eq!((&p, &r), (&original.0, &original.1));
    for id in [1, 2] { assert!(p[&id].last_window_transfers().is_none()); }
}

mod material {
    use super::*;
    use fs_couple::cavity_phs::piston::material::{MaterialGasPiston, graph as material_graph};
    use fs_evidence::ValidityDomain;
    use fs_ga::Vec3;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_material::state_point::MaterialPropertySelection;
    use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body, ResolvedBody};
    use fs_qty::{Density, Temperature};

    fn body(density: f64) -> ResolvedBody {
        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim { key: PropertyKey::new("density", Density::DIMS),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, density), (400.0, density)], dims: Density::DIMS },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
            observations: vec![], provenance: Provenance { source: "synthetic wall specimen".into(),
                license: "CC0-1.0".into(), artifact: None } }).unwrap();
        let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic wall".into(),
            phase: "solid".into(), process: "specified".into(), revision: 0 }, claims, vec![]).unwrap();
        let point = QueryPoint::new().with("T", 300.0).unwrap();
        resolve_body(&[RegionBinding { key: "wall", material: &card, point: &point,
            geometry: VolumeMoments::cuboid(Vec3::new(0.1, 0.1, 0.01)).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly }], 1, || false).unwrap()
    }
    fn material_owners(law: &GasReaction, density: f64)
        -> (BTreeMap<u64, MaterialGasPiston>, BTreeMap<u64, ClosedGasParcel>)
    {
        let (p, r) = owners(law);
        (p.into_iter().map(|(id, p)| (id, MaterialGasPiston::new(p.gas().clone(),
            body(if id == 1 { density } else { 10000.0 }), p.area_m2(), p.velocity_m_s()).unwrap())).collect(), r)
    }
    #[test]
    fn material_path_matches_shared_dynamics_and_retains_exact_source_records() {
        let law = law(8000.0, 1000.0);
        let (p, r) = material_owners(&law, 8000.0);
        let numeric = p.iter().map(|(&id, p)| (id, p.piston().clone())).collect();
        let reference = reactive::step_candidate(&numeric, &r, &links(), 5e-5, &endpoints(), &drives(5e-5),
            &bindings(&law), solve(), || false).unwrap();
        let actual = material_graph::reactive_step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5),
            &bindings(&law), solve(), || false).unwrap();
        assert_eq!(actual.step(), reference.step());
        assert_eq!(actual.reactions(), reference.reactions());
        assert_eq!(actual.reservoirs(), reference.reservoirs());
        for id in [1, 2] {
            assert_eq!(actual.pistons()[&id].piston(), &reference.pistons()[&id]);
            assert_eq!(actual.pistons()[&id].body(), p[&id].body());
            assert_eq!(actual.pistons()[&id].piston().mass_kg(), p[&id].body().inertia().mass());
        }
        let reference = material_graph::step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5), solve(), || false).unwrap();
        let empty = material_graph::reactive_step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5), &[], solve(), || false).unwrap();
        assert_eq!(empty.pistons(), reference.pistons()); assert_eq!(empty.step(), reference.step());
    }
    fn run(density: f64, rate: f64, count: u32) -> (BTreeMap<u64, MaterialGasPiston>, BTreeMap<u64, ClosedGasParcel>) {
        let law = law(8000.0, rate);
        let (mut p, mut r) = material_owners(&law, density);
        let original = p.clone();
        let dt = 0.0016 / f64::from(count);
        for _ in 0..count {
            let next = material_graph::reactive_step_candidate(&p, &r, &links(), dt, &endpoints(), &drives(dt),
                &bindings(&law), solve(), || false).unwrap();
            assert!(next.step().total_energy_residual_j.abs() <= solve().total_energy_tolerance_j);
            let (np, nr, _, _) = next.into_parts(); p = np; r = nr;
            for id in [1, 2] { assert_eq!(p[&id].body(), original[&id].body()); }
        }
        (p, r)
    }
    #[test]
    fn changing_only_density_and_rate_data_propagates_through_reacting_flow_and_motion() {
        let (light, light_r) = run(1000.0, 1000.0, 64);
        let (dense, dense_r) = run(8000.0, 1000.0, 64);
        let (faster, _) = run(8000.0, 2000.0, 64);
        assert!(light[&1].piston().velocity_m_s() > dense[&1].piston().velocity_m_s());
        assert_ne!(light[&1].gas().state().state().temperature, dense[&1].gas().state().state().temperature);
        assert_ne!(light[&1].gas().state().state().sound_speed, dense[&1].gas().state().state().sound_speed);
        assert_ne!(light_r[&3].mass_kg(), dense_r[&3].mass_kg());
        assert!(amount(faster[&1].gas(), "B") > amount(dense[&1].gas(), "B"));
        assert!(faster[&1].piston().velocity_m_s() > dense[&1].piston().velocity_m_s());
    }
    #[test]
    fn material_reacting_trajectory_refines_at_the_declared_first_order() {
        let values = [32, 64, 128].map(|n| run(8000.0, 1000.0, n).0[&1].piston().velocity_m_s());
        let ratio = (values[0]-values[1]).abs()/(values[1]-values[2]).abs();
        assert!((1.8..2.2).contains(&ratio), "{values:?}, ratio {ratio}");
    }
    #[test]
    fn cancellation_after_material_rebinding_and_arriving_source_conflicts_are_atomic() {
        let law = law(8000.0, 1000.0);
        let (p, r) = material_owners(&law, 8000.0);
        let original = (p.clone(), r.clone());
        let mut polls = 0;
        material_graph::reactive_step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5),
            &bindings(&law), solve(), || { polls += 1; false }).unwrap();
        let mut calls = 0;
        assert_eq!(material_graph::reactive_step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5),
            &bindings(&law), solve(), || { calls += 1; calls == polls }), Err(MovingGasError::Cancelled));
        // The transported A/I sources, not just the originally inert chamber,
        // must match the selected reaction's original source catalog.
        let conflicting = super::law(9000.0, 1000.0);
        let b = bindings(&conflicting);
        assert!(matches!(material_graph::reactive_step_candidate(&p, &r, &links(), 5e-5, &endpoints(), &drives(5e-5),
            &b, solve(), || false), Err(MovingGasError::Reaction { .. })));
        assert_eq!((&p, &r), (&original.0, &original.1));
    }
}
