use super::*;
use crate::bridge::tests::{birth, config, material, temperatures};

fn delivered<'a>(id: u64, length_m: f64, temperature_k: f64, law: &'a ThermoMaxwell) -> CellBirth<'a> {
    CellBirth { id, length_m, mass_kg: 1200.0*1e-6*length_m, temperature_k, law }
}
fn near(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol, "{a:e} != {b:e}"); }

#[test]
fn birth_continues_the_current_tip_without_aging_old_or_new_material() {
    let law = material(1e8, vec![(9e7, 0.2)]);
    let seed = birth(8, Support::CantileverLeft, &law);
    let (loaded, _) = seed.advance_controlled(0.4, &temperatures(&seed, 350.0), 9.81,
        &[0.0;8], config(), &mut || false).unwrap();
    let old_tip = loaded.deformed_position_m(0.02, 0.0).unwrap();
    let old_energy = loaded.stored_energy_j();
    let (grown, report) = loaded.append_cantilever_controlled(
        &[delivered(50, 0.005, 500.0, &law)], config(), &mut || false).unwrap();
    assert_eq!(&grown.cells()[..8], loaded.cells());
    assert_eq!(grown.time_s(), 0.4);
    assert_eq!(report.previous_tip_m, old_tip);
    assert_eq!(report.added_cells, 1);
    assert_eq!(grown.stored_energy_j(), old_energy);
    let added = &grown.cells()[8];
    assert_eq!(added.section().time_s(), 0.4);
    assert_eq!(added.section().layers()[&7].birth_time_s, 0.4);
    for point in added.section().layers()[&7].material_points() {
        assert_eq!(point.temperature_k(), 500.0);
        assert_eq!(point.stress_pa(), 0.0);
        assert!(point.branch_strains().iter().all(|&q| q == 0.0));
    }
    assert!(added.slope != 0.0);
    near(grown.deformed_position_m(0.02, 0.0).unwrap()[1], old_tip[1], 1e-15);
    near(report.tip_m[1]-old_tip[1], 0.005*added.slope, 1e-15);
    near(grown.mass_kg()-loaded.mass_kg(), report.incoming_mass_kg, 1e-16);
    near(report.incoming_reference_volume_m3, 5e-9, 1e-20);
}

#[test]
fn enlarged_elastic_load_field_matches_independent_full_span_and_refines() {
    let law = material(1e8, vec![]);
    let ei = 1e8*0.001*0.001_f64.powi(3)/12.0;
    let q = 1200.0*1e-6*9.81;
    let mut errors = Vec::new();
    for n in [4_u32, 8, 16] {
        let seed = birth(n, Support::CantileverLeft, &law);
        let (loaded, r0) = seed.advance_controlled(0.1, &temperatures(&seed, 300.0), 9.81,
            &vec![0.0;n as usize], config(), &mut || false).unwrap();
        near(r0.vertical_reactions_n[0], q*0.02, 1e-12);
        let additions: Vec<_> = (0..n).map(|i| delivered(100+u64::from(i), 0.02/f64::from(n), 300.0, &law)).collect();
        let (grown, _) = loaded.append_cantilever_controlled(&additions, config(), &mut || false).unwrap();
        let (solved, report) = grown.advance_controlled(0.2, &temperatures(&grown, 300.0), 9.81,
            &vec![0.0;2*n as usize], config(), &mut || false).unwrap();
        let length = solved.length_m();
        let tip = solved.deformed_position_m(length, 0.0).unwrap()[1];
        errors.push((tip+q*length.powi(4)/(8.0*ei)).abs());
        near(report.vertical_reactions_n[0], q*length, 1e-12);
        near(report.clamp_reaction_nm, q*length*length/2.0, 1e-12);
        near(report.energy_j.residual, 0.0, 1e-12);
        near(report.load_work_residual_j, 0.0, 1e-12);
        // Build the same final reference mesh directly, independent of birth history.
        let mesh: Vec<_> = grown.cells().iter().map(|c| (c.id(), c.length_m())).collect();
        let whole = Bridge::new(&mesh, Support::CantileverLeft, grown.budget).unwrap()
            .deposit_bonded_layer(7, grown.cells[0].section.layers()[&7].geometry, 1200.0, &law,
                &mesh.iter().map(|&(id,_)| (id,300.0)).collect(), &mut || false).unwrap();
        let (whole, _) = whole.advance_controlled(0.2, &temperatures(&whole,300.0),9.81,
            &vec![0.0;mesh.len()], config(), &mut || false).unwrap();
        near(tip, whole.deformed_position_m(length,0.0).unwrap()[1],1e-14);
    }
    assert!(errors[0]/errors[1] > 3.9 && errors[1]/errors[2] > 3.9, "{errors:?}");
}

#[test]
fn later_hot_material_has_its_own_thermal_birth_reference() {
    let law = material(1e8, vec![]);
    let seed = birth(4, Support::CantileverLeft, &law);
    let (hot, _) = seed.advance_controlled(1.0,&temperatures(&seed,500.0),0.0,&[0.0;4],
        config(),&mut || false).unwrap();
    let (grown, _) = hot.append_cantilever_controlled(&[delivered(20,0.005,500.0,&law)],
        config(),&mut || false).unwrap();
    near(grown.cells()[0].section().axial_strain(),0.002,1e-14);
    assert_eq!(grown.cells()[4].section().axial_strain(),0.0);
    let (cold, _) = grown.advance_controlled(2.0,&temperatures(&grown,300.0),0.0,&[0.0;5],
        config(),&mut || false).unwrap();
    near(cold.cells()[0].section().axial_strain(),0.0,1e-14);
    near(cold.cells()[4].section().axial_strain(),-0.002,1e-14);
    near(cold.deformed_position_m(0.025,0.0).unwrap()[0],0.025-0.005*0.002,1e-14);
    for c in cold.cells() { for p in c.section().layers()[&7].material_points() {
        near(p.stress_pa(),0.0,1e-6);
    } }
}

#[test]
fn batch_and_same_time_segmented_birth_replay_identical_histories() {
    let law = material(1e8,vec![(9e7,0.2)]);
    let seed = birth(4,Support::CantileverLeft,&law);
    let (old,_) = seed.advance_controlled(0.1,&temperatures(&seed,300.0),9.81,&[0.0;4],config(),&mut || false).unwrap();
    let cells = [delivered(10,0.002,400.0,&law),delivered(11,0.003,450.0,&law)];
    let (batch,_) = old.append_cantilever_controlled(&cells,config(),&mut || false).unwrap();
    let (prefix,_) = old.append_cantilever_controlled(&cells[..1],config(),&mut || false).unwrap();
    let (split,_) = prefix.clone().append_cantilever_controlled(&cells[1..],config(),&mut || false).unwrap();
    assert_eq!(batch,split);
    let expected = batch.advance_controlled(0.2,&temperatures(&batch,350.0),9.81,
        &[0.0;6],config(),&mut || false).unwrap().0;
    let actual = split.advance_controlled(0.2,&temperatures(&split,350.0),9.81,
        &[0.0;6],config(),&mut || false).unwrap().0;
    assert_eq!(expected,actual);
}

#[test]
fn future_stations_cannot_support_the_seed_before_their_material_is_born() {
    use crate::bridge::supports::{Config as Supports, Stop};
    use fs_solver::orthant::GramOrthantConfig;
    let law = material(1e8,vec![]);
    let seed = birth(8,Support::CantileverLeft,&law);
    let support = [Stop { id: 80,reference_x_m: 0.025,minimum_deflection_m: 0.0 }];
    let cfg = Supports { contact: GramOrthantConfig { max_variables: 4,max_columns: 64,
        max_setup_products: 4096,max_sweeps: 4096,gradient_tolerance: 1e-13 },
        gap_tolerance_m: 1e-11,complementarity_tolerance_j: 1e-13 };
    assert!(seed.advance_supported_controlled(0.1,&temperatures(&seed,300.0),9.81,&[0.0;8],
        &support,config(),cfg,&mut || false).is_err());
    let (grown,_) = seed.append_cantilever_controlled(&[delivered(20,0.005,300.0,&law)],config(),&mut || false).unwrap();
    let (next,r) = grown.advance_supported_controlled(0.1,&temperatures(&grown,300.0),9.81,&[0.0;9],
        &support,config(),cfg,&mut || false).unwrap();
    assert!(r.reactions[0].upward_force_n > 0.0);
    near(next.deformed_position_m(0.025,0.0).unwrap()[1],0.0,1e-11);
}

#[test]
fn invalid_or_interrupted_batches_publish_no_material_or_history() {
    let law = material(1e8,vec![(9e7,0.2)]);
    let seed = birth(4,Support::CantileverLeft,&law);
    let good = delivered(10,0.003,400.0,&law);
    for bad in [CellBirth { id: 0,..good },CellBirth { id: 11,length_m: 0.0,..good },
        CellBirth { id: 11,mass_kg: f64::NAN,..good },CellBirth { id: 11,temperature_k: 600.0,..good }] {
        assert!(seed.append_cantilever_controlled(&[good,bad],config(),&mut || false).is_err());
    }
    let mut limited = seed.clone(); limited.budget.max_cells = 4;
    assert!(limited.append_cantilever_controlled(&[good],config(),&mut || false).is_err());
    limited.budget.max_cells = 8; limited.budget.max_branch_states = limited.branch_states;
    assert!(limited.append_cantilever_controlled(&[good],config(),&mut || false).is_err());
    let other = CellBirth { id: 11,..good };
    for limit in [0,3,6,7] {
        let mut polls = 0;
        assert!(seed.append_cantilever_controlled(&[good,other],config(),&mut || {
            polls+=1; polls > limit
        }).is_err());
        assert_eq!(seed.cells().len(),4); assert_eq!(seed.time_s(),0.0);
    }
    let span = birth(4,Support::SimplySupported,&law);
    assert!(span.append_cantilever_controlled(&[good],config(),&mut || false).is_err());
    let (same,report) = seed.append_cantilever_controlled(&[],config(),&mut || false).unwrap();
    assert_eq!(same,seed); assert_eq!(report.incoming_mass_kg,0.0);
}
