use super::*;
use fs_material::thermovisco::{Limits, ThermalKnot};

pub(super) fn material(e0: f64, branches: Vec<(f64, f64)>) -> ThermoMaxwell {
    ThermoMaxwell::new(e0, branches, [300.0, 500.0].into_iter().map(|temperature_k| ThermalKnot {
        temperature_k, modulus_scale: 1.0, relaxation_scale: 1.0,
        free_strain: (temperature_k-300.0)*1e-5,
    }).collect(), Limits { max_strain: 0.03, max_strain_rate_s: 10.0 }).unwrap()
}
fn geometry() -> LayerGeometry {
    LayerGeometry { area_m2: 1e-6, centroid_z_m: 0.0, thickness_m: 0.001 }
}
pub(super) fn config() -> Config {
    Config { section: deposition::Config { force_tolerance_n: 1e-7, moment_tolerance_nm: 1e-9,
        energy_tolerance_j: 1e-8, max_strain: 0.03, max_rotation_rad: 0.1 },
        maximum_deflection_m: 0.001 }
}
fn empty(n: u32, support: Support) -> Bridge<'static> {
    Bridge::new(&(0..n).map(|i| (u64::from(i), 0.02/f64::from(n))).collect::<Vec<_>>(), support,
        Budget { max_cells: 128, max_layers_per_cell: 8, max_branch_states: 4096 }).unwrap()
}
pub(super) fn birth<'a>(n: u32, support: Support, law: &'a ThermoMaxwell) -> Bridge<'a> {
    // Empty sections borrow no law, so their covariant lifetime can be shortened.
    let bridge: Bridge<'a> = empty(n, support);
    bridge.deposit_bonded_layer(7, geometry(), 1200.0, law,
        &(0..n).map(|i| (u64::from(i), 300.0)).collect(), &mut || false).unwrap()
}
pub(super) fn temperatures(bridge: &Bridge<'_>, t: f64) -> BTreeMap<(u64, u64), f64> {
    bridge.cells().iter().flat_map(|cell| cell.section().layers().keys()
        .map(move |&layer| ((cell.id(), layer), t))).collect()
}
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() < tolerance, "{a:e} != {b:e}"); }

// G1: independently known elastic beam solutions; spatial error is not a time error.
#[test]
fn elastic_gravity_sag_refines_quadratically_for_both_supports() {
    let law = material(1e8, vec![]);
    let ei = 1e8*0.001*0.001_f64.powi(3)/12.0;
    let q = 1200.0*1e-6*9.81;
    for support in [Support::SimplySupported, Support::CantileverLeft] {
        let mut errors = Vec::new();
        for n in [4, 8, 16] {
            let bridge = birth(n, support, &law);
            let (next, report) = bridge.advance_controlled(1.0, &temperatures(&bridge,300.0), 9.81,
                &vec![0.0; n as usize], config(), &mut || false).unwrap();
            let l = next.length_m();
            let (x, exact) = match support {
                Support::SimplySupported => (l/2.0, -5.0*q*l.powi(4)/(384.0*ei)),
                Support::CantileverLeft => (l, -q*l.powi(4)/(8.0*ei)),
            };
            let value = next.deformed_position_m(x, 0.0).unwrap()[1];
            errors.push((value-exact).abs());
            near(report.vertical_reactions_n.iter().sum(), q*l, 1e-12);
            near(report.maximum_sag_m, -value, 1e-12);
            near(report.load_work_j, report.energy_j.mechanical_work, 1e-12);
            near(report.energy_j.residual, 0.0, 1e-12);
            assert_eq!(report.material_point_updates, 2*n as usize);
        }
        for pair in errors.windows(2) { assert!((3.99..4.01).contains(&(pair[0]/pair[1]))); }
    }
}

#[test]
fn nonuniform_signed_loads_preserve_global_reaction_and_work_balance() {
    let law = material(1e8, vec![]);
    let bridge = Bridge::new(&[(40,0.003),(10,0.007),(20,0.01)], Support::SimplySupported,
        Budget { max_cells: 3,max_layers_per_cell: 2,max_branch_states: 8 }).unwrap()
        .deposit_bonded_layer(7,geometry(),1200.0,&law,
            &BTreeMap::from([(10,300.0),(20,300.0),(40,300.0)]),&mut || false).unwrap();
    let loads = [0.01,-0.003,0.02];
    let (next, r) = bridge.advance_controlled(1.0,&temperatures(&bridge,300.0),0.0,
        &loads,config(),&mut || false).unwrap();
    let weight = 0.01*0.003-0.003*0.007+0.02*0.01;
    let moment = 0.01*0.003*0.0015-0.003*0.007*0.0065+0.02*0.01*0.015;
    near(r.vertical_reactions_n.iter().sum(),weight,1e-14);
    near(r.vertical_reactions_n[1]*0.02,moment,1e-14);
    near(next.deformed_position_m(next.length_m(),0.0).unwrap()[1],0.0,1e-12);
    near(r.load_work_residual_j,0.0,1e-12);
}

// G1: exact standard-linear-solid creep compliance under a step in dead load.
#[test]
fn maxwell_creep_refines_in_time_without_resetting_spring_history() {
    let law = material(1e7, vec![(9e7,0.2)]);
    let elastic = material(1e8, vec![]);
    let baseline = birth(8,Support::SimplySupported,&elastic);
    let (_, r) = baseline.advance_controlled(1.0,&temperatures(&baseline,300.0),9.81,
        &[0.0;8],config(),&mut || false).unwrap();
    let compliance = (1.0-0.9*(-1.0_f64).exp())/1e7;
    let exact = r.maximum_sag_m*1e8*compliance;
    let mut errors = Vec::new();
    for n in [20,40,80] {
        let mut state = birth(8,Support::SimplySupported,&law);
        let t = temperatures(&state,300.0);
        let mut sag = 0.0;
        for i in 1..=n {
            let (next, r) = state.advance_controlled(2.0*f64::from(i)/f64::from(n),&t,9.81,
                &[0.0;8],config(),&mut || false).unwrap();
            assert!(r.maximum_sag_m >= sag);
            sag = r.maximum_sag_m; state = next;
        }
        errors.push((sag-exact).abs());
        assert!(state.cells()[0].section().layers()[&7].material_points()[0].branch_strains()[0] != 0.0);
    }
    assert!(errors[0]/errors[1] > 1.9 && errors[1]/errors[2] > 1.9);
}

#[test]
fn near_maxwell_fluid_limit_retains_elastic_jump_plus_viscous_creep() {
    let law = material(100.0,vec![(1e8,0.1)]);
    let elastic = material(1e8,vec![]);
    let base = birth(8,Support::SimplySupported,&elastic);
    let (_, r) = base.advance_controlled(1.0,&temperatures(&base,300.0),9.81,
        &[0.0;8],config(),&mut || false).unwrap();
    let mut state = birth(8,Support::SimplySupported,&law);
    let t = temperatures(&state,300.0);
    for i in 1..=100 {
        state = state.advance_controlled(f64::from(i)/100.0,&t,9.81,&[0.0;8],
            config(),&mut || false).unwrap().0;
    }
    let sag = -state.deformed_position_m(state.length_m()/2.0,0.0).unwrap()[1];
    // E_eq remains positive: this is a time window below the retardation time,
    // not an admitted zero-stiffness liquid or large-deformation strand model.
    near(sag/(11.0*r.maximum_sag_m),1.0,2e-5);
}

#[test]
fn differential_thermal_contraction_curves_the_spatial_span() {
    let law = material(1e8,vec![]);
    let mut bridge: Bridge<'_> = empty(8,Support::SimplySupported);
    for (layer,z) in [(7,-0.00025),(8,0.00025)] {
        bridge = bridge.deposit_bonded_layer(layer,LayerGeometry { area_m2: 5e-7,
            centroid_z_m: z,thickness_m: 0.0005 },1200.0,&law,
            &(0..8).map(|id| (id,500.0)).collect(),&mut || false).unwrap();
    }
    let t = bridge.cells().iter().flat_map(|c| [((c.id(),7),300.0),((c.id(),8),400.0)]).collect();
    let (cooled,r) = bridge.advance_controlled(1.0,&t,0.0,&[0.0;8],config(),&mut || false).unwrap();
    // Equal-thickness bilayer: kappa = -3 Delta(eigenstrain)/(2 total depth).
    near(cooled.cells()[0].section().curvature_per_m(),-1.5,1e-10);
    near(cooled.deformed_position_m(0.01,0.0).unwrap()[1],1.5*0.02_f64.powi(2)/8.0,1e-12);
    near(r.load_work_j,0.0,1e-14);
    near(r.energy_j.residual,0.0,1e-12);
}

// G0/G5: full-span birth does not erase old stress and same-grid replay agrees.
#[test]
fn deposited_mass_and_stress_free_birth_preserve_existing_history_and_shape() {
    let law = material(1e8,vec![(2e8,1.0)]);
    let first = birth(8,Support::SimplySupported,&law);
    let t = temperatures(&first,300.0);
    let loaded = first.advance_controlled(1.0,&t,9.81,&[0.0;8],config(),&mut || false).unwrap().0;
    let born = loaded.deposit_bonded_layer(8,LayerGeometry { centroid_z_m: 0.001,..geometry() },
        1200.0,&law,&(0..8).map(|id| (id,400.0)).collect(),&mut || false).unwrap();
    near(born.mass_kg(),2.0*loaded.mass_kg(),1e-14);
    for (old,new) in loaded.cells().iter().zip(born.cells()) {
        assert_eq!(old.section().layers()[&7],new.section().layers()[&7]);
        for p in new.section().layers()[&8].material_points() { assert_eq!(p.stress_pa(),0.0); }
    }
    assert_eq!(loaded.deformed_position_m(0.01,0.0),born.deformed_position_m(0.01,0.0));
    let t = temperatures(&born,350.0);
    let direct = born.advance_controlled(2.0,&t,9.81,&[0.0;8],config(),&mut || false).unwrap().0;
    let checkpoint = born.clone();
    let replay = checkpoint.advance_controlled(2.0,&t,9.81,&[0.0;8],config(),&mut || false).unwrap().0;
    assert_eq!(direct,replay);
}

// G4: a late section/material/geometry refusal cannot publish a partial beam.
#[test]
fn invalid_data_global_limits_and_cancellation_preserve_the_original() {
    let law = material(1e8,vec![]);
    let state = birth(8,Support::SimplySupported,&law);
    let before = state.clone();
    let t = temperatures(&state,300.0);
    for allowed in [0,8,20,40] {
        let mut count = 0;
        assert!(state.advance_controlled(1.0,&t,9.81,&[0.0;8],config(),&mut || {
            count += 1; count > allowed
        }).is_err());
        assert_eq!(state,before);
    }
    let mut bad_t = t.clone(); bad_t.insert((7,7),600.0);
    assert!(state.advance_controlled(1.0,&bad_t,9.81,&[0.0;8],config(),&mut || false).is_err());
    assert!(state.advance_controlled(1.0,&t,9.81,&[0.0;8],
        Config { maximum_deflection_m: 1e-20,..config() },&mut || false).is_err());
    assert!(empty(8,Support::SimplySupported).advance_controlled(1.0,&BTreeMap::new(),9.81,
        &[0.0;8],config(),&mut || false).is_err());
    assert!(state.deposit_bonded_layer(8,geometry(),1200.0,&law,
        &(0..8).map(|id| (id,300.0)).collect(),&mut || false).is_err());
    assert_eq!(state,before);
}
