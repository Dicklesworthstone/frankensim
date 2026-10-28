//! Numerical geometry feedback on the original manufactured NASA-9 gas owners.
use super::*;
use crate::sublimation::gas::{tests as fixture, GasSurfaceBinding};
use crate::sublimation::gas::chamber::{GasChamberSolve, advance::advance_chamber_adaptive};
use fs_thermal::sublimation::{SolidSurfaceState, SublimationLaw};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn controls() -> GasChamberAdvanceSolve {
    GasChamberAdvanceSolve { step: GasChamberSolve { surface: fixture::controls(), maximum_surfaces: 8 },
        absolute_mass_tolerance_kg: 1e-14, absolute_temperature_tolerance_k: 1e-5,
        relative_tolerance: 0.0, minimum_dt_s: 1e-12, maximum_dt_s: 0.025, maximum_trials: 10000 }
}
fn geometry(shape: RetreatGeometry) -> MovingChamberGeometry {
    MovingChamberGeometry { surfaces: BTreeMap::from([(1, shape)]),
        maximum_solid_volume_fraction: 0.001, absolute_extent_tolerance_m: 1e-10 }
}
fn surface<'a>(id: u64, phase: &'a dyn SublimationLaw, binding: &'a GasSurfaceBinding,
    chart: (f64,f64), t: f64, shape: RetreatGeometry,
    caloric: &'a dyn Fn(f64) -> Result<f64, fs_thermal::sublimation::SublimationError>,
) -> GasChamberSurface<'a> {
    let mut solid = fixture::solid(t, chart);
    solid.mass_kg = 1e-6;
    solid.area_m2 = shape.at_mass(solid.mass_kg).unwrap().area_m2;
    GasChamberSurface { id, law: phase, binding, solid, temperature_from_enthalpy: caloric }
}

// Independent simultaneous ODE in [solid mass, solid T, vapor mass, gas T].
// The fixture has cp_B=4R/M_B, cv_A=2.5R/M_A, cv_B=3R/M_B and h_s=cp_B*T-L.
fn derivative(y: [f64;4], carrier_mass: f64) -> [f64;4] {
    let [m, ts, mv, tg] = y;
    let ra = R/0.02897;
    let rb = R/0.018015;
    let area = 4.0*core::f64::consts::PI*(3.0*m/(4.0*core::f64::consts::PI*917.0)).powf(2.0/3.0);
    let ps = 1000.0*(200000.0/rb*(1.0/400.0-1.0/ts)).exp();
    let pv = mv*rb*tg/0.01;
    let outgoing = area*0.1*ps/(2.0*core::f64::consts::PI*rb*ts).sqrt();
    let incoming = area*0.1*pv/(2.0*core::f64::consts::PI*rb*tg).sqrt();
    let dm = outgoing-incoming;
    let eo = 3.5*rb*ts;
    let ei = 3.5*rb*tg;
    let energy = dm*eo+incoming*(eo-ei);
    [-dm, (dm*(0.5*rb*ts-200000.0)-incoming*(eo-ei))/(m*4.0*rb), dm,
        (energy-dm*3.0*rb*tg)/(carrier_mass*2.5*ra+mv*3.0*rb)]
}
fn reference(carrier_mass: f64) -> [f64;4] {
    let mut y = [1e-6, 400.0, 0.0, 400.0];
    let dt = 0.2/4096.0;
    for _ in 0..4096 {
        let a = derivative(y, carrier_mass);
        let b = derivative(std::array::from_fn(|i| y[i]+0.5*dt*a[i]), carrier_mass);
        let c = derivative(std::array::from_fn(|i| y[i]+0.5*dt*b[i]), carrier_mass);
        let d = derivative(std::array::from_fn(|i| y[i]+dt*c[i]), carrier_mass);
        y = std::array::from_fn(|i| y[i]+dt*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0);
    }
    y
}

#[test]
fn shrinking_sphere_changes_the_trajectory_and_refines_against_simultaneous_rates() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let invert = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let shape = RetreatGeometry::sphere(917.0).unwrap();
    let rows = [surface(1, &phase, &binding, chart, 400.0, shape, &invert)];
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let exact = reference(gas.mass_kg());
    let mut loose = controls();
    loose.absolute_mass_tolerance_kg = 1e-9;
    loose.absolute_temperature_tolerance_k = 0.05;
    let mut loose_geometry = geometry(shape);
    loose_geometry.absolute_extent_tolerance_m = 1e-5;
    let a = advance_chamber_adaptive_moving(&rows, &loose_geometry, None, &gas, &catalog,
        0.2, loose, || false).unwrap();
    let b = advance_chamber_adaptive_moving(&rows, &geometry(shape), None, &gas, &catalog,
        0.2, controls(), || false).unwrap();
    let temperature = |state: &MovingChamberCandidate<'_>| state.history().audit().surfaces[0].surface_temperature_k[1];
    assert!((temperature(&b)-exact[1]).abs() < 0.1*(temperature(&a)-exact[1]).abs());
    assert!((b.geometry()[&1].final_geometry.mass_kg-exact[0]).abs()
        < 0.1*(a.geometry()[&1].final_geometry.mass_kg-exact[0]).abs());
    let fixed = advance_chamber_adaptive(&rows, &gas, &catalog, 0.2, controls(), || false).unwrap();
    assert!(b.geometry()[&1].final_geometry.mass_kg > fixed.endpoint().surfaces()[0].solid.mass_kg + 1e-9);
    let g = b.geometry()[&1];
    assert!(g.final_geometry.area_m2 < g.initial.area_m2);
    assert_eq!(b.history().endpoint().surfaces()[0].solid.area_m2, g.final_geometry.area_m2);
    fixture::near(g.final_geometry.mass_kg, 917.0*4.0*core::f64::consts::PI*g.final_geometry.extent_m.powi(3)/3.0, 1e-20);
    fixture::near(g.normal_recession_m, b.history().audit().surfaces[0].recession_m, 1e-15);
    fixture::near(b.history().audit().energy_residual_j, 0.0, 1e-7);
    fixture::near(b.history().audit().mass_residual_kg, 0.0, 1e-14);
    assert!(b.history().audit().maximum_accepted_error_ratio <= 1.0);
}

#[test]
fn planar_constant_area_retains_original_physical_states_exactly() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let invert = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let shape = RetreatGeometry::layer(917.0, 1e-5).unwrap();
    let rows = [surface(1, &phase, &binding, chart, 400.0, shape, &invert)];
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let mut geom = geometry(shape);
    geom.absolute_extent_tolerance_m = 1.0;
    let original = advance_chamber_adaptive(&rows, &gas, &catalog, 0.01, controls(), || false).unwrap();
    let moving = advance_chamber_adaptive_moving(&rows, &geom, None, &gas, &catalog, 0.01, controls(), || false).unwrap();
    assert_eq!(moving.gas(), original.gas());
    assert_eq!(moving.history().endpoint().surfaces()[0].solid, original.endpoint().surfaces()[0].solid);
    assert_eq!(moving.history().audit().trials, original.audit().trials);
    fixture::near(moving.geometry()[&1].normal_recession_m,
        (rows[0].solid.mass_kg-moving.geometry()[&1].final_geometry.mass_kg)/(917.0*1e-5), 1e-15);
}

#[test]
fn mixed_shapes_growth_and_reflection_share_one_order_independent_gas_history() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let invert = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let a = RetreatGeometry::sphere(917.0).unwrap();
    let b = RetreatGeometry::cylinder(917.0, 0.001).unwrap();
    let rows = [surface(1, &phase, &binding, chart, 400.0, a, &invert),
        surface(2, &phase, &binding, chart, 300.0, b, &invert)];
    let geom = MovingChamberGeometry { surfaces: BTreeMap::from([(1,a),(2,b)]), ..geometry(a) };
    let law = ThermalAccommodation::new(&catalog, &[(SpeciesId::new("A").unwrap(), 0.1),
        (fixture::species(), 0.2)], [200.0,900.0]).unwrap();
    let laws = BTreeMap::from([(1,law.clone()),(2,law)]);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let forward = advance_chamber_adaptive_moving(&rows, &geom, Some(&laws), &gas, &catalog,
        0.01, controls(), || false).unwrap();
    let reversed = advance_chamber_adaptive_moving(&[rows[1],rows[0]], &geom, Some(&laws),
        &gas, &catalog, 0.01, controls(), || false).unwrap();
    assert_eq!(forward, reversed);
    assert!(forward.geometry()[&1].normal_recession_m > 0.0);
    assert!(forward.geometry()[&2].normal_recession_m < 0.0);
    assert!(forward.geometry()[&2].final_geometry.area_m2 > forward.geometry()[&2].initial.area_m2);
    assert!(forward.history().audit().surfaces.iter().any(|row| row.reflected_heat_into_gas_j != 0.0));
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(forward.gas().species_mass_kg(&carrier).unwrap(), gas.species_mass_kg(&carrier).unwrap(), 1e-14);
    fixture::near(forward.history().audit().gas_stream_energy_residual_j, 0.0, 1e-7);
    fixture::near(forward.history().audit().solid_stream_energy_residual_j, 0.0, 1e-7);
}

#[test]
fn continuation_keeps_actual_area_and_cancellation_publishes_no_geometry_or_gas() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let invert = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let shape = RetreatGeometry::sphere(917.0).unwrap();
    let rows = [surface(1, &phase, &binding, chart, 400.0, shape, &invert)];
    let geom = MovingChamberGeometry { absolute_extent_tolerance_m: 1.0, ..geometry(shape) };
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let saved = gas.clone();
    let controls = GasChamberAdvanceSolve { maximum_dt_s: 0.005, absolute_mass_tolerance_kg: 1.0,
        absolute_temperature_tolerance_k: 100.0, ..controls() };
    let whole = advance_chamber_adaptive_moving(&rows,&geom,None,&gas,&catalog,0.01,controls,|| false).unwrap();
    let first = advance_chamber_adaptive_moving(&rows,&geom,None,&gas,&catalog,0.005,controls,|| false).unwrap();
    let continued = [GasChamberSurface { solid: first.history().endpoint().surfaces()[0].solid, ..rows[0] }];
    let second = advance_chamber_adaptive_moving(&continued,&geom,None,first.gas(),&catalog,0.005,controls,|| false).unwrap();
    assert_eq!(second.gas(), whole.gas());
    assert_eq!(second.geometry()[&1].final_geometry, whole.geometry()[&1].final_geometry);
    fixture::near(first.geometry()[&1].normal_recession_m+second.geometry()[&1].normal_recession_m,
        whole.geometry()[&1].normal_recession_m,1e-15);
    let mut calls = 0;
    advance_chamber_adaptive_moving(&rows,&geom,None,&gas,&catalog,0.01,controls,|| { calls+=1;false }).unwrap();
    for stop in [1,calls/2,calls] {
        let mut count = 0;
        assert!(advance_chamber_adaptive_moving(&rows,&geom,None,&gas,&catalog,0.01,controls,
            || { count+=1;count==stop }).is_err());
        assert_eq!(gas,saved);
        assert_eq!(rows[0].solid.mass_kg,1e-6);
    }
}

#[test]
fn stale_area_density_missing_shapes_and_nondilute_volume_refuse_even_at_zero_time() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog,&binding);
    let invert = |h| fixture::temperature(h,chart);
    let phase = fixture::phase();
    let shape = RetreatGeometry::sphere(917.0).unwrap();
    let row = surface(1,&phase,&binding,chart,400.0,shape,&invert);
    let gas = fixture::gas(&catalog,400.0,0.0);
    let geom = geometry(shape);
    let stale = GasChamberSurface { solid: SolidSurfaceState { area_m2: 1e-4, ..row.solid }, ..row };
    assert!(advance_chamber_adaptive_moving(&[stale],&geom,None,&gas,&catalog,0.0,controls(),|| false).is_err());
    let missing = MovingChamberGeometry { surfaces: BTreeMap::new(), ..geom.clone() };
    assert!(advance_chamber_adaptive_moving(&[row],&missing,None,&gas,&catalog,0.0,controls(),|| false).is_err());
    let dense = MovingChamberGeometry { maximum_solid_volume_fraction: 1e-10, ..geom.clone() };
    assert!(advance_chamber_adaptive_moving(&[row],&dense,None,&gas,&catalog,0.0,controls(),|| false).is_err());
    let other = RetreatGeometry::sphere(1000.0).unwrap();
    let wrong = surface(1,&phase,&binding,chart,400.0,other,&invert);
    assert!(advance_chamber_adaptive_moving(&[wrong],&geometry(other),None,&gas,&catalog,0.0,controls(),|| false).is_err());
    let empty = advance_chamber_adaptive_moving(&[row],&geom,None,&gas,&catalog,0.0,controls(),|| false).unwrap();
    assert_eq!(empty.gas(),&gas);
    assert_eq!(empty.geometry()[&1].initial,empty.geometry()[&1].final_geometry);
}
