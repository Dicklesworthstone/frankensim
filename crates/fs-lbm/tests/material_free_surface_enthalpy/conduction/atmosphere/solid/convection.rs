//! Native finite-solid convection: an independent integral/film heat oracle.
use super::*;
use fs_evidence::{Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, StatisticalCertificate};
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceBodyConvectionQuery, FreeSurfaceBodyGasConvection,
};
use fs_lbm::material::free_surface::enthalpy::conduction::surface::convection::FreeSurfaceHeatTransferCoefficient;
use fs_material::thermal::conduction::LayerConductionConfig;
use std::{cell::RefCell, sync::Arc};

#[path = "convection/radiation.rs"]
mod coupled_radiation;

#[derive(Debug)]
struct Film {
    h: f64,
    maximum_bulk_k: Option<f64>,
    calls: RefCell<Vec<[f64; 3]>>,
}
impl Film {
    fn new(h: f64) -> Self {
        Self { h, maximum_bulk_k: None, calls: RefCell::new(Vec::new()) }
    }
}
impl FreeSurfaceBodyGasConvection for Film {
    fn coefficient(&self, query: FreeSurfaceBodyConvectionQuery<'_>)
        -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        let bulk = query.bulk.state();
        let film = query.film.state();
        assert_eq!(query.body_id, 42);
        assert_eq!(query.face.wall, [-1, 0, 3]);
        assert_eq!(query.outward_normal, [1.0, 0.0, 0.0]);
        assert_eq!(bulk.pressure.to_bits(), film.pressure.to_bits());
        near(film.temperature, 0.5 * query.surface_temperature_k + 0.5 * bulk.temperature, 1e-13);
        self.calls.borrow_mut().push([query.surface_temperature_k, bulk.temperature, film.pressure]);
        if self.maximum_bulk_k.is_some_and(|limit| bulk.temperature > limit) {
            return Err(FreeSurfaceEnthalpyError::Convection("manufactured late finite-gas refusal".into()));
        }
        if query.surface_temperature_k == bulk.temperature { return Ok(None); }
        Ok(Some(Evidence {
            value: FreeSurfaceHeatTransferCoefficient::new(self.h),
            qoi: self.h,
            numerical: NumericalCertificate::exact(self.h),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: Default::default(),
            provenance: ProvenanceHash(0x534f_4c49_4446_494c),
            adjoint_ref: None,
        }))
    }
}
fn solve() -> LayerConductionConfig {
    LayerConductionConfig { max_iterations: 96, absolute_heat_tolerance_w: 1e-12, relative_heat_tolerance: 1e-12 }
}
fn face() -> FreeSurfaceBodyGasFace {
    FreeSurfaceBodyGasFace { distances_m: [0.0005, 0.0], ..solid_face() }
}
// Solve (scale / 2d) [(Tc-200)^2-(Ts-200)^2] = h(Ts-Tg)
// independently of the product's inverse-integral root, in either direction.
fn surface(cell: f64, gas: f64, scale: f64, h: f64) -> f64 {
    let b = h * 0.0005 / scale;
    let c = (cell - 200.0).powi(2) + 2.0 * b * (gas - 200.0);
    200.0 + c / ((b * b + c).sqrt() + b)
}

#[test]
fn integrated_solid_surface_closes_the_same_finite_gas_and_liquid_in_both_directions() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.01, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    for initial_gas in [350.0, 480.0] {
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let film = Film::new(5000.0);
            let mut gas = parcel(initial_gas, 1e-7, 100.0);
            let exposed = gas_contact(gas.model(), gas.point()).face;
            let wall = body::wall(&solid, &chart, &context, 450.0)
                .with_gas_faces(vec![face()])
                .with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap();
            let mut bodies = [wall.clone()];
            let report = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[exposed],
                gas_controls(), false, gate).unwrap().conduction;
            let (mut tl, mut ts, mut tg) = (400.0, 450.0, initial_gas);
            let (mut ql, mut qs) = (0.0, 0.0);
            let r = 8.314_462_618_153_24 / 0.028;
            let mass = 100000.0 * 1e-7 / (r * initial_gas);
            for _ in 0..4 {
                let kg = 2.6 * (0.5 + (100000.0 * tg / initial_gas - 1000.0) / 249000.0);
                let liquid_heat = 0.00025 * 1e-6 / (0.0005 / (tl - 200.0) + 0.0005 / kg) * (tg - tl);
                let solid_heat = 0.00025 * 1e-6 * film.h * (tg - surface(ts, tg, 0.01, film.h));
                tl += liquid_heat / 0.002;
                ts += solid_heat / 0.001;
                tg = (((2.0_f64 + 0.001 * tg).powi(2) - 0.002 * (liquid_heat + solid_heat) / (mass * r)).sqrt() - 2.0) / 0.001;
                ql += liquid_heat;
                qs += solid_heat;
            }
            assert_eq!(report.thermal.substeps, 4);
            let exchange = report.body_environment.as_ref().unwrap();
            let endpoint = &exchange.thermal_endpoint[0];
            let boundary = endpoint.convection.as_ref().unwrap();
            near(gas.state().state().temperature, tg, 2e-10);
            near(endpoint.solid_temperature_k, ts, 2e-11);
            near(boundary.surface_temperature_k, surface(ts, tg, 0.01, film.h), 2e-11);
            near(endpoint.gas.state().temperature, 0.5 * boundary.surface_temperature_k + 0.5 * tg, 2e-11);
            near(exchange.heat_to_solids_j, qs, 2e-10);
            near(report.environment.as_ref().unwrap().heat_to_liquid_j, ql, 2e-10);
            near(report.bodies.as_ref().unwrap().heat_from_gas_j, qs, 2e-10);
            near(boundary.solid_power_to_surface_w + endpoint.heat_into_solid_w, boundary.surface_heat_residual_w, 1e-13);
            assert!(report.finite_gas.unwrap().combined_thermal_balance_defect_j.abs() < 1e-9);
            assert!(Arc::ptr_eq(&boundary.source_line, &exchange.last_applied[0].convection.as_ref().unwrap().source_line));
            let source = Arc::clone(&boundary.source_line);
            let again = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[exposed],
                gas_controls(), false, gate).unwrap().conduction;
            assert!(Arc::ptr_eq(&source, &again.body_environment.unwrap().thermal_endpoint[0].convection.as_ref().unwrap().source_line));
            assert!(film.calls.borrow().iter().any(|q| q[2] != 100000.0));
        });
    }
}

#[test]
fn convection_coefficient_changes_native_motion_through_the_shared_gas() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.01, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    let run = |h| with_exposed_liquid(&liquid, true, |thermal, gate| {
        let film = Film::new(h);
        let mut gas = parcel(350.0, 1e-7, 100.0);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
            .with_gas_faces(vec![face()])
            .with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap()];
        thermal.step_with_finite_gas_environment(
            &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[exposed],
            gas_controls(), false, gate).unwrap();
        (thermal.flow().lattice().grid().state_bits(), gas.state().state().temperature)
    });
    let weak = run(50.0);
    let strong = run(5000.0);
    assert_ne!(weak.0, strong.0);
    assert!(strong.1 > weak.1);
}

#[test]
fn geometry_surface_budget_and_late_provider_refusals_roll_back_every_owner() {
    let liquid = conducting_card(1.0);
    let solid = radiation::emissive_wall(1.0, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let film = Film::new(5000.0);
        let late = Film { maximum_bulk_k: Some(350.001), ..Film::new(5000.0) };
        let mut gas = parcel(350.0, 1e-7, 100.0);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let original = body::wall(&solid, &chart, &context, 450.0);
        assert!(original.clone().with_gas_convection(&film, [300.0, 501.0], solve(), gate).is_err());
        let mut exhausted = solve();
        exhausted.max_iterations = 1;
        exhausted.absolute_heat_tolerance_w = 0.0;
        exhausted.relative_heat_tolerance = 0.0;
        let mut unequal = radiation::radiant_face(450.0);
        unequal.area_m2 *= 0.5;
        let radiating = original.clone().with_radiation_faces(
            vec![unequal], MaterialPropertySelection::SingleClaimOnly).unwrap();
        let before = thermal.total_enthalpy_j().unwrap();
        let bits = thermal.flow().lattice().grid().state_bits();
        let gas_before = gas.state().clone();
        let bodies = [
            original.clone().with_gas_faces(vec![solid_face()]).with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap(),
            radiating.with_gas_faces(vec![face()]).with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap(),
            original.clone().with_gas_faces(vec![face()]).with_gas_convection(&film, [300.0, 500.0], exhausted, gate).unwrap(),
            original.clone().with_gas_faces(vec![face()]).with_gas_convection(&late, [300.0, 500.0], solve(), gate).unwrap(),
        ];
        for (index, body) in bodies.into_iter().enumerate() {
            let source = body.material().clone();
            let mut trial = [body];
            let error = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut trial, &[], &[], &mut gas, &[exposed],
                gas_controls(), false, gate).unwrap_err();
            if index == 2 { assert!(matches!(error, FreeSurfaceEnthalpyError::BodyConvectionBalance { .. })); }
            if index == 3 { assert!(matches!(error, FreeSurfaceEnthalpyError::Convection(_))); }
            assert_eq!(gas.state(), &gas_before);
            assert_eq!(trial[0].specific_enthalpy_j_kg(), original.specific_enthalpy_j_kg());
            assert_eq!(trial[0].material(), &source);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.accepted_steps(), 0);
        }
        assert!(late.calls.borrow().iter().any(|q| q[1] > 350.001));
    });
}

#[test]
fn exact_equilibrium_admits_no_coefficient_and_distinct_radiation_keeps_its_source() {
    let liquid = conducting_card(1.0);
    let solid = radiation::emissive_wall(1.0, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let film = Film::new(5000.0);
        let mut gas = parcel(450.0, 1e-7, 100.0);
        let mut other = radiation::radiant_face(450.0);
        other.wall[2] = 4;
        other.gas[2] = 4;
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
            .with_radiation_faces(vec![other], MaterialPropertySelection::SingleClaimOnly).unwrap()
            .with_gas_faces(vec![face()])
            .with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap()];
        let report = thermal.step_with_finite_gas_environment(
            &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[],
            gas_controls(), false, gate).unwrap().conduction;
        let exchange = report.body_environment.unwrap();
        let sample = &exchange.thermal_endpoint[0];
        assert_eq!(sample.heat_into_solid_w, 0.0);
        assert_eq!(sample.conductance_w_k, 0.0);
        let boundary = sample.convection.as_ref().unwrap();
        assert_eq!(boundary.coefficient, None);
        assert_eq!(boundary.surface_temperature_k, 450.0);
        assert_eq!(boundary.source_line.knots[0].1.receipt.selected,
            bodies[0].material().property("thermal_conductivity").unwrap().answer().receipt.selected);
        assert_eq!(report.body_radiation.unwrap().heat_to_solids_j, 0.0);
    });
}
