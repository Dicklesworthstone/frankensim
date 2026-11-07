//! One conducting surface, separate gas/radiant powers, and finite native stores.
use super::*;

const SIGMA: f64 = 5.670_374_419e-8;
const AREA: f64 = 1e-6;

fn emissivity(t: f64) -> f64 {
    0.2 + 0.003 * (t - 300.0)
}

// Independent temperature bisection with the explicit affine-k primitive and
// fourth-power radiation. Product code solves in the material power coordinate.
fn shared_surface(cell: f64, gas: f64, radiant: f64, h: f64) -> f64 {
    let (mut low, mut high) = (300.0_f64, 500.0_f64);
    for _ in 0..100 {
        let t = 0.5 * low + 0.5 * high;
        let solid = AREA * 0.25 / 0.0005
            * ((cell - 200.0).powi(2) - (t - 200.0).powi(2));
        let gas_power = AREA * h * (t - gas);
        let radiation = AREA * emissivity(t) * SIGMA * (t.powi(4) - radiant.powi(4));
        let residual = solid - gas_power - radiation;
        if residual == 0.0 { return t; }
        if residual > 0.0 { low = t; } else { high = t; }
    }
    0.5 * low + 0.5 * high
}

#[test]
fn one_surface_closes_independent_gas_radiant_and_three_store_energy_in_both_directions() {
    let liquid = conducting_card(1.0);
    let solid = radiation::emissive_wall(1.0, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    for (initial_gas, surroundings) in [(350.0, 900.0), (480.0, 300.0), (450.0, 1000.0)] {
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let film = Film::new(5000.0);
            let mut gas = parcel(initial_gas, 1e-7, 100.0);
            let exposed = gas_contact(gas.model(), gas.point()).face;
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
                .with_radiation_faces(vec![radiation::radiant_face(surroundings)],
                    MaterialPropertySelection::SingleClaimOnly).unwrap()
                .with_gas_faces(vec![face()])
                .with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap()];
            let before_liquid = thermal.total_enthalpy_j().unwrap();
            let before_solid = bodies[0].enthalpy_j();
            let before_gas = gas.state().specific_internal_energy_j_kg();
            let report = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[exposed],
                gas_controls(), false, gate).unwrap().conduction;

            let (mut tl, mut tc, mut tg) = (400.0, 450.0, initial_gas);
            let (mut ql_total, mut qg_total, mut qr_total) = (0.0, 0.0, 0.0);
            let r = 8.314_462_618_153_24 / 0.028;
            let mass = 100000.0 * 1e-7 / (r * initial_gas);
            for _ in 0..4 {
                let kg = 2.6 * (0.5 + (100000.0 * tg / initial_gas - 1000.0) / 249000.0);
                let ql = 0.00025 * AREA / (0.0005 / (tl - 200.0) + 0.0005 / kg) * (tg - tl);
                let ts = shared_surface(tc, tg, surroundings, film.h);
                let qg = 0.00025 * AREA * film.h * (tg - ts);
                let qr = 0.00025 * AREA * emissivity(ts) * SIGMA
                    * (surroundings.powi(4) - ts.powi(4));
                tl += ql / 0.002;
                tc += (qg + qr) / 0.001;
                tg = (((2.0_f64 + 0.001 * tg).powi(2)
                    - 0.002 * (ql + qg) / (mass * r)).sqrt() - 2.0) / 0.001;
                ql_total += ql;
                qg_total += qg;
                qr_total += qr;
            }
            assert_eq!(report.thermal.substeps, 4);
            let exchange = report.body_environment.as_ref().unwrap();
            let sample = &exchange.thermal_endpoint[0];
            let boundary = sample.convection.as_ref().unwrap();
            let radiation = &report.body_radiation.as_ref().unwrap().thermal_endpoint[0];
            let ts = shared_surface(tc, tg, surroundings, film.h);
            near(sample.solid_temperature_k, tc, 2e-10);
            near(gas.state().state().temperature, tg, 2e-10);
            near(boundary.surface_temperature_k, ts, 2e-10);
            assert_eq!(boundary.radiation.as_ref().unwrap(), radiation);
            assert_eq!(radiation.surface_temperature_k, boundary.surface_temperature_k);
            near(radiation.emissivity, emissivity(ts), 2e-12);
            near(sample.gas.state().temperature, 0.5 * (ts + tg), 2e-12);
            near(boundary.solid_power_to_surface_w + sample.heat_into_solid_w
                + radiation.heat_into_solid_w, boundary.surface_heat_residual_w, 1e-13);
            near(exchange.heat_to_solids_j, qg_total, 2e-10);
            near(report.body_radiation.as_ref().unwrap().heat_to_solids_j, qr_total, 2e-10);
            near(report.environment.as_ref().unwrap().heat_to_liquid_j, ql_total, 2e-10);
            near(bodies[0].enthalpy_j() - before_solid, qg_total + qr_total, 2e-10);
            let gas_delta = gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - before_gas);
            near(gas_delta + qg_total + ql_total, 0.0, 2e-9);
            near(thermal.total_enthalpy_j().unwrap() - before_liquid
                + bodies[0].enthalpy_j() - before_solid + gas_delta, qr_total, 2e-9);
            assert!(report.finite_gas.as_ref().unwrap().combined_thermal_balance_defect_j.abs() < 1e-9);

            // Independent elimination of a three-leg passive star. Gas capacity
            // must include the radiant relay even when cell and gas T coincide.
            let material_g = AREA / 0.0005 * 0.25 * (tc + ts - 400.0);
            let gas_g = AREA * film.h;
            let radiant_g = AREA * emissivity(ts) * SIGMA
                * (ts + surroundings) * (ts * ts + surroundings * surroundings);
            let total = material_g + gas_g + radiant_g;
            near(sample.conductance_w_k, material_g * gas_g / total, 2e-10);
            near(radiation.conductance_w_k, material_g * radiant_g / total, 2e-10);
            near(boundary.gas_radiation_conductance_w_k, gas_g * radiant_g / total, 2e-10);
            near(sample.gas_decay_conductance_w_k,
                sample.conductance_w_k + boundary.gas_radiation_conductance_w_k, 1e-13);
            assert!(sample.gas_decay_conductance_w_k > sample.conductance_w_k);
            if initial_gas == 450.0 {
                assert!(boundary.surface_temperature_k > sample.solid_temperature_k);
                assert!(sample.heat_into_solid_w < 0.0);
                assert!(radiation.heat_into_solid_w > 0.0);
            } else {
                assert_ne!(sample.heat_into_solid_w.signum(), radiation.heat_into_solid_w.signum());
            }
        });
    }
}

#[test]
fn emitting_surface_source_and_late_gas_refusals_roll_back_the_shared_transaction() {
    let liquid = conducting_card(1.0);
    for limited_source in [true, false] {
        let solid = radiation::emissive_wall(1.0, if limited_source { 450.0001 } else { 500.0 });
        let chart = body::wall_chart(&solid);
        let context = point(450.0);
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let film = Film { maximum_bulk_k: if limited_source { None } else { Some(450.0001) },
                ..Film::new(5000.0) };
            let mut gas = parcel(450.0, 1e-7, 100.0);
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
                .with_radiation_faces(vec![radiation::radiant_face(1000.0)],
                    MaterialPropertySelection::SingleClaimOnly).unwrap()
                .with_gas_faces(vec![face()])
                .with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap()];
            let before = thermal.total_enthalpy_j().unwrap();
            let populations = thermal.flow().lattice().grid().state_bits();
            let gas_before = gas.state().clone();
            let wall_before = bodies[0].specific_enthalpy_j_kg();
            let source_before = bodies[0].material().clone();
            let error = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[],
                gas_controls(), false, gate).unwrap_err();
            if !limited_source {
                assert!(matches!(error, FreeSurfaceEnthalpyError::Convection(_)));
                assert!(film.calls.borrow().iter().any(|q| q[1] > 450.0001));
            }
            assert_eq!(gas.state(), &gas_before);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
            assert_eq!(bodies[0].material(), &source_before);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
            assert_eq!(thermal.flow().lattice().grid().state_bits(), populations);
            assert_eq!(thermal.accepted_steps(), 0);
        });
    }
}

#[test]
fn matching_radiation_and_natural_equilibrium_do_not_create_heat_or_gas_decay() {
    let liquid = conducting_card(1.0);
    let solid = radiation::emissive_wall(1.0, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let film = Film::new(5000.0);
        let mut gas = parcel(450.0, 1e-7, 100.0);
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)
            .with_radiation_faces(vec![radiation::radiant_face(450.0)],
                MaterialPropertySelection::SingleClaimOnly).unwrap()
            .with_gas_faces(vec![face()])
            .with_gas_convection(&film, [300.0, 500.0], solve(), gate).unwrap()];
        let report = thermal.step_with_finite_gas_environment(
            &[], &[], controls(), &mut bodies, &[], &[], &mut gas, &[],
            gas_controls(), false, gate).unwrap().conduction;
        let sample = &report.body_environment.as_ref().unwrap().thermal_endpoint[0];
        let boundary = sample.convection.as_ref().unwrap();
        assert_eq!(sample.heat_into_solid_w, 0.0);
        assert_eq!(sample.gas_decay_conductance_w_k, 0.0);
        assert_eq!(boundary.surface_temperature_k, 450.0);
        assert_eq!(boundary.coefficient, None);
        assert_eq!(boundary.solid_power_to_surface_w, 0.0);
        assert_eq!(boundary.radiation.as_ref().unwrap().heat_into_solid_w, 0.0);
        assert_eq!(report.body_radiation.unwrap().heat_to_solids_j, 0.0);
    });
}
