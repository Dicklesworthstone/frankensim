//! Source-local bulk viscosity on the original held-stress thermal interval.

use super::*;
use crate::d3q19::CollisionModel3;
use crate::material::free_surface::enthalpy::FreeSurfaceCellHeat;
use crate::material::free_surface::enthalpy::conduction::FreeSurfaceConductionControls;
use std::collections::BTreeMap;

const DT: f64 = 0.0005;

// Independent dimensional closure. At this dx/dt/rho, tau_s=1/2+1.5mu
// and tau_b=1/2+2.25zeta. The manufactured mu is affine in T and gamma,
// so its positive shear root is available without the production iteration.
fn response(temperature: f64, shear_numerator: f64, trace_numerator: f64) -> [f64; 4] {
    let a = 0.5 + 1.5 * viscosity_at(temperature, 0.0, 0.02);
    let rhs = shear_numerator / DT;
    let shear = 2.0 * rhs / (a + (a * a + 0.12 * rhs).sqrt());
    let bulk = sourced_bulk_at(temperature, shear);
    let divergence = trace_numerator / ((0.5 + 2.25 * bulk) * DT);
    let power = viscosity_at(temperature, shear, 0.02) * shear * shear
        + bulk * divergence * divergence;
    [shear, divergence, bulk, power]
}

fn integrated_heat(initial_t: f64, shear_numerator: f64, trace_numerator: f64) -> f64 {
    let evaluate = |temperature| {
        let power = response(temperature, shear_numerator, trace_numerator)[3];
        [40.0 / DT + power / 2_000_000.0, power]
    };
    let step = DT / 256.0;
    let mut temperature = initial_t;
    let mut heat_density = 0.0;
    for _ in 0..256 {
        let a = evaluate(temperature);
        let b = evaluate(temperature + 0.5 * step * a[0]);
        let c = evaluate(temperature + 0.5 * step * b[0]);
        let d = evaluate(temperature + step * c[0]);
        temperature += step * (a[0] + 2.0 * b[0] + 2.0 * c[0] + d[0]) / 6.0;
        heat_density += step * (a[1] + 2.0 * b[1] + 2.0 * c[1] + d[1]) / 6.0;
    }
    heat_density
}

#[test]
fn source_local_bulk_feedback_refreshes_trace_relaxation_with_the_same_thermal_sample() {
    let card = sourced_bulk_card(false);
    let law = sourced_bulk_law(&card);
    let surface = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
        law.reference_state().material().clone(), LiquidViscosityBasis::Dynamic, "T", "rate",
    ).unwrap();
    let caloric = chart(&card, 1000.0);
    let (points, tets, _) = cube([0.0; 3]);
    let h: Vec<_> = points.iter().copied().map(affine_h).collect();
    let velocities: Vec<_> = points.iter().map(|point| {
        let r = sub(*point, [0.004; 3]);
        [0.01 + 20.0 * r[1] + 5.0 * r[0], -0.004 + 5.0 * r[1], 0.002 + 5.0 * r[2]]
    }).collect();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 0.001, DT, 0.1)
            .unwrap().with_bulk_viscosity_source("bulk_eta", 1.2).unwrap();
        let (mut liquid, transfer) = scaling.free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
            domain(), [0.0; 3], source(&points, &tets, &velocities, &h),
            &law, &caloric, bounds(), Default::default(), solve(), cx, &gate,
        ).unwrap();
        let mut held = BTreeMap::new();
        let mut authored = Vec::new();
        let mut expected_heat = 0.0;
        for component in &transfer.mechanical.components {
            for &cell in &component.target_cells {
                let material = liquid.material_at(cell).unwrap().material();
                let temperature = material.state_coordinate("T").unwrap();
                let shear_numerator = 20.0 * (0.5 + 1.5 * viscosity_at(temperature, 20.0, 0.02)) * DT;
                let trace_numerator = 15.0 * (0.5 + 2.25 * sourced_bulk_at(temperature, 20.0)) * DT;
                let mass = liquid.flow().region_mass_kg(|x, y, z| [x, y, z] == cell).unwrap();
                expected_heat += integrated_heat(temperature, shear_numerator, trace_numerator) * mass / 1000.0;
                held.insert(cell, (shear_numerator, trace_numerator, temperature));
                authored.push(FreeSurfaceCellHeat { coordinates: cell, energy_j: mass * 2000.0 * 40.0 });
            }
        }
        let controls = FreeSurfaceConductionControls {
            maximum_faces: 1, maximum_substeps: 34, maximum_dt_s: DT / 32.0,
            energy_tolerance_j: 2e-9,
        };
        let report = liquid.step_with_conduction_and_viscous_feedback(
            &[], &authored, controls, &gate,
        ).unwrap();
        let feedback = report.viscous_feedback.as_ref().unwrap();
        assert!((feedback.heat_j / expected_heat - 1.0).abs() < 0.003);
        for sample in &feedback.thermal_endpoint {
            let (a, b, initial_t) = held[&sample.coordinates];
            let temperature = sample.material.material().state_coordinate("T").unwrap();
            let [shear, divergence, bulk, power] = response(temperature, a, b);
            close(sample.shear_rate_per_s, shear, 2e-8);
            close(sample.divergence_magnitude_per_s, divergence, 2e-8);
            close(sample.bulk_viscosity_pa_s, bulk, 2e-10);
            close(sample.bulk_collision_tau, 0.5 + 2.25 * bulk, 2e-10);
            close(sample.bulk_power_w_m3, bulk * divergence * divergence, 2e-7);
            close(sample.power_w_m3, power, 2e-7);
            assert!(temperature > initial_t + 39.9);
            assert!(divergence < 14.0);
            assert_eq!(sample.material.material().property("bulk_eta").unwrap().answer().receipt.selected,
                law.reference_state().material().property("bulk_eta").unwrap().answer().receipt.selected);
        }
        // Advection/conversion is a later source epoch. Every accepted wet cell,
        // including reconstruction ghosts, has its own admitted endpoint rate.
        for location in liquid.flow().lattice().mass_snapshot().keys() {
            let Some(cell) = FreeSurface3::transport_coordinates(*location) else { continue };
            let bulk = liquid.material_at(cell).unwrap().material().property("bulk_eta").unwrap().value_si();
            let CollisionModel3::CentralMomentBulk { bulk_rate, .. } = liquid.flow().lattice().collision_model_at(cell).unwrap() else { panic!("lost bulk collision") };
            close(bulk_rate, 1.0 / (0.5 + 2.25 * bulk), 2e-13);
        }
        close(report.flow.viscous_heat_j, feedback.heat_j, 0.0);
        close(report.flow.balance_defect_j, 0.0, 2e-9);
        let bits = liquid.flow().lattice().grid().state_bits();
        let enthalpy = liquid.total_enthalpy_j().unwrap();
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(liquid.step_with_conduction_and_viscous_feedback(&[], &[], controls, &cancelled).is_err());
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), enthalpy);
    });
}
