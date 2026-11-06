//! Source derivatives, finite legs and continuous homogeneous admission.

use super::*;
use crate::eos::peng_robinson::{
    PengRobinsonComponentV1, PengRobinsonMixtureV1, PengRobinsonPhaseV1,
};
use fs_qty::{Density, MolarMass, Pressure};

fn methane() -> PengRobinsonMixtureV1 {
    PengRobinsonMixtureV1::new(
        vec![PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
            molar_mass: MolarMass::new(0.016043),
        }],
        vec![1.0], vec![],
    ).unwrap()
}

fn at_volume(model: &PengRobinsonMixtureV1, t: f64, v: f64) -> PengRobinsonDensityCaloricEvaluationV1 {
    model.evaluate_density_caloric(
        Temperature::new(t), Density::new(model.molar_mass().value() / v),
    ).unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual / expected - 1.0).abs() <= tolerance,
        "actual={actual}, expected={expected}, relative tolerance={tolerance}");
}

#[test]
fn analytic_pressure_hessian_and_volume_energy_tangent_match_original_source() {
    let model = methane();
    let t = 300.0;
    for v in [0.0003, 0.001, 0.01] {
        let state = at_volume(&model, t, v);
        let (ptt, ptv) = state.pressure_temperature_second_derivatives();
        let ht = 0.001;
        let hv = v * 1e-5;
        let ptt_difference = (at_volume(&model, t + ht, v).pressure_derivatives().0
            - at_volume(&model, t - ht, v).pressure_derivatives().0) / (2.0 * ht);
        let ptv_difference = (at_volume(&model, t, v + hv).pressure_derivatives().0
            - at_volume(&model, t, v - hv).pressure_derivatives().0) / (2.0 * hv);
        let uv_difference = (at_volume(&model, t, v + hv).residual_internal_energy_j_mol()
            - at_volume(&model, t, v - hv).residual_internal_energy_j_mol()) / (2.0 * hv);
        relative(ptt, ptt_difference, 2e-7);
        relative(ptv, ptv_difference, 2e-7);
        relative(state.internal_energy_volume_derivative_pa(), uv_difference, 2e-7);
        relative(state.internal_energy_volume_derivative_pa(),
            t * state.pressure_derivatives().0 - state.pressure().value(), 1e-12);
    }
    let v = 1e16;
    let dilute = at_volume(&model, t, v);
    let dense = at_volume(&model, t, 0.001);
    let b = dense.covolume_m3_per_mol();
    let expected_factor = dense.internal_energy_volume_derivative_pa()
        * (0.001 * 0.001 + 2.0 * b * 0.001 - b * b);
    relative(dilute.internal_energy_volume_derivative_pa() * v * v, expected_factor, 1e-13);
    assert!(dilute.internal_energy_volume_derivative_pa() > 0.0);
    // Subtracting two ideal-sized terms loses this original residual.
    let cancelled = t * dilute.pressure_derivatives().0 - dilute.pressure().value();
    assert!((cancelled - dilute.internal_energy_volume_derivative_pa()).abs()
        > 0.1 * dilute.internal_energy_volume_derivative_pa());
}

#[test]
fn finite_coordinate_legs_retain_tiny_volume_work_under_finite_heating() {
    let model = methane();
    let rho0 = 10.0_f64;
    let rho1 = rho0.next_up();
    let t0 = 300.0;
    let t1 = 350.0;
    let reference = model.prepare_density_caloric_change(
        Temperature::new(t0), Density::new(rho0), Density::new(rho1),
    ).unwrap();
    let change = reference.evaluate(Temperature::new(t1)).unwrap();
    let (et, ev) = change.coordinate_internal_energy_changes_j_mol();
    let (st, sv) = change.coordinate_entropy_changes_j_mol_k();
    assert!(et.abs() > 1e12 * ev.abs() && ev != 0.0);
    assert!(st.abs() > 1e12 * sv.abs() && sv != 0.0);
    let dv = -model.molar_mass().value() * ((rho1 - rho0) / rho0) / rho1;
    let initial = reference.initial_state();
    let final_state = change.state();
    relative(ev, 0.5 * (initial.internal_energy_volume_derivative_pa()
        + final_state.internal_energy_volume_derivative_pa()) * dv, 2e-13);
    relative(sv + change.ideal_isochoric_entropy_density_correction_j_mol_k(),
        0.5 * (initial.pressure_derivatives().0 + final_state.pressure_derivatives().0) * dv,
        2e-13);
    assert_eq!(et + ev, change.residual_internal_energy_change_j_mol());
    let reverse = model.prepare_density_caloric_change(
        Temperature::new(t1), Density::new(rho1), Density::new(rho0),
    ).unwrap().evaluate(Temperature::new(t0)).unwrap();
    let (rt, rv) = reverse.coordinate_internal_energy_changes_j_mol();
    relative(rt, -et, 1e-13);
    relative(rv, -ev, 1e-13);
    let isochoric = model.prepare_density_caloric_change(
        Temperature::new(t0), Density::new(rho0), Density::new(rho0),
    ).unwrap().evaluate(Temperature::new(t1)).unwrap();
    assert_eq!(isochoric.coordinate_internal_energy_changes_j_mol().1, 0.0);
    assert_eq!(isochoric.coordinate_entropy_changes_j_mol_k().1, 0.0);
    assert_eq!(isochoric.ideal_isochoric_entropy_density_correction_j_mol_k(), 0.0);
}

#[test]
fn rectangle_admission_bounds_original_residual_cv_throughout_stable_region() {
    let model = methane();
    let reference = model.prepare_density_caloric_change(
        Temperature::new(300.0), Density::new(1.0), Density::new(50.0),
    ).unwrap();
    let bound = reference.admit_stable_rectangle(
        [Temperature::new(280.0), Temperature::new(450.0)],
    ).unwrap();
    for t in [280.0, 300.0, 350.0, 400.0, 450.0] {
        for rho in [1.0, 10.0, 20.0, 35.0, 50.0] {
            let state = model.evaluate_density_caloric(
                Temperature::new(t), Density::new(rho),
            ).unwrap();
            assert!(state.residual_cv_j_mol_k() >= bound);
            assert!(state.pressure().value() > 0.0);
            assert!(state.pressure_derivatives().1 < 0.0);
        }
    }
    let reversed = reference.admit_stable_rectangle(
        [Temperature::new(450.0), Temperature::new(280.0)],
    ).unwrap();
    assert_eq!(bound, reversed);
    assert!(reference.admit_stable_rectangle(
        [Temperature::new(310.0), Temperature::new(450.0)],
    ).is_err());
}

#[test]
fn stable_corners_do_not_admit_a_rectangle_crossing_the_spinodal() {
    let model = methane();
    let roots = model.evaluate_pt(Temperature::new(150.0), Pressure::new(1e6)).unwrap();
    let liquid = Density::new(roots.select(PengRobinsonPhaseV1::Liquid).mass_density_kg_per_m3());
    let vapor = Density::new(roots.select(PengRobinsonPhaseV1::Vapor).mass_density_kg_per_m3());
    assert!(liquid.value() > 5.0 * vapor.value());
    for t in [150.0, 151.0] {
        for density in [liquid, vapor] {
            assert!(model.evaluate_density_caloric(Temperature::new(t), density).is_ok());
        }
    }
    let reference = model.prepare_density_caloric_change(
        Temperature::new(150.0), liquid, vapor,
    ).unwrap();
    assert_eq!(reference.admit_stable_rectangle(
        [Temperature::new(150.0), Temperature::new(151.0)],
    ), Err(PengRobinsonErrorV1::MechanicallyUnstable));
}

#[test]
fn rectangle_refuses_nonsmooth_alpha_crossing_even_with_admissible_endpoints() {
    let model = methane();
    let reference = model.prepare_density_caloric_change(
        Temperature::new(300.0), Density::new(1.0), Density::new(2.0),
    ).unwrap();
    assert!(reference.evaluate(Temperature::new(3000.0)).is_ok());
    assert!(reference.admit_stable_rectangle(
        [Temperature::new(300.0), Temperature::new(3000.0)],
    ).is_err());
}
