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

#[test]
fn mechanically_stable_rectangle_still_requires_positive_pressure_at_largest_volume() {
    let model = methane();
    let t = 150.0;
    // Obtain b from the original native EOS rather than copying a model constant.
    let b = at_volume(&model, t, 0.001).covolume_m3_per_mol();
    let initial = at_volume(&model, t, 1.4 * b);
    let final_density = Density::new(model.molar_mass().value() / (2.0 * b));
    let reference = model
        .prepare_density_caloric_change(initial.temperature(), initial.density(), final_density)
        .unwrap();

    // Recover this same state's attraction from its original pressure. On
    // 1.4 <= v/b <= 2 the spinodal threshold decreases to f(2)=49/6,
    // so alpha<f(2) admits mechanical stability on the COMPLETE interval.
    let v0 = initial.molar_volume_m3_per_mol();
    let d0 = v0 * v0 + 2.0 * b * v0 - b * b;
    let attraction = (R * t / (v0 - b) - initial.pressure().value()) * d0;
    let alpha = attraction / (R * t * b);
    assert!(alpha > 8.1 && alpha < 49.0 / 6.0);
    assert!(initial.pressure().value() > 0.0);
    assert!(initial.pressure_derivatives().1 < 0.0);

    // Despite that negative pressure slope everywhere, p crosses zero:
    // about +15.57 MPa at 1.4b and -7.59 MPa at 2b.
    let v1 = model.molar_mass().value() / final_density.value();
    let d1 = v1 * v1 + 2.0 * b * v1 - b * b;
    let final_pressure = R * t / (v1 - b) - attraction / d1;
    let final_slope =
        -R * t / ((v1 - b) * (v1 - b)) + attraction * (2.0 * v1 + 2.0 * b) / (d1 * d1);
    assert!(final_pressure < -7.0e6);
    assert!(final_slope < 0.0);
    assert_eq!(
        reference.admit_stable_rectangle([Temperature::new(t), Temperature::new(t)]),
        Err(PengRobinsonErrorV1::Unrepresentable {
            what: "rectangle cannot guarantee positive pressure",
        }),
    );
}

#[test]
fn negative_residual_cv_uses_the_cold_dense_corner_without_rejecting_the_rectangle() {
    let model = PengRobinsonMixtureV1::new(
        vec![PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(5.1953),
            critical_pressure: Pressure::new(227460.0),
            acentric_factor: -0.385,
            molar_mass: MolarMass::new(0.004002602),
        }],
        vec![1.0],
        vec![],
    )
    .unwrap();
    let reference = model
        .prepare_density_caloric_change(
            Temperature::new(300.0),
            Density::new(1.0),
            Density::new(50.0),
        )
        .unwrap();
    let bound = reference
        .admit_stable_rectangle([Temperature::new(300.0), Temperature::new(500.0)])
        .unwrap();
    let cold_dense_cv = model
        .evaluate_density_caloric(Temperature::new(300.0), Density::new(50.0))
        .unwrap()
        .residual_cv_j_mol_k();
    // Here kappa=-0.259138992 and B<0. The minimum switches to the
    // LOW-temperature, HIGH-density corner; residual Cv need not be positive.
    assert!(bound < 0.0 && bound <= cold_dense_cv);
    assert!((cold_dense_cv + 0.09778055534009236).abs() < 1e-12);
    assert!((bound - cold_dense_cv).abs() < 1e-12);
    // An explicitly supplied positive ideal Cv still admits total calorics.
    assert!(1.5 * R + bound > 0.0);
    for t in [300.0, 350.0, 400.0, 450.0, 500.0] {
        for rho in [1.0, 10.0, 20.0, 35.0, 50.0] {
            let state = model
                .evaluate_density_caloric(Temperature::new(t), Density::new(rho))
                .unwrap();
            assert!(state.residual_cv_j_mol_k() < 0.0);
            assert!(state.residual_cv_j_mol_k() >= bound);
            assert!(state.pressure().value() > 0.0);
            assert!(state.pressure_derivatives().1 < 0.0);
        }
    }
}
