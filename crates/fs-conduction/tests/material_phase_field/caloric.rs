//! Original Cp(T) reaches the assembled phase transport flux and Jacobian.

use super::*;
use fs_material::phase::HeatCapacityKnot;

const COLD: f64 = 320.0;
const HOT: f64 = 580.0;

fn source_curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    // Cp(T)=T, hence h=(T^2-300^2)/2 on the original source interval.
    let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        0.0,
        &[300.0, 600.0].map(|temperature_k| HeatCapacityKnot {
            temperature_k,
            specific_heat_capacity_j_kg_k: temperature_k,
            bulk_density_kg_m3: 1000.0,
        }),
        1e6,
    )
    .unwrap();
    assert_eq!(curve.knots().len(), 2);
    curve
}

fn original_enthalpy(temperature: f64) -> f64 {
    0.5 * (temperature * temperature - 300.0 * 300.0)
}

// One hot vertex has coordinate x with density 3*(1-x)^2. Set
// u=T/COLD, r=HOT/COLD and x=(u^2-1)/(r^2-1). These independent
// polynomial primitives integrate k=T/300 and lambda_j*dk/dh exactly;
// each cold vertex has conditional coordinate (1-x)/3.
fn reference_moments() -> (f64, [f64; 4]) {
    let r = HOT / COLD;
    let r2 = r * r;
    let r4 = r2 * r2;
    let span = r2 - 1.0;
    let value_primitive =
        |u: f64| r4 * u.powi(3) / 3.0 - 2.0 * r2 * u.powi(5) / 5.0 + u.powi(7) / 7.0;
    let all_primitive = |u: f64| r4 * u - 2.0 * r2 * u.powi(3) / 3.0 + u.powi(5) / 5.0;
    let hot_primitive = |u: f64| {
        u.powi(7) / 7.0 - (2.0 * r2 + 1.0) * u.powi(5) / 5.0 + (r4 + 2.0 * r2) * u.powi(3) / 3.0
            - r4 * u
    };
    let mean = (COLD / 300.0) * 6.0 / span.powi(3) * (value_primitive(r) - value_primitive(1.0));
    let all = 6.0 / (300.0 * COLD * span.powi(3)) * (all_primitive(r) - all_primitive(1.0));
    let hot = 6.0 / (300.0 * COLD * span.powi(4)) * (hot_primitive(r) - hot_primitive(1.0));
    let cold = (all - hot) / 3.0;
    assert!(mean > COLD / 300.0 && mean < HOT / 300.0);
    assert!(cold > 0.0 && hot > 0.0);
    (mean, [cold, hot, cold, cold])
}

#[test]
fn original_caloric_law_drives_source_integrated_element_flux_and_assembled_tangent() {
    let mesh = mesh();
    let card = conductivity_card(PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        // The interior source knot must survive the quadratic h(T) mapping.
        knots: vec![(300.0, 1.0), (450.0, 1.5), (600.0, 2.0)],
        dims: KD,
    });
    let curve = source_curve(&card);
    let source = prepare(&card, &curve, [0.0, original_enthalpy(600.0)]);
    assert_eq!(
        source.knots(),
        &[[0.0, 1.0], [56_250.0, 1.5], [135_000.0, 2.0]]
    );
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let h = [COLD, HOT, COLD, COLD].map(original_enthalpy);
    let (mean_k, conductivity_derivatives) = reference_moments();
    let jump = HOT - COLD;
    let dt = 0.1;
    let zero_heat = ScalarField::Uniform(0.0);
    // The consistent P1 source load is (sum(q)+q_i)/120. This zero-net
    // manufactured source balances the independent physical conduction load,
    // so the actual nonlinear/adjoint owner exposes its Jacobian at this field.
    let q = ScalarField::nodal(
        "balanced source-Cp transport",
        4,
        vec![-20.0 * mean_k * jump, 20.0 * mean_k * jump, 0.0, 0.0],
    )
    .unwrap();
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        assert!(std::ptr::eq(
            law.conductivity_source().caloric_curve(),
            &curve
        ));
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let flux = law.element_heat_flux(cx, &mesh, &h).unwrap();
        near(flux[0][0], -mean_k * jump, 2e-11);
        near(flux[0][1], 0.0, 1e-13);
        near(flux[0][2], 0.0, 1e-13);
        let endpoint = solver
            .evaluate_endpoint(
                cx,
                law.problem(&mesh, &boundary, &zero_heat),
                None,
                &h,
                &h,
                dt,
                config(),
                None,
            )
            .unwrap();
        let laplace_temperature = [-jump, jump, 0.0, 0.0];
        for i in 0..4 {
            near(
                endpoint.residual_j[i],
                dt * mean_k * laplace_temperature[i] / 6.0,
                2e-12,
            );
        }
        near(endpoint.energy_residual_j, 0.0, 1e-12);

        let response = solver
            .linearize_step(
                cx,
                law.problem(&mesh, &boundary, &q),
                None,
                &h,
                dt,
                config(),
            )
            .unwrap();
        for (actual, expected) in response.primal().specific_enthalpy_j_kg.iter().zip(h) {
            near(*actual, expected, 1e-7);
        }
        let laplace = [
            [3.0, -1.0, -1.0, -1.0],
            [-1.0, 1.0, 0.0, 0.0],
            [-1.0, 0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0, 1.0],
        ];
        let nodal_temperature = [COLD, HOT, COLD, COLD];
        for column in 0..4 {
            let mut direction = [0.0; 4];
            direction[column] = 1.0;
            let action = response.apply_jacobian(cx, &direction).unwrap();
            for row in 0..4 {
                let expected = if row == column { MASS } else { 0.0 }
                    + dt / 6.0
                        * (conductivity_derivatives[column] * laplace_temperature[row]
                            + mean_k * laplace[row][column] / nodal_temperature[column]);
                near(action[row], expected, 2e-12);
            }
        }
        near(response.primal().source_w, 0.0, 1e-12);
        near(response.primal().stored_energy_change_j, 0.0, 1e-9);
        near(response.primal().energy_residual_j, 0.0, 1e-10);
    });
}

#[test]
fn uniform_smooth_source_knot_keeps_its_classical_caloric_adjoint() {
    let mesh = mesh();
    let card = conductivity_card(PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        // Both original slopes are exactly 0.01 W/(m K^2). Reconstructing
        // dk/dh from separate h chords makes them differ by one f64 ulp.
        knots: vec![(300.0, 1.0), (400.0, 2.0), (600.0, 4.0)],
        dims: KD,
    });
    let curve = source_curve(&card);
    let source = prepare(&card, &curve, [0.0, original_enthalpy(600.0)]);
    assert!(source.knots().contains(&[35_000.0, 2.0]));
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let law = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let response = solver
            .linearize_step(
                cx,
                law.problem(&mesh, &boundary, &q),
                None,
                &[35_000.0; 4],
                0.1,
                config(),
            )
            .unwrap();
        for &temperature in &response.primal().temperature {
            near(temperature, 400.0, 1e-12);
        }
        near(response.primal().energy_residual_j, 0.0, 1e-12);
        let action = response.apply_jacobian(cx, &[1.0, 0.0, 0.0, 0.0]).unwrap();
        // Uniform T removes the d<k>/dh*grad(T) term. The remaining exact
        // source tangent is k/Cp=2/400 times the ordinary element Laplacian.
        for (row, laplace) in [3.0, -1.0, -1.0, -1.0].into_iter().enumerate() {
            let expected = if row == 0 { MASS } else { 0.0 } + 0.1 * (2.0 / 400.0) * laplace / 6.0;
            near(action[row], expected, 1e-13);
        }
    });
}
