//! Reversible narrow-beam potential oracles; synthetic sources, not measurements.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::ReferenceSurface;
use fs_material::laminate::thermal::beam::PreparedLaminateBeam;
use fs_material::laminate::thermal::beam::thermoelastic::{
    LaminateBeamThermoelasticError, ResolvedLaminateBeamThermoelastic,
};
use fs_material::laminate::thermal::source::{
    LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
};
use fs_material::laminate::thermal::thermoelastic::AffineLaminaThermoelastic;
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES as CTE3, ScalarAdmissibility,
    ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

const WIDTH: f64 = 0.02;
#[derive(Clone, Copy)]
struct Ply {
    h: f64,
    angle: f64,
    isotropic: bool,
    e: [[f64; 2]; 2],
    nu: [f64; 2],
    g: [f64; 2],
    alpha: [[f64; 2]; 2],
}
fn ply() -> Ply {
    Ply {
        h: 0.008,
        angle: 0.0,
        isotropic: true,
        e: [[70e9, -30e6], [35e9, 20e6]],
        nu: [0.28, 0.0004],
        g: [20e9, 0.0],
        alpha: [[2e-5, 1e-7], [3e-5, -5e-8]],
    }
}
fn quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", quantity(), t).unwrap()
}
fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(
        name,
        dims,
        if dims == Pressure::DIMS {
            ScalarAdmissibility::StrictlyPositive
        } else {
            ScalarAdmissibility::Finite
        },
    )
    .unwrap()
}
fn source(plies: &[Ply], reference: ReferenceSurface) -> LaminateThermalSource {
    let inputs = plies
        .iter()
        .enumerate()
        .map(|(layer, p)| {
            let mut values = if p.isotropic {
                vec![
                    ("E", Pressure::DIMS, p.e[0]),
                    ("nu", Dims::NONE, p.nu),
                    (CTE, INVERSE_TEMPERATURE_DIMS, p.alpha[0]),
                ]
            } else {
                vec![
                    ("young_modulus_1", Pressure::DIMS, p.e[0]),
                    ("young_modulus_2", Pressure::DIMS, p.e[1]),
                    ("poisson_ratio_12", Dims::NONE, p.nu),
                    ("shear_modulus_12", Pressure::DIMS, p.g),
                    ("shear_modulus_31", Pressure::DIMS, p.g),
                    ("shear_modulus_23", Pressure::DIMS, p.g),
                    (CTE3[0], INVERSE_TEMPERATURE_DIMS, p.alpha[0]),
                    (CTE3[1], INVERSE_TEMPERATURE_DIMS, p.alpha[1]),
                    (CTE3[2], INVERSE_TEMPERATURE_DIMS, [0.0; 2]),
                ]
            };
            values.push((
                DENSITY_PROPERTY,
                Density::DIMS,
                [2700.0 + layer as f64 * 1000.0, 0.0],
            ));
            let mut claims = ClaimSet::new();
            for (name, dims, [value, slope]) in values {
                let density = name == DENSITY_PROPERTY;
                claims
                    .insert_claim(PropertyClaim {
                        key: PropertyKey::new(name, dims),
                        value: if density {
                            PropertyValue::Scalar { value, dims }
                        } else {
                            PropertyValue::Curve {
                                abscissa: "T".into(),
                                abscissa_dims: Temperature::DIMS,
                                knots: vec![
                                    (280.0, value - 20.0 * slope),
                                    (300.0, value),
                                    (360.0, value + 60.0 * slope),
                                ],
                                dims,
                            }
                        },
                        interpolation: if density {
                            InterpolationPolicy::ConstantWithinValidity
                        } else {
                            InterpolationPolicy::LinearInside
                        },
                        validity: ValidityDomain::unconstrained().with_quantity(
                            "T",
                            quantity(),
                            if density { 300.0 } else { 280.0 },
                            if density { 300.0 } else { 360.0 },
                        ),
                        uncertainty: UncertaintyModel::Unstated,
                        observations: vec![],
                        provenance: Provenance {
                            source: "authored free-sided beam potential".into(),
                            license: "test-only".into(),
                            artifact: None,
                        },
                    })
                    .unwrap();
            }
            LaminaSourceInput {
                card: Arc::new(
                    MaterialCard::assemble(
                        MaterialStateId {
                            chemistry: format!("beam layer {layer}"),
                            phase: "solid".into(),
                            process: "numerical thermodynamic oracle".into(),
                            revision: 0,
                        },
                        claims,
                        vec![],
                    )
                    .unwrap(),
                ),
                reference: point(300.0),
                thickness_m: p.h,
                orientation_rad: p.angle,
                elasticity: if p.isotropic {
                    LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungPoisson {
                        young_modulus: requirement("E", Pressure::DIMS),
                        poisson_ratio: requirement("nu", Dims::NONE),
                    })
                } else {
                    LaminaElasticSource::OrthotropicPlane12
                },
                expansion: if p.isotropic {
                    LaminaExpansionSource::Isotropic
                } else {
                    LaminaExpansionSource::Orthotropic
                },
                selection: MaterialPropertySelection::SingleClaimOnly,
            }
        })
        .collect();
    LaminateThermalSource::try_new(inputs, reference, || false).unwrap()
}
fn laws(source: &LaminateThermalSource) -> Vec<AffineLaminaThermoelastic> {
    source
        .layers()
        .iter()
        .map(|p| AffineLaminaThermoelastic::new(p, [280.0, 360.0], 0.02).unwrap())
        .collect()
}
fn prepare(
    source: &LaminateThermalSource,
    laws: &[AffineLaminaThermoelastic],
    t: &[f64],
    axis: f64,
) -> ResolvedLaminateBeamThermoelastic {
    let current = source.at_temperatures(t, || false).unwrap();
    let beam = PreparedLaminateBeam::try_new(Arc::clone(current.prepared()), WIDTH, axis, || false)
        .unwrap();
    ResolvedLaminateBeamThermoelastic::try_new(laws, &current, &beam, || false).unwrap()
}
fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_relaxed_isotropic_capacity_is_the_uniaxial_potential_not_fixed_lateral_strain() {
    let mut p = ply();
    p.e[0][1] = 0.0;
    p.nu[1] = 0.0;
    p.alpha[0][1] = 0.0;
    // An explicitly isotropic family stays planar at arbitrary fixed orientation.
    p.angle = 0.71;
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let prepared = prepare(&source, &laws, &[330.0], 0.0);
    let state = prepared.state(0.002, 0.0, || false).unwrap();
    let (e, alpha, volume, t) = (p.e[0][0], p.alpha[0][0], WIDTH * p.h, 330.0);
    let elastic = 0.002 - alpha * (t - 300.0);
    close(
        state.free_energy_j_m,
        0.5 * volume * e * elastic * elastic,
        1e-11,
        2e-13,
    );
    close(
        state.entropy_j_m_k,
        volume * e * alpha * elastic,
        1e-12,
        2e-13,
    );
    close(
        state.heat_capacity_j_m_k,
        -t * volume * e * alpha * alpha,
        1e-11,
        2e-13,
    );
    close(
        state.generalized_force_temperature_derivative[0],
        -volume * e * alpha,
        1e-10,
        2e-13,
    );
    close(
        state.internal_energy_j_m,
        state.free_energy_j_m + t * state.entropy_j_m_k,
        1e-11,
        2e-13,
    );
    let fixed_lateral_capacity = -2.0 * t * volume * e * alpha * alpha / (1.0 - p.nu[0]);
    assert!((fixed_lateral_capacity - state.heat_capacity_j_m_k).abs() > 1.0);
    assert!(prepared.layer_minimum_heat_capacity_j_m_k(0).unwrap() <= state.heat_capacity_j_m_k);
}

#[test]
fn g1_variable_source_potential_derivatives_and_maxwell_reciprocity_include_relaxation() {
    // 90-degree orthotropy exercises material-axis rotation and distinct alpha2.
    let mut p = ply();
    p.isotropic = false;
    p.angle = core::f64::consts::FRAC_PI_2;
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let (t, eps, kappa, axis, dt) = (332.0, 0.003, 0.06, 0.001, 0.025);
    let prepared = prepare(&source, &laws, &[t], axis);
    let at = prepared.state(eps, kappa, || false).unwrap();
    let low = prepare(&source, &laws, &[t - dt], axis)
        .state(eps, kappa, || false)
        .unwrap();
    let high = prepare(&source, &laws, &[t + dt], axis)
        .state(eps, kappa, || false)
        .unwrap();
    close(
        at.entropy_j_m_k,
        -(high.free_energy_j_m - low.free_energy_j_m) / (2.0 * dt),
        1e-8,
        2e-7,
    );
    close(
        at.heat_capacity_j_m_k,
        -t * (high.free_energy_j_m - 2.0 * at.free_energy_j_m + low.free_energy_j_m) / (dt * dt),
        2e-6,
        2e-6,
    );
    close(
        at.generalized_force_temperature_derivative[0],
        (high.axial_force_n - low.axial_force_n) / (2.0 * dt),
        1e-7,
        2e-7,
    );
    close(
        at.generalized_force_temperature_derivative[1],
        (high.bending_moment_n_m - low.bending_moment_n_m) / (2.0 * dt),
        1e-9,
        2e-7,
    );
    for (coordinate, delta) in [(0, 1e-6), (1, 1e-5)] {
        let mut lo = [eps, kappa];
        let mut hi = lo;
        lo[coordinate] -= delta;
        hi[coordinate] += delta;
        let entropy_lo = prepared
            .state(lo[0], lo[1], || false)
            .unwrap()
            .entropy_j_m_k;
        let entropy_hi = prepared
            .state(hi[0], hi[1], || false)
            .unwrap()
            .entropy_j_m_k;
        close(
            (entropy_hi - entropy_lo) / (2.0 * delta),
            -at.generalized_force_temperature_derivative[coordinate],
            1e-7,
            2e-9,
        );
    }
    // Exact original native strain energy remains the same mechanical potential.
    let native = prepared.beam().evaluate(eps, kappa, || false).unwrap();
    close(at.free_energy_j_m, native.strain_energy_j_m, 1e-11, 2e-13);
    close(at.axial_force_n, native.axial_force_n, 1e-10, 2e-13);
    close(
        at.bending_moment_n_m,
        native.bending_moment_n_m,
        1e-11,
        2e-13,
    );
}

#[test]
fn g3_bimetal_layer_entropy_retains_temperature_weights_and_axis_covariance() {
    let first = ply();
    let mut second = ply();
    second.h = 0.004;
    second.e[0] = [110e9, 15e6];
    second.alpha[0] = [4e-5, -1e-7];
    let source = source(&[first, second], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let prepared = prepare(&source, &laws, &[315.0, 345.0], 0.0);
    let [eps, kappa] = prepared.beam().free_generalized_strain();
    let state = prepared.state(eps, kappa, || false).unwrap();
    close(state.axial_force_n, 0.0, 1e-8, 0.0);
    close(state.bending_moment_n_m, 0.0, 1e-10, 0.0);
    assert!(state.free_energy_j_m > 0.0);
    let layers: Vec<_> = (0..2)
        .map(|i| prepared.layer_state(i, eps, kappa).unwrap())
        .collect();
    close(
        state.internal_energy_j_m,
        state.free_energy_j_m + 315.0 * layers[0].entropy_j_m_k + 345.0 * layers[1].entropy_j_m_k,
        1e-11,
        2e-13,
    );
    let delta = 0.002;
    let shifted = PreparedLaminateBeam::try_new(
        Arc::clone(prepared.source().prepared()),
        WIDTH,
        delta,
        || false,
    )
    .unwrap();
    let shifted =
        ResolvedLaminateBeamThermoelastic::try_new(&laws, prepared.source(), &shifted, || false)
            .unwrap();
    let other = shifted.state(eps - delta * kappa, kappa, || false).unwrap();
    close(other.free_energy_j_m, state.free_energy_j_m, 1e-11, 2e-13);
    close(other.entropy_j_m_k, state.entropy_j_m_k, 1e-12, 2e-13);
    close(
        other.heat_capacity_j_m_k,
        state.heat_capacity_j_m_k,
        1e-11,
        2e-13,
    );
    for (layer, original) in layers.iter().enumerate() {
        let moved = shifted
            .layer_state(layer, eps - delta * kappa, kappa)
            .unwrap();
        close(
            moved.generalized_force_temperature_derivative[1],
            original.generalized_force_temperature_derivative[1]
                + delta * original.generalized_force_temperature_derivative[0],
            1e-11,
            2e-13,
        );
        assert!(
            prepared.layer_minimum_heat_capacity_j_m_k(layer).unwrap()
                <= original.heat_capacity_j_m_k
        );
    }
}

#[test]
fn g4_accidental_off_axis_alignment_does_not_admit_a_temperature_path() {
    let mut p = ply();
    p.isotropic = false;
    p.angle = core::f64::consts::FRAC_PI_4;
    p.e = [[1e9, -1e6], [1e9, 1e6]];
    p.nu = [0.25, 0.0];
    p.g = [4e8, 0.0];
    p.alpha = [[1e-5, 0.0]; 2];
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let current = source.at_temperatures(&[300.0], || false).unwrap();
    let beam = PreparedLaminateBeam::try_new(Arc::clone(current.prepared()), WIDTH, 0.0, || false)
        .unwrap();
    assert!(ResolvedLaminateBeamThermoelastic::try_new(&laws, &current, &beam, || false).is_err());
    let hot = source.at_temperatures(&[320.0], || false).unwrap();
    assert!(
        PreparedLaminateBeam::try_new(Arc::clone(hot.prepared()), WIDTH, 0.0, || false).is_err()
    );
}

#[test]
fn g4_exact_source_identity_domain_and_cancellation_are_required() {
    let source = source(&[ply()], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let prepared = prepare(&source, &laws, &[320.0], 0.0);
    let repeated = source.at_temperatures(&[320.0], || false).unwrap();
    assert_eq!(repeated.prepared().as_ref(), prepared.beam().source());
    assert!(
        ResolvedLaminateBeamThermoelastic::try_new(&laws, &repeated, prepared.beam(), || false)
            .is_err()
    );
    assert!(
        ResolvedLaminateBeamThermoelastic::try_new(
            &laws,
            &prepared.source().clone(),
            prepared.beam(),
            || false
        )
        .is_ok()
    );
    let foreign = super_source(&[ply()]);
    assert!(
        ResolvedLaminateBeamThermoelastic::try_new(
            &foreign,
            prepared.source(),
            prepared.beam(),
            || false
        )
        .is_err()
    );
    assert!(matches!(
        ResolvedLaminateBeamThermoelastic::try_new(
            &laws,
            prepared.source(),
            prepared.beam(),
            || true
        ),
        Err(LaminateBeamThermoelasticError::Cancelled)
    ));
    assert!(matches!(
        prepared.state(0.0, 0.0, || true),
        Err(LaminateBeamThermoelasticError::Cancelled)
    ));
    assert!(prepared.layer_state(1, 0.0, 0.0).is_err());
    assert!(prepared.layer_state(0, 0.03, 0.0).is_err());
    assert!(prepared.layer_state(0, f64::NAN, 0.0).is_err());
    assert!(laws[0].at_temperature(360.1).is_err());
    assert!(prepared.state(0.002, 0.0, || false).is_ok());
}
fn super_source(plies: &[Ply]) -> Vec<AffineLaminaThermoelastic> {
    laws(&source(plies, ReferenceSurface::Midsurface))
}
