//! Fixed spatial approximation: changing E/G must not silently change inertia.

use super::*;
use fs_material::section_temperature::ConcentricSectionTemperatureLaw;
use fs_qty::Temperature;

const RADII: [f64; 2] = [0.008, 0.02];
const YOUNG_COLD: [f64; 2] = [2.0e9, 7.0e9];
const SHEAR_COLD: [f64; 2] = [3.0e6, 9.0e6];
const ETA_COLD: [[f64; 3]; 2] = [[2.0e5, 4.0e6, 3.0e5], [7.0e5, 2.0e6, 1.0e5]];
const HOT_T: [f64; 2] = [350.0, 380.0];

fn area(layer: usize) -> f64 {
    core::f64::consts::PI
        * (RADII[layer].powi(2) - if layer == 0 { 0.0 } else { RADII[layer - 1].powi(2) })
}

fn inertia(layer: usize) -> f64 {
    core::f64::consts::PI / 4.0
        * (RADII[layer].powi(4) - if layer == 0 { 0.0 } else { RADII[layer - 1].powi(4) })
}

fn relative_temperature(temperature: f64) -> f64 {
    (temperature - 300.0) / 100.0
}

fn close(actual: f64, expected: f64) {
    let tolerance = 5.0e-12 * expected.abs().max(1.0e-12);
    assert!((actual - expected).abs() <= tolerance, "{actual:e} != {expected:e}");
}

fn sections() -> (FrameViscousSection, FrameViscousSection) {
    let mut cards = Vec::new();
    let mut layers = Vec::new();
    for layer in 0..2 {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, dims, cold, hot_factor) in [
            ("density", Density::DIMS, RHO, 1.0),
            ("young_modulus", Pressure::DIMS, YOUNG_COLD[layer], 0.6),
            ("shear_modulus", Pressure::DIMS, SHEAR_COLD[layer], 2.8),
            ("axial_eta", DynViscosity::DIMS, ETA_COLD[layer][0], 2.0),
            ("bending_eta", DynViscosity::DIMS, ETA_COLD[layer][1], 0.4),
            ("shear_eta", DynViscosity::DIMS, ETA_COLD[layer][2], 3.0),
        ] {
            let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
            if name.ends_with("_eta") {
                validity = validity.with("omega", BAND.0, BAND.1);
            }
            let constant = name == "density";
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if constant {
                    PropertyValue::Scalar { value: cold, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, cold), (400.0, cold * hot_factor)], dims,
                    }
                },
                validity,
                interpolation: if constant { InterpolationPolicy::ConstantWithinValidity }
                    else { InterpolationPolicy::LinearInside },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured fixed-basis Timoshenko law".into(),
                    license: "test-only".into(), artifact: None,
                },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(
                name, dims, ScalarAdmissibility::Finite,
            ).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: format!("manufactured layer {layer}"), phase: "solid".into(),
            process: "affine thermal mechanical fixture".into(), revision: 0,
        }, claims, vec![]).unwrap();
        let point = QueryPoint::new().with("T", 300.0).unwrap().with("omega", 1.0).unwrap();
        let material = resolve_material_state_point(
            &card, &point, &requirements, MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        layers.push(ConcentricMaterialLayer { outer_radius_m: RADII[layer], material });
        cards.push(card);
    }
    let elastic = ConcentricElasticSection::try_new(0.0, layers).unwrap();
    let law = ConcentricSectionTemperatureLaw::try_new(
        elastic.clone(), &[&cards[0], &cards[1]], &["T", "T"],
    ).unwrap();
    let cold = FrameViscousSection::try_new_with_shear(
        elastic.with_shear_modulus("shear_modulus", SHEAR_CORRECTION).unwrap()
            .with_transverse_viscosity("shear_eta").unwrap(),
        "axial_eta", "bending_eta",
    ).unwrap();
    let hot = cold.with_material_section(law.at_temperatures_k(&HOT_T).unwrap()).unwrap();
    (cold, hot)
}

fn reference_fractions(length: f64) -> (f64, f64) {
    let bending = (0..2).map(|i| YOUNG_COLD[i] * inertia(i)).sum::<f64>();
    let shear = SHEAR_CORRECTION * (0..2).map(|i| SHEAR_COLD[i] * area(i)).sum::<f64>();
    let phi = 12.0 * bending / (shear * length * length);
    (1.0 / (1.0 + phi), phi / (1.0 + phi))
}

#[test]
fn updated_layer_stiffness_and_loss_use_original_timoshenko_strains() {
    let (cold, hot) = sections();
    for length in [0.05_f64, 0.7, 3.0] {
        let basis = cold.reference_kinematics(length).unwrap();
        let (beta, shear_fraction) = reference_fractions(length);
        for x in [[0.001 * length, -0.002 * length, 0.0007 * length],
            [0.0, 0.0, 0.003 * length], [0.0; 3]] {
            let axial = x[0] / length;
            let shear = -shear_fraction * x[1] / (2.0 * length);
            // Two-point Gauss exactly integrates the affine curvature squared.
            let curvature = [-1.0_f64 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()]
                .map(|xi| -(3.0 * beta * xi * x[1] + x[2]) / length.powi(2));
            let curvature_square = 0.5 * (curvature[0].powi(2) + curvature[1].powi(2));
            for (section, temperatures) in [(&cold, [300.0; 2]), (&hot, HOT_T)] {
                let damping = basis.damping_n_s_m(section).unwrap();
                let port_work = core::array::from_fn(|i| damping[i] * x[i] * x[i]);
                let mut total_heat = 0.0;
                for layer in 0..2 {
                    let t = relative_temperature(temperatures[layer]);
                    let young = YOUNG_COLD[layer] * (1.0 - 0.4 * t);
                    let modulus = SHEAR_COLD[layer] * (1.0 + 1.8 * t);
                    let stiffness = basis.layer_stiffness_n_m(section, layer).unwrap();
                    let actual = 0.5 * stiffness.iter().zip(x).map(|(k, x)| k * x * x).sum::<f64>();
                    let expected = 0.5 * length * (young * area(layer) * axial.powi(2)
                        + young * inertia(layer) * curvature_square
                        + SHEAR_CORRECTION * modulus * area(layer) * shear.powi(2));
                    close(actual, expected);
                    let heat = basis.layer_heat_fractions(section, layer).unwrap()
                        .partition_energy_j(port_work).unwrap();
                    close(heat[0], length * ETA_COLD[layer][0] * (1.0 + t) * area(layer) * axial.powi(2));
                    close(heat[1], length * ETA_COLD[layer][1] * (1.0 - 0.6 * t) * inertia(layer) * curvature_square);
                    close(heat[2], length * SHEAR_CORRECTION * ETA_COLD[layer][2] * (1.0 + 2.0 * t)
                        * area(layer) * shear.powi(2));
                    total_heat += heat.iter().sum::<f64>();
                }
                close(total_heat, port_work.iter().sum());
            }
        }
        // Re-solving the static interpolation at hot E/G is a different
        // approximation; a thermal refresh must not silently change this basis.
        let frozen = basis.layer_stiffness_n_m(&hot, 0).unwrap()[1];
        let refitted = hot.reference_kinematics(length).unwrap()
            .layer_stiffness_n_m(&hot, 0).unwrap()[1];
        assert!((frozen - refitted).abs() > 1.0e-4 * frozen.abs());
    }
}

#[test]
fn fixed_basis_parameter_work_matches_constitutive_energy_change() {
    let (cold, hot) = sections();
    let length = 0.13_f64;
    let basis = cold.reference_kinematics(length).unwrap();
    let (beta, shear_fraction) = reference_fractions(length);
    let x = [0.0001, -0.0003, 0.0002];
    let axial_square = (x[0] / length).powi(2);
    let curvature_square = (3.0 * (beta * x[1]).powi(2) + x[2] * x[2]) / length.powi(4);
    let shear_square = (shear_fraction * x[1] / (2.0 * length)).powi(2);
    for layer in 0..2 {
        let old = basis.layer_stiffness_n_m(&cold, layer).unwrap();
        let new = basis.layer_stiffness_n_m(&hot, layer).unwrap();
        let port_change = 0.5 * (0..3).map(|j| (new[j] - old[j]) * x[j] * x[j]).sum::<f64>();
        let t = relative_temperature(HOT_T[layer]);
        let constitutive_change = 0.5 * length * (
            -0.4 * t * YOUNG_COLD[layer] * (area(layer) * axial_square + inertia(layer) * curvature_square)
            + 1.8 * t * SHEAR_COLD[layer] * SHEAR_CORRECTION * area(layer) * shear_square
        );
        close(port_change, constitutive_change);
    }
}

#[test]
fn energy_averaged_curvature_does_not_replace_the_peak_strain_bound() {
    let (section, _) = sections();
    let length = 0.2_f64;
    let basis = section.reference_kinematics(length).unwrap();
    let (beta, shear_fraction) = reference_fractions(length);
    for x in [[0.0002, 0.0007, 0.0003], [0.0, -0.0007, 0.0003], [0.0, 0.0, 0.0003]] {
        let measured = basis.strain_measures(x).unwrap();
        close(measured.axial_strain, x[0] / length);
        close(measured.transverse_shear_strain, -shear_fraction * x[1] / (2.0 * length));
        close(measured.rms_curvature_per_m,
            (3.0 * (beta * x[1]).powi(2) + x[2] * x[2]).sqrt() / length.powi(2));
        close(measured.maximum_absolute_curvature_per_m,
            (3.0 * beta * x[1].abs() + x[2].abs()) / length.powi(2));
        assert!(measured.maximum_absolute_curvature_per_m >= measured.rms_curvature_per_m);
    }
    assert!(basis.strain_measures([f64::NAN, 0.0, 0.0]).is_err());
}
