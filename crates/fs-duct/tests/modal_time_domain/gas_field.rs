//! G1/G3: original primitive sources in the complete modal weak equations.
use super::*;
use fs_duct::modal::time_domain::{ModalMolecularBudget, assemble_lossless_field};
use fs_matdb::QueryPoint;
use fs_material::gas_field::{AxialGasField, PrimitiveGasState};
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium, SpeciesRelaxation};
use fs_material::thermochemical_gas::packs::from_packs;
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

// These existing authored species/model/transport packs retain their immutable
// receipts. They are synthetic analytic data, not measured-gas validation.
#[path = "../../../fs-couple/src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const PRESSURE: f64 = 101_325.0;

fn primitive(temperature_k: f64, fraction_a: f64) -> PrimitiveGasState {
    PrimitiveGasState {
        temperature_k,
        composition: Composition::new(
            CompositionBasis::MoleFraction,
            vec![fraction_a, 1.0 - fraction_a],
        )
        .unwrap(),
    }
}

// The fixture's two species, declared non-relaxing (zero reservoirs). A present
// species without relaxation data is refused, so the empty molecular law names
// them; it gives a molecular state with no relaxation processes.
fn no_relaxation() -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        ["A", "B"]
            .into_iter()
            .map(|name| SpeciesRelaxation {
                species: fs_thermochem::SpeciesId::new(name).unwrap(),
                processes: Vec::new(),
            })
            .collect(),
    )
    .unwrap()
}

fn field(
    stations: Vec<f64>,
    temperatures: &[f64],
    fractions: &[f64],
    cellwise: bool,
) -> AxialGasField {
    let a = gas_support::packs("A", 0.028, 3.5, 200.0);
    let b = gas_support::packs("B", 0.044, 4.5, 200.0);
    let model = from_packs(
        &[a.binding(), b.binding()],
        primitive(300.0, 0.5).composition,
        None,
    )
    .unwrap();
    let states = temperatures
        .iter()
        .zip(fractions)
        .map(|(&t, &f)| primitive(t, f))
        .collect();
    let constructor = if cellwise {
        AxialGasField::new_piecewise_constant
    } else {
        AxialGasField::new
    };
    constructor(
        stations,
        PRESSURE,
        states,
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

fn same_operators(
    a: &fs_duct::modal::time_domain::ModalGalerkin,
    b: &fs_duct::modal::time_domain::ModalGalerkin,
) {
    assert_eq!(a.source().mass(), b.source().mass());
    assert_eq!(a.source().stiffness(), b.source().stiffness());
    assert_eq!(a.source().inverse_mass(), b.source().inverse_mass());
    assert_eq!(a.system().structure(), b.system().structure());
    assert_eq!(
        a.source().inlet_reference_impedance(),
        b.source().inlet_reference_impedance()
    );
    assert_eq!(
        a.source().mouth_reference_impedance(),
        b.source().mouth_reference_impedance()
    );
}

#[test]
fn primitive_uniform_and_cellwise_fields_retain_exact_original_operators_and_sources() {
    let length = 0.08;
    let profile = RadiusProfile::exponential(0.03, 3.0, length).unwrap();
    let uniform = field(vec![0.0, length], &[300.0, 300.0], &[0.8, 0.8], false);
    let regional = field(
        vec![0.0, length / 2.0, length],
        &[300.0, 300.0, 550.0],
        &[0.8, 0.8, 0.3],
        true,
    );
    let cold = *uniform.state_at(0.0).unwrap().state();
    let warm = *regional.state_at(length).unwrap().state();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease, ModalEnd::Ports] {
        let expected =
            assemble_lossless(&profile, &cold, ModeCount(2), 8, end, 124, || false).unwrap();
        let actual =
            assemble_lossless_field(&profile, &uniform, ModeCount(2), 8, end, 124, || false)
                .unwrap();
        same_operators(&actual, &expected);
        assert_eq!(actual.source().gas_field(), Some(&uniform));
        assert_eq!(actual.source().gas_at_element(7), Some(&cold));
        assert_eq!(
            actual.source().gas_at_quadrature(7, 1),
            uniform.uniform_state()
        );

        let gases = two_regions(cold, warm, 8);
        let expected =
            assemble_lossless_regional(&profile, &gases, ModeCount(2), end, 124, || false).unwrap();
        let actual =
            assemble_lossless_field(&profile, &regional, ModeCount(2), 8, end, 124, || false)
                .unwrap();
        same_operators(&actual, &expected);
        assert_eq!(actual.source().gas_field(), Some(&regional));
        assert_eq!(actual.source().regional_gases(), None);
        assert_eq!(actual.source().gas_at_element(0), None);
        assert_eq!(actual.source().outlet_gas(), &warm);
        for element in 0..8 {
            for point in 0..2 {
                let x = actual
                    .source()
                    .gas_quadrature_position_m(element, point)
                    .unwrap();
                let sample = actual.source().gas_at_quadrature(element, point).unwrap();
                assert_eq!(sample, &regional.state_at(x).unwrap());
                assert_eq!(sample.source_packs(), regional.model().source_packs());
            }
        }
        assert!(actual.source().gas_at_quadrature(8, 0).is_none());
        assert!(actual.source().gas_at_quadrature(0, 2).is_none());
        assert!(actual.source().gas_quadrature_position_m(8, 0).is_none());
    }
}

// Independent ideal-mixture expressions for the authored two-species fixture.
// Each NASA-9 cp/R is its supplied constant plus 0.001*T.
fn analytic_gas(t: f64, fraction_a: f64) -> (f64, f64) {
    let mass = 0.028 * fraction_a + 0.044 * (1.0 - fraction_a);
    let cp_r = 4.5 - fraction_a + 0.001 * t;
    let gamma = cp_r / (cp_r - 1.0);
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / mass;
    (PRESSURE / (r * t), gamma * r * t)
}

#[test]
fn continuous_primitive_field_integrates_physical_compression_and_transverse_stiffness() {
    let (length, radius) = (0.08, 0.03);
    let (ta, tb, fa, fb) = (260.0, 800.0, 0.99, 0.01);
    let original = field(vec![0.0, length], &[ta, tb], &[fa, fb], false);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    let area = core::f64::consts::PI * radius * radius;
    let (ma, mb) = (0.044 - 0.016 * fa, 0.044 - 0.016 * fb);
    let (ca, cb) = (4.5 - fa + 0.001 * ta, 4.5 - fb + 0.001 * tb);
    let quotient = (tb - ta) / (mb - ma);
    let inverse_density_integral = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * length / PRESSURE
        * (quotient + (ta - quotient * ma) * det::ln(mb / ma) / (mb - ma));
    let compressibility_integral = length / PRESSURE * (1.0 - det::ln(cb / ca) / (cb - ca));
    let expected_mass = area * compressibility_integral;
    let expected_gradient = area * inverse_density_integral;
    let mut errors = Vec::new();
    for elements in [2, 4, 8] {
        let image = assemble_lossless_field(
            &profile,
            &original,
            ModeCount(2),
            elements,
            ModalEnd::Ports,
            124,
            || false,
        )
        .unwrap();
        let source = image.source();
        let n = source.coordinates();
        let mut mass = 0.0;
        let mut gradient = 0.0;
        let mut transverse = 0.0;
        for row in 0..=elements {
            for column in 0..=elements {
                mass += source.mass()[(2 * row) * n + 2 * column];
                gradient += source.nodes_m()[row]
                    * source.stiffness()[(2 * row) * n + 2 * column]
                    * source.nodes_m()[column];
                transverse += source.stiffness()[(2 * row + 1) * n + 2 * column + 1];
            }
            if row < elements {
                for point in 0..2 {
                    let x = source.gas_quadrature_position_m(row, point).unwrap();
                    let sample = source.gas_at_quadrature(row, point).unwrap().state();
                    let (rho, c2) =
                        analytic_gas(ta + (tb - ta) * x / length, fa + (fb - fa) * x / length);
                    assert!((sample.density / rho - 1.0).abs() < 2e-14);
                    assert!((sample.sound_speed.powi(2) / c2 - 1.0).abs() < 2e-14);
                }
            }
        }
        let expected_transverse =
            (source.radial_eigenvalues()[1] / radius).powi(2) * expected_gradient;
        let defect = [
            (mass / expected_mass - 1.0).abs(),
            (gradient / expected_gradient - 1.0).abs(),
            (transverse / expected_transverse - 1.0).abs(),
        ];
        errors.push(defect);
        assert!((gradient / expected_gradient - transverse / expected_transverse).abs() < 5e-13);
    }
    for axis in 0..3 {
        assert!(
            errors[0][axis] > 1e-7,
            "live Gauss defect required: {errors:?}"
        );
        assert!(
            errors[1][axis] < 0.1 * errors[0][axis] && errors[2][axis] < 0.1 * errors[1][axis],
            "{errors:?}"
        );
        assert!(errors[2][axis] < 1e-6, "{errors:?}");
    }
}

// Independent pressure/flow ODE, integrated from the actual far-end condition.
// Analytic coefficients above do not query the production gas or FE owners.
fn gradient_impedance(
    radius: f64,
    length: f64,
    gamma: f64,
    omega: f64,
    reference: f64,
    end: ModalEnd,
    steps: usize,
) -> C64 {
    let area = core::f64::consts::PI * radius * radius;
    let rate = |x: f64, y: [C64; 2]| {
        let (rho, c2) = analytic_gas(285.0 + 365.0 * x / length, 0.9 - 0.65 * x / length);
        [
            C64::new(0.0, omega * rho / (area * reference)) * y[1],
            C64::new(
                0.0,
                reference * area / (omega * rho) * (omega * omega / c2 - (gamma / radius).powi(2)),
            ) * y[0],
        ]
    };
    let add =
        |a: [C64; 2], b: [C64; 2], scale: f64| [a[0] + b[0].scale(scale), a[1] + b[1].scale(scale)];
    let mut y = match end {
        ModalEnd::Rigid => [C64::ONE, C64::ZERO],
        ModalEnd::PressureRelease => [C64::ZERO, C64::ONE],
        ModalEnd::Ports => unreachable!(),
    };
    let h = -length / steps as f64;
    for step in 0..steps {
        let x = length + step as f64 * h;
        let a = rate(x, y);
        let b = rate(x + 0.5 * h, add(y, a, 0.5 * h));
        let c = rate(x + 0.5 * h, add(y, b, 0.5 * h));
        let d = rate(x + h, add(y, c, h));
        for i in 0..2 {
            y[i] = y[i] + (a[i] + b[i].scale(2.0) + c[i].scale(2.0) + d[i]).scale(h / 6.0);
        }
    }
    (y[0] / y[1]).scale(reference)
}

#[test]
fn primitive_gradient_retains_local_cutoff_and_converges_on_the_complete_modal_response() {
    let (length, radius, omega) = (0.04, 0.05, TAU * 4500.0);
    let original = field(vec![0.0, length], &[285.0, 650.0], &[0.9, 0.25], false);
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    for end in [ModalEnd::Rigid, ModalEnd::PressureRelease] {
        let mut errors = Vec::new();
        for elements in [4, 8, 16] {
            let image = assemble_lossless_field(
                &profile,
                &original,
                ModeCount(2),
                elements,
                end,
                124,
                || false,
            )
            .unwrap();
            let source = image.source();
            let zref = source.inlet_reference_impedance();
            let gamma = source.radial_eigenvalues()[1];
            assert!(omega > source.gas().sound_speed * gamma / radius);
            assert!(omega < source.outlet_gas().sound_speed * gamma / radius);
            let model =
                HarmonicModel::from_ode(image.system(), harmonic_budget(), || false).unwrap();
            let factor = model.factor(omega, || false).unwrap();
            let mut actual = vec![C64::ZERO; 4];
            let mut expected = vec![C64::ZERO; 4];
            for mode in 0..2 {
                let mut drive = vec![C64::ZERO; 2];
                drive[mode] = C64::ONE;
                let response = factor.response(&drive, || false).unwrap();
                actual[mode] = response.output[0];
                actual[2 + mode] = response.output[1];
                assert_eq!(response.dissipated_power_w, 0.0);
                expected[mode * 2 + mode] = gradient_impedance(
                    radius,
                    length,
                    source.radial_eigenvalues()[mode],
                    omega,
                    zref,
                    end,
                    4096,
                );
                let coarse = gradient_impedance(
                    radius,
                    length,
                    source.radial_eigenvalues()[mode],
                    omega,
                    zref,
                    end,
                    2048,
                );
                assert!((expected[mode * 2 + mode] - coarse).abs() < 2e-10 * zref);
            }
            assert!(actual[1].abs() + actual[2].abs() < 1e-10 * zref);
            errors.push(distance(
                &reflection(&actual, zref, 2),
                &reflection(&expected, zref, 2),
            ));
        }
        assert!(errors[0] > 1e-5, "{errors:?}");
        assert!(
            errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 0.025, "{errors:?}");
    }
}

#[test]
fn primitive_field_two_face_work_is_reciprocal_and_lossless() {
    let length = 0.08;
    let original = field(vec![0.0, length], &[300.0, 600.0], &[0.8, 0.2], false);
    let profile = RadiusProfile::exponential(0.03, 4.0, length).unwrap();
    let image = assemble_lossless_field(
        &profile,
        &original,
        ModeCount(2),
        8,
        ModalEnd::Ports,
        124,
        || false,
    )
    .unwrap();
    let model = HarmonicModel::from_ode(image.system(), harmonic_budget(), || false).unwrap();
    let factor = model.factor(TAU * 2500.0, || false).unwrap();
    let mut matrix = [C64::ZERO; 16];
    for column in 0..4 {
        let mut drive = [C64::ZERO; 4];
        drive[column] = C64::ONE;
        let response = factor.response(&drive, || false).unwrap();
        for row in 0..4 {
            matrix[row * 4 + column] = response.output[row];
        }
    }
    let reference = image.source().inlet_reference_impedance();
    for i in 0..4 {
        for j in 0..4 {
            assert!((matrix[i * 4 + j] - matrix[j * 4 + i]).abs() < 2e-9 * reference);
        }
    }
    assert!(
        matrix[1].abs() > 1e-3 * reference,
        "flare must couple radial modes"
    );
    let drive = [
        C64::new(0.1, 0.2),
        C64::new(-0.3, 0.1),
        C64::new(0.2, -0.4),
        C64::new(0.1, 0.3),
    ];
    let response = factor.response(&drive, || false).unwrap();
    let powers: Vec<_> = response
        .output
        .iter()
        .zip(drive)
        .map(|(&p, u)| 0.5 * (p * u.conj()).re)
        .collect();
    assert!(powers.iter().any(|p| p.abs() > 1e-3 * reference));
    assert!(powers.iter().sum::<f64>().abs() < 1e-9 * reference);
    assert_eq!(response.dissipated_power_w, 0.0);
    assert!(response.mean_stored_energy_j > 0.0);
}

#[test]
fn primitive_field_refuses_unresolved_interfaces_capacity_cancellation_and_uniform_molecular_shortcut()
 {
    let length = 0.08;
    let profile = RadiusProfile::exponential(0.03, 0.0, length).unwrap();
    let original = field(
        vec![0.0, length / 3.0, length],
        &[300.0, 400.0, 500.0],
        &[0.8, 0.5, 0.2],
        true,
    );
    let make = |elements, max_states| {
        assemble_lossless_field(
            &profile,
            &original,
            ModeCount(2),
            elements,
            ModalEnd::Ports,
            max_states,
            || false,
        )
    };
    assert!(matches!(make(8, 124), Err(ModalTimeError::Invalid(_))));
    assert!(make(12, 124).is_ok());
    assert!(matches!(make(12, 4), Err(ModalTimeError::Invalid(_))));
    let other = RadiusProfile::exponential(0.03, 0.0, length * 0.99).unwrap();
    assert!(matches!(
        assemble_lossless_field(
            &other,
            &original,
            ModeCount(1),
            12,
            ModalEnd::Rigid,
            124,
            || false
        ),
        Err(ModalTimeError::Invalid(_))
    ));
    let mut all_polls = 0;
    assemble_lossless_field(
        &profile,
        &original,
        ModeCount(2),
        12,
        ModalEnd::Ports,
        124,
        || {
            all_polls += 1;
            false
        },
    )
    .unwrap();
    for limit in [1, 7, all_polls] {
        let mut polls = 0;
        let result = assemble_lossless_field(
            &profile,
            &original,
            ModeCount(2),
            12,
            ModalEnd::Ports,
            124,
            || {
                polls += 1;
                polls == limit
            },
        );
        assert!(matches!(
            result,
            Err(ModalTimeError::Duct(DuctError::Cancelled { .. }))
        ));
        assert_eq!(polls, limit);
    }
    assert!(make(12, 124).is_ok());
    // Even an empty molecular declaration may not relabel a spatial gas field
    // as one uniform gas; that separate physical extension remains explicit.
    let inlet = original.state_at(0.0).unwrap();
    let medium = RelaxingMedium::from_thermochemical(&inlet).unwrap();
    let molecular = no_relaxation().evaluate(&medium).unwrap();
    let result = make(12, 124).unwrap().with_molecular_relaxation(
        molecular,
        ModalMolecularBudget {
            maximum_states: 124,
            maximum_payload_bytes: 16 * 1024 * 1024,
            maximum_preparation_work: 200_000_000,
        },
        || false,
    );
    assert!(matches!(result, Err(ModalTimeError::Invalid(_))));
    let uniform = field(vec![0.0, length], &[300.0, 300.0], &[0.8, 0.8], false);
    let medium = RelaxingMedium::from_thermochemical(uniform.uniform_state().unwrap()).unwrap();
    let molecular = no_relaxation().evaluate(&medium).unwrap();
    assert!(
        assemble_lossless_field(
            &profile,
            &uniform,
            ModeCount(2),
            4,
            ModalEnd::Ports,
            124,
            || false,
        )
        .unwrap()
        .with_molecular_relaxation(
            molecular,
            ModalMolecularBudget {
                maximum_states: 124,
                maximum_payload_bytes: 16 * 1024 * 1024,
                maximum_preparation_work: 200_000_000,
            },
            || false,
        )
        .is_ok()
    );
}
