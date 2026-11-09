use super::finite_layers::debye;
use super::*;
use fs_em::plane_wave::layers::{DielectricLayer, IncoherentLayerCascade, PlanarDielectricLayers};
use fs_math::det;

const OMEGA: f64 = 1e9;

fn k() -> f64 {
    OMEGA * det::sqrt(E) * det::sqrt(M)
}

/// G1: independent geometric sums give a thick transparent plate's two-surface
/// reflection. Its interior angle must retain Snell refraction for TE and TM.
#[test]
fn bare_glass_retains_all_incoherent_internal_reflections() {
    let air = material(E, M);
    let glass = material(2.25 * E, M);
    let front = PlanarDielectricLayers::try_new(
        &air, DielectricResponse::Lossless(&glass), OMEGA, &[],
    ).unwrap();
    let back = PlanarDielectricLayers::try_new(
        &glass, DielectricResponse::Lossless(&air), OMEGA, &[],
    ).unwrap();
    for cosine in [1.0_f64, 0.7, 0.2] {
        let glass_cosine = (1.0 - (1.0 - cosine * cosine) / 2.25).sqrt();
        for polarization in [Polarization::Te, Polarization::Tm] {
            let first = front.prepare_illumination(cosine, polarization).unwrap();
            let second = back.prepare_illumination(glass_cosine, polarization).unwrap();
            let cascade = IncoherentLayerCascade::try_new(&[&first, &second]).unwrap();
            let (yi, yg) = match polarization {
                Polarization::Te => (cosine, 1.5 * glass_cosine),
                Polarization::Tm => (cosine.recip(), 1.5 / glass_cosine),
            };
            let r = ((yi - yg) / (yi + yg)).powi(2);
            let reflected = 2.0 * r / (1.0 + r);
            let transmitted = (1.0 - r) / (1.0 + r);
            let matrix = cascade.power_scattering();
            for side in 0..2 {
                near(matrix[side][side], reflected, 8e-14);
                near(matrix[1 - side][side], transmitted, 8e-14);
            }
            let output = cascade.evaluate([2.0, 0.4], &mut []).unwrap();
            near(output.outgoing_w_m2[0], 2.0 * reflected + 0.4 * transmitted, 2e-13);
            near(output.outgoing_w_m2[1], 2.0 * transmitted + 0.4 * reflected, 2e-13);
            assert_eq!(output.absorbed_w_m2, 0.0);
            assert_eq!(cascade.layer_range(0), Some(0..0));
            assert_eq!(cascade.layer_range(1), Some(0..0));
            assert_eq!(cascade.layer_range(2), None);
            assert_eq!(cascade.angular_frequency_rad_s(), OMEGA);
            assert_eq!(cascade.polarization(), polarization);
            assert!(std::ptr::eq(cascade.blocks()[0], &first));
        }
    }
}

/// G1/G3: exact coherent coating fields are retained, but the thick transparent
/// substrate's phase is averaged. Independently sweep one full round-trip phase
/// in the original all-coherent stack and compare both ports and each coating.
#[test]
fn absorbing_coatings_match_coherent_substrate_phase_average_and_native_heat() {
    let left = material(E, M);
    let substrate = material(4.0 * E, M);
    let right = material(2.25 * E, M);
    let film_a = debye(3.0 * E, 4.0 * E, 0.03, 0.001);
    let film_b = debye(1.5 * E, 0.4 * E, 0.03, 0.001);
    let coating_a = DielectricLayer { response: DielectricResponse::Debye(&film_a), thickness_m: 0.13 / k() };
    let coating_b = DielectricLayer { response: DielectricResponse::Debye(&film_b), thickness_m: 0.21 / k() };
    let front = PlanarDielectricLayers::try_new(
        &left, DielectricResponse::Lossless(&substrate), OMEGA, &[coating_a],
    ).unwrap();
    let back = PlanarDielectricLayers::try_new(
        &substrate, DielectricResponse::Lossless(&right), OMEGA, &[coating_b],
    ).unwrap();
    let cosine = 0.73_f64;
    let q = (4.0 - (1.0 - cosine * cosine)).sqrt();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let first = front.prepare_illumination(cosine, polarization).unwrap();
        let second = back.prepare_illumination(q / 2.0, polarization).unwrap();
        let cascade = IncoherentLayerCascade::try_new(&[&first, &second]).unwrap();
        let mut heat = [0.0; 2];
        let power = cascade.evaluate([1.7, 0.6], &mut heat).unwrap();
        assert!(heat.iter().all(|x| *x > 0.0));
        assert_eq!(cascade.layer_range(0), Some(0..1));
        assert_eq!(cascade.layer_range(1), Some(1..2));
        near(power.absorbed_w_m2, heat[0] + heat[1], 2e-15);
        near(cascade.power_scattering()[0][1], cascade.power_scattering()[1][0], 2e-14);

        let mut expected = [0.0; 4];
        let samples = 512;
        for sample in 0..samples {
            let phase = (32.0 + (sample as f64 + 0.5) / samples as f64) * core::f64::consts::PI;
            let all = PlanarDielectricLayers::try_new(
                &left, DielectricResponse::Lossless(&right), OMEGA,
                &[
                    coating_a,
                    DielectricLayer {
                        response: DielectricResponse::Lossless(&substrate),
                        thickness_m: phase / (q * k()),
                    },
                    coating_b,
                ],
            ).unwrap();
            let coherent = all.prepare_illumination(cosine, polarization).unwrap();
            let y = coherent.port_admittance_s();
            let mut coherent_heat = [0.0; 3];
            for side in 0..2 {
                let mut source = [C64::ZERO; 2];
                source[side] = C64::from_re((2.0 * [1.7, 0.6][side] / y[side]).sqrt());
                let result = coherent.evaluate(source, &mut coherent_heat).unwrap();
                expected[0] += result.power.outgoing_w_m2[0] / samples as f64;
                expected[1] += result.power.outgoing_w_m2[1] / samples as f64;
                expected[2] += coherent_heat[0] / samples as f64;
                expected[3] += coherent_heat[2] / samples as f64;
                assert_eq!(coherent_heat[1], 0.0);
            }
        }
        for (actual, reference) in [power.outgoing_w_m2[0], power.outgoing_w_m2[1], heat[0], heat[1]]
            .into_iter().zip(expected)
        {
            near(actual, reference, 2e-12);
        }

        // A single native block must retain exactly its independent face
        // illumination, rather than dropping interference within the coating.
        let one = IncoherentLayerCascade::try_new(&[&first]).unwrap();
        let y = first.port_admittance_s();
        let input = [C64::from_re((3.4 / y[0]).sqrt()), C64::from_re((1.2 / y[1]).sqrt())];
        let mut native_heat = [0.0];
        let native = first.evaluate_independent_sources(input, &mut native_heat).unwrap();
        let mut one_heat = [0.0];
        let result = one.evaluate([1.7, 0.6], &mut one_heat).unwrap();
        near(result.outgoing_w_m2[0], native.outgoing_w_m2[0], 3e-14);
        near(result.outgoing_w_m2[1], native.outgoing_w_m2[1], 3e-14);
        near(one_heat[0], native_heat[0], 3e-14);
        let mut unchanged = [11.0; 2];
        assert!(cascade.evaluate([-1.0, 0.0], &mut unchanged).is_err());
        assert_eq!(unchanged, [11.0; 2]);
        assert!(cascade.evaluate([f64::NAN, 0.0], &mut unchanged).is_err());
        assert_eq!(unchanged, [11.0; 2]);
        assert!(cascade.evaluate([1.0, 0.0], &mut []).is_err());
        assert_eq!(cascade.evaluate([0.0; 2], &mut heat).unwrap().absorbed_w_m2, 0.0);
        assert_eq!(heat, [0.0; 2]);
    }
}

/// G1: independent lossless intensity interfaces add their R/T transport
/// resistances. Equal media inside repeated native sources remain connected.
#[test]
fn multiple_incoherent_connections_obey_lossless_transport_resistance() {
    let air = material(E, M);
    let glass = material(4.0 * E, M);
    let into = PlanarDielectricLayers::try_new(&air, DielectricResponse::Lossless(&glass), OMEGA, &[]).unwrap();
    let out = PlanarDielectricLayers::try_new(&glass, DielectricResponse::Lossless(&air), OMEGA, &[]).unwrap();
    let a = into.prepare_illumination(1.0, Polarization::Te).unwrap();
    let b = out.prepare_illumination(1.0, Polarization::Te).unwrap();
    let cascade = IncoherentLayerCascade::try_new(&[&a, &b, &a, &b]).unwrap();
    let expected_transmission = 1.0 / (1.0 + 4.0 * (1.0 / 9.0) / (8.0 / 9.0));
    let matrix = cascade.power_scattering();
    for side in 0..2 {
        near(matrix[side][side], 1.0 - expected_transmission, 2e-14);
        near(matrix[1 - side][side], expected_transmission, 2e-14);
    }
}

/// G1: both native reflectances round to one, yet their positive transmissions
/// resolve the finite geometric sum. The denominator cannot be 1-R1*R2.
#[test]
fn nearly_perfect_interfaces_preserve_weak_escape_without_false_cavity_poles() {
    let air = material(E, M);
    let dense = material(1e40 * E, M);
    let into = PlanarDielectricLayers::try_new(&air, DielectricResponse::Lossless(&dense), OMEGA, &[]).unwrap();
    let out = PlanarDielectricLayers::try_new(&dense, DielectricResponse::Lossless(&air), OMEGA, &[]).unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let a = into.prepare_illumination(1.0, polarization).unwrap();
        let b = out.prepare_illumination(1.0, polarization).unwrap();
        assert_eq!(a.tangential_e_scattering()[1][1].norm_sq(), 1.0);
        assert_eq!(b.tangential_e_scattering()[0][0].norm_sq(), 1.0);
        let cascade = IncoherentLayerCascade::try_new(&[&a, &b]).unwrap();
        let matrix = cascade.power_scattering();
        let expected = 2.0 / (1e20 + 1e-20);
        for transmission in [matrix[0][1], matrix[1][0]] {
            assert!(transmission > 0.0);
            near(transmission, expected, 2e-33);
        }
        near(cascade.evaluate([1.0, 0.0], &mut []).unwrap().outgoing_w_m2[1], expected, 2e-33);
    }
}

/// G0/G3: opaque internal barriers can be coherent blocks, but evanescent
/// connecting modes, wrong sources, changed directions and mixed frequencies
/// cannot be silently interpreted as ordinary propagating intensity ports.
#[test]
fn opaque_limits_and_original_connection_admission_are_explicit() {
    let high = material(4.0 * E, M);
    let low = material(E, M);
    let barrier = PlanarDielectricLayers::try_new(
        &high, DielectricResponse::Lossless(&high), OMEGA,
        &[DielectricLayer { response: DielectricResponse::Lossless(&low), thickness_m: 1000.0 / k() }],
    ).unwrap();
    let native = barrier.prepare_illumination(0.5, Polarization::Te).unwrap();
    let cascade = IncoherentLayerCascade::try_new(&[&native, &native]).unwrap();
    let mut heat = [0.0; 2];
    let response = cascade.evaluate([0.7, 0.2], &mut heat).unwrap();
    near(response.outgoing_w_m2[0], 0.7, 3e-14);
    near(response.outgoing_w_m2[1], 0.2, 3e-14);
    assert_eq!(heat, [0.0; 2]);
    assert_eq!(cascade.power_scattering()[0][1], 0.0);
    assert_eq!(cascade.power_scattering()[1][0], 0.0);
    let evanescent_exit = PlanarDielectricLayers::try_new(&high, DielectricResponse::Lossless(&low), OMEGA, &[]).unwrap();
    assert!(evanescent_exit.prepare_illumination(0.5, Polarization::Te).is_err());

    let identity = PlanarDielectricLayers::try_new(&high, DielectricResponse::Lossless(&high), OMEGA, &[]).unwrap();
    let good = identity.prepare_illumination(0.5, Polarization::Te).unwrap();
    let angle = identity.prepare_illumination(0.6, Polarization::Te).unwrap();
    let polarization = identity.prepare_illumination(0.5, Polarization::Tm).unwrap();
    for wrong in [&angle, &polarization] {
        assert!(IncoherentLayerCascade::try_new(&[&good, wrong]).is_err());
    }
    let equal_numbers = material(4.0 * E, M);
    let wrong_source = PlanarDielectricLayers::try_new(&equal_numbers, DielectricResponse::Lossless(&high), OMEGA, &[]).unwrap();
    let wrong_source = wrong_source.prepare_illumination(0.5, Polarization::Te).unwrap();
    assert!(IncoherentLayerCascade::try_new(&[&good, &wrong_source]).is_err());
    let other_frequency = PlanarDielectricLayers::try_new(&high, DielectricResponse::Lossless(&high), 1.1 * OMEGA, &[]).unwrap();
    let other_frequency = other_frequency.prepare_illumination(0.5, Polarization::Te).unwrap();
    assert!(IncoherentLayerCascade::try_new(&[&good, &other_frequency]).is_err());
    assert!(IncoherentLayerCascade::try_new(&[]).is_err());
    assert!(IncoherentLayerCascade::try_new(&[&good; 66]).is_err());
    assert!(IncoherentLayerCascade::try_new(&[&native; 65]).is_err());
    let repeated = IncoherentLayerCascade::try_new(&[&good; 65]).unwrap();
    let output = repeated.evaluate([1.0, 0.0], &mut []).unwrap();
    near(output.outgoing_w_m2[1], 1.0, 2e-14);
    assert_eq!(output.outgoing_w_m2[0], 0.0);
}
