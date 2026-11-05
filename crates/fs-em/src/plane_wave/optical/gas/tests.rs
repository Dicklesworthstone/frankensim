//! G1/G3 authored molecular controls, not measured gas qualification.
use super::*;
use crate::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumExposureLayerLimits, AngularSpectrumExposureLimits,
        AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY,
    },
    gas::GasSpec,
    gas_scattering::{
        AVOGADRO_PER_MOL, MolecularRayleigh, MolecularRayleighMixture,
        composition::MolecularCauchyDispersion,
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_qty::{
    Dims,
    semantic::{Composition, CompositionBasis, SemanticType},
};

const MU: f64 = 1.256_637_062_12e-6;
const R: [f64; 2] = [0.00015, 0.0006];
const B: [f64; 2] = [2e-6, 8e-6];
const C: [f64; 2] = [1e-7, 3e-7];

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:.3e}"
    );
}
fn complex_near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:?} != {b:?}, tolerance {tolerance:.3e}"
    );
}
fn phase(angle: f64) -> C64 {
    let (s, c) = angle.sin_cos();
    C64::new(c, s)
}
fn claim(name: &str, dims: Dims, value: f64, validity: ValidityDomain) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored molecular Maxwell control".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn snapshot(
    properties: &[(&str, Dims, f64)],
    point: QueryPoint,
    validity: ValidityDomain,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(claim(name, dims, value, validity.clone()))
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored molecular optical control".into(),
            phase: "frozen homogeneous gas".into(),
            process: "independent wave reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn magnetic(mu: f64) -> ResolvedMaterialStatePoint {
    snapshot(
        &[(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            mu,
        )],
        QueryPoint::new(),
        ValidityDomain::unconstrained(),
    )
}
fn port() -> ResolvedMaterialStatePoint {
    snapshot(
        &[
            (
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                MU,
            ),
            (
                ELECTRIC_PERMITTIVITY_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
                1.0 / (MU * VACUUM_SPEED_M_S.powi(2)),
            ),
        ],
        QueryPoint::new(),
        ValidityDomain::unconstrained(),
    )
}
fn fixture() -> (RayleighGasState, MolecularRefractionModel) {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let optics: [MolecularRayleigh; 2] = core::array::from_fn(|i| {
        let identity = claim(
            "authored molecular reference refractivity",
            Dims::NONE,
            R[i],
            ValidityDomain::unconstrained(),
        )
        .content_hash();
        MolecularRayleigh::try_new(R[i], 2.5e25, 550.0, 0.03, identity).unwrap()
    });
    let mixture = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optics,
        None,
    )
    .unwrap();
    let state = mixture.evaluate(300.0, 100_000.0).unwrap();
    let model = MolecularRefractionModel::from_mixture(&mixture, state.gas()).unwrap();
    let laws: [MolecularCauchyDispersion; 2] = core::array::from_fn(|i| {
        let identity = claim(
            "authored molecular Cauchy source",
            Dims::NONE,
            B[i],
            ValidityDomain::unconstrained(),
        )
        .content_hash();
        MolecularCauchyDispersion::try_new_in_band(optics[i], B[i], C[i], [300.0, 2500.0], identity)
            .unwrap()
    });
    (state, model.with_dispersion(&laws).unwrap())
}

// Independent molecular equations; no production index/dispersion helper.
fn molecular_reference(numbers: [f64; 2], lambda_nm: f64, dispersive: bool) -> (f64, f64) {
    let x = (lambda_nm * 1e-3).powi(-2);
    let xr = 0.55_f64.powi(-2);
    let mut l = 0.0;
    let mut dl = 0.0;
    for i in 0..2 {
        let (b, c) = if dispersive { (B[i], C[i]) } else { (0.0, 0.0) };
        let n = 1.0 + R[i] + b * (x - xr) + c * (x * x - xr * xr);
        let dn = (-2.0 * b * x - 4.0 * c * x * x) / lambda_nm;
        l += numbers[i] / 2.5e25 * (n * n - 1.0) / (n * n + 2.0);
        dl += numbers[i] / 2.5e25 * 6.0 * n * dn / (n * n + 2.0).powi(2);
    }
    let n = ((1.0 + 2.0 * l) / (1.0 - l)).sqrt();
    (n, n - lambda_nm * 1.5 * dl / (n * (1.0 - l).powi(2)))
}
fn admittance(p: Polarization, n: f64, q: f64, mu: f64) -> f64 {
    match p {
        Polarization::Te => q / (VACUUM_SPEED_M_S * mu),
        Polarization::Tm => n * n / (q * VACUUM_SPEED_M_S * mu),
    }
}
fn beam_limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 32,
        maximum_workspace_bytes: 2_000_000,
        maximum_fft_work: 100_000,
    }
}

#[test]
fn g1_original_gas_phase_group_and_explicit_mu_reach_maxwell() {
    let (gas, model) = fixture();
    let total = gas.gas().density / gas.molar_mass_kg_mol() * AVOGADRO_PER_MOL;
    let expected = molecular_reference([0.75 * total, 0.25 * total], 550.0, false).0;
    let numbers = [1.7e25, 0.6e25];
    let masses: [f64; 2] = core::array::from_fn(|i| {
        numbers[i] / AVOGADRO_PER_MOL * model.components()[i].molar_mass_kg_mol()
    });
    for mu in [MU, 2.3 * MU] {
        let magnetic = magnetic(mu);
        for lambda in [380.0, 550.0, 780.0] {
            let omega = FREQUENCY_WAVELENGTH_NM / lambda;
            let response =
                MolecularGasOpticalResponse::try_lossless_fixed_composition(&gas, &magnetic, omega)
                    .unwrap();
            near(response.phase_index(), expected, 5e-16);
            assert_eq!(response.group_index(), response.phase_index());
            assert!(
                matches!(response.source(), MolecularGasOpticalSource::FixedComposition(s) if core::ptr::eq(s, &gas))
            );
            assert!(core::ptr::eq(response.material(), &magnetic));
            assert_eq!(response.permeability_h_m(), mu);
            assert_eq!(response.complex_index(omega).unwrap().im, 0.0);
            near(
                response.permittivity_f_m(omega).unwrap().re * mu * VACUUM_SPEED_M_S.powi(2),
                expected * expected,
                1e-15,
            );
        }
        for lambda in [300.0, 401.25, 1550.0, 2500.0] {
            let omega = FREQUENCY_WAVELENGTH_NM / lambda;
            let (n, ng) = molecular_reference(numbers, lambda, true);
            for (densities, basis) in [
                (&numbers[..], PartialDensityBasis::MoleculesPerM3),
                (&masses[..], PartialDensityBasis::MassKgPerM3),
            ] {
                let response = MolecularGasOpticalResponse::try_lossless_partial_densities(
                    &model, densities, basis, &magnetic, omega,
                )
                .unwrap();
                near(response.phase_index(), n, 7e-16);
                near(response.group_index(), ng, 7e-16);
                assert!(ng > n);
                let epsilon = response.permittivity_f_m(omega).unwrap();
                assert_eq!(epsilon.im, 0.0);
                near(epsilon.re * mu * VACUUM_SPEED_M_S.powi(2), n * n, 2e-15);
                assert!((epsilon.re * mu * VACUUM_SPEED_M_S.powi(2) - ng * ng).abs() > 1e-6);
                for other in [omega.next_up(), omega.next_down(), 0.0, f64::NAN] {
                    assert!(response.permittivity_f_m(other).is_err());
                }
            }
        }
    }
}

#[test]
fn g3_original_domains_magnetic_conditions_and_source_identity_are_retained() {
    let (gas, model) = fixture();
    let magnetic = magnetic(MU);
    let numbers = [1.7e25, 0.6e25];
    let omega = FREQUENCY_WAVELENGTH_NM / 550.0;
    for (kind, value) in [
        (QuantityKind::AbsoluteTemperature, gas.gas().temperature),
        (QuantityKind::Pressure, gas.gas().pressure),
        (QuantityKind::MassConcentration, gas.gas().density),
    ] {
        let quantity = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        for factor in [1.0, 1.1] {
            let selected = snapshot(
                &[(
                    MAGNETIC_PERMEABILITY_PROPERTY,
                    MAGNETIC_PERMEABILITY_DIMS,
                    MU,
                )],
                QueryPoint::new()
                    .with_quantity("condition", quantity, factor * value)
                    .unwrap(),
                ValidityDomain::unconstrained().with_quantity(
                    "condition",
                    quantity,
                    0.5 * value,
                    2.0 * value,
                ),
            );
            assert_eq!(
                MolecularGasOpticalResponse::try_lossless_fixed_composition(&gas, &selected, omega)
                    .is_ok(),
                factor == 1.0
            );
            assert!(
                MolecularGasOpticalResponse::try_lossless_partial_densities(
                    &model,
                    &numbers,
                    PartialDensityBasis::MoleculesPerM3,
                    &selected,
                    omega
                )
                .is_err()
            );
        }
    }
    let unknown = snapshot(
        &[(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            MU,
        )],
        QueryPoint::new().with("unmapped", 1.0).unwrap(),
        ValidityDomain::unconstrained(),
    );
    assert!(
        MolecularGasOpticalResponse::try_lossless_fixed_composition(&gas, &unknown, omega).is_err()
    );
    let missing = snapshot(
        &[(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            1e-11,
        )],
        QueryPoint::new(),
        ValidityDomain::unconstrained(),
    );
    assert!(
        MolecularGasOpticalResponse::try_lossless_fixed_composition(&gas, &missing, omega).is_err()
    );
    for bad in [
        0.0,
        f64::NAN,
        (FREQUENCY_WAVELENGTH_NM / 780.0).next_down(),
        (FREQUENCY_WAVELENGTH_NM / 380.0).next_up(),
    ] {
        assert!(
            MolecularGasOpticalResponse::try_lossless_fixed_composition(&gas, &magnetic, bad)
                .is_err()
        );
    }
    for bad in [
        (FREQUENCY_WAVELENGTH_NM / 2500.0).next_down(),
        (FREQUENCY_WAVELENGTH_NM / 300.0).next_up(),
    ] {
        assert!(
            MolecularGasOpticalResponse::try_lossless_partial_densities(
                &model,
                &numbers,
                PartialDensityBasis::MoleculesPerM3,
                &magnetic,
                bad
            )
            .is_err()
        );
    }
    for bad in [
        &[][..],
        &[1.0][..],
        &[1.0, 2.0, 3.0][..],
        &[-1.0, 0.0][..],
        &[f64::NAN, 0.0][..],
        &[1e30, 0.0][..],
    ] {
        assert!(
            MolecularGasOpticalResponse::try_lossless_partial_densities(
                &model,
                bad,
                PartialDensityBasis::MoleculesPerM3,
                &magnetic,
                omega
            )
            .is_err()
        );
    }
    let vacuum = MolecularGasOpticalResponse::try_lossless_partial_densities(
        &model,
        &[0.0; 2],
        PartialDensityBasis::MoleculesPerM3,
        &magnetic,
        omega,
    )
    .unwrap();
    assert_eq!((vacuum.phase_index(), vacuum.group_index()), (1.0, 1.0));
    let a = MolecularGasOpticalResponse::try_lossless_partial_densities(
        &model,
        &numbers,
        PartialDensityBasis::MoleculesPerM3,
        &magnetic,
        omega,
    )
    .unwrap();
    let b = MolecularGasOpticalResponse::try_lossless_partial_densities(
        &model,
        &numbers,
        PartialDensityBasis::MoleculesPerM3,
        &magnetic,
        FREQUENCY_WAVELENGTH_NM / 1550.0,
    )
    .unwrap();
    assert!(a.same_source(b));
    let changed_densities = [numbers[0] * 1.1, numbers[1]];
    let changed = MolecularGasOpticalResponse::try_lossless_partial_densities(
        &model,
        &changed_densities,
        PartialDensityBasis::MoleculesPerM3,
        &magnetic,
        omega,
    )
    .unwrap();
    assert!(!a.same_source(changed));
    let foreign_model = model.clone();
    let foreign = MolecularGasOpticalResponse::try_lossless_partial_densities(
        &foreign_model,
        &numbers,
        PartialDensityBasis::MoleculesPerM3,
        &magnetic,
        omega,
    )
    .unwrap();
    assert!(!a.same_source(foreign));
}

#[test]
fn g1_molecular_interfaces_and_films_match_independent_maxwell_fields() {
    let (_, model) = fixture();
    let port = port();
    let magnetic = magnetic(1.6 * MU);
    let numbers = [1.7e25, 0.6e25];
    for lambda in [401.25, 1550.0, 2500.0] {
        let omega = FREQUENCY_WAVELENGTH_NM / lambda;
        let law = MolecularGasOpticalResponse::try_lossless_partial_densities(
            &model,
            &numbers,
            PartialDensityBasis::MoleculesPerM3,
            &magnetic,
            omega,
        )
        .unwrap();
        let n = molecular_reference(numbers, lambda, true).0;
        let interface =
            PlanarDielectricInterface::try_new(&port, DielectricResponse::MolecularGasOptical(law))
                .unwrap();
        let reverse = PlanarDielectricInterface::try_new_with_incident_response(
            DielectricResponse::MolecularGasOptical(law),
            DielectricResponse::Lossless(&port),
        )
        .unwrap();
        for ci in [1.0_f64, 0.8, 0.35] {
            let q = (n * n - (1.0 - ci * ci)).sqrt();
            let kz = omega / VACUUM_SPEED_M_S * q;
            for polarization in [Polarization::Te, Polarization::Tm] {
                let yi = admittance(polarization, 1.0, ci, MU);
                let yg = admittance(polarization, n, q, 1.6 * MU);
                let r = (yi - yg) / (yi + yg);
                let actual = interface.response(omega, ci, polarization).unwrap();
                complex_near(actual.reflected_tangential_e, C64::from_re(r), 4e-15);
                complex_near(
                    actual.transmitted_tangential_e,
                    C64::from_re(1.0 + r),
                    4e-15,
                );
                near(
                    actual.transmitted_normal_wavenumber_rad_m.re / (omega / VACUUM_SPEED_M_S),
                    q,
                    2e-15,
                );
                near(
                    actual.reflected_power_fraction + actual.transmitted_power_fraction,
                    1.0,
                    3e-15,
                );
                complex_near(
                    reverse
                        .response(omega, q / n, polarization)
                        .unwrap()
                        .reflected_tangential_e,
                    C64::from_re(-r),
                    5e-15,
                );
                for delta in [
                    0.0,
                    core::f64::consts::FRAC_PI_2,
                    core::f64::consts::PI,
                    7.3,
                ] {
                    let thickness = delta / kz;
                    let stack = PlanarDielectricLayers::try_new(
                        &port,
                        DielectricResponse::Lossless(&port),
                        omega,
                        &[DielectricLayer {
                            response: DielectricResponse::MolecularGasOptical(law),
                            thickness_m: thickness,
                        }],
                    )
                    .unwrap();
                    let p = phase(-delta);
                    let denominator = C64::ONE - (p * p).scale(r * r);
                    let expected_r = (C64::ONE - p * p).scale(r) / denominator;
                    let expected_t = p.scale(1.0 - r * r) / denominator;
                    let actual = stack.response(ci, polarization).unwrap();
                    complex_near(actual.reflected_tangential_e, expected_r, 6e-15);
                    complex_near(actual.transmitted_tangential_e, expected_t, 1e-14);
                    assert_eq!(actual.absorbed_power_fraction, 0.0);
                    let prepared = stack.prepare_layer_field(0, ci, polarization).unwrap();
                    let a = C64::from_re(1.0 + r) / denominator;
                    for z in [0.0, 0.5 * thickness, thickness] {
                        let forward = a * phase(-kz * z);
                        let backward = (a * p * p * phase(kz * z)).scale(-r);
                        let actual = prepared.field(z).unwrap();
                        complex_near(actual.electric_v_per_m, forward + backward, 2e-14);
                        complex_near(
                            actual.conjugate_magnetic_a_per_m,
                            (forward - backward).scale(yg),
                            8e-17,
                        );
                    }
                }
            }
        }
        let ci = 0.01;
        let gamma = (n * n * (1.0 - ci * ci) - 1.0).sqrt();
        assert!(gamma > 0.0);
        for polarization in [Polarization::Te, Polarization::Tm] {
            let result = reverse.response(omega, ci, polarization).unwrap();
            near(result.reflected_power_fraction, 1.0, 3e-15);
            assert_eq!(result.transmitted_power_fraction, 0.0);
            near(
                result.transmitted_normal_wavenumber_rad_m.im / (omega / VACUUM_SPEED_M_S),
                -gamma,
                2e-14,
            );
        }
    }
}

#[test]
fn g1_original_gas_multimode_beams_preserve_phase_power_and_zero_true_heat() {
    let (_, model) = fixture();
    let port = port();
    let magnetic = magnetic(MU);
    let numbers = [1.7e25, 0.6e25];
    let modes = [
        (2_usize, C64::new(0.2, -0.1), C64::new(0.1, 0.3)),
        (3, C64::new(0.8, 0.2), C64::new(-0.5, 0.4)),
    ];
    for lambda in [401.25, 1550.0] {
        let omega = FREQUENCY_WAVELENGTH_NM / lambda;
        let k0 = omega / VACUUM_SPEED_M_S;
        let wavelength_m = lambda * 1e-9;
        let law = MolecularGasOpticalResponse::try_lossless_partial_densities(
            &model,
            &numbers,
            PartialDensityBasis::MoleculesPerM3,
            &magnetic,
            omega,
        )
        .unwrap();
        let n = molecular_reference(numbers, lambda, true).0;
        let grid = AngularSpectrumGrid {
            width: 16,
            height: 1,
            spacing_m: [0.625 * wavelength_m, wavelength_m],
        };
        let input: Vec<_> = (0..grid.width)
            .map(|j| {
                let mut e = [C64::ZERO; 2];
                for &(bin, ex, ey) in &modes {
                    let p =
                        phase(-core::f64::consts::TAU * bin as f64 * j as f64 / grid.width as f64);
                    e[0] = e[0] + ex * p;
                    e[1] = e[1] + ey * p;
                }
                e
            })
            .collect();
        let stack = PlanarDielectricLayers::try_new(
            &port,
            DielectricResponse::MolecularGasOptical(law),
            omega,
            &[],
        )
        .unwrap();
        let beam =
            PlanarAngularSpectrum::prepare(&stack, grid, &input, beam_limits(), || false).unwrap();
        let area = grid.width as f64 * grid.spacing_m[0] * grid.spacing_m[1];
        let mut expected_power = 0.0;
        for &(bin, ex, ey) in &modes {
            let kx = core::f64::consts::TAU * bin as f64 / (grid.width as f64 * grid.spacing_m[0]);
            let qi = (1.0 - (kx / k0).powi(2)).sqrt();
            let qt = (n * n - (kx / k0).powi(2)).sqrt();
            for (polarization, amplitude) in [(Polarization::Tm, ex), (Polarization::Te, ey)] {
                let yi = admittance(polarization, 1.0, qi, MU);
                let yt = admittance(polarization, n, qt, MU);
                let t = 2.0 * yi / (yi + yt);
                expected_power += 0.5 * area * yt * t * t * amplitude.norm_sq();
            }
        }
        near(
            beam.power().transmitted_w,
            expected_power,
            7e-15 * expected_power,
        );
        let z = 11.0 * wavelength_m;
        let field = beam
            .field(AngularSpectrumSide::Transmitted, z, || false)
            .unwrap();
        for j in 0..grid.width {
            let x = j as f64 * grid.spacing_m[0];
            let mut expected = [C64::ZERO; 3];
            for &(bin, ex, ey) in &modes {
                let kx =
                    core::f64::consts::TAU * bin as f64 / (grid.width as f64 * grid.spacing_m[0]);
                let qi = (1.0 - (kx / k0).powi(2)).sqrt();
                let qt = (n * n - (kx / k0).powi(2)).sqrt();
                let ex = ex.scale(2.0 * qt / (qt + n * n * qi));
                let ey = ey.scale(2.0 * qi / (qi + qt));
                let p = phase(-kx * x - k0 * qt * z);
                let e = [ex, ey, ex.scale(-kx / (k0 * qt))];
                for c in 0..3 {
                    expected[c] = expected[c] + e[c] * p;
                }
            }
            for c in 0..3 {
                complex_near(field.electric_v_per_m()[j][c], expected[c], 9e-14);
            }
        }
        let scattering = model
            .at_wavelength(lambda, PartialDensityBasis::MoleculesPerM3)
            .unwrap()
            .scattering(&numbers)
            .unwrap();
        assert!(scattering.scattering_per_m > 0.0);
        assert_eq!(
            beam.transmitted_absorbed_power_w(0.0, f64::INFINITY, || false)
                .unwrap(),
            0.0
        );
        let heat = beam
            .transmitted_absorption_field(z, beam_limits(), || false)
            .unwrap();
        assert!(heat.absorbed_power_density_w_m3().iter().all(|q| *q == 0.0));
    }
}

#[test]
fn g3_prepared_frequencies_share_original_material_but_no_causal_reservoir() {
    let (_, model) = fixture();
    let port = port();
    let magnetic = magnetic(MU);
    let numbers = [1.7e25, 0.6e25];
    let make_stack = |lambda| {
        let omega = FREQUENCY_WAVELENGTH_NM / lambda;
        let law = MolecularGasOpticalResponse::try_lossless_partial_densities(
            &model,
            &numbers,
            PartialDensityBasis::MoleculesPerM3,
            &magnetic,
            omega,
        )
        .unwrap();
        PlanarDielectricLayers::try_new(
            &port,
            DielectricResponse::Lossless(&port),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::MolecularGasOptical(law),
                thickness_m: 2e-6,
            }],
        )
        .unwrap()
    };
    let stacks = [make_stack(500.0), make_stack(1550.0)];
    assert!(stacks[0].same_sources_and_geometry(&stacks[1]));
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1e-6; 2],
    };
    let beams: Vec<_> = stacks
        .iter()
        .map(|stack| {
            PlanarAngularSpectrum::prepare(
                stack,
                grid,
                &[[C64::ONE, C64::ZERO]],
                beam_limits(),
                || false,
            )
            .unwrap()
        })
        .collect();
    let exposure = AngularSpectrumExposure::try_new(
        &beams[0],
        &beams,
        AngularSpectrumExposureLimits {
            maximum_lines: 2,
            maximum_point_line_pairs: 8,
            maximum_point_frequency_pairs: 8,
            maximum_workspace_bytes: 2_000_000,
            maximum_fft_work: 100_000,
        },
        || false,
    )
    .unwrap();
    assert_eq!(exposure.frequency_count(), 2);
    let result = exposure.prepare_layer(
        0,
        AngularSpectrumExposureLayerLimits {
            maximum_frequency_groups: 2,
            maximum_frequency_layer_mode_pairs: 8,
            maximum_point_line_pairs: 8,
            maximum_point_reservoir_frequency_pairs: 100,
            maximum_workspace_bytes: 2_000_000,
            maximum_fft_work: 100_000,
        },
        || false,
    );
    match result {
        Err(error) => assert!(
            error
                .to_string()
                .contains("original causal material reservoir")
        ),
        Ok(_) => panic!("finite-band gas dispersion is not a causal reservoir"),
    }
}
