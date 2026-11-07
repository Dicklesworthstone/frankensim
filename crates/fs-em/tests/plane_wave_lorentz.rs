#![cfg(feature = "plane-wave")]
//! G1/G3 resonant dielectric scattering from synthetic source cards.
//! These regressions cover constitutive coupling and coherent wave accounting;
//! they do not establish measured-material accuracy or transient Maxwell behavior.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{
    QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

type Entry = (&'static str, QuantitySpec, f64);

fn angular_frequency() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular),
        ValueForm::Static,
    ))
}

fn card(entries: &[Entry], specimen: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, quantity, value) in entries {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic Lorentz wave fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: specimen.into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[cfg(test)]
mod integration {
    use super::*;
    use fs_em::plane_wave::{
        DielectricResponse, PlanarDielectricInterface, Polarization,
        layers::{DielectricLayer, PlanarDielectricLayers},
    };
    use fs_matdb::QueryPoint;
    use fs_material::{
        conductor::{ElectricalPropertyBasis, ResolvedConductorState},
        field_storage::{
            ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
            MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
            dielectric::lorentz::{LorentzPoleProperties, ResolvedLorentzMaterial},
        },
        state_point::{
            MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
            ScalarPropertyRequirement, resolve_material_state_point,
        },
    };
    use fs_math::c64::C64;
    use fs_qty::Dims;

    const E: f64 = 1e-11;
    const M: f64 = 1e-6;
    const NATURAL: f64 = 1e9;
    const SIGMA: f64 = 0.02;

    fn source(entries: &[Entry], specimen: &str) -> ResolvedMaterialStatePoint {
        let requirements = entries
            .iter()
            .map(|&(name, quantity, _)| {
                ScalarPropertyRequirement::try_with_quantity(
                    name,
                    quantity,
                    ScalarAdmissibility::NonNegative,
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        resolve_material_state_point(
            &card(entries, specimen),
            &QueryPoint::new(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }

    fn outside() -> ResolvedMaterialStatePoint {
        source(
            &[
                (EPS_NAME, QuantitySpec::dimensional(EPS), E),
                (MU_NAME, QuantitySpec::dimensional(MU), M),
            ],
            "synthetic surrounding dielectric",
        )
    }

    fn resonant_source(damping: f64, specimen: &str) -> ResolvedMaterialStatePoint {
        let basis = ElectricalPropertyBasis::Conductivity;
        source(
            &[
                ("epsilon_inf", QuantitySpec::dimensional(EPS), 2.0 * E),
                ("delta", QuantitySpec::dimensional(EPS), 4.0 * E),
                ("omega0", angular_frequency(), NATURAL),
                (
                    "gamma",
                    QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                    damping,
                ),
                (MU_NAME, QuantitySpec::dimensional(MU), M),
                (
                    basis.property_name(),
                    QuantitySpec::dimensional(basis.dimensions()),
                    SIGMA,
                ),
            ],
            specimen,
        )
    }

    fn lorentz(source: ResolvedMaterialStatePoint) -> ResolvedLorentzMaterial {
        ResolvedLorentzMaterial::try_from_material(
            source,
            "epsilon_inf",
            &[LorentzPoleProperties {
                permittivity_increment: "delta".into(),
                angular_frequency: "omega0".into(),
                damping_rate: "gamma".into(),
            }],
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }

    fn complex_close(actual: C64, expected: C64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:?} != {expected:?}"
        );
    }

    #[test]
    fn resonance_and_independent_ohmic_loss_preserve_coherent_passive_scattering() {
        let ambient = outside();
        let material = resonant_source(0.03 * NATURAL, "resonant specimen A");
        let polarization = lorentz(material.clone());
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            material,
            ElectricalPropertyBasis::Conductivity,
        )
        .unwrap();
        assert_eq!(
            polarization.material().identity(),
            conductor.material().identity()
        );
        for ratio in [0.5, 0.99, 1.0, 1.15, 2.0] {
            let omega = ratio * NATURAL;
            let bare_epsilon = polarization.permittivity_f_m(omega).unwrap();
            if ratio == 1.15 {
                assert!(bare_epsilon.re < 0.0);
            }
            for (response, expected_epsilon) in [
                (DielectricResponse::Lorentz(&polarization), bare_epsilon),
                (
                    DielectricResponse::LorentzWithOhmic {
                        dielectric: &polarization,
                        conductor: &conductor,
                    },
                    C64::new(bare_epsilon.re, bare_epsilon.im - SIGMA / omega),
                ),
            ] {
                let interface = PlanarDielectricInterface::try_new(&ambient, response).unwrap();
                let layer = DielectricLayer {
                    response,
                    thickness_m: 0.12,
                };
                let film = PlanarDielectricLayers::try_new(
                    &ambient,
                    DielectricResponse::Lossless(&ambient),
                    omega,
                    &[layer],
                )
                .unwrap();
                let split = PlanarDielectricLayers::try_new(
                    &ambient,
                    DielectricResponse::Lossless(&ambient),
                    omega,
                    &[DielectricLayer {
                        response,
                        thickness_m: 0.04,
                    }; 3],
                )
                .unwrap();
                for mode in [Polarization::Te, Polarization::Tm] {
                    let half_space = interface.response(omega, 0.7, mode).unwrap();
                    complex_close(
                        half_space.transmitted_permittivity_f_m,
                        expected_epsilon,
                        3e-14 * expected_epsilon.abs(),
                    );
                    assert!(half_space.transmitted_normal_wavenumber_rad_m.im <= 0.0);
                    assert!(half_space.reflected_power_fraction >= 0.0);
                    assert!(half_space.transmitted_power_fraction >= 0.0);
                    // This is flux entering the lossy half-space, not a claim
                    // that this power reaches an observer at infinite depth.
                    close(
                        half_space.reflected_power_fraction + half_space.transmitted_power_fraction,
                        1.0,
                        3e-13,
                    );
                    let full = film.response(0.7, mode).unwrap();
                    let subdivided = split.response(0.7, mode).unwrap();
                    assert!(full.absorbed_power_fraction > 0.0);
                    assert!(full.reflected_power_fraction >= 0.0);
                    assert!(full.transmitted_power_fraction >= 0.0);
                    close(
                        full.reflected_power_fraction
                            + full.transmitted_power_fraction
                            + full.absorbed_power_fraction,
                        1.0,
                        3e-13,
                    );
                    complex_close(
                        full.reflected_tangential_e,
                        subdivided.reflected_tangential_e,
                        2e-12,
                    );
                    complex_close(
                        full.transmitted_tangential_e,
                        subdivided.transmitted_tangential_e,
                        2e-12,
                    );
                    close(
                        full.absorbed_power_fraction,
                        subdivided.absorbed_power_fraction,
                        2e-12,
                    );
                }
            }
        }
    }

    #[test]
    fn equal_coefficients_from_different_complete_material_sources_cannot_be_combined() {
        let ambient = outside();
        let polarization = lorentz(resonant_source(0.03 * NATURAL, "specimen A"));
        let other = lorentz(resonant_source(0.03 * NATURAL, "specimen B"));
        assert_eq!(
            polarization.permittivity_f_m(NATURAL).unwrap(),
            other.permittivity_f_m(NATURAL).unwrap()
        );
        assert_ne!(
            polarization.material().identity(),
            other.material().identity()
        );
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            other.material().clone(),
            ElectricalPropertyBasis::Conductivity,
        )
        .unwrap();
        let mismatched = DielectricResponse::LorentzWithOhmic {
            dielectric: &polarization,
            conductor: &conductor,
        };
        assert!(PlanarDielectricInterface::try_new(&ambient, mismatched).is_err());
        assert!(
            PlanarDielectricLayers::try_new(
                &ambient,
                DielectricResponse::Lossless(&ambient),
                NATURAL,
                &[DielectricLayer {
                    response: mismatched,
                    thickness_m: 0.12
                }],
            )
            .is_err()
        );
    }

    #[test]
    fn lossless_negative_permittivity_reflects_and_tunnels_without_invented_heat() {
        let ambient = outside();
        let polarization = lorentz(resonant_source(0.0, "lossless resonant specimen"));
        let response = DielectricResponse::Lorentz(&polarization);
        let interface = PlanarDielectricInterface::try_new(&ambient, response).unwrap();
        let omega = 1.15 * NATURAL;
        assert!(polarization.permittivity_f_m(omega).unwrap().re < 0.0);
        let film = PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            omega,
            &[DielectricLayer {
                response,
                thickness_m: 0.12,
            }],
        )
        .unwrap();
        for mode in [Polarization::Te, Polarization::Tm] {
            let half_space = interface.response(omega, 1.0, mode).unwrap();
            close(half_space.reflected_power_fraction, 1.0, 2e-13);
            close(half_space.transmitted_power_fraction, 0.0, 2e-13);
            assert!(half_space.transmitted_normal_wavenumber_rad_m.im < 0.0);
            let finite = film.response(1.0, mode).unwrap();
            assert!(
                finite.transmitted_power_fraction > 0.0 && finite.transmitted_power_fraction < 1.0
            );
            close(finite.absorbed_power_fraction, 0.0, 3e-13);
            close(
                finite.reflected_power_fraction + finite.transmitted_power_fraction,
                1.0,
                3e-13,
            );
            assert!(interface.response(NATURAL, 1.0, mode).is_err());
        }
        assert!(
            PlanarDielectricLayers::try_new(
                &ambient,
                DielectricResponse::Lossless(&ambient),
                NATURAL,
                &[DielectricLayer {
                    response,
                    thickness_m: 0.12
                }],
            )
            .is_err()
        );
    }
}
