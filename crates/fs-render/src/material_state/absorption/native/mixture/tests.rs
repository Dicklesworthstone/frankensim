//! Independent constituent spectra and temperatures drive real observations.

use super::super::tests::{
    LENGTH, budget, camera, close, direction, geometry, origin, palette, planck, settings, with_cx,
};
use super::*;
use crate::material_state::absorption::tests::{complex_fixture, density_key, point};
use crate::volumes::transport::CellPaletteBudget;
use crate::volumes::transport::camera::VolumeCameraBackground;
use fs_material::state_point::{
    MaterialPropertySelection, SpectralComplexIndexStatePoint, resolve_material_state_point,
};
use std::sync::Arc;

fn recipe<'a>(
    cx: &Cx<'_>,
    components: &'a [MaterialAbsorptionComponent<'a>],
) -> MaterialAbsorptionMixture<'a> {
    MaterialAbsorptionMixture::try_new(cx, components, OpticalMixtureBudget { max_components: 8 })
        .unwrap()
}

fn recipes<'a>(
    cx: &Cx<'_>,
    values: Vec<MaterialAbsorptionMixture<'a>>,
    indices: &[u32],
) -> CellPalette<MaterialAbsorptionMixture<'a>> {
    let budget = CellPaletteBudget {
        max_models: values.len(),
        max_cells: indices.len(),
    };
    CellPalette::try_new(cx, values, Arc::from(indices), budget).unwrap()
}

#[test]
fn g1_original_constituent_laws_and_temperatures_drive_camera_and_deposited_power() {
    with_cx(|cx, _| {
        let (card, cold, cold_optics) = complex_fixture([0.01, 0.03]);
        let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
        let hot_optics = cold_optics.at_query_point(&card, &point(400.0)).unwrap();
        let components = [
            MaterialAbsorptionComponent {
                material: MaterialAbsorptionBinding::try_spectral(
                    &cold,
                    &cold_optics,
                    &density_key(),
                )
                .unwrap(),
                mass_weight: 0.25,
            },
            MaterialAbsorptionComponent {
                material: MaterialAbsorptionBinding::try_spectral(
                    &hot,
                    &hot_optics,
                    &density_key(),
                )
                .unwrap(),
                mass_weight: 0.75,
            },
        ];
        let mixture = recipe(cx, &components);
        assert!(core::ptr::eq(mixture.components(), components.as_slice()));
        let density = [900.0];
        let geometry = geometry(cx, &density, 16);
        let wavelengths = [8000.0, 10000.0, 12000.0];
        let field = MaterialAbsorptionMixtureField::try_new_mixtures(
            cx,
            &geometry,
            recipes(cx, vec![mixture], &[0]),
            wavelengths,
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        let observed = field
            .transport(cx, origin(), direction(), f64::INFINITY)
            .unwrap();
        let image = field.render_camera(cx, &camera(), settings()).unwrap();
        let power = [2.0, 3.0, 4.0];
        let deposited = field
            .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &power, 1)
            .unwrap();
        assert_eq!(deposited.medium_identity, field.content_identity());
        assert_eq!(deposited.wavelengths_nm, wavelengths);
        assert_eq!(deposited.cells.len(), 1);
        let illuminated = field
            .render_camera_with_background(cx, &camera(), settings(), |_, _, _, nodes| {
                assert_eq!(nodes, &wavelengths);
                Ok(VolumeCameraBackground {
                    radiance: power,
                    distance_m: f64::INFINITY,
                })
            })
            .unwrap();
        for (lane, (wavelength, k)) in wavelengths.into_iter().zip([0.01, 0.02, 0.03]).enumerate() {
            let a = 4.0 * core::f64::consts::PI * k / (wavelength * 1e-9 * 1000.0);
            let b = 4.0 * core::f64::consts::PI * 2.0 * k / (wavelength * 1e-9 * 800.0);
            let total = 0.25 * a + 0.75 * b;
            let source = (0.25 * a * planck(300.0, wavelength)
                + 0.75 * b * planck(400.0, wavelength))
                / total;
            let tau = 900.0 * total * LENGTH;
            let transmission = (-tau).exp();
            let emitted = source * -(-tau).exp_m1();
            close(observed.transmittance[lane], transmission);
            close(observed.emission_w_m2_sr_nm[lane], emitted);
            assert_eq!(
                image.pixels()[0].emitted_radiance_w_m2_sr_nm()[lane],
                observed.emission_w_m2_sr_nm[lane]
            );
            close(
                illuminated.radiance_w_m2_sr_nm()[0][lane],
                emitted + transmission * power[lane],
            );
            close(
                deposited.absorbed_power_w[lane],
                power[lane] * -(-tau).exp_m1(),
            );
            close(
                deposited.transmitted_power_w[lane],
                power[lane] * transmission,
            );
            close(
                deposited.absorbed_power_w[lane] + deposited.transmitted_power_w[lane],
                power[lane],
            );
            assert_eq!(deposited.scattered_power_w[lane], 0.0);
            if lane == 1 {
                assert!(
                    (source / planck(375.0, wavelength) - 1.0).abs() > 0.05,
                    "a mass-weighted temperature would replace independent thermal sources"
                );
            }
        }
    });
}

#[test]
fn g1_mass_loadings_are_not_normalized_and_common_temperature_reuses_native_source() {
    with_cx(|cx, _| {
        let (_, a, oa) = complex_fixture([0.01, 0.03]);
        let (_, b, ob) = complex_fixture([0.02, 0.06]);
        let a = MaterialAbsorptionBinding::try_spectral(&a, &oa, &density_key()).unwrap();
        let b = MaterialAbsorptionBinding::try_spectral(&b, &ob, &density_key()).unwrap();
        let components = [
            MaterialAbsorptionComponent {
                material: a,
                mass_weight: 2.0,
            },
            MaterialAbsorptionComponent {
                material: b,
                mass_weight: 3.0,
            },
        ];
        let one = [MaterialAbsorptionComponent {
            material: a,
            mass_weight: 1.0,
        }];
        let density = [1000.0];
        let geometry = geometry(cx, &density, 16);
        let field = MaterialAbsorptionMixtureField::try_new_mixtures(
            cx,
            &geometry,
            recipes(cx, vec![recipe(cx, &components)], &[0]),
            [10000.0],
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        let tau = 8.0 * 4.0 * core::f64::consts::PI * 0.02;
        let result = field
            .transport(cx, origin(), direction(), f64::INFINITY)
            .unwrap();
        close(result.transmittance[0], (-tau).exp());
        close(
            result.emission_w_m2_sr_nm[0],
            planck(300.0, 10000.0) * -(-tau).exp_m1(),
        );
        let one = MaterialAbsorptionMixtureField::try_new_mixtures(
            cx,
            &geometry,
            recipes(cx, vec![recipe(cx, &one)], &[0]),
            [10000.0],
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        let original = MaterialAbsorptionField::try_new(
            cx,
            &geometry,
            palette(cx, vec![a], &[0]),
            [10000.0],
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        assert_eq!(
            one.transport(cx, origin(), direction(), f64::INFINITY),
            original.transport(cx, origin(), direction(), f64::INFINITY)
        );
    });
}

#[test]
fn g0_zero_weight_sources_retain_admission_and_budget_but_do_not_require_temperature() {
    with_cx(|cx, _| {
        let (card, accepted, optical) = complex_fixture([0.01, 0.03]);
        let conditions = point(300.0)
            .with_quantity(
                "other_temperature",
                accepted.axis_quantities()["thermal"],
                350.0,
            )
            .unwrap();
        let requirements = accepted
            .properties()
            .iter()
            .map(|p| p.requirement().clone())
            .collect::<Vec<_>>();
        let ambiguous = resolve_material_state_point(
            &card,
            &conditions,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let ambiguous_optics = SpectralComplexIndexStatePoint::try_from_resolved(
            ambiguous.clone(),
            optical.bindings(),
        )
        .unwrap();
        let material =
            MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key()).unwrap();
        let no_temperature =
            MaterialAbsorptionBinding::try_spectral(&ambiguous, &ambiguous_optics, &density_key())
                .unwrap();
        assert_eq!(no_temperature.accepted_temperature_k(), None);
        let components = [
            MaterialAbsorptionComponent {
                material,
                mass_weight: 1.0,
            },
            MaterialAbsorptionComponent {
                material: no_temperature,
                mass_weight: 0.0,
            },
        ];
        let mixture = recipe(cx, &components);
        let density = [1000.0];
        let geometry = geometry(cx, &density, 16);
        let field = MaterialAbsorptionMixtureField::try_new_mixtures(
            cx,
            &geometry,
            recipes(cx, vec![mixture], &[0]),
            [10000.0],
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        assert!(!field.is_vacuum());
        assert_eq!(
            field.materials().values()[0].components()[1]
                .material
                .identity(),
            no_temperature.identity()
        );
        let insufficient = MaterialAbsorptionPreparationBudget {
            max_models: 1,
            max_spectral_evaluations: 1,
        };
        assert!(matches!(
            MaterialAbsorptionMixtureField::try_new_mixtures(
                cx,
                &geometry,
                recipes(cx, vec![mixture], &[0]),
                [10000.0],
                MaterialAbsorptionEmission::None,
                insufficient,
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Budget
            ))
        ));
        let active = [
            components[0],
            MaterialAbsorptionComponent {
                mass_weight: 1.0,
                ..components[1]
            },
        ];
        assert!(matches!(
            MaterialAbsorptionMixtureField::try_new_mixtures(
                cx,
                &geometry,
                recipes(cx, vec![recipe(cx, &active)], &[0]),
                [10000.0],
                MaterialAbsorptionEmission::AcceptedLte,
                budget(),
            ),
            Err(MaterialAbsorptionFieldError::ComponentTemperature {
                model: 0,
                component: 1
            })
        ));
        let zero = [MaterialAbsorptionComponent {
            material: no_temperature,
            mass_weight: 0.0,
        }];
        let vacuum = MaterialAbsorptionMixtureField::try_new_mixtures(
            cx,
            &geometry,
            recipes(cx, vec![recipe(cx, &zero)], &[0]),
            [10000.0],
            MaterialAbsorptionEmission::AcceptedLte,
            budget(),
        )
        .unwrap();
        assert!(vacuum.is_vacuum());
        let power = vacuum
            .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &[3.0], 0)
            .unwrap();
        assert_eq!(power.transmitted_power_w, [3.0]);
        assert!(power.cells.is_empty());
        let reversed = [components[1], components[0]];
        assert_ne!(
            recipe(cx, &reversed).content_identity(),
            mixture.content_identity()
        );
    });
}

#[test]
fn g0_mixture_weights_source_support_range_and_cancellation_refuse_atomically() {
    with_cx(|cx, gate| {
        let (_, accepted, optical) = complex_fixture([0.01, 0.03]);
        let material =
            MaterialAbsorptionBinding::try_spectral(&accepted, &optical, &density_key()).unwrap();
        for weight in [-1.0, f64::NAN, f64::INFINITY] {
            let components = [MaterialAbsorptionComponent {
                material,
                mass_weight: weight,
            }];
            assert!(matches!(
                MaterialAbsorptionMixture::try_new(
                    cx,
                    &components,
                    OpticalMixtureBudget { max_components: 8 }
                ),
                Err(MaterialAbsorptionFieldError::Transport(
                    VolumeTransportError::Mixture
                ))
            ));
        }
        let mut visible_bindings = optical.bindings().to_vec();
        visible_bindings[0].wavelength_nm = 380.0;
        visible_bindings[1].wavelength_nm = 780.0;
        let visible =
            SpectralComplexIndexStatePoint::try_from_resolved(accepted.clone(), &visible_bindings)
                .unwrap();
        let unsupported = [
            MaterialAbsorptionComponent {
                material,
                mass_weight: 1.0,
            },
            MaterialAbsorptionComponent {
                material: MaterialAbsorptionBinding::try_spectral(
                    &accepted,
                    &visible,
                    &density_key(),
                )
                .unwrap(),
                mass_weight: 0.0,
            },
        ];
        assert!(matches!(
            MaterialAbsorptionMixture::try_new(
                cx,
                &unsupported,
                OpticalMixtureBudget { max_components: 8 }
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Wavelength
            ))
        ));
        let oversized = [MaterialAbsorptionComponent {
            material,
            mass_weight: 0.0,
        }; 65];
        assert!(matches!(
            MaterialAbsorptionMixture::try_new(
                cx,
                &oversized,
                OpticalMixtureBudget { max_components: 65 }
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Budget
            ))
        ));
        let (_, tiny, tiny_optics) = complex_fixture([1e-300; 2]);
        let components = [MaterialAbsorptionComponent {
            material: MaterialAbsorptionBinding::try_spectral(&tiny, &tiny_optics, &density_key())
                .unwrap(),
            mass_weight: 1e-100,
        }];
        let mixture = recipe(cx, &components);
        let density = [1000.0];
        let geometry = geometry(cx, &density, 16);
        assert!(matches!(
            MaterialAbsorptionMixtureField::try_new_mixtures(
                cx,
                &geometry,
                recipes(cx, vec![mixture], &[0]),
                [10000.0],
                MaterialAbsorptionEmission::None,
                budget(),
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Unrepresentable
            ))
        ));
        gate.request();
        assert!(matches!(
            MaterialAbsorptionMixture::try_new(
                cx,
                &components,
                OpticalMixtureBudget { max_components: 8 }
            ),
            Err(MaterialAbsorptionFieldError::Transport(
                VolumeTransportError::Cancelled
            ))
        ));
    });
}
