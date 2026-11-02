//! Original native caloric states drive nonisothermal reciprocal observations.
//! Numerical optical laws below are authored fixtures, not material calibration.
use super::*;
use fs_blake3::ContentHash;
use fs_render::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack,
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource,
};
use fs_render::dielectric::{CauchyIor, DielectricCoatingStack, DielectricIor};
use fs_render::thermal_surface::heating::feedback::SurfaceConductorLaw;
use fs_render::thermal_surface::heating::{
    CoatingOpticsModel, NativeCoatingEmission, NativeCoatingEmissionError, StatefulCoatingLayer,
    SubstrateOpticsModel, SurfaceCoatingLaw, SurfaceHeatingError, SurfaceHeatingReceiver,
};
use std::sync::Arc;

const SOURCE: ContentHash = ContentHash([249; 32]);
fn authority() -> ConductorSource {
    ConductorSource::try_new(ContentHash([247; 32]), ConductorDataStatus::Representative).unwrap()
}
fn film(n: f64, k: f64, thickness_nm: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        std::array::from_fn(|i| CoatingIorSample::try_new(380.0 + 50.0 * i as f64, n, k).unwrap()),
        thickness_nm * 1e-9,
        authority(),
    )
    .unwrap()
}
fn varying_film(i: usize, temperature: f64) -> CoherentAbsorbingCoating {
    let d = temperature - 1400.0;
    if i == 0 {
        film(1.6 + 0.00015 * d, 0.05 + 0.00005 * d, 75.0)
    } else {
        film(2.1 - 0.0001 * d, 0.12 + 0.00007 * d, 110.0)
    }
}
fn substrate(temperature: f64) -> ConductorOptics {
    let d = temperature - 1400.0;
    ConductorOptics::try_new(
        std::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.4 + 0.0001 * d, 3.0 - 0.0002 * d)
                .unwrap()
        }),
        authority(),
    )
    .unwrap()
}
fn receiver(id: u64) -> SurfaceHeatingReceiver {
    SurfaceHeatingReceiver {
        thermal_cell: CellId(id),
        material_id: MaterialId(7),
        temperature_bounds_k: (1400.0, 3400.0),
    }
}
fn other_medium() -> DielectricIor {
    DielectricIor::Cauchy(CauchyIor::try_constant(1.5).unwrap())
}
struct HotSensible;
impl ThermalLaw for HotSensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (0.0..=2000.0).contains(&h) {
            Ok(1400.0 + h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static HOT: HotSensible = HotSensible;
struct HotLatent;
impl ThermalLaw for HotLatent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (10.0..=110.0).contains(&h) {
            Ok(2200.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([248; 32])
    }
}
static HOT_LATENT: HotLatent = HotLatent;

// Independent direct SI Planck formula; chosen optical-band fixtures are far
// from its overflow/underflow limits. Production uses its original log owner.
fn blackbody(temperature: f64, wavelength_nm: f64) -> f64 {
    let wavelength = wavelength_nm * 1e-9;
    let h = 6.62607015e-34;
    let c = 299792458.0;
    let k = 1.380649e-23;
    2.0 * h * c * c / wavelength.powi(5) / (h * c / (wavelength * k * temperature)).exp_m1() * 1e-9
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 6e-12 * b.abs().max(1e-100),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn independent_native_film_and_substrate_states_determine_emitted_stokes() {
    let mut cells: Vec<_> = (0..3)
        .map(|i| cell(11 + i as u64, (i + 1) as f64, 1.0, [20.0, 5.0, 10.0][i]))
        .collect();
    for (cell, h) in cells.iter_mut().zip([500.0, 1100.0, 700.0]) {
        cell.specific_h_j_kg = h;
    }
    let native = thermal(cells, &HOT, vec![]);
    let films: Vec<_> = (0..2)
        .map(|i| StatefulCoatingLayer {
            receiver: receiver(11 + i as u64),
            optics: CoatingOpticsModel::Law(Arc::new(
                SurfaceCoatingLaw::try_temperature(
                    authority(),
                    &[
                        (1400.0, varying_film(i, 1400.0)),
                        (3400.0, varying_film(i, 3400.0)),
                    ],
                    2,
                )
                .unwrap(),
            )),
        })
        .collect();
    let metal = SubstrateOpticsModel::Law(Arc::new(
        SurfaceConductorLaw::try_temperature(
            authority(),
            &[(1400.0, substrate(1400.0)), (3400.0, substrate(3400.0))],
            2,
        )
        .unwrap(),
    ));
    let prepare = |thermal: &DepositionSimulation<'_>| {
        NativeCoatingEmission::try_opaque(thermal, &films, receiver(13), &metal, SOURCE, 6, || {
            Ok(())
        })
        .unwrap()
    };
    let observation = prepare(&native);
    assert_eq!(observation.time_s(), 0.0);
    assert_eq!(observation.identity(), prepare(&native).identity());
    let saved = observation.emitter().evaluate(550.0, 0.7).unwrap().stokes();
    let energy = native.energy_j();
    let (heated, _) = native.advance_to(0.4, config(0.1), || Ok(())).unwrap();
    let later = prepare(&heated);
    assert_eq!(later.time_s(), 0.4);
    assert_ne!(later.identity(), observation.identity());
    assert_ne!(
        later.emitter().evaluate(550.0, 0.7).unwrap().stokes(),
        saved
    );
    assert_eq!(
        observation.emitter().evaluate(550.0, 0.7).unwrap().stokes(),
        saved
    );
    assert_eq!(native.energy_j(), energy);
    assert_eq!(native.time_s(), 0.0);
    for (thermal, snapshot) in [(&native, &observation), (&heated, &later)] {
        let temperatures = std::array::from_fn::<_, 3, _>(|i| {
            thermal.temperature_k(CellId(11 + i as u64)).unwrap()
        });
        let optics = [
            varying_film(0, temperatures[0]),
            varying_film(1, temperatures[1]),
        ];
        let layers = optics.each_ref().map(CoatingLayer::Absorbing);
        let stack = CoherentCoatingStack::try_new(&layers).unwrap();
        for lambda in [480.0, 680.0] {
            for cosine in [0.3, 0.8] {
                let partition = stack
                    .layer_power_balance(&substrate(temperatures[2]), lambda, 1.0, cosine)
                    .unwrap();
                let result = snapshot.emitter().evaluate(lambda, cosine).unwrap();
                let mut expected = [0.0; 2];
                for i in 0..2 {
                    assert_eq!(
                        snapshot.emitter().film_temperature_k(i),
                        Some(temperatures[i])
                    );
                    for p in 0..2 {
                        let radiance = partition.layer_absorbed()[i][p]
                            * blackbody(temperatures[i], lambda)
                            / 2.0;
                        near(result.film_radiance_w_m2_sr_nm()[i][p], radiance);
                        expected[p] += radiance;
                    }
                }
                for p in 0..2 {
                    let radiance =
                        partition.substrate_absorbed[p] * blackbody(temperatures[2], lambda) / 2.0;
                    near(result.substrate_radiance_w_m2_sr_nm().unwrap()[p], radiance);
                    expected[p] += radiance;
                }
                near(result.stokes().i, expected[0] + expected[1]);
                near(result.stokes().q, expected[0] - expected[1]);
                assert_eq!((result.stokes().u, result.stokes().v), (0.0, 0.0));
            }
        }
    }
}

fn latent_film(h: f64) -> CoherentAbsorbingCoating {
    film(1.4 + 0.01 * (h - 10.0), 0.02 + 0.0078 * (h - 10.0), 90.0)
}
#[test]
fn equal_temperature_latent_states_preserve_original_optical_enthalpy_and_shared_rows() {
    let mut cells = vec![cell(11, 1.0, 1.0, 4.0), cell(12, 1.0, 1.0, 2.0)];
    cells[0].specific_h_j_kg = 30.0;
    cells[1].specific_h_j_kg = 90.0;
    let native = thermal(cells, &HOT_LATENT, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([248; 32]),
            (2200.0, 2200.0),
            &[(10.0, latent_film(10.0)), (110.0, latent_film(110.0))],
            2,
        )
        .unwrap(),
    );
    let films: Vec<_> = [11, 12, 11]
        .into_iter()
        .map(|id| StatefulCoatingLayer {
            receiver: receiver(id),
            optics: CoatingOpticsModel::Law(Arc::clone(&law)),
        })
        .collect();
    let prepare = |thermal: &DepositionSimulation<'_>, rows| {
        NativeCoatingEmission::try_transmitting(
            thermal,
            &films,
            other_medium(),
            true,
            SOURCE,
            rows,
            || Ok(()),
        )
    };
    assert!(matches!(
        prepare(&native, 1),
        Err(NativeCoatingEmissionError::Binding(
            SurfaceHeatingError::Budget
        ))
    ));
    let snapshot = prepare(&native, 2).unwrap();
    let (heated, _) = native.advance_to(0.5, config(0.125), || Ok(())).unwrap();
    let updated = prepare(&heated, 2).unwrap();
    for id in [11, 12] {
        assert_eq!(
            native.temperature_k(CellId(id)).unwrap(),
            heated.temperature_k(CellId(id)).unwrap()
        );
    }
    assert_ne!(snapshot.emitter().identity(), updated.emitter().identity());
    for (thermal, prepared) in [(&native, &snapshot), (&heated, &updated)] {
        let optics =
            [11, 12, 11].map(|id| latent_film(thermal.cells()[&CellId(id)].specific_h_j_kg));
        let layers = optics.each_ref().map(CoatingLayer::Absorbing);
        let original = DielectricCoatingStack::try_mixed(&layers)
            .unwrap()
            .layer_power_balance(550.0, true, 0.7, 1.0, 1.5)
            .unwrap();
        let observed = prepared.emitter().evaluate(550.0, 0.7).unwrap();
        assert_eq!(observed.substrate_radiance_w_m2_sr_nm(), None);
        for i in 0..3 {
            for p in 0..2 {
                near(
                    observed.film_radiance_w_m2_sr_nm()[i][p],
                    original.layer_absorbed()[i][p] * blackbody(2200.0, 550.0) / 2.0,
                );
            }
        }
        near(
            observed.stokes().i,
            blackbody(2200.0, 550.0) * (original.absorbed[0] + original.absorbed[1]) / 2.0,
        );
    }
}

#[test]
fn native_emission_domain_identity_and_cancellation_refusals_never_change_heat() {
    let mut cell = cell(11, 1.0, 1.0, 0.0);
    cell.specific_h_j_kg = 30.0;
    let native = thermal(vec![cell.clone()], &HOT_LATENT, vec![]);
    let law = Arc::new(
        SurfaceCoatingLaw::try_enthalpy(
            authority(),
            ContentHash([248; 32]),
            (2200.0, 2200.0),
            &[(10.0, latent_film(10.0)), (110.0, latent_film(110.0))],
            2,
        )
        .unwrap(),
    );
    let valid = StatefulCoatingLayer {
        receiver: receiver(11),
        optics: CoatingOpticsModel::Law(law),
    };
    for mutation in 0..3 {
        let mut wrong = valid.clone();
        match mutation {
            0 => wrong.receiver.thermal_cell = CellId(99),
            1 => wrong.receiver.material_id = MaterialId(99),
            _ => wrong.receiver.temperature_bounds_k = (1400.0, 2199.0),
        }
        assert!(
            NativeCoatingEmission::try_transmitting(
                &native,
                &[wrong],
                other_medium(),
                false,
                SOURCE,
                2,
                || Ok(())
            )
            .is_err()
        );
    }
    let mut other = cell;
    other.specific_h_j_kg = 800.0;
    let wrong_law = thermal(vec![other], &HOT, vec![]);
    assert_eq!(wrong_law.temperature_k(CellId(11)).unwrap(), 2200.0);
    assert!(
        NativeCoatingEmission::try_transmitting(
            &wrong_law,
            &[valid.clone()],
            other_medium(),
            false,
            SOURCE,
            2,
            || Ok(())
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        NativeCoatingEmission::try_transmitting(
            &native,
            &[valid.clone()],
            other_medium(),
            false,
            SOURCE,
            2,
            || {
                polls += 1;
                if polls == 3 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        ),
        Err(NativeCoatingEmissionError::Binding(
            SurfaceHeatingError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert_eq!(native.cells()[&CellId(11)].specific_h_j_kg, 30.0);
    assert_eq!(native.time_s(), 0.0);
    let accepted = NativeCoatingEmission::try_transmitting(
        &native,
        &[valid],
        other_medium(),
        false,
        SOURCE,
        2,
        || Ok(()),
    )
    .unwrap();
    assert!(accepted.emitter().radiance_w_m2_sr_nm(550.0, 0.7).unwrap() > 0.0);
    let empty = thermal(vec![], &HOT, vec![]);
    let bare = NativeCoatingEmission::try_transmitting(
        &empty,
        &[],
        other_medium(),
        false,
        SOURCE,
        0,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(bare.emitter().radiance_w_m2_sr_nm(550.0, 0.7).unwrap(), 0.0);
    assert!(bare.emitter().evaluate(1000.0, 0.7).is_err());
}
