use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoatingIorSample, CoherentPeriodicCoating, ConductorDataStatus,
    ConductorSource,
};
use crate::dielectric::CauchyIor;

fn source(tag: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative).unwrap()
}

fn opaque(layers: &[CoatingLayer<'_>], substrate: ConductorOptics) -> CoatingThermalEmission {
    let temperatures = vec![1500.0; layers.len()];
    CoatingThermalEmission::try_opaque(
        layers,
        &temperatures,
        substrate,
        1800.0,
        ContentHash([91; 32]),
    )
    .unwrap()
}

#[test]
fn smooth_binding_preserves_original_bulk_knots_source_and_gap_policy() {
    let representative = ConductorOptics::representative_tungsten();
    let emission = opaque(&[], representative);
    assert!(emission.matches_smooth_conductor(&representative));
    // Numerically identical reconstruction still has different provenance.
    let bulk = ConductorOptics::try_new(representative.samples(), representative.source()).unwrap();
    assert!(!emission.matches_smooth_conductor(&bulk));
    let emission = opaque(&[], bulk);
    assert!(emission.matches_smooth_conductor(&bulk));
    assert!(
        !emission.matches_smooth_conductor(
            &ConductorOptics::try_new(bulk.samples(), source(92)).unwrap()
        )
    );
    let changed_status = ConductorSource::try_new(
        bulk.source().identity(),
        ConductorDataStatus::CallerAssertedMeasured,
    )
    .unwrap();
    assert!(!emission.matches_smooth_conductor(
        &ConductorOptics::try_new(bulk.samples(), changed_status).unwrap()
    ));
    let changed_gap =
        ConductorOptics::try_new_in_band(bulk.samples(), 200.0, bulk.source()).unwrap();
    assert!(!emission.matches_smooth_conductor(&changed_gap));
    let mut changed_knots = bulk.samples();
    let original = changed_knots[4];
    changed_knots[4] = crate::conductor::ConductorIorSample::try_new(
        original.wavelength_nm(),
        original.eta() + 0.01,
        original.k(),
    )
    .unwrap();
    assert!(!emission.matches_smooth_conductor(
        &ConductorOptics::try_new(changed_knots, bulk.source()).unwrap()
    ));
}

#[test]
fn smooth_binding_matches_exact_single_film_and_explicit_substrate_override() {
    let bulk = ConductorOptics::representative_tungsten();
    let film = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.7).unwrap(),
        80.0e-9,
        source(93),
    )
    .unwrap();
    let changed_film =
        CoherentDielectricCoating::try_new(film.ior(), film.thickness_m(), source(94)).unwrap();
    let coated = bulk.with_coating(film);
    let emission = opaque(
        &[CoatingLayer::Dielectric(&film)],
        bulk.with_coating(changed_film),
    );
    assert!(emission.matches_smooth_conductor(&coated));
    assert!(!emission.matches_smooth_conductor(&bulk));
    assert!(!emission.matches_smooth_conductor(&bulk.with_coating(changed_film)));
    assert!(opaque(&[], coated).matches_smooth_conductor(&bulk));
    assert!(!opaque(&[], coated).matches_smooth_conductor(&coated));

    let absorbing = CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|index| {
            CoatingIorSample::try_new(380.0 + 50.0 * index as f64, 1.7, 0.2).unwrap()
        }),
        80.0e-9,
        source(93),
    )
    .unwrap();
    let emission = opaque(&[CoatingLayer::Absorbing(&absorbing)], bulk);
    assert!(emission.matches_smooth_conductor(&bulk.with_absorbing_coating(absorbing)));
    assert!(!emission.matches_smooth_conductor(&coated));
    let sample = emission.evaluate(550.0, 0.6).unwrap();
    for (pol, retained) in sample
        .polarized_radiance_w_m2_sr_nm()
        .into_iter()
        .enumerate()
    {
        let original = sample.film_radiance_w_m2_sr_nm()[0][pol]
            + sample.substrate_radiance_w_m2_sr_nm().unwrap()[pol];
        assert!((retained - original).abs() <= 2.0 * f64::EPSILON * original);
    }
    let two_films = opaque(
        &[
            CoatingLayer::Dielectric(&film),
            CoatingLayer::Absorbing(&absorbing),
        ],
        bulk,
    );
    assert!(!two_films.matches_smooth_conductor(&bulk.with_absorbing_coating(absorbing)));

    let periodic = CoherentPeriodicCoating::try_new(
        &[
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.7).unwrap(), 0.0, 80.0e-9)
                .unwrap(),
        ],
        1,
        source(93),
    )
    .unwrap();
    assert!(
        !opaque(&[CoatingLayer::Dielectric(&film)], bulk)
            .matches_smooth_conductor(&bulk.with_periodic_coating(periodic))
    );
    let transmitting = CoatingThermalEmission::try_transmitting(
        &[CoatingLayer::Absorbing(&absorbing)],
        &[1500.0],
        CauchyIor::try_constant(1.0).unwrap().into(),
        false,
        ContentHash([91; 32]),
    )
    .unwrap();
    assert!(!transmitting.matches_smooth_conductor(&bulk.with_absorbing_coating(absorbing)));
}

#[test]
fn retained_thermal_channels_preserve_the_weak_mode_of_a_brewster_bragg_stack() {
    // At this transverse wave number every n=1/n=2 film boundary is at
    // Brewster incidence for p. The quarter-wave sequence rejects s strongly,
    // while p still reaches the hot substrate. Original positive layer work
    // resolves s even after its contribution vanishes from rounded I and Q.
    let wavelength = 550.0;
    let cosine = 1.0 / 5.0_f64.sqrt();
    let low = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.0).unwrap(),
        wavelength * 1e-9 / (4.0 * cosine),
        source(95),
    )
    .unwrap();
    let high = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(2.0).unwrap(),
        wavelength * 1e-9 / (16.0 * cosine),
        source(96),
    )
    .unwrap();
    let layers = core::array::from_fn::<_, 32, _>(|index| {
        CoatingLayer::Dielectric(if index % 2 == 0 { &high } else { &low })
    });
    let emitter = opaque(&layers, ConductorOptics::representative_tungsten());
    let sample = emitter.evaluate(wavelength, cosine).unwrap();
    let [s, p] = sample.polarized_radiance_w_m2_sr_nm();
    assert!(s > 0.0 && p > 0.0);
    assert!(s / p < 1e-18, "s={s:e}, p={p:e}");
    assert_eq!(0.5 * (sample.stokes().i + sample.stokes().q), 0.0);
    assert_eq!(
        sample.polarized_radiance_w_m2_sr_nm(),
        sample.substrate_radiance_w_m2_sr_nm().unwrap()
    );
    assert!(
        sample
            .film_radiance_w_m2_sr_nm()
            .iter()
            .all(|pair| *pair == [0.0; 2])
    );
}
