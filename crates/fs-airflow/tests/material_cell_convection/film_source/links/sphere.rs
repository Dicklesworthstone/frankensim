//! G1/G3: one source law and one original thermal link exchange sphere heat.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::links::observer::SourceLinkConvectionObserver;

fn sphere_link() -> SourceConvectiveLink {
    SourceConvectiveLink {
        property_temperature: ConvectionPropertyTemperature::Reservoir,
        law: SourceCellConvectionLaw::ForcedSphere { speed_m_s: 0.1 },
        ..link()
    }
}

fn h(wall: f64, bulk: f64) -> f64 {
    let mu = 0.001 * (1.0 - 0.01 * (bulk - 300.0));
    let wall_mu = 0.001 * (1.0 - 0.01 * (wall - 300.0));
    let k = 0.5 + 0.002 * (bulk - 300.0);
    let re = 1000.0 * 0.1 * 0.1 / mu;
    let pr = mu * FLUID_CP / k;
    (2.0 + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0)) * pr.powf(0.4) * (mu / wall_mu).powf(0.25))
        * k
        / 0.1
}

fn endpoint(dt: f64) -> [f64; 2] {
    let (mut lo, mut hi) = (294.0, 310.0);
    for _ in 0..80 {
        let wall = 0.5 * lo + 0.5 * hi;
        let bulk = 290.0 + CP / FLUID_CP * (310.0 - wall);
        let q = LINK_AREA * h(wall, bulk) * (wall - bulk);
        if MASS * CP * (wall - 310.0) + dt * q > 0.0 {
            hi = wall;
        } else {
            lo = wall;
        }
    }
    let wall = 0.5 * lo + 0.5 * hi;
    [wall, 290.0 + CP / FLUID_CP * (310.0 - wall)]
}

#[test]
fn g1_g3_sphere_link_requeries_both_endpoints_and_closes_original_energy() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    for reversed in [false, true] {
        let thermal = pair(&chart, &caloric, reversed, 0.0);
        let saved = thermal.cells().clone();
        let dt = 0.2;
        let expected = endpoint(dt);
        let step = BackwardEulerStep::new(&thermal, dt, || Ok(())).unwrap();
        let increments: Vec<_> = step
            .cell_ids()
            .iter()
            .map(|&id| {
                if id == IDS[0] {
                    MASS * CP * (expected[0] - 310.0)
                } else {
                    MASS * FLUID_CP * (expected[1] - 290.0)
                }
            })
            .collect();
        let trial = step
            .evaluate(
                &increments,
                |state| {
                    let source = evaluate_source_link_convection(
                        state,
                        &[sphere_link()],
                        1,
                        |_, _, t| law.at_temperature(t),
                        || false,
                    )?;
                    Ok::<_, SourceCellConvectionError<MaterialStatePointError>>(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: source.conductances_w_k().clone(),
                        payload: source,
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, result, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        let sample = &result.samples()[0];
        near(sample.property_temperature_k, expected[1], 1e-12);
        near(
            sample.surface.as_ref().unwrap().thermal.temperature_k(),
            expected[0],
            1e-12,
        );
        near(
            sample.conductance_w_k,
            LINK_AREA * h(expected[0], expected[1]),
            1e-9,
        );
        near(
            sample.heat_into_fluid_w,
            MASS * FLUID_CP * (expected[1] - 290.0) / dt,
            1e-8,
        );
        for (id, temperature) in IDS.into_iter().zip(expected) {
            near(next.temperature_k(id).unwrap(), temperature, 1e-10);
        }
        near(next.energy_j(), thermal.energy_j(), 1e-8);
        near(audit.external_heat_j, 0.0, 1e-12);
        assert_eq!(thermal.cells(), &saved);
        assert_eq!(
            sample.material.material().card_identity(),
            card.content_hash()
        );
        assert_eq!(
            sample
                .surface
                .as_ref()
                .unwrap()
                .material
                .material()
                .card_identity(),
            card.content_hash()
        );
    }
}

#[test]
fn g1_native_sphere_link_uses_the_existing_timestep_and_conservation_owner() {
    let chart = curve();
    let (card, caloric) = fluid_source(1.0);
    let law = liquid_law(&card);
    let thermal = pair(&chart, &caloric, false, 0.0);
    let bindings = [sphere_link()];
    let observer = SourceLinkConvectionObserver::new(
        &bindings,
        1,
        |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
            law.at_temperature(t)
        },
    );
    let (next, accepted, audit) = thermal
        .advance_observed(
            &observer,
            0.1,
            StepConfig {
                max_dt_s: 0.005,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
    assert!(next.temperature_k(IDS[0]).unwrap() < 310.0);
    assert!(next.temperature_k(IDS[1]).unwrap() > 290.0);
    near(next.energy_j(), thermal.energy_j(), 1e-8);
    near(audit.external_heat_j, 0.0, 1e-12);
    let sample = &accepted.sample().unwrap().samples()[0];
    near(
        sample.surface.as_ref().unwrap().thermal.temperature_k(),
        next.temperature_k(IDS[0]).unwrap(),
        1e-12,
    );
    near(
        sample.thermal.temperature_k(),
        next.temperature_k(IDS[1]).unwrap(),
        1e-12,
    );
    assert_eq!(next.links()[&5].conductance_w_k, 0.0);
    assert!(observer.sample().is_none());
}
