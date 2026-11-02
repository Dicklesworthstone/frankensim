//! G1/G3: independently solved gap intensities drive both original caloric stores.
use super::*;
use fs_render::dielectric::{DielectricCoatingStack, DielectricIorSample, TabulatedIor};
use fs_render::thermal_surface::heating::coating_exchange::transmitting_pair::*;

const AREA: f64 = 0.004;
const BATHS: [f64; 2] = [430.0, 240.0];
const REVERSE: [bool; 2] = [false, true];

pub(super) fn ir_film(n: f64, k: f64, thickness_m: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new_in_band(
        std::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, n, k).unwrap()
        }),
        thickness_m,
        2000.0,
        authority(),
    )
    .unwrap()
}

fn temperature_film(side: usize, t: f64) -> CoherentAbsorbingCoating {
    match side {
        0 => ir_film(
            1.4 + 0.0003 * (t - 250.0),
            0.08 + 0.0007 * (t - 250.0),
            2e-6,
        ),
        _ => ir_film(
            1.9 + 0.0002 * (t - 250.0),
            0.18 + 0.0005 * (t - 250.0),
            1.2e-6,
        ),
    }
}

fn fixed_film(side: usize) -> CoherentAbsorbingCoating {
    match side {
        0 => ir_film(2.3, 0.22, 1.4e-6),
        _ => ir_film(1.2, 0.035, 3e-6),
    }
}

pub(super) fn vacuum() -> DielectricIor {
    TabulatedIor::try_new_in_band(
        std::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(4000.0 + 2000.0 * i as f64, 1.0).unwrap()
        }),
        2000.0,
    )
    .unwrap()
    .into()
}

pub(super) fn pair() -> CoatingTransmittingPair {
    CoatingTransmittingPair {
        pair_id: 17,
        surfaces: std::array::from_fn(|side| CoatingTransmittingSurface {
            source_identity: ContentHash([246 + side as u8; 32]),
            receiver: SurfaceHeatingReceiver {
                thermal_cell: CellId(11 + side as u64),
                material_id: MaterialId(7),
                temperature_bounds_k: (250.0, 450.0),
            },
            stack: IsothermalCoatingStack::Transmitting {
                films: vec![
                    CoatingOpticsModel::Law(Arc::new(
                        SurfaceCoatingLaw::try_temperature(
                            authority(),
                            &[
                                (250.0, temperature_film(side, 250.0)),
                                (450.0, temperature_film(side, 450.0)),
                            ],
                            2,
                        )
                        .unwrap(),
                    )),
                    CoatingOpticsModel::Fixed(fixed_film(side)),
                ],
                opposite_ior: vacuum(),
                reverse: REVERSE[side],
            },
        }),
        area_m2: AREA,
        baths: BATHS.map(|temperature_k| CoatingExchangeBath::Blackbody { temperature_k }),
        quadrature: CoatingExchangeQuadrature {
            wavelength_bins: 8,
            cosine_bins: 4,
        },
    }
}

// The original coherent fields supply each port's R/T/A. This oracle does not
// use CoatingThermalEmission, any exchange response, or its reduced six-edge
// coefficients. A direct two-unknown intensity balance determines gap light;
// absorbed incoming power minus escaped thermal emission gives each body heat.
// Kirchhoff applies separately to reciprocal ports and polarization modes:
// Greffet et al., Phys. Rev. X 8, 021008 (2018), equations (5)--(6).
fn absorbed_power(temperatures: [f64; 2]) -> [f64; 4] {
    let films = std::array::from_fn::<_, 2, _>(|side| {
        [temperature_film(side, temperatures[side]), fixed_film(side)]
    });
    let layers = films
        .each_ref()
        .map(|f| f.each_ref().map(CoatingLayer::Absorbing));
    let stacks = layers
        .each_ref()
        .map(|l| DielectricCoatingStack::try_mixed_in_band(l, (4000.0, 20000.0)).unwrap());
    let mut watts = [0.0; 4];
    for w in 0..8 {
        let lambda = 4000.0 + (w as f64 + 0.5) * 2000.0;
        let body = temperatures.map(|t| planck(t, lambda));
        let bath = BATHS.map(|t| planck(t, lambda));
        for angular in 0..4 {
            let mu = (angular as f64 + 0.5) / 4.0;
            let gap = std::array::from_fn::<_, 2, _>(|side| {
                stacks[side]
                    .layer_power_balance(lambda, REVERSE[side], mu, 1.0, 1.0)
                    .unwrap()
            });
            let outer = std::array::from_fn::<_, 2, _>(|side| {
                stacks[side]
                    .layer_power_balance(lambda, !REVERSE[side], mu, 1.0, 1.0)
                    .unwrap()
            });
            for polarization in 0..2 {
                let sign = if polarization == 0 { 1.0 } else { -1.0 };
                let reflection = gap.each_ref().map(|work| {
                    work.fresnel.reflection.0[0][0] + sign * work.fresnel.reflection.0[0][1]
                });
                let transmission = gap.each_ref().map(|work| {
                    work.fresnel.transmission.0[0][0] + sign * work.fresnel.transmission.0[0][1]
                });
                let exterior_reflection = outer.each_ref().map(|work| {
                    work.fresnel.reflection.0[0][0] + sign * work.fresnel.reflection.0[0][1]
                });
                let exterior_transmission = outer.each_ref().map(|work| {
                    work.fresnel.transmission.0[0][0] + sign * work.fresnel.transmission.0[0][1]
                });
                let absorption = gap.each_ref().map(|work| work.absorbed[polarization]);
                let exterior_absorption = outer.each_ref().map(|work| work.absorbed[polarization]);
                // X = T1_outer*BL + R1_gap*Y + A1_gap*B1;
                // Y = T2_outer*BR + R2_gap*X + A2_gap*B2.
                let launched = std::array::from_fn::<_, 2, _>(|side| {
                    exterior_transmission[side] * bath[side] + absorption[side] * body[side]
                });
                let rightward = (launched[0] + reflection[0] * launched[1])
                    / (1.0 - reflection[0] * reflection[1]);
                let leftward = launched[1] + reflection[1] * rightward;
                let arriving = [leftward, rightward];
                let body_heat = std::array::from_fn::<_, 2, _>(|side| {
                    exterior_absorption[side] * bath[side] + absorption[side] * arriving[side]
                        - (exterior_absorption[side] + absorption[side]) * body[side]
                });
                let bath_heat = std::array::from_fn::<_, 2, _>(|side| {
                    exterior_reflection[side] * bath[side]
                        + transmission[side] * arriving[side]
                        + exterior_absorption[side] * body[side]
                        - bath[side]
                });
                // SI Planck above includes both polarizations. Their /2 and
                // the azimuth integral 2*pi leave pi per polarization here.
                let measure = std::f64::consts::PI * AREA * mu * 2000.0 / 4.0;
                for (i, heat) in body_heat.into_iter().chain(bath_heat).enumerate() {
                    watts[i] += measure * heat;
                }
            }
        }
    }
    watts
}

#[test]
fn independently_solved_transmitting_gap_and_two_baths_refine_with_native_conduction() {
    let mut first = specimen(340.0, 0.001);
    first.boundary.source_w = 0.02;
    let mut second = specimen(280.0, 0.002);
    second.id = CellId(12);
    let native = thermal(
        vec![first, second],
        &HOT,
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 99,
                a: CellId(11),
                b: CellId(12),
                area_m2: 0.01,
                conductance_w_k: 0.001,
            }),
        }],
    );
    let owner = NativeCoatingTransmittingExchange::try_new(
        &native,
        &[pair()],
        NativeCoatingTransmittingBudget {
            max_pairs: 1,
            max_total_optical_knots: 4,
            max_source_evaluations: 10_000,
            max_optical_evaluations: 1_000_000,
            max_channel_spectral_nodes: 6 * 8,
        },
        || Ok(()),
    )
    .unwrap();
    let expected = absorbed_power([340.0, 280.0]);
    near(expected.iter().sum(), 0.0, 2e-12);
    let power = native.power_balance_observed(&owner, || Ok(())).unwrap();
    let conduction = 0.001 * (340.0 - 280.0);
    near(
        power.net_power_w[&CellId(11)],
        expected[0] + 0.02 - conduction,
        3e-12,
    );
    near(
        power.net_power_w[&CellId(12)],
        expected[1] + conduction,
        3e-12,
    );
    near(
        power.external_power_w,
        -expected[2] - expected[3] + 0.02,
        3e-12,
    );
    assert_eq!(owner.work().optical_evaluations, 4 * 8 * 4);
    assert!(owner.heat_transfers()[0].power_w > 0.0);
    assert!(owner.rates(17).unwrap().bath_transfer.outward_power_w > 0.0);

    let duration = 0.04;
    let step = duration / 128.0;
    let derivative = |x: [f64; 3]| {
        let heat = absorbed_power([x[0], x[1]]);
        let conduction = 0.001 * (x[0] - x[1]);
        [
            (heat[0] + 0.02 - conduction) / 0.01,
            (heat[1] + conduction) / 0.02,
            -heat[2] - heat[3] + 0.02,
        ]
    };
    let mut reference = [340.0, 280.0, 0.0];
    for _ in 0..128 {
        let a = derivative(reference);
        let b = derivative(std::array::from_fn(|i| reference[i] + step * 0.5 * a[i]));
        let c = derivative(std::array::from_fn(|i| reference[i] + step * 0.5 * b[i]));
        let d = derivative(std::array::from_fn(|i| reference[i] + step * c[i]));
        for i in 0..3 {
            reference[i] += step * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    for dt in [0.01, 0.005, 0.0025] {
        let (after, observed, audit) = native
            .advance_observed(&owner, duration, config(dt), || Ok(()))
            .unwrap();
        errors.push(
            (after.temperature_k(CellId(11)).unwrap() - reference[0]).abs()
                + (after.temperature_k(CellId(12)).unwrap() - reference[1]).abs()
                + (audit.external_heat_j - reference[2]).abs(),
        );
        near(
            after.energy_j() - native.energy_j(),
            audit.external_heat_j,
            2e-11,
        );
        near(
            observed.energy().thermal_external_j + 0.02 * duration,
            audit.external_heat_j,
            2e-12,
        );
        near(observed.energy().balance_residual_j, 0.0, 1e-13);
        assert!(observed.energy().by_pair_bath_transfer_j[&17] > 0.0);
        assert!(observed.work().optical_evaluations > 4 * 8 * 4);
        for side in 0..2 {
            let current = after.temperature_k(CellId(11 + side as u64)).unwrap();
            for film in 0..2 {
                assert_eq!(
                    observed.emitter(17, side).unwrap().film_temperature_k(film),
                    Some(current)
                );
            }
        }
    }
    assert!(
        errors[1] < 0.65 * errors[0] && errors[2] < 0.65 * errors[1],
        "{errors:?}"
    );
    assert_eq!(owner.time_s(), 0.0);
    assert_eq!(native.time_s(), 0.0);
}
