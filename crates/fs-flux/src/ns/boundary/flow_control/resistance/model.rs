//! Finite reservoir laws on an arbitrary cached native port mobility.
//!
//! Mobility may have pressure-shift nullspaces or be identically zero. Only
//! the actual coupled pressure system is factored; no datum is manufactured.

use fs_la::factor::lu;
use fs_math::dd::Dd;

/// One external reservoir and linear resistance attached to a prepared port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesResistancePort {
    /// Reservoir pressure increment relative to the original port traction,
    /// divided by the SAME reference density [m²/s²]. Signed and finite.
    pub reservoir_pressure_m2_s2: f64,
    /// Nonnegative external resistance divided by reference density [1/s].
    /// The port law is p = p_reservoir + resistance * Q_outward, where Q is
    /// volume flow per unit depth [m²/s]. Zero prescribes reservoir pressure.
    pub resistance_s_inv: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum ResistanceResponseError {
    Invalid(&'static str),
    Singular,
}

/// Solve Q = Q0 + M*p together with p = p_reservoir + R*Q.
/// The caller still checks the original continuum and integrated-flow laws.
pub(super) fn solve_pressure_law(
    mobility: &[f64],
    baseline_flows: &[f64],
    loads: &[StokesResistancePort],
) -> Result<Vec<f64>, ResistanceResponseError> {
    let count = loads.len();
    let entries = count.checked_mul(count).ok_or(ResistanceResponseError::Invalid(
        "port-response dimensions overflow",
    ))?;
    if count == 0
        || baseline_flows.len() != count
        || mobility.len() != entries
        || mobility.iter().chain(baseline_flows).any(|x| !x.is_finite())
        || loads.iter().any(|load| {
            !load.reservoir_pressure_m2_s2.is_finite()
                || !load.resistance_s_inv.is_finite()
                || load.resistance_s_inv < 0.0
        })
    {
        return Err(ResistanceResponseError::Invalid(
            "finite complete native response and passive reservoir loads required",
        ));
    }
    let mut matrix = vec![0.0; mobility.len()];
    let mut rhs = vec![0.0; count];
    let mut row_scales = vec![0.0_f64; count];
    for (i, load) in loads.iter().enumerate() {
        let resistance = Dd::from_f64(load.resistance_s_inv);
        rhs[i] = (Dd::from_f64(load.reservoir_pressure_m2_s2)
            + resistance * Dd::from_f64(baseline_flows[i]))
        .to_f64();
        for j in 0..count {
            let entry = (Dd::from_f64(if i == j { 1.0 } else { 0.0 })
                - resistance * Dd::from_f64(mobility[i * count + j]))
            .to_f64();
            matrix[i * count + j] = entry;
            row_scales[i] = row_scales[i].max(entry.abs());
        }
    }
    if matrix.iter().chain(&rhs).any(|x| !x.is_finite()) {
        return Err(ResistanceResponseError::Invalid(
            "unrepresentable external hydraulic response",
        ));
    }
    if row_scales.iter().any(|&scale| scale == 0.0) {
        return Err(ResistanceResponseError::Singular);
    }
    let mut scaled = matrix.clone();
    let mut pressures = rhs.clone();
    for i in 0..count {
        pressures[i] /= row_scales[i];
        for j in 0..count {
            scaled[i * count + j] /= row_scales[i];
        }
    }
    let factor = lu(&scaled, count).map_err(|_| ResistanceResponseError::Singular)?;
    factor.solve(&mut pressures);
    for _ in 0..2 {
        let mut correction = Vec::with_capacity(count);
        for i in 0..count {
            let mut residual = Dd::from_f64(rhs[i]);
            for (j, &pressure) in pressures.iter().enumerate() {
                residual =
                    residual - Dd::from_f64(matrix[i * count + j]) * Dd::from_f64(pressure);
            }
            correction.push(residual.to_f64() / row_scales[i]);
        }
        factor.solve(&mut correction);
        for (pressure, delta) in pressures.iter_mut().zip(correction) {
            *pressure += delta;
        }
    }
    if pressures.iter().any(|p| !p.is_finite()) {
        return Err(ResistanceResponseError::Invalid(
            "unrepresentable external port pressure",
        ));
    }
    Ok(pressures)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-9 * expected.abs().max(1.0),
            "{actual} versus {expected}",
        );
    }

    fn flow(mobility: &[f64], baseline: &[f64], pressure: &[f64]) -> Vec<f64> {
        let count = baseline.len();
        (0..count)
            .map(|i| {
                pressure.iter().enumerate().fold(
                    Dd::from_f64(baseline[i]),
                    |sum, (j, &p)| {
                        sum + Dd::from_f64(mobility[i * count + j]) * Dd::from_f64(p)
                    },
                ).to_f64()
            })
            .collect()
    }

    #[test]
    fn singular_mobility_retains_common_pressure_and_series_resistances() {
        let conductance = 0.25;
        let mobility = [-conductance, conductance, conductance, -conductance];
        let baseline = [-0.3, 0.3];
        for (left, right) in [(0.0, 0.0), (0.0, 7.0), (3.0, 7.0)] {
            for drive in [-2.0, 0.0, 2.0] {
                let q = (0.3 + conductance * drive)
                    / (1.0 + conductance * (left + right));
                for shift in [0.0, 17.0, 100_000.0] {
                    let loads = [
                        StokesResistancePort {
                            reservoir_pressure_m2_s2: shift + drive,
                            resistance_s_inv: left,
                        },
                        StokesResistancePort {
                            reservoir_pressure_m2_s2: shift,
                            resistance_s_inv: right,
                        },
                    ];
                    let p = solve_pressure_law(&mobility, &baseline, &loads).unwrap();
                    let actual = flow(&mobility, &baseline, &p);
                    near(actual[0], -q);
                    near(actual[1], q);
                    near(actual.iter().sum(), 0.0);
                    near(p[0], shift + drive - left * q);
                    near(p[1], shift + right * q);
                    let external_loss = left * actual[0].powi(2)
                        + right * actual[1].powi(2);
                    let fluid_power = (p[0] - p[1]) * actual[1];
                    near(drive * actual[1], external_loss + fluid_power);
                    near(fluid_power, q * (q - 0.3) / conductance);
                }
            }
        }
    }

    #[test]
    fn zero_mobility_needs_no_dummy_datum_or_flow_factorization() {
        for q in [0.0, -0.4, 0.4] {
            for resistance in [0.0, 3.0, 10_000.0] {
                let p = solve_pressure_law(
                    &[0.0],
                    &[q],
                    &[StokesResistancePort {
                        reservoir_pressure_m2_s2: 2.0,
                        resistance_s_inv: resistance,
                    }],
                ).unwrap();
                near(p[0], 2.0 + resistance * q);
            }
        }
    }

    #[test]
    fn disconnected_components_keep_independent_pressure_shifts_and_conservation() {
        let mobility = [
            -0.5, 0.5, 0.0, 0.0,
            0.5, -0.5, 0.0, 0.0,
            0.0, 0.0, -0.25, 0.25,
            0.0, 0.0, 0.25, -0.25,
        ];
        let baseline = [-0.3, 0.3, 0.1, -0.1];
        let mut loads = [
            StokesResistancePort { reservoir_pressure_m2_s2: 2.0, resistance_s_inv: 3.0 },
            StokesResistancePort { reservoir_pressure_m2_s2: -1.0, resistance_s_inv: 7.0 },
            StokesResistancePort { reservoir_pressure_m2_s2: 1.0, resistance_s_inv: 0.0 },
            StokesResistancePort { reservoir_pressure_m2_s2: -2.0, resistance_s_inv: 4.0 },
        ];
        let p = solve_pressure_law(&mobility, &baseline, &loads).unwrap();
        let q = flow(&mobility, &baseline, &p);
        near(q[0] + q[1], 0.0);
        near(q[2] + q[3], 0.0);
        near(q[1], (0.3 + 0.5 * 3.0) / (1.0 + 0.5 * 10.0));
        near(q[3], (-0.1 + 0.25 * 3.0) / (1.0 + 0.25 * 4.0));
        for load in &mut loads[..2] {
            load.reservoir_pressure_m2_s2 += 19.0;
        }
        let shifted = solve_pressure_law(&mobility, &baseline, &loads).unwrap();
        let shifted_q = flow(&mobility, &baseline, &shifted);
        for i in 0..4 {
            near(shifted_q[i], q[i]);
            near(shifted[i] - p[i], if i < 2 { 19.0 } else { 0.0 });
        }
    }

    #[test]
    fn only_the_actual_finite_load_system_must_be_nonsingular() {
        let mut load = StokesResistancePort {
            reservoir_pressure_m2_s2: 2.0,
            resistance_s_inv: 1.0,
        };
        assert_eq!(
            solve_pressure_law(&[1.0], &[0.0], &[load]),
            Err(ResistanceResponseError::Singular),
        );
        for resistance in [-1.0, f64::NAN, f64::INFINITY] {
            load.resistance_s_inv = resistance;
            assert!(matches!(
                solve_pressure_law(&[0.0], &[0.0], &[load]),
                Err(ResistanceResponseError::Invalid(_)),
            ));
        }
        load.resistance_s_inv = 1e308;
        assert!(solve_pressure_law(&[-1e308], &[0.0], &[load]).is_err());
        load.resistance_s_inv = 0.0;
        assert!(solve_pressure_law(&[f64::NAN], &[0.0], &[load]).is_err());
        assert!(solve_pressure_law(&[], &[0.0], &[load]).is_err());
        assert!(solve_pressure_law(&[0.0], &[], &[load]).is_err());
        assert!(solve_pressure_law(&[], &[], &[]).is_err());
    }
}
