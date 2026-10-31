//! Scalar midpoint mechanics shared by the original-material sealed-tube owner.
//! The pressure callback owns the gas; this module supplies no replacement EOS.

use core::f64::consts::PI;

#[derive(Clone, Copy)]
pub(crate) struct Input {
    pub radius: f64,
    pub length: f64,
    pub initial: f64,
    pub gas_volume: f64,
    pub viscosity: f64,
    pub density: f64,
    pub capillary_pressure: f64,
    pub inlet_pressure: f64,
    pub gravity: f64,
    pub duration: f64,
}
#[derive(Clone, Copy)]
pub(crate) struct Controls {
    pub iterations: usize,
    pub pressure_tolerance: f64,
    pub length_tolerance: f64,
}
pub(crate) struct Pressure<T> {
    pub value: f64,
    /// Derivative with respect to FINAL gas volume, with old state/heat fixed.
    pub derivative: f64,
    pub payload: T,
}
pub(crate) struct Motion<T> {
    pub length: f64,
    pub volume: f64,
    pub inlet_work: f64,
    pub wall_energy: f64,
    pub gravity_energy: f64,
    pub dissipation: f64,
    pub pressure_residual: f64,
    pub iterations: usize,
    pub pressure: Pressure<T>,
}
#[derive(Debug, PartialEq)]
pub(crate) enum Error<E> {
    Invalid(&'static str),
    Gas(E),
    Cancelled,
    Boundary,
    NonMonotone,
    NotConverged { iterations: usize, residual: f64 },
}

/// Safeguarded Newton with a geometrical bracket and bounded evaluations.
/// The admitted interval has nonnegative hydraulic-plus-gravity derivative.
/// The original stable midpoint gas law must have strictly negative dP/dV;
/// together these establish a strictly increasing mechanical residual.
pub(crate) fn solve<T, E, C: FnMut() -> bool>(
    input: Input,
    controls: Controls,
    mut gas: impl FnMut(f64, &mut C) -> Result<Pressure<T>, E>,
    cancelled: &mut C,
) -> Result<Motion<T>, Error<E>> {
    let i = input;
    if ![
        i.radius,
        i.length,
        i.gas_volume,
        i.viscosity,
        i.density,
        i.duration,
        controls.pressure_tolerance,
        controls.length_tolerance,
    ]
    .into_iter()
    .all(|x| x.is_finite() && x > 0.0)
        || ![i.initial, i.capillary_pressure, i.inlet_pressure, i.gravity]
            .into_iter()
            .all(f64::is_finite)
        || i.initial < 0.0
        || i.initial > i.length
        || i.inlet_pressure < 0.0
        || !(1..=128).contains(&controls.iterations)
    {
        return Err(Error::Invalid(
            "finite physical state, positive tolerances and 1..=128 iterations required",
        ));
    }
    let area = PI * i.radius * i.radius;
    let resistance = 8.0 * i.viscosity / (i.radius * i.radius * i.duration);
    let gravity = i.density * i.gravity;
    if ![area, resistance]
        .into_iter()
        .all(|x| x.is_finite() && x > 0.0)
        || !gravity.is_finite()
    {
        return Err(Error::Invalid(
            "capillary geometry/resistance is unrepresentable",
        ));
    }
    // Downward gravity is admitted only where its negative derivative is
    // offset by the liquid resistance alone; gas stiffness then adds margin.
    let mut lower = if gravity < 0.0 {
        (-0.5 * gravity / resistance).next_up()
    } else {
        0.0
    };
    let mut upper = i.length;
    if lower > i.initial || lower >= upper || resistance * lower + 0.5 * gravity < 0.0 {
        return Err(Error::NonMonotone);
    }
    let mut x = i.initial;
    let mut last_residual = 0.0;
    for iteration in 0..controls.iterations {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        let displacement = x - i.initial;
        let volume = area * displacement;
        if displacement != 0.0 && volume == 0.0 {
            return Err(Error::Invalid("displaced volume is unrepresentable"));
        }
        let gas_volume = i.gas_volume - volume;
        if !gas_volume.is_finite() || gas_volume <= 0.0 {
            return Err(Error::Boundary);
        }
        if volume != 0.0 && gas_volume == i.gas_volume {
            return Err(Error::Invalid(
                "gas volume cannot represent the meniscus displacement",
            ));
        }
        let pressure = gas(gas_volume, cancelled).map_err(Error::Gas)?;
        if !pressure.value.is_finite()
            || pressure.value <= 0.0
            || !pressure.derivative.is_finite()
            || pressure.derivative >= 0.0
        {
            return Err(Error::NonMonotone);
        }
        let midpoint = i.initial + 0.5 * displacement;
        let viscous_pressure = resistance * midpoint * displacement;
        let residual = viscous_pressure
            - ((i.inlet_pressure - pressure.value) + i.capillary_pressure - gravity * midpoint);
        let tangent = resistance * x + 0.5 * gravity - area * pressure.derivative;
        if !residual.is_finite() || !tangent.is_finite() || tangent <= 0.0 {
            return Err(Error::Invalid(
                "mechanical residual or positive tangent is unrepresentable",
            ));
        }
        last_residual = residual;
        let correction = residual / tangent;
        if residual.abs() <= controls.pressure_tolerance
            && (residual == 0.0 || correction.abs() <= controls.length_tolerance)
        {
            let inlet_work = i.inlet_pressure * volume;
            let wall_energy = -i.capillary_pressure * volume;
            let gravity_energy = gravity * midpoint * volume;
            let dissipation = viscous_pressure * volume;
            if ![volume, inlet_work, wall_energy, gravity_energy, dissipation]
                .into_iter()
                .all(f64::is_finite)
                || dissipation < 0.0
            {
                return Err(Error::Invalid("mechanical work is unrepresentable"));
            }
            if cancelled() {
                return Err(Error::Cancelled);
            }
            return Ok(Motion {
                length: x,
                volume,
                inlet_work,
                wall_energy,
                gravity_energy,
                dissipation,
                pressure_residual: residual,
                iterations: iteration + 1,
                pressure,
            });
        }
        if residual < 0.0 {
            lower = x;
        } else {
            upper = x;
        }
        if lower == upper {
            return Err(Error::Boundary);
        }
        let newton = x - correction;
        let next = if newton.is_finite() && newton > lower && newton < upper {
            newton
        } else {
            lower + 0.5 * (upper - lower)
        };
        if next == x {
            return Err(Error::NotConverged {
                iterations: iteration + 1,
                residual,
            });
        }
        x = next;
    }
    Err(Error::NotConverged {
        iterations: controls.iterations,
        residual: last_residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn input() -> Input {
        Input {
            radius: 1e-4,
            length: 0.1,
            initial: 0.02,
            gas_volume: PI * 1e-8 * 0.08,
            viscosity: 1.0,
            density: 1000.0,
            capillary_pressure: 300.0,
            inlet_pressure: 100000.0,
            gravity: 0.0,
            duration: 1.0,
        }
    }
    fn controls() -> Controls {
        Controls {
            iterations: 80,
            pressure_tolerance: 1e-7,
            length_tolerance: 1e-13,
        }
    }
    // Independent exact calorically-perfect implicit-midpoint oracle.
    fn pressure(i: Input, heat: f64, v: f64) -> Pressure<()> {
        let numerator = 2.0 * 100000.0 * i.gas_volume + 0.4 * heat;
        let denominator = 1.4 * v + 0.6 * i.gas_volume;
        Pressure {
            value: numerator / denominator,
            derivative: -1.4 * numerator / denominator.powi(2),
            payload: (),
        }
    }
    #[test]
    fn manufactured_compression_and_drainage_close_the_same_discrete_work() {
        for displacement in [0.001, -0.001] {
            for heat in [0.0, 1e-5] {
                let mut i = input();
                let x = i.initial + displacement;
                let midpoint = 0.5 * (i.initial + x);
                let dv = PI * i.radius.powi(2) * displacement;
                let p = pressure(i, heat, i.gas_volume - dv).value;
                i.inlet_pressure = p - i.capillary_pressure
                    + 8.0 * i.viscosity * midpoint * displacement / (i.radius.powi(2) * i.duration);
                let step = solve(
                    i,
                    controls(),
                    |v, _| Ok::<_, ()>(pressure(i, heat, v)),
                    &mut || false,
                )
                .unwrap();
                assert!((step.length - x).abs() < 1e-13);
                let gas_storage = heat + step.pressure.value * step.volume;
                let residual = step.inlet_work + heat
                    - gas_storage
                    - step.wall_energy
                    - step.gravity_energy
                    - step.dissipation;
                assert!(residual.abs() < 1e-15, "{residual}");
                assert!(step.dissipation > 0.0);
            }
        }
    }
    #[test]
    fn adiabatic_motion_converges_to_continuous_gas_compression_at_second_order() {
        let base = input();
        let target = 0.025;
        let inlet = 125000.0;
        let area = PI * base.radius.powi(2);
        let rate = base.radius.powi(2) / (8.0 * base.viscosity);
        let time_integrand = |x: f64| {
            let volume = base.gas_volume - area * (x - base.initial);
            let pressure = 100000.0 * (base.gas_volume / volume).powf(1.4);
            x / (rate * (inlet + base.capillary_pressure - pressure))
        };
        // Independent continuous solution by finely resolved travel-time integral.
        let n = 8192;
        let h = (target - base.initial) / n as f64;
        let mut sum = time_integrand(base.initial) + time_integrand(target);
        for k in 1..n {
            sum +=
                (if k % 2 == 0 { 2.0 } else { 4.0 }) * time_integrand(base.initial + k as f64 * h);
        }
        let time = h * sum / 3.0;
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let mut i = Input {
                inlet_pressure: inlet,
                duration: time / steps as f64,
                ..base
            };
            let mut initial_pressure = 100000.0;
            for _ in 0..steps {
                let v0 = i.gas_volume;
                let p0 = initial_pressure;
                let step = solve(
                    i,
                    controls(),
                    |v, _| {
                        let denominator = 1.4 * v + 0.6 * v0;
                        Ok::<_, ()>(Pressure {
                            value: 2.0 * p0 * v0 / denominator,
                            derivative: -2.8 * p0 * v0 / denominator.powi(2),
                            payload: (),
                        })
                    },
                    &mut || false,
                )
                .unwrap();
                let v1 = v0 - step.volume;
                initial_pressure = (step.pressure.value * (v0 + v1) - p0 * v0) / v1;
                i.initial = step.length;
                i.gas_volume = v1;
            }
            errors.push((i.initial - target).abs());
        }
        assert!(errors[2] < 1e-7, "{errors:?}");
        for pair in errors.windows(2) {
            assert!(
                pair[0] / pair[1] > 3.8 && pair[0] / pair[1] < 4.2,
                "{errors:?}"
            );
        }
    }
    #[test]
    fn cancellation_source_refusal_and_nonmonotone_gravity_are_explicit() {
        assert!(matches!(
            solve(
                input(),
                controls(),
                |v, _| Ok::<_, ()>(pressure(input(), 0.0, v)),
                &mut || true
            ),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            solve(
                input(),
                controls(),
                |_, _| Err::<Pressure<()>, _>("source"),
                &mut || false
            ),
            Err(Error::Gas("source"))
        ));
        let mut i = input();
        i.gravity = -1e9;
        assert!(matches!(
            solve(
                i,
                controls(),
                |v, _| Ok::<_, ()>(pressure(i, 0.0, v)),
                &mut || false
            ),
            Err(Error::NonMonotone)
        ));
    }
}
