//! Shared source-profile integration. Only preparation performs this work.
use super::*;

impl AtmosphereQuadrature {
    pub(super) fn admit(self) -> Result<(), Error> {
        if !positive(self.relative_tolerance)
            || self.relative_tolerance > 1e-2
            || !(2..=128).contains(&self.maximum_panels)
        {
            return Err(Error::Invalid("atmosphere quadrature controls"));
        }
        Ok(())
    }
}

pub(super) fn average<C, F>(
    gas: &FrozenGas,
    left: f64,
    right: f64,
    control: AtmosphereQuadrature,
    point: &mut F,
    cancelled: &mut C,
) -> Result<Conserved, Error>
where
    C: FnMut() -> bool,
    F: FnMut(f64, &mut C) -> Result<(Conserved, Primitive), Error>,
{
    let mut panels = 1;
    let mut coarse = integrate(left, right, panels, point, cancelled)?;
    while panels < control.maximum_panels {
        panels = (2 * panels).min(control.maximum_panels);
        let fine = integrate(left, right, panels, point, cancelled)?;
        let density_tolerance = control.relative_tolerance * fine.density.abs();
        let energy_tolerance = control.relative_tolerance * fine.energy.abs()
            + fine.density * gas.domain().absolute_energy_tolerance;
        if (fine.density - coarse.density).abs() <= density_tolerance
            && (fine.energy - coarse.energy).abs() <= energy_tolerance
        {
            return Ok(fine);
        }
        coarse = fine;
    }
    Err(Error::Invalid("atmosphere quadrature budget exhausted"))
}

fn integrate<C, F>(
    left: f64,
    right: f64,
    panels: usize,
    point: &mut F,
    cancelled: &mut C,
) -> Result<Conserved, Error>
where
    C: FnMut() -> bool,
    F: FnMut(f64, &mut C) -> Result<(Conserved, Primitive), Error>,
{
    const NODES: [f64; 4] = [
        -0.861_136_311_594_052_6,
        -0.339_981_043_584_856_3,
        0.339_981_043_584_856_3,
        0.861_136_311_594_052_6,
    ];
    const WEIGHTS: [f64; 4] = [
        0.347_854_845_137_453_9,
        0.652_145_154_862_546_1,
        0.652_145_154_862_546_1,
        0.347_854_845_137_453_9,
    ];
    let mut sums = [0.0; 2];
    let mut corrections = [0.0; 2];
    let width = (right - left) / panels as f64;
    for panel in 0..panels {
        let middle = left + (panel as f64 + 0.5) * width;
        for (node, weight) in NODES.into_iter().zip(WEIGHTS) {
            let (q, _) = point(middle + 0.5 * width * node, cancelled)?;
            let factor = weight / (2 * panels) as f64;
            for (field, value) in [q.density, q.energy].into_iter().enumerate() {
                let y = factor * value - corrections[field];
                let next = sums[field] + y;
                corrections[field] = (next - sums[field]) - y;
                sums[field] = next;
            }
        }
    }
    Ok(Conserved {
        density: sums[0],
        momentum: 0.0,
        energy: sums[1],
    })
}
