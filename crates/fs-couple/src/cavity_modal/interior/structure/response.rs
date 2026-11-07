//! Original-coordinate harmonic storage and reciprocal interface work.
//!
//! The acoustic coordinate is pressure potential, whose velocity is pressure.
//! These contractions retain the assembled inverse-density/compliance metrics;
//! they never infer energy or heat from a response peak or a work residual.

use fs_math::{c64::C64, dd::Dd};

/// Cycle-averaged storage of the original structural and fluid fields [J].
/// All inputs are peak phasors. The physical quadratic energy therefore has
/// factor 1/4. Small signed contraction roundoff is retained, not clipped.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StructuralInteriorEnergy {
    pub structural_kinetic_j: f64,
    pub structural_strain_j: f64,
    pub fluid_kinetic_j: f64,
    pub fluid_compression_j: f64,
}

impl StructuralInteriorEnergy {
    #[must_use]
    pub fn structural_j(self) -> f64 {
        self.structural_kinetic_j + self.structural_strain_j
    }

    #[must_use]
    pub fn acoustic_j(self) -> f64 {
        self.fluid_kinetic_j + self.fluid_compression_j
    }

    #[must_use]
    pub fn total_j(self) -> f64 {
        self.structural_j() + self.acoustic_j()
    }
}

#[derive(Debug, PartialEq)]
pub(super) enum ResponseError {
    Invalid(&'static str),
    Cancelled,
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResponseError> {
    if cancelled() {
        Err(ResponseError::Cancelled)
    } else {
        Ok(())
    }
}

/// The same original B supplies both force B^T p and displaced flow B v.
/// Double-double accumulation preserves weak signed interface components.
pub(super) fn interface(
    coupling: &[f64],
    velocity: &[C64],
    pressure: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, Vec<C64>), ResponseError> {
    poll(cancelled)?;
    let ns = velocity.len();
    let na = pressure.len();
    if ns.checked_mul(na) != Some(coupling.len())
        || coupling.iter().any(|x| !x.is_finite())
        || velocity.iter().chain(pressure).any(|&z| !finite(z))
    {
        return Err(ResponseError::Invalid(
            "finite original interface coordinates required",
        ));
    }
    let mut force_re = vec![Dd::from_f64(0.0); ns];
    let mut force_im = vec![Dd::from_f64(0.0); ns];
    let mut flow = Vec::with_capacity(na);
    for a in 0..na {
        poll(cancelled)?;
        let (mut re, mut im) = (Dd::from_f64(0.0), Dd::from_f64(0.0));
        for s in 0..ns {
            let b = Dd::from_f64(coupling[a * ns + s]);
            force_re[s] = force_re[s] + b * Dd::from_f64(pressure[a].re);
            force_im[s] = force_im[s] + b * Dd::from_f64(pressure[a].im);
            re = re + b * Dd::from_f64(velocity[s].re);
            im = im + b * Dd::from_f64(velocity[s].im);
        }
        flow.push(C64::new(re.to_f64(), im.to_f64()));
    }
    let force: Vec<_> = force_re
        .into_iter()
        .zip(force_im)
        .map(|(re, im)| C64::new(re.to_f64(), im.to_f64()))
        .collect();
    if force.iter().chain(&flow).any(|&z| !finite(z)) {
        return Err(ResponseError::Invalid(
            "interface force or flow is unrepresentable",
        ));
    }
    poll(cancelled)?;
    Ok((force, flow))
}

/// Original peak effort/flow work Re(e^H f)/2; never a transfer-row projection.
pub(super) fn work(effort: &[C64], flow: &[C64]) -> f64 {
    if effort.len() != flow.len() {
        return f64::NAN;
    }
    let mut value = Dd::from_f64(0.0);
    for (e, f) in effort.iter().zip(flow) {
        value = value
            + Dd::from_f64(e.re) * Dd::from_f64(f.re)
            + Dd::from_f64(e.im) * Dd::from_f64(f.im);
    }
    (Dd::from_f64(0.5) * value).to_f64()
}

/// Both physical chain-rule terms, including an authored source derivative.
pub(super) fn work_derivative(
    effort: &[C64],
    effort_derivative: Option<&[C64]>,
    flow: &[C64],
    flow_derivative: &[C64],
) -> f64 {
    if effort.len() != flow.len()
        || flow.len() != flow_derivative.len()
        || effort_derivative.is_some_and(|d| d.len() != effort.len())
    {
        return f64::NAN;
    }
    let mut value = Dd::from_f64(0.0);
    for i in 0..flow.len() {
        let (e, f, df) = (effort[i], flow[i], flow_derivative[i]);
        let de = effort_derivative.map_or(C64::ZERO, |d| d[i]);
        value = value
            + Dd::from_f64(de.re) * Dd::from_f64(f.re)
            + Dd::from_f64(de.im) * Dd::from_f64(f.im)
            + Dd::from_f64(e.re) * Dd::from_f64(df.re)
            + Dd::from_f64(e.im) * Dd::from_f64(df.im);
    }
    (Dd::from_f64(0.5) * value).to_f64()
}

#[allow(clippy::too_many_arguments)]
fn block_form(
    matrix: &[f64],
    n: usize,
    offset: usize,
    left: &[C64],
    right: &[C64],
    symmetric_part: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, f64), ResponseError> {
    if n.checked_mul(n) != Some(matrix.len())
        || left.len() != right.len()
        || offset.checked_add(left.len()).is_none_or(|end| end > n)
    {
        return Err(ResponseError::Invalid(
            "original physical block extent required",
        ));
    }
    let mut value = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for (i, &x) in left.iter().enumerate() {
        poll(cancelled)?;
        for (j, &y) in right.iter().enumerate() {
            let a = matrix[(offset + i) * n + offset + j];
            let a = if symmetric_part {
                0.5 * a + 0.5 * matrix[(offset + j) * n + offset + i]
            } else {
                a
            };
            value = value
                + Dd::from_f64(a)
                    * (Dd::from_f64(x.re) * Dd::from_f64(y.re)
                        + Dd::from_f64(x.im) * Dd::from_f64(y.im));
            absolute += a.abs() * x.abs() * y.abs();
        }
    }
    let value = value.to_f64();
    if !value.is_finite() || !absolute.is_finite() {
        return Err(ResponseError::Invalid(
            "physical storage or work is unrepresentable",
        ));
    }
    Ok((value, absolute))
}

fn storage(
    matrix: &[f64],
    n: usize,
    offset: usize,
    field: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ResponseError> {
    let (value, absolute) = block_form(matrix, n, offset, field, field, false, cancelled)?;
    if value < -64.0 * f64::EPSILON * absolute {
        return Err(ResponseError::Invalid(
            "negative original physical storage form",
        ));
    }
    Ok(0.25 * value)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn energy(
    mass: &[f64],
    stiffness: &[f64],
    structural_displacement: &[C64],
    structural_velocity: &[C64],
    pressure_potential: &[C64],
    pressure: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StructuralInteriorEnergy, ResponseError> {
    poll(cancelled)?;
    let ns = structural_displacement.len();
    let n = ns
        .checked_add(pressure.len())
        .ok_or(ResponseError::Invalid("coupled physical extent overflow"))?;
    if structural_velocity.len() != ns || pressure_potential.len() != pressure.len() {
        return Err(ResponseError::Invalid(
            "complete original physical fields required",
        ));
    }
    let result = StructuralInteriorEnergy {
        structural_kinetic_j: storage(mass, n, 0, structural_velocity, cancelled)?,
        structural_strain_j: storage(stiffness, n, 0, structural_displacement, cancelled)?,
        fluid_kinetic_j: storage(stiffness, n, ns, pressure_potential, cancelled)?,
        fluid_compression_j: storage(mass, n, ns, pressure, cancelled)?,
    };
    if !result.total_j().is_finite() {
        return Err(ResponseError::Invalid(
            "total original physical storage is unrepresentable",
        ));
    }
    poll(cancelled)?;
    Ok(result)
}

pub(super) fn dissipation(
    damping: &[f64],
    n: usize,
    offset: usize,
    velocity: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ResponseError> {
    Ok(0.5 * block_form(damping, n, offset, velocity, velocity, true, cancelled)?.0)
}

pub(super) fn dissipation_derivative(
    damping: &[f64],
    n: usize,
    offset: usize,
    velocity: &[C64],
    derivative: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ResponseError> {
    // d(v^H D v/2) = Re(v_w^H D v), with D=(C+C^T)/2 fixed.
    Ok(block_form(damping, n, offset, derivative, velocity, true, cancelled)?.0)
}

#[cfg(test)]
#[path = "response/tests.rs"]
mod tests;
