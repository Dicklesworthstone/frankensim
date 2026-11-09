//! Independent derivatives of original displacement, momentum and memory rows.

use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn check(
    p: &Prepared<'_>,
    direction: &PreparedDirection<'_, '_>,
    omega: f64,
    state: &State,
    base: &State,
    input: &[C64],
    frequency: Option<(&State, &State)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let n = p.source.coordinates;
    let [dm, dc, dk] = &direction.blocks;
    let mut relative = 0.0_f64;
    for i in 0..n {
        poll(cancelled)?;
        let mut row = Row::new();
        row.add(state.q[i], C64::new(0.0, -omega));
        row.add(state.velocity[i], C64::from_re(-1.0));
        if let Some((parameter, _)) = frequency {
            row.add(parameter.q[i], C64::new(0.0, -1.0));
        }
        relative = relative.max(row.relative()?);
        let mut row = Row::new();
        row.add(input[i], C64::from_re(-1.0));
        for j in 0..n {
            let ij = i * n + j;
            row.add(state.q[j], C64::from_re(p.source.stiffness[ij]));
            row.add(
                state.velocity[j],
                C64::new(p.source.damping[ij], -omega * p.source.mass[ij]),
            );
            row.add(base.q[j], C64::from_re(dk[ij]));
            row.add(base.velocity[j], C64::new(dc[ij], -omega * dm[ij]));
            if let Some((parameter, original)) = frequency {
                // d/domega of s*M*v_theta and s*M_theta*v.
                row.add(parameter.velocity[j], C64::new(0.0, -p.source.mass[ij]));
                row.add(original.velocity[j], C64::new(0.0, -dm[ij]));
            }
        }
        for (j, (branch, derivative)) in p.branches.iter().zip(&direction.branches).enumerate() {
            row.add(state.stress[j], C64::from_re(branch.strain[i]));
            row.add(base.stress[j], C64::from_re(derivative.projection[i]));
        }
        relative = relative.max(row.relative()?);
    }
    for (j, (branch, derivative)) in p.branches.iter().zip(&direction.branches).enumerate() {
        poll(cancelled)?;
        let mut row = Row::new();
        row.add(
            state.stress[j],
            C64::new(branch.relaxation_time_s.recip(), -omega),
        );
        row.add(
            base.stress[j],
            C64::from_re(direction.native.branch_rate_derivatives()[j]),
        );
        for i in 0..n {
            row.add(
                state.velocity[i],
                C64::from_re(-branch.stiffness * branch.strain[i]),
            );
            row.add(
                base.velocity[i],
                C64::from_re(-derivative.stiffness * branch.strain[i]),
            );
            row.add(
                base.velocity[i],
                C64::from_re(-branch.stiffness * derivative.projection[i]),
            );
        }
        if let Some((parameter, _)) = frequency {
            row.add(parameter.stress[j], C64::new(0.0, -1.0));
        }
        relative = relative.max(row.relative()?);
    }
    let tolerance = p.limits.physical_relative_residual_tolerance;
    if relative > tolerance {
        return Err(Error::PhysicalResidual {
            relative,
            tolerance,
        });
    }
    poll(cancelled)?;
    Ok(relative)
}
