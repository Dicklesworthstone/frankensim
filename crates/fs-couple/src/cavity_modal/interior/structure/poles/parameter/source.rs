//! One physical derivative assembly shared by quadratic and memory poles.

use super::response::ResponseError;

/// One original physical direction per caller-declared unit of theta.
/// All blocks are complete row-major real arrays. Structural blocks have
/// ns² entries, acoustic/wall blocks na², and coupling na*ns. Zero blocks must
/// be explicit. Matrix directions are signed, not new positive material laws.
#[derive(Clone, Copy, Debug)]
pub struct StructuralInteriorPoleDirection<'d> {
    /// dM_structure/dtheta, symmetric but not necessarily positive.
    pub structural_mass: &'d [f64],
    /// dC_structure/dtheta; all viscous and gyroscopic entries are retained.
    pub structural_damping: &'d [f64],
    /// dK_structure/dtheta, symmetric in the original displacement basis.
    pub structural_stiffness: &'d [f64],
    /// dM_acoustic/dtheta in the original pressure-potential basis.
    pub acoustic_mass: &'d [f64],
    /// dK_acoustic/dtheta in the original pressure-potential basis.
    pub acoustic_stiffness: &'d [f64],
    /// Derivative of the symmetric original relative-wall admittance matrix.
    pub wall_damping: &'d [f64],
    /// dB/dtheta, pressure rows by structural columns. Geometry directions
    /// must include all trace, normal and area-measure derivatives explicitly.
    pub coupling: &'d [f64],
}

impl StructuralInteriorPoleDirection<'_> {
    pub(in crate::cavity_modal::interior::structure) fn assemble(
        &self,
        ns: usize,
        na: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Vec<f64>; 3], ResponseError> {
        poll(cancelled)?;
        for (matrix, size, symmetric) in [
            (self.structural_mass, ns, true),
            (self.structural_damping, ns, false),
            (self.structural_stiffness, ns, true),
            (self.acoustic_mass, na, true),
            (self.acoustic_stiffness, na, true),
            (self.wall_damping, na, true),
        ] {
            if size.checked_mul(size) != Some(matrix.len()) {
                return Err(ResponseError::Invalid(
                    "parameter direction must cover every original block coordinate",
                ));
            }
            for i in 0..size {
                poll(cancelled)?;
                if (0..size).any(|j| {
                    !matrix[i * size + j].is_finite()
                        || symmetric && matrix[i * size + j] != matrix[j * size + i]
                }) {
                    return Err(ResponseError::Invalid(
                        "finite directions with symmetric mass, stiffness and wall blocks required",
                    ));
                }
            }
        }
        if na.checked_mul(ns) != Some(self.coupling.len()) {
            return Err(ResponseError::Invalid(
                "one complete original dB pressure-by-structure direction required",
            ));
        }
        let n = ns.checked_add(na).ok_or(ResponseError::Invalid(
            "coupled derivative coordinate extent overflow",
        ))?;
        let entries = n.checked_mul(n).ok_or(ResponseError::Invalid(
            "coupled derivative matrix extent overflow",
        ))?;
        let mut blocks = [vec![0.0; entries], vec![0.0; entries], vec![0.0; entries]];
        for (order, offset, size, source) in [
            (0, 0, ns, self.structural_mass),
            (1, 0, ns, self.structural_damping),
            (2, 0, ns, self.structural_stiffness),
            (0, ns, na, self.acoustic_mass),
            (1, ns, na, self.wall_damping),
            (2, ns, na, self.acoustic_stiffness),
        ] {
            for i in 0..size {
                poll(cancelled)?;
                blocks[order][(offset + i) * n + offset..(offset + i) * n + offset + size]
                    .copy_from_slice(&source[i * size..(i + 1) * size]);
            }
        }
        for a in 0..na {
            poll(cancelled)?;
            for s in 0..ns {
                let value = self.coupling[a * ns + s];
                if !value.is_finite() {
                    return Err(ResponseError::Invalid(
                        "finite original interface derivative required",
                    ));
                }
                blocks[1][(ns + a) * n + s] = value;
                blocks[1][s * n + ns + a] = -value;
            }
        }
        poll(cancelled)?;
        Ok(blocks)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResponseError> {
    if cancelled() {
        Err(ResponseError::Cancelled)
    } else {
        Ok(())
    }
}
