//! Harmonic motion of supported three-dimensional elastic assemblies.
//!
//! Prescribed coordinates must still be present in the supplied assembly.
//! Assemble with those coordinates unconstrained, then identify them by their
//! original global DOF labels here. This retains consistent cross inertia and
//! viscosity, which cannot be recovered after ordinary support elimination.

use super::TetElasticAssembly;
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseControls, HarmonicResponseError, PreparedSupportedHarmonicResponse,
    SupportedHarmonicResponse,
};
use fs_sparse::Csr;

/// One frequency factor on the complete represented solid coordinate space.
/// P1/P2 tetrahedra and hexahedra use the same native assembly carrier.
pub struct SupportedSolidHarmonicOperator<'a> {
    source: &'a TetElasticAssembly,
    free_dofs: Vec<usize>,
    support_dofs: Vec<usize>,
    response: PreparedSupportedHarmonicResponse,
}

impl SupportedSolidHarmonicOperator<'_> {
    /// Original immutable assembly, including its material/geometry metadata.
    pub fn source(&self) -> &TetElasticAssembly {
        self.source
    }

    /// Original global labels of the coordinates solved by force balance.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Original global labels of the prescribed coordinates, in source order.
    pub fn support_dofs(&self) -> &[usize] {
        &self.support_dofs
    }

    /// Original global labels of every returned displacement/velocity entry.
    pub fn represented_dofs(&self) -> &[usize] {
        &self.source.free_dofs
    }

    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.response.angular_frequency_rad_s()
    }

    /// Solve with peak phasors under `exp(-i omega t)`.
    ///
    /// Applied forces [N] follow `represented_dofs()`, including supported rows;
    /// prescribed displacements [m] follow `support_dofs()`. Returned reactions
    /// [N] are the forces BY the supports ON the solid, in support order. The
    /// full response retains original coordinates and independently evaluates
    /// applied power, support power and full viscous dissipation [W].
    ///
    /// This is a particular steady harmonic response of the original linear
    /// pencil, not transient evolution or a modal truncation. Coordinates
    /// already absent from the assembly remain stationary; their reactions and
    /// their moving-support coupling cannot be inferred. No heat is deposited.
    pub fn solve(
        &self,
        applied_force_n: &[C64],
        support_displacement_m: &[C64],
        cx: &Cx<'_>,
    ) -> Result<SupportedHarmonicResponse, HarmonicResponseError> {
        self.response
            .solve(applied_force_n, support_displacement_m, cx)
    }
}

impl TetElasticAssembly {
    /// Prepare prescribed harmonic support motion without discarding M_fc/C_fc.
    ///
    /// `support_dofs` are strictly increasing ORIGINAL GLOBAL labels still in
    /// `self.free_dofs`, not reduced row numbers. Empty supports admit ordinary
    /// forcing; all represented coordinates may also be prescribed. To retain
    /// every original support reaction, assemble the source without eliminating
    /// its supports first. A coordinate absent from this assembly is refused.
    ///
    /// Optional `damping_ns_m` is the full physical viscosity in the assembly's
    /// represented order, including support cross terms. None means explicitly
    /// zero viscosity. Source M/K and the symmetric part of C undergo the shared
    /// physical admission; no modal loss, mass lumping or artificial damping is
    /// introduced. The skew part of C is retained as nondissipative gyroscopy.
    /// Native stiffness skew within the diagonal-scaled assembly roundoff bound
    /// is symmetrized in conversion scratch only; larger skew is refused. The
    /// original CSR matrices are never rewritten.
    ///
    /// This bounded dense factor is for complete small/medium assemblies. The
    /// control envelope covers CSR-to-dense scratch as well as the shared
    /// factor and output. Cancellation publishes no partial factor or response.
    pub fn prepare_supported_harmonic_response(
        &self,
        support_dofs: &[usize],
        damping_ns_m: Option<&Csr>,
        angular_frequency_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<SupportedSolidHarmonicOperator<'_>, HarmonicResponseError> {
        cx.checkpoint()
            .map_err(|_| HarmonicResponseError::Cancelled)?;
        let n = self.free_dofs.len();
        if n == 0 || self.free_dofs.windows(2).any(|w| w[0] >= w[1]) {
            return Err(HarmonicResponseError::Invalid(
                "solid represented DOFs must be nonempty and strictly increasing",
            ));
        }
        if support_dofs.len() > n || support_dofs.windows(2).any(|w| w[0] >= w[1]) {
            return Err(HarmonicResponseError::Invalid(
                "solid support labels must be strictly increasing and represented",
            ));
        }
        if n > controls.maximum_coordinates {
            return Err(HarmonicResponseError::Budget {
                quantity: "supported solid coordinates",
                requested: n,
                cap: controls.maximum_coordinates,
            });
        }
        let scratch_bytes = n
            .checked_mul(n)
            .and_then(|square| square.checked_mul(24))
            .and_then(|bytes| n.checked_mul(32).and_then(|maps| bytes.checked_add(maps)))
            .ok_or(HarmonicResponseError::Invalid(
                "supported solid harmonic payload overflow",
            ))?;
        let payload = controls
            .maximum_payload_bytes
            .checked_sub(scratch_bytes)
            .ok_or(HarmonicResponseError::Budget {
                quantity: "supported solid conversion bytes",
                requested: scratch_bytes,
                cap: controls.maximum_payload_bytes,
            })?;
        let mut support_rows = Vec::with_capacity(support_dofs.len());
        for &label in support_dofs {
            support_rows.push(self.free_dofs.binary_search(&label).map_err(|_| {
                HarmonicResponseError::Invalid("solid support was already eliminated or is absent")
            })?);
        }
        let mass = dense(&self.mass, n, cx)?;
        let mut stiffness = dense(&self.stiffness, n, cx)?;
        symmetrize_native_stiffness(&mut stiffness, n, cx)?;
        let damping = match damping_ns_m {
            Some(matrix) => dense(matrix, n, cx)?,
            None => vec![0.0; n * n],
        };
        let response = PreparedSupportedHarmonicResponse::try_new(
            &mass,
            &damping,
            &stiffness,
            n,
            &support_rows,
            angular_frequency_rad_s,
            HarmonicResponseControls {
                maximum_payload_bytes: payload,
                ..controls
            },
            cx,
        )?;
        let free_dofs = response
            .free_dofs()
            .iter()
            .map(|&row| self.free_dofs[row])
            .collect();
        Ok(SupportedSolidHarmonicOperator {
            source: self,
            free_dofs,
            support_dofs: support_dofs.to_vec(),
            response,
        })
    }
}

// Native B^T D B assembly evaluates both triangles independently. Match the
// frame pencil's bounded policy: remove only stiffness-scaled roundoff skew,
// preserve every coupling, and leave all original source entries untouched.
fn symmetrize_native_stiffness(
    matrix: &mut [f64],
    n: usize,
    cx: &Cx<'_>,
) -> Result<(), HarmonicResponseError> {
    for row in 0..n {
        cx.checkpoint()
            .map_err(|_| HarmonicResponseError::Cancelled)?;
        for col in 0..=row {
            let (a, b) = (matrix[row * n + col], matrix[col * n + row]);
            if !a.is_finite() || !b.is_finite() {
                return Err(HarmonicResponseError::Invalid("nonfinite solid stiffness"));
            }
            if a == b {
                continue;
            }
            let si = matrix[row * n + row].abs().sqrt();
            let sj = matrix[col * n + col].abs().sqrt();
            let skew = (a / si / sj - b / si / sj).abs();
            if !skew.is_finite() || skew > 64.0 * n as f64 * f64::EPSILON {
                return Err(HarmonicResponseError::Invalid(
                    "solid stiffness asymmetry exceeds assembly roundoff",
                ));
            }
            let mean = 0.5 * a + 0.5 * b;
            matrix[row * n + col] = mean;
            matrix[col * n + row] = mean;
        }
    }
    Ok(())
}

fn dense(matrix: &Csr, n: usize, cx: &Cx<'_>) -> Result<Vec<f64>, HarmonicResponseError> {
    if matrix.nrows() != n || matrix.ncols() != n {
        return Err(HarmonicResponseError::Invalid(
            "solid matrix dimensions must match all represented coordinates",
        ));
    }
    let mut values = vec![0.0; n * n];
    for row in 0..n {
        cx.checkpoint()
            .map_err(|_| HarmonicResponseError::Cancelled)?;
        let (columns, entries) = matrix.row(row);
        for (&column, &entry) in columns.iter().zip(entries) {
            values[row * n + column] = entry;
        }
    }
    Ok(values)
}
