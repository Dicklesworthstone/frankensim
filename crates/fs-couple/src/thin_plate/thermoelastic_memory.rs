//! Full thermoelastic memory on an original DKT modal basis.
//!
//! Each spatial row is sqrt(dA) * f^T * B * Phi, with delta D = f f^T.
//! QR compresses these rows without a numerical rank cutoff. Every retained
//! thermal diffusion mode then scales the SAME full spatial Gram. The existing
//! pHS Maxwell owner supplies storage, heat, exact linear flow and nonlinear
//! composition; no thermal zeta is added to the mechanical resistance.
//!
//! A finite spectrum's error bound is harmonic and band-limited. It is not a
//! time-domain trajectory bound, nor does this constructor band-limit forcing
//! or the nonlinear response. Geometry and elastic-mode truncation are separate.

use super::{
    AcousticRealizeError, AnisotropicPlateThermal, PlateMesh, ThinPlate, VkBody, map_plate,
    plane_stress_section,
};
use fs_exec::CancelGate;
use fs_material::visco::slab_ted::{SlabPronySpectrum, SlabThermoelasticBending};
use fs_math::det;

/// Explicit thermal accuracy and preparation resources. No tolerance or
/// material scale is inferred. These limits add no work to ordinary stepping.
#[derive(Clone, Copy, Debug)]
pub struct ThermoelasticMemoryOptions {
    /// Harmonic approximation band [rad/s], with 0 <= lower < upper.
    pub band_rad_s: (f64, f64),
    /// Absolute scalar response error per unit full thermal relaxation strength.
    /// This is not a relative error in a nearly zero loss component.
    pub normalized_tolerance: f64,
    /// Maximum retained through-thickness thermal diffusion terms (1..=64).
    pub max_thermal_terms: usize,
    /// Maximum appended scalar memory coordinates (0..=256).
    pub max_memory_states: usize,
    /// Bounds DKT projection, original field reads and dense QR arithmetic.
    pub max_projection_terms: usize,
}

pub(super) struct PreparedPlateThermalMemory<'a> {
    pub law: SlabThermoelasticBending,
    pub spectrum: &'a SlabPronySpectrum,
    pub options: ThermoelasticMemoryOptions,
    pub gate: &'a CancelGate,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

impl VkBody {
    /// Prepare the original sampled nonlinear plate with full thermal memory.
    ///
    /// Normal/shear thermal expansion, every retained cross-mode term and
    /// thermal stiffness dispersion remain in the shared physical state.
    /// The capacity at fixed in-plane strain is C_sigma - T alpha^T Q alpha;
    /// it must be positive. Elastic Q is the isothermal equilibrium stiffness.
    /// Conservative prestress is kept in the original mechanical storage.
    ///
    /// The returned spectrum records the omitted diffusion tail and coefficient
    /// error over the declared band. Nominal frequencies and a conservative
    /// instantaneous linear-stiffness bound must fit that band. This is not a
    /// bandlimit on arbitrary forcing or nonlinear motion, a spatial error
    /// bound, or an experimental qualification of the supplied constants.
    ///
    /// Homogeneous frozen-state material, insulated faces, through-thickness
    /// Fourier conduction and linear thermoelastic constitutive coupling are
    /// assumed. The nonlinear membrane storage remains the original owner;
    /// its finite-amplitude membrane thermal loss is not modeled. This explicit
    /// entry retains the complete requested sampled basis under the original
    /// modal-window and shared render-capacity admission.
    ///
    /// # Errors
    /// Duplicate thermal/authored loss, invalid original source/viscosity data,
    /// nonpositive coupled capacity, exhausted thermal/projection/state budgets,
    /// cancellation, or original FE, pHS and render admission failures.
    pub fn from_plate_with_thermoelastic_memory(
        plate: ThinPlate,
        thermal: AnisotropicPlateThermal,
        options: ThermoelasticMemoryOptions,
        gate: &CancelGate,
    ) -> Result<(Self, SlabPronySpectrum), AcousticRealizeError> {
        poll(gate)?;
        if options.max_memory_states > 256 || options.max_projection_terms == 0 {
            return Err(invalid(
                "thermal memory requires explicit bounded projection and state allowances",
            ));
        }
        let section = plane_stress_section(plate)?;
        let law = thermal.resolve_coupled(plate, &section)?;
        let spectrum = law
            .lower_to_prony(
                options.band_rad_s,
                options.normalized_tolerance,
                options.max_thermal_terms,
            )
            .map_err(|_| {
                invalid("thermal diffusion spectrum cannot satisfy its accuracy or term allowance")
            })?;
        let prepared = PreparedPlateThermalMemory {
            law,
            spectrum: &spectrum,
            options,
            gate,
        };
        let body = Self::from_sampled_fe_with_slab(plate, false, None, Some(&prepared))?;
        poll(gate)?;
        Ok((body, spectrum))
    }
}

/// Project the original finite-element fields, with one explicit coordinate
/// scale per column. This is reusable by linear and nonlinear consumers; the
/// caller retains the mechanical mass/stiffness and port maps of that basis.
///
/// Returned projections use interleaved [q0, p0, q1, p1, ...] coordinates and
/// have zero momentum entries. All rows, including weak independent rows, are
/// retained except exactly zero rows. QR is an orthogonal change of the
/// homogeneous spatial thermal coordinates, not an eigenvalue truncation.
/// Its floating-point projection error is separate from the scalar spectrum
/// bound. No material source authority is created by this numeric image.
#[allow(clippy::too_many_arguments)]
pub(crate) fn project_thermal_relaxation(
    mesh: &PlateMesh,
    dof_map: &[Option<usize>],
    modes: &[&[f64]],
    coordinate_scale: &[f64],
    law: SlabThermoelasticBending,
    spectrum: &SlabPronySpectrum,
    max_memory_states: usize,
    max_projection_terms: usize,
    gate: &CancelGate,
) -> Result<Vec<fs_phs::RelaxationBranch>, AcousticRealizeError> {
    poll(gate)?;
    let n = modes.len();
    if n == 0
        || n > 64
        || mesh.tris.is_empty()
        || mesh.nodes.len().checked_mul(3) != Some(dof_map.len())
        || coordinate_scale.len() != n
        || coordinate_scale.iter().any(|s| !s.is_finite() || *s <= 0.0)
        || max_memory_states > 256
        || max_projection_terms == 0
        || spectrum.terms().is_empty()
        || spectrum.terms().len() > 64
        || spectrum.terms()[0].relaxation_time_s().to_bits()
            != law.first_relaxation_time_s().to_bits()
    {
        return Err(invalid(
            "thermal projection needs one original bounded basis and matching slab spectrum",
        ));
    }
    let free = modes[0].len();
    let count = mesh
        .tris
        .len()
        .checked_mul(3)
        .ok_or_else(|| invalid("thermal quadrature count overflows"))?
        .max(n);
    let entries = count
        .checked_mul(n)
        .ok_or_else(|| invalid("thermal projection extent overflows"))?;
    let work = entries
        .checked_mul(8 * n + 64)
        .and_then(|w| free.checked_mul(n).and_then(|fields| w.checked_add(fields)))
        .ok_or_else(|| invalid("thermal projection work overflows"))?;
    if work > max_projection_terms {
        return Err(invalid(
            "thermal projection exceeds its explicit work allowance",
        ));
    }
    if modes
        .iter()
        .any(|v| v.len() != free || v.iter().any(|x| !x.is_finite()))
        || dof_map.iter().flatten().any(|&i| i >= free)
        || mesh.tris.iter().flatten().any(|&i| i >= mesh.nodes.len())
    {
        return Err(invalid(
            "thermal projection needs finite matching original modal fields and DOFs",
        ));
    }
    let factor = law.relaxation_factor();
    let mut rows = vec![0.0; entries];
    let mut row = 0;
    for (element, tri) in mesh.tris.iter().enumerate() {
        poll(gate)?;
        let (curvature, twice_area) = fs_plate::dkt_curvature_quadrature(
            &tri.map(|i| mesh.nodes[i].0),
            &tri.map(|i| mesh.nodes[i].1),
            element,
        )
        .map_err(map_plate)?;
        let weight = twice_area / 6.0;
        if !weight.is_finite() || weight <= 0.0 {
            return Err(invalid(
                "thermal projection requires positive original DKT quadrature",
            ));
        }
        let scale = det::sqrt(weight);
        let dofs: [Option<usize>; 9] = core::array::from_fn(|j| dof_map[3 * tri[j / 3] + j % 3]);
        for b in curvature {
            for (column, (phi, &unit)) in modes.iter().zip(coordinate_scale).enumerate() {
                rows[row * n + column] = scale * thermal_drive(factor, &b, &dofs, phi, unit);
            }
            row += 1;
        }
    }
    if rows.iter().any(|v| !v.is_finite()) {
        return Err(invalid("thermal curvature projection is unrepresentable"));
    }
    let factors = compress_rows(&mut rows, count, n, gate)?;
    let states = factors
        .len()
        .checked_mul(spectrum.terms().len())
        .ok_or_else(|| invalid("thermal memory count overflows"))?;
    if states > max_memory_states {
        return Err(invalid(
            "full thermal memory exceeds its explicit state allowance",
        ));
    }
    let mut branches = Vec::with_capacity(states);
    for term in spectrum.terms() {
        poll(gate)?;
        for factor in &factors {
            let mut projection = vec![0.0; 2 * n];
            for i in 0..n {
                projection[2 * i] = factor[i];
            }
            branches.push(fs_phs::RelaxationBranch {
                projection,
                stiffness: term.weight(),
                relaxation_time_s: term.relaxation_time_s(),
            });
        }
    }
    poll(gate)?;
    Ok(branches)
}

// Retain the material factor until AFTER each original strain projection.
// Precontracting f^T B into rounded coefficients can erase a weak channel that
// becomes observable when the larger strains cancel.
fn thermal_drive(
    factor: [f64; 3],
    curvature: &[[f64; 9]; 3],
    dofs: &[Option<usize>; 9],
    mode: &[f64],
    coordinate_scale: f64,
) -> f64 {
    let strain: [f64; 3] =
        core::array::from_fn(|k| {
            compensated_sum((0..9).map(|j| {
                dofs[j].map_or(0.0, |dof| curvature[k][j] * (mode[dof] * coordinate_scale))
            }))
        });
    compensated_sum((0..3).map(|k| factor[k] * strain[k]))
}

// Neumaier compensation for the fixed three- and nine-term source reductions.
fn compensated_sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

// Preserve weak columns across physical units before the existing QR owner.
// A power-of-two scale has no rounding unless it loses subnormal bits; reject
// that case explicitly instead of erasing a represented source coefficient.
#[allow(clippy::float_cmp)]
pub(crate) fn compress_rows(
    rows: &mut [f64],
    count: usize,
    n: usize,
    gate: &CancelGate,
) -> Result<Vec<Vec<f64>>, AcousticRealizeError> {
    let mut units = vec![1.0; n];
    for j in 0..n {
        poll(gate)?;
        let maximum = (0..count)
            .map(|i| rows[i * n + j].abs())
            .fold(0.0_f64, f64::max);
        if maximum == 0.0 {
            continue;
        }
        units[j] = f64::from_bits((maximum.to_bits() & 0x7ff0_0000_0000_0000).max(1_u64 << 52));
        for i in 0..count {
            let original = rows[i * n + j];
            let scaled = original / units[j];
            if !scaled.is_finite() || scaled * units[j] != original {
                return Err(invalid(
                    "thermal QR scaling loses an original source coefficient",
                ));
            }
            rows[i * n + j] = scaled;
        }
    }
    let qr = fs_la::factor::qr(rows, count, n);
    let mut factors = Vec::with_capacity(n);
    for i in 0..n {
        poll(gate)?;
        let row: Vec<f64> = (0..n).map(|j| qr.r(i, j) * units[j]).collect();
        if row.iter().any(|v| !v.is_finite()) {
            return Err(invalid("thermal QR image is unrepresentable"));
        }
        if row.iter().any(|v| *v != 0.0) {
            factors.push(row);
        }
    }
    Ok(factors)
}

/// A conservative Gershgorin upper bound on the represented instantaneous
/// mechanical stiffness, in the original mass-normalized coordinates. This
/// admits the frozen linear body only; quartic motion is not band-limited.
pub(super) fn admit_memory_band(
    omegas: &[f64],
    branches: &[fs_phs::RelaxationBranch],
    band: (f64, f64),
) -> Result<(), AcousticRealizeError> {
    if omegas
        .iter()
        .any(|w| !w.is_finite() || *w < band.0 || *w > band.1)
    {
        return Err(invalid(
            "thermal approximation band excludes a retained nominal mode",
        ));
    }
    let mut upper: Vec<f64> = omegas.iter().map(|w| (w * w).next_up()).collect();
    for arm in branches {
        let root = det::sqrt(arm.stiffness);
        let values: Vec<f64> = (0..omegas.len())
            .map(|i| (root * arm.projection[2 * i]).abs())
            .collect();
        let norm = values.iter().fold(
            0.0_f64,
            |sum, &v| if v == 0.0 { sum } else { (sum + v).next_up() },
        );
        for (bound, &value) in upper.iter_mut().zip(&values) {
            if value != 0.0 {
                *bound = (*bound + (value * norm).next_up()).next_up();
            }
        }
    }
    let limit = (band.1 * band.1).next_down();
    if upper.iter().any(|v| !v.is_finite() || *v > limit) {
        return Err(invalid(
            "thermal stiffening exceeds the declared frozen linear approximation band",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
