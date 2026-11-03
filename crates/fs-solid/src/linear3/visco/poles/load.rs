//! Original power-conjugate spring/damper/inertance boundary coordinates.
use super::{ViscoPoleError, cap, checkpoint, finite, invalid, poll};
use fs_exec::Cx;
use fs_la::eigen::{MAX_EIGEN_UNPOLLED_SCALAR_WORK, admit_jacobi_eigh, jacobi_eigh};
use fs_math::c64::C64;

/// Passive reciprocal physical boundary, in the supplied port units.
/// No coefficient is fitted, diagonalized, clipped or converted to a modal zeta.
#[derive(Clone, Copy, Debug)]
pub struct ViscoPolePortLoad<'a> {
    /// Number of retained physical ports; all supplied columns are kept.
    pub ports: usize,
    /// Row-major free_displacements*ports B. Opposing nodal force is -B*effort,
    /// port displacement is B^T*u, and flow is B^T*u'. This fixes the power pair.
    pub force_map: &'a [f64],
    /// Symmetric PSD ports*ports K, multiplying port displacement.
    /// For force/velocity ports its unit is N/m.
    pub stiffness: &'a [f64],
    /// Symmetric PSD ports*ports C, multiplying port flow (N*s/m mechanically).
    pub damping: &'a [f64],
    /// Symmetric PSD ports*ports M, multiplying flow derivative (kg mechanically).
    pub inertance: &'a [f64],
    /// Declared real-frequency validity [rad/s], conservatively applied to |s|.
    /// Complex continuation of this stated law is not experimental validation.
    pub band_rad_s: (f64, f64),
}

/// Additional load resources, separate from the original source/solver limits.
#[derive(Clone, Copy, Debug)]
pub struct ViscoPolePortBudget {
    /// Maximum supplied ports; redundant/null columns are not discarded.
    pub maximum_ports: usize,
    /// Retained original load, projection/PSD scratch and port observations.
    /// Excludes allocator overhead and the separately capped source pencil.
    pub maximum_payload_bytes: usize,
    /// Conservative PSD and projection scalar-work envelope, admitted before
    /// allocation. The shared Jacobi per-call unpolled cap also applies.
    pub maximum_preparation_work: u64,
}

/// Load observations for a mode under exp(s*t), with its arbitrary amplitude.
/// These are complex-envelope quadratic forms, not finite-cycle time averages.
#[derive(Clone, Debug)]
pub struct ViscoPolePortResponse {
    /// Original displacement B^T*u, in the caller's port displacement units.
    pub displacement: Vec<C64>,
    /// Original flow s*B^T*u.
    pub flow: Vec<C64>,
    /// Opposing effort (K+s*C+s²*M)*B^T*u.
    pub effort: Vec<C64>,
    /// (d^H*K*d + flow^H*M*flow)/2 [J], real part; roundoff stays signed.
    pub storage_j: f64,
    /// Re(flow^H*C*flow) [W], the boundary's own dissipative channel.
    pub dissipation_w: f64,
    /// Re(flow^H*effort) [W], power into the whole boundary law.
    pub absorbed_power_w: f64,
    /// Absorbed power - 2*Re(s)*storage - dissipation [W].
    pub signed_power_defect_w: f64,
}

pub(super) struct PreparedPortLoad {
    coordinates: usize,
    ports: usize,
    force_map: Vec<f64>,
    coefficients: [Vec<f64>; 3],
    pub(super) band_rad_s: (f64, f64),
}

impl PreparedPortLoad {
    pub(super) fn port_count(&self) -> usize {
        self.ports
    }

    pub(super) fn new(
        coordinates: usize,
        load: ViscoPolePortLoad<'_>,
        budget: ViscoPolePortBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ViscoPoleError> {
        checkpoint(cx)?;
        let p = load.ports;
        if p == 0 {
            return Err(invalid(
                "physical boundary load requires a nonempty port map",
            ));
        }
        cap("load ports", p as u128, budget.maximum_ports as u128)?;
        let np = coordinates
            .checked_mul(p)
            .ok_or_else(|| invalid("boundary force-map extent overflow"))?;
        let pp = p
            .checked_mul(p)
            .ok_or_else(|| invalid("boundary coefficient extent overflow"))?;
        let matrices = [load.stiffness, load.damping, load.inertance];
        if load.force_map.len() != np
            || load.force_map.iter().any(|x| !x.is_finite())
            || matrices
                .iter()
                .any(|a| a.len() != pp || a.iter().any(|x| !x.is_finite()))
            || !load.band_rad_s.0.is_finite()
            || load.band_rad_s.0 < 0.0
            || !load.band_rad_s.1.is_finite()
            || load.band_rad_s.1 < load.band_rad_s.0
        {
            return Err(invalid(
                "finite complete boundary matrices, force map and ordered band required",
            ));
        }
        cap(
            "load numerical payload bytes",
            24 * np as u128 + 80 * pp as u128 + 256 * p as u128,
            budget.maximum_payload_bytes as u128,
        )?;
        let cube = (p as u128).saturating_pow(3);
        let psd_work = 720u128.saturating_mul(cube) + 256 * pp as u128;
        cap(
            "load PSD unpolled work",
            psd_work,
            MAX_EIGEN_UNPOLLED_SCALAR_WORK as u128,
        )?;
        cap(
            "load preparation work",
            3 * psd_work + 6 * np as u128 * (coordinates as u128 + p as u128),
            budget.maximum_preparation_work as u128,
        )?;
        admit_jacobi_eigh(p)
            .map_err(|_| invalid("boundary PSD workspace exceeds shared capacity"))?;
        for (a, name) in matrices.iter().zip(["stiffness", "damping", "inertance"]) {
            admit_coefficient(a, p, name, cx)?;
        }
        checkpoint(cx)?;
        Ok(Self {
            coordinates,
            ports: p,
            force_map: load.force_map.to_vec(),
            coefficients: matrices.map(<[f64]>::to_vec),
            band_rad_s: load.band_rad_s,
        })
    }

    pub(super) fn add_to_pencil(
        &self,
        pencil: &mut [Vec<f64>; 3],
        dimension: usize,
        cx: &Cx<'_>,
    ) -> Result<(), ViscoPoleError> {
        let (n, p) = (self.coordinates, self.ports);
        let mut product = vec![0.0; n * p];
        for (coefficient, source) in pencil.iter_mut().zip(&self.coefficients) {
            for i in 0..n {
                checkpoint(cx)?;
                for b in 0..p {
                    let mut sum = 0.0;
                    for a in 0..p {
                        sum = self.force_map[i * p + a].mul_add(source[a * p + b], sum);
                    }
                    product[i * p + b] = sum;
                }
            }
            for i in 0..n {
                checkpoint(cx)?;
                for j in i..n {
                    let mut value = 0.0;
                    for b in 0..p {
                        value = product[i * p + b].mul_add(self.force_map[j * p + b], value);
                    }
                    if !value.is_finite() {
                        return Err(invalid("physical boundary projection overflow"));
                    }
                    // No arithmetic changes in the identically zero-load path.
                    if value != 0.0 {
                        coefficient[i * dimension + j] += value;
                        if i != j {
                            coefficient[j * dimension + i] += value;
                        }
                    }
                }
            }
        }
        checkpoint(cx)?;
        Ok(())
    }

    pub(super) fn observe(
        &self,
        rate: C64,
        motion: &[C64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<ViscoPolePortResponse, ViscoPoleError> {
        let (n, p) = (self.coordinates, self.ports);
        let mut displacement = vec![C64::ZERO; p];
        for (a, d) in displacement.iter_mut().enumerate() {
            poll(cancelled)?;
            for (i, &u) in motion.iter().enumerate().take(n) {
                *d = *d + u.scale(self.force_map[i * p + a]);
            }
        }
        let flow: Vec<_> = displacement.iter().map(|&d| rate * d).collect();
        let mut effort = vec![C64::ZERO; p];
        let mut storage = 0.0;
        let mut loss = 0.0;
        let mut power = 0.0;
        for a in 0..p {
            poll(cancelled)?;
            let mut elastic = C64::ZERO;
            let mut viscous = C64::ZERO;
            let mut momentum = C64::ZERO;
            for b in 0..p {
                elastic = elastic + displacement[b].scale(self.coefficients[0][a * p + b]);
                viscous = viscous + flow[b].scale(self.coefficients[1][a * p + b]);
                momentum = momentum + flow[b].scale(self.coefficients[2][a * p + b]);
            }
            effort[a] = elastic + viscous + rate * momentum;
            storage += 0.5 * (inner(displacement[a], elastic) + inner(flow[a], momentum));
            loss += inner(flow[a], viscous);
            power += inner(flow[a], effort[a]);
        }
        let defect = power - 2.0 * rate.re * storage - loss;
        if displacement
            .iter()
            .chain(&flow)
            .chain(&effort)
            .any(|&z| !finite(z))
            || [storage, loss, power, defect]
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(invalid("original boundary modal observation overflow"));
        }
        poll(cancelled)?;
        Ok(ViscoPolePortResponse {
            displacement,
            flow,
            effort,
            storage_j: storage,
            dissipation_w: loss,
            absorbed_power_w: power,
            signed_power_defect_w: defect,
        })
    }
}

fn inner(a: C64, b: C64) -> f64 {
    a.re.mul_add(b.re, a.im * b.im)
}

pub(super) fn admit_coefficient(
    a: &[f64],
    n: usize,
    name: &'static str,
    cx: &Cx<'_>,
) -> Result<(), ViscoPoleError> {
    checkpoint(cx)?;
    let refusal = |value| ViscoPoleError::LoadNotPassive {
        coefficient: name,
        minimum_scaled_eigenvalue: value,
    };
    let mut roots = Vec::with_capacity(n);
    for i in 0..n {
        checkpoint(cx)?;
        for j in 0..i {
            if a[i * n + j] != a[j * n + i] {
                return Err(invalid("boundary coefficients must be exactly symmetric"));
            }
        }
        if a[i * n + i] < 0.0 || a[i * n + i] == 0.0 && (0..n).any(|j| a[i * n + j] != 0.0) {
            return Err(refusal(a[i * n + i]));
        }
        roots.push(a[i * n + i].sqrt());
    }
    if (0..n).all(|i| (0..n).all(|j| i == j || a[i * n + j] == 0.0)) {
        return Ok(());
    }
    // Congruence gives scale-independent numerical admission. The ORIGINAL
    // matrices survive unchanged, including their nullspaces and mutual terms.
    let mut normalized = vec![0.0; n * n];
    let allowance = 64.0 * n as f64 * f64::EPSILON;
    for i in 0..n {
        checkpoint(cx)?;
        if roots[i] == 0.0 {
            continue;
        }
        normalized[i * n + i] = 1.0;
        for j in 0..i {
            if roots[j] == 0.0 {
                continue;
            }
            let value = (a[i * n + j] / roots[i]) / roots[j];
            if !value.is_finite() || value.abs() > 1.0 + allowance {
                return Err(refusal(1.0 - value.abs()));
            }
            normalized[i * n + j] = value;
            normalized[j * n + i] = value;
        }
    }
    let (values, _) = jacobi_eigh(&normalized, n);
    checkpoint(cx)?;
    let minimum = values.iter().copied().fold(f64::INFINITY, f64::min);
    if !minimum.is_finite() || minimum < -allowance {
        return Err(refusal(minimum));
    }
    Ok(())
}
