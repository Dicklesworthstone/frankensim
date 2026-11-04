//! Physical effort/flow loads for the original coupled harmonic pencil.
//!
//! A real force map B defines both reaction force B*effort and flow B^T*v.
//! The complete load Z is evaluated by its physical owner before this call.
//! Nothing here fits poles, selects material data, or advances a load history.

use super::{
    HarmonicResponseError, finite, invalid, poll, real_inner, require_psd, scale, scale_complex,
};
use crate::quadratic_fraction;
use fs_exec::Cx;
use fs_la::eigen::admit_jacobi_eigh;
use fs_math::c64::C64;

/// Additional limits for physical-port projection and pointwise admission.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPortControls {
    /// Maximum supplied ports. No port is discarded to satisfy this cap.
    pub maximum_ports: usize,
    /// Maximum `n*p*(n+p)` multiply-add terms in B*Z*B^T.
    /// Pointwise PSD admission has the existing fs-la bounded Jacobi work cap.
    pub maximum_projection_terms: u64,
    /// Conservative additional numerical payload, including source copies,
    /// projection and Hermitian-PSD scratch. Excludes allocator metadata.
    pub maximum_payload_bytes: usize,
}

/// Original physical load variables, in the caller's declared port units.
#[derive(Clone, Debug)]
pub struct HarmonicPortResponse {
    /// Peak conjugate flow amplitudes `B^T v`.
    pub flow: Vec<C64>,
    /// Opposing effort amplitudes `Z flow`; structural reaction is `-B effort`.
    pub effort: Vec<C64>,
    /// `Re(flow^H effort)/2` [W] for power-conjugate physical ports.
    /// This is net absorbed boundary work, not an inferred split into heat,
    /// radiation or stored energy. Signed roundoff is preserved.
    pub average_absorbed_power_w: f64,
}

#[derive(Clone, Debug)]
pub(super) struct PreparedPortImpedance {
    pub(super) coordinates: usize,
    pub(super) ports: usize,
    pub(super) force_map: Vec<f64>,
    pub(super) impedance: Vec<C64>,
}

pub(super) struct PortScaling {
    // B_ia / 2^(coordinate_i + port_a) has magnitude < 2. Port changes
    // are numerical congruences only; all returned flows/efforts use source units.
    pub(super) exponents: Vec<i32>,
    pub(super) minimum_time_exponent: i32,
}

impl PreparedPortImpedance {
    pub(super) fn new(
        n: usize,
        force_map: &[f64],
        ports: usize,
        impedance: &[C64],
        controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<Self, HarmonicResponseError> {
        poll(cx)?;
        if ports == 0 {
            return Err(invalid("physical load requires at least one port"));
        }
        cap("load ports", ports, controls.maximum_ports)?;
        let np = n
            .checked_mul(ports)
            .ok_or(invalid("force map extent overflow"))?;
        let pp = ports
            .checked_mul(ports)
            .ok_or(invalid("impedance extent overflow"))?;
        let real_dimension = ports
            .checked_mul(2)
            .ok_or(invalid("Hermitian extent overflow"))?;
        if force_map.len() != np
            || impedance.len() != pp
            || force_map.iter().any(|x| !x.is_finite())
            || impedance.iter().any(|z| !finite(*z))
        {
            return Err(invalid(
                "finite n*ports force map and ports*ports impedance required",
            ));
        }
        let terms = np
            .checked_mul(
                n.checked_add(ports)
                    .ok_or(invalid("projection work overflow"))?,
            )
            .and_then(|x| u64::try_from(x).ok())
            .ok_or(invalid("projection work overflow"))?;
        if terms > controls.maximum_projection_terms {
            return Err(invalid("physical port projection exceeds its work budget"));
        }
        let bytes = np
            .checked_mul(64)
            .and_then(|x| pp.checked_mul(256).and_then(|y| x.checked_add(y)))
            .and_then(|x| ports.checked_mul(128).and_then(|y| x.checked_add(y)))
            .ok_or(invalid("physical port payload overflow"))?;
        cap(
            "load numerical payload bytes",
            bytes,
            controls.maximum_payload_bytes,
        )?;
        admit_jacobi_eigh(real_dimension)
            .map_err(|_| invalid("load Hermitian PSD admission exceeds fs-la work cap"))?;

        // The Hermitian part governs work even for nonreciprocal loads. Its
        // real representation [Re(H), -Im(H); Im(H), Re(H)] is symmetric and
        // duplicates the eigenvalues; reuse the original numerical PSD owner.
        // Midpoints avoid overflow while preserving zero dissipative channels.
        let mut hermitian = vec![0.0; real_dimension * real_dimension];
        for a in 0..ports {
            poll(cx)?;
            for b in 0..=a {
                let za = impedance[a * ports + b];
                let zb = impedance[b * ports + a];
                let re = f64::midpoint(za.re, zb.re);
                let im = if a == b {
                    0.0
                } else {
                    f64::midpoint(za.im, -zb.im)
                };
                for offset in [0, ports] {
                    hermitian[(a + offset) * real_dimension + b + offset] = re;
                    hermitian[(b + offset) * real_dimension + a + offset] = re;
                }
                hermitian[a * real_dimension + b + ports] = -im;
                hermitian[(b + ports) * real_dimension + a] = -im;
                hermitian[b * real_dimension + a + ports] = im;
                hermitian[(a + ports) * real_dimension + b] = im;
            }
        }
        require_psd(&hermitian, real_dimension, "load Hermitian part", cx)?;
        poll(cx)?;
        Ok(Self {
            coordinates: n,
            ports,
            force_map: force_map.to_vec(),
            impedance: impedance.to_vec(),
        })
    }

    pub(super) fn scaling(
        &self,
        coordinates: &[i32],
        cx: &Cx<'_>,
    ) -> Result<PortScaling, HarmonicResponseError> {
        let p = self.ports;
        let mut exponents = vec![i32::MIN; p];
        for i in 0..self.coordinates {
            poll(cx)?;
            for a in 0..p {
                let b = self.force_map[i * p + a];
                if b != 0.0 {
                    exponents[a] = exponents[a].max(quadratic_fraction(b).1 - coordinates[i]);
                }
            }
        }
        let mut minimum_time_exponent = i32::MIN;
        let headroom = 2 * p.next_power_of_two().ilog2() as i32 + 2;
        for a in 0..p {
            poll(cx)?;
            if exponents[a] == i32::MIN {
                continue;
            }
            for b in 0..p {
                if exponents[b] == i32::MIN {
                    continue;
                }
                let z = self.impedance[a * p + b];
                for value in [z.re, z.im] {
                    if value != 0.0 {
                        minimum_time_exponent = minimum_time_exponent.max(
                            exponents[a] + exponents[b] + quadratic_fraction(value).1 + headroom,
                        );
                    }
                }
            }
        }
        Ok(PortScaling {
            exponents,
            minimum_time_exponent,
        })
    }

    #[allow(clippy::too_many_arguments)] // original pencil's binary scaling
    pub(super) fn add_to_operator(
        &self,
        operator: &mut [C64],
        coordinates: &[i32],
        time_exponent: i32,
        frequency: f64,
        scaling: &PortScaling,
        cx: &Cx<'_>,
    ) -> Result<(), HarmonicResponseError> {
        let (n, p) = (self.coordinates, self.ports);
        let mut b = vec![0.0; n * p];
        let mut z = vec![C64::ZERO; p * p];
        for i in 0..n {
            poll(cx)?;
            for a in 0..p {
                if scaling.exponents[a] != i32::MIN {
                    b[i * p + a] = scale(
                        self.force_map[i * p + a],
                        -coordinates[i] - scaling.exponents[a],
                    )?;
                }
            }
        }
        for a in 0..p {
            poll(cx)?;
            if scaling.exponents[a] == i32::MIN {
                continue;
            }
            for c in 0..p {
                if scaling.exponents[c] != i32::MIN {
                    z[a * p + c] = scale_complex(
                        self.impedance[a * p + c],
                        scaling.exponents[a] + scaling.exponents[c] - time_exponent,
                    )?;
                }
            }
        }
        // O(n*p^2 + n^2*p), not O(n^2*p^2). No inverse of either the bare
        // pencil or Z is needed, so bare poles and lossless load channels work.
        let mut bz = vec![C64::ZERO; n * p];
        for i in 0..n {
            poll(cx)?;
            for a in 0..p {
                for c in 0..p {
                    bz[i * p + c] = bz[i * p + c] + z[a * p + c].scale(b[i * p + a]);
                }
            }
        }
        for i in 0..n {
            poll(cx)?;
            for j in 0..n {
                let mut value = C64::ZERO;
                for a in 0..p {
                    value = value + bz[i * p + a].scale(b[j * p + a]);
                }
                operator[i * n + j] =
                    operator[i * n + j] + C64::new(value.im * frequency, -value.re * frequency);
                if !finite(operator[i * n + j]) {
                    return Err(invalid("loaded dynamic stiffness is not representable"));
                }
            }
        }
        Ok(())
    }

    pub(super) fn response(
        &self,
        velocity: &[C64],
        cx: &Cx<'_>,
    ) -> Result<HarmonicPortResponse, HarmonicResponseError> {
        let p = self.ports;
        let mut flow = vec![C64::ZERO; p];
        for (i, &v) in velocity.iter().enumerate() {
            poll(cx)?;
            for a in 0..p {
                flow[a] = flow[a] + v.scale(self.force_map[i * p + a]);
            }
        }
        let mut effort = vec![C64::ZERO; p];
        let mut power = 0.0;
        for a in 0..p {
            poll(cx)?;
            for b in 0..p {
                effort[a] = effort[a] + self.impedance[a * p + b] * flow[b];
            }
            power += 0.5 * real_inner(flow[a], effort[a]);
        }
        if flow.iter().chain(&effort).any(|z| !finite(*z)) || !power.is_finite() {
            return Err(invalid(
                "physical port response or work is not representable",
            ));
        }
        Ok(HarmonicPortResponse {
            flow,
            effort,
            average_absorbed_power_w: power,
        })
    }
}

fn cap(quantity: &'static str, requested: usize, cap: usize) -> Result<(), HarmonicResponseError> {
    if requested > cap {
        Err(HarmonicResponseError::Budget {
            quantity,
            requested,
            cap,
        })
    } else {
        Ok(())
    }
}
