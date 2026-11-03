//! Original linear boundary states without an invented energy realization.
use super::super::checkpoint;
use super::{
    PreparedViscoPoles, ViscoPoleDynamicBudget, ViscoPoleError, cap, finite, invalid, poll,
};
use fs_exec::Cx;
use fs_math::c64::C64;

/// Original real continuous impedance: z'=A z+B v, effort=C z+D v+E v'.
///
/// The source owns state/port units, its realization, real-frequency support,
/// stability and passivity evidence. Matrices may be copied directly from an
/// original rational source's real state-space view; no pole fit, balancing or
/// sample-clock transformation occurs here. Every original state, including
/// dark states and both real coordinates of conjugate pairs, survives.
///
/// This interface does not infer storage or heat from signed boundary work.
/// It can analyze a declared active or unstable source but does not certify it
/// as passive radiation. An exact physical pHS load should use its own API.
#[derive(Clone, Copy, Debug)]
pub struct ViscoPoleStateSpaceLoad<'a> {
    pub states: usize,
    pub ports: usize,
    /// Row-major states*states, states*ports, ports*states.
    pub a: &'a [f64],
    pub b: &'a [f64],
    pub c: &'a [f64],
    /// Original row-major ports*ports direct and derivative terms.
    pub d: &'a [f64],
    pub e: &'a [f64],
    /// Row-major free_displacements*ports. v=force_map^T*u', and the
    /// opposing body force is -force_map*effort.
    pub force_map: &'a [f64],
    /// Declared continuation support [rad/s], conservatively applied to |s|.
    pub band_rad_s: (f64, f64),
}

/// Original boundary phasors and signed complex-envelope port work.
/// No load storage or dissipation is claimed without a physical witness.
#[derive(Clone, Debug)]
pub struct ViscoPoleStateSpaceResponse {
    pub state: Vec<C64>,
    pub displacement: Vec<C64>,
    pub flow: Vec<C64>,
    pub effort: Vec<C64>,
    /// Re(conj(flow_p)*effort_p), positive into the original boundary.
    /// Individual mutual-port exchanges may be negative.
    pub absorbed_port_power_w: Vec<f64>,
    pub absorbed_power_w: f64,
    /// Original body energy/heat only; excludes the unknown boundary stores.
    pub body_storage_j: f64,
    pub body_dissipation_w: f64,
    /// 2 Re(s) H_body + D_body + signed power into the boundary.
    pub body_decay_power_residual_w: f64,
}

pub(super) struct PreparedStateSpaceLoad {
    coordinates: usize,
    states: usize,
    ports: usize,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
    e: Vec<f64>,
    force_map: Vec<f64>,
    pub(super) band_rad_s: (f64, f64),
}

impl PreparedStateSpaceLoad {
    pub(super) fn port_count(&self) -> usize {
        self.ports
    }

    pub(super) fn new(
        coordinates: usize,
        load: ViscoPoleStateSpaceLoad<'_>,
        budget: ViscoPoleDynamicBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ViscoPoleError> {
        checkpoint(cx)?;
        let (k, p) = (load.states, load.ports);
        if p == 0 {
            return Err(invalid("state-space boundary requires physical ports"));
        }
        cap(
            "original boundary states",
            k as u128,
            budget.maximum_states as u128,
        )?;
        cap(
            "original boundary ports",
            p as u128,
            budget.maximum_ports as u128,
        )?;
        let kk = k
            .checked_mul(k)
            .ok_or_else(|| invalid("boundary state extent overflow"))?;
        let kp = k
            .checked_mul(p)
            .ok_or_else(|| invalid("boundary port extent overflow"))?;
        let pp = p
            .checked_mul(p)
            .ok_or_else(|| invalid("boundary direct extent overflow"))?;
        let mp = coordinates
            .checked_mul(p)
            .ok_or_else(|| invalid("boundary force map overflow"))?;
        if load.a.len() != kk
            || load.b.len() != kp
            || load.c.len() != kp
            || load.d.len() != pp
            || load.e.len() != pp
            || load.force_map.len() != mp
            || !load.band_rad_s.0.is_finite()
            || load.band_rad_s.0 < 0.0
            || !load.band_rad_s.1.is_finite()
            || load.band_rad_s.1 < load.band_rad_s.0
        {
            return Err(invalid(
                "complete original state-space matrices, map and ordered band required",
            ));
        }
        let (m, k128, p128) = (coordinates as u128, k as u128, p as u128);
        cap(
            "state-space boundary payload bytes",
            8 * k128 * k128
                + 16 * k128 * p128
                + 16 * p128 * p128
                + 32 * m * p128
                + 256 * (m + k128 + p128),
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "state-space boundary preparation work",
            8 * m * p128 * (m + p128 + k128) + 4 * k128 * k128 + 64 * (m + k128 + p128),
            budget.maximum_preparation_work as u128,
        )?;
        for matrix in [load.a, load.b, load.c, load.d, load.e, load.force_map] {
            for row in matrix.chunks(64) {
                checkpoint(cx)?;
                if row.iter().any(|x| !x.is_finite()) {
                    return Err(invalid("original state-space boundary is not finite"));
                }
            }
        }
        Ok(Self {
            coordinates,
            states: k,
            ports: p,
            a: load.a.to_vec(),
            b: load.b.to_vec(),
            c: load.c.to_vec(),
            d: load.d.to_vec(),
            e: load.e.to_vec(),
            force_map: load.force_map.to_vec(),
            band_rad_s: load.band_rad_s,
        })
    }

    pub(super) fn add_to_pencil(
        &self,
        a: &mut [Vec<f64>; 3],
        n: usize,
        at: usize,
        cx: &Cx<'_>,
    ) -> Result<(), ViscoPoleError> {
        let (m, k, p) = (self.coordinates, self.states, self.ports);
        // Original host map projects D and E without taking either inverse.
        let mut direct = vec![0.0; m * p];
        let mut derivative = vec![0.0; m * p];
        for i in 0..m {
            checkpoint(cx)?;
            for q in 0..p {
                direct[i * p + q] = (0..p).fold(0.0, |sum, r| {
                    self.force_map[i * p + r].mul_add(self.d[r * p + q], sum)
                });
                derivative[i * p + q] = (0..p).fold(0.0, |sum, r| {
                    self.force_map[i * p + r].mul_add(self.e[r * p + q], sum)
                });
            }
            for j in 0..m {
                a[1][i * n + j] += (0..p).fold(0.0, |sum, q| {
                    direct[i * p + q].mul_add(self.force_map[j * p + q], sum)
                });
                a[2][i * n + j] += (0..p).fold(0.0, |sum, q| {
                    derivative[i * p + q].mul_add(self.force_map[j * p + q], sum)
                });
            }
            for j in 0..k {
                a[0][i * n + at + j] = (0..p).fold(0.0, |sum, q| {
                    self.force_map[i * p + q].mul_add(self.c[q * k + j], sum)
                });
                a[1][(at + j) * n + i] = -(0..p).fold(0.0, |sum, q| {
                    self.b[j * p + q].mul_add(self.force_map[i * p + q], sum)
                });
            }
        }
        for i in 0..k {
            checkpoint(cx)?;
            a[1][(at + i) * n + at + i] = 1.0;
            for j in 0..k {
                a[0][(at + i) * n + at + j] = -self.a[i * k + j];
            }
        }
        Ok(())
    }

    pub(super) fn observe(
        &self,
        prepared: &PreparedViscoPoles<'_>,
        rate: C64,
        vector: &[C64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<ViscoPoleStateSpaceResponse, ViscoPoleError> {
        let (m, k, p) = (self.coordinates, self.states, self.ports);
        let state = &vector[prepared.body_dimension()..];
        let mut displacement = vec![C64::ZERO; p];
        let mut flow = vec![C64::ZERO; p];
        let mut effort = vec![C64::ZERO; p];
        let mut absorbed_port_power_w = vec![0.0; p];
        for q in 0..p {
            poll(cancelled)?;
            displacement[q] = (0..m).fold(C64::ZERO, |sum, i| {
                sum + vector[i].scale(self.force_map[i * p + q])
            });
            flow[q] = rate * displacement[q];
        }
        for q in 0..p {
            poll(cancelled)?;
            effort[q] = (0..k).fold(C64::ZERO, |sum, j| sum + state[j].scale(self.c[q * k + j]));
            for r in 0..p {
                effort[q] = effort[q]
                    + flow[r].scale(self.d[q * p + r])
                    + (rate * flow[r]).scale(self.e[q * p + r]);
            }
            absorbed_port_power_w[q] = (flow[q].conj() * effort[q]).re;
        }
        let absorbed_power_w = absorbed_port_power_w.iter().sum::<f64>();
        let (body_storage_j, body_dissipation_w) =
            super::dynamic::body_energy(prepared, rate, vector, cancelled)?;
        let body_decay_power_residual_w =
            2.0 * rate.re * body_storage_j + body_dissipation_w + absorbed_power_w;
        if state
            .iter()
            .chain(&displacement)
            .chain(&flow)
            .chain(&effort)
            .any(|&z| !finite(z))
            || absorbed_port_power_w.iter().any(|x| !x.is_finite())
            || [
                absorbed_power_w,
                body_storage_j,
                body_dissipation_w,
                body_decay_power_residual_w,
            ]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "original state-space boundary observation is not representable",
            ));
        }
        Ok(ViscoPoleStateSpaceResponse {
            state: state.to_vec(),
            displacement,
            flow,
            effort,
            absorbed_port_power_w,
            absorbed_power_w,
            body_storage_j,
            body_dissipation_w,
            body_decay_power_residual_w,
        })
    }
}
