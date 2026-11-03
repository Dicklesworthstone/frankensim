//! Joint poles with the original independent states of a physical pHS load.
use super::{PreparedViscoPoles, ViscoPoleError, cap, checkpoint, finite, invalid, poll};
use fs_exec::Cx;
use fs_la::eigen::{MAX_EIGEN_UNPOLLED_SCALAR_WORK, admit_jacobi_eigh};
use fs_math::c64::C64;
use fs_phs::PortHamiltonian;

/// A physical impedance law z'=(J-R)Qz+Gv, effort=G^T*Qz.
/// Every state and full mutual entry survives, including unobservable states.
/// No fitting, reduction, sampled-impedance continuation or modal damping fit.
#[derive(Clone, Copy)]
pub struct ViscoPoleDynamicLoad<'a> {
    /// Admitted pHS owner with its exact quadratic storage law, H=z^T*Qz/2.
    /// Nonlinear storage is refused, not frozen at an unstated operating point.
    pub system: &'a PortHamiltonian,
    /// Row-major free_displacements*load_ports B. The load receives B^T*u'
    /// and exerts opposing nodal force -B*effort. Units must be power conjugate.
    pub force_map: &'a [f64],
    /// Declared validity [rad/s], conservatively applied to |s|. This declaration
    /// does not validate the complex continuation experimentally.
    pub band_rad_s: (f64, f64),
}

/// Additional original-load preparation resources; total pencil/solver limits
/// remain in ViscoPoleBudget and the shared numerical solver options.
#[derive(Clone, Copy, Debug)]
pub struct ViscoPoleDynamicBudget {
    pub maximum_states: usize,
    pub maximum_ports: usize,
    /// Load matrices, projections, PSD scratch and physical observations.
    /// Excludes allocator overhead and the separately capped joint pencil.
    pub maximum_payload_bytes: usize,
    pub maximum_preparation_work: u64,
}

/// Complex-envelope quadratic forms under exp(s*t), at arbitrary modal scale.
/// They are not cycle averages and contain no inferred or fitted heat.
#[derive(Clone, Debug)]
pub struct ViscoPoleDynamicResponse {
    /// Original state order of the source pHS load, also in physical_vector's suffix.
    pub state: Vec<C64>,
    pub displacement: Vec<C64>,
    pub flow: Vec<C64>,
    pub effort: Vec<C64>,
    /// Re(z^H*Qz)/2 [J].
    pub storage_j: f64,
    /// Re((Qz)^H*R*(Qz)) [W], independently evaluated from the original R.
    pub dissipation_w: f64,
    /// Re(flow^H*effort) [W], positive INTO the dynamic load.
    pub absorbed_power_w: f64,
    /// absorbed_power - 2*Re(s)*storage - dissipation [W].
    pub signed_power_defect_w: f64,
    /// Original body equilibrium, consistent inertia and Maxwell storage [J].
    pub body_storage_j: f64,
    /// Original Kelvin viscosity and every Maxwell relaxation loss [W].
    pub body_dissipation_w: f64,
    /// 2*Re(s)*(body_storage+load_storage)+body_loss+load_loss [W].
    /// This is observed independently; it never defines either loss channel.
    pub joint_decay_power_defect_w: f64,
}

pub(super) struct PreparedDynamicLoad {
    coordinates: usize,
    states: usize,
    ports: usize,
    force_map: Vec<f64>,
    q: Vec<f64>,
    r: Vec<f64>,
    g: Vec<f64>,
    generator: Vec<f64>,
    pub(super) band_rad_s: (f64, f64),
}

impl PreparedDynamicLoad {
    pub(super) fn new(
        coordinates: usize,
        load: ViscoPoleDynamicLoad<'_>,
        budget: ViscoPoleDynamicBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, ViscoPoleError> {
        checkpoint(cx)?;
        let (k, p) = (load.system.state_dim(), load.system.port_dim());
        if k == 0 || p == 0 {
            return Err(invalid(
                "independent load requires states and physical ports",
            ));
        }
        cap(
            "dynamic load states",
            k as u128,
            budget.maximum_states as u128,
        )?;
        cap(
            "dynamic load ports",
            p as u128,
            budget.maximum_ports as u128,
        )?;
        let kk = k
            .checked_mul(k)
            .ok_or_else(|| invalid("load state extent overflow"))?;
        let kp = k
            .checked_mul(p)
            .ok_or_else(|| invalid("load port extent overflow"))?;
        let mp = coordinates
            .checked_mul(p)
            .ok_or_else(|| invalid("force map extent overflow"))?;
        if load.force_map.len() != mp
            || load.force_map.iter().any(|v| !v.is_finite())
            || !load.band_rad_s.0.is_finite()
            || load.band_rad_s.0 < 0.0
            || !load.band_rad_s.1.is_finite()
            || load.band_rad_s.1 < load.band_rad_s.0
        {
            return Err(invalid(
                "finite complete dynamic force map and ordered band required",
            ));
        }
        cap(
            "dynamic load numerical payload bytes",
            96 * kk as u128
                + 32 * kp as u128
                + 32 * mp as u128
                + 16 * k as u128 * coordinates as u128
                + 128 * (k as u128 + p as u128),
            budget.maximum_payload_bytes as u128,
        )?;
        let psd_work = 720 * (k as u128).pow(3) + 256 * kk as u128;
        cap(
            "dynamic load PSD unpolled work",
            psd_work,
            MAX_EIGEN_UNPOLLED_SCALAR_WORK as u128,
        )?;
        cap(
            "dynamic load preparation work",
            2 * psd_work
                + 4 * (k as u128).pow(3)
                + 4 * p as u128 * kk as u128
                + 4 * coordinates as u128 * kk as u128
                + 4 * coordinates as u128 * kp as u128,
            budget.maximum_preparation_work as u128,
        )?;
        admit_jacobi_eigh(k)
            .map_err(|_| invalid("dynamic load PSD workspace exceeds shared capacity"))?;
        let part = load
            .system
            .storage()
            .quadratic_part()
            .ok_or_else(|| invalid("dynamic load has no declared quadratic storage"))?;
        if !part.exact || part.matrix.len() != kk || part.matrix.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "dynamic load must supply its complete finite quadratic law",
            ));
        }
        let (j, r, g) = load.system.structure();
        super::load::admit_coefficient(&part.matrix, k, "dynamic storage", cx)?;
        super::load::admit_coefficient(r, k, "dynamic resistance", cx)?;
        let mut generator = vec![0.0; kk];
        for i in 0..k {
            checkpoint(cx)?;
            for col in 0..k {
                generator[i * k + col] = (0..k).fold(0.0, |sum, row| {
                    (j[i * k + row] - r[i * k + row]).mul_add(part.matrix[row * k + col], sum)
                });
            }
        }
        if generator.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "original dynamic load generator is not representable",
            ));
        }
        checkpoint(cx)?;
        Ok(Self {
            coordinates,
            states: k,
            ports: p,
            force_map: load.force_map.to_vec(),
            q: part.matrix,
            r: r.to_vec(),
            g: g.to_vec(),
            generator,
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
        // B*G^T is used in both directions, preserving the power-conjugate map.
        let mut coupling = vec![0.0; m * k];
        for i in 0..m {
            checkpoint(cx)?;
            for j in 0..k {
                coupling[i * k + j] = (0..p).fold(0.0, |sum, port| {
                    self.force_map[i * p + port].mul_add(self.g[j * p + port], sum)
                });
            }
            for j in 0..k {
                a[0][i * n + at + j] = (0..k).fold(0.0, |sum, row| {
                    coupling[i * k + row].mul_add(self.q[row * k + j], sum)
                });
                a[1][(at + j) * n + i] = -coupling[i * k + j];
            }
        }
        for i in 0..k {
            checkpoint(cx)?;
            a[1][(at + i) * n + at + i] = 1.0;
            for j in 0..k {
                a[0][(at + i) * n + at + j] = -self.generator[i * k + j];
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
    ) -> Result<ViscoPoleDynamicResponse, ViscoPoleError> {
        let (m, k, p) = (self.coordinates, self.states, self.ports);
        let state = &vector[prepared.body_dimension()..];
        let mut gradient = vec![C64::ZERO; k];
        let mut storage = 0.0;
        for i in 0..k {
            poll(cancelled)?;
            gradient[i] = (0..k).fold(C64::ZERO, |sum, j| sum + state[j].scale(self.q[i * k + j]));
            storage += 0.5 * inner(state[i], gradient[i]);
        }
        let mut loss = 0.0;
        for i in 0..k {
            poll(cancelled)?;
            let resistance = (0..k).fold(C64::ZERO, |sum, j| {
                sum + gradient[j].scale(self.r[i * k + j])
            });
            loss += inner(gradient[i], resistance);
        }
        let mut displacement = vec![C64::ZERO; p];
        let mut flow = vec![C64::ZERO; p];
        let mut effort = vec![C64::ZERO; p];
        let mut power = 0.0;
        for port in 0..p {
            poll(cancelled)?;
            displacement[port] = (0..m).fold(C64::ZERO, |sum, i| {
                sum + vector[i].scale(self.force_map[i * p + port])
            });
            flow[port] = rate * displacement[port];
            effort[port] = (0..k).fold(C64::ZERO, |sum, i| {
                sum + gradient[i].scale(self.g[i * p + port])
            });
            power += inner(flow[port], effort[port]);
        }
        let (body_storage, body_loss) = body_energy(prepared, rate, vector, cancelled)?;
        let defect = power - 2.0 * rate.re * storage - loss;
        let joint_defect = 2.0 * rate.re * (body_storage + storage) + body_loss + loss;
        if state
            .iter()
            .chain(&gradient)
            .chain(&displacement)
            .chain(&flow)
            .chain(&effort)
            .any(|&v| !finite(v))
            || [
                storage,
                loss,
                power,
                defect,
                body_storage,
                body_loss,
                joint_defect,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "original dynamic load observation is not representable",
            ));
        }
        Ok(ViscoPoleDynamicResponse {
            state: state.to_vec(),
            displacement,
            flow,
            effort,
            storage_j: storage,
            dissipation_w: loss,
            absorbed_power_w: power,
            signed_power_defect_w: defect,
            body_storage_j: body_storage,
            body_dissipation_w: body_loss,
            joint_decay_power_defect_w: joint_defect,
        })
    }
}

fn inner(a: C64, b: C64) -> f64 {
    (a.conj() * b).re
}

fn body_energy(
    prepared: &PreparedViscoPoles<'_>,
    rate: C64,
    vector: &[C64],
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<(f64, f64), ViscoPoleError> {
    let m = prepared.free_dofs().len();
    let n = prepared.dimension;
    let mut storage = 0.0;
    let mut loss = 0.0;
    // Dynamic coupling has no host-host entries. These blocks remain the
    // original equilibrium stiffness, Kelvin viscosity and consistent mass.
    for i in 0..m {
        poll(cancelled)?;
        for j in 0..m {
            let at = i * n + j;
            storage += 0.5
                * (inner(vector[i], vector[j]) * prepared.coefficients[0][at]
                    + inner(rate * vector[i], rate * vector[j]) * prepared.coefficients[2][at]);
            loss += inner(rate * vector[i], rate * vector[j]) * prepared.coefficients[1][at];
        }
    }
    for (index, site) in prepared.sites.iter().enumerate() {
        poll(cancelled)?;
        let element = &prepared.source.elements[site.element];
        let branch = &element.law.branches()[site.branch];
        let q = &vector[m + 6 * index..m + 6 * (index + 1)];
        let mut energy_twice = 0.0;
        for i in 0..6 {
            for j in 0..6 {
                energy_twice += inner(q[i], q[j]) * branch.c_pa[i][j];
            }
        }
        energy_twice *= element.points[site.point].volume_m3;
        storage += 0.5 * energy_twice;
        loss += energy_twice / branch.tau_s;
    }
    Ok((storage, loss))
}
