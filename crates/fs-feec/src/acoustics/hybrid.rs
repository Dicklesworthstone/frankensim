//! Power-preserving modal connection of original Helmholtz and resolved LNS regions.
//!
//! The owners use different dual coordinates. Helmholtz imposes inward velocity
//! psi_H Q_H/A_H and measures integral(psi_H p)/A_H. LNS imposes total normal
//! traction psi_L P_L and measures Q_L=-integral(psi_L v.n). Their interface is
//! therefore P_H=T^T P_L, Q_L=-T Q_H, where
//! T[l,h]=integral(psi_L[l] psi_H[h])/A_H[h]. Equating the two modal vectors
//! directly is generally wrong, including for nonuniform shapes on equal areas.
//!
//! Original source traces determine T. Coincident triangles, opposing outward
//! normals and a full-rank retained pairing are required; vertex numbering and
//! winding may differ. Exact P1 triangle moments retain off-diagonal overlap.
//! This is a finite-dimensional Galerkin interface, not pointwise conforming
//! pressure/velocity matching of the full P2 trace. Mode and mesh convergence,
//! and the placement of the interface, remain separate physical obligations.
//!
//! LNS port temperature is insulated; total normal traction couples to acoustic
//! pressure and tangential traction is free. A cut through an important viscous
//! or thermal layer is a model-form approximation, even when port work balances.
//! See COMSOL 6.3, Acoustic--Thermoviscous Acoustic Boundary, which describes
//! these normal-stress/normal-acceleration and adiabatic interface assumptions.
//! https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_multiphysics_couplings.14.12.html
//!
//! Each original loaded FE operator is factored once. A second, at most 128 by
//! 128, interface factor joins their power scattering matrices; no impedance
//! inversion or monolithic combined volume matrix is formed. Arbitrary waves,
//! original fields, and independently integrated work reuse those factors.
//! All fields use exp(+i omega t). Real-reference power coordinates also retain
//! reactive traces; they are not conventional evanescent traveling waves.

use super::helmholtz::boundary_layer::BoundaryLayerPower;
use super::helmholtz::multiport::{
    InteriorMultiport, InteriorPortResponse, PhysicalTrace, PreparedInteriorScattering,
};
use super::helmholtz::{AcousticFeError, Medium};
use super::thermoviscous::{
    PreparedThermoviscousScattering, ThermoviscousError, ThermoviscousFe, ThermoviscousField,
    ThermoviscousMedium, ThermoviscousPower,
};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, det};
use std::collections::{BTreeMap, BTreeSet};

/// Original-owner and interface failures; no mode is silently discarded.
#[derive(Debug, Clone, PartialEq)]
pub enum HybridError {
    /// Invalid indices, incompatible physical traces, or unresolved pairing.
    InvalidInput(&'static str),
    /// Cancellation before a complete prepared object or response exists.
    Cancelled,
    /// Original Helmholtz owner refused.
    Helmholtz(AcousticFeError),
    /// Original resolved velocity/pressure/temperature owner refused.
    Thermoviscous(ThermoviscousError),
    /// Original small complex LU refused a singular interface or cancellation.
    Factor(LuComplexError),
    /// Finite inputs overflowed in coordinate conversion or interface closure.
    NonFinite,
}
impl std::fmt::Display for HybridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FS-FEEC-HYBRID: {self:?}")
    }
}
impl std::error::Error for HybridError {}
impl From<AcousticFeError> for HybridError {
    fn from(e: AcousticFeError) -> Self {
        Self::Helmholtz(e)
    }
}
impl From<ThermoviscousError> for HybridError {
    fn from(e: ThermoviscousError) -> Self {
        Self::Thermoviscous(e)
    }
}
impl From<LuComplexError> for HybridError {
    fn from(e: LuComplexError) -> Self {
        if e == LuComplexError::Cancelled {
            Self::Cancelled
        } else {
            Self::Factor(e)
        }
    }
}

fn poll(cancelled: &dyn Fn() -> bool) -> Result<(), HybridError> {
    if cancelled() {
        Err(HybridError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Source-bound interface pairing. Read-only after construction by the hybrid
/// owner; callers cannot substitute a coincident-area but different trace.
#[derive(Clone, Debug)]
pub struct ModalInterface {
    helmholtz_modes: Vec<usize>,
    thermoviscous_modes: Vec<usize>,
    flow_map: Vec<f64>,
}
impl ModalInterface {
    /// Original Helmholtz mode indices, in the pairing's column order.
    #[must_use]
    pub fn helmholtz_modes(&self) -> &[usize] {
        &self.helmholtz_modes
    }
    /// Original LNS mode indices, in the pairing's row order.
    #[must_use]
    pub fn thermoviscous_modes(&self) -> &[usize] {
        &self.thermoviscous_modes
    }
    /// Row-major T, with Q_L=-T Q_H and P_H=T^T P_L. Dimensionless.
    #[must_use]
    pub fn flow_map(&self) -> &[f64] {
        &self.flow_map
    }
}

type FaceKey = [[u64; 3]; 3];
type FaceModes<'a> = BTreeMap<FaceKey, Vec<(usize, &'a PhysicalTrace)>>;

fn collect_faces<'a>(modes: &'a [Vec<PhysicalTrace>], selected: &[usize]) -> FaceModes<'a> {
    let mut faces = BTreeMap::new();
    for (local, &global) in selected.iter().enumerate() {
        for trace in &modes[global] {
            faces
                .entry(trace.vertices)
                .or_insert_with(Vec::new)
                .push((local, trace));
        }
    }
    faces
}

/// Immutable original owners connected at explicitly selected spatial modes.
pub struct HybridInterior<'a, 'c> {
    helmholtz: &'a InteriorMultiport,
    thermoviscous: &'a ThermoviscousFe<'c>,
    interface: ModalInterface,
    internal: Vec<usize>,
    external: Vec<usize>,
}

impl<'a, 'c> HybridInterior<'a, 'c> {
    /// Pair modes on the same physical triangulation. The coordinate rank must
    /// match and remain nonsingular; unselected modes remain external in the
    /// order Helmholtz then LNS, each preserving its original mode order.
    /// Only geometry-identical faces are admitted; nonmatching mortar meshes
    /// and approximate geometric welding are not silently inferred.
    pub fn new(
        helmholtz: &'a InteriorMultiport,
        thermoviscous: &'a ThermoviscousFe<'c>,
        helmholtz_modes: &[usize],
        thermoviscous_modes: &[usize],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Self, HybridError> {
        poll(cancelled)?;
        let m = helmholtz_modes.len();
        let nh = helmholtz.port_mode_count();
        let nl = thermoviscous.port_mode_count();
        let valid = |ids: &[usize], cap: usize| {
            ids.iter().all(|&i| i < cap)
                && ids.iter().copied().collect::<BTreeSet<_>>().len() == ids.len()
        };
        if m == 0
            || thermoviscous_modes.len() != m
            || !valid(helmholtz_modes, nh)
            || !valid(thermoviscous_modes, nl)
            || 2 * m == nh + nl
        {
            return Err(HybridError::InvalidInput(
                "interface indices/rank or no external ports",
            ));
        }
        let hfaces = collect_faces(&helmholtz.physical_modes, helmholtz_modes);
        let lfaces = collect_faces(&thermoviscous.physical_modes, thermoviscous_modes);
        if hfaces.keys().ne(lfaces.keys()) {
            return Err(HybridError::InvalidInput(
                "interface triangulations are not coincident",
            ));
        }
        let mut flow_map = vec![0.0; m * m];
        for (face, hs) in &hfaces {
            poll(cancelled)?;
            for &(h, ht) in hs {
                for &(l, lt) in &lfaces[face] {
                    let normal_dot: f64 = ht
                        .outward_normal
                        .iter()
                        .zip(lt.outward_normal)
                        .map(|(a, b)| a * b)
                        .sum();
                    if (normal_dot + 1.0).abs() > 1e-12
                        || (ht.area_m2 - lt.area_m2).abs() > 1e-12 * ht.area_m2
                    {
                        return Err(HybridError::InvalidInput(
                            "interface outward normals/areas disagree",
                        ));
                    }
                    let moment = (ht.values.iter().sum::<f64>() * lt.values.iter().sum::<f64>()
                        + ht.values
                            .iter()
                            .zip(lt.values)
                            .map(|(a, b)| a * b)
                            .sum::<f64>())
                        * ht.area_m2
                        / 12.0;
                    flow_map[l * m + h] += moment / helmholtz.port_areas_m2()[helmholtz_modes[h]];
                }
            }
        }
        if flow_map.iter().any(|x| !x.is_finite()) {
            return Err(HybridError::NonFinite);
        }
        // This only admits the retained trace rank; it never inverts the map
        // in the physical solve or replaces zero/evanescent modes by a cutoff.
        let mut scaled = flow_map.clone();
        for row in 0..m {
            let scale = scaled[row * m..(row + 1) * m]
                .iter()
                .map(|x| x.abs())
                .fold(0.0_f64, f64::max);
            if scale == 0.0 {
                return Err(HybridError::InvalidInput(
                    "unresolved interface trace pairing",
                ));
            }
            for value in &mut scaled[row * m..(row + 1) * m] {
                *value /= scale;
            }
        }
        lu_complex_with_cancel(
            &scaled.iter().copied().map(C64::from_re).collect::<Vec<_>>(),
            m,
            cancelled,
        )?;
        let internal = helmholtz_modes
            .iter()
            .copied()
            .chain(thermoviscous_modes.iter().map(|&p| nh + p))
            .collect::<Vec<_>>();
        let external = (0..nh + nl)
            .filter(|p| !internal.contains(p))
            .collect::<Vec<_>>();
        poll(cancelled)?;
        Ok(Self {
            helmholtz,
            thermoviscous,
            interface: ModalInterface {
                helmholtz_modes: helmholtz_modes.to_vec(),
                thermoviscous_modes: thermoviscous_modes.to_vec(),
                flow_map,
            },
            internal,
            external,
        })
    }

    /// Actual source-bound cross-trace pairing.
    #[must_use]
    pub fn interface(&self) -> &ModalInterface {
        &self.interface
    }
    /// Original equilibrium pressure-region medium.
    #[must_use]
    pub fn helmholtz_medium(&self) -> Medium {
        self.helmholtz.medium()
    }
    /// Original resolved-region equilibrium/transport.
    #[must_use]
    pub fn thermoviscous_medium(&self) -> ThermoviscousMedium {
        self.thermoviscous.medium()
    }
    /// Number of unconnected coordinates, with no automatic truncation.
    #[must_use]
    pub fn external_port_count(&self) -> usize {
        self.external.len()
    }

    /// Retain each loaded volume factor and one small interface factor. cp
    /// modifies ONLY the LNS owner; no molecular dispersion is silently added
    /// to the Helmholtz region. Use equilibrium real cp for a common ideal gas.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        cp: C64,
        helmholtz_reference: &[f64],
        thermoviscous_reference: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedHybrid<'_, 'c>, HybridError> {
        poll(cancelled)?;
        let h = self
            .helmholtz
            .prepare_scattering(frequency_hz, helmholtz_reference, cancelled)?;
        self.prepare_loaded(frequency_hz, cp, h, thermoviscous_reference, cancelled)
    }

    /// Use source-evaluated complex cp in LNS and bulk modulus K in Helmholtz,
    /// both in exp(+iwt). This carries caloric relaxation through both regions
    /// without adding a second compressibility loss to the resolved LNS owner.
    /// Helmholtz boundary-layer overlays refuse; resolved LNS conductive solids
    /// remain supported by its original joint velocity/pressure/heat equations.
    pub fn prepare_with_bulk_modulus(
        &self,
        frequency_hz: f64,
        cp: C64,
        bulk_modulus_pa: C64,
        helmholtz_reference: &[f64],
        thermoviscous_reference: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedHybrid<'_, 'c>, HybridError> {
        poll(cancelled)?;
        let h = self.helmholtz.prepare_scattering_with_bulk_modulus(
            frequency_hz,
            bulk_modulus_pa,
            helmholtz_reference,
            cancelled,
        )?;
        self.prepare_loaded(frequency_hz, cp, h, thermoviscous_reference, cancelled)
    }

    fn prepare_loaded<'p>(
        &'p self,
        frequency_hz: f64,
        cp: C64,
        h: PreparedInteriorScattering<'p>,
        thermoviscous_reference: &[f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<PreparedHybrid<'p, 'c>, HybridError> {
        let helmholtz_reference = h.reference_resistance();
        let l = self.thermoviscous.prepare_scattering(
            frequency_hz,
            cp,
            thermoviscous_reference,
            &mut || cancelled(),
        )?;
        let sh = h.scattering(cancelled)?;
        let sl = l.scattering(&mut || cancelled())?;
        let nh = helmholtz_reference.len();
        let nl = thermoviscous_reference.len();
        let n = nh + nl;
        let mut s = vec![C64::ZERO; n * n];
        for i in 0..nh {
            s[i * n..i * n + nh].copy_from_slice(&sh[i * nh..(i + 1) * nh]);
        }
        for i in 0..nl {
            s[(nh + i) * n + nh..(nh + i + 1) * n].copy_from_slice(&sl[i * nl..(i + 1) * nl]);
        }
        let reference = helmholtz_reference
            .iter()
            .chain(thermoviscous_reference)
            .copied()
            .collect::<Vec<_>>();
        let join = junction_scattering(&self.interface, &reference, nh, cancelled)?;
        let ni = self.internal.len();
        let mut a = vec![C64::ZERO; ni * ni];
        for i in 0..ni {
            poll(cancelled)?;
            for j in 0..ni {
                let product = (0..ni).fold(C64::ZERO, |sum, k| {
                    sum + s[self.internal[k] * n + self.internal[j]].scale(join[i * ni + k])
                });
                a[i * ni + j] = C64::from_re(if i == j { 1.0 } else { 0.0 }) - product;
            }
        }
        if a.iter().any(|&z| !finite(z)) {
            return Err(HybridError::NonFinite);
        }
        let factor = lu_complex_with_cancel(&a, ni, cancelled)?;
        poll(cancelled)?;
        Ok(PreparedHybrid {
            h,
            l,
            h_owner: self.helmholtz,
            frequency_hz,
            interface: self.interface.clone(),
            internal: self.internal.clone(),
            external: self.external.clone(),
            reference,
            s,
            join,
            factor,
            nh,
        })
    }
}

// Real lossless transformer in power coordinates. K=sqrt(R_L) T/sqrt(R_H),
// (I+K^T K)a_H=(K^T K-I)b_H+2K^T b_L, a_L=b_L-K(a_H-b_H).
fn junction_scattering(
    interface: &ModalInterface,
    reference: &[f64],
    nh: usize,
    cancelled: &dyn Fn() -> bool,
) -> Result<Vec<f64>, HybridError> {
    let m = interface.helmholtz_modes.len();
    let mut k = vec![0.0; m * m];
    for l in 0..m {
        for h in 0..m {
            k[l * m + h] = interface.flow_map[l * m + h]
                * det::sqrt(reference[nh + interface.thermoviscous_modes[l]])
                / det::sqrt(reference[interface.helmholtz_modes[h]]);
        }
    }
    let mut a = vec![C64::ZERO; m * m];
    for i in 0..m {
        for j in 0..m {
            a[i * m + j] = C64::from_re(
                (if i == j { 1.0 } else { 0.0 })
                    + (0..m).map(|l| k[l * m + i] * k[l * m + j]).sum::<f64>(),
            );
        }
    }
    if a.iter().any(|&z| !finite(z)) {
        return Err(HybridError::NonFinite);
    }
    let factor = lu_complex_with_cancel(&a, m, cancelled)?;
    let mut join = vec![0.0; 4 * m * m];
    for col in 0..2 * m {
        poll(cancelled)?;
        let mut rhs = (0..m)
            .map(|i| {
                if col < m {
                    a[i * m + col] - C64::from_re(if i == col { 2.0 } else { 0.0 })
                } else {
                    C64::from_re(2.0 * k[(col - m) * m + i])
                }
            })
            .collect::<Vec<_>>();
        factor.solve_with_cancel(&mut rhs, cancelled)?;
        for row in 0..m {
            join[row * 2 * m + col] = rhs[row].re;
            join[(m + row) * 2 * m + col] = (if col == m + row { 1.0 } else { 0.0 })
                + (if col < m { k[row * m + col] } else { 0.0 })
                - (0..m).map(|j| k[row * m + j] * rhs[j].re).sum::<f64>();
        }
    }
    if join.iter().any(|x| !x.is_finite()) {
        return Err(HybridError::NonFinite);
    }
    Ok(join)
}

/// All per-frequency factorizations, reusable for scattering and actual fields.
pub struct PreparedHybrid<'a, 'c> {
    h: PreparedInteriorScattering<'a>,
    l: PreparedThermoviscousScattering<'a, 'c>,
    h_owner: &'a InteriorMultiport,
    frequency_hz: f64,
    interface: ModalInterface,
    internal: Vec<usize>,
    external: Vec<usize>,
    reference: Vec<f64>,
    s: Vec<C64>,
    join: Vec<f64>,
    factor: LuComplex,
    nh: usize,
}

/// Independent acoustic work on both sides and physical losses [W].
#[derive(Debug, Clone, Copy)]
pub struct HybridPower {
    pub external_input_w: f64,
    pub helmholtz_interface_w: f64,
    pub thermoviscous_interface_w: f64,
    pub helmholtz_wall: BoundaryLayerPower,
    /// Molecular compression work in the original Helmholtz volume [W].
    pub helmholtz_bulk_w: f64,
    pub thermoviscous: ThermoviscousPower,
}
impl HybridPower {
    /// Signed interface work defect, with inward orientation on both owners.
    #[must_use]
    pub fn interface_residual_w(self) -> f64 {
        self.helmholtz_interface_w + self.thermoviscous_interface_w
    }
    /// Physical volume/wall dissipation, without clipping or a fitted loss term.
    #[must_use]
    pub fn dissipated_w(self) -> f64 {
        self.helmholtz_wall.viscous_w
            + self.helmholtz_wall.thermal_w
            + self.helmholtz_bulk_w
            + self.thermoviscous.dissipated_w()
    }
    #[must_use]
    pub fn balance_residual_w(self) -> f64 {
        self.external_input_w - self.dissipated_w()
    }
}

/// Original fields plus external coordinates and explicit interface defects.
#[derive(Debug, Clone)]
pub struct HybridResponse {
    pub helmholtz: InteriorPortResponse,
    pub thermoviscous: ThermoviscousField,
    pub pressure_pa: Vec<C64>,
    pub inward_flow_m3_s: Vec<C64>,
    pub outgoing: Vec<C64>,
    /// P_H-T^T P_L, in selected Helmholtz mode order [Pa].
    pub interface_pressure_residual_pa: Vec<C64>,
    /// Q_L+T Q_H, in selected LNS mode order [m^3/s].
    pub interface_flow_residual_m3_s: Vec<C64>,
    pub power: HybridPower,
}

impl PreparedHybrid<'_, '_> {
    /// Positive real references on external coordinates, in output order.
    #[must_use]
    pub fn external_reference_resistance(&self) -> Vec<f64> {
        self.external.iter().map(|&p| self.reference[p]).collect()
    }

    fn complete_incident(
        &self,
        incident: &[C64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<C64>, HybridError> {
        poll(cancelled)?;
        if incident.len() != self.external.len() || incident.iter().any(|&z| !finite(z)) {
            return Err(HybridError::InvalidInput(
                "external incident wave coordinates",
            ));
        }
        let n = self.reference.len();
        let ni = self.internal.len();
        let b = self
            .internal
            .iter()
            .map(|&row| {
                self.external
                    .iter()
                    .zip(incident)
                    .fold(C64::ZERO, |sum, (&col, &a)| sum + self.s[row * n + col] * a)
            })
            .collect::<Vec<_>>();
        let mut a = (0..ni)
            .map(|i| (0..ni).fold(C64::ZERO, |sum, j| sum + b[j].scale(self.join[i * ni + j])))
            .collect::<Vec<_>>();
        self.factor.solve_with_cancel(&mut a, cancelled)?;
        let mut all = vec![C64::ZERO; n];
        for (&p, &v) in self
            .external
            .iter()
            .zip(incident)
            .chain(self.internal.iter().zip(&a))
        {
            all[p] = v;
        }
        if all.iter().any(|&z| !finite(z)) {
            return Err(HybridError::NonFinite);
        }
        poll(cancelled)?;
        Ok(all)
    }

    /// Full external power scattering, retaining all mode mixing and using
    /// only the already prepared volume and interface factors.
    pub fn scattering(&self, cancelled: &dyn Fn() -> bool) -> Result<Vec<C64>, HybridError> {
        let m = self.external.len();
        let n = self.reference.len();
        let mut s = vec![C64::ZERO; m * m];
        let mut a = vec![C64::ZERO; m];
        for col in 0..m {
            a[col] = C64::ONE;
            let all = self.complete_incident(&a, cancelled)?;
            a[col] = C64::ZERO;
            for (row, &p) in self.external.iter().enumerate() {
                s[row * m + col] =
                    (0..n).fold(C64::ZERO, |sum, j| sum + self.s[p * n + j] * all[j]);
            }
        }
        if s.iter().any(|&z| !finite(z)) {
            return Err(HybridError::NonFinite);
        }
        poll(cancelled)?;
        Ok(s)
    }

    /// Resolve actual original fields for simultaneous incident waves. The
    /// work defect uses independently recovered P/Q and original volume/surface
    /// losses, rather than defining dissipation as the scattering norm defect.
    pub fn response(
        &self,
        incident: &[C64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<HybridResponse, HybridError> {
        let all = self.complete_incident(incident, cancelled)?;
        let h = self.h.response(&all[..self.nh], cancelled)?;
        let l = self.l.response(&all[self.nh..], &mut || cancelled())?;
        let pressure = h
            .pressure_pa
            .iter()
            .chain(&l.pressure_drive_pa)
            .copied()
            .collect::<Vec<_>>();
        let flow = h
            .inward_flow_m3_s
            .iter()
            .chain(&l.inward_flow_m3_s)
            .copied()
            .collect::<Vec<_>>();
        let pressure_pa = self
            .external
            .iter()
            .map(|&p| pressure[p])
            .collect::<Vec<_>>();
        let inward_flow_m3_s = self.external.iter().map(|&p| flow[p]).collect::<Vec<_>>();
        let outgoing = self
            .external
            .iter()
            .map(|&p| {
                (pressure[p] - flow[p].scale(self.reference[p]))
                    .scale(0.5 / det::sqrt(self.reference[p]))
            })
            .collect::<Vec<_>>();
        let work = |ids: &[usize]| {
            ids.iter()
                .map(|&p| 0.5 * (pressure[p].conj() * flow[p]).re)
                .sum::<f64>()
        };
        let m = self.interface.helmholtz_modes.len();
        let power = HybridPower {
            external_input_w: work(&self.external),
            helmholtz_interface_w: work(&self.internal[..m]),
            thermoviscous_interface_w: work(&self.internal[m..]),
            helmholtz_wall: self
                .h_owner
                .boundary_layer_power(self.frequency_hz, &h.amplitudes)?,
            helmholtz_bulk_w: self.h.bulk_dissipation_w(&h.amplitudes, cancelled)?,
            thermoviscous: l.power,
        };
        let interface_pressure_residual_pa = (0..m)
            .map(|i| {
                pressure[self.internal[i]]
                    - (0..m).fold(C64::ZERO, |sum, j| {
                        sum + pressure[self.internal[m + j]]
                            .scale(self.interface.flow_map[j * m + i])
                    })
            })
            .collect::<Vec<_>>();
        let interface_flow_residual_m3_s = (0..m)
            .map(|i| {
                flow[self.internal[m + i]]
                    + (0..m).fold(C64::ZERO, |sum, j| {
                        sum + flow[self.internal[j]].scale(self.interface.flow_map[i * m + j])
                    })
            })
            .collect::<Vec<_>>();
        if outgoing
            .iter()
            .chain(&interface_pressure_residual_pa)
            .chain(&interface_flow_residual_m3_s)
            .any(|&z| !finite(z))
            || !power.external_input_w.is_finite()
            || !power.interface_residual_w().is_finite()
        {
            return Err(HybridError::NonFinite);
        }
        poll(cancelled)?;
        Ok(HybridResponse {
            helmholtz: h,
            thermoviscous: l,
            pressure_pa,
            inward_flow_m3_s,
            outgoing,
            interface_pressure_residual_pa,
            interface_flow_residual_m3_s,
            power,
        })
    }
}

#[cfg(test)]
mod tests;
