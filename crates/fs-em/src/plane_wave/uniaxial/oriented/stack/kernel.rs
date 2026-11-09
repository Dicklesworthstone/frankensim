//! Stable layer elimination; no inverse evanescent propagation is formed.
use super::{Frequency, ResolvedLayer, Result, bad, poll};
use super::super::kernel::{FieldSum, Input as ModeInput, WaveGeometry};
use super::super::{Incidence, slab::SlabWave};
use super::super::super::kernel::UniaxialWaveKind;
use super::super::slab::kernel::{
    combine, crystal_basis, exit_basis_from_incidence, finite, flux, positive, solve_boundary, tangential,
};
use super::super::spectrum::{UniaxialSpectrumField, rotate};
use crate::plane_wave::angular_spectrum::modes::propagation;
use fs_math::{c64::C64, dd::Dd};

const C0: f64 = 299_792_458.0;
type Vector = [C64; 3];
type Matrix = [[C64; 2]; 2];
type Trace = [[C64; 4]; 2]; // columns are modes, rows are Ex,Ey,Hx,Hy

/// Physical point selection. Every coordinate is global; both adjacent layers
/// admit a shared face, allowing their distinct normal E/H to be inspected.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StackRegion { Incident, Reflected, TotalIncidentSide, Layer(usize), Transmitted }

/// Four actual eigenmode fields with their original decaying-side anchors.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialLayerField {
    pub entrance_z_m: f64,
    pub exit_z_m: f64,
    pub thickness_m: f64,
    pub forward: [SlabWave; 2],
    pub backward: [SlabWave; 2],
}

/// Complete coherent scattering snapshot, retaining every internal layer.
/// Transmitted fields are referenced at the actual final face.
#[derive(Clone, Debug, PartialEq)]
pub struct UniaxialStackResponse {
    pub thickness_m: f64,
    pub tangential_wavenumber_rad_m: f64,
    pub incident_normal_wavenumber_rad_m: f64,
    pub exit_normal_wavenumber_rad_m: C64,
    pub exit_wave_kind: UniaxialWaveKind,
    pub incident_electric_v_m: Vector,
    pub incident_magnetic_a_m: Vector,
    pub reflected_electric_v_m: Vector,
    pub reflected_magnetic_a_m: Vector,
    pub transmitted_electric_v_m: Vector,
    pub transmitted_magnetic_a_m: Vector,
    pub layers: Vec<UniaxialLayerField>,
    pub incident_normal_power_w_m2: f64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
}

impl UniaxialLayerField {
    fn field_sum(&self, point: [f64; 3]) -> Result<FieldSum> {
        if !(self.entrance_z_m..=self.exit_z_m).contains(&point[2]) {
            return Err(bad("uniaxial stack point is outside its selected layer"));
        }
        // Preserve original thickness exactly at each authored face, even when
        // cumulative global coordinates incur one rounding at layer addition.
        // Elsewhere retain the global subtraction residual in the local depth.
        let depth = if point[2] == self.exit_z_m {
            Dd::from_f64(self.thickness_m)
        } else {
            Dd::from_f64(point[2]) - Dd::from_f64(self.entrance_z_m)
        };
        FieldSum::at_local_point([point[0], point[1]], depth).map_err(bad)
    }
    fn add(&self, sum: &mut FieldSum, tangential: f64, axis: [f64; 2]) -> Result<()> {
        for (waves, anchor) in [
            (&self.forward, 0.0),
            (&self.backward, self.thickness_m),
        ] {
            for wave in waves {
                // Keep the signed GLOBAL root. The backward distance is
                // depth-thickness, so its evanescence decays toward the front.
                sum.add_wave(&wave.electric_v_m, &wave.magnetic_a_m, WaveGeometry {
                    tangential_wavenumber_rad_m: tangential,
                    incidence_plane_axis: axis,
                    normal_wavenumber_rad_m: wave.normal_wavenumber_rad_m,
                    anchor_z_m: anchor,
                }).map_err(bad)?;
            }
        }
        Ok(())
    }
    fn face(&self, end: bool) -> Result<UniaxialSpectrumField> {
        let z = if end { self.exit_z_m } else { self.entrance_z_m };
        let mut sum = self.field_sum([0.0, 0.0, z])?;
        self.add(&mut sum, 0.0, [1.0, 0.0])?;
        finish_field(sum)
    }
}

impl UniaxialStackResponse {
    /// All Cartesian E/H and cycle-mean Poynting at [x,y,z] [m], with no
    /// allocation, material reevaluation or boundary solve. Relative phases,
    /// face distances and coherent products remain compensated until output.
    pub fn field_at(&self, region: StackRegion, point: [f64; 3]) -> Result<UniaxialSpectrumField> {
        let mut sum = self.field_sum(region, point)?;
        self.accumulate(&mut sum, region, [1.0,0.0])?;
        finish_field(sum)
    }
    pub(super) fn field_sum(&self, region: StackRegion, point: [f64; 3]) -> Result<FieldSum> {
        if point.iter().any(|v| !v.is_finite()) {
            return Err(bad("finite uniaxial stack field point required"));
        }
        if let StackRegion::Layer(index) = region {
            return self.layers.get(index).ok_or_else(|| bad("unknown uniaxial stack layer"))?
                .field_sum(point);
        }
        if if region == StackRegion::Transmitted { point[2] < self.thickness_m } else { point[2] > 0.0 } {
            return Err(bad("uniaxial stack point is outside its selected exterior"));
        }
        FieldSum::at_point(point).map_err(bad)
    }
    /// The caller admits one point before summing modes of this SAME stack.
    /// Interior waves share the selected layer's compensated local coordinate.
    pub(super) fn accumulate(
        &self, sum: &mut FieldSum, region: StackRegion, axis: [f64; 2],
    ) -> Result<()> {
        if let StackRegion::Layer(index) = region {
            return self.layers.get(index).ok_or_else(|| bad("unknown uniaxial stack layer"))?
                .add(sum, self.tangential_wavenumber_rad_m, axis);
        }
        let mut add = |e: &Vector, h: &Vector, k, anchor| -> Result<()> {
            sum.add_wave(e, h, WaveGeometry {
                tangential_wavenumber_rad_m: self.tangential_wavenumber_rad_m,
                incidence_plane_axis: axis,
                normal_wavenumber_rad_m: k,
                anchor_z_m: anchor,
            }).map_err(bad)
        };
        if matches!(region, StackRegion::Incident | StackRegion::TotalIncidentSide) {
            add(&self.incident_electric_v_m, &self.incident_magnetic_a_m,
                C64::from_re(self.incident_normal_wavenumber_rad_m), 0.0)?;
        }
        if matches!(region, StackRegion::Reflected | StackRegion::TotalIncidentSide) {
            add(&self.reflected_electric_v_m, &self.reflected_magnetic_a_m,
                C64::from_re(-self.incident_normal_wavenumber_rad_m), 0.0)?;
        }
        if region == StackRegion::Transmitted {
            add(&self.transmitted_electric_v_m, &self.transmitted_magnetic_a_m,
                self.exit_normal_wavenumber_rad_m, self.thickness_m)?;
        }
        Ok(())
    }
    pub(super) fn rotate(&mut self, axis: [f64; 2]) {
        for vector in [
            &mut self.incident_electric_v_m, &mut self.incident_magnetic_a_m,
            &mut self.reflected_electric_v_m, &mut self.reflected_magnetic_a_m,
            &mut self.transmitted_electric_v_m, &mut self.transmitted_magnetic_a_m,
        ] { *vector = rotate(*vector, axis); }
        for layer in &mut self.layers {
            for wave in layer.forward.iter_mut().chain(&mut layer.backward) {
                wave.electric_v_m = rotate(wave.electric_v_m, axis);
                wave.magnetic_a_m = rotate(wave.magnetic_a_m, axis);
            }
        }
    }
}

pub(super) fn finish_field(sum: FieldSum) -> Result<UniaxialSpectrumField> {
    let field = sum.finish().map_err(bad)?;
    Ok(UniaxialSpectrumField {
        electric_v_m: field.electric_v_m,
        magnetic_a_m: field.magnetic_a_m,
        poynting_w_m2: field.poynting_w_m2,
    })
}

struct LayerWork {
    field: UniaxialLayerField,
    forward_phase: [C64; 2],
    backward_phase: [C64; 2],
    reflection: Matrix,
    through: Matrix,
}

pub(super) fn workspace_bytes(layers: usize, modes: usize) -> Result<usize> {
    let retained = layers.checked_mul(core::mem::size_of::<UniaxialLayerField>())
        .and_then(|n| n.checked_add(core::mem::size_of::<UniaxialStackResponse>()))
        .and_then(|n| n.checked_add(core::mem::size_of::<[f64;2]>()))
        .and_then(|n| n.checked_mul(modes));
    let scratch = layers.checked_mul(core::mem::size_of::<LayerWork>() + core::mem::size_of::<ResolvedLayer>());
    retained.and_then(|n| n.checked_add(scratch?))
        .and_then(|n| n.checked_add(core::mem::size_of::<FieldSum>() + core::mem::size_of::<UniaxialSpectrumField>()))
        .ok_or_else(|| bad("uniaxial stack workspace size overflow"))
}

fn apply(a: Matrix, b: [C64;2]) -> [C64;2] {
    a.map(|row| row[0]*b[0]+row[1]*b[1])
}
fn trace(waves: [SlabWave;2]) -> Trace {
    waves.map(|wave| tangential(wave.electric_v_m,wave.magnetic_a_m))
}
fn interface(back: [SlabWave;2], next: Trace, forward: [SlabWave;2]) -> Result<(Matrix,Matrix)> {
    let back=trace(back);
    let matrix=core::array::from_fn(|row| [back[0][row],back[1][row],-next[0][row],-next[1][row]]);
    let mut reflection=[[C64::ZERO;2];2];
    let mut through=reflection;
    for (column,rhs) in trace(forward).into_iter().enumerate() {
        let x=solve_boundary(matrix,rhs.map(|v| -v)).map_err(bad)?;
        for row in 0..2 { reflection[row][column]=x[row]; through[row][column]=x[row+2]; }
    }
    Ok((reflection,through))
}
fn front(work: &LayerWork) -> Trace {
    let forward=trace(work.field.forward);
    let backward=trace(work.field.backward);
    // Noncommuting order is B Q R P: backward right-anchor load R acts on
    // forward fields arriving through P, then Q returns them to the left.
    core::array::from_fn(|column| core::array::from_fn(|row| {
        forward[column][row]
            + backward[0][row]*work.backward_phase[0]*work.reflection[0][column]*work.forward_phase[column]
            + backward[1][row]*work.backward_phase[1]*work.reflection[1][column]*work.forward_phase[column]
    }))
}
fn weight(waves: [SlabWave;2], coefficients: [C64;2]) -> [SlabWave;2] {
    core::array::from_fn(|i| SlabWave {
        electric_v_m: waves[i].electric_v_m.map(|z| z*coefficients[i]),
        magnetic_a_m: waves[i].magnetic_a_m.map(|z| z*coefficients[i]),
        ..waves[i]
    })
}

pub(super) fn solve(
    data: &Frequency, incidence: Incidence, electric: [C64;2], plane: [f64;2],
    cancelled: &mut impl FnMut()->bool,
) -> Result<UniaxialStackResponse> {
    poll(cancelled)?;
    if data.layers.is_empty() || electric.iter().any(|&z| !finite(z)) {
        return Err(bad("finite nonzero stack illumination and layers required"));
    }
    let amplitude=positive(electric.iter().fold(0.0_f64,|m,z| m.max(z.re.abs()).max(z.im.abs()))).map_err(bad)?;
    let electric=electric.map(|z| C64::new(z.re/amplitude,z.im/amplitude));
    let mut work=Vec::new();
    work.try_reserve_exact(data.layers.len()).map_err(|_| bad("uniaxial stack workspace allocation failed"))?;
    let mut first=None;
    for layer in &data.layers {
        poll(cancelled)?;
        let [x,y,z]=layer.axis;
        let (basis,forward,backward)=crystal_basis(ModeInput {
            omega:data.omega,ni:data.ni,no:layer.no,ne:layer.ne,mui:data.mui,mut_:layer.mu,
            axis:[x*plane[0]+y*plane[1],y*plane[0]-x*plane[1],z],
            // Every layer retains the SAME original incident index and
            // coordinate. Replacing cosine by rounded p loses shallow flux.
            incidence,
            electric,
        }).map_err(bad)?;
        if first.is_none() { first=Some(basis); }
        let mut phases=[[C64::ZERO;2];2];
        for j in 0..2 {
            phases[0][j]=propagation(forward[j].normal_wavenumber_rad_m,layer.thickness)
                .ok_or_else(|| bad("uniaxial stack forward phase is unrepresentable"))?;
            phases[1][j]=propagation(-backward[j].normal_wavenumber_rad_m,layer.thickness)
                .ok_or_else(|| bad("uniaxial stack backward phase is unrepresentable"))?;
        }
        work.push(LayerWork {
            field:UniaxialLayerField {entrance_z_m:layer.entrance,exit_z_m:layer.exit,
                thickness_m:layer.thickness,forward,backward},
            forward_phase:phases[0],backward_phase:phases[1],
            reflection:[[C64::ZERO;2];2],through:[[C64::ZERO;2];2],
        });
    }
    let basis=first.ok_or_else(|| bad("empty uniaxial stack"))?;
    let kscale=positive((data.omega/C0)*basis.scale).map_err(bad)?;
    let hscale=positive(basis.scale/(data.mui*C0)).map_err(bad)?;
    let (p,q,n)=(basis.p,basis.qi,basis.i);
    let normal=positive(kscale*q).map_err(bad)?;
    let tangential_k=kscale*p;
    if !tangential_k.is_finite() || p!=0.0 && tangential_k==0.0 {
        return Err(bad("uniaxial stack tangential wavenumber is unrepresentable"));
    }
    let incident_e=[electric[0],electric[1],-electric[0].scale(p)/C64::from_re(q)];
    let incident_h=[-incident_e[1].scale(q*hscale),
        (incident_e[0].scale(q)-incident_e[2].scale(p)).scale(hscale),incident_e[1].scale(p*hscale)];
    let incoming=positive(flux(incident_e,incident_h)).map_err(bad)?;
    let reflected=[
        SlabWave {wave_kind:UniaxialWaveKind::Propagating,normal_wavenumber_rad_m:C64::from_re(-normal),
            electric_v_m:[C64::ZERO,C64::ONE,C64::ZERO],
            magnetic_a_m:[C64::from_re(q*hscale),C64::ZERO,C64::from_re(p*hscale)]},
        SlabWave {wave_kind:UniaxialWaveKind::Propagating,normal_wavenumber_rad_m:C64::from_re(-normal),
            electric_v_m:[C64::from_re(q/n),C64::ZERO,C64::from_re(p/n)],
            magnetic_a_m:[C64::ZERO,C64::from_re(-n*hscale),C64::ZERO]},
    ];
    let exit=exit_basis_from_incidence(data.omega,data.ni,incidence,data.nt,data.mut_).map_err(bad)?;
    let mut next=trace(exit);
    for layer in work.iter_mut().rev() {
        poll(cancelled)?;
        (layer.reflection,layer.through)=interface(layer.field.backward,next,layer.field.forward)?;
        next=front(layer);
    }
    let r=trace(reflected);
    let entrance=core::array::from_fn(|row| [r[0][row],r[1][row],-next[0][row],-next[1][row]]);
    let x=solve_boundary(entrance,tangential(incident_e,incident_h).map(|z| -z)).map_err(bad)?;
    let (re,rh)=combine(reflected,[x[0],x[1]]);
    let mut forward=[x[2],x[3]];
    let mut layers=Vec::new();
    layers.try_reserve_exact(work.len()).map_err(|_| bad("uniaxial stack field allocation failed"))?;
    for layer in work {
        poll(cancelled)?;
        let arriving=core::array::from_fn(|i| forward[i]*layer.forward_phase[i]);
        let back=apply(layer.reflection,arriving);
        layers.push(UniaxialLayerField {forward:weight(layer.field.forward,forward),
            backward:weight(layer.field.backward,back),..layer.field});
        forward=apply(layer.through,arriving);
    }
    let (te,th)=combine(exit,forward);
    let reflected_power=-flux(re,rh)/incoming;
    let transmission=flux(te,th)/incoming;
    let transmitted_power=if exit[0].wave_kind==UniaxialWaveKind::Propagating {
        positive(transmission).map_err(bad)?
    } else {
        if !transmission.is_finite() || transmission.abs()>4096.0*f64::EPSILON {
            return Err(bad("uniaxial stack nonpropagating exit has unresolved normal power"));
        }
        0.0
    };
    if !reflected_power.is_finite() || reflected_power<0.0
        || (reflected_power+transmitted_power-1.0).abs()>4096.0*f64::EPSILON {
        return Err(bad("uniaxial stack Maxwell normal-power balance failed"));
    }
    let mut result=UniaxialStackResponse {
        thickness_m:data.layers.last().unwrap().exit,tangential_wavenumber_rad_m:tangential_k,
        incident_normal_wavenumber_rad_m:normal,exit_normal_wavenumber_rad_m:exit[0].normal_wavenumber_rad_m,
        exit_wave_kind:exit[0].wave_kind,incident_electric_v_m:incident_e,incident_magnetic_a_m:incident_h,
        reflected_electric_v_m:re,reflected_magnetic_a_m:rh,transmitted_electric_v_m:te,transmitted_magnetic_a_m:th,
        layers,incident_normal_power_w_m2:incoming,reflected_power_fraction:reflected_power,
        transmitted_power_fraction:transmitted_power,
    };
    boundary_gate(&result,positive(n*hscale).map_err(bad)?,cancelled)?;
    for vector in [&mut result.incident_electric_v_m,&mut result.incident_magnetic_a_m,
        &mut result.reflected_electric_v_m,&mut result.reflected_magnetic_a_m,
        &mut result.transmitted_electric_v_m,&mut result.transmitted_magnetic_a_m] {
        *vector=vector.map(|z| z.scale(amplitude));
        if vector.iter().any(|&z| !finite(z)) {return Err(bad("uniaxial stack dimensional fields are unrepresentable"));}
    }
    for layer in &mut result.layers {
        poll(cancelled)?;
        for wave in layer.forward.iter_mut().chain(&mut layer.backward) {
            wave.electric_v_m=wave.electric_v_m.map(|z| z.scale(amplitude));
            wave.magnetic_a_m=wave.magnetic_a_m.map(|z| z.scale(amplitude));
            if wave.electric_v_m.iter().chain(&wave.magnetic_a_m).any(|&z| !finite(z)) {
                return Err(bad("uniaxial stack dimensional fields are unrepresentable"));
            }
        }
    }
    result.incident_normal_power_w_m2=positive((incoming*amplitude)*amplitude).map_err(bad)?;
    poll(cancelled)?;
    Ok(result)
}

fn boundary_gate(r: &UniaxialStackResponse, admittance: f64, cancelled: &mut impl FnMut()->bool) -> Result<()> {
    let mut left=r.field_at(StackRegion::TotalIncidentSide,[0.0;3])?;
    let continuity=|a: &UniaxialSpectrumField,b: &UniaxialSpectrumField| -> Result<()> {
        let es=[a,b].iter().fold(1.0_f64,|m,f| m.max(f.electric_v_m[0].abs()).max(f.electric_v_m[1].abs()));
        let hs=[a,b].iter().fold(admittance,|m,f| m.max(f.magnetic_a_m[0].abs()).max(f.magnetic_a_m[1].abs()));
        for j in 0..2 {
            if (a.electric_v_m[j]-b.electric_v_m[j]).abs()>4096.0*f64::EPSILON*es
                || (a.magnetic_a_m[j]-b.magnetic_a_m[j]).abs()>4096.0*f64::EPSILON*hs {
                return Err(bad("uniaxial stack tangential Maxwell continuity is unresolved"));
            }
        }
        Ok(())
    };
    for layer in &r.layers {
        poll(cancelled)?;
        let front=layer.face(false)?;
        let back=layer.face(true)?;
        continuity(&left,&front)?;
        for field in [&front,&back] {
            let ratio=field.poynting_w_m2[2]/r.incident_normal_power_w_m2;
            if !ratio.is_finite() || (ratio-r.transmitted_power_fraction).abs()>8192.0*f64::EPSILON {
                return Err(bad("uniaxial stack internal normal-power balance is unresolved"));
            }
        }
        left=back;
    }
    continuity(&left,&r.field_at(StackRegion::Transmitted,[0.0,0.0,r.thickness_m])?)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod anchored_tests;
