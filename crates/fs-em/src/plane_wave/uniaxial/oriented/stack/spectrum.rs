//! Original sparse illumination through the complete coupled crystal stack.
use super::{Frequency, OrientedUniaxialStack, Result, StackRegion, UniaxialStackResponse, bad, kernel, poll};
use super::super::{Incidence, spectrum::{
    FieldSum, UniaxialIncidentMode, UniaxialSpectrumField, UniaxialSpectrumLimits,
}};
use fs_math::c64::C64;

struct Mode { axis: [f64;2], response: UniaxialStackResponse }

/// Original material laws are resolved once per preparation. Each nonzero
/// authored direction is then solved against all ORIGINAL fixed global layer
/// axes. All four waves of the selected layer and every incident mode are
/// coherently summed before the physical Poynting vector is formed.
pub struct PreparedUniaxialStackSpectrum<'a,'s,'m> {
    source: &'a OrientedUniaxialStack<'s,'m>,
    omega_rad_s: f64,
    modes: Vec<Mode>,
    workspace_bytes: usize,
}

impl<'a,'s,'m> PreparedUniaxialStackSpectrum<'a,'s,'m> {
    /// O(modes*layers) preparation, O(modes) point queries, and no query
    /// allocations or repeated source/eigensystem/boundary evaluation.
    pub fn prepare(
        source: &'a OrientedUniaxialStack<'s,'m>, omega_rad_s: f64,
        input: &[UniaxialIncidentMode], limits: UniaxialSpectrumLimits,
        mut cancelled: impl FnMut()->bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let workspace_bytes=admit(source.layers.len(),input.len(),limits)?;
        let data=source.resolve_frequency(omega_rad_s,&mut cancelled)?;
        let modes=prepare_modes(&data,input,&mut cancelled)?;
        Ok(Self {source,omega_rad_s:data.omega,modes,workspace_bytes})
    }
    pub const fn source(&self) -> &'a OrientedUniaxialStack<'s,'m> { self.source }
    pub const fn angular_frequency_rad_s(&self) -> f64 { self.omega_rad_s }
    pub fn retained_mode_count(&self) -> usize { self.modes.len() }
    pub const fn workspace_bytes(&self) -> usize { self.workspace_bytes }

    /// Physical global coordinates [m]. Layer selection retains normal-field
    /// discontinuities at each material face; tangential fields are continuous.
    pub fn field_at(
        &self, region: StackRegion, point: [f64;3], mut cancelled: impl FnMut()->bool,
    ) -> Result<UniaxialSpectrumField> {
        poll(&mut cancelled)?;
        // A completely dark authored spectrum still validates its region.
        if self.modes.is_empty() {
            let supported=match region {
                StackRegion::Layer(index) => {
                    let mut entrance=0.0;
                    let mut interval=None;
                    for (j,layer) in self.source.layers.iter().enumerate() {
                        let exit=entrance+layer.thickness_m;
                        if j==index {interval=Some((entrance,exit)); break;}
                        entrance=exit;
                    }
                    interval.is_some_and(|(a,b)| (a..=b).contains(&point[2]))
                }
                StackRegion::Transmitted => point[2]>=self.source.thickness_m,
                _=>point[2]<=0.0,
            };
            if point.iter().any(|v| !v.is_finite()) || !supported {
                return Err(bad("uniaxial stack point is outside its selected region"));
            }
        }
        field_at(&self.modes,region,point,&mut cancelled)
    }
}

fn admit(layers: usize,modes: usize,limits: UniaxialSpectrumLimits) -> Result<usize> {
    if modes>limits.maximum_modes {
        return Err(bad("uniaxial stack spectrum exceeds its mode/work budget"));
    }
    let bytes=kernel::workspace_bytes(layers,modes)?;
    if bytes>limits.maximum_workspace_bytes {
        return Err(bad("uniaxial stack spectrum exceeds its workspace budget"));
    }
    Ok(bytes)
}
fn prepare_modes(
    data: &Frequency,input: &[UniaxialIncidentMode],cancelled: &mut impl FnMut()->bool,
) -> Result<Vec<Mode>> {
    let mut modes=Vec::new();
    modes.try_reserve_exact(input.len()).map_err(|_| bad("uniaxial stack mode allocation failed"))?;
    for authored in input {
        poll(cancelled)?;
        if !authored.tangential_phase_index.is_finite()
            || authored.electric_v_m.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err(bad("uniaxial stack incident modes must be finite"));
        }
        let axis=super::super::super::kernel::unit_axis(authored.incidence_plane_axis).map_err(bad)?;
        if authored.electric_v_m==[C64::ZERO;2] {continue;}
        let [ex,ey]=authored.electric_v_m;
        let electric=[ex.scale(axis[0])+ey.scale(axis[1]),ey.scale(axis[0])-ex.scale(axis[1])];
        let mut response=kernel::solve(data,Incidence::TangentialIndex(authored.tangential_phase_index),
            electric,axis,cancelled)?;
        response.rotate(axis);
        modes.push(Mode {axis,response});
    }
    poll(cancelled)?;
    Ok(modes)
}
fn field_at(
    modes: &[Mode],region: StackRegion,point: [f64;3],cancelled: &mut impl FnMut()->bool,
) -> Result<UniaxialSpectrumField> {
    poll(cancelled)?;
    let mut sum=FieldSum::new();
    for mode in modes {
        poll(cancelled)?;
        mode.response.accumulate(&mut sum,region,point,mode.axis)?;
    }
    sum.finish()
}

#[cfg(test)]
mod tests;
