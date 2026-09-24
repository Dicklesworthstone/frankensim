//! Distributed wide-tube wall losses on the existing characteristic network.
//!
//! The first-order Zwikker--Kosten telegraph correction is
//! Z'_loss=(1-i)*sqrt(2*rho*mu*omega)/(A*a),
//! Y'_loss=(1-i)*A*(gamma-1)*sqrt(2*mu*omega/rho)/(rho*c*c*a*sqrt(Pr)).
//! This is the same wide-boundary-layer model as fs-duct, NOT an all-shear
//! Bessel or nonlinear gas model. Both viscous and thermal shear must be >=10.
//! fs-phs owns positive Foster identification; fs-vfit owns every retained
//! impedance/admittance state and its midpoint pressure/flow work. No new time
//! integrator, attenuation multiplier, fitted audio gain or delayed force.
//!
//! A physical cell is P/4--Z/2--P/4--Y--P/4--Z/2--P/4. Integer propagation
//! intervals partition the ORIGINAL section transit exactly; no extra delay is
//! introduced by the three instantaneous load nodes. Cold response checks use
//! the actual integer intervals and bilinear loads, not their analogue ideal.

use super::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{RelaxationImpedanceSpec, RelaxationTerm};
use fs_vfit::waveguide::network::admittance::{AdmittanceTerm, RelaxationAdmittanceSpec};

/// Explicit physical section and numerical resolution; unselected sections stay
/// exactly as supplied. Indices refer to the ORIGINAL section list.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscothermalSectionSpec {
    /// Original section index; duplicates refuse.
    pub section: usize,
    /// Number of spatial cells, in 1..=128. Each needs four positive delays.
    pub cells: usize,
    /// Positive low edge of the claimed wide-tube band [Hz].
    pub minimum_frequency_hz: f64,
    /// High edge [Hz], strictly above the low edge; not an output low-pass.
    pub maximum_frequency_hz: f64,
}

/// Physical geometry/gas plus the existing owner's checked eight-arm spectrum.
/// One unit-length spectrum is scaled by each represented cell length.
#[derive(Clone, Debug)]
pub struct WideTubeLoss {
    radius: f64,
    density: f64,
    speed: f64,
    viscosity: f64,
    gamma: f64,
    prandtl: f64,
    dt: f64,
    band: [f64; 2],
    series_gain: f64,
    shunt_gain: f64,
    terms: Vec<(f64, f64)>,
    complex_error: f64,
    resistance_error: f64,
}

/// Mapping and independently checked discretization of ONE supplied section.
#[derive(Clone, Debug)]
pub struct ViscothermalSection {
    /// Original source address and declared physical/numerical band.
    pub source: ViscothermalSectionSpec,
    /// Original endpoints remain unchanged; inserted nodes are appended.
    pub original_nodes: [usize; 2],
    /// Requested physical length [m], not an added radiation end correction.
    pub requested_length_m: f64,
    /// Original integer-transit represented length [m].
    pub represented_length_m: f64,
    /// Half-open range of replacement sections in the lowered graph.
    pub section_range: [usize; 2],
    /// Half-open range of new viscous/thermal load nodes.
    pub node_range: [usize; 2],
    /// Same original total one-way transit, split among the replacement lines.
    pub one_way_samples: usize,
    /// Maximum sampled absolute S11/S22/S21 discrepancy against the homogeneous
    /// wide-tube telegraph section at its REPRESENTED length; fixed gate 0.03.
    pub max_scattering_error: f64,
    /// Physical spectrum, including its separate numerical response checks.
    pub loss: WideTubeLoss,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn checkpoint(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}
fn wave_error(e: fs_vfit::waveguide::WaveguideError) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(e.to_string())
}

impl WideTubeLoss {
    /// Fit only the universal sqrt(omega) kernel, never a pressure waveform.
    /// Poles use a fixed 64-fold guard band and eight existing Foster arms.
    /// Check 65 in-band frequencies: analogue and actual bilinear responses must
    /// each meet 5% relative complex AND real-loss error. Positive fallback
    /// coefficients from the shared fitter confer no admission on their own.
    ///
    /// These are sampled numerical checks of a first-order model, not a physical
    /// validation or all-frequency error bound. No coefficients are retuned.
    ///
    /// # Errors
    /// Invalid gas/geometry/band, narrow boundary layers outside the wide-tube
    /// approximation, unresolved clocks/fits, allocation or cancellation.
    pub fn new(radius: f64, gas: &GasState, dt: f64, band: [f64; 2], gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        checkpoint(gate)?;
        if ![radius, gas.density, gas.sound_speed, gas.dynamic_viscosity, gas.prandtl, dt, band[0], band[1]]
            .iter().all(|x| x.is_finite() && *x > 0.0)
            || !gas.gamma.is_finite() || gas.gamma <= 1.0 || band[1] <= band[0]
            || dt * band[1] > 0.05
        { return Err(invalid("viscothermal loss needs positive geometry/gas/band, gamma>1 and dt*fmax<=0.05")); }
        let rv = radius * det::sqrt(core::f64::consts::TAU * band[0] * gas.density / gas.dynamic_viscosity);
        if !rv.is_finite() || rv.min(rv * det::sqrt(gas.prandtl)) < 10.0 {
            return Err(invalid("wide-tube viscothermal loss requires viscous AND thermal shear >=10 throughout its band"));
        }
        let area = core::f64::consts::PI * radius * radius;
        let series_gain = det::sqrt(2.0 * gas.density * gas.dynamic_viscosity) / (area * radius);
        let shunt_gain = area / (gas.density * gas.sound_speed * gas.sound_speed)
            * (gas.gamma - 1.0) * det::sqrt(2.0 * gas.dynamic_viscosity / gas.density)
            / (radius * det::sqrt(gas.prandtl));
        let low = core::f64::consts::TAU * band[0] / 64.0;
        let high = core::f64::consts::TAU * band[1] * 64.0;
        if ![area, series_gain, shunt_gain, low, high, high*high].iter().all(|x| x.is_finite() && *x > 0.0) {
            return Err(invalid("viscothermal coefficients or identification range overflowed"));
        }
        let terms = fs_phs::foster_sqrt_omega_terms(1.0, low, high, 8)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let mut result = Self { radius, density:gas.density, speed:gas.sound_speed,
            viscosity:gas.dynamic_viscosity, gamma:gas.gamma, prandtl:gas.prandtl,
            dt, band, series_gain, shunt_gain, terms, complex_error:0.0, resistance_error:0.0 };
        for i in 0..=64 {
            checkpoint(gate)?;
            let f = band[0] * det::pow(band[1]/band[0], f64::from(i)/64.0);
            let omega = core::f64::consts::TAU * f;
            let target = C64::new(det::sqrt(omega), -det::sqrt(omega));
            for discrete in [false, true] {
                let value = result.kernel(omega, discrete);
                let complex = (value-target).abs()/target.abs();
                let real = (value.re-target.re).abs()/target.re;
                if !complex.is_finite() || !real.is_finite() || complex>0.05 || real>0.05 {
                    return Err(invalid("viscothermal Foster spectrum exceeds a fixed 5% sampled complex/real-loss error gate"));
                }
                result.complex_error=result.complex_error.max(complex);
                result.resistance_error=result.resistance_error.max(real);
            }
        }
        Ok(result)
    }
    fn kernel(&self, mut omega: f64, discrete: bool) -> C64 {
        if discrete {
            let x=omega*self.dt*0.5;
            omega=2.0*(det::sin(x)/det::cos(x))/self.dt;
        }
        let s=C64::new(0.0,-omega);
        self.terms.iter().fold(C64::ZERO, |v,&(g,p)| v+(s/(s+C64::from_re(p))).scale(g))
    }
    /// Unit-length excess series impedance [Pa s/m^4] and thermal shunt
    /// admittance [m^2/(Pa s)], under exp(-i omega t). No inviscid storage included.
    ///
    /// # Errors
    /// Query outside the admitted positive band.
    pub fn excess_at(&self, frequency_hz: f64, discrete: bool) -> Result<(C64,C64),AcousticRealizeError> {
        if !frequency_hz.is_finite() || frequency_hz<self.band[0] || frequency_hz>self.band[1] {
            return Err(invalid("viscothermal response query is outside its admitted band"));
        }
        let kernel=self.kernel(core::f64::consts::TAU*frequency_hz,discrete);
        Ok((kernel.scale(self.series_gain),kernel.scale(self.shunt_gain)))
    }
    /// Series loss for a positive represented length. No lossless inertance is
    /// added: that storage already belongs to the characteristic lines.
    pub fn series(&self, length: f64) -> Result<RelaxationImpedanceSpec,AcousticRealizeError> {
        if !length.is_finite() || length<=0.0 {return Err(invalid("viscothermal load length must be positive"));}
        let terms:Vec<_>=self.terms.iter().map(|&(g,p)|RelaxationTerm {
            resistance_pa_s_m3:g*self.series_gain*length,rate_per_s:p,
        }).collect();
        RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
            resistance_pa_s_m3:0.0,inertance_pa_s2_m3:0.0,compliance_m3_pa:None,
        },&terms).map_err(wave_error)
    }
    /// Thermal loss for a positive represented length. No adiabatic compliance
    /// is added: characteristic propagation already accounts for it.
    pub fn shunt(&self, length: f64) -> Result<RelaxationAdmittanceSpec,AcousticRealizeError> {
        if !length.is_finite() || length<=0.0 {return Err(invalid("viscothermal load length must be positive"));}
        let terms:Vec<_>=self.terms.iter().map(|&(g,p)|AdmittanceTerm {
            conductance_m3_pa_s:g*self.shunt_gain*length,rate_per_s:p,
        }).collect();
        RelaxationAdmittanceSpec::new(0.0,0.0,&terms).map_err(wave_error)
    }
    /// Physical radius [m].
    #[must_use] pub const fn radius_m(&self)->f64 {self.radius}
    /// Original viscosity [Pa s], not a dimensionless damping knob.
    #[must_use] pub const fn dynamic_viscosity_pa_s(&self)->f64 {self.viscosity}
    /// Original heat-capacity ratio.
    #[must_use] pub const fn gamma(&self)->f64 {self.gamma}
    /// Original gas-derived Prandtl number.
    #[must_use] pub const fn prandtl(&self)->f64 {self.prandtl}
    /// Numerical response discrepancy, excluding wide-model physical error.
    #[must_use] pub const fn max_complex_error(&self)->f64 {self.complex_error}
    /// Separate numerical real-loss discrepancy.
    #[must_use] pub const fn max_real_loss_error(&self)->f64 {self.resistance_error}
}

/// Expand selected sections into positive-delay cells with local viscous and
/// thermal storage. Original node addresses, loads and section endpoints remain
/// intact. No physical state exists yet, so cancellation publishes no partial
/// graph. Empty selection returns the original graph unchanged.
///
/// # Errors
/// Duplicate/absent selections, more than 1024 total cells, unresolvable original
/// transit/quarter intervals, >0.5 rad per cell, failed response or existing
/// network topology/memory/load admission. No numerical limits are relaxed.
pub fn with_viscothermal_sections(mut graph: TubeNetworkSpec, gas:&GasState, dt:f64,
    selections:&[ViscothermalSectionSpec],gate:&CancelGate)
    -> Result<(TubeNetworkSpec,Vec<ViscothermalSection>),AcousticRealizeError>
{
    checkpoint(gate)?;
    if graph.sound_speed_m_s.to_bits()!=gas.sound_speed.to_bits() {
        return Err(invalid("viscothermal medium must equal the network propagation medium"));
    }
    let mut seen=std::collections::BTreeSet::new();
    let mut total=0usize;
    for s in selections {
        if s.section>=graph.sections.len() || !seen.insert(s.section) || !(1..=128).contains(&s.cells) {
            return Err(invalid("viscothermal selections need unique existing sections and 1..=128 cells"));
        }
        total=total.checked_add(s.cells).ok_or_else(||invalid("viscothermal cell count overflow"))?;
    }
    if total>1024 {return Err(invalid("viscothermal expansion exceeds 1024 total cells"));}
    if selections.is_empty() {return Ok((graph,Vec::new()));}
    // Bound the mandatory records before allocating; the final owner additionally
    // admits all wave/energy arrays, load states and lowering scratch.
    let minimum=total.checked_mul(4*core::mem::size_of::<NetworkNode>()+4*core::mem::size_of::<TubeSection>())
        .ok_or_else(||invalid("viscothermal record size overflow"))?;
    if minimum>graph.max_wave_memory_bytes {return Err(invalid("viscothermal records exceed network memory budget"));}
    let original=core::mem::take(&mut graph.sections);
    let mut reports=Vec::new();
    for (index,section) in original.into_iter().enumerate() {
        checkpoint(gate)?;
        let Some(selection)=selections.iter().find(|s|s.section==index) else {graph.sections.push(section);continue;};
        let realized=section.realize(gas.sound_speed,gas.density,dt)?;
        let n=realized.one_way_samples;
        if n<4*selection.cells {return Err(invalid("each viscothermal cell needs four positive propagation intervals"));}
        let loss=WideTubeLoss::new(section.radius_m,gas,dt,
            [selection.minimum_frequency_hz,selection.maximum_frequency_hz],gate)?;
        let section_start=graph.sections.len();let node_start=graph.nodes.len();
        let intervals=4*selection.cells;
        // Rounded cumulative boundaries preserve total transit exactly. Each
        // generated delay is at least one; the checked response uses THESE lengths.
        let boundary=|k:usize| ((n as u128*k as u128+intervals as u128/2)/intervals as u128) as usize;
        let mut left=section.nodes[0];
        for cell in 0..selection.cells {
            checkpoint(gate)?;
            let steps=boundary(4*cell+4)-boundary(4*cell);
            if core::f64::consts::TAU*selection.maximum_frequency_hz*dt*steps as f64>0.5 {
                return Err(invalid("viscothermal spatial cell exceeds 0.5 rad in the declared band; refine the supplied cells/clock"));
            }
            let length=steps as f64*(gas.sound_speed*dt);
            let start=graph.nodes.len();
            let series=loss.series(0.5*length)?;
            graph.nodes.extend([NetworkNode::Series {load:series},
                NetworkNode::ShuntAdmittance {load:loss.shunt(length)?},NetworkNode::Series {load:series}]);
            for quarter in 0..4 {
                let k=4*cell+quarter;
                let last=cell+1==selection.cells && quarter==3;
                let right=if quarter<3 {start+quarter} else if last {section.nodes[1]} else {
                    // Adjacent propagation intervals meet at a zero-volume
                    // junction. It adds no storage, forcing or sample delay.
                    let junction=graph.nodes.len();graph.nodes.push(NetworkNode::Junction);junction
                };
                graph.sections.push(TubeSection {nodes:[left,right],
                    length_m:(boundary(k+1)-boundary(k)) as f64*(gas.sound_speed*dt),
                    radius_m:section.radius_m,max_length_error_m:0.0});
                left=right;
            }
        }
        let mut report=ViscothermalSection {source:*selection,original_nodes:section.nodes,
            requested_length_m:section.length_m,represented_length_m:realized.represented_length_m,
            section_range:[section_start,graph.sections.len()],node_range:[node_start,graph.nodes.len()],
            one_way_samples:n,max_scattering_error:0.0,loss};
        for i in 0..=64 {
            checkpoint(gate)?;
            let f=selection.minimum_frequency_hz*det::pow(selection.maximum_frequency_hz/selection.minimum_frequency_hz,f64::from(i)/64.0);
            let error=section_error(&report,f,&boundary);
            if !error.is_finite() || error>0.03 {return Err(invalid("viscothermal section exceeds fixed 0.03 sampled scattering discrepancy"));}
            report.max_scattering_error=report.max_scattering_error.max(error);
        }
        reports.push(report);
    }
    // Ask the actual owner to admit topology, port scattering and load histories.
    // This temporary zero-state object is cold admission, not a second solver.
    let mut segments=Vec::new();
    for s in &graph.sections {let r=s.realize(gas.sound_speed,gas.density,dt)?;
        segments.push(fs_vfit::waveguide::network::NetworkSegment {nodes:s.nodes,
            one_way_samples:r.one_way_samples,impedance_pa_s_m3:r.impedance_pa_s_m3});}
    fs_vfit::waveguide::network::WaveguideNetwork::new(&graph.nodes,&segments,dt,graph.max_wave_memory_bytes)
        .map_err(wave_error)?;
    checkpoint(gate)?;
    Ok((graph,reports))
}

type Matrix=[C64;4];
fn mul(a:Matrix,b:Matrix)->Matrix {
    [a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3]]
}
fn scattering(m:Matrix)->[C64;3] {
    let d=m[0]+m[1]+m[2]+m[3];
    [(m[0]+m[1]-m[2]-m[3])/d,(-m[0]+m[1]-m[2]+m[3])/d,C64::from_re(2.0)/d]
}
fn section_error(r:&ViscothermalSection,f:f64,boundary:&impl Fn(usize)->usize)->f64 {
    let w=core::f64::consts::TAU*f;let loss=&r.loss;
    let area=core::f64::consts::PI*loss.radius*loss.radius;
    let zc=loss.density*loss.speed/area;
    let ideal=C64::new(det::sqrt(w),-det::sqrt(w));
    let z=C64::new(0.0,-w*loss.density/area)+ideal.scale(loss.series_gain);
    let y=C64::new(0.0,-w*area/(loss.density*loss.speed*loss.speed))+ideal.scale(loss.shunt_gain);
    let propagation=(z*y).sqrt();let g=propagation.scale(r.represented_length_m);
    let ep=det::exp(g.re);let em=det::exp(-g.re);
    let ch=C64::new(0.5*(ep+em)*det::cos(g.im),0.5*(ep-em)*det::sin(g.im));
    let sh=C64::new(0.5*(ep-em)*det::cos(g.im),0.5*(ep+em)*det::sin(g.im));
    let target=scattering([ch,z*sh/propagation.scale(zc),y*sh/propagation.scale(1.0/zc),ch]);
    let mut m=[C64::ONE,C64::ZERO,C64::ZERO,C64::ONE];
    let kernel=loss.kernel(w,true);
    for cell in 0..r.source.cells {
        let length=(boundary(4*cell+4)-boundary(4*cell)) as f64*(loss.speed*loss.dt);
        for q in 0..4 {
            let k=4*cell+q;let angle=w*loss.dt*(boundary(k+1)-boundary(k)) as f64;
            let c=C64::from_re(det::cos(angle));let s=C64::new(0.0,-det::sin(angle));
            m=mul(m,[c,s,s,c]);
            if q==0 || q==2 {m=mul(m,[C64::ONE,kernel.scale(loss.series_gain*length/(2.0*zc)),C64::ZERO,C64::ONE]);}
            if q==1 {m=mul(m,[C64::ONE,C64::ZERO,kernel.scale(loss.shunt_gain*length*zc),C64::ONE]);}
        }
    }
    let actual=scattering(m);
    if actual.iter().chain(&target).any(|x| !x.re.is_finite() || !x.im.is_finite()) { return f64::INFINITY; }
    actual.iter().zip(target).fold(0.0_f64,|worst,(a,b)|worst.max((*a-b).abs()))
}
