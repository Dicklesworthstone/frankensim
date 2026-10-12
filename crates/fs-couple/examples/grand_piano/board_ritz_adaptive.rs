//! Adaptive Ritz enrichment against the supplied source bridge-transfer matrix.
//! Each sweep checks every reciprocal port pair at statics, every distinct source
//! pole, adjacent source-frequency midpoints, and explicit seeds. The finite grid
//! is not a receiver, mesh, or continuous-band accuracy certificate.

use fs_la::eigen_complex::lu_complex;
use fs_math::c64::C64;
use super::{AdaptiveResponseReport, RitzBasis, RitzOptions, dot, normalize, projected_forms, remove, TAU};

// Cold preparation only. Count source-coordinate visits and conservatively
// weighted complex dense work, not wall time or a hard floating-operation count.
// This separate budget admits useful all-key sweeps; the fixed snapshot path's
// original 250-million budget is unchanged.
const MAX_RESPONSE_WORK: u64 = 8_000_000_000;
const ENRICHMENT_POINTS: usize = 16;

struct Work(u64);
impl Work {
    fn charge(&mut self, amount: u64) -> Result<(), String> {
        self.0 = self.0.checked_add(amount).ok_or("adaptive Ritz work estimate overflow")?;
        if self.0 > MAX_RESPONSE_WORK {
            return Err(format!("adaptive bridge response work estimate {} exceeds {MAX_RESPONSE_WORK}; use fewer source modes, retained modes or admitted bridge ports", self.0));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug)]
struct Frequency { hz: f64, omega: f64, squared: f64, scale: f64 }

fn frequency(hz: f64, lambda_max: f64, damping_max: f64) -> Result<Frequency, String> {
    let omega = TAU * hz;
    let squared = omega * omega;
    if !omega.is_finite() || !squared.is_finite() || (hz > 0.0 && squared == 0.0) {
        return Err("adaptive Ritz angular frequency is not representable".into());
    }
    let scale = lambda_max.max(squared).max(omega * damping_max);
    if !scale.is_finite() || scale <= 0.0 {
        return Err("adaptive Ritz dynamic stiffness is not representable".into());
    }
    Ok(Frequency { hz, omega, squared, scale })
}

fn grid(lambda: &[f64], damping: &[f64], seeds: &[f64]) -> Result<Vec<Frequency>, String> {
    let max_k = lambda[lambda.len()-1];
    let max_c = damping.iter().copied().fold(0.0, f64::max);
    let mut points = vec![frequency(0.0, max_k, max_c)?];
    for (&k, &c) in lambda.iter().zip(damping) {
        // These points ARE the source poles, even when sqrt(k)/TAU round-trips
        // to a slightly different float. An undamped pole has no finite target;
        // neither a damping floor nor a nearby substitute is admitted.
        if c == 0.0 {
            return Err("adaptive bridge response grid contains an exact undamped source pole".into());
        }
        let omega = k.sqrt();
        let mut point = frequency(omega/TAU, max_k, max_c)?;
        point.omega = omega;
        point.squared = k;
        point.scale = max_k.max(k).max(omega*max_c);
        if !point.scale.is_finite() {
            return Err("adaptive Ritz source-pole stiffness is not representable".into());
        }
        if let Some(last)=points.last() {
            if last.squared==k {continue;}
            if last.hz>=point.hz {
                return Err("distinct adaptive source poles have no distinguishable representable Hz grid".into());
            }
        }
        points.push(point);
    }
    let midpoints = points.windows(2).map(|pair|
        frequency(f64::midpoint(pair[0].hz,pair[1].hz),max_k,max_c))
        .collect::<Result<Vec<_>,_>>()?;
    points.extend(midpoints);
    for &hz in seeds { points.push(frequency(hz,max_k,max_c)?); }
    points.sort_by(|a,b| a.hz.total_cmp(&b.hz));
    points.dedup_by(|a,b| a.hz == b.hz);
    Ok(points)
}

fn finite(z: C64) -> bool { z.re.is_finite() && z.im.is_finite() }

fn relative_error(got: C64, expected: C64) -> Result<f64,String> {
    error_at_scale(got,expected,expected.abs())
}

fn error_at_scale(got: C64, expected: C64, reference: f64) -> Result<f64,String> {
    let absolute=(got-expected).abs();
    if !absolute.is_finite() || !reference.is_finite() {
        return Err("adaptive compliance error magnitude is not representable".into());
    }
    if reference==0.0 {
        return Ok(if absolute==0.0 {0.0} else {f64::INFINITY});
    }
    let relative=absolute/reference;
    if !relative.is_finite() {
        return Err("adaptive relative compliance error is not representable".into());
    }
    Ok(relative)
}

/// A cross transfer may vanish or change sign at an antiresonance. Its own
/// magnitude is therefore not a suitable denominator. Driving-point magnitudes
/// give a force-unit invariant scale without inventing a compliance floor.
fn transfer_error(got: C64, expected: C64, source_i: C64, source_j: C64)
    -> Result<f64,String> {
    let first=source_i.abs();let second=source_j.abs();
    if !first.is_finite() || !second.is_finite() {
        return Err("adaptive compliance error magnitude is not representable".into());
    }
    // Taking roots separately avoids overflowing or underflowing |Hii|*|Hjj|.
    error_at_scale(got,expected,first.sqrt()*second.sqrt())
}

/// Both source and reduced compliances are multiplied by the same positive
/// frequency-local stiffness scale. Port rescaling also cancels in relative
/// error, avoiding arbitrary force-unit or very small bridge-amplitude floors.
struct Responses {
    points: Vec<Frequency>,
    inverse: Vec<C64>,
    compliance: Vec<C64>,
}
impl Responses {
    fn new(lambda: &[f64], damping: &[f64], ports: &[Vec<f64>], seeds: &[f64],
        work: &mut Work) -> Result<Self,String> {
        let points = grid(lambda,damping,seeds)?;
        let n=lambda.len();let p=ports.len();
        work.charge(points.len() as u64*n as u64*(8+4*p as u64))?;
        let mut inverse=Vec::with_capacity(points.len()*n);
        let mut compliance=Vec::with_capacity(points.len()*p);
        for point in &points {
            let start=inverse.len();
            for i in 0..n {
                let a=(lambda[i]-point.squared)/point.scale;
                let b=-(point.omega*damping[i])/point.scale;
                let value=C64::ONE/C64::new(a,b);
                if !finite(value) {
                    return Err("adaptive source response is singular or not representable".into());
                }
                inverse.push(value);
            }
            for port in ports {
                let mut value=C64::ZERO;
                for i in 0..n { value=value+inverse[start+i].scale(port[i]*port[i]); }
                if !finite(value) { return Err("adaptive source compliance overflow".into()); }
                compliance.push(value);
            }
        }
        Ok(Self {points,inverse,compliance})
    }
}

#[derive(Clone, Copy)]
struct Miss { error: f64, frequency: usize, port: usize }

fn sweep(responses: &Responses, ports: &[Vec<f64>], columns: &[Vec<f64>],
    k: &[f64], c: &[f64], tolerance: f64, work: &mut Work)
    -> Result<(f64,f64,Vec<Miss>),String> {
    let n=ports[0].len();let r=columns.len();let p=ports.len();
    let rank=r as u64;
    let pairs=p as u64*(p as u64-1)/2;
    work.charge(n as u64*rank*p as u64 + responses.points.len() as u64*
        (2*rank*rank*rank+4*p as u64*rank*rank+8*rank*rank+
            pairs*(4*n as u64+4*rank+32)+32*p as u64))?;
    let forces: Vec<Vec<f64>>=ports.iter().map(|port|
        columns.iter().map(|q|dot(port,q)).collect()).collect();
    let mut matrix=vec![C64::ZERO;r*r];let mut rhs=vec![C64::ZERO;r];
    let mut misses: Vec<Miss>=Vec::with_capacity(ENRICHMENT_POINTS);
    let mut maximum_driving=0.0_f64;let mut maximum_transfer=0.0_f64;
    for (fi,point) in responses.points.iter().enumerate() {
        for row in 0..r {for col in 0..r {
            matrix[row*r+col]=C64::new(
                (k[row*r+col]-if row==col {point.squared} else {0.0})/point.scale,
                -(point.omega*c[row*r+col])/point.scale);
        }}
        if matrix.iter().any(|&z|!finite(z)) {
            return Err("adaptive reduced dynamic stiffness overflow".into());
        }
        let factor=lu_complex(&matrix,r).map_err(|_|
            "adaptive reduced bridge response is singular; no damping floor was added")?;
        for (pi,force) in forces.iter().enumerate() {
            for (value,&b) in rhs.iter_mut().zip(force) { *value=C64::from_re(b); }
            factor.solve(&mut rhs);
            let mut got=C64::ZERO;
            for (&b,&q) in force.iter().zip(&rhs) { got=got+q.scale(b); }
            if !finite(got) { return Err("adaptive reduced compliance overflow".into()); }
            let expected=responses.compliance[fi*p+pi];
            let relative=relative_error(got,expected)?;
            maximum_driving=maximum_driving.max(relative);
            let mut drive_error=relative;
            // Reciprocity permits checking one triangle. Reuse this driven
            // response for every output port; no extra solve or all-frequency
            // transfer-matrix cache is needed.
            for pj in pi+1..p {
                let mut cross_got=C64::ZERO;let mut cross_expected=C64::ZERO;
                for (&b,&q) in forces[pj].iter().zip(&rhs) {
                    cross_got=cross_got+q.scale(b);
                }
                for i in 0..n {
                    cross_expected=cross_expected+
                        responses.inverse[fi*n+i].scale(ports[pi][i]*ports[pj][i]);
                }
                if !finite(cross_got) || !finite(cross_expected) {
                    return Err("adaptive bridge transfer overflow".into());
                }
                drive_error=drive_error.max(transfer_error(cross_got,cross_expected,
                    expected,responses.compliance[fi*p+pj])?);
            }
            maximum_transfer=maximum_transfer.max(drive_error);
            // Aggregate output-port misses before enrichment: a driven source
            // response is useful once per frequency, even if many outputs fail.
            if drive_error>tolerance {
                let position=misses.iter().position(|old|drive_error>old.error).unwrap_or(misses.len());
                if position<ENRICHMENT_POINTS {
                    misses.insert(position,Miss {error:drive_error,frequency:fi,port:pi});
                    if misses.len()>ENRICHMENT_POINTS { misses.pop(); }
                }
            }
        }
    }
    Ok((maximum_driving,maximum_transfer,misses))
}

fn append_response(responses: &Responses, point: usize, port: &[f64], columns: &mut Vec<Vec<f64>>,
    snapshots: &mut Vec<Vec<f64>>, limit: usize, keep: usize, work: &mut Work)
    -> Result<(),String> {
    let n=port.len();
    work.charge(n as u64*(24+12*limit as u64))?;
    let response=&responses.inverse[point*n..(point+1)*n];
    let mut real: Vec<_>=response.iter().zip(port).map(|(z,b)|z.re*b).collect();
    let mut imaginary: Vec<_>=response.iter().zip(port).map(|(z,b)|z.im*b).collect();
    for v in [&real,&imaginary] {
        let mut v=v.clone();
        if normalize(&mut v)? { snapshots.push(v); }
    }
    // Joint scaling preserves the relative size of real/imaginary residuals.
    // Prefer the stronger direction if only one rank slot remains.
    let scale=real.iter().chain(&imaginary).map(|x|x.abs()).fold(0.0,f64::max);
    if scale==0.0 { return Ok(()); }
    for v in [&mut real,&mut imaginary] {
        for x in v.iter_mut() { *x/=scale; }
        v[..keep].fill(0.0);
        for _ in 0..2 {for q in columns.iter() {remove(v,q);}}
    }
    if dot(&imaginary,&imaginary)>dot(&real,&real) { std::mem::swap(&mut real,&mut imaginary); }
    let threshold=(64.0*f64::EPSILON*n as f64).powi(2);
    for mut v in [real,imaginary] {
        if columns.len()==limit { break; }
        for _ in 0..2 {for q in columns.iter() {remove(&mut v,q);}}
        if dot(&v,&v)>threshold {
            normalize(&mut v)?;
            columns.push(v);
        }
    }
    Ok(())
}

pub(super) fn build(lambda: &[f64], damping: &[f64], ports: &[Vec<f64>], options: &RitzOptions)
    -> Result<RitzBasis,String> {
    let tolerance=options.transfer_tolerance.ok_or("missing adaptive response tolerance")?;
    let n=lambda.len();let target=options.max_modes.min(n);let keep=options.keep_low_modes;
    let mut work=Work(0);
    let normalized_ports: Vec<Vec<f64>>=ports.iter().map(|port| {
        let scale=port.iter().map(|x|x.abs()).fold(0.0,f64::max);
        port.iter().map(|x|if scale==0.0 {0.0} else {*x/scale}).collect()
    }).collect();
    let responses=Responses::new(lambda,damping,&normalized_ports,&options.sample_hz,&mut work)?;
    let mut columns=Vec::with_capacity(target);
    for low in 0..keep {let mut q=vec![0.0;n];q[low]=1.0;columns.push(q);}
    let mut snapshots=Vec::new();
    // A small initial space leaves room for high-band errors to choose the
    // remaining directions, instead of spending the budget on low-band seeds.
    let initial=(keep+4).min(target);
    let seed_indices: Vec<_>=responses.points.iter().enumerate().filter_map(|(i,point)|
        (point.hz==0.0 || options.sample_hz.contains(&point.hz)).then_some(i)).collect();
    'seed: for port in &normalized_ports {for &point in &seed_indices {
        if columns.len()==initial {break 'seed;}
        append_response(&responses,point,port,&mut columns,&mut snapshots,initial,keep,&mut work)?;
    }}
    if columns.is_empty() {
        return Err("adaptive Ritz reduction has no nonzero bridge response or protected low mode".into());
    }
    let mut rounds=0usize;
    loop {
        work.charge(4*n as u64*columns.len() as u64*columns.len() as u64)?;
        let (stiffness,physical_damping)=projected_forms(lambda,damping,&columns)?;
        let (driving_error,error,misses)=sweep(&responses,&normalized_ports,&columns,
            &stiffness,&physical_damping,tolerance,&mut work)?;
        if error<=tolerance {
            work.charge(snapshots.len() as u64*n as u64*(2*columns.len() as u64+2))?;
            let mut snapshot_error=0.0_f64;
            for snapshot in &snapshots {
                let mut residual=snapshot.clone();
                for q in &columns {remove(&mut residual,q);}
                snapshot_error=snapshot_error.max(dot(&residual,&residual).sqrt());
            }
            return Ok(RitzBasis {columns,stiffness,damping:physical_damping,
                max_relative_snapshot_error:snapshot_error,snapshot_count:snapshots.len(),
                adaptive_response:Some(AdaptiveResponseReport {tolerance,
                    max_relative_driving_point_error:driving_error,
                    max_normalized_transfer_error:error,
                    evaluation_hz:responses.points.iter().map(|point|point.hz).collect(),
                    port_count:ports.len(),enrichment_rounds:rounds})});
        }
        let worst=misses.first().ok_or("adaptive bridge response error is not representable")?;
        if columns.len()==target {
            return Err(format!("adaptive bridge transfer tolerance {tolerance} unmet: normalized error {error} at {} Hz, driven port {}, with {target} retained modes",
                responses.points[worst.frequency].hz,worst.port));
        }
        let previous=columns.len();
        for miss in misses {
            append_response(&responses,miss.frequency,&normalized_ports[miss.port],
                &mut columns,&mut snapshots,target,keep,&mut work)?;
            if columns.len()==target {break;}
        }
        if columns.len()==previous {
            return Err(format!("adaptive bridge transfer tolerance {tolerance} unmet: normalized error {error}; response enrichment has no independent direction"));
        }
        rounds+=1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::bridge_basis;

    fn source(frequencies: &[f64]) -> (Vec<f64>,Vec<f64>) {
        let k:Vec<_>=frequencies.iter().map(|hz|(TAU*hz).powi(2)).collect();
        let c=k.iter().map(|value|0.02*value.sqrt()).collect();
        (k,c)
    }
    fn compliance(k: &[f64], c: &[f64], port: &[f64], hz: f64) -> C64 {
        transfer(k,c,port,port,hz)
    }
    fn transfer(k: &[f64], c: &[f64], drive: &[f64], receiver: &[f64], hz: f64) -> C64 {
        let omega=TAU*hz;
        drive.iter().zip(receiver).enumerate().fold(C64::ZERO,|sum,(i,(a,b))|
            sum+C64::from_re(a*b)/C64::new(k[i]-omega*omega,-omega*c[i]))
    }
    fn reduced_compliance(basis: &RitzBasis, port: &[f64], hz: f64) -> C64 {
        reduced_transfer(basis,port,port,hz)
    }
    fn reduced_transfer(basis: &RitzBasis, drive: &[f64], receiver: &[f64], hz: f64) -> C64 {
        let r=basis.columns.len();let omega=TAU*hz;
        let force:Vec<_>=basis.columns.iter().map(|q|dot(drive,q)).collect();
        let output:Vec<_>=basis.columns.iter().map(|q|dot(receiver,q)).collect();
        let matrix:Vec<_>=(0..r*r).map(|i|C64::new(
            basis.stiffness[i]-if i/r==i%r {omega*omega} else {0.0},-omega*basis.damping[i])).collect();
        let mut rhs:Vec<_>=force.iter().map(|&b|C64::from_re(b)).collect();
        lu_complex(&matrix,r).unwrap().solve(&mut rhs);
        rhs.iter().zip(output).fold(C64::ZERO,|sum,(&q,b)|sum+q.scale(b))
    }

    #[test]
    fn tiny_seed_snapshot_error_cannot_hide_missed_source_resonances() {
        let (k,c)=source(&[100.0,1000.0,1100.0]);let ports=vec![vec![1.0;3]];
        let fixed=bridge_basis(&k,&c,&ports,&RitzOptions::parse("2,1,100").unwrap()).unwrap();
        assert!(fixed.max_relative_snapshot_error<2e-5);
        for hz in [1000.0,1100.0] {
            let exact=compliance(&k,&c,&ports[0],hz);
            assert!((reduced_compliance(&fixed,&ports[0],hz)-exact).abs()/exact.abs()>0.95);
        }
        let error=bridge_basis(&k,&c,&ports,&RitzOptions::parse("adaptive,2,1,0.01,100").unwrap()).unwrap_err();
        assert!(error.contains("tolerance") && error.contains("unmet"),"{error}");
    }

    #[test]
    fn high_band_enrichment_checks_every_transfer_and_preserves_physics() {
        let (k,mut c)=source(&[100.0,350.0,900.0,1000.0,1100.0,2200.0,2400.0]);
        for (i,value) in c.iter_mut().enumerate() { *value*=1.0+0.17*i as f64; }
        // The second port owns two high modes absent from the first. The final
        // source mode is unforced and needs no invented filler direction.
        let ports=vec![vec![1.,0.3,0.7,1.,0.,0.,0.],vec![0.2,0.,0.,0.,1.,0.9,0.]];
        let options=RitzOptions::parse("adaptive,6,1,0.000001,100").unwrap();
        let basis=bridge_basis(&k,&c,&ports,&options).unwrap();
        let report=basis.adaptive_response.as_ref().unwrap();
        assert_eq!(report.port_count,2);assert!(report.enrichment_rounds>0);
        assert!(basis.columns.len()<=6 && basis.columns.len()<k.len());
        assert_eq!(basis.columns[0],vec![1.,0.,0.,0.,0.,0.,0.]);
        assert!(basis.columns[1..].iter().all(|q|q[0]==0.0));
        assert!(report.max_relative_driving_point_error<=report.tolerance);
        assert!(report.max_normalized_transfer_error<=report.tolerance);
        for &hz in &report.evaluation_hz {for drive in &ports {for receiver in &ports {
            let exact=transfer(&k,&c,drive,receiver,hz);
            let scale=compliance(&k,&c,drive,hz).abs().sqrt()*
                compliance(&k,&c,receiver,hz).abs().sqrt();
            let relative=(reduced_transfer(&basis,drive,receiver,hz)-exact).abs()/scale;
            assert!(relative<=report.tolerance,"{hz}Hz: {relative}");
        }}}
        for (i,a) in basis.columns.iter().enumerate() {for (j,b) in basis.columns.iter().enumerate() {
            assert!((dot(a,b)-if i==j {1.0} else {0.0}).abs()<2e-13);
            let r=basis.columns.len();
            assert_eq!(basis.stiffness[i*r+j],basis.stiffness[j*r+i]);
            assert_eq!(basis.damping[i*r+j],basis.damping[j*r+i]);
            let expected:f64=a.iter().zip(b).zip(&c).map(|((x,y),c)|x*y*c).sum();
            assert!((basis.damping[i*r+j]-expected).abs()<1e-11);
        }}
        assert_eq!(basis.columns,bridge_basis(&k,&c,&ports,&options).unwrap().columns);
    }

    #[test]
    fn diagonal_pass_cannot_hide_signed_cross_bridge_antiresonance() {
        let (k,_)=source(&[100.,100.,200.,200.]);
        let c=vec![0.02*TAU*100.;4];
        let strong=0.95_f64.sqrt();let weak=0.05_f64.sqrt();
        // The first drive's static/harmonic seeds span only modes 0 and 2.
        // Each other drive loses exactly 5% of its SELF compliance at every
        // frequency. Between the two resonances, however, their omitted cross
        // contributions add while self contributions nearly cancel.
        let ports=vec![vec![1.,0.,1.,0.],vec![strong,weak,strong,weak],
            vec![strong,weak,strong,-weak],vec![0.;4]];
        let mut options=RitzOptions {max_modes:2,keep_low_modes:0,
            sample_hz:vec![100.*2.5_f64.sqrt()],transfer_tolerance:Some(0.06)};
        let mut work=Work(0);
        let responses=Responses::new(&k,&c,&ports,&options.sample_hz,&mut work).unwrap();
        let columns=vec![vec![1.,0.,0.,0.],vec![0.,0.,1.,0.]];
        let (projected_k,projected_c)=projected_forms(&k,&c,&columns).unwrap();
        let (diagonal,cross,misses)=sweep(&responses,&ports,&columns,
            &projected_k,&projected_c,0.06,&mut work).unwrap();
        assert!((diagonal-0.05).abs()<2e-12,"{diagonal}");
        assert!(cross>2.0,"cross error {cross} hidden by diagonal error {diagonal}");
        for (i,miss) in misses.iter().enumerate() {
            assert!(misses[..i].iter().all(|other|
                other.frequency!=miss.frequency || other.port!=miss.port));
        }
        let refusal=bridge_basis(&k,&c,&ports,&options).unwrap_err();
        assert!(refusal.contains("bridge transfer tolerance") && refusal.contains("unmet"),"{refusal}");
        options.max_modes=4;
        let recovered=bridge_basis(&k,&c,&ports,&options).unwrap();
        let report=recovered.adaptive_response.as_ref().unwrap();
        assert!(report.max_relative_driving_point_error<1e-11);
        assert!(report.max_normalized_transfer_error<1e-11);
        let probe=options.sample_hz[0];
        let exact=transfer(&k,&c,&ports[1],&ports[2],probe);
        let got=reduced_transfer(&recovered,&ports[1],&ports[2],probe);
        assert!((got-exact).abs()<1e-11*exact.abs());
        // Signed force-unit changes must not bypass the same admission.
        let scaled:Vec<Vec<_>>=ports.iter().zip([-1e100,1e-100,-3e75,1.])
            .map(|(port,scale)|port.iter().map(|b|b*scale).collect()).collect();
        options.max_modes=2;
        assert!(bridge_basis(&k,&c,&scaled,&options).is_err());
        options.max_modes=4;
        assert!(bridge_basis(&k,&c,&scaled,&options).unwrap().adaptive_response
            .unwrap().max_normalized_transfer_error<1e-11);
    }

    #[test]
    fn cross_transfer_scale_preserves_zero_ports_and_extreme_finite_units() {
        assert_eq!(transfer_error(C64::ZERO,C64::ZERO,C64::ZERO,C64::ONE).unwrap(),0.0);
        assert_eq!(transfer_error(C64::ONE,C64::ZERO,C64::ZERO,C64::ONE).unwrap(),f64::INFINITY);
        assert_eq!(transfer_error(C64::ONE,C64::ZERO,C64::ONE,C64::ONE).unwrap(),1.0);
        for magnitude in [1e-200,1e200] {
            let value=C64::from_re(magnitude);
            assert!((transfer_error(value,C64::ZERO,value,value).unwrap()-1.0).abs()<1e-14);
        }
        let huge=C64::new(f64::MAX,f64::MAX);
        assert!(transfer_error(C64::ZERO,C64::ZERO,huge,C64::ONE).is_err());
    }

    #[test]
    fn adaptive_parser_poles_grid_and_work_limits_are_explicit() {
        let options=RitzOptions::parse("adaptive,32,4,0.05,100,1500").unwrap();
        assert_eq!(options.transfer_tolerance,Some(0.05));
        assert!(RitzOptions::parse("32,4,100,1500").unwrap().transfer_tolerance.is_none());
        for text in ["adaptive,32,4,0,100","adaptive,32,4,1,100",
            "adaptive,32,4,NaN,100","adaptive,32,4,0.1","adaptive,32,4,0.1,100,90"] {
            assert!(RitzOptions::parse(text).is_err(),"{text}");
        }
        let (k,c)=source(&[100.0,1000.0,1100.0]);
        let points=grid(&k,&c,&[75.,100.]).unwrap();
        for expected in [0.,50.,75.,100.,550.,1000.,1050.,1100.] {
            assert!(points.iter().any(|point|(point.hz-expected).abs()<1e-10),"{expected}");
        }
        let error=bridge_basis(&k,&[0.,c[1],c[2]],&[vec![1.;3]],
            &RitzOptions::parse("adaptive,3,1,0.01,100").unwrap()).unwrap_err();
        assert!(error.contains("exact undamped source pole"));
        let mut work=Work(MAX_RESPONSE_WORK-1);work.charge(1).unwrap();
        assert!(work.charge(1).unwrap_err().contains("work estimate"));
    }

    #[test]
    fn unrepresentable_error_norms_and_colliding_source_poles_cannot_pass() {
        let huge=C64::new(f64::MAX,f64::MAX);
        assert!(finite(huge));
        assert!(relative_error(C64::ZERO,huge).is_err());
        assert!(relative_error(C64::from_re(f64::MAX),C64::from_re(-f64::MAX)).is_err());
        assert!(relative_error(C64::ONE,C64::from_re(f64::MIN_POSITIVE*f64::EPSILON)).is_err());
        assert_eq!(relative_error(C64::ZERO,C64::ZERO).unwrap(),0.0);
        assert_eq!(relative_error(C64::ONE,C64::ZERO).unwrap(),f64::INFINITY);
        let collided=grid(&[1.0,1.0_f64.next_up()],&[0.01,0.01],&[0.1]).unwrap_err();
        assert!(collided.contains("distinguishable representable Hz"),"{collided}");
    }
}
