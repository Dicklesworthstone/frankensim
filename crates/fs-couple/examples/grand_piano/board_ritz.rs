//! Cold bridge-response Ritz spaces in an existing mass-normalized modal basis.
//! No eigenproblem is solved here. The source model has M=I, positive diagonal
//! K and nonnegative diagonal C. Projection retains the complete K and C;
//! a caller may subsequently diagonalize K through the existing modal owner.
//! The snapshot error is a relative modal displacement projection norm, not
//! a transfer-function error bound, source-spectrum certificate or convergence
//! claim. Source eigenvalue certificates do not transfer to Ritz frequencies.

use std::f64::consts::TAU;

pub const MAX_SOURCE_MODES: usize = 512;
pub const MAX_RITZ_MODES: usize = 128;
pub const MAX_SAMPLE_FREQUENCIES: usize = 16;
const MAX_PORTS: usize = 88;
/// Conservative source-coordinate visits in snapshot orthogonalization,
/// residual measurement and projected-matrix construction, checked before
/// allocating snapshots. This is a work admission bound, not a timing claim.
const MAX_SCALAR_WORK: usize = 250_000_000;

#[derive(Clone, Debug)]
pub struct RitzOptions {
    pub max_modes: usize,
    pub keep_low_modes: usize,
    /// Positive, strictly increasing harmonic samples; statics are implicit.
    pub sample_hz: Vec<f64>,
}

impl RitzOptions {
    /// `max_modes,keep_low_modes,frequency_hz,...`, with 1..=16 samples.
    pub fn parse(text: &str) -> Result<Self, String> {
        let fields: Vec<_> = text.split(',').take(MAX_SAMPLE_FREQUENCIES + 3).collect();
        if !(3..=MAX_SAMPLE_FREQUENCIES + 2).contains(&fields.len()) {
            return Err("board reduction needs max_modes,keep_low_modes and 1..=16 frequencies in Hz".into());
        }
        let options = Self {
            max_modes: fields[0].trim().parse().map_err(|_| "invalid Ritz mode budget")?,
            keep_low_modes: fields[1].trim().parse().map_err(|_| "invalid protected low-mode count")?,
            sample_hz: fields[2..].iter().map(|s| s.trim().parse::<f64>()
                .map_err(|_| "invalid Ritz sample frequency".to_string())).collect::<Result<_,_>>()?,
        };
        options.validate()?;
        Ok(options)
    }

    pub(crate) fn validate(&self) -> Result<(), String> {
        if !(1..=MAX_RITZ_MODES).contains(&self.max_modes) || self.keep_low_modes > self.max_modes
            || !(1..=MAX_SAMPLE_FREQUENCIES).contains(&self.sample_hz.len())
            || self.sample_hz.iter().any(|f| !f.is_finite() || *f <= 0.0)
            || self.sample_hz.windows(2).any(|w| w[0] >= w[1]) {
            return Err("Ritz reduction requires 1..=128 modes, keep_low_modes <= max_modes, and 1..=16 finite positive strictly increasing frequencies".into());
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct RitzBasis {
    /// Columns in source-modal coordinates. Protected low columns are exact
    /// canonical vectors; every other column is zero in those coordinates.
    pub columns: Vec<Vec<f64>>,
    /// Full row-major matrices in the returned orthonormal basis.
    pub stiffness: Vec<f64>,
    pub damping: Vec<f64>,
    /// Maximum ||x-Q Q^T x|| / ||x|| over nonzero real-valued snapshots.
    /// Mass-normalized displacement only; no actual transfer-error certificate.
    pub max_relative_snapshot_error: f64,
    /// Nonzero static, harmonic-real and harmonic-imaginary snapshots, before
    /// rank selection. Duplicate directions still count as supplied snapshots.
    pub snapshot_count: usize,
}

fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x,y)| x*y).sum() }

/// Normalize without forming an overflowing physical snapshot norm.
fn normalize(v: &mut [f64]) -> Result<bool, String> {
    if v.iter().any(|x| !x.is_finite()) { return Err("nonfinite Ritz response snapshot".into()); }
    let scale = v.iter().map(|x| x.abs()).fold(0.0, f64::max);
    if scale == 0.0 { return Ok(false); }
    let norm = v.iter().map(|x| (x/scale).powi(2)).sum::<f64>().sqrt();
    for x in v { *x = (*x/scale)/norm; }
    Ok(true)
}

fn remove(v: &mut [f64], q: &[f64]) {
    let coefficient = dot(v,q);
    for (x,b) in v.iter_mut().zip(q) { *x -= coefficient*b; }
}

/// Build a bounded space from selected bridge forces in the complete supplied
/// source-modal model. `lambda` is ascending angular-frequency squared;
/// `damping[i] = 2*zeta_i*sqrt(lambda[i])` is the physical viscous rate.
/// Each port contains one force coefficient per source mode. The caller must
/// retain the original slice count/authority separately and use the SAME final
/// nodal transformation for bridge mechanics, surface motion and acoustics.
pub fn bridge_basis(lambda: &[f64], damping: &[f64], ports: &[Vec<f64>], options: &RitzOptions)
    -> Result<RitzBasis, String> {
    options.validate()?;
    let n = lambda.len();
    if !(1..=MAX_SOURCE_MODES).contains(&n) || damping.len() != n || options.keep_low_modes > n
        || lambda.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || lambda.windows(2).any(|w| w[0] > w[1])
        || damping.iter().any(|x| !x.is_finite() || *x < 0.0)
        || !(1..=MAX_PORTS).contains(&ports.len())
        || ports.iter().any(|p| p.len() != n || p.iter().any(|x| !x.is_finite())) {
        return Err("Ritz reduction needs 1..=512 ascending positive source modes, matching nonnegative damping, 1..=88 finite bridge vectors and an admitted low-mode count".into());
    }
    let target = options.max_modes.min(n);
    let extra = target-options.keep_low_modes;
    let snapshot_bound = ports.len()*(1+2*options.sample_hz.len());
    // Five visits per candidate/added column: twice dot+subtract, then norm.
    // Final direct projection uses two per column. Remaining terms bound
    // snapshot normalization, pivot reorthogonalization and both dense forms.
    let work = n.checked_mul(snapshot_bound.checked_mul(5*extra+2*target+10)
        .and_then(|v| v.checked_add(6*target*target)).ok_or("Ritz work estimate overflow")?)
        .ok_or("Ritz work estimate overflow")?;
    if work > MAX_SCALAR_WORK {
        return Err(format!("Ritz snapshot work estimate {work} exceeds {MAX_SCALAR_WORK}; use fewer explicit source modes, retained modes, bridge ports or frequency samples"));
    }
    // Refuse exact undamped poles even when a particular port vanishes there.
    // Adding a damping floor would change the supplied mechanical model.
    for &hz in &options.sample_hz {
        let omega = TAU*hz;
        let squared = omega*omega;
        if !squared.is_finite() || squared == 0.0 {
            return Err("Ritz sample angular frequency is not representable".into());
        }
        if lambda.iter().zip(damping).any(|(&k,&c)| k == squared && c == 0.0) {
            return Err(format!("Ritz sample {hz} Hz is an exact undamped source pole"));
        }
    }
    let mut snapshots = Vec::with_capacity(snapshot_bound);
    let mut push = |mut v: Vec<f64>| -> Result<(), String> {
        if normalize(&mut v)? { snapshots.push(v); }
        Ok(())
    };
    for port in ports {
        push(port.iter().zip(lambda).map(|(b,k)| b/k).collect())?;
        for &hz in &options.sample_hz {
            let omega = TAU*hz;
            let mut real = vec![0.0;n];
            let mut imaginary = vec![0.0;n];
            for i in 0..n {
                let a = lambda[i]-omega*omega;
                let b = omega*damping[i];
                let scale = a.abs().max(b.abs());
                if !scale.is_finite() || scale == 0.0 {
                    return Err("Ritz harmonic response denominator is not representable".into());
                }
                let ar = a/scale;let br = b/scale;
                let denominator = ar*ar+br*br;
                // e^(-i omega t): inverse of (lambda-omega^2-i omega*c).
                real[i] = (port[i]*(ar/denominator))/scale;
                imaginary[i] = (port[i]*(br/denominator))/scale;
            }
            push(real)?;push(imaginary)?;
        }
    }
    let mut residuals = snapshots.clone();
    for v in &mut residuals { v[..options.keep_low_modes].fill(0.0); }
    let mut selected = vec![false;snapshots.len()];
    let mut columns = Vec::with_capacity(target);
    for i in 0..options.keep_low_modes {
        let mut q = vec![0.0;n];q[i] = 1.0;columns.push(q);
    }
    // Relative to each FULL normalized snapshot, before removing low modes.
    // The final measured residual includes every numerically dependent tail.
    let rank_tolerance = 64.0*f64::EPSILON*n as f64;
    while columns.len() < target {
        let mut pivot = None;let mut largest = rank_tolerance*rank_tolerance;
        for (i,v) in residuals.iter().enumerate() {
            if selected[i] { continue; }
            let squared = dot(v,v);
            if squared > largest { largest = squared;pivot = Some(i); }
        }
        let Some(pivot) = pivot else { break; };
        selected[pivot] = true;
        let mut q = residuals[pivot].clone();
        for _ in 0..2 { for old in &columns { remove(&mut q,old); } }
        // At most one pivot attempt per available column, including a final
        // numerically dependent strongest candidate: preserve the work bound.
        if dot(&q,&q) <= rank_tolerance*rank_tolerance { break; }
        normalize(&mut q)?;
        for (i,v) in residuals.iter_mut().enumerate() {
            if !selected[i] { for _ in 0..2 { remove(v,&q); } }
        }
        columns.push(q);
    }
    if columns.is_empty() { return Err("Ritz reduction has no nonzero bridge response or protected low mode".into()); }
    let mut max_relative_snapshot_error = 0.0_f64;
    for snapshot in &snapshots {
        let mut residual = snapshot.clone();
        for q in &columns {
            let coefficient = dot(snapshot,q);
            for (x,b) in residual.iter_mut().zip(q) { *x -= coefficient*b; }
        }
        let error = (dot(&residual,&residual)/dot(snapshot,snapshot)).sqrt();
        max_relative_snapshot_error = max_relative_snapshot_error.max(error);
    }
    let r = columns.len();
    let mut stiffness = vec![0.0;r*r];let mut projected_damping = vec![0.0;r*r];
    for row in 0..r { for col in row..r {
        let mut k = 0.0;let mut c = 0.0;
        for i in 0..n {
            let product = columns[row][i]*columns[col][i];
            k += lambda[i]*product;c += damping[i]*product;
        }
        if !k.is_finite() || !c.is_finite() {
            return Err("projected Ritz stiffness or damping overflow".into());
        }
        stiffness[row*r+col] = k;stiffness[col*r+row] = k;
        projected_damping[row*r+col] = c;projected_damping[col*r+row] = c;
    } }
    Ok(RitzBasis { columns, stiffness, damping: projected_damping,
        max_relative_snapshot_error, snapshot_count: snapshots.len() })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(max_modes: usize, keep_low_modes: usize) -> (Vec<f64>,Vec<f64>,Vec<Vec<f64>>,RitzOptions) {
        (vec![1.,4.,9.,16.,25.,36.],vec![0.11,0.07,0.43,0.2,0.31,0.6],
            vec![vec![0.7,-0.4,1.2,0.3,-0.8,0.9]],
            RitzOptions { max_modes,keep_low_modes,sample_hz:vec![0.31] })
    }
    fn coordinates(columns: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
        columns.iter().map(|q|dot(q,x)).collect()
    }
    fn quadratic(a: &[f64], x: &[f64]) -> f64 {
        (0..x.len()).flat_map(|r| (0..x.len()).map(move |c|x[r]*a[r*x.len()+c]*x[c])).sum()
    }

    #[test]
    fn static_and_harmonic_bridge_responses_satisfy_the_projected_physical_equations() {
        let (lambda,c,ports,options) = fixture(4,1);
        let basis = bridge_basis(&lambda,&c,&ports,&options).unwrap();
        assert_eq!(basis.columns.len(),4);assert_eq!(basis.snapshot_count,3);
        assert!(basis.max_relative_snapshot_error < 1e-13);
        let r = basis.columns.len();
        let force = coordinates(&basis.columns,&ports[0]);
        let static_x: Vec<_> = ports[0].iter().zip(&lambda).map(|(b,k)| b/k).collect();
        let static_q = coordinates(&basis.columns,&static_x);
        let omega = TAU*options.sample_hz[0];
        // Independent direct complex division, safe on this small SI fixture.
        let real: Vec<_> = (0..lambda.len()).map(|i| {
            let a=lambda[i]-omega*omega;let b=omega*c[i];
            ports[0][i]*a/(a*a+b*b)
        }).collect();
        let imaginary: Vec<_> = (0..lambda.len()).map(|i| {
            let a=lambda[i]-omega*omega;let b=omega*c[i];
            ports[0][i]*b/(a*a+b*b)
        }).collect();
        let qr = coordinates(&basis.columns,&real);let qi = coordinates(&basis.columns,&imaginary);
        for row in 0..r {
            let static_force: f64 = (0..r).map(|col|basis.stiffness[row*r+col]*static_q[col]).sum();
            let dynamic_real: f64 = (0..r).map(|col| {
                let a=basis.stiffness[row*r+col]-if row==col {omega*omega} else {0.};
                a*qr[col]+omega*basis.damping[row*r+col]*qi[col]
            }).sum();
            let dynamic_imaginary: f64 = (0..r).map(|col| {
                let a=basis.stiffness[row*r+col]-if row==col {omega*omega} else {0.};
                a*qi[col]-omega*basis.damping[row*r+col]*qr[col]
            }).sum();
            assert!((static_force-force[row]).abs()<2e-12);
            assert!((dynamic_real-force[row]).abs()<2e-12);
            assert!(dynamic_imaginary.abs()<2e-12);
        }
    }

    #[test]
    fn protected_low_modes_and_full_energy_forms_survive_rank_selection() {
        let (lambda,c,ports,options) = fixture(4,2);
        let basis = bridge_basis(&lambda,&c,&ports,&options).unwrap();
        assert_eq!(basis.columns.len(),4);
        for low in 0..2 {
            assert_eq!(basis.columns[low],(0..6).map(|i|if i==low {1.} else {0.}).collect::<Vec<_>>());
            assert!(basis.columns[2..].iter().all(|q|q[low]==0.));
        }
        for (i,a) in basis.columns.iter().enumerate() { for (j,b) in basis.columns.iter().enumerate() {
            assert!((dot(a,b)-if i==j {1.} else {0.}).abs()<1e-14);
        } }
        let q = [0.3,-0.2,0.7,-0.4];
        let full: Vec<_> = (0..6).map(|i|basis.columns.iter().zip(q).map(|(column,x)|column[i]*x).sum::<f64>()).collect();
        assert!((dot(&full,&full)-dot(&q,&q)).abs()<1e-14);
        let expected_k: f64 = full.iter().zip(&lambda).map(|(x,k)|x*x*k).sum();
        let expected_c: f64 = full.iter().zip(&c).map(|(x,c)|x*x*c).sum();
        assert!((quadratic(&basis.stiffness,&q)-expected_k).abs()<1e-13);
        assert!((quadratic(&basis.damping,&q)-expected_c).abs()<1e-14);
        assert!(basis.damping[2*4+3].abs()>1e-6,"heterogeneous physical damping must not become diagonal");
        assert!(expected_c>=0.);
        for i in 0..4 {for j in 0..4 {
            assert_eq!(basis.stiffness[i*4+j],basis.stiffness[j*4+i]);
            assert_eq!(basis.damping[i*4+j],basis.damping[j*4+i]);
        }}
    }

    #[test]
    fn pivots_are_repeatable_with_stable_ties_and_no_invented_rank() {
        let lambda=[1.,4.,9.,16.];let c=[0.;4];
        let ports=vec![vec![1.,0.,0.,0.],vec![0.,4.,0.,0.]];
        let options=RitzOptions::parse("1,0,0.1").unwrap();
        let a=bridge_basis(&lambda,&c,&ports,&options).unwrap();
        let b=bridge_basis(&lambda,&c,&ports,&options).unwrap();
        assert_eq!(a.columns,vec![vec![1.,0.,0.,0.]]);
        assert_eq!(a.columns,b.columns);assert_eq!(a.stiffness,b.stiffness);
        assert_eq!(a.damping,b.damping);assert_eq!(a.max_relative_snapshot_error,1.);
        let rank=bridge_basis(&lambda,&c,&ports[..1],&RitzOptions::parse("4,0,0.1").unwrap()).unwrap();
        assert_eq!(rank.columns.len(),1,"do not fill unused capacity with arbitrary modes");
        let (lambda,c,ports,options)=fixture(2,1);
        let compressed=bridge_basis(&lambda,&c,&ports,&options).unwrap();
        assert!(compressed.max_relative_snapshot_error>1e-3 && compressed.max_relative_snapshot_error<=1.+1e-14);
    }

    #[test]
    fn parser_and_budget_admission_are_bounded_and_explicit() {
        let a=RitzOptions::parse("128,64,100,1000,2200").unwrap();
        assert_eq!(a.max_modes,128);assert_eq!(a.keep_low_modes,64);
        for text in ["", "4,1", "0,0,100", "129,0,100", "4,5,100", "4,-1,100",
            "4,0,0", "4,0,NaN", "4,0,inf", "4,0,200,100", "4,0,100,100"] {
            assert!(RitzOptions::parse(text).is_err(),"{text}");
        }
        let too_many=format!("128,0,{}",(1..=17).map(|n|n.to_string()).collect::<Vec<_>>().join(","));
        assert!(RitzOptions::parse(&too_many).is_err());
        let lambda: Vec<_>=(1..=MAX_SOURCE_MODES).map(|i|i as f64).collect();
        let c=vec![0.1;lambda.len()];let ports=vec![vec![1.;lambda.len()];MAX_PORTS];
        let expensive=RitzOptions {max_modes:128,keep_low_modes:0,sample_hz:(1..=16).map(f64::from).collect()};
        assert!(bridge_basis(&lambda,&c,&ports,&expensive).unwrap_err().contains("work estimate"));
        let admitted=bridge_basis(&lambda,&c,&ports[..1],&RitzOptions::parse("1,1,0.1").unwrap()).unwrap();
        assert_eq!(admitted.columns[0].len(),MAX_SOURCE_MODES);
        let too_many=vec![1.;MAX_SOURCE_MODES+1];
        assert!(bridge_basis(&too_many,&too_many,&[too_many.clone()],&a).is_err());
        let full=bridge_basis(&vec![1.;128],&vec![0.1;128],&[vec![0.;128]],
            &RitzOptions::parse("128,128,0.1").unwrap()).unwrap();
        assert_eq!(full.columns.len(),128);assert_eq!(full.snapshot_count,0);
        assert_eq!(full.max_relative_snapshot_error,0.);
    }

    #[test]
    fn exact_undamped_poles_and_invalid_physics_refuse_without_added_loss() {
        let options=RitzOptions::parse("1,0,1").unwrap();
        let lambda=[TAU*TAU];let ports=vec![vec![1.]];
        assert!(bridge_basis(&lambda,&[0.],&ports,&options).unwrap_err().contains("exact undamped source pole"));
        let damped=bridge_basis(&lambda,&[0.003],&ports,&options).unwrap();
        assert_eq!(damped.damping,vec![0.003]);
        assert_eq!(damped.stiffness,lambda.to_vec());
        for bad in [vec![-1.],vec![f64::NAN],vec![f64::INFINITY]] {
            assert!(bridge_basis(&lambda,&bad,&ports,&options).is_err());
        }
        assert!(bridge_basis(&[-1.],&[0.1],&ports,&options).is_err());
        assert!(bridge_basis(&lambda,&[0.1],&[],&options).is_err());
        assert!(bridge_basis(&lambda,&[0.1],&vec![vec![1.];89],&options).is_err());
        assert!(bridge_basis(&lambda,&[0.1],&[vec![0.]],&options).is_err());
    }
}
