//! Orthonormal acoustic ports in the complete mass-normalized board basis.
//!
//! Rows are Q^T. The acoustic velocity is Q^T v and its reaction is Q f, so
//! mechanical power is unchanged. Orthonormality also makes acoustic-port
//! kinetic energy the exact selected part of the full board kinetic energy.
//! This is a coordinate map; only comparison with the complete physical load
//! can establish whether a truncated acoustic space is adequate.
use fs_math::c64::C64;
use fs_couple::render::plate::impact::radiation::{MAX_PORTS, fit::Fit};

#[derive(Clone, Debug)]
pub struct PortBasis {
    board_ports: usize,
    vectors: Vec<Vec<f64>>,
    pub(super) orthogonality_error: f64,
}

impl PortBasis {
    /// Admit a supplied basis without normalizing, deleting or reordering rows.
    pub fn new(board_ports: usize, vectors: Vec<Vec<f64>>) -> Result<Self, String> {
        let ports = vectors.len();
        if !(1..=crate::linear::MAX_BOARD_MODES).contains(&board_ports)
            || !(1..=32).contains(&ports) || ports > board_ports
            || vectors.iter().any(|row| row.len() != board_ports
                || row.iter().any(|v| !v.is_finite())) {
            return Err("radiation basis requires 1..32 complete rows over 1..128 board coordinates".into());
        }
        let tolerance = 256. * f64::EPSILON * (board_ports + ports + 1) as f64;
        let mut orthogonality_error = 0.0_f64;
        for (i, a) in vectors.iter().enumerate() {
            for (j, b) in vectors.iter().enumerate().take(i + 1) {
                let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
                let error = (dot - if i == j { 1. } else { 0. }).abs();
                if !error.is_finite() || error > tolerance {
                    return Err("radiation basis is not orthonormal in the loaded board energy coordinates".into());
                }
                orthogonality_error = orthogonality_error.max(error);
            }
        }
        Ok(Self { board_ports, vectors, orthogonality_error })
    }

    pub fn board_ports(&self) -> usize { self.board_ports }
    pub fn ports(&self) -> usize { self.vectors.len() }
    pub fn vectors(&self) -> &[Vec<f64>] { &self.vectors }

    /// Q^T Z Q, retaining signed complex cross-port terms.
    pub fn project_impedance(&self, impedance: &[C64]) -> Result<Vec<C64>, String> {
        let n = self.board_ports;
        let r = self.ports();
        self.check_matrix(impedance, n)?;
        let mut right = vec![C64::ZERO; n * r];
        for i in 0..n {
            for (a, row) in self.vectors.iter().enumerate() {
                for (j, &q) in row.iter().enumerate() {
                    right[i * r + a] = right[i * r + a] + impedance[i * n + j].scale(q);
                }
            }
        }
        let mut result = vec![C64::ZERO; r * r];
        for (a, row) in self.vectors.iter().enumerate() {
            for b in 0..r {
                for (i, &q) in row.iter().enumerate() {
                    result[a * r + b] = result[a * r + b] + right[i * r + b].scale(q);
                }
            }
        }
        self.check_matrix(&result, r)?;
        Ok(result)
    }

    /// Q Z Q^T in the original complete board coordinates. This lift preserves
    /// reciprocity and positive-realness, but cannot restore an omitted field.
    pub fn lift_impedance(&self, impedance: &[C64]) -> Result<Vec<C64>, String> {
        let n = self.board_ports;
        let r = self.ports();
        self.check_matrix(impedance, r)?;
        let mut left = vec![C64::ZERO; n * r];
        for i in 0..n {
            for b in 0..r {
                for (a, row) in self.vectors.iter().enumerate() {
                    left[i * r + b] = left[i * r + b] + impedance[a * r + b].scale(row[i]);
                }
            }
        }
        let mut result = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                for (b, row) in self.vectors.iter().enumerate() {
                    result[i * n + j] = result[i * n + j] + left[i * r + b].scale(row[j]);
                }
            }
        }
        self.check_matrix(&result, n)?;
        Ok(result)
    }

    fn check_matrix(&self, matrix: &[C64], n: usize) -> Result<(), String> {
        if matrix.len() != n * n
            || matrix.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err("radiation projection requires a complete finite impedance matrix".into());
        }
        Ok(())
    }
}

/// Select radiating combinations from geometry alone, before seeing any BEM
/// frequencies or receiver positions. Rows of Q^T are Euclidean-orthonormal in
/// the bank's mass-normalized board coordinates. The area-weighted Gram picks
/// the dominant normal-velocity fields, not a new structural eigenspace.
/// Each geometry owner supplies its existing panel ceiling.
pub(crate) fn geometry_basis(weights:&[Vec<f64>],areas:&[f64],maximum:usize,
    maximum_panels:usize)->Result<(PortBasis,f64),String> {
    let n=weights.len();let panels=areas.len();
    if !(1..=crate::linear::MAX_BOARD_MODES).contains(&n)
        || !(1..=MAX_PORTS).contains(&maximum)
        || !(1..=fs_bem::helmholtz::MAX_DENSE_PANELS).contains(&maximum_panels)
        || !(1..=maximum_panels).contains(&panels)
        || areas.iter().any(|a|!a.is_finite() || *a<=0.)
        || weights.iter().any(|row|row.len()!=panels || row.iter().any(|v|!v.is_finite())) {
        return Err("radiation port selection needs complete finite board fields, positive panel areas and 1..32 ports".into());
    }
    let mut weighted=vec![vec![0.;panels];n];let mut scale=0.0_f64;
    for (out,row) in weighted.iter_mut().zip(weights) {
        for ((value,&shape),&area) in out.iter_mut().zip(row).zip(areas) {
            *value=shape*area.sqrt();
            if !value.is_finite() || (shape!=0. && *value==0.) {
                return Err("radiation surface weighting exceeds finite numerical range".into());
            }
            scale=scale.max(value.abs());
        }
    }
    if scale==0. {
        // A physically stationary exterior has exactly zero load. Keep one
        // named coordinate for the shared fitter's existing zero-load image.
        let mut row=vec![0.;n];row[0]=1.;
        return Ok((PortBasis::new(n,vec![row])?,0.));
    }
    for value in weighted.iter_mut().flatten() {
        let old=*value;*value/=scale;
        if old!=0. && *value==0. {
            return Err("radiation surface scale hides a nonzero field".into());
        }
    }
    let mut gram=vec![0.;n*n];let mut identity=vec![0.;n*n];
    for i in 0..n {
        identity[i*n+i]=1.;
        for j in 0..=i {
            let value=weighted[i].iter().zip(&weighted[j]).map(|(a,b)|a*b).sum();
            gram[i*n+j]=value;gram[j*n+i]=value;
        }
    }
    let trace=(0..n).map(|i|gram[i*n+i]).sum::<f64>();
    let mut pairs=fs_modal::eigh_gen_dense(&gram,&identity,n).map_err(|e|e.to_string())?;
    let tolerance=256.*f64::EPSILON*n as f64*trace;
    if !trace.is_finite() || trace<=0. || pairs.len()!=n
        || pairs.iter().any(|p|!p.lambda.is_finite() || p.lambda< -tolerance
            || !p.residual.is_finite() || p.residual>tolerance
            || p.phi.len()!=n || p.phi.iter().any(|v|!v.is_finite())) {
        return Err("radiation surface Gram eigenproblem did not resolve its positive semidefinite field".into());
    }
    // Stable ordering retains the eigensolver's order for equal eigenvalues.
    // No pressure sample or fit failure changes this geometry-only selection.
    pairs.sort_by(|a,b|b.lambda.total_cmp(&a.lambda));
    let mut vectors=Vec::new();
    for pair in pairs.into_iter().filter(|p|p.lambda>tolerance).take(maximum) {
        let mut row=pair.phi;let mut pivot=0;
        for i in 1..n {if row[i].abs()>row[pivot].abs() {pivot=i;}}
        if row[pivot]<0. {for value in &mut row {*value= -*value;}}
        vectors.push(row);
    }
    if vectors.is_empty() {return Err("radiation surface has no resolved nonzero port".into());}
    let basis=PortBasis::new(n,vectors)?;
    let mut error=0.0_f64;let mut reference=0.0_f64;
    let mut projected=vec![0.;n];
    for panel in 0..panels {
        projected.fill(0.);
        for row in basis.vectors() {
            let amplitude=row.iter().zip(&weighted).map(|(q,w)|q*w[panel]).sum::<f64>();
            for (value,q) in projected.iter_mut().zip(row) {*value+=q*amplitude;}
        }
        for (value,row) in projected.iter().zip(&weighted) {
            reference=reference.hypot(row[panel]);error=error.hypot(value-row[panel]);
        }
    }
    let surface_rms_error=error/reference;
    if !surface_rms_error.is_finite() {return Err("nonfinite radiation surface projection error".into());}
    Ok((basis,surface_rms_error))
}

/// Apply the same whole-matrix tolerances as the existing shared fit to the
/// lifted response. Both even training and odd held-out samples are checked;
/// no held-out value selects the geometry basis or any fit coefficient.
pub(crate) fn validate_complete_fit(omega:&[f64],data:&[Vec<C64>],basis:&PortBasis,
    fitted:&mut Fit)->Result<(),String> {
    let ports=basis.board_ports();
    if data.len()!=omega.len() || data.iter().any(|row|row.len()!=ports*ports
        || row.iter().any(|z|!z.re.is_finite() || !z.im.is_finite())) {
        return Err("projected radiation validation requires the complete original BEM matrices".into());
    }
    let mut peak=0.0_f64;let mut rms=0.0_f64;let mut rp=0.0_f64;let mut rr=0.0_f64;
    for offset in 0..2 {
        let mut max_reference=0.0_f64;let mut max_error=0.0_f64;
        let mut references=0.0_f64;let mut errors=0.0_f64;
        let mut real_max=0.0_f64;let mut real_error=0.0_f64;
        let mut real_references=0.0_f64;let mut real_errors=0.0_f64;
        for f in (offset..omega.len()).step_by(2) {
            let lifted=basis.lift_impedance(&fitted.model.impedance(omega[f])?)?;
            let error:Vec<_>=lifted.iter().zip(&data[f]).map(|(a,b)|*a-*b).collect();
            let reference=data[f].iter().fold(0.0_f64,|s,z|s.hypot(z.abs()));
            let delta=error.iter().fold(0.0_f64,|s,z|s.hypot(z.abs()));
            max_reference=max_reference.max(reference);max_error=max_error.max(delta);
            references=references.hypot(reference);errors=errors.hypot(delta);
            let mut h=0.0_f64;let mut dh=0.0_f64;
            for i in 0..ports {for j in 0..ports {
                h=h.hypot((data[f][i*ports+j]+data[f][j*ports+i].conj()).scale(0.5).abs());
                dh=dh.hypot((error[i*ports+j]+error[j*ports+i].conj()).scale(0.5).abs());
            }}
            real_max=real_max.max(h);real_error=real_error.max(dh);
            real_references=real_references.hypot(h);real_errors=real_errors.hypot(dh);
        }
        let ratio=|error:f64,reference:f64|if error==0. {0.} else {error/reference};
        peak=peak.max(ratio(max_error,max_reference));rms=rms.max(ratio(errors,references));
        rp=rp.max(ratio(real_error,real_max));rr=rr.max(ratio(real_errors,real_references));
    }
    if ![peak,rms,rp,rr].iter().all(|x|x.is_finite()) || peak>0.05 || rms>0.02 || rp>0.10 || rr>0.05 {
        return Err(format!("projected BEM load fit refused on the complete board: complex peak/RMS={peak:.5}/{rms:.5}, power peak/RMS={rp:.5}/{rr:.5}; limits .05/.02/.10/.05; supply more radiation ports or retain the complete load"));
    }
    fitted.peak_error=peak;fitted.rms_error=rms;
    fitted.resistance_peak_error=rp;fitted.resistance_rms_error=rr;
    Ok(())
}
