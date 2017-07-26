use  std::f64;

///robust sum
pub fn robust_sum(e:Vec<f64>, f:Vec<f64>) -> Vec<f64> {
	unimplemented!();
	return linear_expansion_sum(e, f)
}

///linear expansion sum
fn linear_expansion_sum(e:Vec<f64>, f:Vec<f64>) -> Vec<f64> {
	let ne = e.len();
	let nf = f.len();
	if ne == 1 && nf == 1 {
		return scalar_scalar(e[0], f[0]);
	}
	let n = ne + nf;
	let mut g = vec![0.0; n];
	let mut count = 0;
	let mut eptr  = 0;
	let mut fptr  = 0;
	let mut ei    = e[eptr];
	let mut ea    = ei.abs();
	let mut fi    = f[fptr];
	let mut fa    = fi.abs();
    let (mut a, mut b): (f64, f64);

	if ea < fa {
		b = ei;
		eptr += 1;
		if eptr < ne {
			ei = e[eptr];
			ea = ei.abs();
		}
	} else {
		b = fi;
		fptr += 1;
		if fptr < nf {
			fi = f[fptr];
			fa = fi.abs();
		}
	}
	if (eptr < ne && ea < fa) || (fptr >= nf) {
		a = ei;
		eptr += 1;
		if eptr < ne {
			ei = e[eptr];
			ea = ei.abs();
		}
	} else {
		a = fi;
		fptr += 1;
		if fptr < nf {
			fi = f[fptr];
			fa = fi.abs();
		}
	}
	let mut x = a + b;
	let mut bv = x - a;
	let mut y = b - bv;
	let mut q0 = y;
	let mut q1 = x;
	let (mut _x,mut _bv,mut _av,mut _br,mut _ar): (f64 ,f64 ,f64 ,f64 ,f64 );
	while eptr < ne && fptr < nf {
		if ea < fa {
			a = ei;
			eptr += 1;
			if eptr < ne {
				ei = e[eptr];
				ea = ei.abs();
			}
		} else {
			a = fi;
			fptr += 1;
			if fptr < nf {
				fi = f[fptr];
				fa = fi.abs();
			}
		}
		b = q0;
		x = a + b;
		bv = x - a;
		y = b - bv;
		if y != 0.0 {
			g[count] = y;
			count += 1;
		}
		_x = q1 + x;
		_bv = _x - q1;
		_av = _x - _bv;
		_br = x - _bv;
		_ar = q1 - _av;
		q0 = _ar + _br;
		q1 = _x;
	}
	while eptr < ne {
		a = ei;
		b = q0;
		x = a + b;
		bv = x - a;
		y = b - bv;
		if y != 0.0 {
			g[count] = y;
			count += 1;
		}
		_x = q1 + x;
		_bv = _x - q1;
		_av = _x - _bv;
		_br = x - _bv;
		_ar = q1 - _av;
		q0 = _ar + _br;
		q1 = _x;
		eptr += 1;
		if eptr < ne {
			ei = e[eptr];
		}
	}
	while fptr < nf {
		a = fi;
		b = q0;
		x = a + b;
		bv = x - a;
		y = b - bv;
		if y != 0.0 {
			g[count] = y;
			count += 1;
		}
		_x = q1 + x;
		_bv = _x - q1;
		_av = _x - _bv;
		_br = x - _bv;
		_ar = q1 - _av;
		q0 = _ar + _br;
		q1 = _x;
		fptr += 1;
		if fptr < nf {
			fi = f[fptr];
		}
	}
	if q0 != 0.0 {
		g[count] = q0;
		count += 1;
	}
	if q1 != 0.0 {
		g[count] = q1;
		count += 1;
	}
	if count == 0 {
		g[count] = 0.0;
		count += 1;
	}

    let mut gv = Vec::with_capacity(count);
    for i in 0..count {
        gv.push(g[i])
    }

	return gv
}

///scalar sum: easy case: add two scalars
fn scalar_scalar(a:f64, b:f64) -> Vec<f64> {
	let  x  = a + b;
	let  bv = x - a;
	let  av = x - bv;
	let  br = b - bv;
	let  ar = a - av;
	let  y  = ar + br;
	if y != 0.0 {
		return vec!(y, x)
	}
	return vec!(x)
}


#[cfg(test)]
mod robust_sum_test {
    use super::robust_sum;
    fn compare(a: &[f64], b: &[f64]) -> bool {
        a == b
    }
    #[test]
    fn test_robust_sum() {
        let a = robust_sum(vec![1.0, 64.0], vec![-1e-64, 1e64]);
        let b = vec![1.0, 1e+64];
		println!("{:?}", a)
//        assert!(compare(a.as_slice(), b.as_slice()));
    }
}

