///robust sum
pub fn robust_sum(e: &[f64], f: &[f64]) -> Vec<f64> {
    linear_expansion(e, f)
}

///linear expansion sum
fn linear_expansion(e: &[f64], f: &[f64]) -> Vec<f64> {
    let ne = e.len();
    let nf = f.len();
    if ne == 1 && nf == 1 {
        return scalar_scalar(e[0], f[0]);
    }
    let n = ne + nf;
    let mut count = 0usize;
    let mut eptr  = 0usize;
    let mut fptr  = 0usize;
    let mut g = vec![0f64; n];
    let mut ei = e[eptr];
    let mut ea = ei.abs();
    let mut fi = f[fptr];
    let mut fa = fi.abs();
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
    let (mut _x, mut _bv, mut _av, mut _br, mut _ar): (f64, f64, f64, f64, f64);
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

    g[0..count].to_vec()
}

///scalar sum: easy case: add two scalars
fn scalar_scalar(a: f64, b: f64) -> Vec<f64> {
    let x = a + b;
    let bv = x - a;
    let av = x - bv;
    let br = b - bv;
    let ar = a - av;
    let y = ar + br;
    if y != 0.0 {
        return vec!(y, x);
    }
    vec!(x)
}


#[cfg(test)]
mod robust_sum_test {
    extern crate rand;
    extern crate validate_robust_seq;

    use self::rand::random;
    use super::robust_sum;
    use self::validate_robust_seq::validate_sequence as validate;


    fn rnd() -> f64 {
        random::<f64>()
    }

    fn cmp(va: Vec<f64>, vb: Vec<f64>) -> bool {
        (va.len() == vb.len()) && va.iter().zip(vb).all(|(a, b)| f64_same_val(*a, b))
    }

    fn f64_same_val(a: f64, b: f64) -> bool {
        a == b
        //(a.is_nan() && b.is_nan()) || (a == b)
    }

    #[test]
    fn test_robust_sum() {
        assert_eq!(robust_sum(&vec!(1., 64.0), &vec!(-1e-64, 1e64)), [-1e-64, 65., 1e64]);
        assert_eq!(robust_sum(&vec!(0.), &vec!(0.)), [0.]);
        assert_eq!(robust_sum(&vec!(0.), &vec!(1.)), [1.]);
        assert_eq!(robust_sum(&vec!(1., 1e64), &vec!(1e-64, 2.)), [1e-64, 3., 1e64]);
        assert_eq!(robust_sum(&vec!(1.), &vec!(1e-64, 1e-16)), [1e-64, 1e-16, 1.]);

        for i in -10..(10 + 1) {
            for j in -10..(10 + 1) {
                assert_eq!(robust_sum(&vec!(i as f64), &vec!(j as f64)), [(i + j) as f64]);
            }
        }

        assert!(validate(
            &robust_sum(
                &vec!(5.711861227349496e-133, 1e-116),
                &vec!(5.711861227349496e-133, 1e-116),
            )
        ));


        let mut nois = Vec::<f64>::with_capacity(10);
        let mut expect = Vec::<f64>::with_capacity(10);
        for i in 0..10 {
            nois.push(2f64.powi(-1000 + 53 * i));
            expect.push(2f64.powi(-999 + 53 * i));
        }
        let x = robust_sum(&nois, &nois);
        assert!(cmp(x.clone(), expect));
        assert!(validate(&x));

        assert!(cmp(robust_sum(&vec!(0.), &vec!(1., 1e64)), vec!(1., 1e64)));

        let mut s = vec![0.0];
        for _ in 0..1000 {
            s = robust_sum(&s, &vec![rnd() * 2f64.powf(rnd() * 1800. - 900.)]);
            assert!(validate(&s))
        }
    }
}

