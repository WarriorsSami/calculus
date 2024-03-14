use calculus_lib::{solve_polynomial, EPS};

fn main() {
    // bernoulli method for finding maximal root of a polynomial
    let n = 5;
    let a = [1.0, 5.0, 0.0, 0.0, 0.0, -5.0];
    let mut y = vec![0.0; n];

    const ITMAX: usize = 6;

    y[0] = n as f64;
    y[1] = -a[1] / a[0];

    for i in 2..n {
        y[i] = (1..=(i - 1)).fold(0.0, |acc, j| acc - y[j] * a[i - j] / a[0]);
        y[i] -= i as f64 * a[i] / a[0];
    }

    let mut i = 0;
    let mut m = 0;
    let mut x;

    loop {
        y.push((1..=n).fold(0.0, |acc, j| acc - y[n + i - j] * a[j] / a[0]));
        x = y[n + i] / y[n + i - 1];

        println!("y_{} = {}", i, y[n + i]);
        println!("x_{} = {}", i, x);

        m += 1;
        i += 1;

        let p_x = solve_polynomial(a.to_vec(), x);
        println!("p(x_{}): {}", i, p_x);

        if p_x < EPS || m > ITMAX {
            break;
        }
    }

    match m > ITMAX {
        true => println!("No convergence in {} iterations", ITMAX),
        false => println!("Maximal root: {}", x),
    }
}
