use calculus_lib::{Solution, EPS};

fn main() {
    // bairstow method for finding roots of a polynomial
    let mut n = 4;
    let mut a = [1.0, 1.0, -10.0, -34.0, -26.0];

    let p0 = 0.1;
    let q0 = 0.0;

    let mut p = p0;
    let mut q = q0;

    let mut xs: Vec<Solution> = vec![];

    while n >= 3 {
        let mut b: Vec<f64> = vec![0.0; n + 1];

        loop {
            b[0] = a[0];
            b[1] = a[1] - p * b[0];

            for i in 2..=n {
                b[i] = a[i] - p * b[i - 1] - q * b[i - 2];
            }

            let mut c = vec![0.0; n];
            c[0] = b[0];
            c[1] = b[1] - p * c[0];

            for i in 2..=n - 1 {
                c[i] = b[i] - p * c[i - 1] - q * c[i - 2];
            }

            let d = c[n - 2] * c[n - 2] - c[n - 3] * c[n - 1] + c[n - 3] * b[n - 1];
            let dp = (-b[n - 1]) * c[n - 2] + b[n] * c[n - 3];
            let dq = (-b[n]) * c[n - 2] + b[n - 1] * c[n - 1] - b[n - 1] * b[n - 1];

            p -= dp / d;
            q -= dq / d;

            if f64::max(b[n - 1].abs(), (b[n] + p * b[n - 1]).abs()) < EPS {
                break;
            }
        }

        // solve the quadratic equation x^2 + px + q = 0
        let (x1, x2) = calculus_lib::solve_quadratic_equation(1.0, p, q);
        xs.push(x1);
        xs.push(x2);

        n -= 2;

        a[..(n + 1)].copy_from_slice(&b[..(n + 1)]);
    }

    match n {
        2 => {
            let (x1, x2) = calculus_lib::solve_quadratic_equation(a[0], a[1], a[2]);
            xs.push(x1);
            xs.push(x2);
        }
        1 => {
            let x = calculus_lib::solve_linear_equation(a[0], a[1]);
            xs.push(x);
        }
        _ => {}
    }

    for (idx, x) in xs.iter().enumerate() {
        match x {
            Solution::Left(x) => println!("x{} = {}", idx + 1, x),
            Solution::Right(x) => println!("x{} = {}", idx + 1, x),
        }
    }
}
