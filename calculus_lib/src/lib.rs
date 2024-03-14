use either::Either;
use num::complex::Complex64;

pub fn newton_secant(x_prev: f32, x_curr: f32, f: fn(f32) -> f32) -> f32 {
    x_curr - f(x_curr) * (x_curr - x_prev) / (f(x_curr) - f(x_prev))
}

pub const EPS: f64 = 1e-3;

pub type Solution = Either<f64, Complex64>;

pub fn solve_linear_equation(a: f64, b: f64) -> Solution {
    Either::Left(-b / a)
}

pub fn solve_quadratic_equation(a: f64, b: f64, c: f64) -> (Solution, Solution) {
    let d = b * b - 4.0 * a * c;

    match d {
        d if d > 0.0 => {
            let x1 = (-b + d.sqrt()) / (2.0 * a);
            let x2 = (-b - d.sqrt()) / (2.0 * a);
            (Either::Left(x1), Either::Left(x2))
        }
        d if d < 0.0 => {
            let x1 = Complex64::new(-b / (2.0 * a), (-d).sqrt() / (2.0 * a));
            let x2 = Complex64::new(-b / (2.0 * a), -(-d).sqrt() / (2.0 * a));
            (Either::Right(x1), Either::Right(x2))
        }
        _ => {
            let x = -b / (2.0 * a);
            (Either::Left(x), Either::Left(x))
        }
    }
}

// p(x) = a[0] * x^n + a[1] * x^(n-1) + ... + a[n-1] * x + a[n]
pub fn solve_polynomial(a: Vec<f64>, x: f64) -> f64 {
    a.iter().enumerate().fold(0.0, |acc, (i, &ai)| {
        acc + ai * x.powi((a.len() - 1 - i) as i32)
    })
}
