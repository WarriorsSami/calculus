use calculus_lib::{Matrix, EPS};

fn main() {
    let A = Matrix::new(vec![
        vec![-1.0, -0.1, -0.1],
        vec![-0.2, -1.0, -0.1],
        vec![-0.2, -0.2, -1.0],
    ]);

    let b = [0.0, -1.2, -1.3, -1.4];

    // Seidel-Gauss method
    let mut x = vec![0.0; A.n() + 1];

    let q = (1..=A.n())
        .map(|i| {
            (1..=A.n())
                .map(|j| if i == j { 0.0 } else { A.at(i, j).abs() })
                .sum::<f64>()
                / A.at(i, i).abs()
        })
        .fold(f64::NEG_INFINITY, f64::max);

    if q < 1.0 {
        let mut y = vec![0.0; A.n() + 1];
        for i in 1..=A.n() {
            let ay = (1..=i - 1).map(|j| A.at(i, j) * y[j]).sum::<f64>();
            let ax = (i + 1..=A.n()).map(|j| A.at(i, j) * x[j]).sum::<f64>();
            y[i] = (b[i] - ay - ax) / A.at(i, i);
        }

        let max_xy = (1..=A.n())
            .map(|i| (x[i] - y[i]).abs())
            .fold(f64::NEG_INFINITY, f64::max);

        // compute the lowest natural number k such that q^k / (1 - q) * max_xy < eps
        let k = ((1.0 - q) * EPS / max_xy).log(q).ceil() as i32;

        for _ in 2..=k {
            for i in 1..=A.n() {
                x[i] = y[i];
                let ay = (1..=i - 1).map(|j| A.at(i, j) * y[j]).sum::<f64>();
                let ax = (i + 1..=A.n()).map(|j| A.at(i, j) * x[j]).sum::<f64>();
                y[i] = (b[i] - ay - ax) / A.at(i, i);
            }
        }

        println!(
            "The solution after {} iterations is: {:?}",
            k,
            x.iter().skip(1).collect::<Vec<&f64>>()
        );
    } else {
        println!("The matrix A is not convergent");
    }
}
