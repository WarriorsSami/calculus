use calculus_lib::Matrix;

fn main() {
    let mut A = Matrix::new(vec![
        vec![2.0, 2.0, 3.0],
        vec![4.0, 5.0, 6.0],
        vec![1.0, 2.0, 2.0],
    ]);

    let b = [0.0, 7.0, 15.0, 5.0];

    // LR Doolittle factorization
    if A.at(1, 1) == 0.0 {
        println!("Cannot factorize matrix A");
    } else {
        for i in 2..=A.n() {
            A.set(i, 1, A.at(i, 1) / A.at(1, 1));
        }

        for k in 2..A.n() {
            for j in k..=A.n() {
                let mut sum = 0.0;
                for s in 1..k {
                    sum += A.at(k, s) * A.at(s, j);
                }
                A.set(k, j, A.at(k, j) - sum);
            }
            if A.at(k, k) == 0.0 {
                println!("Cannot factorize matrix A");
                return;
            } else {
                for i in (k + 1)..=A.n() {
                    let mut sum = 0.0;
                    for s in 1..k {
                        sum += A.at(i, s) * A.at(s, k);
                    }
                    A.set(i, k, (A.at(i, k) - sum) / A.at(k, k));
                }
            }
        }

        let mut sum = 0.0;
        for s in 1..A.n() {
            sum += A.at(A.n(), s) * A.at(s, A.n());
        }
        A.set(A.n(), A.n(), A.at(A.n(), A.n()) - sum);
    }

    println!("{}", A);

    // Solve system Ax = b

    let mut y = vec![0.0; A.n() + 1];
    let mut x = vec![0.0; A.n() + 1];

    y[1] = b[1];
    for i in 2..=A.n() {
        let mut sum = 0.0;
        for j in 1..i {
            sum += A.at(i, j) * y[j];
        }
        y[i] = b[i] - sum;
    }

    x[A.n()] = y[A.n()] / A.at(A.n(), A.n());
    for i in (1..A.n()).rev() {
        let mut sum = 0.0;
        for j in i + 1..=A.n() {
            sum += A.at(i, j) * x[j];
        }
        x[i] = (y[i] - sum) / A.at(i, i);
    }

    println!("Solution:");

    for i in 1..=A.n() {
        println!("x[{}] = {}", i, x[i]);
    }
}
